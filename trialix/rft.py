"""Participant-level RFT using the fitted model's whitened residuals.

pyrft is optional and imported only when RFT is requested. Neither the observed
statistic map nor the preprocessing kernel is used to estimate null smoothness.
"""

from dataclasses import dataclass
import logging

import numpy as np
from scipy import special, stats

logger = logging.getLogger("trialix")


def require_pyrft():
    """Fail before fitting if the optional RFT backend is unavailable."""
    try:
        import pyrft
        from pyrft import estimate_smoothness_residuals, threshold_map
    except ImportError as exc:
        raise RuntimeError(
            "RFT requires pyrft 0.2 with its NIfTI adapter. From the Trialix "
            "checkout, install it with: python -m pip install -e '.[rft]'"
        ) from exc
    return pyrft, estimate_smoothness_residuals, threshold_map


@dataclass
class RFTContext:
    """One residual smoothness estimate shared by a run's T contrasts."""

    mask_img: object
    residual_df: float
    smoothness: dict
    preparation: dict


def prepare_rft_context(flm):
    """Assemble a full whitened residual stack in fitted-mask voxel order.

    Nilearn groups AR fits by estimated coefficients. Reconstruct the residual
    matrix using labels_ rather than assuming dictionary order is voxel order.
    Public flm.residuals is deliberately not used: the model is fitted in a
    whitened space. pyrft standardizes the raw whitened stack by its RSS/df.
    """
    _, estimate_smoothness, _ = require_pyrft()
    if any(len(getattr(flm, key, [])) != 1
           for key in ('design_matrices_', 'labels_', 'results_')):
        raise ValueError("Trialix RFT currently requires exactly one fitted run.")
    if getattr(flm, 'minimize_memory', True):
        raise ValueError("RFT needs full regression results (minimize_memory=False).")
    design = flm.design_matrices_[0]
    labels = np.asarray(flm.labels_[0])
    results = flm.results_[0]
    n_scans = len(design)
    if labels.ndim != 1 or not results or n_scans < 4:
        raise ValueError("RFT requires voxel labels, regression results, and at least four scans.")
    whitened = np.empty((n_scans, labels.size), dtype=float)
    assigned = np.zeros(labels.size, dtype=bool)
    degrees_of_freedom = []
    for label, result in results.items():
        selected = labels == label
        if not np.any(selected):
            continue
        residuals = getattr(result, 'whitened_residuals', None)
        if residuals is None:
            raise ValueError("Fitted regression results do not expose whitened_residuals.")
        residuals = np.asarray(residuals, dtype=float)
        if residuals.shape != (n_scans, int(selected.sum())):
            raise ValueError("Whitened residual shape does not match the fitted design/voxel labels.")
        whitened[:, selected] = residuals
        assigned[selected] = True
        degrees_of_freedom.append(float(result.df_residuals))
    if not np.all(assigned) or not degrees_of_freedom:
        raise ValueError("Some fitted voxels have no whitened residuals.")
    df = degrees_of_freedom[0]
    if (not np.all(np.isfinite(degrees_of_freedom)) or df <= 1 or df > n_scans
            or not np.allclose(degrees_of_freedom, df, rtol=0, atol=1e-8)):
        raise ValueError("RFT requires common finite residual degrees of freedom > 1 across AR groups.")
    mask_img = flm.masker_.mask_img_
    mask = np.asarray(mask_img.dataobj) != 0
    if mask.ndim != 3 or int(mask.sum()) != labels.size:
        raise ValueError("Fitted mask and regression voxel order are inconsistent.")
    residual_img = flm.masker_.inverse_transform(whitened)
    smoothness = estimate_smoothness(
        np.asarray(residual_img.dataobj, dtype=float), mask,
        residual_df=df, standardized=False, n_observations=n_scans,
    )
    preparation = {
        'Source': 'FirstLevelModel.results_[0][label].whitened_residuals',
        'NoiseModel': flm.noise_model,
        'Whitened': True,
        'Standardization': 'full whitened residual stack / sqrt(RSS / residual_df)',
        'ObservationCount': int(n_scans),
        'ResidualDegreesOfFreedom': df,
        'DesignRank': int(n_scans - df),
        'AutoregressiveGroupCount': len(degrees_of_freedom),
        'SearchMaskSource': 'FirstLevelModel.masker_.mask_img_',
        'SearchMaskVoxelCount': int(mask.sum()),
        'DerivativeBoundary': 'complete forward 2x2x2 cubes inside fitted mask',
    }
    logger.info("RFT residual smoothness: df=%s, scans=%s, FWHM(vox)=%s, derivative voxels=%s",
                df, n_scans, smoothness['fwhm_vox'], smoothness['derivative_voxels'])
    return RFTContext(mask_img, df, smoothness, preparation)


def threshold_rft_contrast(outputs, context, alpha):
    """Threshold the native T field and retain original Z scores for output.

    Two-sided height FWE uses alpha/2 per tail. No cluster inference, map-only
    smoothness approximation, or automatic Bonferroni fallback is performed.
    """
    pyrft, _, threshold_map = require_pyrft()
    if outputs.get('stat') is None or outputs.get('z_score') is None:
        raise ValueError("RFT requires both native T and Z-score contrast maps.")
    native = outputs['stat']
    z_img = outputs['z_score']
    if native.shape != z_img.shape or not np.allclose(native.affine, z_img.affine,
                                                     rtol=0, atol=1e-5):
        raise ValueError("Native statistic and Z-score maps do not share the same grid.")
    result = threshold_map(
        native, context.mask_img, stat='T', df=context.residual_df,
        alpha=alpha, tail='two-sided', method='rft',
        fwhm=context.smoothness['fwhm_vox'], fwhm_units='vox',
    )
    z_data = np.asarray(z_img.dataobj, dtype=float)
    mask = np.asarray(context.mask_img.dataobj) != 0
    if not np.all(np.isfinite(z_data[mask])):
        raise ValueError("Z-score map has nonfinite values inside the fitted search mask.")
    survives = np.asarray(result.thresholded_map.dataobj) != 0
    header = z_img.header.copy()
    header.set_data_dtype(np.float64)
    thresholded_z = z_img.__class__(np.where(survives, z_data, 0.), z_img.affine, header)
    # Same marginal T-to-Z transform as the model's Z output, for report display.
    z_threshold = float(-special.ndtri_exp(stats.t.logsf(result.threshold, context.residual_df)))
    if not np.isfinite(z_threshold):
        raise ValueError("RFT T threshold cannot be represented as a finite Z threshold.")
    metadata = dict(result.metadata)
    metadata['smoothness'] = dict(context.smoothness)
    metadata.update({
        'Backend': 'pyrft', 'BackendVersion': pyrft.__version__,
        'ResidualPreparation': dict(context.preparation),
        'SourceStatistic': 'T', 'OutputStatistic': 'Z',
        'TThreshold': float(result.threshold), 'ZThreshold': z_threshold,
        'Selection': 'retain original Z scores where native abs(T) >= TThreshold',
        'Scope': 'two-sided peak/height FWE per contrast within the fitted model mask',
    })
    logger.info("RFT alpha=%s: |T| >= %.6f (df=%s), equivalent |Z| >= %.6f; resels=%s",
                alpha, result.threshold, context.residual_df, z_threshold, metadata['resels'])
    return thresholded_z, z_threshold, metadata
