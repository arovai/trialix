"""First-level model building, fitting, and contrast computation."""

import logging

from nilearn.glm import threshold_stats_img
from nilearn.glm.first_level import FirstLevelModel

logger = logging.getLogger("trialix")


def build_first_level_model(t_r, subject_label=None, model_params=None,
                            debug=False):
    """
    Build a nilearn FirstLevelModel with the given parameters.

    Parameters
    ----------
    t_r : float
        Repetition time in seconds.
    subject_label : str or None
        Subject identifier.
    model_params : dict or None
        Model parameters (hrf_model, drift_model, drift_order, high_pass,
        smoothing_fwhm). Missing keys use nilearn defaults.
    debug : bool
        If True, increase nilearn verbosity.

    Returns
    -------
    flm : FirstLevelModel
        Configured (but not yet fitted) FirstLevelModel.
    """
    if model_params is None:
        model_params = {}

    # Set verbose level: 0=silent, 1=warnings, 2=info (only used for debug)
    verbose = 2 if debug else 0

    flm = FirstLevelModel(
        t_r=t_r,
        slice_time_ref=0.0,
        hrf_model=model_params.get('hrf_model', 'glover'),
        drift_model=model_params.get('drift_model', 'cosine'),
        high_pass=model_params.get('high_pass', 0.01),
        drift_order=model_params.get('drift_order', 1),
        smoothing_fwhm=model_params.get('smoothing_fwhm'),
        noise_model='ar1',
        standardize=False,
        signal_scaling=0,
        minimize_memory=False,
        verbose=verbose,
        n_jobs=1,
        subject_label=subject_label,
    )

    logger.info(f"FirstLevelModel configured: TR={t_r}s, "
                f"HRF={flm.hrf_model}, drift={flm.drift_model}, "
                f"high_pass={flm.high_pass}, "
                f"smoothing_fwhm={flm.smoothing_fwhm}")

    return flm


def fit_model(flm, bold_img, events_df, confounds_df=None, design_matrix=None):
    """
    Fit the FirstLevelModel.

    Parameters
    ----------
    flm : FirstLevelModel
        Configured model.
    bold_img : str or Nifti1Image
        Path to (or loaded) preprocessed BOLD data.
    events_df : pd.DataFrame
        Events with onset, duration, trial_type columns.
    confounds_df : pd.DataFrame or None
        Confound regressors.
    design_matrix : pd.DataFrame or None
        Complete design for a physically cropped image. When supplied, events
        and confounds are already incorporated and are not passed to fit.

    Returns
    -------
    flm : FirstLevelModel
        Fitted model.
    """
    logger.info("Fitting FirstLevelModel...")
    if design_matrix is None:
        flm.fit(bold_img, events=events_df, confounds=confounds_df)
    else:
        flm.fit(bold_img, design_matrices=design_matrix)
    logger.info("Model fitting complete.")

    # Log design matrix columns
    dm = flm.design_matrices_[0]
    logger.info(f"Design matrix columns: {list(dm.columns)}")

    return flm


def compute_contrasts(flm, contrasts):
    """
    Compute contrasts and return all output types.

    Parameters
    ----------
    flm : FirstLevelModel
        Fitted model.
    contrasts : list of str
        Contrast definition strings (e.g., "cond1 - cond2").

    Returns
    -------
    results : dict
        {contrast_name: {output_type: Nifti1Image}} dictionary.
        Output types: z_score, stat, p_value, effect_size, effect_variance.
    """
    results = {}

    for contrast_def in contrasts:
        logger.info(f"Computing contrast: '{contrast_def}'")

        outputs = flm.compute_contrast(contrast_def, output_type='all')
        results[contrast_def] = outputs

        logger.info(f"  -> Got maps: {list(outputs.keys())}")

    return results


def threshold_contrast_map(z_map, alpha=0.05, height_control=None):
    """
    Apply statistical thresholding to a z-score map.

    Parameters
    ----------
    z_map : Nifti1Image
        Z-score statistical map.
    alpha : float
        P-value threshold.
    height_control : str or None
        Correction method: 'fdr', 'bonferroni', or None (uncorrected).

    Returns
    -------
    thresholded_map : Nifti1Image
        Thresholded z-score map.
    threshold_value : float
        The z-scale threshold that was applied.
    """
    logger.info(f"Thresholding z-map: alpha={alpha}, "
                f"height_control={height_control}")

    thresholded_map, threshold_value = threshold_stats_img(
        z_map,
        alpha=alpha,
        height_control=height_control,
        cluster_threshold=0,
        two_sided=True,
    )

    logger.info(f"  -> Applied z-threshold: {threshold_value:.4f}")

    return thresholded_map, threshold_value


def generate_report(flm, contrasts, alpha=0.001, height_control='fpr'):
    """
    Generate an HTML report for the fitted model.

    Parameters
    ----------
    flm : FirstLevelModel
        Fitted model.
    contrasts : list of str
        Contrast definitions to include in the report.
    alpha : float
        Significance level for thresholding.
    height_control : str or None
        False-positive control method ('fpr', 'fdr', 'bonferroni', or None).

    Returns
    -------
    report : HTMLReport
        Nilearn HTML report object.
    """
    logger.info("Generating HTML report...")
    report = flm.generate_report(
        contrasts=contrasts,
        alpha=alpha,
        height_control=height_control,
        two_sided=True,
        plot_type='glass',
    )
    logger.info("Report generation complete.")
    return report


def get_default_contrasts(flm):
    """
    Generate identity contrasts for each condition in the design matrix.

    Excludes drift/constant columns from contrast generation.

    Parameters
    ----------
    flm : FirstLevelModel
        Fitted model.

    Returns
    -------
    contrasts : list of str
        One contrast string per condition column.
    """
    dm = flm.design_matrices_[0]
    columns = list(dm.columns)

    # Exclude drift and constant columns
    skip_prefixes = ('drift_', 'cosine_', 'constant')
    conditions = [
        col for col in columns
        if not any(col.lower().startswith(p) for p in skip_prefixes)
    ]

    logger.info(f"Auto-generated identity contrasts: {conditions}")

    return conditions
