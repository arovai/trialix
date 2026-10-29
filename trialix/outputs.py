"""Output path construction and file saving for Trialix."""

import json
import logging
import re
from pathlib import Path

import nibabel as nib

from trialix import __version__

logger = logging.getLogger("trialix")

# Map nilearn output_type keys to BIDS desc entity values
OUTPUT_TYPE_MAP = {
    'z_score': 'zScore',
    'stat': 'tStat',
    'p_value': 'pValue',
    'effect_size': 'effectSize',
    'effect_variance': 'effectVariance',
}


def sanitize_contrast_name(contrast_str):
    """
    Convert a contrast string to a valid BIDS entity value.

    Examples
    --------
    >>> sanitize_contrast_name("cond1 - cond2")
    'cond1MinusCond2'
    >>> sanitize_contrast_name("left_hand + right_hand")
    'leftHandPlusRightHand'
    """
    name = contrast_str.strip()
    name = name.replace(' - ', 'Minus')
    name = name.replace(' + ', 'Plus')
    name = name.replace(' * ', 'Times')
    name = name.replace(' / ', 'Div')
    name = name.replace('-', 'Minus')
    name = name.replace('+', 'Plus')
    name = name.replace('*', 'Times')
    name = name.replace('/', 'Div')
    name = name.replace(' ', '')
    # Remove any remaining non-alphanumeric characters (except underscore)
    name = re.sub(r'[^a-zA-Z0-9]', '', name)
    return name


def build_bids_prefix(subject, session=None, task=None, space=None):
    """
    Build the BIDS filename prefix.

    Returns
    -------
    prefix : str
        e.g., "sub-01_ses-pre_task-motor_space-MNI152NLin2009cAsym"
    """
    parts = [f"sub-{subject}"]
    if session:
        parts.append(f"ses-{session}")
    if task:
        parts.append(f"task-{task}")
    if space:
        parts.append(f"space-{space}")
    return "_".join(parts)


def build_output_dir(output_dir, subject, session=None):
    """
    Build and create the BIDS-compliant output directory for a subject.

    Returns
    -------
    out_dir : Path
        e.g., OUTPUT_DIR/sub-01/ses-pre/ or OUTPUT_DIR/sub-01/
    """
    out_dir = Path(output_dir) / f"sub-{subject}"
    if session:
        out_dir = out_dir / f"ses-{session}"
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def save_contrast_maps(contrast_outputs, output_dir, subject, session, task,
                       space, contrast_name, metadata=None):
    """
    Save all output maps for a single contrast.

    Parameters
    ----------
    contrast_outputs : dict
        {output_type: Nifti1Image} from compute_contrast(output_type='all')
    output_dir : Path
        Base output directory.
    subject : str
        Subject label.
    session : str or None
        Session label.
    task : str
        Task label.
    space : str
        Space label.
    contrast_name : str
        Original contrast string.

    Returns
    -------
    saved_files : list of Path
        Paths to saved files.
    """
    out_dir = build_output_dir(output_dir, subject, session)
    prefix = build_bids_prefix(subject, session, task, space)
    con_label = sanitize_contrast_name(contrast_name)

    saved = []

    for output_type, img in contrast_outputs.items():
        desc = OUTPUT_TYPE_MAP.get(output_type, output_type)
        fname = f"{prefix}_contrast-{con_label}_desc-{desc}_statmap.nii.gz"
        fpath = out_dir / fname

        nib.save(img, str(fpath))
        _write_map_metadata(fpath, metadata, Contrast=contrast_name)
        logger.info(f"  Saved: {fpath.name}")
        saved.append(fpath)

    return saved


def save_thresholded_map(thresholded_img, output_dir, subject, session, task,
                         space, contrast_name, alpha, height_control, metadata=None):
    """
    Save a thresholded z-score map with threshold info in the filename.

    Parameters
    ----------
    thresholded_img : Nifti1Image
        Thresholded z-score map.
    output_dir : Path
        Base output directory.
    subject, session, task, space : str
        BIDS entities.
    contrast_name : str
        Original contrast string.
    alpha : float
        Alpha threshold used.
    height_control : str or None
        Correction method used.

    Returns
    -------
    fpath : Path
        Path to saved file.
    """
    out_dir = build_output_dir(output_dir, subject, session)
    prefix = build_bids_prefix(subject, session, task, space)
    con_label = sanitize_contrast_name(contrast_name)

    # Build descriptive label for thresholding parameters
    alpha_str = str(alpha).replace('.', 'p')
    hc_str = height_control if height_control else "unc"
    desc = f"zScoreThresh{hc_str.capitalize()}Alpha{alpha_str}"

    fname = f"{prefix}_contrast-{con_label}_desc-{desc}_statmap.nii.gz"
    fpath = out_dir / fname

    nib.save(thresholded_img, str(fpath))
    _write_map_metadata(fpath, metadata, Contrast=contrast_name,
                        Alpha=alpha, HeightControl=height_control)
    logger.info(f"  Saved thresholded: {fpath.name}")
    return fpath


def _write_map_metadata(fpath, metadata, **details):
    """Write provenance alongside each unthresholded or thresholded map."""
    if metadata is None:
        return
    sidecar = fpath.with_name(fpath.name.removesuffix('.nii.gz') + '.json')
    with sidecar.open('w') as stream:
        json.dump({**metadata, **details}, stream, indent=2, allow_nan=False)


def save_report(report, output_dir, subject, session, task, space):
    """
    Save the HTML report.

    Parameters
    ----------
    report : HTMLReport
        Nilearn report object.
    output_dir : Path
        Base output directory.
    subject, session, task, space : str
        BIDS entities.

    Returns
    -------
    fpath : Path
        Path to saved HTML file.
    """
    out_dir = build_output_dir(output_dir, subject, session)
    prefix = build_bids_prefix(subject, session, task, space)
    fname = f"{prefix}_report.html"
    fpath = out_dir / fname

    report.save_as_html(str(fpath))
    logger.info(f"  Saved report: {fpath.name}")
    return fpath


def write_dataset_description(output_dir):
    """
    Write a BIDS dataset_description.json to the output directory.

    Parameters
    ----------
    output_dir : Path
        Output directory root.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    desc = {
        "Name": "Trialix first-level GLM results",
        "BIDSVersion": "1.8.0",
        "PipelineDescription": {
            "Name": "trialix",
            "Version": __version__,
            "Description": (
                "First-level fMRI analysis for trial-based designs "
                "using nilearn."
            ),
        },
        "GeneratedBy": [
            {
                "Name": "trialix",
                "Version": __version__,
            }
        ],
    }

    fpath = output_dir / "dataset_description.json"
    with open(fpath, 'w') as f:
        json.dump(desc, f, indent=2)

    logger.info(f"Wrote {fpath}")
