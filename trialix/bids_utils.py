"""BIDS data discovery and validation utilities for Trialix."""

import logging
import re
from pathlib import Path

import pandas as pd
from bids import BIDSLayout

logger = logging.getLogger("trialix")


def setup_bids_layout(input_dir, derivatives_dict=None, is_raw=True):
    """
    Create BIDSLayout object(s) for data discovery.

    Parameters
    ----------
    input_dir : Path
        Path to BIDS rawdata or preprocessed data.
    derivatives_dict : dict or None
        Dict of name -> path for derivatives (from --derivatives parsing).
    is_raw : bool
        If True, input_dir is rawdata. If False, it is preprocessed data.

    Returns
    -------
    raw_layout : BIDSLayout or None
        Layout for rawdata (only in Mode 1).
    preproc_layout : BIDSLayout
        Layout for preprocessed data.
    """
    input_dir = Path(input_dir).resolve()

    if is_raw:
        logger.info(f"Setting up BIDS layout from rawdata: {input_dir}")
        raw_layout = BIDSLayout(str(input_dir), validate=False)

        preproc_layout = None
        if derivatives_dict and 'preproc' in derivatives_dict:
            preproc_path = Path(derivatives_dict['preproc']).resolve()
            logger.info(f"Setting up derivatives layout: {preproc_path}")
            preproc_layout = BIDSLayout(
                str(preproc_path), validate=False, is_derivative=True
            )

        return raw_layout, preproc_layout
    else:
        logger.info(f"Setting up BIDS layout from derivatives: {input_dir}")
        preproc_layout = BIDSLayout(
            str(input_dir), validate=False, is_derivative=True
        )
        return None, preproc_layout


def get_subjects(layout, participant_label=None):
    """
    Get list of subjects to process, validating requested labels.

    Parameters
    ----------
    layout : BIDSLayout
        BIDS layout for subject discovery.
    participant_label : list of str or None
        Subjects to process. None means all available subjects.

    Returns
    -------
    subjects : list of str
        Validated subject labels (without 'sub-' prefix).
    """
    available = layout.get_subjects()
    logger.debug(f"Available subjects in layout: {available}")

    if participant_label is None:
        return sorted(available)

    missing = [s for s in participant_label if s not in available]
    if missing:
        raise ValueError(
            f"Requested subjects not found in dataset: {missing}. "
            f"Available subjects: {sorted(available)}"
        )

    return sorted(participant_label)


def validate_task(layout, task):
    """Validate that the given task exists in the BIDS layout."""
    available_tasks = layout.get_tasks()
    if task not in available_tasks:
        raise ValueError(
            f"Task '{task}' not found in dataset. "
            f"Available tasks: {available_tasks}"
        )
    logger.debug(f"Task '{task}' validated (available: {available_tasks})")


def get_sessions(layout, session=None):
    """
    Get list of sessions to process.

    Parameters
    ----------
    layout : BIDSLayout
        BIDS layout for session discovery.
    session : str or None
        Specific session to process. None means all sessions.

    Returns
    -------
    sessions : list of str or [None]
        Session labels, or [None] if dataset has no sessions.
    """
    available = layout.get_sessions()

    if not available:
        return [None]

    if session is not None:
        if session not in available:
            raise ValueError(
                f"Session '{session}' not found. "
                f"Available sessions: {available}"
            )
        return [session]

    return sorted(available)


def find_bold_file(preproc_layout, subject, task, session=None,
                   space="MNI152NLin2009cAsym"):
    """
    Find the preprocessed BOLD NIfTI file.

    Parameters
    ----------
    preproc_layout : BIDSLayout
        Layout for preprocessed data (fMRIPrep output).
    subject : str
        Subject label.
    task : str
        Task label.
    session : str or None
        Session label.
    space : str
        Template space.

    Returns
    -------
    bold_path : str
        Absolute path to the preprocessed BOLD file.
    """
    filters = dict(
        subject=subject,
        task=task,
        space=space,
        desc='preproc',
        suffix='bold',
        extension=['.nii.gz', '.nii'],
        return_type='file',
    )
    if session:
        filters['session'] = session

    files = preproc_layout.get(**filters)

    if not files:
        ses_str = f" ses-{session}" if session else ""
        raise FileNotFoundError(
            f"No preprocessed BOLD file found for sub-{subject}{ses_str} "
            f"task-{task} space-{space}. Check that the space label matches "
            f"the preprocessed data."
        )

    if len(files) > 1:
        logger.warning(
            f"Multiple BOLD files found for sub-{subject}, using: {files[0]}"
        )

    logger.debug(f"BOLD file: {files[0]}")
    return files[0]


def _normalize_entity(value):
    """Normalize entity values for robust string matching."""
    if value is None:
        return None
    return str(value).strip()


def _extract_bids_entities_from_name(file_path):
    """Extract BIDS-like entities from a filename (best effort)."""
    name = Path(file_path).name
    entities = {}

    patterns = {
        'subject': r'(?:^|_)sub-([^_]+)',
        'task': r'(?:^|_)task-([^_]+)',
        'session': r'(?:^|_)ses-([^_]+)',
        'run': r'(?:^|_)run-([^_]+)',
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, name)
        if match:
            entities[key] = match.group(1)

    return entities


def _filter_events_candidates(files, subject, task, session=None, run=None,
                              allow_global_subject=True):
    """Filter candidate event files by required BIDS entities."""
    subject = _normalize_entity(subject)
    task = _normalize_entity(task)
    session = _normalize_entity(session)
    run = _normalize_entity(run)

    filtered = []
    for fpath in files:
        entities = _extract_bids_entities_from_name(fpath)

        if _normalize_entity(entities.get('task')) != task:
            continue

        # Subject can be specific (sub-XX...) or optionally global.
        sub_entity = _normalize_entity(entities.get('subject'))
        if sub_entity is None:
            if not allow_global_subject:
                continue
        elif sub_entity != subject:
            continue

        # If session/run are present in filename, they must match.
        ses_entity = _normalize_entity(entities.get('session'))
        if ses_entity is not None and ses_entity != session:
            continue

        run_entity = _normalize_entity(entities.get('run'))
        if run_entity is not None and run_entity != run:
            continue

        filtered.append(str(fpath))

    return sorted(filtered)


def find_events_file(raw_layout, subject, task, session=None, run=None):
    """
    Find the events TSV file from rawdata.

    Parameters
    ----------
    raw_layout : BIDSLayout
        Layout for rawdata.
    subject : str
        Subject label.
    task : str
        Task label.
    session : str or None
        Session label.
    run : str or int or None
        Run label.

    Returns
    -------
    events_path : str
        Path to the events TSV file.
    """
    files = []

    # Preferred query path when available in pybids.
    if hasattr(raw_layout, 'get_events'):
        try:
            files = raw_layout.get_events(
                subject=subject,
                task=task,
                session=session,
                run=run,
                return_type='file',
            )
        except TypeError:
            # Some pybids versions expose get_events with a narrower signature.
            files = raw_layout.get_events(return_type='file')

    # Fallback query path for compatibility across pybids versions.
    if not files:
        filters = dict(
            task=task,
            suffix='events',
            extension='.tsv',
            return_type='file',
        )
        if session:
            filters['session'] = session
        files = raw_layout.get(**filters)

    matched = _filter_events_candidates(
        files,
        subject=subject,
        task=task,
        session=session,
        run=run,
        allow_global_subject=False,
    )

    # If no participant-specific file is found, allow global task events.
    if not matched:
        logger.debug(
            f"No participant-specific events for sub-{subject}, "
            f"trying global task-level events."
        )
        global_candidates = []
        for fpath in files:
            entities = _extract_bids_entities_from_name(fpath)
            if 'subject' in entities:
                continue
            if _normalize_entity(entities.get('task')) != _normalize_entity(task):
                continue
            ses_entity = _normalize_entity(entities.get('session'))
            if ses_entity is not None and ses_entity != _normalize_entity(session):
                continue
            run_entity = _normalize_entity(entities.get('run'))
            if run_entity is not None and run_entity != _normalize_entity(run):
                continue
            global_candidates.append(str(fpath))
        matched = sorted(global_candidates)

    if not matched:
        ses_str = f" ses-{session}" if session else ""
        run_str = f" run-{run}" if run is not None else ""
        raise FileNotFoundError(
            f"No events file found for sub-{subject}{ses_str}{run_str} "
            f"task-{task} in rawdata."
        )

    if len(matched) > 1:
        ses_str = f" ses-{session}" if session else ""
        run_str = f" run-{run}" if run is not None else ""
        raise ValueError(
            f"Multiple events files matched for sub-{subject}{ses_str}{run_str} "
            f"task-{task}: {matched}. Refine your dataset to keep a unique match."
        )

    logger.debug(f"Events file: {matched[0]}")
    return matched[0]


def find_events_file_in_dir(events_dir, subject, task, session=None, run=None):
    """
    Recursively search an events directory and match by BIDS entities.

    This is used when INPUT_DIR is preprocessed data and --events-file is not
    provided. The directory does not need to be a full BIDS dataset.
    """
    events_dir = Path(events_dir)
    if not events_dir.is_dir():
        raise ValueError(f"Events directory does not exist: {events_dir}")

    candidates = [
        str(p) for p in events_dir.rglob("*_events.tsv") if p.is_file()
    ]

    matched = _filter_events_candidates(
        candidates,
        subject=subject,
        task=task,
        session=session,
        run=run,
        allow_global_subject=False,
    )

    ses_str = f" ses-{session}" if session else ""
    run_str = f" run-{run}" if run is not None else ""

    if not matched:
        raise FileNotFoundError(
            f"No events file found in {events_dir} for "
            f"sub-{subject}{ses_str}{run_str} task-{task}."
        )

    if len(matched) > 1:
        raise ValueError(
            f"Multiple events files matched in {events_dir} for "
            f"sub-{subject}{ses_str}{run_str} task-{task}: {matched}."
        )

    logger.debug(f"Events file from events-dir: {matched[0]}")
    return matched[0]


def find_confounds_file(preproc_layout, subject, task, session=None):
    """
    Find the confounds timeseries TSV file from fMRIPrep output.

    Tries both 'timeseries' and 'regressors' suffixes for compatibility
    across fMRIPrep versions.

    Parameters
    ----------
    preproc_layout : BIDSLayout
        Layout for preprocessed data.
    subject : str
        Subject label.
    task : str
        Task label.
    session : str or None
        Session label.

    Returns
    -------
    confounds_path : str or None
        Path to confounds TSV file, or None if not found.
    """
    base_filters = dict(
        subject=subject,
        task=task,
        desc='confounds',
        extension='.tsv',
        return_type='file',
    )
    if session:
        base_filters['session'] = session

    # Try 'timeseries' suffix (newer fMRIPrep)
    filters = {**base_filters, 'suffix': 'timeseries'}
    files = preproc_layout.get(**filters)
    if files:
        logger.debug(f"Confounds file: {files[0]}")
        return files[0]

    # Try 'regressors' suffix (older fMRIPrep)
    filters = {**base_filters, 'suffix': 'regressors'}
    files = preproc_layout.get(**filters)
    if files:
        logger.debug(f"Confounds file (legacy): {files[0]}")
        return files[0]

    ses_str = f" ses-{session}" if session else ""
    logger.warning(
        f"No confounds file found for sub-{subject}{ses_str} task-{task}"
    )
    return None


def get_repetition_time(preproc_layout, bold_file):
    """
    Extract repetition time (TR) from BOLD file metadata or NIfTI header.

    Parameters
    ----------
    preproc_layout : BIDSLayout
        Layout for preprocessed data.
    bold_file : str
        Path to the BOLD NIfTI file.

    Returns
    -------
    t_r : float
        Repetition time in seconds.
    """
    # Try BIDS JSON sidecar metadata first
    try:
        metadata = preproc_layout.get_metadata(bold_file)
        if 'RepetitionTime' in metadata:
            t_r = float(metadata['RepetitionTime'])
            logger.info(f"TR from BIDS metadata: {t_r}s")
            return t_r
    except Exception:
        logger.debug("Could not read TR from BIDS metadata, trying NIfTI header")

    # Fallback to NIfTI header
    import nibabel as nib
    img = nib.load(bold_file)
    zooms = img.header.get_zooms()
    if len(zooms) >= 4 and float(zooms[3]) > 0:
        t_r = float(zooms[3])
        logger.info(f"TR from NIfTI header: {t_r}s")
        return t_r

    raise ValueError(
        f"Could not determine TR for {bold_file}. "
        f"Ensure BIDS metadata (JSON sidecar) or NIfTI header contains TR."
    )


def load_events(events_file, conditions=None):
    """
    Load events from TSV and optionally filter by trial_type.

    Parameters
    ----------
    events_file : str or Path
        Path to events TSV file.
    conditions : list of str or None
        Trial types to include. None means include all.

    Returns
    -------
    events : pd.DataFrame
        Events DataFrame with onset, duration, trial_type columns.
    """
    events = pd.read_csv(events_file, sep='\t')

    required_cols = ['onset', 'duration', 'trial_type']
    missing = [c for c in required_cols if c not in events.columns]
    if missing:
        raise ValueError(
            f"Events file {events_file} missing required columns: {missing}. "
            f"Found columns: {list(events.columns)}"
        )

    available_types = sorted(events['trial_type'].unique().tolist())
    logger.info(f"Trial types in events file: {available_types}")

    if conditions:
        invalid = [c for c in conditions if c not in available_types]
        if invalid:
            raise ValueError(
                f"Requested conditions not in events file: {invalid}. "
                f"Available trial types: {available_types}"
            )
        events = events[events['trial_type'].isin(conditions)].copy()
        logger.info(f"Filtered events to conditions: {conditions}")

    return events


def load_confounds(confounds_file, regressors=None):
    """
    Load confound regressors from TSV file.

    Parameters
    ----------
    confounds_file : str or Path
        Path to confounds TSV file.
    regressors : list of str or None
        Column names to select. None means no confounds used.

    Returns
    -------
    confounds : pd.DataFrame or None
        Selected confound columns, or None if no regressors.
    """
    if not regressors or confounds_file is None:
        return None

    confounds = pd.read_csv(confounds_file, sep='\t')

    missing = [r for r in regressors if r not in confounds.columns]
    if missing:
        available = list(confounds.columns)
        # Show first 30 columns to avoid overwhelming output
        shown = available[:30]
        suffix = f" ... ({len(available)} total)" if len(available) > 30 else ""
        raise ValueError(
            f"Requested regressors not in confounds file: {missing}. "
            f"Available columns: {shown}{suffix}"
        )

    selected = confounds[regressors].copy()

    # Handle NaN values (common in first row for derivative-based regressors)
    n_nan = int(selected.isna().sum().sum())
    if n_nan > 0:
        logger.warning(
            f"Filling {n_nan} NaN value(s) in confounds with 0 "
            f"(common for derivative columns in first volume)"
        )
        selected = selected.fillna(0)

    return selected
