"""Main pipeline orchestration for Trialix."""

import logging
import time

from trialix import __version__
from trialix.bids_utils import (
    find_bold_file,
    find_confounds_file,
    find_events_file,
    find_events_file_in_dir,
    get_repetition_time,
    get_sessions,
    get_subjects,
    load_confounds,
    load_events,
    setup_bids_layout,
    validate_task,
)
from trialix.cli import Colors, parse_derivatives_arg
from trialix.config import load_config, merge_config, validate_config
from trialix.logging_setup import setup_logging
from trialix.model import (
    build_first_level_model,
    compute_contrasts,
    fit_model,
    generate_report,
    get_default_contrasts,
    threshold_contrast_map,
)
from trialix.outputs import (
    save_contrast_maps,
    save_report,
    save_thresholded_map,
    write_dataset_description,
)

logger = logging.getLogger("trialix")


def _print_banner():
    """Print startup banner."""
    print(f"""
{Colors.BOLD}{Colors.GREEN}╔══════════════════════════════════════════════════════════════════════════════╗
║                          TRIALIX v{__version__}                                       ║
║           First-Level fMRI Analysis for Trial-Based Designs                   ║
╚══════════════════════════════════════════════════════════════════════════════╝{Colors.END}
""")


def _resolve_threshold_pairs(alpha_list, height_control_list):
    """
    Pair alpha values with height_control methods.

    If lists differ in length, pad the shorter one:
    - missing height_control -> None (uncorrected)
    - missing alpha -> reuse last alpha

    Returns
    -------
    pairs : list of (float, str or None)
    """
    if height_control_list is None:
        height_control_list = []

    # Normalize 'none' to None
    hc_normalized = []
    for hc in height_control_list:
        if hc is None or (isinstance(hc, str) and hc.lower() == 'none'):
            hc_normalized.append(None)
        else:
            hc_normalized.append(hc)

    n_alpha = len(alpha_list)
    n_hc = len(hc_normalized)

    pairs = []
    n = max(n_alpha, n_hc)
    for i in range(n):
        a = alpha_list[min(i, n_alpha - 1)]
        hc = hc_normalized[i] if i < n_hc else None
        pairs.append((a, hc))

    return pairs


def run_pipeline(args):
    """
    Main pipeline entry point.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    """
    _print_banner()
    pipeline_start = time.time()

    # ------------------------------------------------------------------
    # 1. Load config and merge
    # ------------------------------------------------------------------
    config_dict = {}
    if args.config:
        logger.info(f"Loading config file: {args.config}")
        config_dict = load_config(args.config)

    config = merge_config(args, config_dict)

    # Setup logging (needs output_dir)
    log = setup_logging(config['output_dir'], debug=config['debug'])
    log.info("=" * 72)
    log.info(f"Trialix v{__version__} starting")
    log.info("=" * 72)

    # Validate merged config
    config = validate_config(config)
    log.info("Configuration validated successfully.")

    # Log configuration summary
    log.info(f"  Input dir:      {config['input_dir']}")
    log.info(f"  Output dir:     {config['output_dir']}")
    log.info(f"  Task:           {config.get('task')}")
    log.info(f"  Session:        {config.get('session', 'all')}")
    log.info(f"  Space:          {config.get('space', 'MNI152NLin2009cAsym')}")
    log.info(f"  Conditions:     {config.get('conditions', 'all')}")
    log.info(f"  Regressors:     {config.get('regressors', 'none')}")
    log.info(f"  Contrasts:      {config.get('contrasts', 'auto')}")
    log.info(f"  Alpha:          {config.get('alpha')}")
    log.info(f"  Height control: {config.get('height_control', 'none')}")
    log.info(f"  HRF model:      {config.get('hrf_model', 'glover')}")
    log.info(f"  Drift model:    {config.get('drift_model', 'cosine')}")
    log.info(f"  High pass:      {config.get('high_pass', 0.01)}")
    log.info(f"  Smoothing FWHM: {config.get('smoothing_fwhm', 'none')}")

    # ------------------------------------------------------------------
    # 2. Determine input mode and setup BIDS layouts
    # ------------------------------------------------------------------
    derivatives_dict = parse_derivatives_arg(config.get('derivatives'))
    has_preproc = 'preproc' in derivatives_dict
    is_raw = has_preproc  # Mode 1 if preproc derivatives specified

    raw_layout, preproc_layout = setup_bids_layout(
        config['input_dir'],
        derivatives_dict=derivatives_dict,
        is_raw=is_raw,
    )

    # Determine which layout to use for subject/task discovery
    discovery_layout = raw_layout if is_raw else preproc_layout

    # Validate task
    task = config['task']
    validate_task(discovery_layout, task)

    # Get space
    space = config.get('space', 'MNI152NLin2009cAsym')

    # ------------------------------------------------------------------
    # 3. Get subjects and sessions to process
    # ------------------------------------------------------------------
    subjects = get_subjects(
        discovery_layout,
        participant_label=config.get('participant_label'),
    )
    log.info(f"Subjects to process: {subjects}")

    sessions = get_sessions(
        discovery_layout,
        session=config.get('session'),
    )
    log.info(f"Sessions to process: {sessions}")

    # ------------------------------------------------------------------
    # 4. Write dataset_description.json
    # ------------------------------------------------------------------
    write_dataset_description(config['output_dir'])

    # ------------------------------------------------------------------
    # 5. Prepare thresholding parameters
    # ------------------------------------------------------------------
    alpha_list = config.get('alpha', [0.05])
    hc_list = config.get('height_control')
    threshold_pairs = _resolve_threshold_pairs(alpha_list, hc_list)
    log.info(f"Thresholding: {threshold_pairs}")

    # ------------------------------------------------------------------
    # 6. Model parameters for nilearn
    # ------------------------------------------------------------------
    model_params = {
        'hrf_model': config.get('hrf_model', 'glover'),
        'drift_model': config.get('drift_model', 'cosine'),
        'drift_order': config.get('drift_order', 1),
        'high_pass': config.get('high_pass', 0.01),
        'smoothing_fwhm': config.get('smoothing_fwhm'),
    }

    # ------------------------------------------------------------------
    # 7. Main processing loop
    # ------------------------------------------------------------------
    n_total = len(subjects) * len(sessions)
    n_done = 0
    n_failed = 0

    for subject in subjects:
        for session in sessions:
            n_done += 1
            ses_str = f" ses-{session}" if session else ""
            log.info("")
            log.info("=" * 72)
            log.info(f"Processing sub-{subject}{ses_str} "
                     f"[{n_done}/{n_total}]")
            log.info("=" * 72)

            try:
                _process_subject_session(
                    config=config,
                    raw_layout=raw_layout,
                    preproc_layout=preproc_layout,
                    is_raw=is_raw,
                    subject=subject,
                    session=session,
                    task=task,
                    space=space,
                    model_params=model_params,
                    threshold_pairs=threshold_pairs,
                )
            except Exception as e:
                n_failed += 1
                log.error(
                    f"FAILED: sub-{subject}{ses_str}: {e}"
                )
                if config['debug']:
                    import traceback
                    log.error(traceback.format_exc())

    # ------------------------------------------------------------------
    # 8. Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - pipeline_start
    log.info("")
    log.info("=" * 72)
    log.info("PIPELINE COMPLETE")
    log.info(f"  Processed: {n_done - n_failed}/{n_total} successfully")
    if n_failed:
        log.warning(f"  Failed:    {n_failed}/{n_total}")
    log.info(f"  Time:      {elapsed:.1f}s")
    log.info(f"  Output:    {config['output_dir']}")
    log.info("=" * 72)


def _process_subject_session(config, raw_layout, preproc_layout, is_raw,
                             subject, session, task, space, model_params,
                             threshold_pairs):
    """
    Run the full first-level analysis for one subject/session.

    Parameters
    ----------
    config : dict
        Merged configuration.
    raw_layout : BIDSLayout or None
        Rawdata layout (Mode 1 only).
    preproc_layout : BIDSLayout
        Preprocessed data layout.
    is_raw : bool
        Whether input is rawdata (Mode 1).
    subject : str
        Subject label.
    session : str or None
        Session label.
    task : str
        Task label.
    space : str
        Template space.
    model_params : dict
        FirstLevelModel parameters.
    threshold_pairs : list of (float, str or None)
        (alpha, height_control) pairs for thresholding.
    """
    ses_str = f" ses-{session}" if session else ""

    # ------------------------------------------------------------------
    # A. Find BOLD file
    # ------------------------------------------------------------------
    logger.info(f"Looking for BOLD file: sub-{subject}{ses_str} "
                f"task-{task} space-{space}")
    bold_file = find_bold_file(preproc_layout, subject, task, session, space)
    logger.info(f"  BOLD: {bold_file}")

    # ------------------------------------------------------------------
    # B. Get TR
    # ------------------------------------------------------------------
    t_r = get_repetition_time(preproc_layout, bold_file)

    # Parse BIDS entities from selected BOLD file to resolve run-specific
    # events files when available.
    bold_entities = preproc_layout.parse_file_entities(bold_file)
    run = bold_entities.get('run')

    # ------------------------------------------------------------------
    # C. Find and load events
    # ------------------------------------------------------------------
    events_file = config.get('events_file')
    if events_file:
        # Mode 2: user-provided events file
        logger.info(f"  Events (user-provided): {events_file}")
    elif is_raw and raw_layout:
        # Mode 1: find events in rawdata
        events_file = find_events_file(
            raw_layout, subject, task, session=session, run=run
        )
        logger.info(f"  Events (rawdata): {events_file}")
    elif config.get('events_dir'):
        events_file = find_events_file_in_dir(
            config['events_dir'], subject, task, session=session, run=run
        )
        logger.info(f"  Events (events-dir): {events_file}")
    else:
        raise FileNotFoundError(
            f"No events file for sub-{subject}{ses_str}. "
            f"Use --events-file, --events-dir, or provide rawdata with "
            f"--derivatives."
        )

    events_df = load_events(events_file, conditions=config.get('conditions'))

    # ------------------------------------------------------------------
    # D. Find and load confounds
    # ------------------------------------------------------------------
    confounds_file = config.get('confounds_file')
    confounds_df = None

    if confounds_file:
        logger.info(f"  Confounds (user-provided): {confounds_file}")
    elif config.get('regressors'):
        # Need confounds for regressors; find from fMRIPrep
        confounds_file = find_confounds_file(
            preproc_layout, subject, task, session
        )
        if confounds_file:
            logger.info(f"  Confounds (fMRIPrep): {confounds_file}")
        else:
            raise FileNotFoundError(
                f"Regressors requested but no confounds file found for "
                f"sub-{subject}{ses_str}."
            )

    confounds_df = load_confounds(confounds_file, config.get('regressors'))
    if confounds_df is not None:
        logger.info(f"  Confound regressors: {list(confounds_df.columns)}")
    else:
        logger.info("  No confound regressors included.")

    # ------------------------------------------------------------------
    # E. Build and fit model
    # ------------------------------------------------------------------
    flm = build_first_level_model(
        t_r=t_r,
        subject_label=subject,
        model_params=model_params,
        debug=config['debug'],
    )

    flm = fit_model(flm, bold_file, events_df, confounds_df)

    # ------------------------------------------------------------------
    # F. Determine contrasts
    # ------------------------------------------------------------------
    contrasts = config.get('contrasts')
    if not contrasts:
        contrasts = get_default_contrasts(flm)
        logger.info(f"Using auto-generated contrasts: {contrasts}")
    else:
        logger.info(f"Using user-defined contrasts: {contrasts}")

    # ------------------------------------------------------------------
    # G. Compute contrasts
    # ------------------------------------------------------------------
    contrast_results = compute_contrasts(flm, contrasts)

    # ------------------------------------------------------------------
    # H. Save outputs
    # ------------------------------------------------------------------
    output_dir = config['output_dir']

    for contrast_name, outputs in contrast_results.items():
        logger.info(f"Saving maps for contrast: '{contrast_name}'")

        # Save unthresholded maps
        save_contrast_maps(
            outputs, output_dir, subject, session, task, space, contrast_name
        )

        # Apply thresholding and save thresholded maps
        z_map = outputs.get('z_score')
        if z_map is not None:
            for alpha, hc in threshold_pairs:
                try:
                    thresholded, threshold_val = threshold_contrast_map(
                        z_map, alpha=alpha, height_control=hc
                    )
                    save_thresholded_map(
                        thresholded, output_dir, subject, session, task,
                        space, contrast_name, alpha, hc
                    )
                except Exception as e:
                    logger.warning(
                        f"Thresholding failed (alpha={alpha}, "
                        f"hc={hc}): {e}"
                    )

    # ------------------------------------------------------------------
    # I. Generate and save report
    # ------------------------------------------------------------------
    try:
        alpha, hc = threshold_pairs[0] if threshold_pairs else (0.001, 'fpr')
        report = generate_report(flm, contrasts, alpha=alpha, height_control=hc)
        save_report(report, output_dir, subject, session, task, space)
    except Exception as e:
        logger.warning(f"Report generation failed: {e}")

    logger.info(f"sub-{subject}{ses_str} complete.")
