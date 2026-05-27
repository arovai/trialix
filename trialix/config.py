"""Configuration loading and merging for Trialix."""

from pathlib import Path

import yaml


# Keys that can appear in both CLI and config (CLI overrides config)
SHARED_KEYS = {
    'participant_label',
    'task',
    'session',
    'space',
    'conditions',
    'regressors',
    'contrasts',
    'alpha',
    'height_control',
}

# Config-only keys passed to nilearn FirstLevelModel
MODEL_KEYS = {
    'hrf_model',
    'drift_model',
    'drift_order',
    'high_pass',
    'smoothing_fwhm',
}

# Defaults for model parameters
MODEL_DEFAULTS = {
    'hrf_model': 'glover',
    'drift_model': 'cosine',
    'drift_order': 1,
    'high_pass': 0.01,
    'smoothing_fwhm': None,
}


def load_config(config_path):
    """
    Load a YAML configuration file.

    Parameters
    ----------
    config_path : str or Path
        Path to the YAML configuration file.

    Returns
    -------
    config : dict
        Parsed configuration dictionary.
    """
    config_path = Path(config_path)
    if not config_path.is_file():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    if config is None:
        config = {}

    return config


def merge_config(args, config=None):
    """
    Merge CLI arguments with config file values.

    CLI values override config values for shared keys.
    Config-only model keys are taken from config with defaults as fallback.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    config : dict or None
        Configuration from YAML file.

    Returns
    -------
    merged : dict
        Unified configuration dictionary.
    """
    if config is None:
        config = {}

    merged = {}

    # Model parameters: config -> defaults
    for key in MODEL_KEYS:
        merged[key] = config.get(key, MODEL_DEFAULTS.get(key))

    # Shared keys: CLI overrides config
    for key in SHARED_KEYS:
        cli_val = getattr(args, key, None)
        config_val = config.get(key)

        if cli_val is not None:
            merged[key] = cli_val
        elif config_val is not None:
            merged[key] = config_val

    # CLI-only values
    merged['input_dir'] = args.input_dir
    merged['output_dir'] = args.output_dir
    merged['analysis_level'] = args.analysis_level
    merged['debug'] = args.debug
    merged['events_file'] = getattr(args, 'events_file', None)
    merged['events_dir'] = getattr(args, 'events_dir', None)
    merged['confounds_file'] = getattr(args, 'confounds_file', None)
    merged['derivatives'] = args.derivatives  # raw list from argparse

    # Ensure alpha has a default
    if merged.get('alpha') is None:
        merged['alpha'] = [0.05]

    return merged


def validate_config(config):
    """
    Validate the merged configuration for consistency.

    Parameters
    ----------
    config : dict
        Merged configuration dictionary.

    Raises
    ------
    ValueError
        If required values are missing or invalid.
    """
    from trialix.cli import parse_derivatives_arg

    # Task is required
    if not config.get('task'):
        raise ValueError(
            "--task is required. Specify via CLI or config file."
        )

    # Parse derivatives to determine input mode
    derivatives_dict = parse_derivatives_arg(config.get('derivatives'))
    has_preproc = 'preproc' in derivatives_dict
    has_events_file = config.get('events_file') is not None
    has_events_dir = config.get('events_dir') is not None

    # Input mode rules:
    # - Rawdata mode: requires --derivatives preproc=PATH.
    # - Preprocessed mode: requires --events-file OR --events-dir.
    if not has_preproc and not has_events_file and not has_events_dir:
        raise ValueError(
            "Either --derivatives preproc=PATH must be provided (when INPUT_DIR "
            "is rawdata), or --events-file/--events-dir must be provided "
            "(when INPUT_DIR is preprocessed data)."
        )

    # Validate input_dir exists
    input_dir = config['input_dir']
    if not input_dir.is_dir():
        raise ValueError(f"INPUT_DIR does not exist: {input_dir}")

    # Validate derivatives paths exist
    for key, path in derivatives_dict.items():
        if not path.is_dir():
            raise ValueError(f"Derivatives path for '{key}' does not exist: {path}")

    # Validate events-file exists if provided
    events_file = config.get('events_file')
    if events_file and not events_file.is_file():
        raise ValueError(f"Events file does not exist: {events_file}")

    # Validate events-dir exists if provided
    events_dir = config.get('events_dir')
    if events_dir and not events_dir.is_dir():
        raise ValueError(f"Events directory does not exist: {events_dir}")

    # Validate confounds-file exists if provided
    confounds_file = config.get('confounds_file')
    if confounds_file and not confounds_file.is_file():
        raise ValueError(f"Confounds file does not exist: {confounds_file}")

    # Validate alpha values
    alpha_list = config.get('alpha', [0.05])
    for a in alpha_list:
        if not 0 < a < 1:
            raise ValueError(f"Alpha must be between 0 and 1, got: {a}")

    # Validate height_control values
    hc_list = config.get('height_control') or []
    valid_hc = {'fdr', 'bonferroni', 'none', None}
    for hc in hc_list:
        if hc not in valid_hc:
            raise ValueError(
                f"Invalid height_control: '{hc}'. "
                f"Must be 'fdr', 'bonferroni', or 'none'."
            )

    # Validate HRF model
    valid_hrf = {
        'glover', 'glover + derivative', 'glover + derivative + dispersion',
        'spm', 'spm + derivative', 'spm + derivative + dispersion',
        'fir',
    }
    hrf = config.get('hrf_model', 'glover')
    if isinstance(hrf, str) and hrf not in valid_hrf:
        raise ValueError(
            f"Invalid hrf_model: '{hrf}'. Valid options: {valid_hrf}"
        )

    # Validate drift model
    valid_drift = {'cosine', 'polynomial', None}
    drift = config.get('drift_model', 'cosine')
    if drift not in valid_drift:
        raise ValueError(
            f"Invalid drift_model: '{drift}'. Valid options: {valid_drift}"
        )

    return config
