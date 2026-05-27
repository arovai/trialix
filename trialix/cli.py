"""
Trialix CLI: Command-line interface for first-level fMRI analysis.

Follows the BIDS-application convention:
    trialix INPUT_DIR OUTPUT_DIR participant [options]
"""

import argparse
import sys
import textwrap
from pathlib import Path

from trialix import __version__


class Colors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


class ColoredHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter with colored section headers."""

    def __init__(self, prog, indent_increment=2, max_help_position=40, width=100):
        super().__init__(prog, indent_increment, max_help_position, width)

    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = f'{Colors.BOLD}Usage:{Colors.END} '
        return super()._format_usage(usage, actions, groups, prefix)

    def start_section(self, heading):
        if heading:
            heading = f'{Colors.BOLD}{Colors.CYAN}{heading}{Colors.END}'
        super().start_section(heading)


def create_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser for trialix."""

    description = textwrap.dedent(f"""
    {Colors.BOLD}{Colors.GREEN}╔══════════════════════════════════════════════════════════════════════════════╗
    ║                          TRIALIX v{__version__}                                      ║
    ║           First-Level fMRI Analysis for Trial-Based Designs                  ║
    ╚══════════════════════════════════════════════════════════════════════════════╝{Colors.END}

    {Colors.BOLD}Description:{Colors.END}
      Trialix performs BIDS-compliant first-level (participant-level) fMRI analysis
      for trial-based experimental designs (block or event-related). It wraps nilearn
      tools to fit GLM models on preprocessed fMRI data (e.g., fMRIPrep output).

    {Colors.BOLD}Workflow:{Colors.END}
      1. Discover input data from BIDS dataset structure
      2. Validate events, confounds, and preprocessed data
      3. Build first-level GLM design matrix
      4. Fit the General Linear Model (FirstLevelModel)
      5. Compute contrasts (z-score, t-stat, p-value, effect size, variance)
      6. Apply statistical thresholding
      7. Generate BIDS-compliant outputs and HTML reports
    """)

    epilog = textwrap.dedent(f"""
    {Colors.BOLD}{Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}
    {Colors.BOLD}EXAMPLES{Colors.END}
    {Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}

    {Colors.BOLD}Basic Usage (rawdata + fMRIPrep derivatives):{Colors.END}

      {Colors.YELLOW}# Process all subjects for a given task{Colors.END}
      trialix /data/rawdata /data/output participant \\
          --task motor --derivatives preproc=/data/derivatives/fmriprep

      {Colors.YELLOW}# Process specific subjects with contrasts{Colors.END}
      trialix /data/rawdata /data/output participant \\
          --task motor --participant-label 01 02 \\
          --derivatives preproc=/data/derivatives/fmriprep \\
          --conditions left_hand right_hand \\
          --contrasts "left_hand - right_hand" "left_hand + right_hand"

    {Colors.BOLD}Using Preprocessed Data Directly:{Colors.END}

      {Colors.YELLOW}# INPUT_DIR is the preprocessed data folder{Colors.END}
      trialix /data/derivatives/fmriprep /data/output participant \\
          --task motor --events-dir /data/rawdata

      {Colors.YELLOW}# Use one global events file for all participants{Colors.END}
      trialix /data/derivatives/fmriprep /data/output participant \
          --task motor --events-file /data/rawdata/task-motor_events.tsv

    {Colors.BOLD}Full Analysis:{Colors.END}

      {Colors.YELLOW}# All options: conditions, regressors, contrasts, thresholding{Colors.END}
      trialix /data/rawdata /data/output participant \\
          --task motor \\
          --derivatives preproc=/data/derivatives/fmriprep \\
          --conditions left_hand right_hand \\
          --regressors trans_x trans_y trans_z rot_x rot_y rot_z \\
          --contrasts "left_hand - right_hand" "left_hand + right_hand" \\
          --alpha 0.05 0.001 \\
          --height-control fdr bonferroni

    {Colors.BOLD}Using a Config File:{Colors.END}

      {Colors.YELLOW}# Use a YAML configuration file{Colors.END}
      trialix /data/rawdata /data/output participant --config config.yaml

      {Colors.YELLOW}# Config file + CLI overrides{Colors.END}
      trialix /data/rawdata /data/output participant \\
          --config config.yaml --participant-label 01

    {Colors.BOLD}{Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}
    {Colors.BOLD}INPUT MODES{Colors.END}
    {Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}

      {Colors.CYAN}Mode 1:{Colors.END}  INPUT_DIR = BIDS rawdata folder
               Requires --derivatives preproc=PATH
               Events files are read from the rawdata folder

      {Colors.CYAN}Mode 2:{Colors.END}  INPUT_DIR = Preprocessed data folder (e.g., fMRIPrep output)
               Requires --events-file PATH or --events-dir PATH
               Preprocessed data is read directly from INPUT_DIR

    {Colors.BOLD}{Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}
    {Colors.BOLD}CONFIG FILE{Colors.END}
    {Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}

      The following parameters are config-file only (not available via CLI)
      and are passed directly to nilearn's FirstLevelModel:

        {Colors.CYAN}hrf_model{Colors.END}       HRF model (default: 'glover')
        {Colors.CYAN}drift_model{Colors.END}     Drift model (default: 'cosine')
        {Colors.CYAN}drift_order{Colors.END}     Drift order for polynomial model (default: 1)
        {Colors.CYAN}high_pass{Colors.END}       High-pass filter cutoff in Hz (default: 0.01)
        {Colors.CYAN}smoothing_fwhm{Colors.END}  Spatial smoothing FWHM in mm (default: None)

    {Colors.BOLD}{Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}
    {Colors.BOLD}MORE INFORMATION{Colors.END}
    {Colors.GREEN}═══════════════════════════════════════════════════════════════════════════════{Colors.END}

      Documentation:  https://github.com/arovai/trialix
      Version:        {__version__}
    """)

    parser = argparse.ArgumentParser(
        prog="trialix",
        description=description,
        epilog=epilog,
        formatter_class=ColoredHelpFormatter,
        add_help=False,
    )

    # =========================================================================
    # REQUIRED ARGUMENTS
    # =========================================================================
    required = parser.add_argument_group(
        f'{Colors.BOLD}Required Arguments{Colors.END}'
    )

    required.add_argument(
        "input_dir",
        type=Path,
        metavar="INPUT_DIR",
        help="Path to BIDS rawdata folder (with --derivatives) or to "
               "preprocessed data folder (with --events-file/--events-dir).",
    )

    required.add_argument(
        "output_dir",
        type=Path,
        metavar="OUTPUT_DIR",
        help="Path to output directory for analysis derivatives.",
    )

    required.add_argument(
        "analysis_level",
        choices=["participant"],
        metavar="{participant}",
        help="Analysis level. Currently only 'participant' is supported.",
    )

    # =========================================================================
    # GENERAL OPTIONS
    # =========================================================================
    general = parser.add_argument_group(
        f'{Colors.BOLD}General Options{Colors.END}'
    )

    general.add_argument(
        "-h", "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )

    general.add_argument(
        "--version",
        action="version",
        version=f"trialix {__version__}",
        help="Show program version and exit.",
    )

    general.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode with increased verbosity.",
    )

    general.add_argument(
        "--config",
        type=Path,
        metavar="FILE",
        help="Path to YAML configuration file. "
             "CLI arguments override config file settings.",
    )

    # =========================================================================
    # DERIVATIVES OPTIONS
    # =========================================================================
    derivatives = parser.add_argument_group(
        f'{Colors.BOLD}Input Derivatives{Colors.END}'
    )

    derivatives.add_argument(
        "--derivatives",
        action="append",
        metavar="KEY=PATH",
        dest="derivatives",
        help="Specify location of BIDS derivatives. Format: key=path "
             "(e.g., preproc=/data/derivatives/fmriprep). Can be specified "
             "multiple times. Use 'preproc=' for the preprocessed data path.",
    )

    # =========================================================================
    # BIDS ENTITY FILTERS
    # =========================================================================
    filters = parser.add_argument_group(
        f'{Colors.BOLD}BIDS Entity Filters{Colors.END}',
        "Filter which data to process based on BIDS entities."
    )

    filters.add_argument(
        "--participant-label",
        metavar="LABEL",
        dest="participant_label",
        nargs='+',
        help="Process one or more participants (without 'sub-' prefix). "
             "Default: all subjects in dataset.",
    )

    filters.add_argument(
        "--task",
        metavar="TASK",
        help="Task to process (without 'task-' prefix). REQUIRED.",
    )

    filters.add_argument(
        "--session",
        metavar="SESSION",
        help="Session to process (without 'ses-' prefix). Default: None.",
    )

    filters.add_argument(
        "--space",
        metavar="SPACE",
        default="MNI152NLin2009cAsym",
        help="Template space for preprocessed data "
             "(default: 'MNI152NLin2009cAsym').",
    )

    # =========================================================================
    # MODELING OPTIONS
    # =========================================================================
    modeling = parser.add_argument_group(
        f'{Colors.BOLD}Modeling Options{Colors.END}',
        "Configure the first-level GLM design and contrasts."
    )

    modeling.add_argument(
        "--events-file",
        type=Path,
        metavar="PATH",
        dest="events_file",
        help="Path to a single global events TSV file. If omitted and "
             "INPUT_DIR is preprocessed data, --events-dir is required.",
    )

    modeling.add_argument(
        "--events-dir",
        type=Path,
        metavar="PATH",
        dest="events_dir",
        help="Directory to recursively search for participant-specific "
             "*_events.tsv files when INPUT_DIR is preprocessed data and "
             "--events-file is not provided.",
    )

    modeling.add_argument(
        "--confounds-file",
        type=Path,
        metavar="PATH",
        dest="confounds_file",
        help="Path to custom confounds TSV file. Default: the one found "
             "in fMRIPrep output.",
    )

    modeling.add_argument(
        "--conditions",
        metavar="COND",
        nargs='+',
        help="Trial types from events.tsv to include in the model. "
             "Default: all trial types found in events file.",
    )

    modeling.add_argument(
        "--regressors",
        metavar="REG",
        nargs='+',
        help="Confound regressors to include (column names from confounds "
             "file, e.g., trans_x trans_y rot_x). Default: none.",
    )

    modeling.add_argument(
        "--contrasts",
        metavar="CON",
        nargs='+',
        help="Contrast definitions as strings (e.g., 'cond1 - cond2'). "
             "Passed directly to nilearn's compute_contrast. "
             "Default: one identity contrast per condition.",
    )

    # =========================================================================
    # STATISTICAL THRESHOLDING
    # =========================================================================
    stats = parser.add_argument_group(
        f'{Colors.BOLD}Statistical Thresholding{Colors.END}',
        "Control statistical thresholding of output maps."
    )

    stats.add_argument(
        "--alpha",
        metavar="ALPHA",
        nargs='+',
        type=float,
        default=[0.05],
        help="Statistical p-value threshold(s). Passed to "
             "nilearn.glm.threshold_stats_img. Default: 0.05.",
    )

    stats.add_argument(
        "--height-control",
        metavar="METHOD",
        nargs='+',
        dest="height_control",
        help="Multiple comparison correction method(s): 'fdr', 'bonferroni', "
             "or 'none'. Paired with --alpha values. Default: none.",
    )

    return parser


def parse_derivatives_arg(derivatives_list):
    """
    Parse --derivatives KEY=PATH arguments into a dictionary.

    Expected format: key=path (e.g., preproc=/data/derivatives/fmriprep)
    """
    if not derivatives_list:
        return {}
    result = {}
    for item in derivatives_list:
        if "=" not in item:
            raise ValueError(
                f"Invalid --derivatives argument: '{item}'. "
                f"Expected format: key=path (e.g., preproc=/data/derivatives/fmriprep)"
            )
        key, path = item.split("=", 1)
        result[key.strip()] = Path(path.strip())
    return result


def main():
    """Main CLI entry point for trialix."""
    parser = create_parser()
    args = parser.parse_args()

    from trialix.pipeline import run_pipeline

    try:
        run_pipeline(args)
    except KeyboardInterrupt:
        print("\nInterrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"\n{Colors.RED}Error: {e}{Colors.END}", file=sys.stderr)
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
