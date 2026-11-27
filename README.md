<div align="center">

# Trialix

**First-level fMRI analysis for trial-based experimental designs**

[Installation](#installation) | [Quick Start](#quick-start) | [Usage](#usage) | [Configuration](#configuration) | [Workflow](#workflow) | [Output Structure](#output-structure) | [Events File Format](#events-file-format) | [Expected fMRIPrep Output Structure](#expected-fmriprep-output-structure) | [License](#license)

</div>

---

Trialix is a BIDS-application CLI tool that performs participant-level (first-level) fMRI analysis for trial-based experimental conditions — block designs, event-related designs, or mixed designs.

At its core, Trialix is a BIDS wrapper for [nilearn](https://nilearn.github.io/). GLM fitting, contrast computation, and the existing FDR/Bonferroni thresholding use Nilearn, in particular [`FirstLevelModel`](https://nilearn.github.io/stable/modules/generated/nilearn.glm.first_level.FirstLevelModel.html). Trialix adds a BIDS-aware layer on top, handling data discovery, confound selection, and BIDS-compliant output generation so that users can run first-level analyses directly from fMRIPrep derivatives with minimal boilerplate. Version **0.2.0** adds optional participant-level height-FWE correction using [pyrft](https://github.com/arovai/pyrft).

> **Note:** Trialix does **not** perform fMRI preprocessing. It expects data preprocessed by [fMRIPrep](https://fmriprep.org/) or a similar tool. Group-level (second-level) analysis is handled by the companion tool **StatCraft**.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/arovai/trialix.git
cd trialix

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install trialix and its dependencies
pip install -e .

# Optional: install the pinned pyrft backend for RFT correction
pip install -e '.[rft]'
```

### Dependencies

- Python ≥ 3.10
- [nilearn](https://nilearn.github.io/) ≥ 0.10
- [pybids](https://bids-standard.github.io/pybids/) ≥ 0.16
- [nibabel](https://nipy.org/nibabel/) ≥ 4.0
- numpy, scipy ≥ 1.9, pandas, pyyaml, matplotlib
- Optional `[rft]` extra: pyrft 0.2.0 with its NIfTI adapter, pinned to
  [commit b635329](https://github.com/arovai/pyrft/commit/b63532940f78999850132dcdc1bec9317b98058d).
  This currently requires Git/network access at installation because pyrft is not
  yet published on PyPI. Ordinary analyses do not require pyrft.

---

## Quick Start

### Basic usage (BIDS rawdata + fMRIPrep derivatives)

```bash
trialix /data/rawdata /data/output participant \
    --task motor \
    --derivatives preproc=/data/derivatives/fmriprep
```

### Preprocessed input with participant-specific events

```bash
trialix /data/derivatives/fmriprep /data/output participant \
  --task motor \
  --events-dir /data/rawdata
```

### Full analysis with contrasts and thresholding

```bash
trialix /data/rawdata /data/output participant \
    --task motor \
    --derivatives preproc=/data/derivatives/fmriprep \
    --participant-label 01 02 03 \
    --conditions left_hand right_hand \
    --regressors trans_x trans_y trans_z rot_x rot_y rot_z \
    --contrasts "left_hand - right_hand" "left_hand + right_hand" \
    --alpha 0.05 0.001 \
    --height-control fdr bonferroni
```

### Using a config file

```bash
trialix /data/rawdata /data/output participant --config config.yaml
```

CLI arguments always override config file values.

---

## Usage

```
trialix INPUT_DIR OUTPUT_DIR participant [options]
```

### Positional Arguments

| Argument | Description |
|---|---|
| `INPUT_DIR` | Path to BIDS rawdata folder (with `--derivatives`) or to preprocessed data folder (with `--events-file`/`--events-dir`). |
| `OUTPUT_DIR` | Path to output directory for analysis derivatives. |
| `participant` | Analysis level (currently only `participant` is supported). |

### Input Modes

**Mode 1: BIDS rawdata + derivatives**
- `INPUT_DIR` points to the BIDS rawdata folder
- Use `--derivatives preproc=PATH` to specify fMRIPrep output
- Events files are discovered automatically from rawdata

**Mode 2: Preprocessed data directly**
- `INPUT_DIR` points to the preprocessed data folder (e.g., fMRIPrep output)
- Use `--events-file PATH` for one global events TSV file, or
- Use `--events-dir PATH` to recursively find participant-specific `*_events.tsv`

### General Options

| Option | Description |
|---|---|
| `--help` | Show help message and exit. |
| `--version` | Show version and exit. |
| `--debug` | Enable debug mode with increased verbosity. |
| `--config FILE` | Path to YAML configuration file. |

### Derivatives

| Option | Description |
|---|---|
| `--derivatives KEY=PATH` | Specify BIDS derivatives location. Use `preproc=PATH` for preprocessed data. Can be specified multiple times. |

### BIDS Entity Filters

| Option | Default | Description |
|---|---|---|
| `--participant-label LABEL [LABEL ...]` | all | Subjects to process (without `sub-` prefix). |
| `--task TASK` | *required* | Task to process (without `task-` prefix). |
| `--session SESSION` | None | Session to process (without `ses-` prefix). |
| `--space SPACE` | `MNI152NLin2009cAsym` | Template space. |

### Modeling Options

| Option | Default | Description |
|---|---|---|
| `--events-file PATH` | auto | Path to a single global events TSV file (Mode 2). |
| `--events-dir PATH` | auto | Recursive search root for participant-specific `*_events.tsv` (Mode 2 when `--events-file` is omitted). |
| `--confounds-file PATH` | auto | Custom confounds TSV file. Default: fMRIPrep output. |
| `--conditions COND [COND ...]` | all | Trial types to include from events.tsv. |
| `--censor-trials TRIAL_TYPE` | disabled | Remove volumes acquired within matching events: `[onset, onset + duration)`. Also configurable as `censor_trials`. |
| `--regressors REG [REG ...]` | none | Confound regressors (column names from confounds file). |
| `--contrasts CON [CON ...]` | auto | Contrast definitions (e.g., `"cond1 - cond2"`). |

### Temporal censoring

**Behavior change:** `--censor-trials` removes matching intervals. The initial
implementation kept them; rerun affected analyses with this corrected version.

Use `--censor-trials consolidatedSleep` to remove volumes whose acquisition times
fall inside `consolidatedSleep` events in the resolved events TSV. The same discovery
rules apply, including participant-, session-, and run-specific events files.
The censoring windows are read before filtering events with `--conditions`.

```bash
trialix /data/derivatives/fmriprep /data/output participant \
    --task motor --events-dir /data/rawdata \
    --censor-trials consolidatedSleep --conditions left_hand right_hand \
    --regressors trans_x trans_y trans_z
```

For an event with onset `10` and duration `120`, the removed interval is
`10 <= time < 130` seconds. Volume `i` is sampled at `i * TR` (zero-based).
All matching windows are removed; overlapping windows never remove a volume twice.
Windows extending beyond the scan affect only existing volumes. Missing trial
types, invalid or nonpositive durations, nonfinite timings, fewer than two
retained volumes, and confound/BOLD length mismatches cause that analysis to fail
with an explanatory log message.

BOLD and every selected confound use exactly the same volume indices. Task
regressors are convolved on the **original timeline**, then the same rows are
selected, preserving HRF contributions from events before a retained window.
Events are not shifted or compressed before convolution. Drift terms and the
constant are generated by Nilearn on a new regular timeline containing only the
retained volumes, so cosine drift count depends on the retained length. Disjoint
windows are concatenated at the original TR for drift generation and AR(1)
fitting; gaps are treated as adjacent samples, with no extra HRF padding.
Any precomputed confound columns (including fMRIPrep cosine columns) are selected
as supplied, rather than recomputed; use `regressors` to choose them deliberately.

`--conditions` controls modeled trial types independently. The censoring trial
type is automatically excluded from the model, even if listed in `--conditions`.
In YAML, use `censor_trials: consolidatedSleep` (or `null` to disable it); a supplied
CLI `--censor-trials` value takes precedence.

Each output statistical map, including thresholded maps, has a JSON sidecar
with `TemporalCensoring` provenance: enabled status, trial type, source events
file, intervals, TR, original/retained/removed counts, zero-based retained volume
indices and original acquisition times, removed indices, the explicit removal
operation, excluded model trial type, task-regressor timing strategy, and drift
settings/columns. Without
censoring the sidecar records `Enabled: false`.

### Statistical Thresholding

| Option | Default | Description |
|---|---|---|
| `--alpha ALPHA [ALPHA ...]` | `0.05` | P-value threshold(s). |
| `--height-control METHOD [METHOD ...]` | none | Correction method: `fdr`, `bonferroni`, `rft` (requires `[rft]` extra), or `none`. |

### Participant-level RFT correction

Install the optional backend from the checkout, then request `rft`:

```bash
python -m pip install -e '.[rft]'
trialix /data/rawdata /data/output participant \
    --task motor --derivatives preproc=/data/derivatives/fmriprep \
    --contrasts "left_hand - right_hand" --height-control rft --alpha 0.05
```

YAML configuration is equivalent:

```yaml
alpha: [0.05]
height_control: [rft]
```

You can request multiple methods, for example `--height-control rft fdr` with
`--alpha 0.05 0.05`. Existing element-wise pairing is retained: a missing alpha
reuses the last value; a missing method means `none`. CLI alpha overrides YAML
only when explicitly supplied; otherwise YAML is respected before the 0.05 fallback.
RFT requires `0 < alpha <= 0.1`. It is opt-in; the default correction is unchanged.

For each fitted participant/session run, Trialix:

1. Reconstructs the complete **whitened residual** stack from Nilearn's AR groups
   in fitted-mask voxel order. It does not use the public residual image, observed
   Z/T-map smoothness, or the preprocessing kernel as a substitute.
2. Obtains the residual degrees of freedom and observation count from the fitted
   model, including any trial censoring. Smoothness is estimated once through pyrft
   and reused across the run's contrasts and alpha values.
3. Computes resels from the fitted model's actual search mask and estimated
   anisotropic FWHM. Derivative estimation uses complete in-mask forward cubes.
4. Applies two-sided **T-field height FWE** correction: alpha/2 per tail, using
   fitted residual df. No automatic Bonferroni fallback is used for `rft`.
5. Saves the original Z scores only where native absolute T exceeds its RFT cutoff.
   Output naming remains compatible with existing downstream consumers, for example
   `*_desc-zScoreThreshRftAlpha0p05_statmap.nii.gz`.

FWE is controlled spatially for each contrast; there is no additional correction
across the collection of requested contrasts.

The integration supports one fitted volumetric run and named one-dimensional
T contrasts. It does not implement F-field correction, joint multi-run RFT,
surfaces, or cluster-extent inference. The full residual stack needs additional
memory while estimating smoothness. Zero residual variance, too few usable scans,
insufficient derivative stencils, unsupported geometry, inconsistent df or backend
errors fail the participant analysis; they do not trigger map-only estimation or
another correction method. An RFT run with any failed participant/session exits
nonzero after processing the remaining participants.

RFT map sidecars include a `Thresholding` block and an `RFT` block with the actual
T/Z cutoffs, full/per-tail alpha, fitted df, residual preparation, observation
count/design rank, mask geometry, FWHM in voxels/mm, resels, backend version and
any diagnostics. Existing `TemporalCensoring` provenance is preserved.

The HTML report uses the **first requested alpha/method pair**. When that method
is RFT, Trialix supplies its equivalent numerical Z cutoff to Nilearn and labels
the report as RFT height-FWE correction. A companion `*_report.json` records the
RFT provenance. Cluster sizes in this report are descriptive; no cluster-FWE
probabilities are claimed. When another method comes first, its report behavior
is unchanged; the RFT NIfTI files and sidecars are still saved.

Numerical pyrft–SPM comparisons are maintained separately in
[pyrft-spm-validation](https://github.com/arovai/pyrft-spm-validation).
Scientific calibration remains deferred. Per `INSTRUCTIONS.md`, this Trialix
integration has only received syntax/static source review; runtime analyses and
integration tests have not been executed.

### Config-File-Only Parameters

These parameters can only be set in the YAML configuration file and are passed directly to nilearn's `FirstLevelModel`:

| Parameter | Default | Description |
|---|---|---|
| `hrf_model` | `glover` | HRF model (`glover`, `spm`, `glover + derivative`, etc.) |
| `drift_model` | `cosine` | Drift model (`cosine`, `polynomial`, or `null`) |
| `drift_order` | `1` | Drift order for polynomial model |
| `high_pass` | `0.01` | High-pass filter cutoff in Hz |
| `smoothing_fwhm` | `null` | Spatial smoothing FWHM in mm |

---

## Configuration

Trialix supports a YAML configuration file. CLI arguments override config values. See [`config_example.yaml`](config_example.yaml) for a complete example.

```yaml
# Required
task: motor

# Optional BIDS filters
participant_label:
  - "01"
  - "02"
space: MNI152NLin2009cAsym

# Modeling
conditions:
  - left_hand
  - right_hand
regressors:
  - trans_x
  - trans_y
  - trans_z
contrasts:
  - "left_hand - right_hand"

# Thresholding
alpha:
  - 0.05
height_control:
  - fdr

# FirstLevelModel parameters (config-only)
hrf_model: glover
drift_model: cosine
high_pass: 0.01
smoothing_fwhm: 6
```

---

## Workflow

1. **Validate inputs** — Check CLI arguments and config file; verify data paths.
2. **Discover data** — Use pybids to find subjects, tasks, sessions, BOLD files, events, and confounds.
3. **Loop over subjects and sessions:**
  - Load preprocessed BOLD data
  - Resolve events TSV
  - Rawdata mode: pybids lookup in rawdata (participant-specific first, then global)
  - Preprocessed mode: `--events-file` (global) or recursive `--events-dir` lookup by BIDS entities
  - Load events TSV (filter by requested conditions)
  - Load confounds TSV (select requested regressors)
  - Extract TR from BIDS metadata or NIfTI header
  - Build `FirstLevelModel` with specified parameters
  - Fit the GLM
  - Compute contrasts (z-score, t-stat, p-value, effect size, effect variance)
  - Apply statistical thresholding (optionally residual-based pyrft RFT height FWE)
  - Save BIDS-compliant NIfTI outputs
  - Generate and save HTML report
4. **Write `dataset_description.json`** in output directory.

---

## Output Structure

```
OUTPUT_DIR/
├── dataset_description.json
├── logs/
│   └── trialix_YYYYMMDD_HHMMSS.log
├── sub-01/
│   ├── [ses-SES/]
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-zScore_statmap.nii.gz
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-tStat_statmap.nii.gz
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-pValue_statmap.nii.gz
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-effectSize_statmap.nii.gz
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-effectVariance_statmap.nii.gz
│   │   ├── sub-01_[ses-SES_]task-TASK_space-SPACE_contrast-CON_desc-zScoreThresh*_statmap.nii.gz
│   │   └── sub-01_[ses-SES_]task-TASK_space-SPACE_report.html
│   └── ...
├── sub-02/
│   └── ...
└── ...
```

---

## Events File Format

The events TSV file must contain at least three columns:

| Column | Type | Description |
|---|---|---|
| `onset` | float | Event onset time in seconds |
| `duration` | float | Event duration in seconds |
| `trial_type` | string | Condition label |

Example:

```tsv
onset	duration	trial_type
0.0	15.0	left_hand
20.0	15.0	right_hand
40.0	15.0	left_hand
60.0	15.0	right_hand
```

---

## Expected fMRIPrep Output Structure

Trialix expects fMRIPrep output organized as:

```
PREPROC_PATH/
├── sub-SUB/
│   ├── [ses-SES/]
│   │   └── func/
│   │       ├── sub-SUB_[ses-SES_]task-TASK_space-SPACE_desc-preproc_bold.nii.gz
│   │       └── sub-SUB_[ses-SES_]task-TASK_desc-confounds_timeseries.tsv
│   └── ...
└── ...
```

---

## Acknowledgments

Trialix is built on top of [nilearn](https://nilearn.github.io/), which provides the core neuroimaging statistical routines used throughout this tool. We gratefully acknowledge the nilearn developers and community for their work in making high-quality, well-tested fMRI analysis methods accessible in Python.

Trialix also relies on [pybids](https://bids-standard.github.io/pybids/) for BIDS dataset querying, [nibabel](https://nipy.org/nibabel/) for neuroimaging file I/O, and [fMRIPrep](https://fmriprep.org/) as the expected preprocessing pipeline.

---

## License

MIT
