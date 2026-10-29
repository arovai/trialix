"""Synchronized volume selection with original-time task regressors."""

import numpy as np
import pandas as pd
from nilearn.glm.first_level import make_first_level_design_matrix
from nilearn.image import index_img, load_img


def prepare_censored_run(bold_img, events, confounds, selection_events,
                         trial_type, flm, events_file):
    """Crop BOLD and regressors identically, rebuilding drifts after selection.

    Task convolution uses the original regular acquisition grid. Drift generation
    uses a new regular grid of retained volumes, concatenating disjoint windows.
    Selection events must be unfiltered, including when --conditions is supplied.
    """
    img = load_img(bold_img)
    if len(img.shape) != 4:
        raise ValueError("Trial censoring requires a 4D BOLD image.")
    n_scans = img.shape[3]
    t_r = float(flm.t_r)
    if not np.isfinite(t_r) or t_r <= 0:
        raise ValueError("Trial censoring requires a finite, positive TR.")
    if confounds is not None and len(confounds) != n_scans:
        raise ValueError(
            f"Confounds have {len(confounds)} rows but BOLD has {n_scans} volumes."
        )
    windows = selection_events.loc[
        selection_events['trial_type'] == trial_type, ['onset', 'duration']
    ]
    if windows.empty:
        raise ValueError(f"Censor trial_type '{trial_type}' not found in {events_file}.")
    try:
        values = windows.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Censor onsets and durations must be numeric.") from exc
    if not np.isfinite(values).all() or (values[:, 1] <= 0).any():
        raise ValueError("Censor onsets must be finite and durations finite and positive.")
    times = np.arange(n_scans) * t_r
    keep = np.zeros(n_scans, dtype=bool)
    intervals = []
    for onset, duration in values:
        end = onset + duration
        if not np.isfinite(end):
            raise ValueError("Censor interval end must be finite.")
        keep |= (times >= onset) & (times < end)
        intervals.append({'Onset': float(onset), 'Duration': float(duration),
                          'End': float(end)})
    indices = np.flatnonzero(keep)
    if len(indices) < 2:
        raise ValueError("Trial censoring must retain at least two BOLD volumes.")

    # Do not compress event timing before HRF convolution: preceding events can
    # still contribute to retained volumes through their haemodynamic response.
    full_design = make_first_level_design_matrix(
        times, events=events, hrf_model=flm.hrf_model, drift_model=None,
        fir_delays=flm.fir_delays, min_onset=flm.min_onset,
    )
    regressors = full_design.drop(columns='constant').iloc[indices].reset_index(drop=True)
    if confounds is not None:
        regressors = pd.concat(
            [regressors, confounds.iloc[indices].reset_index(drop=True)], axis=1
        )
    reduced_times = np.arange(len(indices)) * t_r
    design = make_first_level_design_matrix(
        reduced_times, drift_model=flm.drift_model, drift_order=flm.drift_order,
        high_pass=flm.high_pass, add_regs=regressors,
    )
    metadata = {
        'Enabled': True,
        'TrialType': trial_type,
        'EventsFile': str(events_file),
        'IntervalConvention': 'onset <= volume_index * RepetitionTime < onset + duration',
        'Intervals': intervals,
        'RepetitionTime': t_r,
        'OriginalVolumeCount': int(n_scans),
        'RetainedVolumeCount': int(len(indices)),
        'RemovedVolumeCount': int(n_scans - len(indices)),
        'RetainedVolumeIndices': indices.tolist(),
        'VolumeIndexBase': 0,
        'TaskRegressors': 'HRF convolution on original timeline, then volume selection',
        'DriftTimeline': 'Retained volumes concatenated in acquisition order at original TR',
        'DriftModel': flm.drift_model,
        'HighPass': flm.high_pass,
        'DriftOrder': flm.drift_order,
        'DriftColumns': [c for c in design if c not in regressors.columns and c != 'constant'],
    }
    return index_img(img, indices), design, metadata
