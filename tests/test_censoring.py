import json
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from nilearn.glm.first_level import make_first_level_design_matrix

from trialix.censoring import prepare_censored_run
from trialix.cli import create_parser
from trialix.config import merge_config, validate_config
from trialix.model import build_first_level_model, fit_model
from trialix.outputs import save_contrast_maps, save_thresholded_map


def run_inputs():
    n = 160
    rng = np.random.default_rng(14)
    data = 100 + rng.normal(size=(3, 3, 3, n))
    # Encode volume identity in one voxel to catch selection/indexing errors.
    data[0, 0, 0] = np.arange(n)
    img = nib.Nifti1Image(data, np.eye(4))
    events = pd.DataFrame({
        'onset': [10., 40., 8., 18., 42.],
        'duration': [20., 20., 1., 1., 1.],
        'trial_type': ['focus', 'focus', 'task', 'task', 'task'],
    })
    confounds = pd.DataFrame({'motion': np.arange(n, dtype=float)})
    flm = build_first_level_model(2.)
    return img, events, confounds, flm


@pytest.mark.parametrize('drift', ['cosine', 'polynomial', None])
def test_alignment_hrf_and_reduced_drifts(drift):
    img, events, confounds, flm = run_inputs()
    flm.drift_model = drift
    task_events = events[events.trial_type == 'task']
    cropped, design, metadata = prepare_censored_run(
        img, task_events, confounds, events, 'focus', flm, 'sub-01_events.tsv'
    )
    indices = np.setdiff1d(np.arange(160), np.r_[5:15, 20:30])
    np.testing.assert_array_equal(cropped.get_fdata()[0, 0, 0], indices)
    np.testing.assert_allclose(design.motion, indices)
    original = make_first_level_design_matrix(
        np.arange(160) * 2., events=task_events, drift_model=None
    )
    np.testing.assert_allclose(design.task, original.task.iloc[indices])
    reduced_drifts = make_first_level_design_matrix(
        np.arange(140) * 2., drift_model=drift, drift_order=flm.drift_order,
        high_pass=flm.high_pass,
    )
    for col in reduced_drifts:
        np.testing.assert_allclose(design[col], reduced_drifts[col], atol=1e-12)
    assert metadata['RetainedVolumeIndices'] == indices.tolist()
    assert metadata['OriginalVolumeCount'] == 160
    assert metadata['RetainedVolumeCount'] == 140
    if drift == 'cosine':
        assert len(metadata['DriftColumns']) == 5  # 280s, rather than 320s.
        full = make_first_level_design_matrix(np.arange(160) * 2.)
        assert len([c for c in full if c.startswith('drift_')]) == 6


def test_real_fit_and_contrast():
    img, events, confounds, flm = run_inputs()
    cropped, design, _ = prepare_censored_run(
        img, events[events.trial_type == 'task'], confounds,
        events, 'focus', flm, 'events.tsv',
    )
    flm.mask_img = nib.Nifti1Image(np.ones(img.shape[:3], dtype=np.uint8), img.affine)
    fitted = fit_model(flm, cropped, events, confounds, design)
    assert fitted.design_matrices_[0].shape[0] == cropped.shape[3] == 140
    assert fitted.compute_contrast('task', output_type='z_score').shape == img.shape[:3]


def test_overlapping_windows_clipping_and_half_open_bounds():
    img, events, confounds, flm = run_inputs()
    windows = pd.DataFrame({
        'onset': [-5., 10., 310.], 'duration': [25., 20., 50.],
        'trial_type': ['focus'] * 3,
    })
    _, _, meta = prepare_censored_run(
        img, events[events.trial_type == 'task'], None,
        windows, 'focus', flm, 'events.tsv',
    )
    assert meta['RetainedVolumeIndices'] == list(range(15, 155))


@pytest.mark.parametrize('onset,duration', [(10, 0), (10, -1), (np.nan, 20),
                                           (10, np.inf), ('bad', 20)])
def test_invalid_windows(onset, duration):
    img, events, confounds, flm = run_inputs()
    windows = pd.DataFrame({'onset': [onset], 'duration': [duration],
                            'trial_type': ['focus']})
    with pytest.raises(ValueError, match='Censor'):
        prepare_censored_run(img, events, confounds, windows, 'focus', flm, 'events.tsv')


def test_missing_empty_and_confounds_mismatch():
    img, events, confounds, flm = run_inputs()
    with pytest.raises(ValueError, match='not found'):
        prepare_censored_run(img, events, None, events, 'absent', flm, 'events.tsv')
    windows = pd.DataFrame({'onset': [0.], 'duration': [1000.], 'trial_type': ['focus']})
    with pytest.raises(ValueError, match='at least two'):
        prepare_censored_run(img, events, None, windows, 'focus', flm, 'events.tsv')
    with pytest.raises(ValueError, match='Confounds have'):
        prepare_censored_run(img, events, confounds.iloc[:-1], events, 'focus', flm, 'events.tsv')


def test_late_sleep_window_removed_at_fractional_tr():
    # Mirrors the reported dataset's timing and header precision, without data.
    tr = 2.4600000381469727
    n = 2713
    img = nib.Nifti1Image(np.arange(n, dtype=float).reshape(1, 1, 1, n), np.eye(4))
    events = pd.DataFrame({
        'onset': [10., 1000., 3000., 4900., 5033.16],
        'duration': [4.92, 4.92, 4.92, 4.92, 1650.66],
        'trial_type': ['task'] * 4 + ['consolidatedSleep'],
    })
    flm = build_first_level_model(tr)
    cropped, design, meta = prepare_censored_run(
        img, events, pd.DataFrame({'motion': np.arange(n, dtype=float)}),
        events, 'consolidatedSleep', flm, 'events.tsv',
    )
    assert cropped.shape[-1] == len(design) == meta['RetainedVolumeCount'] == 2046
    assert meta['RemovedVolumeCount'] == 667
    assert meta['RetainedVolumeIndices'] == list(range(2046))
    assert meta['RemovedVolumeIndices'] == list(range(2046, n))
    assert meta['Operation'] == 'remove_matching_intervals'
    assert 'consolidatedSleep' not in design
    np.testing.assert_allclose(design.motion, np.arange(2046))
    original = make_first_level_design_matrix(
        np.arange(n) * tr, events=events[events.trial_type == 'task'], drift_model=None
    )
    np.testing.assert_allclose(design.task, original.task.iloc[:2046], atol=1e-10)
    # Task responses still appear near all four original onsets, including late ones.
    for onset in [10., 1000., 3000., 4900.]:
        window = (design.index >= onset) & (design.index < onset + 30)
        assert design.loc[window, 'task'].max() > .1


def test_removed_only_condition_fails_clearly():
    img, events, _, flm = run_inputs()
    with pytest.raises(ValueError, match='No modeled events remain'):
        prepare_censored_run(img, events[events.trial_type == 'focus'], None,
                             events, 'focus', flm, 'events.tsv')


def test_cli_config_precedence():
    parser = create_parser()
    args = parser.parse_args(['input', 'output', 'participant'])
    assert merge_config(args, {'censor_trials': 'focus'})['censor_trials'] == 'focus'
    args = parser.parse_args(['input', 'output', 'participant', '--censor-trials', 'other'])
    assert merge_config(args, {'censor_trials': 'focus'})['censor_trials'] == 'other'
    for invalid in ['', [], 1]:
        with pytest.raises(ValueError, match='censor_trials'):
            validate_config({'censor_trials': invalid})


def test_confounds_discovery_matches_bold_run():
    from trialix.bids_utils import find_confounds_file

    calls = []

    def get(**filters):
        calls.append(filters)
        return ['run-02_confounds.tsv'] if filters['suffix'] == 'regressors' else []

    path = find_confounds_file(SimpleNamespace(get=get), '01', 'motor', run=2)
    assert path == 'run-02_confounds.tsv'
    assert all(filters['run'] == 2 for filters in calls)


def test_sidecars(tmp_path):
    img, events, _, flm = run_inputs()
    _, _, meta = prepare_censored_run(
        img, events, None, events, 'focus', flm, 'sub-01_events.tsv'
    )
    stat = nib.Nifti1Image(np.ones((3, 3, 3)), np.eye(4))
    metadata = {'TemporalCensoring': meta}
    paths = save_contrast_maps({'z_score': stat}, tmp_path, '01', None,
                              'motor', 'MNI', 'task', metadata=metadata)
    paths.append(save_thresholded_map(stat, tmp_path, '01', None, 'motor',
                                     'MNI', 'task', .05, 'fdr', metadata=metadata))
    for path in paths:
        sidecar = path.with_name(path.name.removesuffix('.nii.gz') + '.json')
        saved = json.loads(sidecar.read_text())
        assert saved['TemporalCensoring'] == meta
        assert saved['Contrast'] == 'task'


def test_uncensored_fit_keeps_existing_path():
    calls = []
    flm = SimpleNamespace(
        fit=lambda *args, **kwargs: calls.append((args, kwargs)),
        design_matrices_=[pd.DataFrame({'task': [1.]})],
    )
    fit_model(flm, 'bold.nii.gz', 'events', 'confounds')
    assert calls == [(('bold.nii.gz',), {'events': 'events', 'confounds': 'confounds'})]


@pytest.mark.parametrize('subject,onset', [('01', 10.), ('02', 40.)])
def test_pipeline_subject_events_and_output_provenance(tmp_path, monkeypatch, subject, onset):
    from trialix import pipeline

    img, events, confounds, flm = run_inputs()
    bold = tmp_path / f'sub-{subject}_bold.nii.gz'
    nib.save(img, bold)
    # A subject's marker is independent of the conditions included in the GLM.
    events = events[events.trial_type == 'task'].copy()
    events.loc[len(events) + 10] = [onset, 20., 'focus']
    event_file = tmp_path / f'sub-{subject}_events.tsv'
    events.to_csv(event_file, sep='\t', index=False)
    confound_file = tmp_path / 'confounds.tsv'
    confounds.to_csv(confound_file, sep='\t', index=False)
    flm.mask_img = nib.Nifti1Image(np.ones(img.shape[:3], dtype=np.uint8), img.affine)
    monkeypatch.setattr(pipeline, 'find_bold_file', lambda *args: bold)
    monkeypatch.setattr(pipeline, 'get_repetition_time', lambda *args: 2.)
    monkeypatch.setattr(pipeline, 'find_events_file_in_dir',
                        lambda root, sub, task, **kwargs: root / f'sub-{sub}_events.tsv')
    monkeypatch.setattr(pipeline, 'build_first_level_model', lambda **kwargs: flm)
    monkeypatch.setattr(pipeline, 'generate_report', lambda *args, **kwargs: None)
    monkeypatch.setattr(pipeline, 'save_report', lambda *args: None)
    monkeypatch.setattr(pipeline, 'threshold_contrast_map',
                        lambda img, **kwargs: (img, 1.))
    config = {
        'events_dir': tmp_path, 'conditions': ['task'], 'censor_trials': 'focus',
        'confounds_file': confound_file, 'regressors': ['motion'],
        'debug': False, 'output_dir': tmp_path / 'output', 'contrasts': ['task'],
    }
    layout = SimpleNamespace(parse_file_entities=lambda path: {'run': 1})
    pipeline._process_subject_session(
        config, None, layout, False, subject, None, 'motor', 'MNI', {}, [(0.05, None)]
    )
    assert list(flm.design_matrices_[0].columns[:2]) == ['task', 'motion']
    assert 'focus' not in flm.design_matrices_[0]
    sidecars = list((tmp_path / 'output').rglob('*.json'))
    assert len(sidecars) == 6  # five contrast outputs plus thresholded z map
    for path in sidecars:
        meta = json.loads(path.read_text())['TemporalCensoring']
        assert meta['EventsFile'] == str(event_file)
        assert meta['RetainedVolumeIndices'] == [i for i in range(160) if not int(onset / 2) <= i < int(onset / 2) + 10]
