# Changelog

## 0.2.0

### Participant-level RFT

- Add opt-in `--height-control rft` and matching YAML configuration.
- Provide an optional `[rft]` extra pinned to pyrft 0.2.0's Git revision.
- Reconstruct full whitened residuals from fitted AR groups and estimate
  anisotropic smoothness once per participant/session run.
- Apply two-sided native T-field height-FWE thresholds using the fitted mask
  and residual df, including after temporal censoring; preserve Z-score outputs.
- Record RFT cutoffs, smoothness, resels, mask geometry and residual preparation
  in JSON sidecars. Use the actual cutoff and RFT title in HTML reports.
- Fail requested RFT processing explicitly rather than silently omitting maps
  or substituting another method. RFT participant failures give a nonzero status.

### Configuration and requirements

- Require Python >= 3.10 and declare SciPy >= 1.9 explicitly.
- Respect YAML alpha when no CLI alpha is supplied; retain default alpha 0.05.
- Accept scalar YAML alpha/method values as well as lists.
- Keep existing correction defaults and non-RFT thresholding algorithms unchanged.

### Verification

Repository instructions allow syntax checks only. Python/TOML/YAML syntax and
whitespace were checked; no Trialix run or unit/integration suite was executed.
Scientific FWE calibration is deferred. Runtime integration remains to be checked
in the user's analysis environment.
