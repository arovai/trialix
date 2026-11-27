# Trialix session memory

## Participant-level RFT integration (2026-11-27)

- Version moves from 0.1.0 to 0.2.0; Python minimum becomes 3.10.
- Add opt-in `--height-control rft` / YAML `height_control: [rft]`.
- Optional `[rft]` extra pins pyrft 0.2.0 to Git commit
  b63532940f78999850132dcdc1bec9317b98058d while PyPI publication is deferred.
- Reconstruct the full whitened residual stack from the fitted run's AR groups
  using labels_ and results_; do not use the public unwhitened residual image.
- Get residual degrees of freedom from fitted regression results, not from the
  original uncensored scan count. Require uniform df across groups and one run.
- Standardize/estimate global roughness through pyrft once per run. Use the
  fitted model mask. Threshold native T contrasts with alpha/2 per tail.
- Retain original Z scores at surviving T voxels to preserve output conventions.
- Record actual T/Z thresholds, residual preparation, widths, resels, masks,
  backend version and censoring provenance in map/report JSON sidecars.
- RFT HTML reports use a precomputed Z cutoff and explicit RFT title; never pass
  an unsupported `rft` control method to Nilearn or substitute an uncorrected map.
- RFT failure marks the participant failed and produces a nonzero final status.
- CLI alpha defaults now defer to YAML before falling back to 0.05.

## Verification restriction

INSTRUCTIONS.md explicitly prohibits running the tool or writing/running unit
or integration tests and allows syntax checks only. No analysis, CLI invocation,
installation/build check, or test suite has been run for this change. Python AST
and TOML syntax parsing and whitespace checks have passed for the implementation.

## Remaining work

- Complete user-facing usage/release documentation.
- User-side runtime checks should cover ordinary and censored AR(1) fits, YAML
  overrides, multiple contrasts/alphas, fixed-cutoff HTML reports, missing extra,
  degenerate residuals and RFT error propagation. These have NOT been executed.
- Scientific FWE calibration remains a separate deferred project.
