## Resubmission

The previous 2.4.0 submission was flagged by the incoming pretest with one
NOTE on the Windows flavor: the `estimate_jsd()` example exceeded 10s elapsed
(38s). The example ran two two-dimensional bootstrap demonstrations whose
per-resample multivariate bandwidth selection is slow on that platform. The
example now demonstrates bootstrapping on a single feature (univariate
bandwidth selection) and grouped estimation as a point estimate; it runs in
well under one second on the test machines below. No other changes relative
to the previous submission.

## Submission

This is an update to the CRAN release `phontrast` 2.3.1.

Version 2.4.0 corrects public release metadata and adds one opt-in feature.
The 2.4.0 GitHub tag initially retained `Version: 2.3.1`, and the package
citation remained frozen at version 2.0.0; the package version, NEWS, tag,
version-specific Zenodo DOI, and version-aware citation are now consistent.
If this update arrives close enough to 2.3.1 to trigger the
days-since-last-update note, the short interval is to correct the public
version and citation metadata that users currently see.

## Changes

- Added the opt-in proportion-standardized Pillai estimates documented in
  `NEWS.md` (`pillai_overlap(proportion_standardized = TRUE)`); existing
  defaults and the original return fields are unchanged.
- Corrected the package version metadata to 2.4.0.
- Replaced the stale citation with a version-aware citation for Grant M. Berry
  and the 2.4.0 Zenodo DOI <doi:10.5281/zenodo.21864533>.
- Added tests that keep the installed citation version and DOI synchronized
  with `DESCRIPTION`.

## Test environments

- Local: Ubuntu 24.04.4, R 4.3.3, `R CMD check --no-manual --timings`
  (offline build environment; the network-dependent CRAN incoming checks
  could not run there)
- R-hub v2 (GitHub Actions): Windows Server 2022, R-devel
  (2026-08-19 r90430 ucrt); Ubuntu, R-devel — both `R CMD check` clean,
  test suite 401 passing / 0 failures
- GitHub Actions: Ubuntu latest, R release, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

The sole local note is environmental: the suggested package `tuneR` is not
installable in the offline build container, so it was unavailable for
checking there. It is available on CRAN's machines, and the previous
submission's Debian pretest ran fully OK against the same Suggests list.
After the example revision above, `estimate_jsd`'s example runs in 0.7s
locally and no example exceeds the check timing thresholds.

The local PDF-manual check was skipped (no TeX installation); all Rd checks,
examples, tests, and vignettes completed successfully.
