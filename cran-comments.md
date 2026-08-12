## Submission

This is an update to the CRAN release `phontrast` 2.3.1.

Version 2.4.0 is being submitted three days after 2.3.1 because it corrects
release metadata that users currently see: the GitHub 2.4.0 tag initially
retained `Version: 2.3.1`, and the package citation remained frozen at version
2.0.0. The tag, package version, NEWS, version-specific Zenodo DOI, and
version-aware package citation are now consistent.

## Changes

- Added the opt-in proportion-standardized Pillai estimates documented in
  `NEWS.md`; existing defaults and the original return fields are unchanged.
- Corrected the package version metadata to 2.4.0.
- Replaced the stale citation with a version-aware citation for Grant M. Berry
  and the 2.4.0 Zenodo DOI <doi:10.5281/zenodo.21864533>.
- Added tests that keep the installed citation version and DOI synchronized
  with `DESCRIPTION`.

## Test environments

- Local: macOS 26.5.2 (aarch64-apple-darwin23), R 4.6.1,
  `R CMD check --as-cran --no-manual`
- GitHub Actions: Ubuntu latest, R release, `R CMD check`

## R CMD check results

0 errors | 0 warnings | 1 note

The sole note is the expected incoming-check note:

```
Days since last update: 3
```

The unusually short interval is necessary to correct the public package
version and citation metadata described above.

The local PDF-manual check could not run because the local TeX installation is
missing `inconsolata.sty`. All Rd checks, examples, tests, vignettes, and the
HTML manual check completed successfully.
