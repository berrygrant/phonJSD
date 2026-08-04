## Resubmission

This is a resubmission. Compared to the previously submitted version (2.1.0):

- I replaced the `\\dontrun{}` example for `extract_mfcc()` with a fast,
  executable example. It creates a short WAV file under `tempdir()`, guards the
  optional `tuneR` dependency with `requireNamespace()`, and removes the file
  after use.
- I unwrapped the short `phontrast()` and `hier_boot_jsd_model()` bootstrap
  examples after measuring them below five seconds. I also reduced the
  illustrative `estimate_jsd()` bootstrap count from 50 to 10, bringing that
  example below five seconds so it could be unwrapped. No examples now use
  `\\dontrun{}` or `\\donttest{}`.
- I removed direct access to `.GlobalEnv` and `.Random.seed`. Seeded internal
  sampling now uses a private deterministic generator, so `eval_seed` remains
  reproducible without modifying the user's workspace or RNG stream.
- Versions 2.2.0 and 2.3.0 added a multivariate-normal density backend and
  distribution-aware plotting, respectively. Version 2.3.1 contains the CRAN
  review fixes above; the previously published GitHub tags remain unchanged.

## Submission

This is a new submission (first time on CRAN). `phontrast` is a rename and
reframing of the previously GitHub-only package `phonJSD` (released through
v1.2.0). Development, testing, and documentation were assisted by AI tools, as
disclosed in the README; all design, analysis, and release decisions were made
by the maintainer.

## Test environments

- Local: macOS 26.5.2 (aarch64-apple-darwin23), R 4.6.1,
  `R CMD check --as-cran --no-manual`

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for the reviewer

* The spell-check flags domain terminology in the Description and
  documentation (e.g. "Bhattacharyya", "Mahalanobis", "Pillai", "MFCCs",
  "sociophonetics", "phonological"). These are standard terms from phonetics
  and multivariate statistics; they are listed in `inst/WORDLIST`.
* The `Description` cites Lin (1991) <doi:10.1109/18.61115> for the
  Jensen-Shannon divergence.
* `ggplot2`, `mgcv`, and `tuneR` are used only conditionally
  (`requireNamespace()`), so they are in Suggests rather than Imports.
