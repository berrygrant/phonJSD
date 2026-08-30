## Submission

This update fixes the test failure reported by the CRAN additional checks
(tests-MKL flavor: R-devel with Intel MKL) for `phontrast` 2.4.0.

One test expected the package's documented "nonsingular" error for a
collinear two-feature design, but the underlying guard decided singularity
by whether `chol()` of the within-class error SSCP threw an error. That is
BLAS-dependent on exactly singular input: reference LAPACK errors where MKL
can return a tiny positive pivot, so under MKL the degenerate design slipped
past the guard and failed later inside `summary.manova()` with a different
message. The guard now decides rank deficiency by R's tolerance-based QR of
the residual matrix -- the same criterion `summary.manova()` applies -- so
the documented error is raised identically on every BLAS build. A
near-collinear regression test pins the behavior. Estimates on well-posed
designs are unchanged, and there are no other changes.

## Test environments

- Local: Ubuntu 24.04.4, R 4.3.3 with reference BLAS/LAPACK, and the full
  test suite additionally run against OpenBLAS (offline build environment;
  the network-dependent CRAN incoming checks could not run there)
- R-hub v2 (GitHub Actions): Windows R-devel; Ubuntu R-devel
- GitHub Actions: Ubuntu latest, R release, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

The sole local note is environmental: the suggested package `tuneR` is not
installable in the offline build container, so it was unavailable for
checking there. It is available on CRAN's machines.

The local PDF-manual check was skipped (no TeX installation); all Rd checks,
examples, tests, and vignettes completed successfully.
