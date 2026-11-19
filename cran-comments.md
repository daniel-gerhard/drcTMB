## Resubmission

This is a resubmission. In this version I have:

* Reduced the check time (NOTE on r-devel-windows-x86_64: overall checktime 12 min): the slow
  tests (bootstrap, profile likelihood and Monte Carlo integration with random effects) are now
  skipped on CRAN with testthat::skip_on_cran(), which reduces the time of the tests by about 70%,
  and the parametric bootstrap test in the vignette uses fewer simulations.

In the previous resubmission I had:

* Shortened and simplified the Description field in DESCRIPTION, removing the words flagged as
  possibly misspelled (names of models, authors and organisations).

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.6.2
* win-builder: Windows, R-devel (2026-11-17 r90926 ucrt)
* R-hub (GitHub Actions): Windows, R-devel; macOS arm64, R-devel; macOS x86_64, R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The installed size of the package is large (about 75 MB on Linux), almost all of it in libs: the
compiled TMB template, as for other packages that use TMB (e.g. glmmTMB). Depending on the
platform, this may give a NOTE on the installed size.

Examples that take a few seconds (model averaging with bootstrap, the EFSA model averaging) are
wrapped in \donttest{}.
