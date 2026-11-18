## Submission

This is the first submission of drcTMB.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.6.2
* win-builder: Windows, R-devel (2026-11-17 r90926 ucrt)
* R-hub (GitHub Actions): Windows, R-devel; macOS arm64, R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

  The possibly misspelled words in DESCRIPTION are names and terms spelled correctly: EFSA
  (European Food Safety Authority), Gompertz, Michaelis-Menten, Kristensen et al., hormesis.

The installed size of the package is large (about 75 MB on Linux), almost all of it in libs: the
compiled TMB template, as for other packages that use TMB (e.g. glmmTMB). Depending on the
platform, this may give a NOTE on the installed size.

Examples that take a few seconds (model averaging with bootstrap, the EFSA model averaging) are
wrapped in \donttest{}.
