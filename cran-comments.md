## Submission

This is the first submission of drcTMB.

## Test environments

* Local: Ubuntu 26.04.1 LTS, R 4.6.2

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The installed size of the package is large (about 75 MB on Linux), almost all of it in libs: the
compiled TMB template, as for other packages that use TMB (e.g. glmmTMB). Depending on the
platform, this may give a NOTE on the installed size.

Examples that take a few seconds (model averaging with bootstrap, the EFSA model averaging) are
wrapped in \donttest{}.
