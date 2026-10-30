# RTMBdist 1.2.0

## Purpose of this submission

This release adds a number of new distributions and copulas, extends support
for one-step-ahead residuals, and fixes several distribution functions and
derivatives that were incorrect. RTMB (>= 2.0) is now required. Random number
streams for `rskewt()`, `rskewnorm()` and `rwrpcauchy()` differ from earlier
versions. NEWS.md has the details.

## Test environments

- GitHub Actions: macOS, Windows and Ubuntu (R release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

RTMBdist has one strong reverse dependency, LaMa, and one package that
suggests it, multiSA.
