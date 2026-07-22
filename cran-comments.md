# RTMBdist 1.0.5

## Reason for submission

RTMBdist was archived on 2026-07-15 following an R-devel regression in S4
method dispatch for `apply()`. This was a change in R itself, not in RTMBdist,
and R-core has since reverted it. The remaining macOS build issue in the
dependency RTMB has also been resolved by its maintainer, and RTMB now checks
cleanly on all CRAN platforms. RTMBdist passes checks on all platforms again.

This version also adds minor stability improvements to `rgmrf()`.

## Test environments

- local: macOS, R 4.x (release)
- macOS builder (r-devel)
- win-builder (r-devel and r-release)

## R CMD check results

0 errors | 0 warnings | 0 notes
