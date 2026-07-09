# RTMBdist 1.0.5

## Reason for submission

This patch release fixes the ERROR reported in the CRAN check results
(additional issue, M1mac): re-building the vignette 'Examples.Rmd' failed on
r-devel (aarch64-apple-darwin) with "invalid 'type' (complex) of argument".

Cause: the internal helpers pmin.ad()/pmax.ad() used
apply(cbind(x, y), 1, min/max), which relied on class-preserving method
dispatch for RTMB's AD types. They have been replaced by dispatch-safe
equivalents, removing the failure mode entirely.

## Note on r-release-macos-x86_64

The ERROR on this flavor ("Package required but not available: 'RTMB'") is
caused by the dependency RTMB not being available there. This is outside the
control of this package.

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked on:

* local macOS (aarch64), R release
* GitHub Actions: macOS (R release), Windows (R release),
  Ubuntu (R devel, release, oldrel-1)
