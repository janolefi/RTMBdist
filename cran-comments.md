# RTMBdist 1.1.0

## Purpose of this submission

This release removes the dependency on 'gamlss.dist', which is scheduled for
archival on 2026-09-27.

RTMBdist previously imported seventeen quantile and random generation functions
from 'gamlss.dist'. All of them are now implemented within RTMBdist itself, and
'gamlss.dist' has been dropped from Imports. The archival therefore no longer
affects this package, nor its strong reverse dependency 'LaMa'.

The ported functions were checked against 'gamlss.dist' 6.1-1 over a grid of
roughly 1750 parameter combinations before and after the change. All results are
unchanged except where noted in NEWS.md, and random number generation is
identical stream-for-stream under a common seed.

## License change

The license has changed from MIT to GPL-2 | GPL-3.

Several of the distributions in this package (BCCG, BCT, BCPE, PE, PE2, GG, GPO,
exGAUS and PARETO) are derived from the corresponding families in 'gamlss.dist',
which is released under GPL-2 | GPL-3. Previously most of that code was reached
through Imports; it is now contained in RTMBdist itself. The package is
therefore released under the same terms as the code it is derived from, matching
the 'gamlss.dist' license field exactly. The derived files carry comments
identifying the 'gamlss.dist' source file each was taken from.

## Test environments

- local: macOS 15 (aarch64-apple-darwin20), R 4.5.3
- win-builder: devel and release
- macOS builder: R release

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

RTMBdist has one strong reverse dependency, 'LaMa', and one package that
suggests it, 'multiSA'. Both were checked against this version and are
unaffected.

## Other changes

Besides removing the dependency, this release fixes a number of argument
handling bugs that were inherited together with the imported functions, most
importantly argument recycling when a single quantile is evaluated against
vectorised parameters. It also adds the Bell and Conway-Maxwell-binomial
distributions and an AD-compatible Lambert W function. NEWS.md has the details.
