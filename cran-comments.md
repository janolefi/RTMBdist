# RTMBdist 1.2.1

## Purpose of this submission

I am aware that this update comes soon after version 1.2.0 (published
2026-10-31), and apologise for that. The update is needed for a paper I am
about to submit: the paper relies on a new feature of `LaMa::MCreport()`, a thin
wrapper around `RTMBdist::mcreport()`, and that feature needs the changes in this
release, so the update to LaMa can only follow once this version is on CRAN.

`mcreport()` gains a `par` argument for the mean of the samples and now uses a
supplied precision matrix `Q` for models without random effects, including one
with additional trailing hyperparameter coordinates. Defaults are unchanged.
NEWS.md has the details.

## Test environments

- GitHub Actions: macOS and Windows (R release), Ubuntu (R devel, release and oldrel-1)
- Local: macOS

## R CMD check results

0 errors | 0 warnings | 1 note

- Days since last update: see above.

## Reverse dependencies

RTMBdist has one strong reverse dependency, LaMa, whose test suite passes with
this version, and one package that suggests it, multiSA, which does not use
`mcreport()`.
