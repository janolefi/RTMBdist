# Temporary workaround for RTMB 2.0: the derivatives of pbeta() and qbeta() in the shapes
# are wrong (often zero) when a shape is an integer or half-integer, and their second
# derivatives are inaccurate within about 1e-4 of such values. This is fixed in the TMB
# development version. Once the fix is on CRAN, require that RTMB version in DESCRIPTION,
# replace every beta_ad(pbeta, ...) by pbeta(...) and beta_ad(qbeta, ...) by
# qbeta(...), and delete this file and tests/testthat/test-pbeta-workaround.R.
#
# f is pbeta or qbeta. In AD context, the shapes are moved by -+d and -+2d, with d about
# 0.0014 (relative for shapes below 1, which keeps them positive), and the four results are
# combined by Richardson extrapolation, which cancels the error up to O(d^4). Values then
# agree with stats::pbeta() to about 1e-11, and first and second derivatives are correct.
# d is irrational so that round shapes never move onto an integer or half-integer.
beta_ad <- function(f, x, shape1, shape2) {
  if (!ad_context()) return(f(x, shape1, shape2))

  e <- sqrt(2) * 1e-3
  d1 <- e * shape1 / (1 + shape1)
  d2 <- e * shape2 / (1 + shape2)
  m <- function(k) (f(x, shape1 + k * d1, shape2 + k * d2) + f(x, shape1 - k * d1, shape2 - k * d2)) / 2

  (4 * m(1) - m(2)) / 3
}
