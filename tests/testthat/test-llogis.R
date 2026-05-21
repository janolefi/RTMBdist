# Tests for the log-logistic distribution
# Support: x > 0; parameters alpha > 0 (scale/median), beta > 0 (shape)
# Note: for beta < 1 the density is unbounded at x = 0, so normalisation
# tests use beta >= 1 to avoid non-finite values at the lower integration boundary.
#

test_that("llogis passes check_continuous_dist (alpha=1, beta=2)", {
  check_continuous_dist(
    dllogis, pllogis, qllogis,
    xs    = c(0.25, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    alpha = 1, beta = 2
  )
})

test_that("llogis passes check_continuous_dist (alpha=3, beta=1.5)", {
  check_continuous_dist(
    dllogis, pllogis, qllogis,
    xs    = c(0.5, 1, 3, 6, 12),
    lower = 0, upper = Inf,
    alpha = 3, beta = 1.5
  )
})

test_that("llogis AD gradient has no NaN", {
  check_ad_gradient(dllogis, rllogis, alpha = 1, beta = 2)
})
