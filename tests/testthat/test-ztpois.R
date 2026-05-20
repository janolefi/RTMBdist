# Tests for the zero-truncated Poisson distribution
# Support starts at 1, so sum_support begins at 1.

test_that("ztpois passes discrete distribution checks (lambda=3)", {
  check_discrete_dist(
    dfun        = dztpois,
    pfun        = pztpois,
    xs_int      = c(1, 2, 3, 5, 8),
    sum_support = 1:100,
    lambda = 3
  )
})

test_that("ztpois passes discrete distribution checks (lambda=8)", {
  check_discrete_dist(
    dfun        = dztpois,
    pfun        = pztpois,
    xs_int      = c(1, 3, 6, 10, 15),
    sum_support = 1:100,
    lambda = 8
  )
})
