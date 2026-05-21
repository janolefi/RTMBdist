# Tests for the zero-inflated Poisson distribution
# Note: unlike zero-inflated *continuous* distributions, d(0) already encodes
# the total probability mass at 0, so sum(d(0:N)) = 1 holds directly.

test_that("zipois passes discrete distribution checks (lambda=3, zeroprob=0.3)", {
  check_discrete_dist(
    dfun        = dzipois,
    pfun        = pzipois,
    xs_int      = c(0, 1, 3, 5, 8),
    sum_support = 0:100,
    lambda = 3, zeroprob = 0.3
  )
})

test_that("zipois passes discrete distribution checks (lambda=5, zeroprob=0.2)", {
  check_discrete_dist(
    dfun        = dzipois,
    pfun        = pzipois,
    xs_int      = c(0, 2, 4, 7, 10),
    sum_support = 0:100,
    lambda = 5, zeroprob = 0.2
  )
})

test_that("zipois AD gradient has no NaN", {
  check_ad_gradient(dzipois,    rzipois,    lambda = 3, zeroprob = 0.2)
})
