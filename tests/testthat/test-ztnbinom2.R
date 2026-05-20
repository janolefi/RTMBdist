# Tests for the zero-truncated negative binomial distribution (mean/size parameterisation)
# Support starts at 1.

test_that("ztnbinom2 passes discrete distribution checks (mu=4, size=2)", {
  check_discrete_dist(
    dfun        = dztnbinom2,
    pfun        = pztnbinom2,
    xs_int      = c(1, 2, 4, 7, 11),
    sum_support = 1:200,
    mu = 4, size = 2
  )
})

test_that("ztnbinom2 passes discrete distribution checks (mu=8, size=5)", {
  check_discrete_dist(
    dfun        = dztnbinom2,
    pfun        = pztnbinom2,
    xs_int      = c(1, 3, 6, 10, 15),
    sum_support = 1:200,
    mu = 8, size = 5
  )
})
