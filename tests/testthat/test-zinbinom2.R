# Tests for the zero-inflated negative binomial distribution (mean/size parameterisation)

test_that("zinbinom2 passes discrete distribution checks (mu=4, size=2, zeroprob=0.2)", {
  check_discrete_dist(
    dfun        = dzinbinom2,
    pfun        = pzinbinom2,
    xs_int      = c(0, 1, 3, 6, 10),
    sum_support = 0:200,
    mu = 4, size = 2, zeroprob = 0.2
  )
})

test_that("zinbinom2 passes discrete distribution checks (mu=8, size=5, zeroprob=0.15)", {
  check_discrete_dist(
    dfun        = dzinbinom2,
    pfun        = pzinbinom2,
    xs_int      = c(0, 2, 6, 10, 15),
    sum_support = 0:200,
    mu = 8, size = 5, zeroprob = 0.15
  )
})
