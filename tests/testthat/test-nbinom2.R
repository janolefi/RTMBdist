# Tests for the negative binomial distribution (mean/size parameterisation)

test_that("nbinom2 passes discrete distribution checks (mu=3, size=2)", {
  check_discrete_dist(
    dfun        = dnbinom2,
    pfun        = pnbinom2,
    xs_int      = c(0, 1, 3, 6, 10),
    sum_support = 0:200,
    mu = 3, size = 2
  )
})

test_that("nbinom2 passes discrete distribution checks (mu=10, size=5)", {
  check_discrete_dist(
    dfun        = dnbinom2,
    pfun        = pnbinom2,
    xs_int      = c(0, 3, 8, 14, 20),
    sum_support = 0:200,
    mu = 10, size = 5
  )
})

test_that("nbinom2 AD gradient has no NaN", {
  check_ad_gradient(dnbinom2,   rnbinom2,   mu = 10, size = 5)
})
