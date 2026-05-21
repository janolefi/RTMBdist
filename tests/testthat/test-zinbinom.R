# Tests for the zero-inflated negative binomial distribution (size/prob parameterisation)

test_that("zinbinom passes discrete distribution checks (size=3, prob=0.4, zeroprob=0.2)", {
  check_discrete_dist(
    dfun        = dzinbinom,
    pfun        = pzinbinom,
    xs_int      = c(0, 1, 3, 6, 10),
    sum_support = 0:200,
    size = 3, prob = 0.4, zeroprob = 0.2
  )
})

test_that("zinbinom passes discrete distribution checks (size=5, prob=0.6, zeroprob=0.3)", {
  check_discrete_dist(
    dfun        = dzinbinom,
    pfun        = pzinbinom,
    xs_int      = c(0, 1, 3, 5, 8),
    sum_support = 0:100,
    size = 5, prob = 0.6, zeroprob = 0.3
  )
})

test_that("zinbinom AD gradient has no NaN", {
  check_ad_gradient(dzinbinom,  rzinbinom,  size = 5, prob = 0.4, zeroprob = 0.2)
})
