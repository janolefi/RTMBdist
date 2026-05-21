# Tests for the zero-truncated negative binomial distribution (size/prob parameterisation)
# Support starts at 1.

test_that("ztnbinom passes discrete distribution checks (size=3, prob=0.4)", {
  check_discrete_dist(
    dfun        = dztnbinom,
    pfun        = pztnbinom,
    xs_int      = c(1, 2, 4, 7, 11),
    sum_support = 1:200,
    size = 3, prob = 0.4
  )
})

test_that("ztnbinom passes discrete distribution checks (size=5, prob=0.6)", {
  check_discrete_dist(
    dfun        = dztnbinom,
    pfun        = pztnbinom,
    xs_int      = c(1, 2, 4, 6, 9),
    sum_support = 1:100,
    size = 5, prob = 0.6
  )
})

test_that("ztnbinom AD gradient has no NaN", {
  check_ad_gradient(dztnbinom,  rztnbinom,  size = 5, prob = 0.4)
})
