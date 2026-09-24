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

test_that("ztnbinom supports simulation and OSA residuals", {
  set.seed(1)
  y <- rztnbinom(40, 2, 0.4)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dztnbinom(y, 2, RTMB::plogis(lp), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lp = qlogis(0.4)), silent = TRUE)
  expect_true(all(obj$simulate()$y >= 1))
  check_osa_cdf(dztnbinom, pztnbinom, c(1, 2, 4, 9), size = 2, prob = 0.4, .discrete = TRUE)
})
