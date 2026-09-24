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

test_that("ztnbinom2 AD gradient has no NaN", {
  check_ad_gradient(dztnbinom2, rztnbinom2, mu = 5, size = 3)
})

test_that("ztnbinom2 supports simulation and OSA residuals", {
  set.seed(1)
  y <- rztnbinom2(40, 3, 2)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dztnbinom2(y, exp(lmu), 2, log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lmu = log(3)), silent = TRUE)
  expect_true(all(obj$simulate()$y >= 1))
  check_osa_cdf(dztnbinom2, pztnbinom2, c(1, 2, 4, 9), mu = 3, size = 2, .discrete = TRUE)
})
