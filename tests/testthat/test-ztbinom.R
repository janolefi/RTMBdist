# Tests for the zero-truncated binomial distribution
# Support starts at 1.

test_that("ztbinom passes discrete distribution checks (size=10, prob=0.4)", {
  check_discrete_dist(
    dfun        = dztbinom,
    pfun        = pztbinom,
    xs_int      = c(1, 3, 5, 7, 10),
    sum_support = 1:10,
    size = 10, prob = 0.4
  )
})

test_that("ztbinom passes discrete distribution checks (size=20, prob=0.3)", {
  check_discrete_dist(
    dfun        = dztbinom,
    pfun        = pztbinom,
    xs_int      = c(1, 4, 8, 12, 18),
    sum_support = 1:20,
    size = 20, prob = 0.3
  )
})

test_that("ztbinom AD gradient has no NaN", {
  check_ad_gradient(dztbinom,   rztbinom,   size = 10, prob = 0.4)
})

test_that("ztbinom supports simulation and OSA residuals", {
  set.seed(1)
  y <- rztbinom(40, 5, 0.3)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dztbinom(y, 5, RTMB::plogis(lp), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lp = qlogis(0.3)), silent = TRUE)
  sim <- obj$simulate()$y
  expect_true(all(sim >= 1 & sim <= 5))
  check_osa_cdf(dztbinom, pztbinom, c(1, 2, 3, 5), size = 5, prob = 0.3, .discrete = TRUE)
})
