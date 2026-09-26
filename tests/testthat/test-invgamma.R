# Tests for the inverse gamma distribution
# If X ~ Gamma(shape, scale), then 1/X ~ InvGamma(shape, rate=1/scale)

test_that("invgamma passes standard distribution checks (shape=2, rate=1)", {
  check_continuous_dist(
    dfun  = dinvgamma,
    pfun  = pinvgamma,
    qfun  = qinvgamma,
    xs    = c(0.2, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    shape = 2, rate = 1
  )
})

test_that("invgamma passes standard distribution checks (shape=3, rate=2)", {
  check_continuous_dist(
    dfun  = dinvgamma,
    pfun  = pinvgamma,
    qfun  = qinvgamma,
    xs    = c(0.1, 0.3, 0.6, 1, 2),
    lower = 0, upper = Inf,
    shape = 3, rate = 2
  )
})

test_that("invgamma AD gradient has no NaN", {
  check_ad_gradient(dinvgamma,  rinvgamma,  shape = 3, rate = 2)
})

test_that("invgamma functions accept scale instead of rate", {
  x <- c(0.3, 1, 2)
  expect_equal(dinvgamma(x, 3, scale = 0.5), dinvgamma(x, 3, rate = 2))
  expect_equal(pinvgamma(x, 3, scale = 0.5), pinvgamma(x, 3, rate = 2))
  expect_equal(qinvgamma(c(0.1, 0.5, 0.9), 3, scale = 0.5), qinvgamma(c(0.1, 0.5, 0.9), 3, rate = 2))
  set.seed(1); a <- rinvgamma(5, 3, scale = 0.5)
  set.seed(1); b <- rinvgamma(5, 3, rate = 2)
  expect_equal(a, b)
  expect_error(pinvgamma(1, 3, rate = 2, scale = 0.5), "not both")
  expect_error(dinvgamma(1, 3), "must be given")
})

test_that("invgamma works with scale as an AD parameter, for simulation and OSA as well", {
  check_ad_cdf(pinvgamma, dinvgamma, c(0.3, 0.8, 2), shape = 3, scale = 0.5)
  set.seed(1)
  y <- rinvgamma(30, 3, scale = 0.5)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dinvgamma(y, exp(lshape), scale = exp(lscale), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lshape = log(3), lscale = log(0.5)), silent = TRUE)
  set.seed(2)
  expect_true(all(obj$simulate()$y > 0))
  res <- RTMB::oneStepPredict(obj, method = "cdf", trace = FALSE)
  expect_equal(res$residual, qnorm(pinvgamma(y, 3, scale = 0.5)), tolerance = 1e-6)
})

test_that("pinvgamma has correct Hessians and finite third derivatives, also at q = scale", {
  q <- c(0.2, 0.8, 2, 5, 10)
  check_ad_cdf_hessian(pinvgamma, q, shape = 3, scale = 2)
  F3 <- RTMB::MakeTape(function(p) sum(log(pinvgamma(q, p[1], scale = p[2]))), c(3, 2))$jacfun()$jacfun()
  expect_true(all(is.finite(F3$jacobian(c(3, 2)))))
})
