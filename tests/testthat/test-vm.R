# Tests for the von Mises distribution (circular, support [-pi, pi])
# No quantile function exists, so the round-trip test is skipped (qfun = NULL)

test_that("vm passes standard distribution checks (mu=0, kappa=1)", {
  check_continuous_dist(
    dfun  = dvm,
    pfun  = pvm,
    qfun  = NULL,
    xs    = c(-2, -1, 0, 1, 2),
    lower = -pi, upper = pi,
    mu = 0, kappa = 1
  )
})

test_that("vm passes standard distribution checks (mu=1, kappa=3)", {
  check_continuous_dist(
    dfun  = dvm,
    pfun  = pvm,
    qfun  = NULL,
    xs    = c(-1, 0, 1, 2, 2.5),
    lower = -pi, upper = pi,
    mu = 1, kappa = 3
  )
})

test_that("vm AD gradient has no NaN", {
  check_ad_gradient(dvm,        rvm,        mu = 1, kappa = 2)
})

test_that("pvm under AD matches pvm outside AD", {
  q <- c(-3, -1, 0.99, 1, 1.01, 2, 3.1, 4.5)
  for (kappa in c(0.5, 5, 50)) {
    check_ad_cdf(pvm, dvm, q, mu = 1, kappa = kappa)
  }
  # very concentrated: the derivative in q at q = mu is accurate to about 1e-6
  check_ad_cdf(pvm, dvm, q, mu = 1, kappa = 500, .tol = 1e-5, .grad_tol = 1e-5)
  check_ad_cdf(pvm, dvm, q, mu = 1, kappa = 3, .args = list(lower.tail = FALSE, log.p = TRUE))
  check_ad_cdf(pvm, dvm, q, mu = 1, kappa = 3, .args = list(from = 0.5))
})

test_that("pvm under AD is differentiable in the origin from", {
  q <- c(-3, -1, 0, 1, 2, 3.1)
  F <- RTMB::MakeTape(function(p) pvm(q, 1, 3, from = p), -2)
  expect_equal(F(-2), pvm(q, 1, 3, from = -2), tolerance = 1e-8)
  # moving the origin changes the CDF by minus the density at the origin
  expect_equal(as.vector(F$jacobian(-2)), rep(-dvm(-2, 1, 3), length(q)), tolerance = 1e-8)
})

test_that("pvm under AD stays correct when mu moves angles across mu +- pi", {
  q <- c(-3, -1, 0, 1, 2, 3.1)
  F <- RTMB::MakeTape(function(p) pvm(q, p[1], p[2]), c(0, 2))
  expect_equal(F(c(3, 2)), pvm(q, 3, 2), tolerance = 1e-8)
  expect_equal(F(c(-2.5, 7)), pvm(q, -2.5, 7), tolerance = 1e-8)
})

test_that("vm supports OSA residuals, with the circle cut at -pi", {
  set.seed(1)
  check_osa_cdf(dvm, function(q, mu, kappa) pvm(q, mu, kappa, from = -pi),
                rvm(30, 1, 2), mu = 1, kappa = 2)
})

test_that("vm OSA residuals are valid for mixtures over different mean directions", {
  # the predictive distribution function averages over the components, which is only
  # the distribution function of the mixture with an origin that is the same for all
  set.seed(1)
  y <- c(rvm(20, -2, 3), rvm(20, 2, 1))
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(log(0.3 * exp(dvm(y, mu1, kappa1, log = TRUE)) +
             0.7 * exp(dvm(y, mu2, kappa2, log = TRUE))))
  }
  par <- list(mu1 = -2, kappa1 = 3, mu2 = 2, kappa2 = 1)
  obj <- RTMB::MakeADFun(fn, par, silent = TRUE)
  res <- RTMB::oneStepPredict(obj, method = "cdf", trace = FALSE)
  Fmix <- 0.3 * pvm(y, -2, 3, from = -pi) + 0.7 * pvm(y, 2, 1, from = -pi)
  expect_equal(res$residual, qnorm(Fmix), tolerance = 1e-6)
})
