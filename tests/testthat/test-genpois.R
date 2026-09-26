# Tests for the generalised Poisson distribution

test_that("genpois passes discrete distribution checks (lambda=3, phi=0.3)", {
  check_discrete_dist(
    dfun        = dgenpois,
    pfun        = pgenpois,
    xs_int      = c(0, 1, 3, 5, 8),
    sum_support = 0:100,
    lambda = 3, phi = 0.3
  )
})

test_that("genpois passes discrete distribution checks (lambda=5, phi=0.5)", {
  check_discrete_dist(
    dfun        = dgenpois,
    pfun        = pgenpois,
    xs_int      = c(0, 2, 4, 7, 11),
    sum_support = 0:200,
    lambda = 5, phi = 0.5
  )
})

test_that("genpois AD gradient has no NaN", {
  check_ad_gradient(dgenpois,   rgenpois,   lambda = 5, phi = 0.5)
})

test_that("qgenpois honours lower.tail and log.p", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qgenpois(p, lambda = 5, phi = 0.5),
               qgenpois(1 - p, lambda = 5, phi = 0.5, lower.tail = FALSE))
  expect_equal(qgenpois(p, lambda = 5, phi = 0.5),
               qgenpois(log(p), lambda = 5, phi = 0.5, log.p = TRUE))
})

test_that("pgenpois recycles lambda and phi like gamlss.dist did", {
  target <- pgenpois(0:5, 3, 0.2)
  expect_equal(pgenpois(0:5, rep(3, 6), rep(0.2, 6)), target)
  expect_equal(pgenpois(0:5, c(3, 3), c(0.2, 0.2)), target)
})

test_that("pgenpois is AD-compatible in its parameters, so OSA works with estimated parameters", {
  q <- c(0, 1, 3, 6)
  check_ad_cdf(pgenpois, dgenpois, q, lambda = 3, phi = 0.2, .dq = FALSE)
  set.seed(1)
  y <- rgenpois(40, 3, 0.2)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dgenpois(y, exp(llambda), phi, log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(llambda = log(3), phi = 0.2), silent = TRUE)
  res <- RTMB::oneStepPredict(obj, method = "cdf", discrete = TRUE, trace = FALSE)
  expect_equal(res$Fx, pgenpois(y, 3, 0.2), tolerance = 1e-8)
})
