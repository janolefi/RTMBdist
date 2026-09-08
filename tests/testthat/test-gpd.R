# Tests for the generalised Pareto distribution

test_that("gpd passes standard distribution checks (xi > 0)", {
  check_continuous_dist(
    dfun  = dgpd,
    pfun  = pgpd,
    qfun  = qgpd,
    xs    = c(0.1, 0.6, 1.5, 4, 12),
    lower = 0, upper = Inf,
    mu = 0, sigma = 1, xi = 0.3
  )
})

test_that("gpd passes standard distribution checks (xi < 0, bounded above)", {
  check_continuous_dist(
    dfun  = dgpd,
    pfun  = pgpd,
    qfun  = qgpd,
    xs    = c(2.1, 3, 4.5, 6, 6.9),
    lower = 2, upper = 7, # support is 2 <= x <= mu - sigma/xi = 7
    mu = 2, sigma = 2, xi = -0.4
  )
})

test_that("gpd passes standard distribution checks (xi = 0)", {
  check_continuous_dist(
    dfun  = dgpd,
    pfun  = pgpd,
    qfun  = qgpd,
    xs    = c(0.1, 0.5, 1.2, 3, 8),
    lower = 0, upper = Inf,
    mu = 0, sigma = 1, xi = 0
  )
})

test_that("gpd AD gradient has no NaN", {
  check_ad_gradient(dgpd, rgpd, mu = 0, sigma = 1, xi = 0.3)
  check_ad_gradient(dgpd, rgpd, mu = 0, sigma = 2, xi = -0.3)
  check_ad_gradient(dgpd, rgpd, mu = 0, sigma = 1, xi = 0)
})

test_that("gpd reduces to the exponential distribution at xi = 0", {
  xs <- c(0, 0.5, 1, 3, 9)
  expect_equal(dgpd(xs, 0, 2, 0), stats::dexp(xs, 1 / 2))
  expect_equal(pgpd(xs, 0, 2, 0), stats::pexp(xs, 1 / 2))
})

test_that("gpd reduces to the Pareto distribution for a matching threshold", {
  # a GPD with mu = sigma / xi is Pareto with shape 1 / xi on x > mu
  xs <- c(1.01, 1.5, 3, 10, 100)
  expect_equal(dgpd(xs, mu = 1, sigma = 1 / 3, xi = 1 / 3), dpareto(xs, mu = 3))
  expect_equal(pgpd(xs, mu = 1, sigma = 1 / 3, xi = 1 / 3), ppareto(xs, mu = 3))
})

test_that("gpd is zero below the threshold and the end points are right", {
  expect_equal(dgpd(c(-3, -1e-8), 0, 2, 0.3), c(0, 0))
  expect_equal(pgpd(c(-3, -1e-8), 0, 2, 0.3), c(0, 0))

  # the threshold itself belongs to the support, as it does for stats::dexp
  expect_equal(dgpd(0, 0, 2, 0.3), 1 / 2)
  expect_equal(pgpd(0, 0, 2, 0.3), 0)

  expect_equal(qgpd(0, 1, 2, 0.3), 1)
  expect_equal(qgpd(1, 1, 2, 0.3), Inf)

  # xi < 0 bounds the support above at mu - sigma / xi
  expect_equal(dgpd(c(5, 6), 0, 2, -0.4), c(0, 0))
  expect_equal(pgpd(c(5, 6), 0, 2, -0.4), c(1, 1))
  expect_equal(qgpd(1, 0, 2, -0.4), 5)
})

test_that("the gpd gradient with respect to xi is exact at xi = 0", {
  set.seed(1)
  x <- rgpd(50, 0, 1, 0.2)
  nll <- function(p) -sum(dgpd(x, 0, exp(p[1]), p[2], log = TRUE))
  F <- RTMB::MakeTape(nll, c(0, 0.2))

  p <- c(0, 0)
  g <- as.vector(F$jacobian(p))
  num <- sapply(seq_along(p), function(i) {
    h <- 1e-5; up <- lo <- p; up[i] <- up[i] + h; lo[i] <- lo[i] - h
    (nll(up) - nll(lo)) / (2 * h)
  })
  expect_false(any(is.nan(g)))
  expect_equal(g, num, tolerance = 1e-6)
  expect_true(abs(g[2]) > 1)
})

test_that("gpd gradients stay finite when observations fall outside the support", {
  set.seed(2)
  x <- rgpd(50, 0, 1, -0.4)
  nll <- function(p) -sum(dgpd(x, p[1], exp(p[2]), p[3], log = TRUE))
  F <- RTMB::MakeTape(nll, c(0, 0, -0.4))

  expect_equal(F(c(0, 0, -1.5)), Inf)   # upper end point pulled below the data
  expect_equal(F(c(0.5, 0, -0.4)), Inf) # threshold pushed above the data
  expect_false(any(is.nan(as.vector(F$jacobian(c(0, 0, -1.5))))))
  expect_false(any(is.nan(as.vector(F$jacobian(c(0.5, 0, -0.4))))))
})

test_that("qgpd honours lower.tail and log.p and recycles its arguments", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qgpd(p, 1, 2, 0.3), qgpd(1 - p, 1, 2, 0.3, lower.tail = FALSE))
  expect_equal(qgpd(p, 1, 2, 0.3), qgpd(log(p), 1, 2, 0.3, log.p = TRUE))
  expect_length(qgpd(0.5, mu = 1:4), 4)
  expect_length(pgpd(1:6, 0, c(1, 2), 0.1), 6)
})

test_that("dgpd and pgpd reject a non-positive scale", {
  expect_error(dgpd(1, 0, -1, 0), "sigma")
  expect_error(pgpd(1, 0, 0, 0), "sigma")
})
