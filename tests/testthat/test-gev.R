# Tests for the generalised extreme value distribution

test_that("gev passes standard distribution checks (xi > 0, heavy tail)", {
  check_continuous_dist(
    dfun  = dgev,
    pfun  = pgev,
    qfun  = qgev,
    xs    = c(-2, -0.5, 0.3, 1.5, 6),
    lower = -5, upper = Inf, # support is x > mu - sigma/xi = -5
    mu = 0, sigma = 1, xi = 0.2
  )
})

test_that("gev passes standard distribution checks (xi < 0, bounded above)", {
  check_continuous_dist(
    dfun  = dgev,
    pfun  = pgev,
    qfun  = qgev,
    xs    = c(-1, 0.5, 2, 3.5, 4.5),
    lower = -Inf, upper = 6, # support is x < mu - sigma/xi = 6
    mu = 1, sigma = 2, xi = -0.4
  )
})

test_that("gev passes standard distribution checks (xi = 0)", {
  check_continuous_dist(
    dfun  = dgev,
    pfun  = pgev,
    qfun  = qgev,
    xs    = c(-1.5, -0.5, 0.5, 2, 5),
    mu = 0, sigma = 1, xi = 0
  )
})

test_that("gev AD gradient has no NaN", {
  check_ad_gradient(dgev, rgev, mu = 0, sigma = 1, xi = 0.3)
  check_ad_gradient(dgev, rgev, mu = 2, sigma = 3, xi = -0.3)
  check_ad_gradient(dgev, rgev, mu = 0, sigma = 1, xi = 0)
})

test_that("gev reduces to the Gumbel distribution at xi = 0", {
  xs <- c(-3, -1, 0, 1, 4, 10)
  expect_equal(dgev(xs, 1, 2, 0), dgumbel(xs, 1, 2))
  expect_equal(pgev(xs, 1, 2, 0), pgumbel(xs, 1, 2))
  expect_equal(qgev(c(0.1, 0.5, 0.9), 1, 2, 0), qgumbel(c(0.1, 0.5, 0.9), 1, 2))
})

test_that("gev is zero outside its support and the end points are right", {
  # xi > 0: bounded below by mu - sigma/xi
  expect_equal(dgev(c(-6, -5), 0, 1, 0.2), c(0, 0))
  expect_equal(pgev(c(-6, -5), 0, 1, 0.2), c(0, 0))
  expect_equal(dgev(-5 + 1e-8, 0, 1, 0.2), 0) # density vanishes at the end point
  expect_equal(qgev(0, 0, 1, 0.2), -5)
  expect_equal(qgev(1, 0, 1, 0.2), Inf)

  # xi < 0: bounded above by mu - sigma/xi
  expect_equal(dgev(c(5, 6), 0, 1, -0.2), c(0, 0))
  expect_equal(pgev(c(5, 6), 0, 1, -0.2), c(1, 1))
  expect_equal(qgev(0, 0, 1, -0.2), -Inf)
  expect_equal(qgev(1, 0, 1, -0.2), 5)
})

test_that("the gev gradient with respect to xi is exact at xi = 0", {
  # the branch-free log(1 + u) / u keeps the shape derivative correct at the
  # value people start an optimiser at; a branch on iszero(xi) would give 0 here
  set.seed(1)
  x <- rgev(50, 0, 1, 0.2)
  nll <- function(p) -sum(dgev(x, p[1], exp(p[2]), p[3], log = TRUE))
  F <- RTMB::MakeTape(nll, c(0, 0, 0.2))

  p <- c(0, 0, 0)
  g <- as.vector(F$jacobian(p))
  num <- sapply(seq_along(p), function(i) {
    h <- 1e-5; up <- lo <- p; up[i] <- up[i] + h; lo[i] <- lo[i] - h
    (nll(up) - nll(lo)) / (2 * h)
  })
  expect_false(any(is.nan(g)))
  expect_equal(g, num, tolerance = 1e-6)
  expect_true(abs(g[3]) > 1) # and in particular it is not zero
})

test_that("gev gradients stay finite when observations fall outside the support", {
  set.seed(2)
  x <- rgev(50, 0, 1, -0.3) # bounded above by 1 / 0.3
  nll <- function(p) -sum(dgev(x, p[1], exp(p[2]), p[3], log = TRUE))
  F <- RTMB::MakeTape(nll, c(0, 0, -0.3))

  # a shape this small puts the upper end point below the largest observation
  expect_equal(F(c(0, 0, -0.9)), Inf)
  expect_false(any(is.nan(as.vector(F$jacobian(c(0, 0, -0.9))))))
  expect_false(any(is.nan(as.vector(F$jacobian(c(-5, 0, 0.6))))))
})

test_that("qgev honours lower.tail and log.p and recycles its arguments", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qgev(p, 1, 2, 0.3), qgev(1 - p, 1, 2, 0.3, lower.tail = FALSE))
  expect_equal(qgev(p, 1, 2, 0.3), qgev(log(p), 1, 2, 0.3, log.p = TRUE))
  expect_length(qgev(0.5, mu = 1:4), 4)
  expect_length(dgev(1:6, mu = c(0, 1), sigma = 1, xi = c(0, 0.2, -0.2)), 6)
})

test_that("dgev and pgev reject a non-positive scale", {
  expect_error(dgev(1, 0, -1, 0), "sigma")
  expect_error(pgev(1, 0, 0, 0), "sigma")
})
