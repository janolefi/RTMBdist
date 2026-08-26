# Tests for the Bell distribution.
# The support is unbounded, so sums over the support use a grid that is wide
# enough for the tail to be negligible at the parameter values tested.

test_that("bell passes discrete distribution checks (theta = 0.5)", {
  check_discrete_dist(
    dfun = dbell, pfun = pbell,
    xs_int = c(0, 1, 2, 5, 10), sum_support = 0:200,
    theta = 0.5
  )
})

test_that("bell passes discrete distribution checks (theta = 2)", {
  check_discrete_dist(
    dfun = dbell, pfun = pbell,
    xs_int = c(0, 5, 15, 30, 60), sum_support = 0:400,
    theta = 2
  )
})

test_that("dbell matches the defining formula", {
  # first Bell numbers, https://oeis.org/A000110
  B <- c(1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975)
  x <- 0:10
  for (theta in c(0.3, 1, 2.5)) {
    expect_equal(dbell(x, theta),
                 exp(1 - exp(theta)) * theta^x * B / factorial(x),
                 tolerance = 1e-12,
                 label = paste("pmf at theta =", theta))
  }
})

test_that("dbell has the right mean and variance", {
  for (theta in c(0.2, 1, 2, 3)) {
    xs <- 0:ceiling(theta * exp(theta) + 25 * sqrt(theta * exp(theta) * (1 + theta)) + 60)
    d <- dbell(xs, theta)
    m <- sum(xs * d)
    expect_equal(m, theta * exp(theta), tolerance = 1e-8,
                 label = paste("mean at theta =", theta))
    expect_equal(sum((xs - m)^2 * d), theta * exp(theta) * (1 + theta), tolerance = 1e-7,
                 label = paste("variance at theta =", theta))
  }
})

test_that("bell approaches the Poisson distribution as theta -> 0", {
  # mean theta * exp(theta) -> theta, and the dispersion index 1 + theta -> 1
  expect_equal(dbell(0:6, 1e-7), dpois(0:6, 1e-7), tolerance = 1e-6)
  expect_equal(dbell(0:6, 1e-10), dpois(0:6, 1e-10), tolerance = 1e-9)
})

test_that("log Bell numbers are exact and survive double overflow", {
  B <- c(1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570,
         4213597, 27644437, 190899322)
  expect_equal(lbell(0:14), log(B), tolerance = 1e-14)

  # B_218 = 6.1e306 is the last Bell number representable in double precision,
  # so a pmf built on the Bell numbers directly returns NaN from B_219 on
  expect_true(is.finite(exp(lbell(218))))
  expect_true(is.infinite(exp(lbell(219))))
  expect_true(all(is.finite(lbell(c(219, 300, 1000, 5000)))))

  # the cached triangle and Dobinski's formula are independent computations
  expect_equal(lbell(3000), RTMBdist:::lbell_dobinski(3000), tolerance = 1e-13)
  expect_equal(lbell(20000), RTMBdist:::lbell_dobinski(20000), tolerance = 1e-13)

  # growing the cache must not disturb entries already in it
  before <- lbell(0:50)
  invisible(lbell(4000))
  expect_identical(lbell(0:50), before)
})

test_that("the pmf still normalises for large theta, where counts exceed 218", {
  # theta = 6 gives a mean of 2421: every observation is past the point at
  # which the Bell numbers overflow
  expect_equal(sum(dbell(0:6000, 6)), 1, tolerance = 1e-10)
})

test_that("bell still works where the Bell numbers must come from Dobinski", {
  skip_on_cran()

  # theta = 8 has mean 23848, so the whole grid sits above the cached triangle
  # and every Bell number is obtained from the saddle-point sum instead
  theta <- 8
  mu <- theta * exp(theta)
  expect_gt(mu, bell_tri_max)

  q <- qbell(c(0.01, 0.5, 0.99), theta)
  expect_true(all(is.finite(q)))
  expect_true(all(diff(q) > 0))
  expect_equal(pbell(q[2], theta), 0.5, tolerance = 1e-3)
  # the median of a distribution this concentrated sits next to the mean
  expect_equal(q[2], mu, tolerance = 1e-3)

  # the pmf still integrates to one and has the right mean over a window that
  # covers the mass
  xs <- seq(round(mu - 12 * sqrt(mu * (1 + theta))),
            round(mu + 12 * sqrt(mu * (1 + theta))))
  d <- dbell(xs, theta)
  expect_equal(sum(d), 1, tolerance = 1e-8)
  expect_equal(sum(xs * d), mu, tolerance = 1e-6)
})

test_that("the vectorised Dobinski sum agrees with per-value evaluation", {
  n <- c(20001, 25000, 60000, 250000)
  expect_equal(RTMBdist:::lbell_dobinski(n),
               vapply(n, function(i) RTMBdist:::lbell_dobinski(i), numeric(1)),
               tolerance = 1e-14)
  # and with the exact triangle, where both are available
  expect_equal(RTMBdist:::lbell_dobinski(c(2, 5, 100, 5000, 20000)),
               lbell(c(2, 5, 100, 5000, 20000)), tolerance = 1e-13)
})

test_that("pbell accumulates the pmf", {
  expect_equal(pbell(0:40, 1.4), cumsum(dbell(0:40, 1.4)), tolerance = 1e-12)
  expect_equal(pbell(-1, 1), 0)
  expect_equal(pbell(0, 1), dbell(0, 1))
  expect_equal(pbell(500, 1), 1, tolerance = 1e-12)
  expect_equal(pbell(0:10, 1, log.p = TRUE), log(pbell(0:10, 1)), tolerance = 1e-12)
})

test_that("qbell inverts pbell", {
  for (theta in c(0.3, 1, 2.5)) {
    xs <- 0:qbell(0.9999, theta)
    expect_equal(qbell(pbell(xs, theta), theta), as.numeric(xs),
                 label = paste("round-trip at theta =", theta))
  }
})

test_that("qbell handles the endpoints and the far tail", {
  expect_equal(qbell(0, 1), 0)
  expect_equal(qbell(1, 1), Inf) # as for qpois(), the support is unbounded
  expect_true(is.na(qbell(NA, 1)))
  # the summation grid is grown automatically until it reaches p
  expect_true(is.finite(qbell(1 - 1e-12, 3)))
  expect_gt(qbell(1 - 1e-12, 3), qbell(0.999, 3))
  expect_equal(qbell(0.3, 1), qbell(0.7, 1, lower.tail = FALSE))
  expect_equal(qbell(log(0.3), 1, log.p = TRUE), qbell(0.3, 1))
})

test_that("qbell returns the smallest x with F(x) >= p", {
  theta <- 1.7
  p <- c(0.001, 0.1, 0.5, 0.9, 0.99, 0.99999)
  ref <- vapply(p, function(pp) {
    k <- 0
    while (pbell(k, theta) < pp * (1 - 64 * .Machine$double.eps)) k <- k + 1
    k
  }, numeric(1))
  expect_equal(qbell(p, theta), ref)
})

test_that("bell handles edge cases", {
  # non-integer and negative values carry no mass
  expect_equal(dbell(-1, 1), 0)
  expect_equal(dbell(1.5, 1), 0)
  expect_equal(dbell(c(-1, 1.5), 1, log = TRUE), c(-Inf, -Inf))
  # x = 0 is written without forming 0 * log(theta)
  expect_equal(dbell(0, 1e-300), exp(-expm1(1e-300)))
  expect_false(is.nan(dbell(0, 1e-300)))
})

test_that("bell recycles x and theta", {
  expect_equal(dbell(c(1, 2), c(0.5, 2)), c(dbell(1, 0.5), dbell(2, 2)))
  expect_equal(dbell(2, c(0.5, 2)), c(dbell(2, 0.5), dbell(2, 2)))
  expect_equal(pbell(c(2, 5), c(0.5, 2)), c(pbell(2, 0.5), pbell(5, 2)))
  expect_equal(qbell(c(0.3, 0.7), c(0.5, 2)), c(qbell(0.3, 0.5), qbell(0.7, 2)))
})

test_that("bell rejects invalid arguments", {
  expect_error(dbell(1, 0), "theta")
  expect_error(dbell(1, -1), "theta")
  expect_error(dbell(1, Inf), "theta")
  expect_error(qbell(1.2, 1), "p must be in")
  expect_error(qbell(-0.1, 1), "p must be in")
})

test_that("rbell matches the pmf", {
  set.seed(42)
  n <- 3e5
  r <- rbell(n, 1.3)
  expect_true(all(r >= 0 & r == floor(r)))
  emp <- as.numeric(table(factor(r, levels = 0:14))) / n
  expect_equal(emp, dbell(0:14, 1.3), tolerance = 0.005)
  # the compound Poisson representation must reproduce the moments
  expect_equal(mean(r), 1.3 * exp(1.3), tolerance = 0.05)
  expect_equal(var(r), 1.3 * exp(1.3) * 2.3, tolerance = 0.5)
})

test_that("rbell recycles theta and handles n given as a vector", {
  set.seed(1)
  r <- rbell(4, c(0.01, 5))
  expect_length(r, 4)
  expect_true(all(r[c(1, 3)] < r[c(2, 4)])) # theta = 0.01 vs theta = 5
  expect_length(rbell(c(9, 9, 9), 1), 3)
})

# ---------------------------------------------------------------------------
# AD behaviour
# ---------------------------------------------------------------------------

test_that("bell AD gradient has no NaN", {
  check_ad_gradient(dbell, rbell, theta = 1.2)
})

test_that("dbell AD gradient matches finite differences and the tape replays", {
  set.seed(42)
  y <- rbell(40, 1.2)
  nll <- function(par) -sum(dbell(y, exp(par), log = TRUE))
  environment(nll) <- environment()
  F <- RTMB::MakeTape(nll, 0.2)

  h <- 1e-6
  for (pt in c(-1.5, 0.2, 1.1, 2)) {
    expect_equal(F(pt), nll(pt), tolerance = 1e-9,
                 label = paste("tape replay at", pt))
    expect_equal(as.numeric(F$jacobian(pt)), (nll(pt + h) - nll(pt - h)) / (2 * h),
                 tolerance = 1e-4, label = paste("gradient at", pt))
  }
})

test_that("pbell is AD-able in theta", {
  q <- c(0, 2, 5)
  f <- function(theta) sum(pbell(q, theta))
  environment(f) <- environment()
  F <- RTMB::MakeTape(f, 1.2)
  expect_equal(F(1.2), f(1.2), tolerance = 1e-12)
  h <- 1e-6
  expect_equal(as.numeric(F$jacobian(1.2)), (f(1.2 + h) - f(1.2 - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("dbell refuses to differentiate with respect to x", {
  # log B_x is only defined on the integers; the smooth extension given by
  # Dobinski's formula disagrees with B_0 and is not integrable below zero,
  # so it would not give valid oneStepGeneric residuals either
  expect_error(RTMB::MakeTape(function(x) dbell(x, 1), 3), "not differentiable")
  expect_error(RTMB::MakeTape(function(q) pbell(q, 1), 3), "not differentiable")
})

test_that("the MLE reproduces the exponential family moment identity", {
  # Bell is a one-parameter exponential family with sufficient statistic x,
  # so the fitted mean must equal the sample mean exactly
  set.seed(42)
  y <- rbell(200, 1.2)
  obj <- RTMB::MakeADFun(function(p) {
    RTMB::getAll(p)
    -sum(dbell(y, exp(logtheta), log = TRUE))
  }, list(logtheta = 0), silent = TRUE)
  o <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(o$convergence, 0)
  theta <- exp(o$par[[1]])
  expect_equal(theta * exp(theta), mean(y), tolerance = 1e-6)
})

test_that("bell supports simulation and OSA residuals via the CDF", {
  skip_on_cran()

  set.seed(123)
  y <- rbell(80, 1.2)
  dat <- list(y = y)

  fn <- function(par) {
    RTMB::getAll(par, dat)
    y <- RTMB::OBS(y)
    -sum(dbell(y, exp(logtheta), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(logtheta = 0), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(opt$convergence, 0)

  set.seed(1)
  sim <- obj$simulate()
  expect_true(all(sim$y >= 0 & sim$y == floor(sim$y)))
  expect_length(sim$y, 80)

  res <- RTMB::oneStepPredict(obj, method = "cdf", discrete = TRUE, trace = FALSE)
  expect_false(any(is.nan(res$residual)))
  expect_lt(abs(mean(res$residual)), 0.4)
  expect_lt(abs(sd(res$residual) - 1), 0.4)
})
