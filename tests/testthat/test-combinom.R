# Tests for the Conway-Maxwell-binomial distribution.
# Support is 0:size, so the pmf sums exactly over a finite grid and the
# q(p(x)) round-trip applies on the integers.

test_that("combinom passes discrete distribution checks (size=12, prob=0.35, nu=1.7)", {
  check_discrete_dist(
    dfun        = dcombinom,
    pfun        = pcombinom,
    xs_int      = c(0, 2, 5, 8, 12),
    sum_support = 0:12,
    size = 12, prob = 0.35, nu = 1.7
  )
})

test_that("combinom passes discrete distribution checks (size=20, prob=0.6, nu=0.5)", {
  check_discrete_dist(
    dfun        = dcombinom,
    pfun        = pcombinom,
    xs_int      = c(0, 5, 10, 15, 20),
    sum_support = 0:20,
    size = 20, prob = 0.6, nu = 0.5
  )
})

test_that("combinom passes discrete distribution checks for negative nu (size=15, nu=-0.5)", {
  # nu < 0 is a valid, strongly over-dispersed regime: the support is finite,
  # so the normalising constant converges for every real nu.
  check_discrete_dist(
    dfun        = dcombinom,
    pfun        = pcombinom,
    xs_int      = c(0, 3, 7, 11, 15),
    sum_support = 0:15,
    size = 15, prob = 0.4, nu = -0.5
  )
})

test_that("combinom reduces to the binomial distribution at nu = 1", {
  expect_equal(dcombinom(0:12, 12, 0.3, 1), dbinom(0:12, 12, 0.3), tolerance = 1e-12)
  expect_equal(pcombinom(0:12, 12, 0.3, 1), pbinom(0:12, 12, 0.3), tolerance = 1e-12)
  expect_equal(qcombinom(c(0.1, 0.5, 0.9), 12, 0.3, 1),
               qbinom(c(0.1, 0.5, 0.9), 12, 0.3), tolerance = 1e-12)
})

test_that("combinom q(p(x)) round-trip returns x", {
  for (par in list(c(12, 0.35, 1.7), c(20, 0.6, 0.5), c(15, 0.4, -0.5))) {
    s <- par[1]; pr <- par[2]; nu <- par[3]
    expect_equal(qcombinom(pcombinom(0:s, s, pr, nu), s, pr, nu), 0:s,
                 label = paste("round-trip at size =", s, "nu =", nu))
  }
})

test_that("combinom handles edge cases", {
  expect_equal(pcombinom(-1, 12, 0.35, 1.7), 0)
  expect_equal(pcombinom(12, 12, 0.35, 1.7), 1, tolerance = 1e-12)
  expect_equal(pcombinom(99, 12, 0.35, 1.7), 1, tolerance = 1e-12)
  # outside the support, and non-integer values, have zero mass outside of AD
  expect_equal(dcombinom(-1, 12, 0.35, 1.7), 0)
  expect_equal(dcombinom(13, 12, 0.35, 1.7), 0)
  expect_equal(dcombinom(3.5, 12, 0.35, 1.7), 0)
  # size = 0 puts all mass on zero
  expect_equal(dcombinom(0, 0, 0.4, 1.7), 1)
  expect_equal(qcombinom(0, 12, 0.35, 1.7), 0)
  expect_equal(qcombinom(1, 12, 0.35, 1.7), 12)
})

test_that("combinom recycles size, prob and nu", {
  expect_equal(dcombinom(c(2, 3), c(5, 10), 0.4, 1.2),
               c(dcombinom(2, 5, 0.4, 1.2), dcombinom(3, 10, 0.4, 1.2)))
  expect_equal(dcombinom(2, 5, 0.4, c(1, 2)),
               c(dcombinom(2, 5, 0.4, 1), dcombinom(2, 5, 0.4, 2)))
  expect_equal(pcombinom(c(2, 3), c(5, 10), c(0.4, 0.6), 1.2),
               c(pcombinom(2, 5, 0.4, 1.2), pcombinom(3, 10, 0.6, 1.2)))
})

test_that("combinom rejects invalid arguments", {
  expect_error(dcombinom(2, -1, 0.4, 1), "size")
  expect_error(dcombinom(2, 5.5, 0.4, 1), "size")
  expect_error(dcombinom(2, 5, 0, 1), "prob")
  expect_error(dcombinom(2, 5, 1, 1), "prob")
  expect_error(dcombinom(2, 5, 0.4, Inf), "nu")
})

test_that("rcombinom matches the pmf", {
  set.seed(42)
  n <- 2e5
  r <- rcombinom(n, 12, 0.35, 1.7)
  expect_true(all(r >= 0 & r <= 12 & r == floor(r)))
  emp <- as.numeric(table(factor(r, levels = 0:12))) / n
  expect_equal(emp, dcombinom(0:12, 12, 0.35, 1.7), tolerance = 0.01)
})

# ---------------------------------------------------------------------------
# AD behaviour
# ---------------------------------------------------------------------------

# check_ad_gradient() is not used here because it tapes *every* parameter
# passed to it, including size; dcombinom deliberately rejects an advector
# size, since the support 0:size cannot depend on a parameter.
test_that("combinom AD gradient has no NaN and matches finite differences", {
  set.seed(42)
  size <- 12
  x <- rcombinom(20, size, 0.35, 1.7)

  nll <- function(par) -sum(dcombinom(x, size, RTMB::plogis(par[1]), par[2], log = TRUE))
  environment(nll) <- environment()

  F <- RTMB::MakeTape(nll, c(0, 1))
  grad <- F$jacobian(c(0, 1))
  expect_false(any(is.nan(grad)), label = "no NaN in AD gradient")

  h <- 1e-6
  fd <- vapply(1:2, function(j) {
    e <- numeric(2); e[j] <- h
    (nll(c(0, 1) + e) - nll(c(0, 1) - e)) / (2 * h)
  }, numeric(1))
  expect_equal(as.numeric(grad), fd, tolerance = 1e-5)
})

# This is the test that catches a density whose tape silently freezes values
# recorded at the taping point (as e.g. an in-tape root finder would).
test_that("combinom tape replays correctly away from the recording point", {
  set.seed(42)
  size <- 12
  x <- rcombinom(20, size, 0.35, 1.7)

  nll <- function(par) -sum(dcombinom(x, size, RTMB::plogis(par[1]), par[2], log = TRUE))
  environment(nll) <- environment()
  F <- RTMB::MakeTape(nll, c(0, 1))

  for (pt in list(c(0, 1), c(-0.6, 1.7), c(1.2, 0.4), c(-1.5, 3), c(0.3, -0.5))) {
    expect_equal(F(pt), nll(pt), tolerance = 1e-10,
                 label = paste("tape replay at", paste(pt, collapse = ", ")))
    gad <- as.numeric(F$jacobian(pt))
    h <- 1e-6
    fd <- vapply(1:2, function(j) {
      e <- numeric(2); e[j] <- h
      (nll(pt + e) - nll(pt - e)) / (2 * h)
    }, numeric(1))
    expect_equal(gad, fd, tolerance = 1e-4,
                 label = paste("gradient at", paste(pt, collapse = ", ")))
  }
})

test_that("combinom tapes when only some arguments are advectors", {
  # The internal accumulators must pick up the advector type from whichever
  # of prob / nu is a parameter, not just from one of them.
  xs <- c(2, 5, 7)
  cases <- list(
    "nu only, d"   = list(f = function(nu) -sum(dcombinom(xs, 12, 0.35, nu, log = TRUE)), pt = 1.7),
    "nu only, p"   = list(f = function(nu) sum(pcombinom(xs, 12, 0.35, nu)), pt = 1.7),
    "prob only, d" = list(f = function(z) -sum(dcombinom(xs, 12, RTMB::plogis(z), 1.7, log = TRUE)), pt = 0.2),
    "prob only, p" = list(f = function(z) sum(pcombinom(xs, 12, RTMB::plogis(z), 1.7)), pt = 0.2),
    "vector size"  = list(f = function(z) -sum(dcombinom(xs, c(8, 12, 15), RTMB::plogis(z), 1.7, log = TRUE)), pt = 0.2),
    "vector nu"    = list(f = function(nu) -sum(dcombinom(xs, 12, 0.35, rep(nu, 3), log = TRUE)), pt = 1.7)
  )
  for (nm in names(cases)) {
    f <- cases[[nm]]$f; pt <- cases[[nm]]$pt
    environment(f) <- environment()
    F <- RTMB::MakeTape(f, pt)
    expect_equal(F(pt), f(pt), tolerance = 1e-10, label = paste("value:", nm))
    expect_false(any(is.nan(F$jacobian(pt))), label = paste("no NaN gradient:", nm))
  }
})

test_that("dcombinom is differentiable with respect to x", {
  # The pmf is written with lgamma, so it extends smoothly in x; this is what
  # makes OSA residuals possible. Compare against the analytic derivative
  # d/dx log f(x) = nu * (digamma(size - x + 1) - digamma(x + 1)) + logit(prob).
  size <- 12; prob <- 0.35; nu <- 1.7
  F <- RTMB::MakeTape(function(x) dcombinom(x, size, prob, nu, log = TRUE), 4)

  for (x0 in c(1, 4, 7, 10)) {
    analytic <- nu * (digamma(size - x0 + 1) - digamma(x0 + 1)) +
      log(prob) - log1p(-prob)
    expect_equal(as.numeric(F$jacobian(x0)), analytic, tolerance = 1e-8,
                 label = paste("d/dx log-pmf at x =", x0))
    expect_equal(F(x0), dcombinom(x0, size, prob, nu, log = TRUE),
                 tolerance = 1e-12, label = paste("value at x =", x0))
  }
})

test_that("pcombinom accepts an advector q and is exact on the integers", {
  size <- 12; prob <- 0.35; nu <- 1.7
  F <- RTMB::MakeTape(function(q) pcombinom(q, size, prob, nu), 4)

  for (q0 in c(0, 4, 7, 12)) {
    expect_equal(F(q0), pcombinom(q0, size, prob, nu), tolerance = 1e-12,
                 label = paste("CDF value at q =", q0))
  }
  # a discrete CDF is a step function, so its derivative in q is zero
  expect_equal(as.numeric(F$jacobian(4)), 0)
})

test_that("combinom supports simulation and OSA residuals", {
  skip_on_cran()

  set.seed(123)
  size <- 12
  y <- rcombinom(60, size, 0.35, 1.7)
  dat <- list(y = y, size = size)

  fn <- function(par) {
    RTMB::getAll(par, dat)
    y <- RTMB::OBS(y)
    -sum(dcombinom(y, size, RTMB::plogis(lp), nu, log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lp = 0, nu = 1), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(opt$convergence, 0)

  # simulation goes through dGenericSim -> rcombinom
  set.seed(1)
  sim <- obj$simulate()
  expect_true(all(sim$y >= 0 & sim$y <= size & sim$y == floor(sim$y)))

  # both OSA methods need d (resp. p) to be AD-able in x
  for (meth in c("oneStepGeneric", "cdf")) {
    res <- RTMB::oneStepPredict(obj, method = meth, discrete = TRUE, trace = FALSE)
    expect_false(any(is.nan(res$residual)), label = paste("no NaN,", meth))
    expect_lt(abs(mean(res$residual)), 0.4)
    expect_lt(abs(sd(res$residual) - 1), 0.4)
  }
})
