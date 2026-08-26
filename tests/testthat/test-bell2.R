# Tests for the Bell distribution parameterised by its mean.
# theta = W(mu) with W the principal branch of the Lambert W function, so the
# checks below are largely about that transformation being exact and AD-able.

test_that("bell2 passes discrete distribution checks (mu = 1)", {
  check_discrete_dist(
    dfun = dbell2, pfun = pbell2,
    xs_int = c(0, 1, 2, 5, 10), sum_support = 0:200,
    mu = 1
  )
})

test_that("bell2 passes discrete distribution checks (mu = 15)", {
  check_discrete_dist(
    dfun = dbell2, pfun = pbell2,
    xs_int = c(0, 5, 15, 30, 60), sum_support = 0:400,
    mu = 15
  )
})

test_that("bell2 agrees with bell at theta = W(mu)", {
  for (mu in c(0.5, 1, 3, 20, 500)) {
    theta <- lambertW(mu)
    expect_equal(dbell2(0:20, mu), dbell(0:20, theta), tolerance = 1e-14,
                 label = paste("pmf at mu =", mu))
    expect_equal(pbell2(0:20, mu), pbell(0:20, theta), tolerance = 1e-14,
                 label = paste("cdf at mu =", mu))
    expect_equal(qbell2(c(0.1, 0.5, 0.9), mu), qbell(c(0.1, 0.5, 0.9), theta),
                 label = paste("quantiles at mu =", mu))
  }
})

test_that("mu really is the mean", {
  for (mu in c(0.5, 3, 20)) {
    xs <- 0:ceiling(mu + 30 * sqrt(mu * (1 + lambertW(mu))) + 60)
    d <- dbell2(xs, mu)
    expect_equal(sum(xs * d), mu, tolerance = 1e-8,
                 label = paste("mean at mu =", mu))
    # variance is mu * (1 + W(mu))
    expect_equal(sum((xs - mu)^2 * d), mu * (1 + lambertW(mu)), tolerance = 1e-7,
                 label = paste("variance at mu =", mu))
  }
})

test_that("bell2 recycles and rejects invalid arguments", {
  expect_equal(dbell2(c(1, 2), c(1, 5)), c(dbell2(1, 1), dbell2(2, 5)))
  expect_error(dbell2(1, 0), "mu")
  expect_error(dbell2(1, -1), "mu")
  expect_error(qbell2(0.5, 0), "mu")
  expect_error(rbell2(1, 0), "mu")
})

test_that("qbell2 inverts pbell2", {
  for (mu in c(0.5, 3, 20)) {
    xs <- 0:qbell2(0.9999, mu)
    expect_equal(qbell2(pbell2(xs, mu), mu), as.numeric(xs),
                 label = paste("round-trip at mu =", mu))
  }
})

test_that("rbell2 matches the pmf and has mean mu", {
  set.seed(42)
  n <- 3e5
  r <- rbell2(n, 4)
  expect_true(all(r >= 0 & r == floor(r)))
  expect_equal(as.numeric(table(factor(r, levels = 0:14))) / n, dbell2(0:14, 4),
               tolerance = 0.005)
  expect_equal(mean(r), 4, tolerance = 0.05)
})

# ---------------------------------------------------------------------------
# Lambert W
# ---------------------------------------------------------------------------

test_that("lambertW solves W exp(W) = x on the principal branch", {
  x <- c(0, 1e-300, 1e-16, 0.1, 0.5, 1, exp(1), 10, 1e3, 1e8, 1e50, 1e300)
  w <- lambertW(x)
  expect_equal(w * exp(w), x, tolerance = 1e-13)
  expect_equal(lambertW(exp(1)), 1, tolerance = 1e-14)
  expect_equal(lambertW(0), 0)
})

test_that("lambertW handles the negative part of the principal branch", {
  x <- c(-exp(-1), -0.36, -0.3, -0.2, -0.1)
  w <- lambertW(x)
  expect_equal(w * exp(w), x, tolerance = 1e-12)
  expect_equal(lambertW(-exp(-1)), -1, tolerance = 1e-8)
  expect_true(all(w >= -1))
  # below -1/e the principal branch is undefined
  expect_warning(v <- lambertW(-0.5), "principal branch")
  expect_true(is.nan(v))
  expect_true(is.nan(lambertW(NA)))
})

test_that("lambertW is differentiable to third order under AD", {
  # the gradient of a Laplace approximation involves the third derivative
  # tensor, so all three orders have to be right
  for (x0 in c(1e-8, 0.5, 2.5, 1e4)) {
    w <- lambertW(x0)
    a1 <- exp(-w) / (1 + w)
    a2 <- -exp(-2 * w) * (2 + w) / (1 + w)^3
    a3 <- exp(-3 * w) * (2 * w^2 + 8 * w + 9) / (1 + w)^5

    F1 <- RTMB::MakeTape(function(z) lambertW(z), x0)
    F2 <- F1$jacfun()
    F3 <- F2$jacfun()
    F4 <- F3$jacfun()

    expect_equal(as.numeric(F1(x0)), w, tolerance = 1e-14, label = paste("W at", x0))
    expect_equal(as.numeric(F2(x0)), a1, tolerance = 1e-10, label = paste("W' at", x0))
    expect_equal(as.numeric(F3(x0)), a2, tolerance = 1e-10, label = paste("W'' at", x0))
    expect_equal(as.numeric(F4(x0)), a3, tolerance = 1e-10, label = paste("W''' at", x0))
  }
})

test_that("the lambertW tape replays away from its recording point", {
  # catches an atomic that froze the values recorded when the tape was built
  F <- RTMB::MakeTape(function(z) sum(lambertW(z)), c(1, 2, 3))
  for (pt in list(c(1, 2, 3), c(0.01, 50, 1e6), c(5, 5, 5))) {
    expect_equal(F(pt), sum(lambertW(pt)), tolerance = 1e-12,
                 label = paste("replay at", paste(pt, collapse = ", ")))
  }
  expect_false(any(is.nan(F$jacobian(c(0.01, 50, 1e6)))))
})

# ---------------------------------------------------------------------------
# AD behaviour
# ---------------------------------------------------------------------------

test_that("bell2 AD gradient has no NaN", {
  check_ad_gradient(dbell2, rbell2, mu = 4)
})

test_that("dbell2 AD gradient matches finite differences and the tape replays", {
  set.seed(42)
  y <- rbell2(40, 4)
  nll <- function(par) -sum(dbell2(y, exp(par), log = TRUE))
  environment(nll) <- environment()
  F <- RTMB::MakeTape(nll, 1.4)

  h <- 1e-6
  for (pt in c(0.5, 1.4, 2.5, 3.5)) {
    expect_equal(F(pt), nll(pt), tolerance = 1e-9,
                 label = paste("tape replay at", pt))
    expect_equal(as.numeric(F$jacobian(pt)), (nll(pt + h) - nll(pt - h)) / (2 * h),
                 tolerance = 1e-4, label = paste("gradient at", pt))
  }
})

test_that("dbell2 works with a per-observation mean", {
  set.seed(42)
  n <- 40
  y <- rbell2(n, 4)
  z <- seq_len(n) / n
  f <- function(b) -sum(dbell2(y, exp(b[1] + b[2] * z), log = TRUE))
  environment(f) <- environment()
  F <- RTMB::MakeTape(f, c(1.4, 0.1))
  g <- as.numeric(F$jacobian(c(1.4, 0.1)))
  expect_false(any(is.nan(g)))

  h <- 1e-6
  fd <- vapply(1:2, function(j) {
    e <- numeric(2); e[j] <- h
    (f(c(1.4, 0.1) + e) - f(c(1.4, 0.1) - e)) / (2 * h)
  }, numeric(1))
  expect_equal(g, fd, tolerance = 1e-4)
})

test_that("the dbell2 MLE of mu is the sample mean", {
  set.seed(42)
  y <- rbell2(200, 4)
  obj <- RTMB::MakeADFun(function(p) {
    RTMB::getAll(p)
    -sum(dbell2(y, exp(logmu), log = TRUE))
  }, list(logmu = log(3)), silent = TRUE)
  o <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(o$convergence, 0)
  expect_equal(exp(o$par[[1]]), mean(y), tolerance = 1e-6)
})

test_that("bell2 works with random effects under the Laplace approximation", {
  skip_on_cran()

  # grouped random effects, so that the random-effect variance is identified
  # separately from the overdispersion the Bell distribution already has
  set.seed(1)
  m <- 60
  rep_n <- 20
  u0 <- rnorm(m, 0, 0.6)
  g <- rep(seq_len(m), each = rep_n)
  y <- rbell2(m * rep_n, exp(1.2 + u0[g]))
  dat <- list(y = y, g = g)

  fn <- function(p) {
    RTMB::getAll(p, dat)
    -sum(dnorm(u, 0, exp(logsd), log = TRUE)) -
      sum(dbell2(y, exp(b + u[g]), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(b = 1, logsd = log(0.5), u = numeric(m)),
                         random = "u", silent = TRUE)

  # the Laplace gradient is where the third derivatives of lambertW enter
  pv <- c(1, log(0.5))
  step <- 1e-5
  fd <- vapply(1:2, function(j) {
    e <- numeric(2); e[j] <- step
    (obj$fn(pv + e) - obj$fn(pv - e)) / (2 * step)
  }, numeric(1))
  expect_equal(as.numeric(obj$gr(pv)), fd, tolerance = 1e-4)

  # the gradient check above is the exactness test; these are loose
  # finite-sample checks that the fit lands in the right neighbourhood
  o <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(o$convergence, 0)
  expect_equal(o$par[[1]], 1.2, tolerance = 0.2)
  expect_equal(exp(o$par[[2]]), 0.6, tolerance = 0.35)

  sdr <- RTMB::sdreport(obj)
  expect_true(all(is.finite(sqrt(diag(sdr$cov.fixed)))))
})
