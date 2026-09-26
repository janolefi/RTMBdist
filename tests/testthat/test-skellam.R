# Tests for the Skellam distribution
# Support is all integers; there is no quantile function.

test_that("skellam passes discrete distribution checks (mu1=2, mu2=3)", {
  check_discrete_dist(
    dfun        = dskellam,
    pfun        = pskellam,
    xs_int      = c(-5, -2, 0, 2, 4),
    sum_support = -100:100,
    mu1 = 2, mu2 = 3
  )
})

test_that("skellam passes discrete distribution checks (mu1=5, mu2=1)", {
  check_discrete_dist(
    dfun        = dskellam,
    pfun        = pskellam,
    xs_int      = c(-2, 0, 3, 6, 9),
    sum_support = -100:100,
    mu1 = 5, mu2 = 1
  )
})

test_that("skellam AD gradient has no NaN", {
  check_ad_gradient(dskellam,   rskellam,   mu1 = 5, mu2 = 2)
})

# independent reference: X = N1 - N2, so P(X <= q) = sum_n P(N2 = n) P(N1 <= q + n),
# with the upper tail summed directly as well
ref_pskellam <- function(q, mu1, mu2, lower.tail = TRUE) {
  n <- 0:(ceiling(mu2 + 40 * sqrt(mu2)) + 60)
  sapply(q, function(k) sum(stats::dpois(n, mu2) * stats::ppois(k + n, mu1, lower.tail = lower.tail)))
}

test_that("pskellam matches the Poisson convolution", {
  for (p in list(c(2, 3), c(0.05, 0.3), c(7, 0.2), c(200, 180), c(1e-3, 1e-3))) {
    m <- p[1] - p[2]; s <- sqrt(p[1] + p[2])
    q <- unique(round(m + c(-6, -3, -1, 0, 1, 3, 6) * s))
    expect_equal(pskellam(q, p[1], p[2]), ref_pskellam(q, p[1], p[2]), tolerance = 1e-9,
                 label = paste("lower tail, mu =", p[1], p[2]))
    expect_equal(pskellam(q, p[1], p[2], lower.tail = FALSE),
                 ref_pskellam(q, p[1], p[2], lower.tail = FALSE), tolerance = 1e-9,
                 label = paste("upper tail, mu =", p[1], p[2]))
  }
})

test_that("pskellam keeps small tail probabilities accurate", {
  # relative accuracy far out in both tails, where 1 - p would lose all digits
  q_lo <- c(-25, -15); q_hi <- c(15, 25)
  expect_equal(pskellam(q_lo, 2, 3) / ref_pskellam(q_lo, 2, 3), c(1, 1), tolerance = 1e-8)
  expect_equal(pskellam(q_hi, 2, 3, lower.tail = FALSE) / ref_pskellam(q_hi, 2, 3, FALSE),
               c(1, 1), tolerance = 1e-8)
  # log.p far out in the tails, compared with a log-sum-exp of the pmf
  q <- -40
  lp <- dskellam(-150:q, 2, 3, log = TRUE)
  expect_equal(pskellam(q, 2, 3, log.p = TRUE), max(lp) + log(sum(exp(lp - max(lp)))),
               tolerance = 1e-10)
  lp <- dskellam(40:150, 2, 3, log = TRUE) # P(X > 39) = P(X >= 40)
  expect_equal(pskellam(-q - 1, 2, 3, lower.tail = FALSE, log.p = TRUE),
               max(lp) + log(sum(exp(lp - max(lp)))), tolerance = 1e-10)
  # where the pmf itself underflows, the tail is zero rather than an error
  expect_equal(suppressWarnings(pskellam(c(-400, 400), 2, 3, log.p = TRUE)), c(-Inf, 0))
})

test_that("pskellam handles edge cases and recycles its arguments", {
  expect_equal(pskellam(c(-Inf, Inf, NA), 2, 3), c(0, 1, NA))
  expect_equal(pskellam(c(-Inf, Inf), 2, 3, lower.tail = FALSE), c(1, 0))
  expect_equal(pskellam(1.7, 2, 3), pskellam(1, 2, 3))
  # symmetry: P(X <= q; mu1, mu2) = P(X >= -q; mu2, mu1)
  q <- -4:4
  expect_equal(pskellam(q, 2, 3), 1 - pskellam(-q - 1, 3, 2), tolerance = 1e-12)
  # parameters per element
  expect_equal(pskellam(c(-1, 2, 0), c(1, 4, 2), c(3, 1, 2)),
               c(pskellam(-1, 1, 3), pskellam(2, 4, 1), pskellam(0, 2, 2)))
  expect_length(pskellam(0:5, c(1, 2), 3), 6)
  expect_error(RTMB::MakeTape(function(x) pskellam(x, 2, 3), 1), "numeric data")
})

test_that("pskellam under AD matches pskellam outside AD", {
  check_ad_cdf(pskellam, dskellam, c(-8, -2, 0, 1, 3, 9), mu1 = 2, mu2 = 3, .dq = FALSE)
  check_ad_cdf(pskellam, dskellam, c(-3, 0, 4), mu1 = 0.3, mu2 = 5, .dq = FALSE,
               .args = list(lower.tail = FALSE, log.p = TRUE))
  # dF/dmu1 = -P(X = q) exactly, from the integral representation
  q <- c(-5, -1, 0, 2, 6)
  J <- RTMB::MakeTape(function(p) pskellam(q, p[1], p[2]), c(2, 3))$jacobian(c(2, 3))
  expect_equal(J[, 1], -dskellam(q, 2, 3), tolerance = 1e-10)
})

test_that("pskellam can be taped in only one of mu1 and mu2", {
  q <- c(-5, -1, 0, 2, 6)
  F1 <- RTMB::MakeTape(function(m) pskellam(q, m, 3), 2)
  F2 <- RTMB::MakeTape(function(m) pskellam(q, 2, m), 3)
  expect_equal(F1(2), pskellam(q, 2, 3))
  expect_equal(F2(3), pskellam(q, 2, 3))
  expect_equal(as.vector(F1$jacobian(2)), -dskellam(q, 2, 3), tolerance = 1e-10)
})

test_that("the pskellam tape stays correct when the other tail would be the smaller one", {
  # the tail integrated directly is chosen at tape time; both are exact, so
  # re-evaluating at parameters that move the mean across q must stay correct
  q <- c(-6, -1, 0, 2, 7)
  F <- RTMB::MakeTape(function(p) pskellam(q, p[1], p[2]), c(2, 3))
  for (p in list(c(9, 1), c(0.5, 8), c(30, 25))) {
    expect_equal(F(p), ref_pskellam(q, p[1], p[2]), tolerance = 1e-9)
  }
})

test_that("pskellam has correct second derivatives, as needed for the Laplace approximation", {
  # several elements in both tails; second derivatives through RTMB's Vectorize() with
  # integrate() were wrong, so this compares with finite differences of values
  q <- c(-6, -4, -1, 0, 3, 5)
  f <- function(p) sum(log(pskellam(q, p[1], p[2])))
  par <- c(2.5, 3.5)
  H <- RTMB::MakeTape(f, par)$jacfun()$jacobian(par)
  h <- 1e-4
  H_fd <- outer(1:2, 1:2, Vectorize(function(i, j) {
    ei <- replace(c(0, 0), i, h); ej <- replace(c(0, 0), j, h)
    (f(par + ei + ej) - f(par + ei - ej) - f(par - ei + ej) + f(par - ei - ej)) / (4 * h^2)
  }))
  expect_equal(H, H_fd, tolerance = 1e-5, ignore_attr = TRUE)
})

test_that("skellam supports OSA residuals, also with estimated parameters", {
  check_osa_cdf(dskellam, pskellam, c(-6, -2, 0, 0, 3, 8), mu1 = 2, mu2 = 3, .discrete = TRUE)

  set.seed(1)
  y <- rskellam(50, 4, 2.5)
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(dskellam(y, exp(lmu1), exp(lmu2), log = TRUE))
  }
  obj <- RTMB::MakeADFun(fn, list(lmu1 = log(4), lmu2 = log(2.5)), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(opt$convergence, 0)
  res <- RTMB::oneStepPredict(obj, method = "cdf", discrete = TRUE, trace = FALSE)
  mu <- exp(opt$par)
  expect_equal(res$Fx, pskellam(y, mu[1], mu[2]), tolerance = 1e-8)
  expect_equal(res$px, dskellam(y, mu[1], mu[2]), tolerance = 1e-8)
  expect_false(any(is.nan(res$residual)))
})
