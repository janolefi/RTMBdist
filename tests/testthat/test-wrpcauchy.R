# Tests for the wrapped Cauchy distribution (circular, support [-pi, pi])
# The distribution function has its origin at mu - pi, so p and q are checked
# on (mu - pi, mu + pi].

test_that("wrpcauchy log=TRUE is consistent with log(density) (mu=0, rho=0.5)", {
  xs <- c(-2, -1, 0, 1, 2)
  expect_equal(
    dwrpcauchy(xs, mu = 0, rho = 0.5, log = TRUE),
    log(dwrpcauchy(xs, mu = 0, rho = 0.5)),
    tolerance = 1e-10
  )
})

test_that("wrpcauchy density integrates to 1 (mu=0, rho=0.5)", {
  result <- integrate(dwrpcauchy, lower = -pi, upper = pi, mu = 0, rho = 0.5)
  expect_equal(result$value, 1, tolerance = 1e-4)
})

test_that("wrpcauchy log=TRUE is consistent with log(density) (mu=1, rho=0.8)", {
  xs <- c(-2, -0.5, 1, 1.5, 2.5)
  expect_equal(
    dwrpcauchy(xs, mu = 1, rho = 0.8, log = TRUE),
    log(dwrpcauchy(xs, mu = 1, rho = 0.8)),
    tolerance = 1e-10
  )
})

test_that("wrpcauchy density integrates to 1 (mu=1, rho=0.8)", {
  result <- integrate(dwrpcauchy, lower = -pi, upper = pi, mu = 1, rho = 0.8)
  expect_equal(result$value, 1, tolerance = 1e-4)
})

test_that("wrpcauchy AD gradient has no NaN", {
  check_ad_gradient(dwrpcauchy, rwrpcauchy, mu = 0, rho = 0.5)
})

test_that("wrpcauchy passes standard distribution checks", {
  check_continuous_dist(
    dfun = dwrpcauchy, pfun = pwrpcauchy, qfun = qwrpcauchy,
    xs = c(-3, -1, 0, 0.5, 3), lower = -pi, upper = pi,
    mu = 0, rho = 0.5
  )
  check_continuous_dist(
    dfun = dwrpcauchy, pfun = pwrpcauchy, qfun = qwrpcauchy,
    xs = c(-2, 0, 1, 1.2, 3.5), lower = 1 - pi, upper = 1 + pi,
    mu = 1, rho = 0.8
  )
})

test_that("pwrpcauchy integrates the density from mu - pi", {
  for (mu in c(0, 0.7, -2)) {
    for (rho in c(0.1, 0.6, 0.95)) {
      xs <- mu + c(-3, -1.5, -0.2, 0.4, 2, 3.1)
      num <- sapply(xs, function(x)
        integrate(dwrpcauchy, mu - pi, x, mu = mu, rho = rho, rel.tol = 1e-10)$value)
      expect_equal(pwrpcauchy(xs, mu, rho), num, tolerance = 1e-8)
    }
  }
})

test_that("pwrpcauchy has its origin at mu - pi and is periodic", {
  mu <- 0.7; rho <- 0.6
  expect_equal(pwrpcauchy(mu, mu, rho), 0.5)
  expect_equal(pwrpcauchy(mu - pi, mu, rho), 0)
  expect_equal(pwrpcauchy(mu + pi, mu, rho), 1)
  xs <- mu + c(-3, -1, 0.5, 2.5)
  expect_equal(pwrpcauchy(xs + 2 * pi, mu, rho), pwrpcauchy(xs, mu, rho))
  expect_equal(pwrpcauchy(xs - 4 * pi, mu, rho), pwrpcauchy(xs, mu, rho))
  # on [-pi, pi], the distribution function drops back to 0 at mu - pi
  expect_gt(pwrpcauchy(mu - pi - 0.01, mu, rho), 0.99)
  expect_lt(pwrpcauchy(mu - pi + 0.01, mu, rho), 0.01)
})

test_that("wrpcauchy with rho = 0 is the circular uniform", {
  xs <- c(-3, -1, 0, 2, 3)
  expect_equal(pwrpcauchy(xs, 0, 0), (xs + pi) / (2 * pi))
  expect_equal(qwrpcauchy(c(0.1, 0.5, 0.9), 0, 0), 2 * pi * c(0.1, 0.5, 0.9) - pi)
})

test_that("qwrpcauchy returns angles in [mu - pi, mu + pi]", {
  mu <- 2; rho <- 0.7
  expect_equal(qwrpcauchy(0, mu, rho), mu - pi)
  expect_equal(qwrpcauchy(1, mu, rho), mu + pi)
  expect_equal(qwrpcauchy(0.5, mu, rho), mu)
  expect_true(all(abs(qwrpcauchy(seq(0, 1, by = 0.05), mu, rho) - mu) <= pi))
})

test_that("wrpcauchy p and q honour lower.tail and log.p and recycle", {
  xs <- c(-3, -1, 0, 1, 3)
  expect_equal(pwrpcauchy(xs, 1, 0.5, log.p = TRUE), log(pwrpcauchy(xs, 1, 0.5)))
  expect_equal(pwrpcauchy(xs, 1, 0.5, lower.tail = FALSE, log.p = TRUE),
               log(1 - pwrpcauchy(xs, 1, 0.5)))
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qwrpcauchy(p, 1, 0.5), qwrpcauchy(1 - p, 1, 0.5, lower.tail = FALSE))
  expect_equal(qwrpcauchy(p, 1, 0.5), qwrpcauchy(log(p), 1, 0.5, log.p = TRUE))
  expect_length(pwrpcauchy(0.5, mu = 0, rho = c(0.1, 0.3, 0.5, 0.7)), 4)
  expect_length(qwrpcauchy(0.5, mu = 1:4, rho = 0.5), 4)
  expect_equal(pwrpcauchy(0.5, mu = 0, rho = c(0.1, 0.7)),
               c(pwrpcauchy(0.5, 0, 0.1), pwrpcauchy(0.5, 0, 0.7)))
})

test_that("pwrpcauchy rejects rho outside [0, 1)", {
  expect_error(pwrpcauchy(1, 0, 1), "rho")
  expect_error(qwrpcauchy(0.5, 0, -0.1), "rho")
  expect_error(qwrpcauchy(1.5, 0, 0.5), "p must")
})

test_that("pwrpcauchy AD gradient matches the density and has no NaN", {
  # dF/dq is the density, also at the cut mu +- pi where tan is only finite
  xs <- c(0.3 - pi, -2, 0, 1, 0.3 + pi)
  F <- RTMB::MakeTape(function(q) pwrpcauchy(q, 0.3, 0.6), xs)
  expect_equal(diag(F$jacobian(xs)), dwrpcauchy(xs, 0.3, 0.6), tolerance = 1e-8)

  G <- RTMB::MakeTape(function(par) sum(pwrpcauchy(xs, par[1], par[2])), c(0.3, 0.6))
  expect_false(any(is.nan(G$jacobian(c(0.3, 0.6)))))
  H <- RTMB::MakeTape(function(par) sum(qwrpcauchy(c(0.1, 0.5, 0.9), par[1], par[2])), c(0.3, 0.6))
  expect_false(any(is.nan(H$jacobian(c(0.3, 0.6)))))
})

test_that("rwrpcauchy draws from the wrapped Cauchy and respects wrap", {
  set.seed(1)
  x <- rwrpcauchy(20000, mu = 2.5, rho = 0.6)
  expect_true(all(x >= -pi & x < pi))
  # uniform probability integral transform
  expect_gt(ks.test(pwrpcauchy(x, 2.5, 0.6), "punif")$p.value, 0.01)
  # mean resultant length of the wrapped Cauchy is rho (standard error 0.004)
  expect_equal(abs(mean(exp(1i * x))), 0.6, tolerance = 0.05)
  expect_equal(Arg(mean(exp(1i * x))), 2.5, tolerance = 0.05)

  y <- rwrpcauchy(1000, mu = 2.5, rho = 0.6, wrap = FALSE)
  expect_true(all(abs(y - 2.5) <= pi))
  expect_length(rwrpcauchy(4, mu = 1:4, rho = 0.5), 4)
  expect_error(rwrpcauchy(5, 0, 1), "rho")
})

test_that("pwrpcauchy with from integrates the density from that origin", {
  for (mu in c(0, 2)) {
    for (from in c(-pi, 0.4, mu + pi, 10)) {
      xs <- from + c(0.3, 1, 2.5, 4, 6)
      num <- sapply(xs, function(x)
        integrate(dwrpcauchy, from, x, mu = mu, rho = 0.7, rel.tol = 1e-10)$value)
      expect_equal(pwrpcauchy(xs, mu, 0.7, from = from), num, tolerance = 1e-8)
      # periodic in q
      expect_equal(pwrpcauchy(xs - 2 * pi, mu, 0.7, from = from),
                   pwrpcauchy(xs, mu, 0.7, from = from), tolerance = 1e-10)
    }
  }
})

test_that("qwrpcauchy with from inverts pwrpcauchy on [from, from + 2 pi]", {
  p <- c(0, 0.01, 0.3, 0.5, 0.8, 0.99, 1)
  for (mu in c(0, 2, -3)) {
    for (from in c(-pi, 0.4, mu - pi, mu + pi, 10, -7)) {
      q <- qwrpcauchy(p, mu, 0.6, from = from)
      expect_equal(q[1], from)
      expect_equal(q[length(q)], from + 2 * pi)
      expect_true(all(diff(q) > 0))
      expect_equal(pwrpcauchy(q[-c(1, length(q))], mu, 0.6, from = from),
                   p[-c(1, length(p))], tolerance = 1e-10)
      xs <- from + c(0.2, 1, 3, 5, 6)
      expect_equal(qwrpcauchy(pwrpcauchy(xs, mu, 0.6, from = from), mu, 0.6, from = from),
                   xs, tolerance = 1e-8)
    }
  }
  expect_equal(qwrpcauchy(p, 1, 0.6, from = -pi),
               qwrpcauchy(1 - p, 1, 0.6, from = -pi, lower.tail = FALSE))
})

test_that("from = mu - pi reproduces the default origin", {
  xs <- c(-3, -1, 0.5, 2, 3)
  expect_equal(pwrpcauchy(xs, 1, 0.5, from = 1 - pi), pwrpcauchy(xs, 1, 0.5))
  expect_equal(pwrpcauchy(xs, 1, 0.5, from = 1 - pi, lower.tail = FALSE),
               pwrpcauchy(xs, 1, 0.5, lower.tail = FALSE))
  p <- c(0.05, 0.5, 0.95)
  expect_equal(qwrpcauchy(p, 1, 0.5, from = 1 - pi), qwrpcauchy(p, 1, 0.5))
})

test_that("a fixed from gives uniform PIT values for a mixture over mu", {
  # with the default origin, each mu cuts the circle elsewhere and the averaged
  # distribution function is not the distribution function of the mixture
  set.seed(1)
  n <- 20000
  mu <- sample(c(0, 2), n, replace = TRUE)
  x <- rwrpcauchy(n, mu, 0.5)
  u_fixed <- 0.5 * pwrpcauchy(x, 0, 0.5, from = -pi) + 0.5 * pwrpcauchy(x, 2, 0.5, from = -pi)
  u_default <- 0.5 * pwrpcauchy(x, 0, 0.5) + 0.5 * pwrpcauchy(x, 2, 0.5)
  expect_gt(ks.test(u_fixed, "punif")$p.value, 0.01)
  expect_lt(ks.test(u_default, "punif")$p.value, 1e-10)
})

test_that("pwrpcauchy with from has correct AD derivatives", {
  xs <- c(-3, -1, 0, 1, 3)
  F <- RTMB::MakeTape(function(q) pwrpcauchy(q, 0.3, 0.6, from = -pi), xs)
  expect_equal(diag(F$jacobian(xs)), dwrpcauchy(xs, 0.3, 0.6), tolerance = 1e-8)
  # the derivative with respect to the origin is minus the density there
  G <- RTMB::MakeTape(function(f) pwrpcauchy(1, 0.3, 0.6, from = f), -2)
  expect_equal(G$jacobian(-2)[1, 1], -dwrpcauchy(-2, 0.3, 0.6), tolerance = 1e-8)
  H <- RTMB::MakeTape(function(par) sum(qwrpcauchy(c(0.1, 0.5, 0.9), par[1], par[2], from = -pi)),
                      c(0.3, 0.6))
  expect_false(any(is.nan(H$jacobian(c(0.3, 0.6)))))
})

test_that("pwrpcauchy with a fixed origin stays correct when the tape is re-evaluated", {
  q <- c(-2, 0, 2)
  F <- RTMB::MakeTape(function(m) pwrpcauchy(q, m, 0.5, from = -pi), 0)
  expect_equal(F(2.5), pwrpcauchy(q, 2.5, 0.5, from = -pi), tolerance = 1e-10)
  expect_equal(F(-2.5), pwrpcauchy(q, -2.5, 0.5, from = -pi), tolerance = 1e-10)
})

test_that("wrpcauchy supports OSA residuals, with the circle cut at -pi", {
  set.seed(1)
  check_osa_cdf(dwrpcauchy, function(q, mu, rho) pwrpcauchy(q, mu, rho, from = -pi),
                rwrpcauchy(30, 1, 0.6), mu = 1, rho = 0.6)
})

test_that("wrpcauchy OSA residuals are valid for mixtures over different mean directions", {
  set.seed(1)
  y <- c(rwrpcauchy(20, -2, 0.7), rwrpcauchy(20, 2, 0.4))
  fn <- function(par) {
    RTMB::getAll(par)
    y <- RTMB::OBS(y)
    -sum(log(0.3 * exp(dwrpcauchy(y, mu1, rho1, log = TRUE)) +
             0.7 * exp(dwrpcauchy(y, mu2, rho2, log = TRUE))))
  }
  par <- list(mu1 = -2, rho1 = 0.7, mu2 = 2, rho2 = 0.4)
  obj <- RTMB::MakeADFun(fn, par, silent = TRUE)
  res <- RTMB::oneStepPredict(obj, method = "cdf", trace = FALSE)
  Fmix <- 0.3 * pwrpcauchy(y, -2, 0.7, from = -pi) + 0.7 * pwrpcauchy(y, 2, 0.4, from = -pi)
  expect_equal(res$residual, qnorm(Fmix), tolerance = 1e-6)
})
