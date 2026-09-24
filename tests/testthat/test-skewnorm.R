# Tests for the skew normal distribution

test_that("skewnorm passes standard distribution checks (xi=0, omega=1, alpha=0)", {
  check_continuous_dist(
    dfun  = dskewnorm,
    pfun  = pskewnorm,
    qfun  = qskewnorm,
    xs    = c(-2, -1, 0, 1, 2),
    xi = 0, omega = 1, alpha = 0
  )
})

test_that("skewnorm passes standard distribution checks (xi=1, omega=2, alpha=3)", {
  check_continuous_dist(
    dfun  = dskewnorm,
    pfun  = pskewnorm,
    qfun  = qskewnorm,
    xs    = c(0.7, 1.6, 2.3, 3.3, 4.9),
    xi = 1, omega = 2, alpha = 3
  )
})

test_that("skewnorm AD gradient has no NaN", {
  check_ad_gradient(dskewnorm,  rskewnorm,  xi = 0, omega = 1, alpha = 2)
})

test_that("pskewnorm under AD matches sn::psn in value and gradient", {
  q   <- c(-3, -1, 0, 0.5, 2, 6)
  par <- c(xi = 0.3, omega = 1.5, alpha = 3)

  F <- RTMB::MakeTape(function(p) pskewnorm(q, xi = p[1], omega = p[2], alpha = p[3]), par)
  expect_equal(F(par), sn::psn(q, 0.3, 1.5, 3), tolerance = 1e-10)

  # tape stays correct at other parameters, including a different skew direction
  expect_equal(F(c(-2, 0.7, -1)), sn::psn(q, -2, 0.7, -1), tolerance = 1e-10)

  h <- 1e-5
  J_fd <- sapply(1:3, function(j) {
    e <- replace(numeric(3), j, h)
    (sn::psn(q, par[1] + e[1], par[2] + e[2], par[3] + e[3]) -
       sn::psn(q, par[1] - e[1], par[2] - e[2], par[3] - e[3])) / (2 * h)
  })
  expect_equal(F$jacobian(par), J_fd, tolerance = 1e-7, ignore_attr = TRUE)

  # derivative with respect to q is the density
  G <- RTMB::MakeTape(function(x) pskewnorm(x, xi = 0.3, omega = 1.5, alpha = 3), q)
  expect_equal(diag(G$jacobian(q)), dskewnorm(q, 0.3, 1.5, 3), tolerance = 1e-10)
})

test_that("pskewnorm under AD handles tails, log.p and recycling", {
  par <- c(xi = 0.3, omega = 1.5, alpha = 3)

  F <- RTMB::MakeTape(function(p)
    pskewnorm(c(-4, 1, 8), xi = p[1], omega = p[2], alpha = p[3],
              lower.tail = FALSE, log.p = TRUE), par)
  expect_equal(F(par), log1p(-sn::psn(c(-4, 1, 8), 0.3, 1.5, 3)), tolerance = 1e-8)

  q <- c(-1, 0, 1, 2)
  G <- RTMB::MakeTape(function(p)
    pskewnorm(q, xi = p[1:2], omega = 1.2, alpha = p[3:6]), c(0, 1, -2, 0, 1, 5))
  expect_equal(G(c(0, 1, -2, 0, 1, 5)), sn::psn(q, c(0, 1, 0, 1), 1.2, c(-2, 0, 1, 5)),
               tolerance = 1e-10)

  H <- RTMB::MakeTape(function(p) pskewnorm(c(-Inf, Inf), xi = p[1], omega = p[2], alpha = p[3]), par)
  expect_equal(H(par), c(0, 1))
})

test_that("pskewnorm under AD is accurate for a narrow density far from zero", {
  xi <- 100; omega <- 0.002; alpha <- 4
  q <- xi + omega * c(-5, 0, 0.5, 1.5, 5)
  z <- (q - xi) / omega
  F <- RTMB::MakeTape(function(p) pskewnorm(q, p[1], p[2], p[3]), c(xi, omega, alpha))
  J <- F$jacobian(c(xi, omega, alpha))
  expect_equal(F(c(xi, omega, alpha)), sn::psn(z, 0, 1, alpha), tolerance = 1e-8)
  # finite differences in xi and omega are inaccurate at this scale, but in a location-scale
  # family dF/dxi = -f(q) and dF/domega = -z f(q), and dF/dalpha is that of the standard case
  expect_equal(J[, 1], -dskewnorm(q, xi, omega, alpha), tolerance = 1e-8)
  expect_equal(J[, 2], -z * dskewnorm(q, xi, omega, alpha), tolerance = 1e-8)
  h <- 1e-5
  expect_equal(J[, 3], (sn::psn(z, 0, 1, alpha + h) - sn::psn(z, 0, 1, alpha - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("skewnorm supports OSA residuals", {
  set.seed(1)
  check_osa_cdf(dskewnorm, pskewnorm, rskewnorm(30, 1, 2, 3), xi = 1, omega = 2, alpha = 3)
})
