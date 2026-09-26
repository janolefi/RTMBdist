# Tests for the ex-Gaussian distribution
# X = normal(mu, sigma) + exponential(lambda)

test_that("exgauss passes standard distribution checks (mu=0, sigma=1, lambda=1)", {
  check_continuous_dist(
    dfun  = dexgauss,
    pfun  = pexgauss,
    qfun  = qexgauss,
    xs    = c(-1, 0, 1, 2, 4),
    lower = -5, upper = 20,
    mu = 0, sigma = 1, lambda = 1
  )
})

test_that("exgauss passes standard distribution checks (mu=2, sigma=0.5, lambda=2)", {
  check_continuous_dist(
    dfun  = dexgauss,
    pfun  = pexgauss,
    qfun  = qexgauss,
    xs    = c(1, 2, 2.5, 3, 4),
    lower = -1, upper = 15,
    mu = 2, sigma = 0.5, lambda = 2
  )
})

test_that("exgauss AD gradient has no NaN", {
  check_ad_gradient(dexgauss,   rexgauss,   mu = 2, sigma = 0.5, lambda = 2)
})

test_that("dexgauss is accurate far in the left tail, where it used to be floored", {
  # reference: the convolution of the normal and the exponential, integrated on the log scale
  conv <- function(x, mu, sigma, lambda) {
    lf <- function(t) stats::dnorm(x - t, mu, sigma, log = TRUE) - lambda * t + log(lambda)
    m <- lf(0)
    m + log(stats::integrate(function(t) exp(lf(t) - m), 0, Inf, rel.tol = 1e-12)$value)
  }
  for (x in c(-40, -10, 0, 3)) {
    expect_equal(dexgauss(x, 0, 1, 2, log = TRUE), conv(x, 0, 1, 2), tolerance = 1e-10)
  }
  f <- function(p) sum(dexgauss(c(-40, -5, 1), p[1], p[2], p[3], log = TRUE))
  par <- c(0.2, 1.1, 2)
  J <- RTMB::MakeTape(f, par)$jacobian(par)
  h <- 1e-6
  J_fd <- sapply(1:3, function(j) { e <- replace(numeric(3), j, h); (f(par + e) - f(par - e)) / (2 * h) })
  expect_equal(as.vector(J), J_fd, tolerance = 1e-6)
})
