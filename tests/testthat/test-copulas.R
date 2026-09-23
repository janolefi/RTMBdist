# Tests for the bivariate copula densities
# The densities are checked against mixed finite differences of the copula
# CDFs, d^2 C / du dv, and for uniform margins.

fd_density <- function(C, u, v, h = 1e-4) {
  (C(u + h, v + h) - C(u + h, v - h) - C(u - h, v + h) + C(u - h, v - h)) / (4 * h^2)
}

margins <- function(logc, at = c(0.05, 0.3, 0.5, 0.77, 0.95)) {
  c(sapply(at, function(v)
      integrate(function(u) exp(logc(u, rep(v, length(u)))), 0, 1, rel.tol = 1e-10)$value),
    sapply(at, function(u)
      integrate(function(v) exp(logc(rep(u, length(v)), v)), 0, 1, rel.tol = 1e-10)$value))
}

u <- c(0.1, 0.3, 0.5, 0.62, 0.9); v <- c(0.2, 0.8, 0.5, 0.35, 0.95)

test_that("cgumbel is the density of Cgumbel", {
  for (theta in c(1.2, 2, 5)) {
    expect_equal(exp(cgumbel(theta)(u, v)), fd_density(Cgumbel(theta), u, v), tolerance = 1e-5)
    expect_equal(margins(cgumbel(theta)), rep(1, 10), tolerance = 1e-6)
  }
  # theta = 1 is independence
  expect_equal(cgumbel(1)(u, v), rep(0, 5))
})

test_that("cfrank is the density of Cfrank, for either sign of theta", {
  for (theta in c(-6, -1, 0.5, 3, 10)) {
    expect_equal(exp(cfrank(theta)(u, v)), fd_density(Cfrank(theta), u, v), tolerance = 1e-5)
    expect_equal(margins(cfrank(theta)), rep(1, 10), tolerance = 1e-6)
  }
})

test_that("cfrank stays accurate for strong dependence", {
  # closed form at u = v, where the density is theta (1 - e^-theta) e^(-2 theta u) / den^2
  theta <- 30; w <- c(0.2, 0.5, 0.9)
  den <- 2 * exp(-theta * w) - exp(-2 * theta * w) - exp(-theta)
  expect_equal(cfrank(theta)(w, w),
               log(theta * (1 - exp(-theta))) - 2 * theta * w - 2 * log(den), tolerance = 1e-12)
})

test_that("cgaussian and cclayton are the densities of their CDFs and have uniform margins", {
  expect_equal(exp(cclayton(2)(u, v)), fd_density(Cclayton(2), u, v), tolerance = 1e-5)
  expect_equal(margins(cclayton(2)), rep(1, 10), tolerance = 1e-6)
  expect_equal(margins(cgaussian(0.6)), rep(1, 10), tolerance = 1e-6)
  expect_equal(margins(cgaussian(-0.4)), rep(1, 10), tolerance = 1e-6)
})

test_that("copula densities have AD gradients without NaN", {
  F <- RTMB::MakeTape(function(p) sum(cgumbel(p[1])(u, v)) + sum(cfrank(p[2])(u, v)) +
                        sum(cclayton(p[3])(u, v)) + sum(cgaussian(p[4])(u, v)),
                      c(2, -3, 1.5, 0.3))
  expect_false(any(is.nan(F$jacobian(c(2, -3, 1.5, 0.3)))))
})
