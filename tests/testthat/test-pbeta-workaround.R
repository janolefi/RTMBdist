# Tests for beta_ad() in R/pbeta-workaround.R, which works around wrong derivatives of
# RTMB's pbeta() and qbeta() at integer and half-integer shapes. Delete with the workaround.

fd_grad <- function(f, p, h = 1e-5) {
  sapply(seq_along(p), function(j) {
    e <- replace(0 * p, j, h)
    (f(p + e) - f(p - e)) / (2 * h)
  })
}
fd_hess <- function(f, p, h = 1e-4) {
  d <- length(p)
  outer(seq_len(d), seq_len(d), Vectorize(function(i, j) {
    ei <- replace(0 * p, i, h); ej <- replace(0 * p, j, h)
    (f(p + ei + ej) - f(p + ei - ej) - f(p - ei + ej) + f(p - ei - ej)) / (4 * h^2)
  }))
}

test_that("beta_ad gives correct values and derivatives at integer and half-integer shapes", {
  for (p in list(c(0.3, 2, 3), c(0.8, 4, 2), c(0.3, 3, 7), c(0.3, 2.5, 4), c(0.1, 0.5, 3),
                 c(0.6, 1.3, 0.7))) {
    ref <- function(p) stats::pbeta(p[1], p[2], p[3])
    F <- RTMB::MakeTape(function(p) beta_ad(pbeta, p[1], p[2], p[3]), p)
    expect_equal(F(p), ref(p), tolerance = 1e-11)
    expect_equal(as.vector(F$jacobian(p)), fd_grad(ref, p), tolerance = 1e-7)
    expect_equal(F$jacfun()$jacobian(p), fd_hess(ref, p), tolerance = 1e-4)

    refq <- function(p) stats::qbeta(p[1], p[2], p[3])
    G <- RTMB::MakeTape(function(p) beta_ad(qbeta, p[1], p[2], p[3]), p)
    expect_equal(G(p), refq(p), tolerance = 1e-10)
    expect_equal(as.vector(G$jacobian(p)), fd_grad(refq, p), tolerance = 1e-6)
  }
  # outside AD it is the plain function
  expect_identical(beta_ad(pbeta, 0.3, 2, 3), stats::pbeta(0.3, 2, 3))
})

test_that("the distribution functions using pbeta have correct gradients at integer shapes", {
  # shapes 4 and 6, where RTMB's pbeta has a zero derivative in shape2
  check_ad_cdf(pbeta2, dbeta2, c(0.2, 0.5, 0.8), mu = 0.4, phi = 10)
  # shape2 = q + 1 is an integer
  check_ad_cdf(pnbinom2, dnbinom2, c(0, 2, 5), mu = 3, size = 2, .dq = FALSE)
})
