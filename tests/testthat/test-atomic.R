# Tests for the atomic tapes in R/atomic.R

test_that("pgamma_ad matches stats::pgamma, in and outside AD context", {
  x <- c(-1, 0, 1e-8, 0.3, 1, 1.999, 2, 2.001, 7, 60, 1e3, Inf)
  for (a in c(0.01, 0.5, 1, 3, 50, 1000)) {
    expect_equal(pgamma_ad(x, a), stats::pgamma(x, a), tolerance = 1e-14)
    F <- RTMB::MakeTape(function(p) pgamma_ad(x, p), a)
    expect_equal(F(a), stats::pgamma(x, a), tolerance = 1e-14)
  }
})

test_that("pgamma_ad has correct second and finite third derivatives, also where RTMB's pgamma has not", {
  # RTMB's pgamma has non-finite higher derivatives in x for x < 1, and at x = 1 exactly
  f <- function(p) pgamma_ad(p[1], p[2])
  for (x in c(0.05, 0.5, 1, 1.5, 2, 3, 20)) {
    for (a in c(0.3, 2, 10)) {
      p <- c(x, a)
      G <- RTMB::MakeTape(f, p)$jacfun()
      H <- G$jacobian(p)
      h <- 1e-5
      H_fd <- sapply(1:2, function(j) {
        e <- replace(c(0, 0), j, h)
        (as.vector(G(p + e)) - as.vector(G(p - e))) / (2 * h)
      })
      expect_equal(H, H_fd, tolerance = 1e-6, ignore_attr = TRUE)
      expect_true(all(is.finite(G$jacfun()$jacobian(p))))
    }
  }
  # at x = 0 the result is 0, and all derivatives stay finite
  F3 <- RTMB::MakeTape(f, c(0, 0.5))$jacfun()$jacfun()
  expect_true(all(is.finite(F3$jacobian(c(0, 0.5)))))
})

test_that("the atomic tape is reused across independent tapes", {
  F1 <- RTMB::MakeTape(function(p) pgamma_ad(1.3, p), 2)
  F2 <- RTMB::MakeTape(function(p) sum(pgamma_ad(c(0.4, 5), p)), 3)
  expect_equal(F1(2.5), stats::pgamma(1.3, 2.5), tolerance = 1e-14)
  expect_equal(F2(3.5), sum(stats::pgamma(c(0.4, 5), 3.5)), tolerance = 1e-14)
})
