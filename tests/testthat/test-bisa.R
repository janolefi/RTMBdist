# Tests for the Birnbaum-Saunders distribution

test_that("bisa passes standard distribution checks", {
  check_continuous_dist(dbisa, pbisa, qbisa, xs = c(0.5, 1, 2, 4), lower = 0,
                        scale = 2, shape = 0.5)
  check_continuous_dist(dbisa, pbisa, qbisa, xs = c(0.01, 0.3, 3, 30), lower = 0,
                        scale = 1.5, shape = 2.2)
})

test_that("dbisa and pbisa match VGAM and extraDistr", {
  # reference values from VGAM::dbisa / pbisa, baked in; extraDistr::dfatigue agrees
  x <- c(0.05, 0.6, 2, 7, 40)
  b <- c(2, 2, 0.5, 3, 3)
  a <- c(0.5, 1.7, 0.5, 0.2, 3)
  logd <- c(-72.1040739200034, -1.05012129611431, -5.19579498189046, -10.692043600063,
            -5.66591950854042)
  cdf <- c(3.01259162906102e-35, 0.226092870770399, 0.99865010196837,
           0.999993625162539, 0.869890602918604)
  expect_equal(dbisa(x, b, a, log = TRUE), logd, tolerance = 1e-12)
  expect_equal(pbisa(x, b, a), cdf, tolerance = 1e-12)
})

test_that("scale is the median and the mean is scale * (1 + shape^2 / 2)", {
  expect_equal(pbisa(2, 2, 0.7), 0.5)
  expect_equal(qbisa(0.5, 2, 0.7), 2)
  m <- integrate(function(x) x * dbisa(x, 2, 0.7), 0, Inf)$value
  expect_equal(m, 2 * (1 + 0.7^2 / 2), tolerance = 1e-7)
})

test_that("both tails are accurate far out", {
  expect_equal(pbisa(1e-3, 2, 0.5), pnorm(2 * sinh(0.5 * log(5e-4)) / 0.5))
  expect_equal(pbisa(500, 2, 0.5, lower.tail = FALSE), pnorm(-2 * sinh(0.5 * log(250)) / 0.5))
  expect_equal(pbisa(500, 2, 0.5, lower.tail = FALSE, log.p = TRUE),
               pnorm(-2 * sinh(0.5 * log(250)) / 0.5, log.p = TRUE))
  # q(p) stays accurate for small p, where the textbook formula cancels
  expect_equal(pbisa(qbisa(1e-20, 2, 0.5), 2, 0.5), 1e-20)
})

test_that("bisa is zero below its support", {
  expect_equal(dbisa(c(-1, 0), 2, 0.5), c(0, 0))
  expect_equal(pbisa(c(-1, 0), 2, 0.5), c(0, 0))
})

test_that("bisa AD derivatives are correct", {
  check_ad_gradient(dbisa, rbisa, scale = 2, shape = 0.5)
  xs <- c(0.5, 1, 2, 4)
  check_ad_cdf(pbisa, dbisa, xs, scale = 2, shape = 0.5)
  check_ad_cdf_hessian(pbisa, xs, scale = 2, shape = 0.5)
  F <- RTMB::MakeTape(function(s) qbisa(c(0.1, 0.5, 0.9), 2, s), 0.5)
  h <- 1e-6
  expect_equal(as.vector(F$jacobian(0.5)),
               (qbisa(c(0.1, 0.5, 0.9), 2, 0.5 + h) - qbisa(c(0.1, 0.5, 0.9), 2, 0.5 - h)) / (2 * h),
               tolerance = 1e-6)
})

test_that("bisa supports OSA residuals", {
  set.seed(1)
  check_osa_cdf(dbisa, pbisa, rbisa(20, 2, 0.5), scale = 2, shape = 0.5)
})

test_that("rbisa draws from the bisa", {
  set.seed(1)
  x <- rbisa(20000, 2, 0.5)
  expect_gt(ks.test(x, pbisa, 2, 0.5)$p.value, 0.01)
  expect_length(rbisa(4, scale = 1:4, 0.5), 4)
})

test_that("bisa recycles and rejects non-positive parameters", {
  expect_equal(dbisa(1, 2, c(0.5, 2)), c(dbisa(1, 2, 0.5), dbisa(1, 2, 2)))
  expect_error(dbisa(1, 0, 0.5), "scale")
  expect_error(pbisa(1, 2, 0), "shape")
  expect_error(rbisa(1, 2, -1), "shape")
})
