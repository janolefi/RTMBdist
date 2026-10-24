# Tests for the GB2 distribution and its special cases Dagum, Singh-Maddala and Lomax.
# Shapes are not integers of 3 or more in the AD tests, because RTMB's pbeta() has a zero
# derivative in its second shape there (RTMB 2.0).

test_that("gb2 and its special cases pass standard distribution checks", {
  check_continuous_dist(dgb2, pgb2, qgb2, xs = c(0.3, 1, 2, 5), lower = 0,
                        scale = 2, shape1 = 3, shape2 = 1.5, shape3 = 0.8)
  check_continuous_dist(dgb2, pgb2, qgb2, xs = c(0.01, 0.5, 4, 20), lower = 0,
                        scale = 0.7, shape1 = 0.9, shape2 = 0.6, shape3 = 2.2)
  check_continuous_dist(ddagum, pdagum, qdagum, xs = c(0.3, 1, 2, 5), lower = 0,
                        scale = 2, shape1 = 3, shape2 = 1.5)
  check_continuous_dist(dsinmad, psinmad, qsinmad, xs = c(0.3, 1, 2, 5), lower = 0,
                        scale = 2, shape1 = 3, shape3 = 0.8)
  check_continuous_dist(dlomax, plomax, qlomax, xs = c(0.3, 1, 2, 5), lower = 0,
                        scale = 2, shape3 = 1.8)
})

test_that("dgb2 and pgb2 match gamlss.dist", {
  # reference values from gamlss.dist::dGB2 / pGB2 with mu = scale, sigma = shape1,
  # nu = shape2, tau = shape3, baked in: gamlss.dist is scheduled for archival
  x <- c(0.05, 0.7, 2, 9, 60, 0.4, 3, 1e3)
  b <- c(0.6, 3, 0.6, 3, 0.6, 3, 1.5, 2)
  a <- c(0.8, 4, 4, 0.8, 0.8, 4, 2.5, 1.2)
  p <- c(0.5, 2.3, 0.5, 2.3, 2.3, 0.5, 1.7, 0.9)
  q <- c(0.7, 3.4, 3.4, 0.7, 3.4, 0.7, 1.3, 2.6)
  logd <- c(0.705969752120814, -8.6250940285255, -15.7094583609617, -3.79940564515269,
            -13.9475097341533, -2.64620642688972, -2.02625746610723, -25.3398446875072)
  cdf <- c(0.280541267759449, 1.3730964302621e-05, 0.999999977712334,
           0.322260324455328, 0.999980036084721, 0.0141875346512904, 0.849214488313934,
           0.999999996829452)
  expect_equal(dgb2(x, b, a, p, q, log = TRUE), logd, tolerance = 1e-12)
  expect_equal(pgb2(x, b, a, p, q), cdf, tolerance = 1e-12)
})

test_that("the special cases match VGAM", {
  # reference values from VGAM::ddagum, dsinmad, dlomax and the p functions, baked in
  x <- c(0.3, 1, 4)
  expect_equal(ddagum(x, 2, 3, 1.5), c(0.002916377, 0.1481481, 0.1047566), tolerance = 1e-6)
  expect_equal(dsinmad(x, 2, 3, 0.8), c(0.02683675, 0.2426871, 0.09196122), tolerance = 1e-6)
  expect_equal(dlomax(x, 2, 0.8), c(0.3110311, 0.192795, 0.05536582), tolerance = 1e-6)
  expect_equal(pdagum(x, 2, 3, 1.5), c(0.0001950813, 0.03703704, 0.8380525), tolerance = 1e-6)
  expect_equal(psinmad(x, 2, 3, 0.8), c(0.002691825, 0.08992333, 0.8275727), tolerance = 1e-6)
  expect_equal(plomax(x, 2, 0.8), c(0.1057855, 0.2770188, 0.5847564), tolerance = 1e-6)
})

test_that("gb2 contains the log-logistic and beta prime distributions", {
  x <- c(0.1, 1, 3, 20)
  expect_equal(dgb2(x, 2, 3, 1, 1), dllogis(x, alpha = 2, beta = 3))
  expect_equal(pgb2(x, 2, 3, 1, 1), pllogis(x, alpha = 2, beta = 3))
  expect_equal(dgb2(x, 1, 1, 2.5, 1.5), dbetaprime(x, 2.5, 1.5))
})

test_that("both tails are accurate far out", {
  # closed forms of the Lomax distribution: P(X > x) = (1 + x / b)^(-q)
  expect_equal(plomax(1e8, 2, 0.8, lower.tail = FALSE), (1 + 5e7)^-0.8)
  expect_equal(plomax(1e-12, 2, 0.8), -expm1(-0.8 * log1p(5e-13)))
  expect_equal(plomax(1e300, 2, 0.8, lower.tail = FALSE, log.p = TRUE),
               -0.8 * log(5e299))
})

test_that("gb2 is zero below its support and handles x = 0", {
  expect_equal(dgb2(c(-1, -1e-10), 2, 3, 1.5, 0.8), c(0, 0))
  expect_equal(pgb2(c(-1, 0), 2, 3, 1.5, 0.8), c(0, 0))
  expect_equal(dlomax(0, 2, 0.8), 0.8 / 2) # a p = 1: the density is finite at 0
})

test_that("gb2 AD derivatives are correct", {
  check_ad_gradient(dgb2, rgb2, scale = 2, shape1 = 3, shape2 = 1.5, shape3 = 0.8)
  check_ad_gradient(dlomax, rlomax, scale = 2, shape3 = 1.8)
  xs <- c(0.3, 1, 2, 5)
  check_ad_cdf(pgb2, dgb2, xs, scale = 2, shape1 = 2.7, shape2 = 1.5, shape3 = 0.8)
  check_ad_cdf(pgb2, dgb2, xs, scale = 2, shape1 = 2.7, shape2 = 1.5, shape3 = 0.8,
               .args = list(lower.tail = FALSE))
  check_ad_cdf_hessian(pgb2, xs, scale = 2, shape1 = 2.7, shape2 = 1.5, shape3 = 0.8)
  # qgb2 is differentiable as well
  F <- RTMB::MakeTape(function(s) qgb2(c(0.1, 0.5, 0.9), 2, s, 1.5, 0.8), 2.7)
  h <- 1e-6
  expect_equal(as.vector(F$jacobian(2.7)),
               (qgb2(c(0.1, 0.5, 0.9), 2, 2.7 + h, 1.5, 0.8) -
                  qgb2(c(0.1, 0.5, 0.9), 2, 2.7 - h, 1.5, 0.8)) / (2 * h), tolerance = 1e-6)
})

test_that("gb2 and its special cases support OSA residuals", {
  set.seed(1)
  check_osa_cdf(dgb2, pgb2, rgb2(20, 2, 3, 1.5, 0.8),
                scale = 2, shape1 = 3, shape2 = 1.5, shape3 = 0.8)
  check_osa_cdf(dlomax, plomax, rlomax(20, 2, 1.8), scale = 2, shape3 = 1.8)
})

test_that("rgb2 draws from the gb2", {
  set.seed(1)
  x <- rgb2(20000, 2, 3, 1.5, 0.8)
  expect_gt(ks.test(x, pgb2, 2, 3, 1.5, 0.8)$p.value, 0.01)
  expect_length(rgb2(4, scale = 1:4, 3, 1.5, 0.8), 4)
})

test_that("gb2 recycles and rejects non-positive parameters", {
  expect_equal(dgb2(1, 2, c(1, 3), 1.5, 0.8), c(dgb2(1, 2, 1, 1.5, 0.8), dgb2(1, 2, 3, 1.5, 0.8)))
  expect_error(dgb2(1, 0, 3, 1.5, 0.8), "scale")
  expect_error(pgb2(1, 2, -3, 1.5, 0.8), "shape1")
  expect_error(qdagum(0.5, 2, 3, 0), "shape2")
  expect_error(rlomax(1, 2, 0), "shape3")
})
