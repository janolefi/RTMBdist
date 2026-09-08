# Tests for the Yule-Simon distribution

test_that("yules passes standard discrete checks", {
  check_discrete_dist(
    dfun        = dyules,
    pfun        = pyules,
    xs_int      = c(1, 2, 5, 20, 100),
    sum_support = 1:5000,
    shape = 2.5
  )
  check_discrete_dist(
    dfun        = dyules,
    pfun        = pyules,
    xs_int      = c(1, 3, 9, 40),
    sum_support = 1:100000, # power-law tail, so the sum converges slowly
    shape = 1.5
  )
})

test_that("yules AD gradient has no NaN", {
  check_ad_gradient(dyules, ryules, shape = 2.5)
  check_ad_gradient(dyules, ryules, shape = 0.8)
})

test_that("dyules and pyules match VGAM", {
  # reference values from VGAM::dyules / VGAM::pyules, baked in so the test
  # needs no dependency on VGAM
  g <- expand.grid(x = c(1, 2, 3, 5, 10, 50), shape = c(0.3, 1, 2.5, 6))
  logd <- c(
    -1.46633706879343, -2.29924619172853, -2.80002147964102,
    -3.44143667311061, -4.32450081512627, -6.40165682214712,
    -0.693147180559945, -1.79175946922805, -2.484906649788,
    -3.40119738166216, -4.70048036579242, -7.84384863815247,
    -0.336472236621213, -1.84054963339749, -2.85215054507597,
    -4.25394909273182, -6.340984341942, -11.6606155237031,
    -0.154150679827258, -2.23359222150709, -3.73766961828337,
    -5.95324333428778, -9.49902194476105, -19.4160425115798)
  cdf <- c(
    0.230769230769231, 0.331103678929766, 0.391912435390696,
    0.466355801132686, 0.558661408878205, 0.723532225406169, 0.5,
    0.666666666666667, 0.75, 0.833333333333333, 0.909090909090909,
    0.980392156862745, 0.714285714285714, 0.873015873015873,
    0.930735930735931, 0.971583971583972, 0.992949734341577,
    0.999827460308127, 0.857142857142857, 0.964285714285714,
    0.988095238095238, 0.997835497835498, 0.999875124875125,
    0.999999969200857)
  expect_equal(dyules(g$x, g$shape, log = TRUE), logd)
  expect_equal(pyules(g$x, g$shape), cdf)
})

test_that("yules support starts at one", {
  expect_equal(dyules(c(-2, 0, 0.5), 2), c(0, 0, 0))
  expect_equal(dyules(c(-2, 0), 2, log = TRUE), c(-Inf, -Inf))
  expect_equal(pyules(c(-2, 0, 0.5), 2), c(0, 0, 0))
  expect_equal(pyules(1, 2), dyules(1, 2)) # the whole mass below 2 sits at 1
  expect_true(dyules(1, 2) > 0)
})

test_that("pyules is the cumulative sum of dyules", {
  expect_equal(pyules(1:60, 2.5), cumsum(dyules(1:60, 2.5)))
  expect_equal(pyules(1:60, 0.7), cumsum(dyules(1:60, 0.7)))
})

test_that("yules is the beta-negative binomial special case it claims to be", {
  k <- 1:60
  expect_equal(dyules(k, 2.5), dbnbinom(k - 1, size = 1, shape1 = 2.5, shape2 = 1))
  # and the Waring with sigma = mu, shifted to start at one
  expect_equal(dyules(k, (3 + 1) / 3), dwaring(k - 1, 3, 3))
})

test_that("yules moments match the closed forms", {
  k <- 1:3000000
  for (rho in c(3, 5)) {
    d <- dyules(k, rho); m <- sum(k * d)
    expect_equal(m, rho / (rho - 1), tolerance = 1e-6)
    expect_equal(sum((k - m)^2 * d), rho^2 / ((rho - 1)^2 * (rho - 2)), tolerance = 1e-5)
  }
})

test_that("ryules draws follow the mass function", {
  set.seed(3)
  x <- ryules(2e4, 2.5)
  expect_true(all(x == floor(x)) && min(x) >= 1)
  e <- c(dyules(1:15, 2.5), 1 - sum(dyules(1:15, 2.5)))
  o <- as.vector(table(factor(pmin(x, 16), levels = 1:16)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
})

test_that("dyules rejects a non-positive shape", {
  expect_error(dyules(1, 0), "shape")
  expect_error(pyules(1, -1), "shape")
})
