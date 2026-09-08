# Tests for the zero-modified beta-binomial distributions.
# betabinom itself has no cdf, so neither do these; pfun = NULL throughout.

test_that("zi/zt/h beta-binomial pass discrete distribution checks", {
  check_discrete_dist(dfun = dzibetabinom, pfun = NULL, xs_int = 0:10,
                      sum_support = 0:10, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.3)
  check_discrete_dist(dfun = dztbetabinom, pfun = NULL, xs_int = 1:10,
                      sum_support = 0:10, size = 10, shape1 = 2, shape2 = 3)
  check_discrete_dist(dfun = dhbetabinom, pfun = NULL, xs_int = 0:10,
                      sum_support = 0:10, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
})

test_that("beta-binomial variants have AD gradients without NaN", {
  check_ad_gradient(dzibetabinom, rzibetabinom, size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.3)
  check_ad_gradient(dztbetabinom, rztbetabinom, size = 10, shape1 = 2, shape2 = 3)
  check_ad_gradient(dhbetabinom,  rhbetabinom,  size = 10, shape1 = 2, shape2 = 3, zeroprob = 0.4)
})

test_that("beta-binomial variants satisfy their defining identities", {
  sz <- 10; a <- 2; b <- 3; zp <- 0.35; x <- 1:10
  p0 <- dbetabinom(0, sz, a, b)
  expect_equal(dzibetabinom(0, sz, a, b, zp), zp + (1 - zp) * p0)
  expect_equal(dzibetabinom(x, sz, a, b, zp), (1 - zp) * dbetabinom(x, sz, a, b))
  expect_equal(dztbetabinom(x, sz, a, b),     dbetabinom(x, sz, a, b) / (1 - p0))
  expect_equal(dhbetabinom(0, sz, a, b, zp),  zp)
  expect_equal(dhbetabinom(x, sz, a, b, zp),  (1 - zp) * dztbetabinom(x, sz, a, b))
  # a hurdle with zeroprob = P(X = 0) is the ordinary beta-binomial
  expect_equal(dhbetabinom(0:sz, sz, a, b, p0), dbetabinom(0:sz, sz, a, b))
  # a hurdle with zeroprob = 0 is the zero-truncated beta-binomial
  expect_equal(dhbetabinom(0:sz, sz, a, b, 0), dztbetabinom(0:sz, sz, a, b))
})

test_that("beta-binomial variants have no mass below zero", {
  expect_equal(dzibetabinom(c(-3, -1), 10, 2, 3, 0.3), c(0, 0))
  expect_equal(dztbetabinom(c(-3, -1), 10, 2, 3), c(0, 0))
  expect_equal(dhbetabinom(c(-3, -1), 10, 2, 3, 0.4), c(0, 0))
})

test_that("hbetabinom handles the boundary values of zeroprob without NaN", {
  expect_equal(dhbetabinom(0:3, 10, 2, 3, 1), c(1, 0, 0, 0))
  for (zp in c(0, 1)) {
    expect_false(any(is.nan(dhbetabinom(c(-1, 0, 1, 5), 10, 2, 3, zp))))
    expect_false(any(is.nan(dhbetabinom(c(-1, 0, 1, 5), 10, 2, 3, zp, log = TRUE))))
  }
})

test_that("beta-binomial variant RNGs have the right zero behaviour", {
  set.seed(42)
  p0 <- dbetabinom(0, 10, 2, 3)
  expect_equal(mean(rzibetabinom(1e5, 10, 2, 3, 0.3) == 0), 0.3 + 0.7 * p0, tolerance = 0.02)
  y <- rztbetabinom(1e4, 10, 2, 3)
  expect_true(all(y >= 1 & y <= 10))
  y <- rhbetabinom(1e5, 10, 2, 3, 0.4)
  expect_equal(mean(y == 0), 0.4, tolerance = 0.02)
  expect_equal(mean(y[y > 0]), sum((1:10) * dztbetabinom(1:10, 10, 2, 3)), tolerance = 0.05)
})
