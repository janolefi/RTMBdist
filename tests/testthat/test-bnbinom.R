# Tests for the beta-negative binomial distribution

test_that("bnbinom passes standard discrete checks", {
  check_discrete_dist(
    dfun        = dbnbinom,
    pfun        = NULL, # no closed-form distribution function
    xs_int      = c(0, 1, 3, 8, 25),
    sum_support = 0:2000,
    size = 3, shape1 = 4.5, shape2 = 2
  )
  check_discrete_dist(
    dfun        = dbnbinom,
    pfun        = NULL,
    xs_int      = c(0, 2, 6, 15, 40),
    sum_support = 0:20000, # heavier tail, so a wider support is needed
    size = 0.5, shape1 = 2.5, shape2 = 1.5
  )
})

test_that("bnbinom AD gradient has no NaN", {
  check_ad_gradient(dbnbinom, rbnbinom, size = 3, shape1 = 4.5, shape2 = 2)
  check_ad_gradient(dbnbinom, rbnbinom, size = 0.7, shape1 = 2, shape2 = 3)
})

test_that("bnbinom mass is zero outside the support", {
  expect_equal(dbnbinom(c(-3, -1, -0.5), 3, 2, 2), c(0, 0, 0))
  expect_equal(dbnbinom(c(-3, -0.5), 3, 2, 2, log = TRUE), c(-Inf, -Inf))
  expect_true(dbnbinom(0, 3, 2, 2) > 0)
})

test_that("bnbinom collapses to the negative binomial as shape1 grows", {
  # X | p ~ NB(size, p) with p ~ Beta(a, b); as a, b -> Inf with a/(a+b) fixed
  # the prior degenerates at that value
  k <- 0:40
  expect_equal(dbnbinom(k, 3, 1e8, 1e8 * (1 / 0.4 - 1)), stats::dnbinom(k, 3, 0.4),
               tolerance = 1e-7)
})

test_that("bnbinom reproduces the mean and variance of its own draws", {
  size <- 3; a <- 4.5; b <- 2
  k <- 0:200000
  d <- dbnbinom(k, size, a, b)
  expect_equal(sum(k * d), size * b / (a - 1), tolerance = 1e-6)
  expect_equal(sum((k - sum(k * d))^2 * d),
               size * b * (size + a - 1) * (b + a - 1) / ((a - 2) * (a - 1)^2),
               tolerance = 1e-4)
})

test_that("rbnbinom draws follow the mass function", {
  set.seed(1)
  x <- rbnbinom(2e4, 3, 4.5, 2)
  expect_true(all(x == floor(x)) && min(x) >= 0)
  e <- c(dbnbinom(0:14, 3, 4.5, 2), 1 - sum(dbnbinom(0:14, 3, 4.5, 2)))
  o <- as.vector(table(factor(pmin(x, 15), levels = 0:15)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
})

test_that("bnbinom refuses OSA residuals and recycles its arguments", {
  expect_length(dbnbinom(0:5, 3, c(2, 4, 6), 2), 6)
  expect_length(rbnbinom(5, 3, c(2, 4), 2), 5)
  expect_error(
    dbnbinom(structure(1, class = "osa"), 3, 2, 2),
    "does not support OSA"
  )
})

test_that("dbnbinom rejects non-positive parameters", {
  expect_error(dbnbinom(1, 0, 2, 2), "size")
  expect_error(dbnbinom(1, 3, -1, 2), "shape1")
  expect_error(dbnbinom(1, 3, 2, 0), "shape")
})
