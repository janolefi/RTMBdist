# Tests for the hurdle (zero-altered) negative binomial distribution

test_that("hnbinom passes discrete distribution checks (size=2, prob=0.4, zeroprob=0.35)", {
  check_discrete_dist(dfun = dhnbinom, pfun = phnbinom, xs_int = 0:12,
                      sum_support = 0:2000, size = 2, prob = 0.4, zeroprob = 0.35)
})

test_that("hnbinom AD gradient has no NaN", {
  check_ad_gradient(dhnbinom, rhnbinom, size = 2, prob = 0.4, zeroprob = 0.35)
})

test_that("zeroprob is exactly the probability of a zero", {
  for (size in c(0.5, 2, 10)) for (prob in c(0.1, 0.5, 0.9)) for (zp in c(0.05, 0.5, 0.95)) {
    expect_equal(dhnbinom(0, size, prob, zp), zp)
  }
})

test_that("the positive part is the rescaled zero-truncated negative binomial", {
  expect_equal(dhnbinom(1:40, 2, 0.4, 0.35), 0.65 * dztnbinom(1:40, 2, 0.4))
  expect_equal(phnbinom(0:40, 2, 0.4, 0.35), cumsum(dhnbinom(0:40, 2, 0.4, 0.35)))
})

test_that("hnbinom reduces to the negative binomial when zeroprob = P(X = 0)", {
  p0 <- dnbinom(0, size = 2, prob = 0.4)
  expect_equal(dhnbinom(0:60, 2, 0.4, p0), dnbinom(0:60, size = 2, prob = 0.4))
  expect_equal(phnbinom(0:60, 2, 0.4, p0), pnbinom(0:60, size = 2, prob = 0.4))
})

test_that("hnbinom handles boundary zeroprob and values outside the support", {
  expect_equal(dhnbinom(0:40, 2, 0.4, 0), dztnbinom(0:40, 2, 0.4))
  expect_equal(dhnbinom(0:2, 2, 0.4, 1), c(1, 0, 0))
  expect_equal(dhnbinom(c(-3, -1), 2, 0.4, 0.35), c(0, 0))
  expect_equal(phnbinom(c(-3, -1), 2, 0.4, 0.35), c(0, 0))
  for (zp in c(0, 1)) expect_false(any(is.nan(dhnbinom(c(-1, 0, 1, 5), 2, 0.4, zp, log = TRUE))))
})

test_that("pztnbinom returns zero below its support rather than NaN", {
  # regression test: pbeta() is NaN for a negative second shape parameter
  expect_equal(pztnbinom(c(-5, -3, -2, -1), 2, 0.4), c(0, 0, 0, 0))
  expect_equal(pztnbinom2(c(-5, -3, -2, -1), 4, 2), c(0, 0, 0, 0))
})

test_that("rhnbinom draws match the intended zero probability and conditional mean", {
  set.seed(42)
  y <- rhnbinom(2e5, 2, 0.4, 0.35)
  expect_true(all(y >= 0 & y == floor(y)))
  expect_equal(mean(y == 0), 0.35, tolerance = 0.02)
  expect_equal(mean(y[y > 0]), sum((1:600) * dztnbinom(1:600, 2, 0.4)), tolerance = 0.05)
})
