# Tests for the reparameterised hurdle negative binomial distribution

test_that("hnbinom2 passes discrete distribution checks (mu=4, size=2, zeroprob=0.35)", {
  check_discrete_dist(dfun = dhnbinom2, pfun = phnbinom2, xs_int = 0:15,
                      sum_support = 0:2000, mu = 4, size = 2, zeroprob = 0.35)
})

test_that("hnbinom2 AD gradient has no NaN", {
  check_ad_gradient(dhnbinom2, rhnbinom2, mu = 4, size = 2, zeroprob = 0.35)
})

test_that("hnbinom2 is hnbinom with prob = size / (size + mu)", {
  mu <- 4; size <- 2; zp <- 0.35
  prob <- size / (size + mu)
  expect_equal(dhnbinom2(0:40, mu, size, zp), dhnbinom(0:40, size, prob, zp))
  expect_equal(phnbinom2(0:40, mu, size, zp), phnbinom(0:40, size, prob, zp))
})

test_that("zeroprob is exactly the probability of a zero", {
  for (mu in c(0.5, 4, 20)) for (size in c(0.5, 2, 10)) for (zp in c(0.05, 0.5, 0.95)) {
    expect_equal(dhnbinom2(0, mu, size, zp), zp)
  }
})

test_that("the positive part is the rescaled zero-truncated negative binomial", {
  expect_equal(dhnbinom2(1:40, 4, 2, 0.35), 0.65 * dztnbinom2(1:40, 4, 2))
  expect_equal(phnbinom2(0:40, 4, 2, 0.35), cumsum(dhnbinom2(0:40, 4, 2, 0.35)))
})

test_that("hnbinom2 reduces to the negative binomial when zeroprob = P(X = 0)", {
  p0 <- dnbinom(0, size = 2, mu = 4)
  expect_equal(dhnbinom2(0:60, 4, 2, p0), dnbinom(0:60, size = 2, mu = 4))
  expect_equal(phnbinom2(0:60, 4, 2, p0), pnbinom(0:60, size = 2, mu = 4))
})

test_that("hnbinom2 handles boundary zeroprob and values outside the support", {
  expect_equal(dhnbinom2(0:40, 4, 2, 0), dztnbinom2(0:40, 4, 2))
  expect_equal(dhnbinom2(0:2, 4, 2, 1), c(1, 0, 0))
  expect_equal(dhnbinom2(c(-3, -1), 4, 2, 0.35), c(0, 0))
  expect_equal(phnbinom2(c(-3, -1), 4, 2, 0.35), c(0, 0))
  for (zp in c(0, 1)) expect_false(any(is.nan(dhnbinom2(c(-1, 0, 1, 5), 4, 2, zp, log = TRUE))))
})

test_that("rhnbinom2 draws match the intended zero probability and conditional mean", {
  set.seed(42)
  y <- rhnbinom2(2e5, 4, 2, 0.35)
  expect_true(all(y >= 0 & y == floor(y)))
  expect_equal(mean(y == 0), 0.35, tolerance = 0.02)
  expect_equal(mean(y[y > 0]), sum((1:600) * dztnbinom2(1:600, 4, 2)), tolerance = 0.05)
})
