# Tests for the hurdle (zero-altered) binomial distribution

test_that("hbinom passes discrete distribution checks (size=15, prob=0.3, zeroprob=0.4)", {
  check_discrete_dist(dfun = dhbinom, pfun = phbinom, xs_int = 0:10,
                      sum_support = 0:15, size = 15, prob = 0.3, zeroprob = 0.4)
})

test_that("hbinom passes discrete distribution checks (size=5, prob=0.7, zeroprob=0.1)", {
  check_discrete_dist(dfun = dhbinom, pfun = phbinom, xs_int = 0:5,
                      sum_support = 0:5, size = 5, prob = 0.7, zeroprob = 0.1)
})

test_that("hbinom AD gradient has no NaN", {
  check_ad_gradient(dhbinom, rhbinom, size = 15, prob = 0.3, zeroprob = 0.4)
})

test_that("zeroprob is exactly the probability of a zero", {
  for (size in c(1, 5, 20)) for (prob in c(0.1, 0.5, 0.9)) for (zp in c(0.05, 0.5, 0.95)) {
    expect_equal(dhbinom(0, size, prob, zp), zp)
  }
})

test_that("the positive part is the rescaled zero-truncated binomial", {
  expect_equal(dhbinom(1:15, 15, 0.3, 0.4), 0.6 * dztbinom(1:15, 15, 0.3))
  expect_equal(phbinom(0:15, 15, 0.3, 0.4), cumsum(dhbinom(0:15, 15, 0.3, 0.4)))
})

test_that("hbinom reduces to the binomial when zeroprob = P(X = 0)", {
  for (prob in c(0.2, 0.5)) {
    p0 <- dbinom(0, 15, prob)
    expect_equal(dhbinom(0:15, 15, prob, p0), dbinom(0:15, 15, prob))
    expect_equal(phbinom(0:15, 15, prob, p0), pbinom(0:15, 15, prob))
  }
})

test_that("hbinom handles boundary zeroprob and values outside the support", {
  expect_equal(dhbinom(0:15, 15, 0.3, 0), dztbinom(0:15, 15, 0.3))
  expect_equal(dhbinom(0:2, 15, 0.3, 1), c(1, 0, 0))
  expect_equal(dhbinom(c(-3, -1), 15, 0.3, 0.4), c(0, 0))
  expect_equal(phbinom(c(-3, -1), 15, 0.3, 0.4), c(0, 0))
  for (zp in c(0, 1)) expect_false(any(is.nan(dhbinom(c(-1, 0, 1, 5), 15, 0.3, zp, log = TRUE))))
})

test_that("rhbinom draws match the intended zero probability and conditional mean", {
  set.seed(42)
  y <- rhbinom(2e5, 15, 0.3, 0.4)
  expect_true(all(y >= 0 & y <= 15 & y == floor(y)))
  expect_equal(mean(y == 0), 0.4, tolerance = 0.02)
  expect_equal(mean(y[y > 0]), sum((1:15) * dztbinom(1:15, 15, 0.3)), tolerance = 0.02)
})
