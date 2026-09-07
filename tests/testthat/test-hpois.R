# Tests for the hurdle (zero-altered) Poisson distribution

test_that("hpois passes discrete distribution checks (lambda=2, zeroprob=0.4)", {
  check_discrete_dist(
    dfun   = dhpois,
    pfun   = phpois,
    xs_int = 0:6,
    lambda = 2, zeroprob = 0.4
  )
})

test_that("hpois passes discrete distribution checks (lambda=8, zeroprob=0.1)", {
  check_discrete_dist(
    dfun   = dhpois,
    pfun   = phpois,
    xs_int = 0:15,
    lambda = 8, zeroprob = 0.1
  )
})

test_that("hpois AD gradient has no NaN", {
  check_ad_gradient(dhpois, rhpois, lambda = 2, zeroprob = 0.4)
})

test_that("zeroprob is exactly the probability of a zero", {
  for (lambda in c(0.2, 1, 5, 20)) {
    for (zeroprob in c(0.05, 0.5, 0.95)) {
      expect_equal(dhpois(0, lambda, zeroprob), zeroprob)
    }
  }
})

test_that("the positive part is the rescaled zero-truncated Poisson", {
  xs <- 1:15
  expect_equal(dhpois(xs, 3, 0.4), 0.6 * dztpois(xs, 3))
  # and the cdf is the cumulative pmf
  expect_equal(phpois(0:20, 3, 0.4), cumsum(dhpois(0:20, 3, 0.4)))
})

test_that("hpois reduces to the Poisson when zeroprob = exp(-lambda)", {
  for (lambda in c(0.5, 2, 6)) {
    expect_equal(dhpois(0:30, lambda, exp(-lambda)), dpois(0:30, lambda))
    expect_equal(phpois(0:30, lambda, exp(-lambda)), ppois(0:30, lambda))
  }
})

test_that("hpois handles the boundary values of zeroprob without NaN", {
  # zeroprob = 0 leaves no mass at zero, i.e. the zero-truncated Poisson
  expect_equal(dhpois(0:20, 2, 0), dztpois(0:20, 2))
  # zeroprob = 1 puts all mass at zero
  expect_equal(dhpois(0:3, 2, 1), c(1, 0, 0, 0))
  for (zeroprob in c(0, 1)) {
    expect_false(any(is.nan(dhpois(c(-1, 0, 1, 5), 2, zeroprob))))
    expect_false(any(is.nan(dhpois(c(-1, 0, 1, 5), 2, zeroprob, log = TRUE))))
  }
})

test_that("hpois has no mass below zero", {
  expect_equal(dhpois(c(-3, -1), 2, 0.4), c(0, 0))
  expect_equal(phpois(c(-3, -1), 2, 0.4), c(0, 0))
})

test_that("rhpois draws have the right zero probability and positive mean", {
  set.seed(42)
  y <- rhpois(2e5, lambda = 2, zeroprob = 0.4)
  expect_true(all(y >= 0 & y == floor(y)))
  expect_equal(mean(y == 0), 0.4, tolerance = 0.02)
  # mean of the zero-truncated Poisson is lambda / (1 - exp(-lambda))
  expect_equal(mean(y[y > 0]), 2 / (1 - exp(-2)), tolerance = 0.02)
})

test_that("hpois recycles a single count against vector parameters", {
  expect_equal(dhpois(2, lambda = c(1, 2, 3), zeroprob = 0.3),
               vapply(c(1, 2, 3), function(l) dhpois(2, l, 0.3), numeric(1)))
  expect_equal(phpois(2, lambda = 2, zeroprob = c(.1, .5, .9)),
               vapply(c(.1, .5, .9), function(z) phpois(2, 2, z), numeric(1)))
})
