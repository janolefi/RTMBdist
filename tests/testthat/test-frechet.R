# Tests for the Frechet distribution

test_that("frechet passes standard distribution checks", {
  check_continuous_dist(
    dfun  = dfrechet,
    pfun  = pfrechet,
    qfun  = qfrechet,
    xs    = c(0.3, 0.7, 1.2, 2.5, 8),
    lower = 0, upper = Inf,
    mu = 0, sigma = 1, alpha = 3
  )
})

test_that("frechet passes standard distribution checks (shifted, heavy tail)", {
  check_continuous_dist(
    dfun  = dfrechet,
    pfun  = pfrechet,
    qfun  = qfrechet,
    xs    = c(2.2, 3, 5, 12, 60),
    lower = 2, upper = Inf,
    mu = 2, sigma = 1.5, alpha = 0.8
  )
})

test_that("frechet AD gradient has no NaN", {
  check_ad_gradient(dfrechet, rfrechet, mu = 0, sigma = 1, alpha = 3)
  check_ad_gradient(dfrechet, rfrechet, mu = 0, sigma = 2, alpha = 0.7)
})

test_that("frechet is the gev with shape 1 / alpha", {
  # Frechet(mu, sigma, alpha) is GEV(mu + sigma, sigma / alpha, 1 / alpha)
  xs <- c(0.2, 1, 3, 20)
  expect_equal(dfrechet(xs, 0, 1, 3), dgev(xs, mu = 1, sigma = 1 / 3, xi = 1 / 3))
  expect_equal(pfrechet(xs, 0, 1, 3), pgev(xs, mu = 1, sigma = 1 / 3, xi = 1 / 3))
  xs2 <- c(2.5, 4, 10)
  expect_equal(dfrechet(xs2, 2, 1.5, 0.8), dgev(xs2, mu = 3.5, sigma = 1.5 / 0.8, xi = 1 / 0.8))
})

test_that("frechet is zero outside its support and the end points are right", {
  expect_equal(dfrechet(c(-1, 0), 0, 1, 3), c(0, 0))
  expect_equal(pfrechet(c(-1, 0), 0, 1, 3), c(0, 0))
  expect_equal(dfrechet(c(1, 2), 2, 1, 3), c(0, 0)) # support is x > mu
  expect_equal(qfrechet(0, 2, 1, 3), 2)
  expect_equal(qfrechet(1, 2, 1, 3), Inf)
})

test_that("frechet gradients stay finite when observations fall outside the support", {
  set.seed(2)
  x <- rfrechet(50, 0, 1, 3)
  nll <- function(p) -sum(dfrechet(x, p[1], exp(p[2]), exp(p[3]), log = TRUE))
  F <- RTMB::MakeTape(nll, c(0, 0, log(3)))

  expect_equal(F(c(1, 0, log(3))), Inf) # location pushed above part of the data
  expect_false(any(is.nan(as.vector(F$jacobian(c(1, 0, log(3)))))))
})

test_that("qfrechet honours lower.tail and log.p", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qfrechet(p, 1, 2, 3), qfrechet(1 - p, 1, 2, 3, lower.tail = FALSE))
  expect_equal(qfrechet(p, 1, 2, 3), qfrechet(log(p), 1, 2, 3, log.p = TRUE))
})

test_that("dfrechet rejects non-positive scale and shape", {
  expect_error(dfrechet(1, 0, -1, 1), "sigma")
  expect_error(dfrechet(1, 0, 1, -2), "alpha")
  expect_error(pfrechet(1, 0, 1, 0), "alpha")
})
