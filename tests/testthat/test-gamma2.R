# Tests for the gamma2 distribution
# (gamma reparameterised by mean and standard deviation)
# shape = mean^2/sd^2, scale = sd^2/mean

test_that("gamma2 passes standard distribution checks (mean=2, sd=1)", {
  # shape = 4, scale = 0.5 — well away from boundary singularities
  check_continuous_dist(
    dfun  = dgamma2,
    pfun  = pgamma2,
    qfun  = qgamma2,
    xs    = c(0.5, 1, 2, 3, 4),
    lower = 0, upper = Inf,
    mean = 2, sd = 1
  )
})

test_that("gamma2 passes standard distribution checks (mean=0.5, sd=0.5)", {
  # shape = 1 (exponential), scale = 0.5
  check_continuous_dist(
    dfun  = dgamma2,
    pfun  = pgamma2,
    qfun  = qgamma2,
    xs    = c(0.1, 0.3, 0.5, 1, 2),
    lower = 0, upper = Inf,
    mean = 0.5, sd = 0.5
  )
})

test_that("gamma2 AD gradient has no NaN", {
  check_ad_gradient(dgamma2, rgamma2, mean = 2, sd = 1)
})

test_that("gamma2 matches dgamma with converted shape and scale", {
  x <- c(0.5, 1, 2, 3, 4)
  mean <- 2; sd <- 1
  shape <- mean^2 / sd^2
  scale <- sd^2 / mean
  expect_equal(
    dgamma2(x, mean = mean, sd = sd),
    dgamma(x, shape = shape, scale = scale),
    tolerance = 1e-10
  )
})
