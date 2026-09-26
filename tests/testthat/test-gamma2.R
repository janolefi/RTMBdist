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

test_that("pgamma2 has correct Hessians and finite third derivatives, also at q = scale", {
  # mean 2 and sd 1 give scale 0.5, so q = 0.5 is where RTMB's pgamma has a NaN Hessian
  q <- c(0.1, 0.5, 0.8, 2, 4)
  check_ad_cdf_hessian(pgamma2, q, mean = 2, sd = 1)
  F3 <- RTMB::MakeTape(function(p) sum(log(pgamma2(q, p[1], p[2]))), c(2, 1))$jacfun()$jacfun()
  expect_true(all(is.finite(F3$jacobian(c(2, 1)))))
})
