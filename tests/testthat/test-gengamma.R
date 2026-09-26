# Tests for the generalised gamma distribution (GAMLSS parameterisation)

test_that("gengamma passes standard distribution checks (mu=1, sigma=0.5, nu=1)", {
  check_continuous_dist(
    dfun  = dgengamma,
    pfun  = pgengamma,
    qfun  = qgengamma,
    xs    = c(0.3, 0.6, 1, 2, 4),
    lower = 0, upper = Inf,
    mu = 1, sigma = 0.5, nu = 1
  )
})

test_that("gengamma passes standard distribution checks (mu=2, sigma=1, nu=2)", {
  check_continuous_dist(
    dfun  = dgengamma,
    pfun  = pgengamma,
    qfun  = qgengamma,
    xs    = c(0.5, 1, 2, 4, 7),
    lower = 0, upper = Inf,
    mu = 2, sigma = 1, nu = 2
  )
})

test_that("gengamma AD gradient has no NaN", {
  check_ad_gradient(dgengamma,  rgengamma,  mu = 2, sigma = 1, nu = 2)
})

test_that("pgengamma has correct Hessians and finite third derivatives", {
  q <- c(0.2, 0.5, 0.8, 3, 6)
  check_ad_cdf_hessian(pgengamma, q, mu = 3, sigma = 0.5, nu = 2)
  F3 <- RTMB::MakeTape(function(p) sum(log(pgengamma(q, p[1], p[2], p[3]))), c(3, 0.5, 2))$jacfun()$jacfun()
  expect_true(all(is.finite(F3$jacobian(c(3, 0.5, 2)))))
})
