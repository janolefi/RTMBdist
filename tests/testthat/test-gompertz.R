# Tests for the Gompertz distribution
# Wikipedia parameterisation: eta (shape) > 0, b (rate) > 0
# Support: x >= 0

test_that("gompertz passes standard distribution checks (eta=1, b=1)", {
  check_continuous_dist(
    dfun  = dgompertz,
    pfun  = pgompertz,
    qfun  = qgompertz,
    xs    = c(0.1, 0.5, 1, 2, 3),
    lower = 0, upper = Inf,
    eta = 1, b = 1
  )
})

test_that("gompertz passes standard distribution checks (eta=2, b=0.5)", {
  check_continuous_dist(
    dfun  = dgompertz,
    pfun  = pgompertz,
    qfun  = qgompertz,
    xs    = c(0.2, 0.8, 1.5, 3, 5),
    lower = 0, upper = Inf,
    eta = 2, b = 0.5
  )
})

test_that("gompertz AD gradient has no NaN", {
  check_ad_gradient(dgompertz, rgompertz, eta = 1, b = 1)
})
