# Tests for the power generalised Weibull distribution

test_that("pgweibull passes standard distribution checks (scale=1, shape=1, powershape=1)", {
  check_continuous_dist(
    dfun  = dpgweibull,
    pfun  = ppgweibull,
    qfun  = qpgweibull,
    xs    = c(0.2, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    scale = 1, shape = 1, powershape = 1
  )
})

test_that("pgweibull passes standard distribution checks (scale=2, shape=2, powershape=3)", {
  check_continuous_dist(
    dfun  = dpgweibull,
    pfun  = ppgweibull,
    qfun  = qpgweibull,
    xs    = c(0.8, 2.1, 3.9, 7.1, 15.0),
    lower = 0, upper = Inf,
    scale = 2, shape = 2, powershape = 3
  )
})

test_that("pgweibull AD gradient has no NaN", {
  check_ad_gradient(dpgweibull, rpgweibull, scale = 2, shape = 2, powershape = 3)
})
