# Tests for the truncated normal distribution

test_that("truncnorm passes standard distribution checks (mean=0, sd=1, min=-2, max=2)", {
  check_continuous_dist(
    dfun  = dtruncnorm,
    pfun  = ptruncnorm,
    qfun  = qtruncnorm,
    xs    = c(-1.5, -0.5, 0, 0.5, 1.5),
    lower = -2, upper = 2,
    mean = 0, sd = 1, min = -2, max = 2
  )
})

test_that("truncnorm passes standard distribution checks (mean=2, sd=1, min=0, max=Inf)", {
  check_continuous_dist(
    dfun  = dtruncnorm,
    pfun  = ptruncnorm,
    qfun  = qtruncnorm,
    xs    = c(0.3, 0.8, 1.5, 2.5, 4.0),
    lower = 0, upper = Inf,
    mean = 2, sd = 1, min = 0
  )
})

test_that("truncnorm AD gradient has no NaN", {
  check_ad_gradient(dtruncnorm, rtruncnorm, mean = 0, sd = 1, min = -3, max = 3)
})
