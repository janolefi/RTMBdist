# Tests for the inverse gamma distribution
# If X ~ Gamma(shape, scale), then 1/X ~ InvGamma(shape, rate=1/scale)

test_that("invgamma passes standard distribution checks (shape=2, rate=1)", {
  check_continuous_dist(
    dfun  = dinvgamma,
    pfun  = pinvgamma,
    qfun  = qinvgamma,
    xs    = c(0.2, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    shape = 2, rate = 1
  )
})

test_that("invgamma passes standard distribution checks (shape=3, rate=2)", {
  check_continuous_dist(
    dfun  = dinvgamma,
    pfun  = pinvgamma,
    qfun  = qinvgamma,
    xs    = c(0.1, 0.3, 0.6, 1, 2),
    lower = 0, upper = Inf,
    shape = 3, rate = 2
  )
})

test_that("invgamma AD gradient has no NaN", {
  check_ad_gradient(dinvgamma,  rinvgamma,  shape = 3, rate = 2)
})
