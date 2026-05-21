# Tests for the inverse chi-squared distribution
# If X ~ Chi2(nu), then 1/X ~ InvChisq(nu)
# scale defaults to 1/df (standard parameterisation)

test_that("invchisq passes standard distribution checks (df=5)", {
  check_continuous_dist(
    dfun  = dinvchisq,
    pfun  = pinvchisq,
    qfun  = qinvchisq,
    xs    = c(0.05, 0.1, 0.2, 0.5, 1),
    lower = 0, upper = Inf,
    df = 5
  )
})

test_that("invchisq passes standard distribution checks (df=10, scale=0.5)", {
  check_continuous_dist(
    dfun  = dinvchisq,
    pfun  = pinvchisq,
    qfun  = qinvchisq,
    xs    = c(0.1, 0.2, 0.4, 0.8, 1.5),
    lower = 0, upper = Inf,
    df = 10, scale = 0.5
  )
})

test_that("invchisq AD gradient has no NaN", {
  check_ad_gradient(dinvchisq,  rinvchisq,  df = 10, scale = 0.5)
})
