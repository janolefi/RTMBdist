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

test_that("pinvchisq has correct Hessians and finite third derivatives", {
  q <- c(0.05, 0.2, 1, 3, 8)
  check_ad_cdf_hessian(pinvchisq, q, df = 5, scale = 0.5)
  F3 <- RTMB::MakeTape(function(p) sum(log(pinvchisq(q, p[1], p[2]))), c(5, 0.5))$jacfun()$jacfun()
  expect_true(all(is.finite(F3$jacobian(c(5, 0.5)))))
})
