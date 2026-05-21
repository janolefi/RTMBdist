# Tests for the truncated t distribution

test_that("trunct passes standard distribution checks (df=5, min=-2, max=2)", {
  check_continuous_dist(
    dfun  = dtrunct,
    pfun  = ptrunct,
    qfun  = qtrunct,
    xs    = c(-1.5, -0.5, 0, 0.5, 1.5),
    lower = -2, upper = 2,
    df = 5, min = -2, max = 2
  )
})

test_that("trunct passes standard distribution checks (df=10, min=0, max=Inf)", {
  check_continuous_dist(
    dfun  = dtrunct,
    pfun  = ptrunct,
    qfun  = qtrunct,
    xs    = c(0.3, 1.0, 2.0, 3.0, 5.0),
    lower = 0, upper = Inf,
    df = 10, min = 0
  )
})

test_that("trunct AD gradient has no NaN", {
  check_ad_gradient(dtrunct,    rtrunct,    df = 5, min = -3, max = 3)
})
