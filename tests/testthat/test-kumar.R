# Tests for the Kumaraswamy distribution
# CDF: pkumar(q, a, b) = 1 - (1 - q^a)^b

test_that("kumar passes standard distribution checks (a=2, b=3)", {
  check_continuous_dist(
    dfun  = dkumar,
    pfun  = pkumar,
    qfun  = qkumar,
    xs    = c(0.1, 0.3, 0.5, 0.7, 0.9),
    lower = 0, upper = 1,
    a = 2, b = 3
  )
})

test_that("kumar passes standard distribution checks (a=0.5, b=2)", {
  check_continuous_dist(
    dfun  = dkumar,
    pfun  = pkumar,
    qfun  = qkumar,
    xs    = c(0.05, 0.2, 0.4, 0.6, 0.85),
    lower = 0, upper = 1,
    a = 0.5, b = 2
  )
})

test_that("kumar AD gradient has no NaN", {
  check_ad_gradient(dkumar, rkumar, a = 2, b = 3)
})
