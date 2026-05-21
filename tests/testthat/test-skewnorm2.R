# Tests for the reparameterised skew normal distribution (mean/sd parameterisation)

test_that("skewnorm2 passes standard distribution checks (mean=0, sd=1, alpha=0)", {
  check_continuous_dist(
    dfun  = dskewnorm2,
    pfun  = pskewnorm2,
    qfun  = qskewnorm2,
    xs    = c(-2, -1, 0, 1, 2),
    mean = 0, sd = 1, alpha = 0
  )
})

test_that("skewnorm2 passes standard distribution checks (mean=2, sd=1.5, alpha=2)", {
  check_continuous_dist(
    dfun  = dskewnorm2,
    pfun  = pskewnorm2,
    qfun  = qskewnorm2,
    xs    = c(-0.2, 0.9, 1.9, 2.9, 4.7),
    mean = 2, sd = 1.5, alpha = 2
  )
})

test_that("skewnorm2 AD gradient has no NaN", {
  check_ad_gradient(dskewnorm2, rskewnorm2, mean = 0, sd = 1, alpha = 2)
})
