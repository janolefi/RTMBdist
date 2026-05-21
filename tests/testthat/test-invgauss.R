# Tests for the inverse Gaussian distribution

test_that("invgauss passes standard distribution checks (mean=1, shape=1)", {
  check_continuous_dist(
    dfun  = dinvgauss,
    pfun  = pinvgauss,
    qfun  = qinvgauss,
    xs    = c(0.3, 0.6, 1, 2, 4),
    lower = 0, upper = Inf,
    mean = 1, shape = 1
  )
})

test_that("invgauss passes standard distribution checks (mean=2, shape=3)", {
  check_continuous_dist(
    dfun  = dinvgauss,
    pfun  = pinvgauss,
    qfun  = qinvgauss,
    xs    = c(0.5, 1, 2, 3, 5),
    lower = 0, upper = Inf,
    mean = 2, shape = 3
  )
})

test_that("invgauss AD gradient has no NaN", {
  check_ad_gradient(dinvgauss, rinvgauss, mean = 2, shape = 3)
})
