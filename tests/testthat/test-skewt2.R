# Tests for the moment-parameterised skew t distribution

test_that("skewt2 passes standard distribution checks (mean=0, sd=1, skew=0, df=5)", {
  check_continuous_dist(
    dfun  = dskewt2,
    pfun  = pskewt2,
    qfun  = qskewt2,
    xs    = c(-2, -1, 0, 1, 2),
    mean = 0, sd = 1, skew = 0, df = 5
  )
})

test_that("skewt2 passes standard distribution checks (mean=1, sd=2, skew=2, df=5)", {
  check_continuous_dist(
    dfun  = dskewt2,
    pfun  = pskewt2,
    qfun  = qskewt2,
    xs    = c(-1.5, -0.3, 0.7, 1.9, 4.5),
    mean = 1, sd = 2, skew = 2, df = 5
  )
})

test_that("skewt2 AD gradient has no NaN", {
  check_ad_gradient(dskewt2,    rskewt2,    mean = 0, sd = 1, skew = 2, df = 5)
})
