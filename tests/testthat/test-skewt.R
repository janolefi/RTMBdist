# Tests for the skew t distribution

test_that("skewt passes standard distribution checks (mu=0, sigma=1, skew=0, df=5)", {
  check_continuous_dist(
    dfun  = dskewt,
    pfun  = pskewt,
    qfun  = qskewt,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, skew = 0, df = 5
  )
})

test_that("skewt passes standard distribution checks (mu=0, sigma=1, skew=2, df=5)", {
  check_continuous_dist(
    dfun  = dskewt,
    pfun  = pskewt,
    qfun  = qskewt,
    xs    = c(-0.4, 0.2, 0.7, 1.3, 2.6),
    mu = 0, sigma = 1, skew = 2, df = 5
  )
})

test_that("skewt AD gradient has no NaN", {
  check_ad_gradient(dskewt,     rskewt,     mu = 0, sigma = 1, skew = 2, df = 5)
})
