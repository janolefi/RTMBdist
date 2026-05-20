# Tests for the location-scale t distribution

test_that("t2 passes standard distribution checks (mu=0, sigma=1, df=5)", {
  check_continuous_dist(
    dfun  = dt2,
    pfun  = pt2,
    qfun  = qt2,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, df = 5
  )
})

test_that("t2 passes standard distribution checks (mu=2, sigma=3, df=10)", {
  check_continuous_dist(
    dfun  = dt2,
    pfun  = pt2,
    qfun  = qt2,
    xs    = c(-4, 0, 2, 5, 8),
    mu = 2, sigma = 3, df = 10
  )
})
