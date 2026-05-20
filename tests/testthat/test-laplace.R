# Tests for the Laplace distribution

test_that("laplace passes standard distribution checks (mu=0, b=1)", {
  check_continuous_dist(
    dfun  = dlaplace,
    pfun  = plaplace,
    qfun  = qlaplace,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, b = 1
  )
})

test_that("laplace passes standard distribution checks (mu=2, b=3)", {
  check_continuous_dist(
    dfun  = dlaplace,
    pfun  = plaplace,
    qfun  = qlaplace,
    xs    = c(-2, 0, 2, 4, 6),
    mu = 2, b = 3
  )
})
