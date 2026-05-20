# Tests for the power exponential distribution (and its variant powerexp2)
# nu = 2 gives the normal, nu = 1 gives the Laplace

test_that("powerexp passes standard distribution checks (mu=0, sigma=1, nu=2)", {
  check_continuous_dist(
    dfun  = dpowerexp,
    pfun  = ppowerexp,
    qfun  = qpowerexp,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, nu = 2
  )
})

test_that("powerexp passes standard distribution checks (mu=1, sigma=2, nu=1)", {
  check_continuous_dist(
    dfun  = dpowerexp,
    pfun  = ppowerexp,
    qfun  = qpowerexp,
    xs    = c(-3, 0, 1, 3, 6),
    mu = 1, sigma = 2, nu = 1
  )
})

test_that("powerexp2 passes standard distribution checks (mu=0, sigma=1, nu=2)", {
  check_continuous_dist(
    dfun  = dpowerexp2,
    pfun  = ppowerexp2,
    qfun  = qpowerexp2,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, nu = 2
  )
})

test_that("powerexp2 passes standard distribution checks (mu=1, sigma=2, nu=1)", {
  check_continuous_dist(
    dfun  = dpowerexp2,
    pfun  = ppowerexp2,
    qfun  = qpowerexp2,
    xs    = c(-3, 0, 1, 3, 6),
    mu = 1, sigma = 2, nu = 1
  )
})
