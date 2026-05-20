# Tests for the Box-Cox Cole and Green distribution

test_that("bccg passes standard distribution checks (mu=1, sigma=0.1, nu=1)", {
  check_continuous_dist(
    dfun  = dbccg,
    pfun  = pbccg,
    qfun  = qbccg,
    xs    = c(0.85, 0.93, 1.0, 1.07, 1.16),
    lower = 0, upper = Inf,
    mu = 1, sigma = 0.1, nu = 1
  )
})

test_that("bccg passes standard distribution checks (mu=3, sigma=0.3, nu=2)", {
  check_continuous_dist(
    dfun  = dbccg,
    pfun  = pbccg,
    qfun  = qbccg,
    xs    = c(1.5, 2.2, 3.0, 4.0, 5.5),
    lower = 0, upper = Inf,
    mu = 3, sigma = 0.3, nu = 2
  )
})
