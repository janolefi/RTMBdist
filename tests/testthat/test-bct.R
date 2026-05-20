# Tests for the Box-Cox t distribution

test_that("bct passes standard distribution checks (mu=5, sigma=0.1, nu=1, tau=5)", {
  check_continuous_dist(
    dfun  = dbct,
    pfun  = pbct,
    qfun  = qbct,
    xs    = c(4.0, 4.6, 5.0, 5.4, 6.0),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.1, nu = 1, tau = 5
  )
})

test_that("bct passes standard distribution checks (mu=5, sigma=0.3, nu=2, tau=10)", {
  check_continuous_dist(
    dfun  = dbct,
    pfun  = pbct,
    qfun  = qbct,
    xs    = c(1.5, 3.0, 5.0, 7.5, 11.0),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.3, nu = 2, tau = 10
  )
})
