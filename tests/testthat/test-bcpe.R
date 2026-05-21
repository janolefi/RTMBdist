# Tests for the Box-Cox Power Exponential distribution

test_that("bcpe passes standard distribution checks (mu=5, sigma=0.1, nu=1, tau=2)", {
  check_continuous_dist(
    dfun  = dbcpe,
    pfun  = pbcpe,
    qfun  = qbcpe,
    xs    = c(4.0, 4.5, 5.0, 5.5, 6.2),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.1, nu = 1, tau = 2
  )
})

test_that("bcpe passes standard distribution checks (mu=5, sigma=0.3, nu=2, tau=1.5)", {
  check_continuous_dist(
    dfun  = dbcpe,
    pfun  = pbcpe,
    qfun  = qbcpe,
    xs    = c(1.5, 3.0, 5.0, 7.5, 11.0),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.3, nu = 2, tau = 1.5
  )
})

test_that("bcpe AD gradient has no NaN", {
  check_ad_gradient(dbcpe,      rbcpe,      mu = 5, sigma = 0.3, nu = 2, tau = 1.5)
})
