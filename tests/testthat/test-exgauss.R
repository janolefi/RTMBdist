# Tests for the ex-Gaussian distribution
# X = normal(mu, sigma) + exponential(lambda)

test_that("exgauss passes standard distribution checks (mu=0, sigma=1, lambda=1)", {
  check_continuous_dist(
    dfun  = dexgauss,
    pfun  = pexgauss,
    qfun  = qexgauss,
    xs    = c(-1, 0, 1, 2, 4),
    lower = -5, upper = 20,
    mu = 0, sigma = 1, lambda = 1
  )
})

test_that("exgauss passes standard distribution checks (mu=2, sigma=0.5, lambda=2)", {
  check_continuous_dist(
    dfun  = dexgauss,
    pfun  = pexgauss,
    qfun  = qexgauss,
    xs    = c(1, 2, 2.5, 3, 4),
    lower = -1, upper = 15,
    mu = 2, sigma = 0.5, lambda = 2
  )
})

test_that("exgauss AD gradient has no NaN", {
  check_ad_gradient(dexgauss,   rexgauss,   mu = 2, sigma = 0.5, lambda = 2)
})
