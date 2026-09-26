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

test_that("powerexp AD gradient has no NaN", {
  check_ad_gradient(dpowerexp,  rpowerexp,  mu = 0, sigma = 1, nu = 2)
})

test_that("ppowerexp and ppowerexp2 have the right gradients at and around q = mu", {
  # the derivative used to be NaN at q = mu, for nu below as well as above 1
  for (nu in c(0.5, 1, 1.5, 2, 5)) {
    for (pd in list(c(ppowerexp, dpowerexp), c(ppowerexp2, dpowerexp2))) {
      pfun <- pd[[1]]; dfun <- pd[[2]]
      check_ad_cdf(pfun, dfun, c(-2, 1 - 1e-2, 1 + 1e-2, 3), mu = 1, sigma = 1.5, nu = nu)
      # at q = mu the density has a cusp for nu <= 1, where finite differences are
      # inaccurate; there F = 1/2 for all sigma and nu, so the gradient is (-f(mu), 0, 0)
      F <- RTMB::MakeTape(function(p) pfun(1, p[1], p[2], p[3]), c(1, 1.5, nu))
      expect_equal(as.vector(F$jacobian(c(1, 1.5, nu))), c(-dfun(1, 1, 1.5, nu), 0, 0),
                   tolerance = 1e-10)
      G <- RTMB::MakeTape(function(x) pfun(x, 1, 1.5, nu), 1)
      expect_equal(as.vector(G$jacobian(1)), dfun(1, 1, 1.5, nu), tolerance = 1e-10)
    }
  }
})
