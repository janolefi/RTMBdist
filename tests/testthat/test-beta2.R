# Tests for the beta2 distribution
# (beta reparameterised by mean mu and concentration phi)
#
# Parameters chosen so shape1 = mu*phi and shape2 = (1-mu)*phi are both > 1,
# which avoids boundary singularities that would complicate numerical integration.

test_that("beta2 passes standard distribution checks (mu=0.4, phi=5)", {
  # shape1 = 2, shape2 = 3
  check_continuous_dist(
    dfun  = dbeta2,
    pfun  = pbeta2,
    qfun  = qbeta2,
    xs    = c(0.1, 0.3, 0.5, 0.7, 0.9),
    lower = 0, upper = 1,
    mu = 0.4, phi = 5
  )
})

test_that("beta2 passes standard distribution checks (mu=0.7, phi=20)", {
  # shape1 = 14, shape2 = 6 — different region of parameter space
  check_continuous_dist(
    dfun  = dbeta2,
    pfun  = pbeta2,
    qfun  = qbeta2,
    xs    = c(0.05, 0.25, 0.5, 0.75, 0.95),
    lower = 0, upper = 1,
    mu = 0.7, phi = 20
  )
})
