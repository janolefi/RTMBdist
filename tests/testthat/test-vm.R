# Tests for the von Mises distribution (circular, support [-pi, pi])
# No quantile function exists, so the round-trip test is skipped (qfun = NULL)

test_that("vm passes standard distribution checks (mu=0, kappa=1)", {
  check_continuous_dist(
    dfun  = dvm,
    pfun  = pvm,
    qfun  = NULL,
    xs    = c(-2, -1, 0, 1, 2),
    lower = -pi, upper = pi,
    mu = 0, kappa = 1
  )
})

test_that("vm passes standard distribution checks (mu=1, kappa=3)", {
  check_continuous_dist(
    dfun  = dvm,
    pfun  = pvm,
    qfun  = NULL,
    xs    = c(-1, 0, 1, 2, 2.5),
    lower = -pi, upper = pi,
    mu = 1, kappa = 3
  )
})
