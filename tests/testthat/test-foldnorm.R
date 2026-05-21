# Tests for the folded normal distribution
# No quantile function exists, so the round-trip test is skipped (qfun = NULL)

test_that("foldnorm passes standard distribution checks (mu=0, sigma=1)", {
  check_continuous_dist(
    dfun  = dfoldnorm,
    pfun  = pfoldnorm,
    qfun  = NULL,
    xs    = c(0.2, 0.5, 1, 2, 3),
    lower = 0, upper = Inf,
    mu = 0, sigma = 1
  )
})

test_that("foldnorm passes standard distribution checks (mu=1, sigma=2)", {
  check_continuous_dist(
    dfun  = dfoldnorm,
    pfun  = pfoldnorm,
    qfun  = NULL,
    xs    = c(0.1, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    mu = 1, sigma = 2
  )
})

test_that("foldnorm AD gradient has no NaN", {
  check_ad_gradient(dfoldnorm,  rfoldnorm,  mu = 1, sigma = 2)
})
