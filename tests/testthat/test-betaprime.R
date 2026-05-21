# Tests for the beta prime distribution
# If X ~ Beta(a, b), then X/(1-X) ~ BetaPrime(a, b)
# CDF: pbetaprime(q, a, b) = pbeta(q/(1+q), a, b)

test_that("betaprime passes standard distribution checks (shape1=2, shape2=3)", {
  check_continuous_dist(
    dfun  = dbetaprime,
    pfun  = pbetaprime,
    qfun  = qbetaprime,
    xs    = c(0.1, 0.5, 1, 2, 4),
    lower = 0, upper = Inf,
    shape1 = 2, shape2 = 3
  )
})

test_that("betaprime passes standard distribution checks (shape1=3, shape2=5)", {
  check_continuous_dist(
    dfun  = dbetaprime,
    pfun  = pbetaprime,
    qfun  = qbetaprime,
    xs    = c(0.1, 0.3, 0.6, 1, 2),
    lower = 0, upper = Inf,
    shape1 = 3, shape2 = 5
  )
})

test_that("pbetaprime matches pbeta via the q/(1+q) transform", {
  q <- c(0.1, 0.5, 1, 2, 4)
  a <- 2; b <- 3
  expect_equal(
    pbetaprime(q, shape1 = a, shape2 = b),
    pbeta(q / (1 + q), shape1 = a, shape2 = b),
    tolerance = 1e-10
  )
})

test_that("betaprime AD gradient has no NaN", {
  check_ad_gradient(dbetaprime, rbetaprime, shape1 = 2, shape2 = 3)
})
