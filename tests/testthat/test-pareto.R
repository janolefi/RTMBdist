# Tests for the Pareto distribution

test_that("pareto passes standard distribution checks (mu=1)", {
  check_continuous_dist(
    dfun  = dpareto,
    pfun  = ppareto,
    qfun  = qpareto,
    xs    = c(1.2, 1.6, 2.5, 4.0, 9.0),
    lower = 1, upper = Inf,
    mu = 1
  )
})

test_that("pareto passes standard distribution checks (mu=5)", {
  check_continuous_dist(
    dfun  = dpareto,
    pfun  = ppareto,
    qfun  = qpareto,
    xs    = c(1.05, 1.15, 1.3, 1.7, 2.4),
    lower = 1, upper = Inf,
    mu = 5
  )
})

test_that("pareto AD gradient has no NaN", {
  check_ad_gradient(dpareto,    rpareto,    mu = 5)
})

test_that("qpareto honours lower.tail and log.p", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qpareto(p, mu = 2), qpareto(1 - p, mu = 2, lower.tail = FALSE))
  expect_equal(qpareto(p, mu = 2), qpareto(log(p), mu = 2, log.p = TRUE))
  # the support starts at 1 and p = 1 is the upper end point
  expect_equal(qpareto(0, mu = 2), 1)
  expect_equal(qpareto(1, mu = 2), Inf)
})
