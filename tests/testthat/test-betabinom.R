# Tests for the beta-binomial distribution
# Support is 0:size; there is no quantile function.

test_that("betabinom passes discrete distribution checks (size=10, shape1=2, shape2=3)", {
  check_discrete_dist(
    dfun        = dbetabinom,
    pfun        = pbetabinom,
    xs_int      = c(0, 2, 5, 8, 10),
    sum_support = 0:10,
    size = 10, shape1 = 2, shape2 = 3
  )
})

test_that("betabinom passes discrete distribution checks (size=20, shape1=0.5, shape2=0.5)", {
  check_discrete_dist(
    dfun        = dbetabinom,
    pfun        = pbetabinom,
    xs_int      = c(0, 5, 10, 15, 20),
    sum_support = 0:20,
    size = 20, shape1 = 0.5, shape2 = 0.5
  )
})

test_that("betabinom AD gradient has no NaN", {
  check_ad_gradient(dbetabinom, rbetabinom, size = 20, shape1 = 2, shape2 = 3)
})

test_that("pbetabinom matches the binomial distribution function mixed over the beta prior", {
  ref <- function(q, lower.tail = TRUE) sapply(q, function(x) stats::integrate(function(p)
    stats::pbinom(x, 10, p, lower.tail = lower.tail) * stats::dbeta(p, 2, 30), 0, 1,
    rel.tol = 1e-12)$value)
  q <- c(0, 1, 3, 6, 9)
  expect_equal(pbetabinom(q, 10, 2, 30), ref(q), tolerance = 1e-10)
  # small upper tails are summed directly rather than taken as 1 - p
  expect_equal(pbetabinom(q, 10, 2, 30, lower.tail = FALSE), ref(q, FALSE), tolerance = 1e-10)
  expect_identical(pbetabinom(c(10, 12), 10, 2, 3), c(1, 1))
  expect_equal(pbetabinom(c(-1, NA), 10, 2, 3), c(0, NA))
  # parameters per element
  expect_equal(pbetabinom(c(2, 5), c(10, 12), c(1, 3), c(2, 4)),
               c(sum(dbetabinom(0:2, 10, 1, 2)), sum(dbetabinom(0:5, 12, 3, 4))), tolerance = 1e-12)
})

test_that("pbetabinom under AD matches pbetabinom outside AD", {
  q <- c(0, 2, 5, 8, 10)
  check_ad_cdf(pbetabinom, dbetabinom, q, shape1 = 2, shape2 = 3, .fixed = list(size = 10), .dq = FALSE)
  # below size, as P(X > size) = 0 has no finite log
  check_ad_cdf(pbetabinom, dbetabinom, q[q < 10], shape1 = 2, shape2 = 3, .fixed = list(size = 10),
               .args = list(lower.tail = FALSE, log.p = TRUE))
  expect_error(RTMB::MakeTape(function(x) pbetabinom(x, 10, 2, 3), 1), "numeric data")
})

test_that("betabinom supports OSA residuals", {
  check_osa_cdf(dbetabinom, pbetabinom, c(0, 2, 5, 8, 10), shape1 = 2, shape2 = 3,
                .fixed = list(size = 10), .discrete = TRUE)
})
