# Tests for the beta-binomial distribution
# Support is 0:size; no p or q functions exist so only tests 1 and 4 apply.

test_that("betabinom passes discrete distribution checks (size=10, shape1=2, shape2=3)", {
  check_discrete_dist(
    dfun        = dbetabinom,
    pfun        = NULL,
    xs_int      = c(0, 2, 5, 8, 10),
    sum_support = 0:10,
    size = 10, shape1 = 2, shape2 = 3
  )
})

test_that("betabinom passes discrete distribution checks (size=20, shape1=0.5, shape2=0.5)", {
  check_discrete_dist(
    dfun        = dbetabinom,
    pfun        = NULL,
    xs_int      = c(0, 5, 10, 15, 20),
    sum_support = 0:20,
    size = 20, shape1 = 0.5, shape2 = 0.5
  )
})

test_that("betabinom AD gradient has no NaN", {
  check_ad_gradient(dbetabinom, rbetabinom, size = 20, shape1 = 2, shape2 = 3)
})
