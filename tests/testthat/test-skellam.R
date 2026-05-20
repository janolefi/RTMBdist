# Tests for the Skellam distribution
# Support is all integers; no p or q functions exist so only tests 1 and 4 apply.

test_that("skellam passes discrete distribution checks (mu1=2, mu2=3)", {
  check_discrete_dist(
    dfun        = dskellam,
    pfun        = NULL,
    xs_int      = c(-5, -2, 0, 2, 4),
    sum_support = -100:100,
    mu1 = 2, mu2 = 3
  )
})

test_that("skellam passes discrete distribution checks (mu1=5, mu2=1)", {
  check_discrete_dist(
    dfun        = dskellam,
    pfun        = NULL,
    xs_int      = c(-2, 0, 3, 6, 9),
    sum_support = -100:100,
    mu1 = 5, mu2 = 1
  )
})
