# Tests for the zero-truncated binomial distribution
# Support starts at 1.

test_that("ztbinom passes discrete distribution checks (size=10, prob=0.4)", {
  check_discrete_dist(
    dfun        = dztbinom,
    pfun        = pztbinom,
    xs_int      = c(1, 3, 5, 7, 10),
    sum_support = 1:10,
    size = 10, prob = 0.4
  )
})

test_that("ztbinom passes discrete distribution checks (size=20, prob=0.3)", {
  check_discrete_dist(
    dfun        = dztbinom,
    pfun        = pztbinom,
    xs_int      = c(1, 4, 8, 12, 18),
    sum_support = 1:20,
    size = 20, prob = 0.3
  )
})
