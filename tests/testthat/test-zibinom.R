# Tests for the zero-inflated binomial distribution

test_that("zibinom passes discrete distribution checks (size=10, prob=0.4, zeroprob=0.2)", {
  check_discrete_dist(
    dfun        = dzibinom,
    pfun        = pzibinom,
    xs_int      = c(0, 2, 4, 7, 10),
    sum_support = 0:10,
    size = 10, prob = 0.4, zeroprob = 0.2
  )
})

test_that("zibinom passes discrete distribution checks (size=20, prob=0.6, zeroprob=0.3)", {
  check_discrete_dist(
    dfun        = dzibinom,
    pfun        = pzibinom,
    xs_int      = c(0, 5, 10, 15, 20),
    sum_support = 0:20,
    size = 20, prob = 0.6, zeroprob = 0.3
  )
})

test_that("zibinom AD gradient has no NaN", {
  check_ad_gradient(dzibinom,   rzibinom,   size = 10, prob = 0.4, zeroprob = 0.2)
})
