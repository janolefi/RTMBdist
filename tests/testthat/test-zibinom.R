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

test_that("pzibinom includes the binomial zeros at q = 0 and matches the mass function", {
  expect_equal(pzibinom(0, 10, 0.4, 0.2), 0.2 + 0.8 * stats::dbinom(0, 10, 0.4))
  expect_equal(pzibinom(0:10, 10, 0.4, 0.2), cumsum(dzibinom(0:10, 10, 0.4, 0.2)))
  expect_equal(pzibinom(c(-2, -1), 10, 0.4, 0.2), c(0, 0))
  check_osa_cdf(dzibinom, pzibinom, c(0, 0, 2, 5, 10), prob = 0.4, zeroprob = 0.2,
                .fixed = list(size = 10), .discrete = TRUE)
})
