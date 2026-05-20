# Tests for the generalised Poisson distribution

test_that("genpois passes discrete distribution checks (lambda=3, phi=0.3)", {
  check_discrete_dist(
    dfun        = dgenpois,
    pfun        = pgenpois,
    xs_int      = c(0, 1, 3, 5, 8),
    sum_support = 0:100,
    lambda = 3, phi = 0.3
  )
})

test_that("genpois passes discrete distribution checks (lambda=5, phi=0.5)", {
  check_discrete_dist(
    dfun        = dgenpois,
    pfun        = pgenpois,
    xs_int      = c(0, 2, 4, 7, 11),
    sum_support = 0:200,
    lambda = 5, phi = 0.5
  )
})
