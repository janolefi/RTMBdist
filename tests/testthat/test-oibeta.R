# Tests for the one-inflated beta distribution
# Point mass at 1; continuous part on (0, 1)

test_that("oibeta passes inflated distribution checks (shape1=2, shape2=3, oneprob=0.2)", {
  check_inflated_dist(
    dfun         = doibeta,
    pfun         = poibeta,
    xs           = c(0.1, 0.3, 0.5, 0.7, 0.9),
    lower        = 0, upper = 1,
    point_masses = 1,
    shape1 = 2, shape2 = 3, oneprob = 0.2
  )
})

test_that("oibeta passes inflated distribution checks (shape1=0.5, shape2=2, oneprob=0.3)", {
  check_inflated_dist(
    dfun         = doibeta,
    pfun         = poibeta,
    xs           = c(0.05, 0.2, 0.5, 0.8, 0.95),
    lower        = 0, upper = 1,
    point_masses = 1,
    shape1 = 0.5, shape2 = 2, oneprob = 0.3
  )
})
