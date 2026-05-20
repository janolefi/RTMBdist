# Tests for the one-inflated beta distribution (mu/phi parameterisation)
# Point mass at 1; continuous part on (0, 1)

test_that("oibeta2 passes inflated distribution checks (mu=0.3, phi=5, oneprob=0.2)", {
  check_inflated_dist(
    dfun         = doibeta2,
    pfun         = poibeta2,
    xs           = c(0.05, 0.2, 0.4, 0.6, 0.85),
    lower        = 0, upper = 1,
    point_masses = 1,
    mu = 0.3, phi = 5, oneprob = 0.2
  )
})

test_that("oibeta2 passes inflated distribution checks (mu=0.7, phi=10, oneprob=0.15)", {
  check_inflated_dist(
    dfun         = doibeta2,
    pfun         = poibeta2,
    xs           = c(0.1, 0.4, 0.6, 0.8, 0.95),
    lower        = 0, upper = 1,
    point_masses = 1,
    mu = 0.7, phi = 10, oneprob = 0.15
  )
})
