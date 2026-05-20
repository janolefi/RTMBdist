# Tests for the zero-one-inflated beta distribution (mu/phi parameterisation)
# Point masses at 0 and 1; continuous part on (0, 1)

test_that("zoibeta2 passes inflated distribution checks (mu=0.4, phi=5, zeroprob=0.15, oneprob=0.1)", {
  check_inflated_dist(
    dfun         = dzoibeta2,
    pfun         = pzoibeta2,
    xs           = c(0.1, 0.3, 0.5, 0.7, 0.9),
    lower        = 0, upper = 1,
    point_masses = c(0, 1),
    mu = 0.4, phi = 5, zeroprob = 0.15, oneprob = 0.1
  )
})

test_that("zoibeta2 passes inflated distribution checks (mu=0.6, phi=8, zeroprob=0.2, oneprob=0.2)", {
  check_inflated_dist(
    dfun         = dzoibeta2,
    pfun         = pzoibeta2,
    xs           = c(0.05, 0.25, 0.5, 0.75, 0.95),
    lower        = 0, upper = 1,
    point_masses = c(0, 1),
    mu = 0.6, phi = 8, zeroprob = 0.2, oneprob = 0.2
  )
})
