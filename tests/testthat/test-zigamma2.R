# Tests for the zero-inflated reparameterised gamma distribution (mean/sd)

test_that("zigamma2 passes zero-inflated distribution checks (mean=2, sd=1, zeroprob=0.2)", {
  check_zeroinfl_dist(
    dfun = dzigamma2,
    pfun = pzigamma2,
    xs   = c(0.5, 1, 2, 3, 5),
    mean = 2, sd = 1, zeroprob = 0.2
  )
})

test_that("zigamma2 passes zero-inflated distribution checks (mean=1, sd=2, zeroprob=0.4)", {
  check_zeroinfl_dist(
    dfun = dzigamma2,
    pfun = pzigamma2,
    xs   = c(0.1, 0.5, 1, 2, 5),
    mean = 1, sd = 2, zeroprob = 0.4
  )
})
