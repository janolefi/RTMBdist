# Tests for the zero-inflated gamma distribution

test_that("zigamma passes zero-inflated distribution checks (shape=2, scale=1, zeroprob=0.3)", {
  check_zeroinfl_dist(
    dfun = dzigamma,
    pfun = pzigamma,
    xs   = c(0.5, 1, 2, 3, 5),
    shape = 2, scale = 1, zeroprob = 0.3
  )
})

test_that("zigamma passes zero-inflated distribution checks (shape=0.5, scale=2, zeroprob=0.1)", {
  check_zeroinfl_dist(
    dfun = dzigamma,
    pfun = pzigamma,
    xs   = c(0.1, 0.5, 1, 3, 6),
    shape = 0.5, scale = 2, zeroprob = 0.1
  )
})

test_that("zigamma AD gradient has no NaN", {
  check_ad_gradient(dzigamma,   rzigamma,   shape = 2, scale = 1, zeroprob = 0.2)
})
