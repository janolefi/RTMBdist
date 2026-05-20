# Tests for the zero-inflated beta distribution

test_that("zibeta passes zero-inflated distribution checks (shape1=2, shape2=3, zeroprob=0.2)", {
  check_zeroinfl_dist(
    dfun = dzibeta,
    pfun = pzibeta,
    xs   = c(0.1, 0.3, 0.5, 0.7, 0.9),
    upper = 1,
    shape1 = 2, shape2 = 3, zeroprob = 0.2
  )
})

test_that("zibeta passes zero-inflated distribution checks (shape1=0.5, shape2=0.5, zeroprob=0.3)", {
  check_zeroinfl_dist(
    dfun = dzibeta,
    pfun = pzibeta,
    xs   = c(0.05, 0.2, 0.5, 0.8, 0.95),
    upper = 1,
    shape1 = 0.5, shape2 = 0.5, zeroprob = 0.3
  )
})
