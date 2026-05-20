# Tests for the zero-inflated Weibull distribution

test_that("ziweibull passes zero-inflated distribution checks (shape=2, scale=1, zeroprob=0.3)", {
  check_zeroinfl_dist(
    dfun = dziweibull,
    pfun = pziweibull,
    xs   = c(0.3, 0.6, 1, 1.5, 2.5),
    shape = 2, scale = 1, zeroprob = 0.3
  )
})

test_that("ziweibull passes zero-inflated distribution checks (shape=0.5, scale=2, zeroprob=0.15)", {
  check_zeroinfl_dist(
    dfun = dziweibull,
    pfun = pziweibull,
    xs   = c(0.2, 0.7, 1.5, 3, 6),
    shape = 0.5, scale = 2, zeroprob = 0.15
  )
})
