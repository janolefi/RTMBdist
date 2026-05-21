# Tests for the zero-inflated inverse Gaussian distribution

test_that("ziinvgauss passes zero-inflated distribution checks (mean=1, shape=1, zeroprob=0.25)", {
  check_zeroinfl_dist(
    dfun = dziinvgauss,
    pfun = pziinvgauss,
    xs   = c(0.3, 0.6, 1, 2, 4),
    mean = 1, shape = 1, zeroprob = 0.25
  )
})

test_that("ziinvgauss passes zero-inflated distribution checks (mean=2, shape=3, zeroprob=0.1)", {
  check_zeroinfl_dist(
    dfun = dziinvgauss,
    pfun = pziinvgauss,
    xs   = c(0.5, 1, 2, 3, 5),
    mean = 2, shape = 3, zeroprob = 0.1
  )
})

test_that("ziinvgauss AD gradient has no NaN", {
  check_ad_gradient(dziinvgauss,rziinvgauss,mean = 2, shape = 3, zeroprob = 0.2)
})
