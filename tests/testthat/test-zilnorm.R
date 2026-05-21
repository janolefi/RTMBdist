# Tests for the zero-inflated log-normal distribution

test_that("zilnorm passes zero-inflated distribution checks (meanlog=0, sdlog=1, zeroprob=0.2)", {
  check_zeroinfl_dist(
    dfun = dzilnorm,
    pfun = pzilnorm,
    xs   = c(0.3, 0.7, 1, 2, 4),
    meanlog = 0, sdlog = 1, zeroprob = 0.2
  )
})

test_that("zilnorm passes zero-inflated distribution checks (meanlog=1, sdlog=0.5, zeroprob=0.35)", {
  check_zeroinfl_dist(
    dfun = dzilnorm,
    pfun = pzilnorm,
    xs   = c(1, 2, 3, 5, 8),
    meanlog = 1, sdlog = 0.5, zeroprob = 0.35
  )
})

test_that("zilnorm AD gradient has no NaN", {
  check_ad_gradient(dzilnorm,   rzilnorm,   meanlog = 0, sdlog = 1, zeroprob = 0.2)
})
