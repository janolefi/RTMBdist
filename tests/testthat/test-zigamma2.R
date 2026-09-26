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

test_that("zigamma2 AD gradient has no NaN", {
  check_ad_gradient(dzigamma2,  rzigamma2,  mean = 2, sd = 1, zeroprob = 0.2)
})

test_that("pzigamma2 matches the gamma distribution with the implied scale", {
  q <- c(0.5, 1, 2, 3, 5)
  # mean 2 and sd 1 give shape 4 and scale 0.5
  expect_equal(pzigamma2(q, 2, 1, 0.2), 0.2 + 0.8 * stats::pgamma(q, shape = 4, scale = 0.5))
  check_ad_cdf(pzigamma2, dzigamma2, q, mean = 2, sd = 1, zeroprob = 0.2)
})
