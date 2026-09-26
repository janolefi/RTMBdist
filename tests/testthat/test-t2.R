# Tests for the location-scale t distribution

test_that("t2 passes standard distribution checks (mu=0, sigma=1, df=5)", {
  check_continuous_dist(
    dfun  = dt2,
    pfun  = pt2,
    qfun  = qt2,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, df = 5
  )
})

test_that("t2 passes standard distribution checks (mu=2, sigma=3, df=10)", {
  check_continuous_dist(
    dfun  = dt2,
    pfun  = pt2,
    qfun  = qt2,
    xs    = c(-4, 0, 2, 5, 8),
    mu = 2, sigma = 3, df = 10
  )
})

test_that("t2 AD gradient has no NaN", {
  check_ad_gradient(dt2,        rt2,        mu = 0, sigma = 1, df = 5)
})

test_that("pt.ad matches stats::pt, including derivatives at and around q = 0", {
  q <- c(-3, -1e-3, -2e-4, 0, 2e-4, 1e-3, 0.4, 5)
  h <- 1e-5
  for (df in c(0.5, 1, 5, 50)) {
    F <- RTMB::MakeTape(function(x) pt.ad(x, df), q)
    expect_equal(F(q), stats::pt(q, df), tolerance = 1e-12)
    # dF/dq is the density, also at q = 0, where it used to be 0
    expect_equal(diag(F$jacobian(q)), stats::dt(q, df), tolerance = 1e-9)
    G <- RTMB::MakeTape(function(d) pt.ad(q, d), df)
    expect_equal(as.vector(G$jacobian(df)), (stats::pt(q, df + h) - stats::pt(q, df - h)) / (2 * h),
                 tolerance = 1e-6)
  }
  # higher derivatives at 0, which the Laplace approximation needs: F'' = 0, F''' = -(df + 1) / df * f(0)
  H1 <- RTMB::MakeTape(function(x) pt.ad(x, 5), 0)$jacfun()
  H2 <- H1$jacfun()
  H3 <- H2$jacfun()
  expect_equal(c(H1(0), H2(0), H3(0)), c(stats::dt(0, 5), 0, -6 / 5 * stats::dt(0, 5)), tolerance = 1e-10)
})

test_that("pt2 has the right gradient at q = mu", {
  check_ad_cdf(pt2, dt2, c(-1, 2, 2 + 1e-4, 5), mu = 2, sigma = 3, df = 10)
})
