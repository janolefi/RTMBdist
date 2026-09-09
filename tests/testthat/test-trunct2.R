# Tests for the truncated location-scale t distribution

test_that("trunct2 passes standard distribution checks (df=5, mu=0, sigma=1, min=-2, max=2)", {
  check_continuous_dist(
    dfun  = dtrunct2,
    pfun  = ptrunct2,
    qfun  = qtrunct2,
    xs    = c(-1.5, -0.5, 0, 0.5, 1.5),
    lower = -2, upper = 2,
    df = 5, mu = 0, sigma = 1, min = -2, max = 2
  )
})

test_that("trunct2 passes standard distribution checks (df=10, mu=2, sigma=1, min=0, max=Inf)", {
  check_continuous_dist(
    dfun  = dtrunct2,
    pfun  = ptrunct2,
    qfun  = qtrunct2,
    xs    = c(0.5, 1.0, 2.0, 3.0, 4.0),
    lower = 0, upper = Inf,
    df = 10, mu = 2, sigma = 1, min = 0
  )
})

test_that("trunct2 AD gradient has no NaN", {
  check_ad_gradient(dtrunct2,   rtrunct2,   df = 5, mu = 0, sigma = 1, min = -3, max = 3)
})

test_that("trunct2 AD gradient is finite when a bound is infinite", {
  set.seed(1)
  x <- rtrunct2(50, df = 5, mu = 2, sigma = 2, min = -1, max = Inf)

  grad <- function(min, max) {
    F <- RTMB::MakeTape(function(par) {
      -sum(dtrunct2(x, df = 5, mu = par[1], sigma = exp(par[2]),
                    min = min, max = max, log = TRUE))
    }, c(2, log(2)))
    as.vector(F$jacobian(c(2, log(2))))
  }

  expect_false(any(is.nan(grad(-1, Inf))))
  expect_false(any(is.nan(grad(-Inf, 1))))
  expect_false(any(is.nan(grad(-Inf, Inf))))

  expect_equal(grad(-1, Inf), grad(-1, 1e8))
  expect_equal(grad(-Inf, Inf), grad(-1e8, 1e8))
})
