# Tests for the truncated normal distribution

test_that("truncnorm passes standard distribution checks (mean=0, sd=1, min=-2, max=2)", {
  check_continuous_dist(
    dfun  = dtruncnorm,
    pfun  = ptruncnorm,
    qfun  = qtruncnorm,
    xs    = c(-1.5, -0.5, 0, 0.5, 1.5),
    lower = -2, upper = 2,
    mean = 0, sd = 1, min = -2, max = 2
  )
})

test_that("truncnorm passes standard distribution checks (mean=2, sd=1, min=0, max=Inf)", {
  check_continuous_dist(
    dfun  = dtruncnorm,
    pfun  = ptruncnorm,
    qfun  = qtruncnorm,
    xs    = c(0.3, 0.8, 1.5, 2.5, 4.0),
    lower = 0, upper = Inf,
    mean = 2, sd = 1, min = 0
  )
})

test_that("truncnorm AD gradient has no NaN", {
  check_ad_gradient(dtruncnorm, rtruncnorm, mean = 0, sd = 1, min = -3, max = 3)
})

test_that("truncnorm AD gradient is finite when a bound is infinite", {
  set.seed(1)
  x <- rtruncnorm(50, mean = 2, sd = 2, min = -1, max = Inf)

  grad <- function(min, max) {
    F <- RTMB::MakeTape(function(par) {
      -sum(dtruncnorm(x, mean = par[1], sd = exp(par[2]),
                      min = min, max = max, log = TRUE))
    }, c(2, log(2)))
    as.vector(F$jacobian(c(2, log(2))))
  }

  expect_false(any(is.nan(grad(-1, Inf))))
  expect_false(any(is.nan(grad(-Inf, 1))))
  expect_false(any(is.nan(grad(-Inf, Inf))))

  # an infinite bound must agree with a far-away finite one
  expect_equal(grad(-1, Inf), grad(-1, 1e8))
  expect_equal(grad(-Inf, Inf), grad(-1e8, 1e8))
})
