# Tests for the mean-parameterised beta-negative binomial distribution

test_that("bnbinom2 passes standard discrete checks", {
  check_discrete_dist(
    dfun        = dbnbinom2,
    pfun        = NULL, # no closed-form distribution function
    xs_int      = c(0, 1, 4, 10, 30),
    sum_support = 0:20000,
    mu = 4, sigma = 0.4, nu = 0.5
  )
  check_discrete_dist(
    dfun        = dbnbinom2,
    pfun        = NULL,
    xs_int      = c(0, 3, 9, 20, 60),
    sum_support = 0:2000,
    mu = 10, sigma = 0.2, nu = 0.4
  )
})

test_that("bnbinom2 AD gradient has no NaN", {
  check_ad_gradient(dbnbinom2, rbnbinom2, mu = 4, sigma = 0.4, nu = 0.5)
  check_ad_gradient(dbnbinom2, rbnbinom2, mu = 12, sigma = 0.8, nu = 1.5)
})

test_that("bnbinom2 is the bnbinom reparameterisation it claims to be", {
  k <- 0:60
  mu <- 4; sigma <- 0.4; nu <- 0.5
  expect_equal(
    dbnbinom2(k, mu, sigma, nu),
    dbnbinom(k, size = 1 / nu, shape1 = 1 / sigma + 1, shape2 = mu * nu / sigma)
  )
  expect_equal(dbnbinom2(c(-2, -0.5), mu, sigma, nu), c(0, 0))
})

test_that("mu is exactly the mean and the variance matches the closed form", {
  k <- 0:2000000
  for (p in list(c(4, 0.4, 0.5), c(10, 0.3, 1.5))) {
    mu <- p[1]; sigma <- p[2]; nu <- p[3]
    d <- dbnbinom2(k, mu, sigma, nu)
    m <- sum(k * d)
    expect_equal(m, mu, tolerance = 1e-6)
    expect_equal(sum((k - m)^2 * d),
                 mu * (sigma + nu) * (mu * nu + 1) / (nu * (1 - sigma)),
                 tolerance = 1e-5)
  }
})

test_that("bnbinom2 collapses to the negative binomial as sigma goes to zero", {
  k <- 0:60
  expect_equal(dbnbinom2(k, 5, 1e-9, 0.5), stats::dnbinom(k, size = 1 / 0.5, mu = 5),
               tolerance = 1e-6)
})

test_that("rbnbinom2 draws follow the mass function", {
  set.seed(2)
  x <- rbnbinom2(2e4, 4, 0.4, 0.5)
  expect_true(all(x == floor(x)) && min(x) >= 0)
  e <- c(dbnbinom2(0:14, 4, 0.4, 0.5), 1 - sum(dbnbinom2(0:14, 4, 0.4, 0.5)))
  o <- as.vector(table(factor(pmin(x, 15), levels = 0:15)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
})

test_that("bnbinom2 recovers its parameters where bnbinom stalls on the ridge", {
  # the point of the reparameterisation: mu is pinned to the mean, so the
  # likelihood no longer has the size / shape2 ridge the original one has
  set.seed(7)
  x <- rbnbinom2(50000, 4, 0.4, 0.5)
  f <- function(p) -sum(dbnbinom2(x, exp(p[1]), exp(p[2]), exp(p[3]), log = TRUE))
  F <- RTMB::MakeTape(f, log(c(4, 0.4, 0.5)))
  o <- nlminb(log(c(1, 1, 1)), F, function(p) as.vector(F$jacobian(p)))
  expect_equal(o$convergence, 0)
  expect_equal(exp(o$par), c(4, 0.4, 0.5), tolerance = 0.15)
})

test_that("bnbinom2 refuses OSA residuals and recycles its arguments", {
  expect_length(dbnbinom2(0:5, c(1, 2), 0.5, 0.5), 6)
  expect_length(rbnbinom2(5, c(2, 4), 0.5, 0.5), 5)
  expect_error(
    dbnbinom2(structure(1, class = "osa"), 4, 0.4, 0.5),
    "does not support OSA"
  )
})

test_that("dbnbinom2 rejects non-positive parameters", {
  expect_error(dbnbinom2(1, 0, 0.4, 0.5), "mu")
  expect_error(dbnbinom2(1, 4, -1, 0.5), "sigma")
  expect_error(dbnbinom2(1, 4, 0.4, 0), "nu")
})
