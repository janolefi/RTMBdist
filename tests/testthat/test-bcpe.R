# Tests for the Box-Cox Power Exponential distribution

test_that("bcpe passes standard distribution checks (mu=5, sigma=0.1, nu=1, tau=2)", {
  check_continuous_dist(
    dfun  = dbcpe,
    pfun  = pbcpe,
    qfun  = qbcpe,
    xs    = c(4.0, 4.5, 5.0, 5.5, 6.2),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.1, nu = 1, tau = 2
  )
})

test_that("bcpe passes standard distribution checks (mu=5, sigma=0.3, nu=2, tau=1.5)", {
  check_continuous_dist(
    dfun  = dbcpe,
    pfun  = pbcpe,
    qfun  = qbcpe,
    xs    = c(1.5, 3.0, 5.0, 7.5, 11.0),
    lower = 0, upper = Inf,
    mu = 5, sigma = 0.3, nu = 2, tau = 1.5
  )
})

test_that("bcpe AD gradient has no NaN", {
  check_ad_gradient(dbcpe,      rbcpe,      mu = 5, sigma = 0.3, nu = 2, tau = 1.5)
})

test_that("pbcpe has the right gradient at and around q = mu", {
  # q = mu gives z = 0, where the derivative used to be NaN
  h <- 1e-5
  for (tau in c(0.8, 2, 4)) {
    check_ad_cdf(pbcpe, dbcpe, c(3, 5 - 0.05, 5 + 0.05, 8), mu = 5, sigma = 0.2, nu = 0.5, tau = tau)
    # at q = mu the density has a cusp for tau <= 1, where finite differences in mu are
    # inaccurate; as mu is a scale parameter, dF/dmu = -f(mu) there
    f <- function(p) pbcpe(5, p[1], p[2], p[3], p[4])
    par <- c(5, 0.2, 0.5, tau)
    J <- as.vector(RTMB::MakeTape(f, par)$jacobian(par))
    expect_equal(J[1], -dbcpe(5, 5, 0.2, 0.5, tau), tolerance = 1e-10)
    J_fd <- sapply(2:4, function(j) { e <- replace(numeric(4), j, h); (f(par + e) - f(par - e)) / (2 * h) })
    expect_equal(J[2:4], J_fd, tolerance = 1e-6)
  }
})
