# Tests for the skew t distribution

test_that("skewt passes standard distribution checks (mu=0, sigma=1, skew=0, df=5)", {
  check_continuous_dist(
    dfun  = dskewt,
    pfun  = pskewt,
    qfun  = qskewt,
    xs    = c(-2, -1, 0, 1, 2),
    mu = 0, sigma = 1, skew = 0, df = 5
  )
})

test_that("skewt passes standard distribution checks (mu=0, sigma=1, skew=2, df=5)", {
  check_continuous_dist(
    dfun  = dskewt,
    pfun  = pskewt,
    qfun  = qskewt,
    xs    = c(-0.4, 0.2, 0.7, 1.3, 2.6),
    mu = 0, sigma = 1, skew = 2, df = 5
  )
})

test_that("skewt AD gradient has no NaN", {
  check_ad_gradient(dskewt,     rskewt,     mu = 0, sigma = 1, skew = 2, df = 5)
})

test_that("pskewt under AD matches pskewt outside AD", {
  # df is not taped here, because sn::pst is only exact for integer df, so finite
  # differences in df are not accurate enough, see the next test
  q <- c(-3, -0.4, 0.2, 0.7, 1.3, 2.6, 15)
  F <- check_ad_cdf(pskewt, dskewt, q, mu = 0.2, sigma = 1.3, skew = 2, .fixed = list(df = 5))
  # tape stays correct at other parameters, including a different skew direction
  expect_equal(F(c(1, 0.5, -3)), pskewt(q, 1, 0.5, -3, 5), tolerance = 1e-8)
  check_ad_cdf(pskewt, dskewt, c(-40, 0, 60), mu = 0, sigma = 1, skew = -1, .fixed = list(df = 2),
               .args = list(lower.tail = FALSE, log.p = TRUE))
})

test_that("pskewt under AD has the right derivative in df", {
  q <- c(-3, 0, 0.2, 2.6)
  ref <- function(df) sapply(q, function(x) stats::integrate(function(y) dskewt(y, 0, 1.3, 2, df),
    -Inf, x, rel.tol = 1e-13, abs.tol = 0)$value)
  F <- RTMB::MakeTape(function(p) pskewt(q, 0, 1.3, 2, p), 3.5)
  h <- 1e-4
  expect_equal(as.vector(F$jacobian(3.5)), (ref(3.5 + h) - ref(3.5 - h)) / (2 * h), tolerance = 1e-6)
})

test_that("skewt supports OSA residuals", {
  set.seed(1)
  check_osa_cdf(dskewt, pskewt, rskewt(30, 1, 2, 3, 5), mu = 1, sigma = 2, skew = 3, df = 5)
})

test_that("the skew derivative of dskewt and pskewt is not zero at skew = 0", {
  # this used to be exactly zero, so that skew could not move away from a start at zero
  x <- c(-1, 0.5, 2)
  h <- 1e-5
  for (f in list(dskewt, pskewt)) {
    g <- function(a) f(x, 0, 1, a, 5)
    F <- RTMB::MakeTape(g, 0)
    expect_equal(as.vector(F$jacobian(0)), (g(h) - g(-h)) / (2 * h), tolerance = 1e-6)
  }
})

test_that("pskewt under AD has correct second derivatives", {
  check_ad_cdf_hessian(pskewt, c(-1, 0.2, 0.5, 2, 4), mu = 0.3, sigma = 1.5, skew = 2,
                       .fixed = list(df = 5))
})

test_that("rskewt recycles vector parameters and matches pskewt", {
  # sn::rst is not vectorised over its parameters, which rskewt used to call
  set.seed(1)
  s <- rep(1:2, 5000)
  x <- rskewt(length(s), c(-3, 3)[s], c(1, 2)[s], c(-2, 4)[s], c(4, Inf)[s])
  expect_length(x, 10000)
  expect_gt(ks.test(x[s == 1], function(q) sn::pst(q, -3, 1, -2, 4))$p.value, 1e-3)
  expect_gt(ks.test(x[s == 2], function(q) sn::psn(q, 3, 2, 4))$p.value, 1e-3)
  expect_length(rskewt(0), 0)
})

test_that("pskewt and qskewt recycle vector skew and df", {
  # sn::pst and sn::qst only take a single alpha and nu
  mu <- c(-1, 0, 2, 0); sigma <- c(1, 2, 0.5, 1); skew <- c(-2, 3, 0, 3); df <- c(4, 7.5, Inf, 7.5)
  q <- c(-1.5, 0.3, 2.2, 4)
  ref <- mapply(function(q, m, o, a, v) sn::pst(q, m, o, a, v), q, mu, sigma, skew, df)
  expect_equal(pskewt(q, mu, sigma, skew, df), ref, tolerance = 1e-12)
  p <- c(0.1, 0.5, 0.8, 0.99)
  ref <- mapply(function(p, m, o, a, v) sn::qst(p, m, o, a, v), p, mu, sigma, skew, df)
  expect_equal(qskewt(p, mu, sigma, skew, df), ref, tolerance = 1e-8)
  expect_equal(pskewt2(qskewt2(p, 1, 2, skew, c(4, 7.5, 30, 7.5)), 1, 2, skew, c(4, 7.5, 30, 7.5)),
               p, tolerance = 1e-7)
})

test_that("pskewt is accurate in both tails", {
  # 1 - p cancels in the upper tail, and sn's far tail approximation (|z| > 30) for
  # non-integer df was off by 0.8% here
  tail_int <- function(q, a, nu) {
    integrate(function(u) sn::dst(q / u, 0, 1, a, nu) * abs(q) / u^2, 0, 1, rel.tol = 1e-13)$value
  }
  expect_equal(pskewt(40, 0, 1, 2, 5.5, lower.tail = FALSE), tail_int(40, 2, 5.5), tolerance = 1e-8)
  expect_equal(pskewt(12, 0, 1, 0, 5, lower.tail = FALSE), tail_int(12, 0, 5), tolerance = 1e-8)
  expect_equal(pskewt(-20, 0, 1, 2, 4, log.p = TRUE), log(tail_int(-20, 2, 4)), tolerance = 1e-8)
  expect_equal(pskewt(c(-20, 1, 40), 1, 2, 2, 5.5, lower.tail = FALSE),
               1 - pskewt(c(-20, 1, 40), 1, 2, 2, 5.5), tolerance = 1e-8)
})
