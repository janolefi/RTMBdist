# Tests for the Waring distribution

test_that("waring passes standard discrete checks", {
  check_discrete_dist(
    dfun        = dwaring,
    pfun        = pwaring,
    xs_int      = c(0, 1, 4, 15, 60),
    sum_support = 0:5000,
    mu = 2, sigma = 0.5
  )
  check_discrete_dist(
    dfun        = dwaring,
    pfun        = pwaring,
    xs_int      = c(0, 2, 8, 30),
    sum_support = 0:100000, # heavier tail
    mu = 1, sigma = 0.9
  )
})

test_that("waring AD gradient has no NaN", {
  check_ad_gradient(dwaring, rwaring, mu = 2, sigma = 0.5)
  check_ad_gradient(dwaring, rwaring, mu = 8, sigma = 1.4)
})

test_that("dwaring and pwaring match gamlss.dist", {
  # reference values from gamlss.dist::dWARING / pWARING, baked in: gamlss.dist
  # is scheduled for archival, so the test must not depend on it
  h <- expand.grid(x = c(0, 1, 3, 10, 40), mu = c(0.4, 2, 9), sigma = c(0.2, 1, 2.2))
  logd <- c(
    -0.287682072451781, -1.79175946922805, -4.00733318523247,
    -8.5762165318657, -16.3227910364791, -0.980829253011727,
    -1.5114575040739, -2.46346631855013, -5.03261420081503,
    -10.8978498198304, -2.14006616349627, -2.28464739230738,
    -2.56513435422462, -3.46667975300553, -6.38729376081168,
    -0.182321556793955, -2.32238772029023, -4.27845024080956,
    -7.17978694885177, -11.0559303427293, -0.693147180559945,
    -1.6094379124341, -2.86220088092947, -5.2040066870768, -8.79815271810719,
    -1.70474809223843, -1.99243016469021, -2.49595648597458,
    -3.79173683955364, -6.54271208537289, -0.117783035656383,
    -2.79193168508291, -4.66963358411169, -7.21863792295076,
    -10.4669688467109, -0.485507815781701, -1.79384063543188,
    -3.23232074972234, -5.51142153116156, -8.64407818016039,
    -1.33828514193353, -1.80828877117927, -2.54037985111975,
    -4.11592578542821, -6.8346034618606)
  cdf <- c(
    0.75, 0.916666666666667, 0.984848484848485, 0.999622926093514,
    0.999999429573144, 0.375, 0.595588235294118, 0.815531475748194,
    0.978260869565218, 0.999845850289802, 0.117647058823529,
    0.219457013574661, 0.384729218247531, 0.713811912867509,
    0.976160275034878, 0.833333333333333, 0.931372549019608,
    0.976430976430976, 0.996038483305037, 0.999680977121502, 0.5, 0.7,
    0.857142857142857, 0.967032967032967, 0.996828752642706,
    0.181818181818182, 0.318181818181818, 0.505494505494506,
    0.785714285714286, 0.964705882352942, 0.888888888888889,
    0.950191570881226, 0.979490646833446, 0.994870400914846,
    0.99921375938535, 0.615384615384615, 0.781704781704782,
    0.893935656647521, 0.969697301487663, 0.995045303278517,
    0.262295081967213, 0.426229508196721, 0.61567231605179, 0.84198880448076,
    0.967386967258838)
  expect_equal(dwaring(h$x, h$mu, h$sigma, log = TRUE), logd)
  expect_equal(pwaring(h$x, h$mu, h$sigma), cdf)
})

test_that("waring mass is zero below zero", {
  expect_equal(dwaring(c(-2, -0.5), 2, 0.5), c(0, 0))
  expect_equal(pwaring(c(-2, -0.5), 2, 0.5), c(0, 0))
  expect_equal(pwaring(0, 2, 0.5), dwaring(0, 2, 0.5))
})

test_that("pwaring is the cumulative sum of dwaring", {
  expect_equal(pwaring(0:60, 2, 0.5), cumsum(dwaring(0:60, 2, 0.5)))
  expect_equal(pwaring(0:60, 7, 1.3), cumsum(dwaring(0:60, 7, 1.3)))
})

test_that("waring is the beta-negative binomial with nu = 1", {
  k <- 0:60
  expect_equal(dwaring(k, 2, 0.5), dbnbinom2(k, 2, 0.5, nu = 1))
  expect_equal(dwaring(k, 2, 0.5),
               dbnbinom(k, size = 1, shape1 = 1 / 0.5 + 1, shape2 = 2 / 0.5))
})

test_that("mu is exactly the mean and the variance matches the closed form", {
  k <- 0:3000000
  for (p in list(c(2, 0.5), c(6, 0.25))) {
    mu <- p[1]; sigma <- p[2]
    d <- dwaring(k, mu, sigma); m <- sum(k * d)
    expect_equal(m, mu, tolerance = 1e-6)
    expect_equal(sum((k - m)^2 * d), mu * (sigma + 1) * (mu + 1) / (1 - sigma),
                 tolerance = 1e-5)
  }
})

test_that("rwaring draws follow the mass function", {
  set.seed(4)
  x <- rwaring(2e4, 2, 0.5)
  expect_true(all(x == floor(x)) && min(x) >= 0)
  e <- c(dwaring(0:15, 2, 0.5), 1 - sum(dwaring(0:15, 2, 0.5)))
  o <- as.vector(table(factor(pmin(x, 16), levels = 0:16)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
})

test_that("dwaring rejects non-positive parameters", {
  expect_error(dwaring(1, 0, 0.5), "mu")
  expect_error(dwaring(1, 2, -1), "sigma")
  expect_error(pwaring(1, 2, 0), "sigma")
})
