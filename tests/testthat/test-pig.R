# Tests for the Poisson-inverse Gaussian distribution
# dpig is dsichel with nu = -1/2, so large counts, gradients, recycling and parameter
# recovery are tested there

test_that("pig passes standard discrete checks", {
  check_discrete_dist(
    dfun        = dpig,
    pfun        = ppig,
    xs_int      = c(0, 1, 4, 10, 30),
    sum_support = 0:2000,
    mu = 4, sigma = 0.5
  )
  check_discrete_dist(
    dfun        = dpig,
    pfun        = ppig,
    xs_int      = c(0, 2, 8, 25, 80),
    sum_support = 0:20000,
    mu = 5, sigma = 2
  )
})

test_that("pig AD gradient has no NaN", {
  check_ad_gradient(dpig, rpig, mu = 4, sigma = 0.5)
  check_ad_gradient(dpig, rpig, mu = 20, sigma = 3)
})

test_that("dpig and ppig match gamlss.dist", {
  # reference values from gamlss.dist::dPIG / pPIG, baked in: gamlss.dist
  # is scheduled for archival, so the test must not depend on it
  h <- expand.grid(x = c(0, 1, 4, 15, 45, 120), mu = c(0.5, 4, 25), sigma = c(0.05, 0.6, 3))
  logd <- c(
    -0.493901531919199, -1.21144379456386, -6.2567748187826, -34.3179768983386,
    -123.300533898765, -351.648949572045, -3.66431913239846, -2.44626088958918,
    -1.72213690009475, -9.05806386413568, -43.2850193256362, -136.786777742091,
    -17.4165738677394, -14.8240795271189, -10.0664871035362, -3.68992656063948,
    -5.9165050121868, -28.1158331890106, -0.441518440112253, -1.36966743529507,
    -5.15492498169116, -17.6839984700507, -48.7174676734076, -123.740008377245,
    -2.34719819293077, -1.83983279058706, -2.19892642511174, -5.37913683162147,
    -12.5190710406804, -28.1302195089493, -7.6129406047167, -6.11105838209108,
    -4.17394767393998, -3.40719262242949, -5.0402212705624, -8.66543773693638,
    -0.333333333333333, -1.71962769445322, -4.50928840059165, -9.69760827283077,
    -19.987489514943, -43.0385341416545, -1.33333333333333, -1.55647688464754,
    -2.82071881842926, -5.12152373688534, -7.97590231780616, -12.5039565597706,
    -3.7627352424815, -3.05249933602076, -2.82528555923557, -3.96529505909846,
    -5.61906870161597, -7.53299101943776)
  cdf <- c(
    0.610240869379575, 0.908007923851427, 0.999748173872166, 0.999999999999998,
    0.999999999999998, 0.999999999999998, 0.0256216102705508, 0.112238462003158,
    0.630302416658077, 0.999930055471762, 1, 1, 2.72946843071263e-08,
    3.92035245583633e-07, 5.71227729389601e-05, 0.0862692710852453,
    0.989996440893919, 0.999999999998408, 0.643059232436591, 0.897250713048999,
    0.997598672276412, 0.999999989109927, 1, 1, 0.095636742888097,
    0.254480726983912, 0.671395237974569, 0.984313937960702, 0.999984995230449,
    0.999999999997272, 0.000494017011706373, 0.00271221887655357,
    0.0337432633916745, 0.383484418234044, 0.876545593723864, 0.996011152447591,
    0.716531310573789, 0.895664138217236, 0.985507531508786, 0.999864813928481,
    0.999999994446557, 1, 0.263597138115727, 0.474474848608308, 0.756274950693755,
    0.946310360110749, 0.995010631512861, 0.999930062971895, 0.0232201407189876,
    0.0704608470393736, 0.251823723500906, 0.616164185928754, 0.857473055923175,
    0.964712221747896)
  expect_equal(dpig(h$x, h$mu, h$sigma, log = TRUE), logd, tolerance = 1e-10)
  expect_equal(ppig(h$x, h$mu, h$sigma), cdf, tolerance = 1e-10)
})

test_that("mu is the mean and the variance is mu + sigma * mu^2", {
  k <- 0:200000
  for (p in list(c(4, 0.5), c(10, 2))) {
    mu <- p[1]; sigma <- p[2]
    d <- dpig(k, mu, sigma)
    m <- sum(k * d)
    expect_equal(m, mu, tolerance = 1e-8)
    expect_equal(sum((k - m)^2 * d), mu + sigma * mu^2, tolerance = 1e-6)
  }
})

test_that("pig collapses to the Poisson as sigma goes to zero", {
  k <- 0:40
  expect_equal(dpig(k, 6, 1e-9), stats::dpois(k, 6), tolerance = 1e-7)
})

test_that("rpig draws follow the mass function", {
  set.seed(3)
  x <- rpig(2e4, 4, 0.5)
  expect_true(all(x == floor(x)) && min(x) >= 0)
  e <- c(dpig(0:14, 4, 0.5), 1 - sum(dpig(0:14, 4, 0.5)))
  o <- as.vector(table(factor(pmin(x, 15), levels = 0:15)))
  expect_gt(suppressWarnings(stats::chisq.test(o, p = e)$p.value), 0.01)
})

test_that("dpig is zero below the support", {
  expect_equal(dpig(c(-3, -1), 4, 0.5), c(0, 0))
  expect_equal(ppig(-1, 4, 0.5), 0)
})

test_that("dpig rejects non-positive parameters", {
  expect_error(dpig(1, 0, 0.5), "mu")
  expect_error(dpig(1, 4, -1), "sigma")
  expect_error(rpig(1, 4, 0), "sigma")
})

test_that("ppig is AD-compatible and gives OSA residuals", {
  q <- c(0, 1, 4, 10, 30)
  check_ad_cdf(ppig, dpig, q, mu = 4, sigma = 0.5, .dq = FALSE)
  check_ad_cdf_hessian(ppig, q, mu = 4, sigma = 0.5)
  check_osa_cdf(dpig, ppig, q, mu = 4, sigma = 0.5, .discrete = TRUE)
})
