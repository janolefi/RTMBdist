# Tests for the half-Cauchy distribution

test_that("halfcauchy passes standard distribution checks", {
  check_continuous_dist(
    dfun = dhalfcauchy, pfun = phalfcauchy, qfun = qhalfcauchy,
    xs = c(0.1, 0.7, 2, 9, 300), lower = 0, upper = Inf,
    sigma = 1
  )
  check_continuous_dist(
    dfun = dhalfcauchy, pfun = phalfcauchy, qfun = qhalfcauchy,
    xs = c(0.5, 3, 12, 80), lower = 0, upper = Inf,
    sigma = 5
  )
})

test_that("halfcauchy AD gradient has no NaN", {
  check_ad_gradient(dhalfcauchy, rhalfcauchy, sigma = 1)
  check_ad_gradient(dhalfcauchy, rhalfcauchy, sigma = 4)
})

test_that("dhalfcauchy and phalfcauchy match their closed forms", {
  xs <- c(0, 0.3, 1, 4, 1000)
  for (s in c(0.3, 1, 20)) {
    expect_equal(dhalfcauchy(xs, s), 2 / (pi * s * (1 + (xs / s)^2)))
    expect_equal(phalfcauchy(xs, s), (2 / pi) * atan(xs / s))
    expect_equal(qhalfcauchy(c(0.1, 0.5, 0.9), s), s * tan(pi * c(0.1, 0.5, 0.9) / 2))
  }
})

test_that("halfcauchy mass sits on the closed positive half line", {
  expect_equal(dhalfcauchy(c(-4, -1, -1e-9), 2), c(0, 0, 0))
  expect_equal(dhalfcauchy(c(-4, -1), 2, log = TRUE), c(-Inf, -Inf))
  expect_equal(phalfcauchy(c(-4, -1), 2), c(0, 0))
  expect_equal(dhalfcauchy(0, 2), 1 / pi) # 2 / (pi * sigma)
  expect_equal(phalfcauchy(0, 2), 0)
  expect_equal(qhalfcauchy(0, 2), 0)
  # tan(pi / 2) is only 1.6e16 in double precision, so this is returned explicitly
  expect_equal(qhalfcauchy(1, 2), Inf)
})

test_that("halfcauchy is the half-t with one degree of freedom", {
  xs <- c(0, 0.2, 1, 3, 10)
  expect_equal(dhalfcauchy(xs, 2), dhalft(xs, 1, 2))
  expect_equal(phalfcauchy(xs, 2), phalft(xs, 1, 2))
  expect_equal(qhalfcauchy(c(0.1, 0.5, 0.9), 2), qhalft(c(0.1, 0.5, 0.9), 1, 2))
})

test_that("qhalfcauchy honours lower.tail and log.p and recycles", {
  p <- c(0.01, 0.1, 0.5, 0.9, 0.99)
  expect_equal(qhalfcauchy(p, 2), qhalfcauchy(1 - p, 2, lower.tail = FALSE))
  expect_equal(qhalfcauchy(p, 2), qhalfcauchy(log(p), 2, log.p = TRUE))
  expect_length(qhalfcauchy(0.5, sigma = 1:4), 4)
  expect_length(dhalfcauchy(1:6, sigma = c(1, 2)), 6)
})

test_that("dhalfcauchy rejects a non-positive scale", {
  expect_error(dhalfcauchy(1, 0), "sigma")
  expect_error(phalfcauchy(1, -1), "sigma")
})
