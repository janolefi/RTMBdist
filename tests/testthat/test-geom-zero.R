# Tests for the zero-modified geometric distributions.
# The base geometric is RTMB's negative binomial with size = 1; RTMBdist does not
# export d/p/q/rgeom, so stats' versions stay visible to users.

test_that("the geometric variants agree with the negative binomial at size = 1", {
  x <- 0:200
  expect_equal(dztgeom(x, 0.3),          dztnbinom(x, 1, 0.3))
  expect_equal(dzigeom(x, 0.3, 0.4),     dzinbinom(x, 1, 0.3, 0.4))
  expect_equal(dhgeom(x, 0.3, 0.4),      dhnbinom(x, 1, 0.3, 0.4))
})

test_that("zigeom, ztgeom and hgeom pass discrete distribution checks", {
  check_discrete_dist(dfun = dzigeom, pfun = pzigeom, xs_int = 0:10,
                      sum_support = 0:5000, prob = 0.3, zeroprob = 0.4)
  check_discrete_dist(dfun = dztgeom, pfun = pztgeom, xs_int = 1:10,
                      sum_support = 0:5000, prob = 0.3)
  check_discrete_dist(dfun = dhgeom, pfun = phgeom, xs_int = 0:10,
                      sum_support = 0:5000, prob = 0.3, zeroprob = 0.4)
})

test_that("the geometric variants satisfy their defining identities", {
  prob <- 0.3; zp <- 0.4; x <- 1:200
  expect_equal(dzigeom(0, prob, zp), zp + (1 - zp) * prob)
  expect_equal(dzigeom(x, prob, zp), (1 - zp) * dgeom(x, prob))
  expect_equal(dztgeom(x, prob),     dgeom(x, prob) / (1 - prob))
  expect_equal(dhgeom(0, prob, zp),  zp)
  expect_equal(dhgeom(x, prob, zp),  (1 - zp) * dztgeom(x, prob))
  # a hurdle with zeroprob = P(X = 0) is the ordinary geometric
  expect_equal(dhgeom(0:200, prob, prob), dgeom(0:200, prob))
  # a hurdle with zeroprob = 0 is the zero-truncated geometric
  expect_equal(dhgeom(0:200, prob, 0), dztgeom(0:200, prob))
})

test_that("geometric variants have AD gradients without NaN", {
  check_ad_gradient(dzigeom, rzigeom, prob = 0.3, zeroprob = 0.4)
  check_ad_gradient(dztgeom, rztgeom, prob = 0.3)
  check_ad_gradient(dhgeom,  rhgeom,  prob = 0.3, zeroprob = 0.4)
})

test_that("geometric variant RNGs have the right zero behaviour", {
  set.seed(42)
  expect_equal(mean(rzigeom(2e5, 0.3, 0.4) == 0), 0.4 + 0.6 * 0.3, tolerance = 0.02)
  expect_true(all(rztgeom(1e4, 0.3) >= 1))
  y <- rhgeom(2e5, 0.3, 0.4)
  expect_equal(mean(y == 0), 0.4, tolerance = 0.02)
  expect_equal(mean(y[y > 0]), sum((1:5000) * dztgeom(1:5000, 0.3)), tolerance = 0.05)
})
