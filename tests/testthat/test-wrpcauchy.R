# Tests for the wrapped Cauchy distribution (circular, support [-pi, pi])
# Only a density function exists; p and q are not implemented.

test_that("wrpcauchy log=TRUE is consistent with log(density) (mu=0, rho=0.5)", {
  xs <- c(-2, -1, 0, 1, 2)
  expect_equal(
    dwrpcauchy(xs, mu = 0, rho = 0.5, log = TRUE),
    log(dwrpcauchy(xs, mu = 0, rho = 0.5)),
    tolerance = 1e-10
  )
})

test_that("wrpcauchy density integrates to 1 (mu=0, rho=0.5)", {
  result <- integrate(dwrpcauchy, lower = -pi, upper = pi, mu = 0, rho = 0.5)
  expect_equal(result$value, 1, tolerance = 1e-4)
})

test_that("wrpcauchy log=TRUE is consistent with log(density) (mu=1, rho=0.8)", {
  xs <- c(-2, -0.5, 1, 1.5, 2.5)
  expect_equal(
    dwrpcauchy(xs, mu = 1, rho = 0.8, log = TRUE),
    log(dwrpcauchy(xs, mu = 1, rho = 0.8)),
    tolerance = 1e-10
  )
})

test_that("wrpcauchy density integrates to 1 (mu=1, rho=0.8)", {
  result <- integrate(dwrpcauchy, lower = -pi, upper = pi, mu = 1, rho = 0.8)
  expect_equal(result$value, 1, tolerance = 1e-4)
})
