# Tests for the Wishart distribution
# x is a positive definite matrix, or a p x p x n array for multiple observations

test_that("wishart log=TRUE is consistent with log(density) (2x2 matrix, nu=5, Sigma=I)", {
  nu    <- 5
  Sigma <- diag(2)
  # single positive definite matrix
  x <- matrix(c(3, 0.5, 0.5, 2), nrow = 2)
  expect_equal(
    dwishart(x, nu = nu, Sigma = Sigma, log = TRUE),
    log(dwishart(x, nu = nu, Sigma = Sigma)),
    tolerance = 1e-10
  )
})

test_that("wishart log=TRUE is consistent with log(density) (3D array of matrices, nu=6)", {
  nu    <- 6
  Sigma <- matrix(c(2, 0.4, 0.4, 1), nrow = 2)
  # three positive definite 2x2 matrices stored in a 2x2x3 array
  x <- array(c(4, 0.8, 0.8, 2,
               3, 0.5, 0.5, 1.5,
               5, 1.0, 1.0, 3), dim = c(2, 2, 3))
  expect_equal(
    dwishart(x, nu = nu, Sigma = Sigma, log = TRUE),
    log(dwishart(x, nu = nu, Sigma = Sigma)),
    tolerance = 1e-10
  )
})
