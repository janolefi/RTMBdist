# Tests for the vmf2 parameterisation of the von Mises-Fisher distribution
# theta = mu * kappa (encodes both direction and concentration)
# x is a matrix of unit vectors

test_that("vmf2 log=TRUE is consistent with log(density) (3D, theta=c(2,0,0))", {
  s <- 1 / sqrt(3)
  x <- rbind(c(1, 0, 0), c(0, 1, 0), c(s, s, s), c(s, -s, s))
  theta <- c(2, 0, 0)  # mu = (1,0,0), kappa = 2
  expect_equal(
    dvmf2(x, theta = theta, log = TRUE),
    log(dvmf2(x, theta = theta)),
    tolerance = 1e-10
  )
})

test_that("vmf2 log=TRUE is consistent with log(density) (2D, theta=c(3,4))", {
  x <- rbind(c(1, 0), c(0, 1), c(-1, 0), c(1/sqrt(2), 1/sqrt(2)))
  theta <- c(3, 4)  # kappa = 5, mu = (0.6, 0.8)
  expect_equal(
    dvmf2(x, theta = theta, log = TRUE),
    log(dvmf2(x, theta = theta)),
    tolerance = 1e-10
  )
})
