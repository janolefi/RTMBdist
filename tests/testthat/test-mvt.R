# Tests for the multivariate t distribution
# x is a matrix with one observation per row

test_that("mvt log=TRUE is consistent with log(density) (2D, mu=c(0,0), Sigma=I, df=5)", {
  x <- rbind(c(-1,  1),
             c( 0,  2),
             c( 1, -1),
             c( 2,  0))
  mu    <- c(0, 0)
  Sigma <- diag(2)
  expect_equal(
    dmvt(x, mu = mu, Sigma = Sigma, df = 5, log = TRUE),
    log(dmvt(x, mu = mu, Sigma = Sigma, df = 5)),
    tolerance = 1e-10
  )
})

test_that("mvt log=TRUE is consistent with log(density) (3D, mu=c(1,2,3), df=10)", {
  x <- rbind(c( 0, 1, 2),
             c( 1, 2, 3),
             c( 2, 3, 4),
             c(-1, 0, 1))
  mu    <- c(1, 2, 3)
  Sigma <- matrix(c(2, 0.5, 0.3,
                    0.5, 1, 0.2,
                    0.3, 0.2, 1.5), nrow = 3)
  expect_equal(
    dmvt(x, mu = mu, Sigma = Sigma, df = 10, log = TRUE),
    log(dmvt(x, mu = mu, Sigma = Sigma, df = 10)),
    tolerance = 1e-10
  )
})
