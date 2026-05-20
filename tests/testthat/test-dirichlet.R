# Tests for the Dirichlet distribution
# x is a matrix with rows on the simplex (each row sums to 1)

test_that("dirichlet log=TRUE is consistent with log(density) (alpha=c(2,3,5))", {
  x <- rbind(c(0.2, 0.3, 0.5),
             c(0.1, 0.5, 0.4),
             c(0.4, 0.4, 0.2))
  alpha <- c(2, 3, 5)
  expect_equal(
    ddirichlet(x, alpha, log = TRUE),
    log(ddirichlet(x, alpha)),
    tolerance = 1e-10
  )
})

test_that("dirichlet log=TRUE is consistent with log(density) (alpha=c(0.5,0.5,0.5))", {
  x <- rbind(c(0.05, 0.15, 0.80),
             c(0.33, 0.33, 0.34),
             c(0.6,  0.2,  0.2))
  alpha <- c(0.5, 0.5, 0.5)
  expect_equal(
    ddirichlet(x, alpha, log = TRUE),
    log(ddirichlet(x, alpha)),
    tolerance = 1e-10
  )
})
