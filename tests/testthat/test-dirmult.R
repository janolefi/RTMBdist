# Tests for the Dirichlet-multinomial distribution
# x is a matrix of non-negative integer count vectors; rows must sum to size

test_that("dirmult log=TRUE is consistent with log(density) (size=20, alpha=c(2,3,5))", {
  x <- rbind(c(4,  6, 10),
             c(2, 10,  8),
             c(8,  8,  4))
  expect_equal(
    ddirmult(x, size = 20, alpha = c(2, 3, 5), log = TRUE),
    log(ddirmult(x, size = 20, alpha = c(2, 3, 5))),
    tolerance = 1e-10
  )
})

test_that("dirmult log=TRUE is consistent with log(density) (size=10, alpha=c(1,1,1,1))", {
  x <- rbind(c(1, 2, 3, 4),
             c(3, 3, 2, 2),
             c(5, 2, 2, 1))
  expect_equal(
    ddirmult(x, size = 10, alpha = c(1, 1, 1, 1), log = TRUE),
    log(ddirmult(x, size = 10, alpha = c(1, 1, 1, 1))),
    tolerance = 1e-10
  )
})
