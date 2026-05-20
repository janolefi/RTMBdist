# Tests for the von Mises-Fisher distribution
# x and mu are matrices of unit vectors (each row has norm 1)

test_that("vmf log=TRUE is consistent with log(density) (3D, kappa=2)", {
  s <- 1 / sqrt(3)
  x <- rbind(c(1, 0, 0), c(0, 1, 0), c(s, s, s), c(s, -s, s))
  mu <- c(1, 0, 0)
  expect_equal(
    dvmf(x, mu = mu, kappa = 2, log = TRUE),
    log(dvmf(x, mu = mu, kappa = 2)),
    tolerance = 1e-10
  )
})

test_that("vmf log=TRUE is consistent with log(density) (2D, kappa=5)", {
  x <- rbind(c(1, 0), c(0, 1), c(-1, 0), c(1/sqrt(2), 1/sqrt(2)))
  mu <- c(1/sqrt(2), 1/sqrt(2))
  expect_equal(
    dvmf(x, mu = mu, kappa = 5, log = TRUE),
    log(dvmf(x, mu = mu, kappa = 5)),
    tolerance = 1e-10
  )
})
