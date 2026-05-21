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

test_that("dirmult AD gradient has no NaN", {
  set.seed(42)
  alpha_true <- c(2, 3, 5)
  size <- 20
  x <- rdirmult(100, size, alpha_true)
  nll <- function(par) {
    -sum(ddirmult(x, size = size, alpha = par, log = TRUE))
  }
  F <- tryCatch(
    RTMB::MakeTape(nll, alpha_true),
    error = function(e) { fail(paste("MakeTape failed:", conditionMessage(e))); NULL }
  )
  if (!is.null(F)) {
    grad <- F$jacobian(alpha_true)
    expect_false(any(is.nan(grad)), label = "no NaN in AD gradient")
  }
})
