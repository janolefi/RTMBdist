test_that("rgmrf survives an indefinite precision matrix", {
  Qbad <- Matrix::Diagonal(5)
  Qbad[1, 1] <- -1
  expect_warning(s <- rgmrf(3, mean = 0, Q = Qbad), "not PD")
  expect_true(all(is.finite(s)))
  expect_equal(dim(s), c(3, 5))
})
