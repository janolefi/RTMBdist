# Tests for mcreport()

make_obj <- function() {
  set.seed(1)
  x <- rnorm(200, 1, 2)
  nll <- function(par) {
    getAll(par)
    -sum(dnorm(x, mu, exp(log_sigma), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll, list(mu = 0, log_sigma = 0), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  obj
}

test_that("supplied Q equal to the Hessian reproduces default samples without random effects", {
  obj <- make_obj()
  H <- obj$he(obj$env$last.par.best)
  set.seed(2); s1 <- suppressMessages(mcreport(obj, 50))
  set.seed(2); s2 <- mcreport(obj, 50, Q = H)
  expect_equal(s1, s2)
})

test_that("augmented Q returns parameter samples with unchanged mean and covariance", {
  obj <- make_obj()
  p_hat <- obj$env$last.par.best
  H <- obj$he(p_hat)
  Q <- Matrix::bdiag(H, Matrix::Matrix(4))
  set.seed(3); s <- mcreport(obj, 1e5, Q = Q, par = c(p_hat, 10))
  expect_named(s, c("mu", "log_sigma"))
  S <- cbind(unlist(s$mu), unlist(s$log_sigma))
  expect_equal(colMeans(S), unname(p_hat), tolerance = 0.01)
  expect_equal(cov(S), solve(H), tolerance = 0.02, ignore_attr = TRUE)
  expect_error(mcreport(obj, 10, Q = Q), "par must have length")
  expect_error(mcreport(obj, 10, Q = H[1, 1, drop = FALSE]), "at least 2 rows")
})

test_that("par overrides the mean", {
  obj <- make_obj()
  H <- obj$he(obj$env$last.par.best)
  s <- mcreport(obj, 1e4, Q = H, par = c(5, -1))
  expect_equal(mean(unlist(s$mu)), 5, tolerance = 0.01)
  expect_equal(mean(unlist(s$log_sigma)), -1, tolerance = 0.01)
})
