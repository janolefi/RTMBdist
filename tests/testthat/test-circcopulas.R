# Tests for the circular-linear copulas cjw() and cfold()
# u is the circular, v the linear margin

# integrates exp(logc) over u for fixed v and over v for fixed u
copula_margins <- function(logc, at = c(0.05, 0.3, 0.5, 0.77, 0.95)) {
  cu <- sapply(at, function(v)
    integrate(function(u) exp(logc(u, rep(v, length(u)))), 0, 1, rel.tol = 1e-10)$value)
  cv <- sapply(at, function(u)
    integrate(function(v) exp(logc(rep(u, length(v)), v)), 0, 1, rel.tol = 1e-10)$value)
  c(cu, cv)
}

test_that("cjw has uniform margins for von Mises and wrapped Cauchy binding densities", {
  for (q in c(1, -1)) {
    expect_equal(copula_margins(cjw(dvm, mu = 0.5, kappa = 3, q = q)), rep(1, 10), tolerance = 1e-8)
    expect_equal(copula_margins(cjw(dwrpcauchy, mu = -1, rho = 0.6, q = q)), rep(1, 10), tolerance = 1e-8)
  }
})

test_that("cjw matches the Johnson-Wehrly density", {
  u <- c(0.1, 0.4, 0.9); v <- c(0.2, 0.7, 0.5)
  expect_equal(cjw(dvm, mu = 1, kappa = 2)(u, v),
               log(2 * pi) + dvm(2 * pi * (u - v), 1, 2, log = TRUE))
  expect_equal(cjw(dwrpcauchy, mu = 1, rho = 0.4, q = -1)(u, v),
               log(2 * pi) + dwrpcauchy(2 * pi * (u + v), 1, 0.4, log = TRUE))
  # zero concentration gives the independence copula
  expect_equal(cjw(dvm, mu = 0, kappa = 0)(u, v), rep(0, 3))
  expect_equal(cjw(dwrpcauchy, mu = 0, rho = 0)(u, v), rep(0, 3))
  # a function name works as well as the function
  expect_equal(cjw("dvm", mu = 1, kappa = 2)(u, v), cjw(dvm, mu = 1, kappa = 2)(u, v))
})

test_that("cjw with dcopula reproduces the joint density of the simulation recipe", {
  # f(theta, s) = 2 pi g(2 pi (F1(theta) - F2(s))) f1(theta) f2(s)
  angle <- c(-2, -0.3, 0.1, 1.5); step <- c(0.2, 1.1, 2.4, 0.6)
  d1 <- dwrpcauchy(angle, 0, 0.5, log = TRUE); p1 <- pwrpcauchy(angle, 0, 0.5)
  d2 <- dweibull(step, 2, 1, log = TRUE); p2 <- pweibull(step, 2, 1)
  expect_equal(
    dcopula(d1, d2, p1, p2, copula = cjw(dwrpcauchy, mu = 0.3, rho = 0.7), log = TRUE),
    log(2 * pi) + dwrpcauchy(2 * pi * (p1 - p2), 0.3, 0.7, log = TRUE) + d1 + d2
  )
})

test_that("cjw simulation recipe has uniform margins", {
  set.seed(1)
  n <- 20000
  v <- runif(n)
  z <- rvm(n, mu = 1, kappa = 4)
  u <- (z / (2 * pi) - v) %% 1
  expect_gt(ks.test(u, "punif")$p.value, 0.01)
  # the density of (u, v) is the copula: z is recovered from u and v
  expect_equal(abs(mean(exp(1i * 2 * pi * (u + v)))), besselI(4, 1) / besselI(4, 0), tolerance = 0.02)
})

test_that("cjw rejects q other than 1 and -1", {
  expect_error(cjw(dvm, mu = 0, kappa = 1, q = 2), "q must")
  expect_error(cjw(dvm, mu = 0, kappa = 1, q = 0.5), "q must")
  expect_error(cjw(dvm, mu = 0, kappa = 1, q = c(1, -1)), "q must")
})

test_that("cfold is a copula and symmetric in the sign of the angle", {
  # at u = 1/2 the folded value is 1, where the linear copula is degenerate in
  # v, so the margins are checked away from it
  for (cop in list(cgaussian(0.6), cgaussian(-0.4), cclayton(2), cgumbel(1.5), cfrank(-3))) {
    folded <- cfold(cop)
    expect_equal(copula_margins(folded, at = c(0.05, 0.3, 0.43, 0.77, 0.95)),
                 rep(1, 10), tolerance = 1e-6)
    u <- c(0.05, 0.2, 0.45); v <- c(0.1, 0.5, 0.9)
    expect_equal(folded(u, v), folded(1 - u, v))
    expect_equal(folded(u, v), cop(2 * u, v))
  }
})

test_that("cfold stays finite at the mean direction and its antipode", {
  folded <- cfold(cgaussian(0.6))
  expect_true(all(is.finite(folded(c(0.5, 1e-17, 1 - 1e-17), c(0.3, 0.3, 0.3)))))
})

test_that("cfold with positive dependence makes long steps straight", {
  # with a margin centred at 0 and the default origin, 1 - |2u - 1| is large
  # for small |angle|, so long steps should be more likely to be straight
  folded <- cfold(cgaussian(0.7))
  u_straight <- pwrpcauchy(0.1, 0, 0.5); u_turn <- pwrpcauchy(2.5, 0, 0.5)
  expect_gt(folded(u_straight, 0.95), folded(u_turn, 0.95))
  expect_lt(folded(u_straight, 0.05), folded(u_turn, 0.05))
  expect_equal(folded(u_turn, 0.3), folded(pwrpcauchy(-2.5, 0, 0.5), 0.3))
})

test_that("circular-linear copula likelihoods have AD gradients without NaN", {
  set.seed(1)
  n <- 50
  v <- runif(n)
  u <- (rwrpcauchy(n, 0, 0.6) / (2 * pi) + v) %% 1
  angle <- qwrpcauchy(u, 0, 0.5); step <- qweibull(v, 2, 1)
  angle[1] <- 0 # exactly at the mean direction, where cfold has its kink

  nll_jw <- function(par) {
    d1 <- dwrpcauchy(angle, par[1], plogis(par[2]), log = TRUE)
    p1 <- pwrpcauchy(angle, par[1], plogis(par[2]))
    d2 <- dweibull(step, exp(par[3]), exp(par[4]), log = TRUE)
    p2 <- 1 - exp(-(step / exp(par[4]))^exp(par[3]))
    -sum(dcopula(d1, d2, p1, p2, copula = cjw(dvm, mu = par[5], kappa = exp(par[6])), log = TRUE))
  }
  par <- c(0.1, 0, log(2), 0, 0.2, log(2))
  F <- RTMB::MakeTape(nll_jw, par)
  expect_false(any(is.nan(F$jacobian(par))))
  expect_equal(F(par), nll_jw(par))

  nll_fold <- function(par) {
    d1 <- dwrpcauchy(angle, 0, plogis(par[1]), log = TRUE)
    p1 <- pwrpcauchy(angle, 0, plogis(par[1]))
    d2 <- dweibull(step, 2, 1, log = TRUE)
    p2 <- 1 - exp(-step^2)
    -sum(dcopula(d1, d2, p1, p2, copula = cfold(cgaussian(tanh(par[2]))), log = TRUE))
  }
  par <- c(0, 0.5)
  G <- RTMB::MakeTape(nll_fold, par)
  expect_false(any(is.nan(G$jacobian(par))))
  expect_equal(G(par), nll_fold(par))
})

test_that("cjw recovers the parameters of the binding density", {
  set.seed(2)
  n <- 3000
  v <- runif(n)
  u <- (rvm(n, mu = 0.8, kappa = 3) / (2 * pi) + v) %% 1
  angle <- qwrpcauchy(u, 0, 0.5); step <- qweibull(v, 2, 1)

  nll <- function(pars) {
    par <- pars$par
    d1 <- dwrpcauchy(angle, 0, plogis(par[1]), log = TRUE)
    p1 <- pwrpcauchy(angle, 0, plogis(par[1]))
    d2 <- dweibull(step, exp(par[2]), 1, log = TRUE)
    p2 <- 1 - exp(-step^exp(par[2]))
    -sum(dcopula(d1, d2, p1, p2, copula = cjw(dvm, mu = par[3], kappa = exp(par[4])), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll, list(par = c(0, 0, 0, 0)), silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  est <- c(plogis(opt$par[1]), exp(opt$par[2]), opt$par[3], exp(opt$par[4]))
  expect_equal(unname(est), c(0.5, 2, 0.8, 3), tolerance = 0.1)
})
