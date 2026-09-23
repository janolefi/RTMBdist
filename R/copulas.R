#' Joint density under a bivariate copula
#'
#' Computes the joint density (or log-density) of a bivariate distribution
#' constructed from two arbitrary margins combined with a specified copula.
#'
#' @param d1,d2 Marginal density values. If \code{log = TRUE}, supply the log-density.
#' If \code{log = FALSE}, supply the raw density.
#' @param p1,p2 Marginal CDF values. Need not be supplied on log scale.
#' @param copula A function of two arguments (\code{u}, \code{v}) returning
#'   the log copula density \eqn{\log c(u,v)}.
#'   You can either construct this yourself or use the copula constructors available (see details)
#' @param log Logical; if \code{TRUE}, return the log joint density. In this case,
#'   \code{d1} and \code{d2} must be on the log scale.
#'
#' @returns Joint density (or log-density) under the bivariate copula.
#'
#' @details
#' The joint density is
#' \deqn{f(x,y) = c(F_1(x), F_2(y)) \, f_1(x) f_2(y),}
#' where \eqn{F_i} are the marginal CDFs, \eqn{f_i} are the marginal densities,
#' and \eqn{c} is the copula density.
#'
#' The marginal densities \code{d1}, \code{d2} and CDFs \code{p1}, \code{p2}
#' must be differentiable for automatic differentiation (AD) to work.
#'
#' Available copula constructors are:
#' - \code{\link{cgaussian}} (Gaussian copula)
#' - \code{\link{cclayton}} (Clayton copula)
#' - \code{\link{cgumbel}} (Gumbel copula)
#' - \code{\link{cfrank}} (Frank copula)
#'
#' Circular-linear copulas, e.g. for turning angles and step lengths, with the
#' circular margin as the first one:
#' - \code{\link{cjw}} (Johnson-Wehrly copula with a circular binding density such as \code{\link{dvm}} or \code{\link{dwrpcauchy}})
#' - \code{\link{cfold}} (folds any of the copulas above into one that is symmetric in the sign of the angle)
#'
#' @seealso [ddcopula()], [dmvcopula()]
#'
#' @export
#'
#' @examples
#' # Normal + Exponential margins with Gaussian copula
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
#'
#' # Normal + Beta margins with Clayton copula
#' x <- c(0.5, 1); y <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
#'
#' # Normal + Beta margins with Gumbel copula
#' x <- c(0.5, 1); y <- c(0.2, 0.4)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
#'
#' # Normal + Exponential margins with Frank copula
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
dcopula <- function(d1, d2, p1, p2, copula = cgaussian(0), log = FALSE) {

  # ensure numeric stability of uniforms
  eps <- .Machine$double.eps
  p1 <- pmin.ad(pmax.ad(p1, eps), 1 - eps)
  p2 <- pmin.ad(pmax.ad(p2, eps), 1 - eps)

  # evaluate log copula density
  logC <- copula(p1, p2)

  # combine with marginal densities
  if (log) {
    logC + d1 + d2
  } else {
    exp(logC) * d1 * d2
  }
}

#' Joint probability under a discrete bivariate copula
#'
#' Computes the joint probability mass function of two \strong{discrete} margins
#' combined with a copula CDF.
#'
#' @param d1,d2 Marginal p.m.f. values at the observed points, \strong{not} on log-scale.
#' @param p1,p2 Marginal CDF values at the observed points.
#' @param Copula A function of two arguments returning the copula CDF.
#' @param log Logical; if \code{TRUE}, return the log joint density. In this case,
#'   \code{d1} and \code{d2} must be on the log scale.
#'
#' @details
#' The joint probability mass function for two discrete margins is
#' \deqn{
#' \Pr(Y_1 = y_1, Y_2 = y_2) =
#' C(F_1(y_1), F_2(y_2))
#' - C(F_1(y_1-1), F_2(y_2))
#' - C(F_1(y_1), F_2(y_2-1))
#' + C(F_1(y_1-1), F_2(y_2-1)),
#' }
#' where \eqn{F_i} are the marginal CDFs, and \eqn{C} is the copula CDF.
#'
#' Available copula CDF constructors are:
#' - \code{\link{Cclayton}} (Clayton copula)
#' - \code{\link{Cgumbel}} (Gumbel copula)
#' - \code{\link{Cfrank}} (Frank copula)
#'
#' @return Joint probability (or log-probability) under chosen copula
#'
#' @seealso [dcopula()], [dmvcopula()]
#'
#' @export
#'
#' @examples
#' x <- c(3,5); y <- c(2,4)
#' d1 <- dpois(x, 4); d2 <- dpois(y, 3)
#' p1 <- ppois(x, 4); p2 <- ppois(y, 3)
#' ddcopula(d1, d2, p1, p2, Copula = Cclayton(2), log = FALSE)
ddcopula <- function(d1, d2, p1, p2, Copula, log = FALSE) {

  # ensure numeric stability of uniforms
  eps <- .Machine$double.eps
  p1 <- pmin.ad(pmax.ad(p1, eps), 1 - eps)
  p2 <- pmin.ad(pmax.ad(p2, eps), 1 - eps)
  p1m1 <- p1 - d1 # F(x_1 - 1)
  p2m1 <- p2 - d2 # F(x_2 - 1)
  p1m1 <- pmin.ad(pmax.ad(p1m1, eps), 1 - eps)
  p2m1 <- pmin.ad(pmax.ad(p2m1, eps), 1 - eps)

  # finite-difference formula
  prob <- Copula(p1, p2) -
    Copula(p1m1, p2) -
    Copula(p1, p2m1) +
    Copula(p1m1, p2m1)

  # ensure numeric stability of result
  prob <- pmin.ad(pmax.ad(prob, eps), 1-eps)

  if(log) {
    return(log(prob))
  }
  return(prob)
}

#' Gaussian copula constructor
#'
#' Returns a function computing the log density of the bivariate Gaussian copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' @param rho Correlation parameter (\eqn{-1 < rho < 1}).
#' @return Function of two arguments (u,v) returning log copula density.
#'
#' The Gaussian copula density is
#' \deqn{
#' c(u,v;\rho) = \frac{1}{\sqrt{1-\rho^2}}
#' \exp\left\{-\frac{1}{2(1-\rho^2)} (z_1^2 - 2 \rho z_1 z_2 + z_2^2) + \frac{1}{2}(z_1^2 + z_2^2) \right\},
#' }
#' where \eqn{z_1 = \Phi^{-1}(u)}, \eqn{z_2 = \Phi^{-1}(v)}, and \eqn{-1 < \rho < 1}.
#'
#' @seealso [cclayton()], [cgumbel()], [cfrank()]
#'
#' @importFrom RTMB qnorm
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cgaussian(0.5), log = TRUE)
cgaussian <- function(rho = 0) {
  function(u, v) {
    # u, v assumed clipped in dcopula()
    z1 <- qnorm(u)
    z2 <- qnorm(v)

    # compute log copula density using Cholesky-style formula
    logdet <- -0.5 * log1p(-rho^2)  # log(1 - rho^2) stable
    exponent <- -0.5 / (1 - rho^2) * (z1^2 - 2 * rho * z1 * z2 + z2^2) +
      0.5 * (z1^2 + z2^2)

    logdet + exponent
  }
}

#' Clayton copula constructors
#'
#' Construct a function that computes the log density or CDF of the bivariate Clayton copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' The Clayton copula density is
#' \deqn{
#' c(u,v;\theta) = (1+\theta) (uv)^{-(1+\theta)}
#' \left( u^{-\theta} + v^{-\theta} - 1 \right)^{-(2\theta+1)/\theta}, \quad \theta > 0.
#' }
#'
#' @seealso [cgaussian()], [cgumbel()], [cfrank()]
#'
#' @param theta Positive dependence parameter (\eqn{\theta > 0}).
#'
#' @return A function of two arguments (u,v) returning log copula \strong{density} (\code{cclayton}) or copula \strong{CDF} (\code{Cclayton}).
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cclayton(2), log = TRUE)
#'
#' # CDF version (for discrete copulas)
#' Cclayton(1.5)(0.5, 0.4)
#' @name cclayton
NULL
#' @rdname cclayton
#' @export
cclayton <- function(theta) {
  function(u, v) {
    log(theta + 1) -
      (theta + 1) * (log(u) + log(v)) -
      ((2 * theta + 1) / theta) * log(u^(-theta) + v^(-theta) - 1)
  }
}

#' @rdname cclayton
#' @export
Cclayton <- function(theta) {
  function(u, v) {
    (u^(-theta) + v^(-theta) - 1)^(-1/theta)
  }
}

#' Gumbel copula constructors
#'
#' Construct functions that compute either the log density or the CDF
#' of the bivariate Gumbel copula, intended for use with \code{\link{dcopula}}.
#'
#' The Gumbel copula density
#'
#' \deqn{
#' c(u,v;\theta) = \exp\Big[-\big((-\log u)^\theta + (-\log v)^\theta\big)^{1/\theta}\Big] \cdot h(u,v;\theta),
#' }
#' where \eqn{h(u,v;\theta)} contains the derivative terms ensuring the function is a density.
#'
#' @seealso [cgaussian()], [cclayton()], [cfrank()]
#'
#' @param theta Dependence parameter (\eqn{\theta >= 1}).
#'
#' @return A function of two arguments \code{(u, v)} returning either
#' the log copula density (\code{cgumbel}) or the copula CDF (\code{Cgumbel}).
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(0.2, 0.4)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dbeta(y, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pbeta(y, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cgumbel(1.5), log = TRUE)
#'
#' # CDF version (for discrete copulas)
#' Cgumbel(1.5)(0.5, 0.4)
#' @name cgumbel
NULL

#' @rdname cgumbel
#' @export
cgumbel <- function(theta) {
  function(u, v) {
    lu <- -log(u); lv <- -log(v)
    A <- (lu^theta + lv^theta)^(1/theta)

    logC <- -A + (theta - 1) * (log(lu) + log(lv)) +
      (-2 + 2/theta) * log(lu^theta + lv^theta) -
      log(u) - log(v)

    term <- 1 + (theta - 1) * (lu * lv)^theta / (lu^theta + lv^theta)^2
    logC + log(term)
  }
}

#' @rdname cgumbel
#' @export
Cgumbel <- function(theta) {
  function(u, v) {
    exp(-((-log(u))^theta + (-log(v))^theta)^(1/theta))
  }
}

#' Frank copula constructor
#'
#' Returns a function computing the log density of the bivariate Frank copula,
#' intended to be used with \code{\link{dcopula}}.
#'
#' The Frank copula density is
#' \deqn{
#' c(u,v;\theta) = \frac{\theta (1-e^{-\theta}) e^{-\theta(u+v)}}
#' {\left[(e^{-\theta u}-1)(e^{-\theta v}-1) + (1 - e^{-\theta}) \right]^2}, \quad \theta \ne 0.
#' }
#'
#' @seealso [cgaussian()], [cclayton()], [cgumbel()]
#'
#' @param theta Dependence parameter (\eqn{\theta \neq 0}).
#'
#' @return A function of two arguments \code{(u, v)} returning either
#' the log copula density (\code{cfrank}) or the copula CDF (\code{Cfrank}).
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2)
#' dcopula(d1, d2, p1, p2, copula = cfrank(2), log = TRUE)
#' @name cfrank
NULL

#' @rdname cfrank
#' @export
cfrank <- function(theta) {
  function(u, v) {
    # u, v assumed clipped in (0,1)
    eu <- exp(-theta * u)
    ev <- exp(-theta * v)
    et  <- exp(-theta)

    log_num <- log(theta * (1 - et)) - theta * (u + v)
    den <- (1 - et) + (eu - 1) * (ev - 1)
    den <- den * den
    log_den <- log(den)

    log_num - log_den
  }
}

#' @rdname cfrank
#' @export
Cfrank <- function(theta) {
  function(u, v) {
    eu <- exp(-theta * u)
    ev <- exp(-theta * v)
    et <- exp(-theta)

    -1/theta * log(1 + (eu - 1) * (ev - 1) / (et - 1))
  }
}


#' Joint density under a multivariate copula
#'
#' Computes the joint density (or log-density) of a distribution
#' constructed from any number of arbitrary margins combined with a specified copula.
#'
#' @param D Matrix of marginal density values of with rows corresponding to observations and columns corresponding to dimensions.
#' If \code{log = TRUE}, supply the log-densities.
#' If \code{log = FALSE}, supply the raw densities
#' @param P Matrix of marginal CDF values of the same dimension as \code{D}. Need not be supplied on log scale.
#' @param copula A function of a matrix argument \code{U} returning
#'   the log copula density \eqn{\log c(u_1, ... u_d)}. The columns of \code{U} correspond to dimensions.
#'   You can either construct this yourself or use the copula constructors available (see details)
#' @param log Logical; if \code{TRUE}, return the log joint density. In this case,
#'   \code{D} must be on the log scale.
#'
#' @returns Joint density (or log-density) under the chosen copula.
#'
#' @details
#' The joint density is
#' \deqn{f(x_1, \dots, x_d) = c(F_1(x_1), \dots, F_d(x_d)) \, f_1(x_1) \dots f_d(x_d),}
#' where \eqn{F_i} are the marginal CDFs, \eqn{f_i} are the marginal densities,
#' and \eqn{c} is the copula density.
#'
#' The marginal densities \code{d_1, \dots, d_d} and CDFs \code{p_1, \dots, p_d}
#' must be differentiable for automatic differentiation (AD) to work.
#'
#' Available multivariate copula constructors are:
#' - \code{\link{cmvgauss}} (Multivariate Gaussian copula)
#' - \code{\link{cgmrf}} (Multivariate Gaussian copula parameterised by precision (inverse correlation) matrix)
#'
#' @seealso [dcopula()], [ddcopula()]
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2); z <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE); d3 <- dbeta(z, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2); p3 <- pbeta(z, 2, 1)
#' R <- matrix(c(1,0.5,0.3,0.5,1,0.4,0.3,0.4,1), nrow = 3)
#'
#' ### Multivariate Gaussian copula
#' ## Based on correlation matrix
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cmvgauss(R), log = TRUE)
#'
#' ## Based on precision matrix
#' Q <- solve(R)
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cgmrf(Q), log = TRUE)
#'
#' ## Parameterisation inside a model
#' # using RTMB::unstructured to get a valid correlation matrix
#' library(RTMB)
#' d <- 5 # dimension
#' cor_func <- unstructured(d)
#' npar <- length(cor_func$parms())
#' R <- cor_func$corr(rep(0.1, npar))
dmvcopula <- function(D, P, copula = cmvgauss(diag(ncol(D))), log = FALSE) {

  # ensure numeric stability of uniforms
  eps <- .Machine$double.eps
  P <- apply(P, 2, function(p) pmin.ad(pmax.ad(p, eps), 1 - eps))

  # evaluate log copula density
  log_c <- copula(P)

  # combine with marginal densities
  if (log) {
    log_c + rowSums(D)
  } else {
    exp(log_c) * apply(D, 1, prod)
  }
}


#' Multivariate Gaussian copula constructor
#'
#' Returns a function computing the log density of the multivariate Gaussian copula,
#' intended to be used with \code{\link{dmvcopula}}.
#'
#' @seealso [cgmrf()]
#'
#' @param R Positive definite correlation matrix (unit diagonal)
#'
#' @return Function with matrix argument \code{U} returning log copula density.
#'
#' @export
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2); z <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE); d3 <- dbeta(z, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2); p3 <- pbeta(z, 2, 1)
#' R <- matrix(c(1,0.5,0.3,0.5,1,0.4,0.3,0.4,1), nrow = 3)
#'
#' ## Based on correlation matrix
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cmvgauss(R), log = TRUE)
#'
#' ## Based on precision matrix
#' Q <- solve(R)
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cgmrf(Q), log = TRUE)
#'
#' ## Parameterisation inside a model
#' # using RTMB::unstructured to get a valid correlation matrix
#' library(RTMB)
#' d <- 5 # dimension
#' cor_func <- unstructured(d)
#' npar <- length(cor_func$parms())
#' R <- cor_func$corr(rep(0.1, npar))
cmvgauss <- function(R) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) {
    stop("R must be a square matrix.")
  }

  if (!ad_context()) {
    if (!isTRUE(all.equal(R, t(R), tolerance = 1e-8))) {
      stop("R must be symmetric.")
    }
    if (inherits(try(chol(R), silent = TRUE), "try-error")) {
      stop("R must be positive definite.")
    }
    if (!isTRUE(all.equal(diag(R), rep(1, nrow(R)), tolerance = 1e-8))) {
      stop("R must have unit diagonal (correlation matrix).")
    }
  }

  d <- nrow(R)
  invR <- solve(R)
  logdetR <- determinant(R, logarithm = TRUE)$modulus
  P <- invR - diag(d)

  function(U) {
    if (ncol(U) != d) stop("Input dimension does not match correlation matrix R.")

    Z <- apply(U, 2, qnorm)
    quadform <- rowSums((Z %*% P) * Z)
    -0.5 * logdetR - 0.5 * quadform
  }
}


#' Multivariate Gaussian copula constructor parameterised by inverse correlation matrix
#'
#' Returns a function computing the log density of the multivariate Gaussian copula,
#' parameterised by the inverse correlation matrix.
#'
#' \strong{Caution:} Parameterising the inverse correlation directly is difficult, as inverting it needs to yield a positive definite matrix with \strong{unit diagonal}.
#' Hence we still advise parameterising the correaltion matrix \code{R} and computing its inverse.
#' This function is useful when you need access to the precision (i.e. inverse correlation) in your likelihood function.
#'
#' @seealso [cmvgauss()]
#'
#' @param Q Inverse of a positive definite correlation matrix with unit diagonal. Can either be sparse or dense matrix.
#'
#' @return Function with matrix argument \code{U} returning log copula density.
#'
#' @export
#' @import RTMB
#'
#' @examples
#' x <- c(0.5, 1); y <- c(1, 2); z <- c(0.2, 0.8)
#' d1 <- dnorm(x, 1, log = TRUE); d2 <- dexp(y, 2, log = TRUE); d3 <- dbeta(z, 2, 1, log = TRUE)
#' p1 <- pnorm(x, 1); p2 <- pexp(y, 2); p3 <- pbeta(z, 2, 1)
#' R <- matrix(c(1,0.5,0.3,0.5,1,0.4,0.3,0.4,1), nrow = 3)
#'
#' ## Based on correlation matrix
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cmvgauss(R), log = TRUE)
#'
#' ## Based on precision matrix
#' Q <- solve(R)
#' dmvcopula(cbind(d1, d2, d3), cbind(p1, p2, p3), copula = cgmrf(Q), log = TRUE)
#'
#' ## Parameterisation inside a model
#' # using RTMB::unstructured to get a valid correlation matrix
#' library(RTMB)
#' d <- 5 # dimension
#' cor_func <- unstructured(d)
#' npar <- length(cor_func$parms())
#' R <- cor_func$corr(rep(0.1, npar))
cgmrf <- function(Q) {
  if (nrow(Q) != ncol(Q)) {
    stop("Q must be a square matrix.")
  }

  d <- nrow(Q)

  ## if Q is dense, convert to sparse matrix
  if (is.matrix(Q)) {
    # Precompute log|Q| and the "precision minus I" matrix
    L <- chol(Q)
    logdetQ <- 2 * sum(log(diag(L)))
    Prec <- Q - AD(diag(d))

    f <- function(U) {
      if (ncol(U) != d) stop("Input dimension does not match correlation matrix R.")
      Z <- apply(U, 2, qnorm)

      quadform <- rowSums((Z %*% Prec) * Z)
      0.5 * logdetQ - 0.5 * quadform
    }

  } else {
    f <- function(U) {
      if (ncol(U) != d) stop("Input dimension does not match correlation matrix R.")
      Z <- apply(U, 2, qnorm)

      logdens <- RTMB::dgmrf(Z, Q = Q, log = TRUE)
      logdens <- logdens + 0.5 * rowSums(Z^2)
      logdens <- logdens + 0.5 * d * log(2 * pi)
      return(logdens)
    }
  }
  return(f)
}
