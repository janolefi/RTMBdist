#' Multivariate t distribution
#'
#' Density and and random generation for the multivariate t distribution
#'
#' @details
#' This implementation of \code{dmvt} allows for automatic differentiation with \code{RTMB}.
#'
#' Note: for \code{df} \eqn{\le 1} the mean is undefined, and for \code{df} \eqn{\le 2} the covariance is infinite.
#' For \code{df} > 2, the covariance is \code{df/(df-2) * Sigma}.
#'
#' @param x vector or matrix of quantiles
#' @param n number of random values to return.
#' @param mu vector or matrix of location parameters (mean if \code{df} > 1)
#' @param Sigma positive definite scale matrix (proportional to the covariance matrix if \code{df} > 2)
#' @param df degrees of freedom; must be positive
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{dmvt} gives the density, \code{rmvt} generates random deviates.
#'
#' @examples
#' # single mu
#' mu <- c(1,2,3)
#' Sigma <- diag(c(1,1,1))
#' df <- 5
#' x <- rmvt(2, mu, Sigma, df)
#' d <- dmvt(x, mu, Sigma, df)
#' # vectorised over mu
#' mu <- rbind(c(1,2,3), c(0, 0.5, 1))
#' x <- rmvt(2, mu, Sigma, df)
#' d <- dmvt(x, mu, Sigma, df)
#' @name mvt
NULL

#' @rdname mvt
#' @export
#' @import RTMB
dmvt <- function(x, mu, Sigma, df, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order

    if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-10))) {
      stop("Sigma must be symmetric.")
    }

    if(any(df <= 0)) {
      stop("Degrees of freedom 'df' must be positive.")
    }
  }

  # potentially escape to RNG or produce error for CDF
  if(inherits(x, "simref")) {
    n <- if (is.matrix(x)) nrow(x) else 1
    x[] <- rmvt(n, mu=mu, Sigma=Sigma, df=df)
    return(0)
  }
  if(inherits(x, "osa")) {
    stop("Multivariate t does not support OSA residuals.")
  }

  # Input reshaping
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)

  d <- ncol(x) # dimension

  if (is.null(dim(mu))){
    if(length(mu) == 1) mu <- matrix(mu, nrow = 1, ncol = d)
    else mu <- matrix(mu, nrow = 1)
  }

  # Dimension checks
  if (ncol(x) != ncol(mu)) {
    stop("x and mu must have the same number of columns (dimension).")
  }
  if (nrow(x) > 1) {
    if (nrow(mu) != 1 && nrow(mu) != nrow(x)) {
      stop("If x is a matrix, mu must have either one row or the same number of rows as x.")
    }
    if (nrow(mu) == 1) {
      mu <- matrix(rep(mu, each = nrow(x)), nrow = nrow(x))
    }
  }

  # Sigma checks
  if (length(dim(Sigma)) != 2 || nrow(Sigma) != d || ncol(Sigma) != d) {
    stop("Sigma must be a square matrix with dimension equal to ncol(x).")
  }

  # Ensure df has correct length
  if (length(df) == 1){
    df <- rep(df, nrow(x))
  } else if(length(df) != nrow(x)) {
    stop("df must be either a scalar or have the same length as the number of rows in x.")
  }


  # Log determinant of Sigma
  logdet <- determinant(Sigma, logarithm=TRUE)$modulus

  # Mahalanobis distances
  z <- x - mu # centered
  quadform <- rowSums((z %*% solve(Sigma)) * z)

  # log density constant
  const <- lgamma((df + d) / 2) - lgamma(df / 2) -
    0.5 * (d * log(df * pi) + logdet)

  logdens <- const - 0.5 * (df + d) * log1p(quadform / df)

  if(log) return(logdens)
  return(exp(logdens))
}
#' @rdname mvt
#' @export
#' @importFrom stats rchisq rnorm
rmvt <- function(n, mu, Sigma, df) {

  # Input reshaping
  if (is.null(dim(mu))) mu <- matrix(mu, nrow = 1)
  d <- ncol(mu) # dimension
  if (nrow(mu) == 1) {
    mu <- matrix(rep(mu, each = n), nrow = n)
  } else if (nrow(mu) != n) {
    stop("If mu has multiple rows, it must have the same number of rows as n.")
  }

  # Sigma checks
  if (length(dim(Sigma)) != 2 || nrow(Sigma) != d || ncol(Sigma) != d) {
    stop("Sigma must be a square matrix with dimension equal to ncol(x).")
  }

  if (!isTRUE(all.equal(Sigma, t(Sigma), tolerance = 1e-10))) {
    stop("Sigma must be symmetric.")
  }
  if(any(df <= 0)) {
    stop("Degrees of freedom 'df' must be positive.")
  }

  # Ensure df has correct length
  if (length(df) == 1){
    df <- rep(df, n)
  } else if(length(df) != n) {
    stop("df must be either a scalar or have the same length as the number of rows in x.")
  }

  # Simulation procedure
  u <- rchisq(n, df=df)

  z <- matrix(rnorm(n * d), nrow = n, ncol = d)
  L <- chol(Sigma)
  y <- z %*% L # y ~ N(0, Sigma)
  x <- y / sqrt(u / df) + mu

  return(x)
}
