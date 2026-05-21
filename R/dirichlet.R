#' Dirichlet distribution
#'
#' Density and and random generation for the Dirichlet distribution.
#'
#' @details
#' This implementation of \code{ddirichlet} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x vector or matrix of quantiles. If \code{x} is a vector, it needs to sum to one.
#' If \code{x} is a matrix, each row should sum to one.
#' @param n number of random values to return.
#' @param alpha vector or matrix of positive shape parameters
#' @param log logical; if \code{TRUE}, densities \eqn{p} are returned as \eqn{\log(p)}.
#'
#' @return
#' \code{ddirichlet} gives the density, \code{rdirichlet} generates random deviates.
#'
#' @examples
#' # single alpha
#' alpha <- c(1,2,3)
#' x <- rdirichlet(1, alpha)
#' d <- ddirichlet(x, alpha)
#' # vectorised over alpha
#' alpha <- rbind(alpha, 2*alpha)
#' x <- rdirichlet(2, alpha)
#  d <- ddirichlet(x, alpha)
#' @name dirichlet
NULL

#' @rdname dirichlet
#' @export
#' @import RTMB
ddirichlet <- function(x, alpha, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
  }

  # potentially escape to RNG or produce error for CDF
  if(inherits(x, "simref")) {
    n <- if (is.matrix(x)) nrow(x) else 1
    x[] <- rdirichlet(n, alpha=alpha)
    return(0)
  }
  if(inherits(x, "osa")) {
    return(ddirichlet_osa(x=x, alpha=alpha, log=log))
  }

  # Check if x and alpha are vectors by checking if they have dimensions
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  if (is.null(dim(alpha))) alpha <- matrix(alpha, nrow = 1)

  # Ensure x and alpha have the same number of columns
  if (ncol(x) != ncol(alpha)) {
    stop("x and alpha must have the same number of columns (categories).")
  }
  # If x is a matrix, check that alpha either has the same number of rows or one
  if(nrow(x) > 1) {
    if(nrow(alpha) != 1 && nrow(alpha) != nrow(x)) {
      stop("If x is a matrix, alpha must have either one row or the same number of rows as x.")
    }
  }
  # If x has more than one row and alpha has only one row, replicate it
  if(nrow(x) > 1 && nrow(alpha) == 1){
    alpha <- matrix(rep(alpha, each = nrow(x)), nrow = nrow(x))
  }

  # Compute log of the multivariate beta function B(alpha) for each row
  log_B_alpha <- rowSums(lgamma(alpha)) - lgamma(rowSums(alpha))

  # Compute log density for each row
  log_density <- rowSums((alpha - 1) * log(x)) - log_B_alpha

  log_density <- as.numeric(log_density)

  # Return result based on the 'log' argument
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
#' @rdname dirichlet
#' @export
#' @importFrom stats rgamma
rdirichlet <- function(n, alpha) {

  if(is.vector(alpha)){
    longalpha <- rep(alpha, n)
    k <- length(alpha)
  } else if(is.matrix(alpha)){
    if(nrow(alpha) != n) stop("Number of rows in alpha must match n.")

    longalpha <- as.vector(t(alpha))
    k <- ncol(alpha)
  }

  x <- rgamma(n * k, shape = longalpha, scale = 1)
  x <- matrix(x, nrow = n, ncol = k, byrow = TRUE)
  x / rowSums(x)
}

ddirichlet_osa <- function(x, alpha, log = FALSE) {
  ## log only
  stopifnot(isTRUE(log))
  ## Matrix case
  if (is.matrix(x)) {
    if (is.matrix(alpha))
      stopifnot(identical(dim(x), dim(alpha)))
    else
      alpha <- matrix(alpha, nrow(x), ncol(x), byrow=TRUE)
    ans <- AD(numeric(nrow(x)))
    for (i in seq_len(nrow(x))) {
      ans[[i]] <- ddirichlet_osa(x[i,], alpha[i,], log=log)
    }
    return(ans)
  }
  ## Vector case
  alpha <- rep(alpha, length.out=length(x))
  ## Permute
  perm <- order(attr(x@keep, "ord")) ## FIXME: Make extractor in osa.R ?
  x <- x[perm]
  alpha <- alpha[perm]
  ## Factorize in successive betas
  sx <- sum(x@x)
  sa <- sum(alpha)
  ## retun value
  ans <- 0
  if (length(x) >= 2) {
    ## Draw first
    sa <- sa - alpha[1]
    ans <- ans + dbeta(x[1], alpha[1], sa, log=TRUE)
    ## Draw the rest, but not the last
    for (i in seq_along(x)[-c(1, length(x))]) {
      sx <- sx - x@x[i-1]
      sa <- sa - alpha[i]
      ## x[i] ~ Scaled Beta, but 'dbeta' doesn't have a scale argument
      xi <- x[i]
      xi@x <- xi@x / sx
      ans <- ans + dbeta(xi, alpha[i], sa, log = TRUE)
      ans <- ans - xi@keep[,1] * log(sx)
    }
  }
  ## Draw last: Always a one-point measure
  if (length(x) >= 1) {
    xi <- x[length(x)]
    xi@x <- 0
    ans <- ans - dbinom(xi, 1, 0, log=TRUE)
  }
  ans
}
