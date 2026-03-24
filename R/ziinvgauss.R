#' Zero-inflated inverse Gaussian distribution
#'
#' Density, distribution function, and random generation for
#' the zero-inflated inverse Gaussian distribution.
#'
#' @details
#' This implementation of \code{zidinvgauss} allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param n number of random values to return
#' @param mean location parameter
#' @param shape shape parameter, must be positive.
#' @param zeroprob zero-probability, must be in \eqn{[0, 1]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dziinvgauss} gives the density, \code{pziinvgauss} gives the distribution function, and \code{rziinvgauss} generates random deviates.
#'
#' @examples
#' x <- rziinvgauss(1, 1, 2, 0.5)
#' d <- dziinvgauss(x, 1, 2, 0.5)
#' p <- pziinvgauss(x, 1, 2, 0.5)
#' @name ziinvgauss
NULL

#' @rdname ziinvgauss
#' @export
dziinvgauss <- function(x, mean = 1, shape = 1, zeroprob = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dziinvgauss", x=x, mean=mean, shape=shape, zeroprob=zeroprob, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dziinvgauss", x=x, mean=mean, shape=shape, zeroprob=zeroprob, log=log))
  }

  # logdens <- numeric(length(x))
  # zero_idx <- (x == 0)
  #
  # # Zero inflation part
  # logdens[zero_idx] <- log(zeroprob)
  # logdens[!zero_idx] <- log(1 - zeroprob) +
  #   dinvgauss(x[!zero_idx], mean = mean, shape = shape, log = TRUE)

  logdens <- dinvgauss(x, mean = mean, shape = shape, log = TRUE)
  logdens <- log_zi(x, logdens, zeroprob)

  if (log) return(logdens)
  return(exp(logdens))
}
#' @rdname ziinvgauss
#' @export
#' @usage pziinvgauss(q, mean = 1, shape = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE)
pziinvgauss <- function(q, mean = 1, shape = 1, zeroprob = 0, lower.tail = TRUE, log.p = FALSE) {

  # p <- numeric(length(q))
  #
  # below_zero <- q < 0
  # is_zero <- q == 0
  # positive <- q > 0
  #
  # p[below_zero] <- 0
  # p[is_zero] <- zeroprob
  # p[positive] <- zeroprob + (1 - zeroprob) * pinvgauss(q[positive], mean=mean, shape=shape)

  # s1 <- 2 * sign(q) - 1 # gives -3 for q < 0, -1 for q == 0, and 1 for q > 0
  # s2 <- sign(3 + s1) # only zero or 1

  # p <- 0.5 * (1 - s1) * s2 * zeroprob +
  #   0.5 * (1 + s1) * s2 * (zeroprob + (1 - zeroprob) * pinvgauss(q, mean, shape))

  p <- iszero(q) * zeroprob +
    ispos_strict(q) * (zeroprob + (1 - zeroprob) * pinvgauss(q, mean, shape))

  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p))
  return(p)
}
#' @rdname ziinvgauss
#' @export
#' @importFrom statmod rinvgauss
rziinvgauss <- function(n, mean = 1, shape = 1, zeroprob = 0) {
  # ensure mean, shape >= 0 zeroprob in [0,1]
  if (any(mean <= 0)) stop("mean must be > 0")
  if (any(shape <= 0)) stop("shape must be > 0")
  if (any(zeroprob < 0 | zeroprob > 1)) stop("zeroprob must be in [0,1]")

  u <- runif(n)
  res <- rep(0, n)
  is_zero <- u < zeroprob
  res[!is_zero] <- rinvgauss(sum(!is_zero), mean = mean, shape = shape)

  u <- runif(n)
  res <- ifelse(u < zeroprob, 0, rinvgauss(n, mean = mean, shape = shape))
  return(res)
}
