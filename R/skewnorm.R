# Vecotorizing sn package functions
psn <- Vectorize(sn::psn)
qsn <- Vectorize(sn::qsn)
rsn <- Vectorize(sn::rsn)

#' Skew normal distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the skew normal distribution.
#'
#' @details
#' This implementation of \code{dskewnorm} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#' See \code{sn::\link[sn]{dsn}} for more details.
#'
#' @seealso [skewnorm2], [skewt], [skewt2]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param xi location parameter
#' @param omega scale parameter, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param ... additional parameters to be passed to the \code{sn} package functions for \code{pskewnorm} and \code{qskewnorm}.
#'
#' @return
#' \code{dskewnorm} gives the density, \code{pskewnorm} gives the distribution function, \code{qskewnorm} gives the quantile function, and \code{rskewnorm} generates random deviates.
#'
#' @examples
#' # alpha is skew parameter
#' x <- rskewnorm(1, alpha = 1)
#' d <- dskewnorm(x, alpha = 1)
#' p <- pskewnorm(x, alpha = 1)
#' q <- qskewnorm(p, alpha = 1)
#' @name skewnorm
NULL

#' @rdname skewnorm
#' @export
#' @importFrom RTMB dnorm
#' @importFrom RTMB pnorm
dskewnorm <- function(x, xi = 0, omega = 1, alpha = 0, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure omega > 0
    if (any(omega <= 0)) stop("omega must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dskewnorm", x=x, xi=xi, omega=omega, alpha=alpha, log=log))
  }
  if(inherits(x, "osa")) {
    # return(dGenericOSA("dskewnorm", x=x, xi=xi, omega=omega, alpha=alpha, log=log))
    stop("Currently, skew normal does not support OSA residuals.")
  }

  z = (x - xi) / omega # standardised observation

  log_normal_density <- RTMB::dnorm(z, log = TRUE)
  log_skew_component <- log(2) - log(omega) + log(1e-300 + RTMB::pnorm(alpha * z))

  log_density = log_normal_density + log_skew_component

  if(log) {
    return(log_density)
  } else{
    return(exp(log_density))
  }
}

#' @rdname skewnorm
#' @export
#' @importFrom sn psn
pskewnorm <- function(q, xi = 0, omega = 1, alpha = 0, ...) {

  if(!ad_context()) {
    # ensure omega > 0
    if (any(omega <= 0)) stop("omega must be strictly positive.")
  }

  psn(x = q, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn qsn
qskewnorm <- function(p, xi = 0, omega = 1, alpha = 0, ...) {

  if(!ad_context()) {
    # ensure omega > 0
    if (any(omega <= 0)) stop("omega must be strictly positive.")
  }

  qsn(p = p, xi = xi, omega = omega, alpha = alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom sn rsn
rskewnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  if (any(omega <= 0)) stop("omega must be strictly positive.")

  # elementwise generation
  if (length(xi) == 1L && length(omega) == 1L && length(alpha) == 1L) {
    return(sn::rsn(n, xi, omega, alpha))
  }

  len <- max(length(xi), length(omega), length(alpha), n)
  xi <- rep_len(xi, len)
  omega <- rep_len(omega, len)
  alpha <- rep_len(alpha, len)

  vapply(seq_len(len), function(i)
    sn::rsn(1, xi[i], omega[i], alpha[i]), numeric(1))
}
