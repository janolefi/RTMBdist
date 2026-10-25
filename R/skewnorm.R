# Vectorising sn package functions
psn <- Vectorize(sn::psn)
qsn <- Vectorize(sn::qsn)

# n standard skew normal draws, vectorised over alpha, via the representation
# delta * |U0| + sqrt(1 - delta^2) * U1 with U0, U1 iid N(0, 1)
rsn_std <- function(n, alpha) {
  alpha <- rep_len(alpha, n)
  delta <- alpha / sqrt(1 + alpha^2)
  delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n)
}

#' Skew normal distribution
#'
#' Density, distribution function, quantile function, and random generation for
#' the skew normal distribution.
#'
#' @details
#' This implementation of \code{dskewnorm} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#' When called during AD taping, \code{pskewnorm} instead integrates the density numerically with the AD-compatible \code{\link[RTMB:ADintegrate]{integrate}} of \code{RTMB}, which makes it AD-compatible as well.
#' See \code{sn::\link[sn]{dsn}} for more details.
#'
#' \deqn{f(x;\,\xi,\omega,\alpha) = \frac{2}{\omega}\,\phi\!\left(\frac{x-\xi}{\omega}\right)\Phi\!\left(\alpha\frac{x-\xi}{\omega}\right),}
#' where \eqn{\phi} and \eqn{\Phi} are the standard normal PDF and CDF.
#'
#' @seealso [skewnorm2], [skewt], [skewt2]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return
#' @param xi location parameter
#' @param omega scale parameter, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
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
    return(dGenericOSA("dskewnorm", x=x, xi=xi, omega=omega, alpha=alpha, log=log))
  }

  z = (x - xi) / omega # standardised observation

  log_normal_density <- RTMB::dnorm(z, log = TRUE)
  log_skew_component <- log(2) - log(omega) + RTMB::pnorm(alpha * z, log.p = TRUE)

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
pskewnorm <- function(q, xi = 0, omega = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE) {

  # centre for numerical integration is the mean, written to allow alpha = +/- Inf
  a <- value_of(alpha)
  centre <- value_of(xi) + value_of(omega) * sqrt(2 / pi) * sign(a) / sqrt(1 + 1 / a^2)
  pars <- list(xi = xi, omega = omega, alpha = alpha)

  if(!ad_context()) {
    # ensure omega > 0
    if (any(omega <= 0)) stop("omega must be strictly positive.")
  } else {
    # sn::psn is not AD-compatible, hence integrate the density numerically
    return(numerical_cdf(dskewnorm, q, pars, centre = centre, scale = omega,
                         lower.tail = lower.tail, log.p = log.p))
  }

  p <- psn(x = q, xi = xi, omega = omega, alpha = alpha)
  # sn is exact for alpha = +/- Inf, where the density is zero on one side
  refine_tails(p, dskewnorm, q, pars, centre = centre, scale = omega,
               lower.tail = lower.tail, log.p = log.p, use = is.finite(alpha))
}

#' @rdname skewnorm
#' @export
#' @importFrom sn qsn
qskewnorm <- function(p, xi = 0, omega = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE) {

  if(!ad_context()) {
    # ensure omega > 0
    if (any(omega <= 0)) stop("omega must be strictly positive.")
  }

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qsn(p = p, xi = xi, omega = omega, alpha = alpha)
}

#' @rdname skewnorm
#' @export
rskewnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  if (any(omega <= 0)) stop("omega must be strictly positive.")

  if (length(n) > 1L) n <- length(n)
  rep_len(xi, n) + rep_len(omega, n) * rsn_std(n, alpha)
}
