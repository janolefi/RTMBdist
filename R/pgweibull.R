## This implementation was written by Julia Dyck
## and is taken from https://github.com/julia-dyck/BWSPsignal/blob/main/R/pgw_functions.R
## it is slightly modified to fit package style and RTMB requirements.

#' Power generalized Weibull distribution
#'
#' @description Survival, hazard, cumulative distribution,
#' density, quantile and sampling function for the power generalized
#' Weibull (PgW) distribution with parameters \code{scale}, \code{shape} and \code{powershape}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param scale positive scale parameter
#' @param shape positive shape parameter
#' @param powershape positive power shape parameter
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @details The survival function of the PgW distribution is:
#' \deqn{
#'     S(x) = \exp \left\{ 1 - \left[ 1 + \left(\frac{x}{\theta}\right)^{\nu}\right]^{\frac{1}{\gamma}} \right\}.
#' }
#' The hazard function is
#' \deqn{
#'\frac{\nu}{\gamma\theta^{\nu}}\cdot x^{\nu-1}\cdot  \left[ 1 + \left(\frac{x}{\theta}\right)^{\nu}\right]^{\frac{1}{\gamma-1}}
#' }
#' The cumulative distribution function is then \eqn{F(x) = 1 - S(x)} and the density function
#' is \eqn{S(x)\cdot h(x)}.
#'
#' If both shape parameters equal 1, the PgW distribution reduces to the exponential distribution
#' (see \code{\link[stats]{dexp}}) with \eqn{\texttt{rate} = 1/\texttt{scale}}
#' If the power shape parameter equals 1, the PgW distribution simplifies to the Weibull distribution
#' (see \code{\link[stats]{dweibull}}) with the same parametrization.
#'
#' @return
#' \code{dpgweibull} gives the density, \code{ppgweibull} gives the distribution function, \code{qpgweibull} gives the quantile function, and \code{rpgweibull} generates random deviates.
#' \code{spgweibull} gives the survival function and \code{hpgweibull} gives the hazard function.
#'
#' @examples
#' x <- rpgweibull(1, 2, 2, 3)
#' d <- dpgweibull(x, 2, 2, 3)
#' p <- ppgweibull(x, 2, 2, 3)
#' q <- qpgweibull(p, 2, 2, 3)
#' s <- spgweibull(x, 2, 2, 3)
#' h <- hpgweibull(x, 2, 2, 3)
#' @name pgweibull
NULL

#' @rdname pgweibull
#' @export
spgweibull <- function(q, scale = 1, shape = 1, powershape = 1, log = FALSE) {

  if(!ad_context()) {
    # ensure scale, shape, powershape > 0
    if (any(scale <= 0)) stop("scale must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(powershape <= 0)) stop("powershape must be strictly positive.")
  }

  # renaming to match the formula
  theta <- scale
  nu <- shape
  gamma <- powershape

  log1pxtn <- log1p((q / theta)^nu)
  log_inner <- log1pxtn / gamma     # stabilised exponent
  Slog <- -expm1(log_inner)         # log survival

  if(log) return(Slog)
  return(exp(Slog))
}

#' @rdname pgweibull
#' @export
hpgweibull <- function(x, scale = 1, shape = 1, powershape = 1, log = FALSE) {

  if(!ad_context()) {
    # ensure scale, shape, powershape > 0
    if (any(scale <= 0)) stop("scale must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(powershape <= 0)) stop("powershape must be strictly positive.")
  }

  # renaming to match the formula
  theta <- scale
  nu <- shape
  gamma <- powershape

  log_hazard_values <- log(nu) - log(gamma) -nu * log(theta) +
    (nu-1) * log(x) + (1/gamma - 1) * log1p((x/theta)^nu)
  if(log == TRUE){
    return(log_hazard_values)
  }
  else{
    hazard_values = exp(log_hazard_values)
    return(hazard_values)
  }
}

#' @rdname pgweibull
#' @export
#' @usage
#' ppgweibull(q, scale = 1, shape = 1, powershape = 1,
#'            lower.tail = TRUE, log.p = FALSE)
ppgweibull <- function(q, scale = 1, shape = 1, powershape = 1,
                       lower.tail = TRUE, log.p = FALSE) {

  p <- 1 - spgweibull(q, scale=scale, shape=shape, powershape=powershape)

  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)

  return(p)
}

#' @rdname pgweibull
#' @export
dpgweibull <- function(x, scale = 1, shape = 1, powershape = 1, log = FALSE) {

  if(!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure scale, shape, powershape > 0
    if (any(scale <= 0)) stop("scale must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(powershape <= 0)) stop("powershape must be strictly positive.")
  }

  # potentially escape to RNG or CDF
  if(inherits(x, "simref")) {
    return(dGenericSim("dpgweibull", x=x, scale=scale, shape=shape, powershape=powershape, log=log))
  }
  if(inherits(x, "osa")) {
    return(dGenericOSA("dpgweibull", x=x, scale=scale, shape=shape, powershape=powershape, log=log))
  }

  # renaming to match the formula
  theta <- scale
  nu <- shape
  gamma <- powershape

  log1pxtn <- log1p((x/theta)^nu)
  log_inner <- log1pxtn / gamma # stabilising double power
  log_S <- -expm1(log_inner)
  log_h <- log(nu) - log(gamma) - nu * log(theta) +
    (nu-1) * log(x) + (1/gamma - 1) * log1pxtn

  logdens <- log_S + log_h

  if(log) return(logdens)
  return(exp(logdens))
}

#' @rdname pgweibull
#' @export
qpgweibull <- function(p, scale = 1, shape = 1, powershape = 1) {

  if(!ad_context()) {
    # ensure scale, shape, powershape > 0
    if (any(scale <= 0)) stop("scale must be strictly positive.")
    if (any(shape <= 0)) stop("shape must be strictly positive.")
    if (any(powershape <= 0)) stop("powershape must be strictly positive.")
  }

  # renaming to match the formula
  theta <- scale
  nu <- shape
  gamma <- powershape

  # inverse cumulative distribution fct.
  theta * ( ( 1 - log1p(-p) )^gamma - 1 )^(1/nu)
}

#' @rdname pgweibull
#' @export
#' @importFrom stats runif
rpgweibull <- function(n, scale = 1, shape = 1, powershape = 1){

  # ensure scale, shape, powershape > 0
  if (any(scale <= 0)) stop("scale must be strictly positive.")
  if (any(shape <= 0)) stop("shape must be strictly positive.")
  if (any(powershape <= 0)) stop("powershape must be strictly positive.")

  u = stats::runif(n) # generate random sample from uniform distribution
  qpgweibull(u, scale = scale, shape = shape, powershape = powershape)
}
