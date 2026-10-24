#' Generalised beta distribution of the second kind and its special cases
#'
#' Density, distribution function, quantile function, and random generation for
#' the generalised beta distribution of the second kind (GB2) and three of its
#' special cases: the Dagum, Singh-Maddala (Burr XII) and Lomax (Pareto II)
#' distributions.
#'
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' The GB2 distribution with scale \eqn{b} and shape parameters \eqn{a}, \eqn{p}
#' and \eqn{q} has density
#' \deqn{f(x) = \frac{a\, x^{ap - 1}}{b^{ap} B(p, q) \left(1 + (x/b)^a\right)^{p + q}}, \quad x > 0,}
#' and its distribution function is a regularised incomplete beta function,
#' \deqn{F(x) = I_z(p, q), \qquad z = \frac{(x/b)^a}{1 + (x/b)^a}.}
#' The special cases fix some of the shapes, and their arguments keep the names
#' of the shapes that remain:
#' \itemize{
#'   \item Dagum: \eqn{q = 1}, \code{ddagum(x, scale, shape1, shape2)}.
#'   \item Singh-Maddala, or Burr type XII: \eqn{p = 1}, \code{dsinmad(x, scale, shape1, shape3)}.
#'   \item Lomax, or Pareto type II: \eqn{a = p = 1}, \code{dlomax(x, scale, shape3)}.
#' }
#' The log-logistic distribution \code{\link{dllogis}} is the GB2 with
#' \eqn{p = q = 1}, and the beta prime distribution \code{\link{dbetaprime}} the
#' GB2 with \eqn{a = b = 1}. The \eqn{k}-th moment exists for \eqn{-ap < k < aq}.
#'
#' The distribution functions are differentiable, so one-step-ahead (OSA)
#' residuals via \code{RTMB::\link[RTMB]{oneStepPredict}} with
#' \code{method = "cdf"} are supported.
#'
#' @references
#' McDonald, J. B. (1984) Some generalized functions for the size distribution of income.
#' Econometrica, 52, 647-663.
#'
#' Kleiber, C. and Kotz, S. (2003) Statistical Size Distributions in Economics and Actuarial Sciences. Wiley.
#'
#' @seealso [llogis], [betaprime], [pareto]
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param scale scale parameter \eqn{b}, must be positive.
#' @param shape1 shape parameter \eqn{a}, must be positive.
#' @param shape2 shape parameter \eqn{p}, must be positive.
#' @param shape3 shape parameter \eqn{q}, must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P[X \le x]}, otherwise \eqn{P[X > x]}.
#'
#' @return
#' \code{d...} gives the density, \code{p...} gives the distribution function,
#' \code{q...} gives the quantile function, and \code{r...} generates random deviates.
#'
#' @examples
#' set.seed(1)
#' x <- rgb2(5, scale = 2, shape1 = 3, shape2 = 1.5, shape3 = 0.8)
#' d <- dgb2(x, 2, 3, 1.5, 0.8)
#' p <- pgb2(x, 2, 3, 1.5, 0.8)
#' q <- qgb2(p, 2, 3, 1.5, 0.8)
#'
#' d <- ddagum(x, 2, 3, 1.5)
#' d <- dsinmad(x, 2, 3, 0.8)
#' d <- dlomax(x, 2, 0.8)
#' @name gb2
NULL

#' @rdname gb2
#' @export
dgb2 <- function(x, scale = 1, shape1, shape2, shape3, log = FALSE) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    check_gb2(scale, shape1, shape2, shape3)
  }

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dgb2", x = x, scale = scale, shape1 = shape1, shape2 = shape2,
                       shape3 = shape3, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dgb2", x = x, scale = scale, shape1 = shape1, shape2 = shape2,
                       shape3 = shape3, log = log))
  }

  # x is kept away from zero, where log(x) would be -Inf, and the density is set to
  # zero below the support by the indicator
  xc <- pmax.ad(x, .Machine$double.xmin)
  t <- shape1 * log(xc / scale) # log((x / b)^a)

  logdens <- log(shape1) - log(xc) + shape2 * t - lbeta(shape2, shape3) -
    (shape2 + shape3) * logspace_add(0, t) +
    log(1 - smaller(x, 0))

  if (log) return(logdens)
  return(exp(logdens))
}

#' @rdname gb2
#' @export
#' @importFrom RTMB plogis
pgb2 <- function(q, scale = 1, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) check_gb2(scale, shape1, shape2, shape3)

  t <- shape1 * log(q * (1 - smaller(q, 0)) / scale) # negative q is treated as 0

  # the upper tail as I_{1 - z}(q, p) rather than 1 - p, which keeps it accurate
  if (lower.tail) {
    p <- beta_ad(pbeta, plogis(t), shape2, shape3)
  } else {
    p <- beta_ad(pbeta, plogis(-t), shape3, shape2)
  }

  if (log.p) p <- log(p)
  p
}

#' @rdname gb2
#' @export
qgb2 <- function(p, scale = 1, shape1, shape2, shape3, lower.tail = TRUE, log.p = FALSE) {

  if (!ad_context()) check_gb2(scale, shape1, shape2, shape3)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  if (!ad_context()) {
    if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  }

  scale * exp(qlogis(beta_ad(qbeta, p, shape2, shape3)) / shape1)
}

#' @rdname gb2
#' @export
#' @importFrom stats rbeta
rgb2 <- function(n, scale = 1, shape1, shape2, shape3) {

  check_gb2(scale, shape1, shape2, shape3)

  scale * exp(qlogis(rbeta(n, shape2, shape3)) / shape1)
}

check_gb2 <- function(scale, shape1, shape2, shape3) {
  if (any(scale <= 0)) stop("scale must be positive.")
  if (any(shape1 <= 0)) stop("shape1 must be positive.")
  if (any(shape2 <= 0)) stop("shape2 must be positive.")
  if (any(shape3 <= 0)) stop("shape3 must be positive.")
}


#' @rdname gb2
#' @export
ddagum <- function(x, scale = 1, shape1, shape2, log = FALSE) {

  if (!ad_context()) simulation_check(as.list(environment()))

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("ddagum", x = x, scale = scale, shape1 = shape1, shape2 = shape2, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("ddagum", x = x, scale = scale, shape1 = shape1, shape2 = shape2, log = log))
  }

  dgb2(x, scale, shape1, shape2, 1, log = log)
}

#' @rdname gb2
#' @export
pdagum <- function(q, scale = 1, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  pgb2(q, scale, shape1, shape2, 1, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
qdagum <- function(p, scale = 1, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
  qgb2(p, scale, shape1, shape2, 1, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
rdagum <- function(n, scale = 1, shape1, shape2) {
  rgb2(n, scale, shape1, shape2, 1)
}


#' @rdname gb2
#' @export
dsinmad <- function(x, scale = 1, shape1, shape3, log = FALSE) {

  if (!ad_context()) simulation_check(as.list(environment()))

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dsinmad", x = x, scale = scale, shape1 = shape1, shape3 = shape3, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dsinmad", x = x, scale = scale, shape1 = shape1, shape3 = shape3, log = log))
  }

  dgb2(x, scale, shape1, 1, shape3, log = log)
}

#' @rdname gb2
#' @export
psinmad <- function(q, scale = 1, shape1, shape3, lower.tail = TRUE, log.p = FALSE) {
  pgb2(q, scale, shape1, 1, shape3, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
qsinmad <- function(p, scale = 1, shape1, shape3, lower.tail = TRUE, log.p = FALSE) {
  qgb2(p, scale, shape1, 1, shape3, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
rsinmad <- function(n, scale = 1, shape1, shape3) {
  rgb2(n, scale, shape1, 1, shape3)
}


#' @rdname gb2
#' @export
dlomax <- function(x, scale = 1, shape3, log = FALSE) {

  if (!ad_context()) simulation_check(as.list(environment()))

  # potentially escape to RNG or CDF
  if (inherits(x, "simref")) {
    return(dGenericSim("dlomax", x = x, scale = scale, shape3 = shape3, log = log))
  }
  if (inherits(x, "osa")) {
    return(dGenericOSA("dlomax", x = x, scale = scale, shape3 = shape3, log = log))
  }

  dgb2(x, scale, 1, 1, shape3, log = log)
}

#' @rdname gb2
#' @export
plomax <- function(q, scale = 1, shape3, lower.tail = TRUE, log.p = FALSE) {
  pgb2(q, scale, 1, 1, shape3, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
qlomax <- function(p, scale = 1, shape3, lower.tail = TRUE, log.p = FALSE) {
  qgb2(p, scale, 1, 1, shape3, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gb2
#' @export
rlomax <- function(n, scale = 1, shape3) {
  rgb2(n, scale, 1, 1, shape3)
}
