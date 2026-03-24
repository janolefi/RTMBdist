#' Moment-parameterised skew t distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the skew t distribution reparameterised so that \code{mean} and \code{sd}
#' correspond to the *true* mean and standard deviation.
#'
#' @details
#'
#' This corresponds to the skew t type 2 distribution in GAMLSS (\code{\link[gamlss.dist]{ST2}}), see pp. 411-412 of Rigby et al. (2019) and the version implemented in the \code{sn} package.
#' However, it is reparameterised in terms of a standard deviation parameter \code{sd} rather than just a scale parameter \code{sigma}. Details of this reparameterisation are given below.
#' This implementation of \code{dskewt} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#' See \code{sn::\link[sn]{dst}} for more details.
#'
#' \strong{Caution:} In a numerial optimisation, the \code{skew} parameter should NEVER be initialised with exactly zero.
#' This will cause the initial and all subsequent derivatives to be exactly zero and hence the parameter will remain at its initial value.
#'
#' For given \code{skew} \eqn{= \alpha} and \code{df} = \eqn{\nu}, define
#' \deqn{
#' \delta = \alpha / \sqrt{1 + \alpha^2}, \qquad
#' b_\nu = \sqrt{\nu / \pi}\, \Gamma((\nu-1)/2)/\Gamma(\nu/2),
#' }
#' then
#' \deqn{
#' E(X) = \mu + \sigma \delta b_\nu,\quad
#' Var(X) = \sigma^2 \left( \frac{\nu}{\nu-2} - \delta^2 b_\nu^2 \right).
#' }
#'
#' @seealso [skewt], [skewnorm], [skewnorm2]
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of random values to return.
#' @param mean mean parameter
#' @param sd standard deviation parameter, must be positive.
#' @param skew skewness parameter, can be positive or negative.
#' @param df degrees of freedom, must be positive.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param tol a scalar value which regulates the accuracy of the result of qsn, measured on the probability scale.
#' @param method an integer value between 0 and 5 which selects the computing method; see ‘Details’ in the \code{\link[sn]{pst}} documentation below for the meaning of these values. If method=0 (default value), an automatic choice is made among the four actual computing methods, depending on the other arguments.
#'
#' @return
#' \code{dskewt2} gives the density, \code{pskewt2} gives the distribution function, \code{qskewt2} gives the quantile function, and \code{rskewt2} generates random deviates.
#'
#' @examples
#' x <- rskewt2(1, 1, 2, 5, 5)
#' d <- dskewt2(x, 1, 2, 5, 5)
#' p <- pskewt2(x, 1, 2, 5, 5)
#' q <- qskewt2(p, 1, 2, 5, 5)
#' @name skewt2
NULL

# internal helper: convert (mean, sd) → (mu, sigma)
.skewt2_internal_params <- function(mean, sd, skew, df) {

  if (!ad_context()) {
    args <- as.list(environment())
    simulation_check(args) # informative error message if likelihood in wrong order
    # ensure sigma, df > 2
    if (any(sd <= 0)) stop("sd must be strictly positive.")
    if (any(df <= 2)) stop("df must be larger than 2 for finite variance.")
  }

  delta <- skew / sqrt(1 + skew * skew)
  b_nu <- exp(0.5 * log(df / pi) + lgamma((df - 1) / 2) - lgamma(df / 2))
  s_corr <- sqrt(df / (df - 2) - delta * delta * b_nu * b_nu)
  mu <- mean - sd * delta * b_nu / s_corr
  sigma <- sd / s_corr
  list(mu = mu, sigma = sigma)
}

#' @rdname skewt2
#' @export
#' @import RTMB
dskewt2 <- function(x, mean = 0, sd = 1, skew = 0, df = 1e2, log = FALSE) {
  p <- .skewt2_internal_params(mean, sd, skew, df)
  dskewt(x, mu = p$mu, sigma = p$sigma, skew = skew, df = df, log = log)
}

#' @rdname skewt2
#' @export
#' @usage
#' pskewt2(q, mean = 0, sd = 1, skew = 0, df = 100,
#'         method = 0, lower.tail = TRUE, log.p = FALSE)
pskewt2 <- function(q, mean = 0, sd = 1, skew = 0, df = 1e2,
                    method = 0, lower.tail = TRUE, log.p = FALSE) {
  p <- .skewt2_internal_params(mean, sd, skew, df)
  pskewt(q, mu = p$mu, sigma = p$sigma, skew = skew, df = df,
         method = method, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname skewt2
#' @export
qskewt2 <- function(p, mean = 0, sd = 1, skew = 0, df = 1e2,
                    tol = 1e-8, method = 0) {
  p2 <- .skewt2_internal_params(mean, sd, skew, df)
  qskewt(p, mu = p2$mu, sigma = p2$sigma, skew = skew, df = df,
         tol = tol, method = method)
}

#' @rdname skewt2
#' @export
rskewt2 <- function(n, mean = 0, sd = 1, skew = 0, df = 1e2) {
  p <- .skewt2_internal_params(mean, sd, skew, df)
  rskewt(n, mu = p$mu, sigma = p$sigma, skew = skew, df = df)
}
