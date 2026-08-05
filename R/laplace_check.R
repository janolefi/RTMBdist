#' Assess the accuracy of the Laplace approximation
#'
#' Diagnoses how well the Laplace approximation used by \code{TMB}/\code{RTMB}
#' approximates the marginal likelihood, via the importance-sampling check of
#' evaluating the integrand relative to the fitted Gaussian.
#'
#' For fixed parameters \eqn{\hat\theta}, the marginal likelihood is
#' \eqn{Z = \int f(y\mid x) f(x)\, dx}. Writing the importance weight
#' \deqn{\tilde\ell(x) = f(y\mid x) f(x) / N(x\mid \hat x, H^{-1}),}
#' the Laplace approximation equals \eqn{\tilde\ell(\hat x)}, while
#' \eqn{E_{x\sim N(\hat x, H^{-1})}[\tilde\ell(x)] = Z} exactly. The two coincide
#' if and only if the joint negative log-likelihood is quadratic in the latent
#' variables. The spread of the log-weights therefore measures the departure from
#' this ideal: a standard deviation near zero indicates a near-exact
#' approximation.
#'
#' @param obj A fitted \code{TMB}/\code{RTMB} object (as returned by
#'   \code{MakeADFun}) with random effects. The model is assumed to have been
#'   optimised, so that \code{obj$env$last.par.best} holds the maximum-likelihood
#'   parameters and the corresponding mode of the latent variables.
#' @param nSamples Number of Monte Carlo samples drawn from the Gaussian
#'   approximation. Default \code{1000}.
#'
#' @return An object of class \code{"laplace_check"}: a list with the log-weights
#'   \code{logw}, their standard deviation \code{sd_logw}, the Laplace and
#'   importance-sampling estimates of the log marginal likelihood
#'   (\code{lZ_laplace}, \code{lZ_is}), their difference \code{log_bias}, and the
#'   relative effective sample size \code{ess_ratio}.
#'
#' @importFrom stats sd
#'
#' @examples
#' # Chicken weight example; taken from RTMB Introduction vignette
#' data(ChickWeight)
#'
#' parameters <- list(
#'   mua=0,          ## Mean slope
#'   sda=1,          ## Std of slopes
#'   mub=0,          ## Mean intercept
#'   sdb=1,          ## Std of intercepts
#'   sdeps=1,        ## Residual Std
#'   a=rep(0, 50),   ## Random slope by chick
#'   b=rep(0, 50)    ## Random intercept by chick
#' )
#'
#' jnll <- function(parms) {
#'   getAll(ChickWeight, parms, warn=FALSE)
#'   ## Optional (enables extra RTMB features)
#'   weight <- OBS(weight)
#'   ## Initialize joint negative log likelihood
#'   nll <- 0
#'   ## Random slopes
#'   nll <- nll - sum(dnorm(a, mean=mua, sd=sda, log=TRUE))
#'   ## Random intercepts
#'   nll <- nll - sum(dnorm(b, mean=mub, sd=sdb, log=TRUE))
#'   ## Data
#'   predWeight <- a[Chick] * Time + b[Chick]
#'   nll <- nll - sum(dnorm(weight, predWeight, sd=sdeps, log=TRUE))
#'   ## Get predicted weight uncertainties
#'   ADREPORT(predWeight)
#'   ## Return
#'   nll
#' }
#'
#' obj <- MakeADFun(jnll, parameters, random=c("a", "b"), silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#'
#' chk <- laplace_check(obj)
#' chk
#' # Laplace approximation exact here: linear Gaussian-Gaussian example
#' @export
laplace_check <- function(obj, nSamples = 1000) {
  if (!is.list(obj) || is.null(obj$env)) {
    stop("`obj` must be a TMB/RTMB object as returned by MakeADFun().")
  }
  if (!is.numeric(nSamples) || length(nSamples) != 1L || nSamples < 1) {
    stop("`nSamples` must be a single positive integer.")
  }
  nSamples <- as.integer(nSamples)

  random_ind <- obj$env$random
  if (is.null(random_ind) || length(random_ind) == 0L) {
    stop("`obj` has no random effects; the Laplace approximation does not apply.")
  }

  p_hat <- obj$env$last.par.best
  if (is.null(p_hat)) {
    stop("`obj$env$last.par.best` is empty; optimise the model before calling.")
  }

  x_hat  <- p_hat[random_ind]
  H      <- obj$env$spHess(p_hat, random = TRUE)  # precision of x | y at the mode
  lZ_lap <- as.numeric(-obj$fn(p_hat[-random_ind]))

  samples <- rgmrf(nSamples, x_hat, H)            # draws from N(x_hat, H^-1)

  logw <- numeric(nSamples)
  par  <- p_hat
  gaussians <- dgmrf(samples, x_hat, H, log = TRUE)
  for (i in seq_len(nSamples)) {
    par[random_ind] <- samples[i, ]
    logw[i] <- -obj$env$f(par)
  }
  logw <- logw - gaussians

  logsumexp <- function (x) {
    max_x <- max(x)
    max_x + log(sum(exp(x - max_x)))
  }

  lZ_is <- logsumexp(logw) - log(nSamples)       # unbiased log marginal likelihood
  w     <- exp(logw - max(logw))
  ess   <- sum(w)^2 / sum(w^2)                     # Kish effective sample size

  structure(
    list(
      logw       = logw,
      nSamples   = nSamples,
      sd_logw    = sd(logw),
      lZ_laplace = lZ_lap,
      lZ_is      = lZ_is,
      log_bias   = lZ_is - lZ_lap,
      ess_ratio  = ess / nSamples
    ),
    class = "laplace_check"
  )
}

#' @param x An object of class \code{"laplace_check"}.
#' @param digits Number of decimal places for the diagnostic values.
#' @param ... Unused.
#' @rdname laplace_check
#' @export
print.laplace_check <- function(x, digits = 4, ...) {
  num <- function(v, d = digits) formatC(v, format = "f", digits = d)
  sci <- function(v) formatC(v, format = "e", digits = 1)

  rel_bias  <- abs(x$log_bias) / abs(x$lZ_laplace)   # bias relative to log-lik magnitude
  lik_ratio <- exp(x$log_bias)                       # IS / Laplace on the likelihood scale

  verdict <- if (x$ess_ratio > 0.5 && rel_bias < 1e-3) {
    "The Laplace approximation appears accurate."
  } else if (x$ess_ratio > 0.2) {
    "Plausible, but the check is only moderately reliable; interpret with caution."
  } else {
    paste("The check is unreliable here (very low effective sample size, e.g.",
          "due to high latent dimension); it neither confirms nor rules out",
          "an accurate approximation.")
  }

  cat("Laplace approximation check (importance sampling)\n")
  cat("=================================================\n")
  cat(sprintf("Monte Carlo samples:         %d\n\n", x$nSamples))

  cat("Log marginal likelihood\n")
  cat(sprintf("  Laplace:                   %s\n", num(x$lZ_laplace, 2)))
  cat(sprintf("  Importance sampling:       %s   (unbiased; should match Laplace)\n",
              num(x$lZ_is, 2)))
  cat(sprintf("  Bias (IS - Laplace):       %s   (log scale; 0 = exact)\n",
              num(x$log_bias)))
  cat(sprintf("    likelihood ratio:        %s   (exp(bias); 1 = no error)\n",
              num(lik_ratio, 3)))
  cat(sprintf("    relative to log-lik:     %s   (0 = exact)\n", sci(rel_bias)))

  cat("\nDiagnostics\n")
  cat(sprintf("  SD of log-weights:         %s   (0 = exact; grows with non-Gaussianity and dimension)\n",
              num(x$sd_logw)))
  cat(sprintf("  Effective sample size:     %s%%   (100%% = ideal; low may reflect high dimension)\n",
              num(100 * x$ess_ratio, 1)))

  cat("\n", verdict, "\n", sep = "")
  invisible(x)
}
