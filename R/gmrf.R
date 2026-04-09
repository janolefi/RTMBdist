#' Sample from a multivariate Gaussian with a sparse precision matrix
#'
#' Draw samples from a multivariate Gaussian distribution specified by a sparse precision matrix.
#' This is numerically efficient for high-dimensional but sparse systems.
#'
#' @param n Number of samples to draw.
#' @param mean Mean vector (or scalar, which will be recycled to match the dimension of \code{Q}).
#' @param Q Sparse precision matrix (\eqn{\Sigma^{-1}}).
#'
#' @returns A matrix of samples with rows corresponding to samples and columns to dimensions.
#' @export
#'
#' @importFrom Matrix Cholesky solve
#' @importFrom stats rnorm
#'
#' @examples
#' rgmrf(3, mean = c(1, 2, 3), Q = Matrix::Diagonal(3))
rgmrf <- function(n, mean = 0, Q) {
  d <- nrow(Q)
  if (length(mean) == 1) mean <- rep(mean, d)
  stopifnot(length(mean) == d)

  # force symmetric
  Q <- Matrix::forceSymmetric(Q)

  ### efficient sampling (taken from RTMBdist:::rgmrf0)

  # safe Cholesky factorization with jitter if needed
  L <- tryCatch(Matrix::Cholesky(Q, super = TRUE, LDL = FALSE),
                error = function(e) NULL,
                warning = function(w) NULL)
  if (is.null(L)) {
    eps <- 1e-08 * mean(diag(Q))
    attempts <- 0
    while (is.null(L) && attempts < 20) {
      Q <- Q + Matrix::Diagonal(x = rep(eps, nrow(Q)))
      L <- tryCatch(Matrix::Cholesky(Q, super = TRUE, LDL = FALSE),
                    error = function(e) NULL)
      eps <- eps * 2
      attempts <- attempts + 1
    }
    warning(paste0("Precision matrix is not PD, adding jitter...\n",
                   "Required ", attempts, " attempts"))
    if (is.null(L)) stop("Matrix still not PD after 20 jitter attempts")
  }

  # L <- Matrix::Cholesky(Q, super = TRUE, LDL = FALSE)
  u <- matrix(rnorm(ncol(L) * n), ncol(L), n)
  u <- Matrix::solve(L, u, system = "Lt")
  u <- Matrix::solve(L, u, system = "Pt")

  # adding mean
  samples <- as.matrix(u) + mean

  # returning transposed (rows = samples, columns = dimensions)
  t(samples)
}
