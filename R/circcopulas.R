#' Johnson-Wehrly circular-linear copula constructor
#'
#' Returns a function computing the log density of the circular-linear copula
#' of Johnson and Wehrly (1978), intended to be used with \code{\link{dcopula}}
#' to join a circular and a linear margin, such as the turning angles and step
#' lengths of an animal track.
#'
#' @details
#' The copula density is
#' \deqn{c(u, v) = 2\pi \, g\bigl(2\pi(u - q v)\bigr),}
#' where \eqn{u} is the distribution function of the circular margin,
#' \eqn{v} that of the linear margin, \eqn{g} is a density on the circle, the
#' \emph{binding density}, and \eqn{q = \pm 1}. The joint density of an angle
#' \eqn{\theta} and a step length \eqn{s} is then
#' \deqn{f(\theta, s) = 2\pi \, g\bigl(2\pi(F_1(\theta) - q F_2(s))\bigr) f_1(\theta) f_2(s).}
#' Because \eqn{g} is periodic, \eqn{c} is a copula for any circular density
#' and any integer \eqn{q \neq 0}; only \eqn{q = \pm 1} is allowed here.
#' Any circular density of this package that allows for automatic
#' differentiation can be used for \eqn{g}, in particular \code{\link{dvm}} and
#' \code{\link{dwrpcauchy}}. Its concentration controls the strength of the
#' dependence, and \code{q = -1} reverses its direction.
#'
#' In \code{\link{dcopula}}, the circular margin must be the \emph{first} one,
#' i.e. \code{d1} and \code{p1} belong to the angle. Its distribution function
#' can be \code{\link{pwrpcauchy}}, which allows for automatic differentiation.
#' Where the circle is cut open for this distribution function only shifts
#' \eqn{u} by a constant, which the location of \eqn{g} absorbs. The likelihood
#' is therefore the same for any origin, but the meaning of the location of
#' \eqn{g} is not.
#'
#' The dependence of this copula is a helix: as \eqn{v} goes from 0 to 1, the
#' angle it favours turns once around the circle. The copula is hence not
#' symmetric in the sign of the angle, \eqn{c(u, v) \neq c(1 - u, v)}. With a
#' turning angle margin centred at 0 and the location of \eqn{g} at 0, for
#' example, medium steps are straight, short steps turn one way and long steps
#' the other. It therefore cannot capture the pattern most common in movement
#' data, where long steps are straight and short steps turn in either
#' direction. For this, see \code{\link{cfold}}.
#'
#' Random pairs from the copula are obtained by drawing \eqn{v} uniformly and
#' \eqn{z} from \eqn{g}, and setting \eqn{u = (z / (2\pi) + q v) \bmod 1};
#' see the examples.
#'
#' @param g circular density function with first argument \code{x} and an
#'   argument \code{log}, such as \code{\link{dvm}} or \code{\link{dwrpcauchy}}.
#' @param ... parameters passed to \code{g}, such as \code{mu} and \code{kappa}
#'   for \code{\link{dvm}} or \code{mu} and \code{rho} for \code{\link{dwrpcauchy}}.
#' @param q direction of the dependence, either \code{1} or \code{-1}.
#'
#' @return Function of two arguments \code{(u, v)} returning the log copula
#'   density, with \code{u} for the circular and \code{v} for the linear margin.
#'
#' @references
#' Johnson, R. A. and Wehrly, T. E. (1978) Some angular-linear distributions and
#' related regression models. Journal of the American Statistical Association,
#' 73, 602-606, doi:10.1080/01621459.1978.10480062.
#'
#' Hodel, F. H. and Fieberg, J. R. (2022) Circular-linear copulae for animal
#' movement data. Methods in Ecology and Evolution, 13,
#' doi:10.1111/2041-210X.13821.
#'
#' @seealso [dcopula()], [cfold()], [vm], [wrpcauchy]
#'
#' @export
#'
#' @examples
#' # turning angles with wrapped Cauchy margin, step lengths with Weibull margin
#' angle <- c(-2, -0.3, 0.1, 1.5); step <- c(0.2, 1.1, 2.4, 0.6)
#' d1 <- dwrpcauchy(angle, 0, 0.5, log = TRUE); p1 <- pwrpcauchy(angle, 0, 0.5)
#' d2 <- dweibull(step, 2, 1, log = TRUE); p2 <- pweibull(step, 2, 1)
#'
#' # von Mises binding density
#' dcopula(d1, d2, p1, p2, copula = cjw(dvm, mu = 0, kappa = 2), log = TRUE)
#'
#' # wrapped Cauchy binding density, dependence in the other direction
#' dcopula(d1, d2, p1, p2, copula = cjw(dwrpcauchy, mu = 0, rho = 0.7, q = -1), log = TRUE)
#'
#' # simulation from the joint distribution with von Mises binding density
#' n <- 1000
#' v <- runif(n)
#' z <- rvm(n, mu = 0, kappa = 2)
#' u <- (z / (2 * pi) + v) %% 1
#' angle <- qwrpcauchy(u, 0, 0.5)
#' step <- qweibull(v, 2, 1)
cjw <- function(g, ..., q = 1) {

  g <- match.fun(g)
  if (length(q) != 1 || !(q %in% c(-1, 1))) stop("q must be 1 or -1.")
  pars <- list(...)

  function(u, v) {
    log(2 * pi) + do.call(g, c(list(2 * pi * (u - q * v)), pars, list(log = TRUE)))
  }
}

#' Folded circular-linear copula constructor
#'
#' Turns a copula for two linear variables into a circular-linear copula that
#' is symmetric in the sign of the angle, intended to be used with
#' \code{\link{dcopula}} to join a turning angle and a step length.
#'
#' @details
#' The copula density is
#' \deqn{c(u, v) = c_0\bigl(1 - |2u - 1|, \, v\bigr),}
#' where \eqn{u} is the distribution function of the circular margin,
#' \eqn{v} that of the linear margin and \eqn{c_0} is the density of the copula
#' passed as \code{copula}, e.g. \code{\link{cgaussian}}, \code{\link{cclayton}},
#' \code{\link{cgumbel}} or \code{\link{cfrank}}. This is a copula for any
#' \eqn{c_0}. It is the rectangular patchwork copula of Hodel and Fieberg (2022),
#' with the linear copula in the rectangle \eqn{u \le 1/2} and its mirror image
#' in the rectangle \eqn{u > 1/2}.
#'
#' The map \eqn{u \mapsto 1 - |2u - 1|} folds the circle at \eqn{u = 1/2}. For a
#' circular margin that is symmetric about its mean direction and whose
#' distribution function is cut open at the antipode of that mean direction, as
#' by default in \code{\link{pwrpcauchy}} and \code{\link{pvm}}, \eqn{u = 1/2} is
#' the mean direction. \eqn{1 - |2u - 1|} is then the distribution function of
#' the angular distance to the mean direction, counted from the antipode: it is
#' 1 for angles at the mean direction and 0 for angles opposite to it. In
#' words, \eqn{c_0} links the straightness of a step to its length. A copula
#' with positive dependence, such as \code{cgaussian(rho)} with \code{rho > 0},
#' makes long steps straight and lets short steps turn in either direction,
#' the pattern most common in movement data. Unlike \code{\link{cjw}}, the
#' resulting copula is symmetric, \eqn{c(u, v) = c(1 - u, v)}.
#'
#' In \code{\link{dcopula}}, the circular margin must be the \emph{first} one,
#' i.e. \code{d1} and \code{p1} belong to the angle.
#'
#' @param copula function of two arguments \code{(u, v)} returning a log
#'   copula density for two linear variables, e.g. \code{cgaussian(0.5)}.
#'
#' @return Function of two arguments \code{(u, v)} returning the log copula
#'   density, with \code{u} for the circular and \code{v} for the linear margin.
#'
#' @references
#' Hodel, F. H. and Fieberg, J. R. (2022) Circular-linear copulae for animal
#' movement data. Methods in Ecology and Evolution, 13,
#' doi:10.1111/2041-210X.13821.
#'
#' Durante, F., Saminger-Platz, S. and Sarkoci, P. (2009) Rectangular patchwork
#' for bivariate copulas and tail dependence. Communications in Statistics -
#' Theory and Methods, 38, 2515-2527, doi:10.1080/03610920802571203.
#'
#' @seealso [dcopula()], [cjw()], [vm], [wrpcauchy]
#'
#' @export
#'
#' @examples
#' # turning angles with von Mises margin, step lengths with Weibull margin
#' angle <- c(-2, -0.3, 0.1, 1.5); step <- c(0.2, 1.1, 2.4, 0.6)
#' d1 <- dvm(angle, 0, 2, log = TRUE); p1 <- pvm(angle, 0, 2)
#' d2 <- dweibull(step, 2, 1, log = TRUE); p2 <- pweibull(step, 2, 1)
#' dcopula(d1, d2, p1, p2, copula = cfold(cgaussian(0.5)), log = TRUE)
#'
#' # wrapped Cauchy margin, which also allows for automatic differentiation
#' d1 <- dwrpcauchy(angle, 0, 0.5, log = TRUE); p1 <- pwrpcauchy(angle, 0, 0.5)
#' dcopula(d1, d2, p1, p2, copula = cfold(cclayton(2)), log = TRUE)
#'
#' # the copula is symmetric in the sign of the angle
#' cop <- cfold(cclayton(2))
#' cop(c(0.2, 0.8), 0.9)
cfold <- function(copula) {
  function(u, v) {
    # angles at the mean direction fold onto 1, where linear copulas are
    # unbounded, so the folded value is kept inside the unit interval
    eps <- .Machine$double.eps
    w <- 1 - abs(2 * u - 1)
    w <- pmin.ad(pmax.ad(w, eps), 1 - eps)
    copula(w, v)
  }
}
