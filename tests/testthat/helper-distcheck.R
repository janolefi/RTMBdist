# Reusable helpers for standard distribution property tests.
# testthat loads every helper-*.R file automatically before running tests,
# so these functions are available in all test files without source().

# Tests covered:
#   1. log = TRUE  — d(x, log=TRUE) == log(d(x))
#   2. round-trip  — q(p(x)) == x            (continuous only)
#   3. lower.tail  — p(x, lower.tail=FALSE) == 1 - p(x)
#   4. normalised  — integral of pdf == 1  /  sum of pmf == 1
#   5. AD gradient — MakeTape(nll) Jacobian contains no NaN


#' Check that the AD gradient of a distribution's NLL is free of NaN (test 5)
#'
#' Simulates n observations from rfun using the supplied parameter values,
#' builds a negative log-likelihood with dfun inside RTMB::MakeTape, evaluates
#' the Jacobian at those same parameter values, and asserts that no gradient
#' component is NaN.
#'
#' @param .dfun density function, e.g. dllogis
#' @param .rfun random generation function, e.g. rllogis
#' @param ...   named distribution parameters at sensible true values,
#'              e.g. alpha = 1, beta = 2
#' @param .n    number of simulated observations (default 20)
#' @param .seed random seed (default 42)
#'
#' @note Parameters use a leading dot so that common distribution parameter
#'   names (df, mu, sigma, ...) cannot partially match the helper's own
#'   formals via R's partial-matching rules.
check_ad_gradient <- function(.dfun, .rfun, ..., .n = 20, .seed = 42) {
  set.seed(.seed)
  params <- list(...)
  x   <- do.call(.rfun, c(list(n = .n), params))
  par <- unlist(params)   # named numeric vector, used as tape starting point

  nll <- function(par) {
    x <- RTMB::OBS(x)
    -sum(do.call(.dfun, c(list(x = x), as.list(par), list(log = TRUE))))
  }

  environment(nll) <- environment() # make sure nll has access to x

  F <- tryCatch(
    RTMB::MakeTape(nll, par),
    error = function(e) {
      fail(paste("MakeTape failed:", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(F)) {
    grad <- F$jacobian(par)
    expect_false(
      any(is.nan(grad)),
      label = "no NaN in AD gradient"
    )
  }
}


#' Run tests 1-4 for a continuous distribution
#'
#' @param dfun  density function, e.g. dbeta2
#' @param pfun  CDF function,     e.g. pbeta2
#' @param qfun  quantile function, e.g. qbeta2; pass NULL to skip the round-trip test
#' @param xs    numeric vector of test points inside the support
#' @param lower lower integration bound (default -Inf)
#' @param upper upper integration bound (default  Inf)
#' @param ...   named distribution parameters, e.g. mu = 0.4, phi = 5
#'
#' @note Do NOT include log, log.p, or lower.tail in ... — they are handled
#'   internally by each check.
check_continuous_dist <- function(dfun, pfun, qfun = NULL,
                                  xs,
                                  lower = -Inf, upper = Inf,
                                  ...) {

  # 1. log = TRUE consistency -------------------------------------------------
  expect_equal(
    dfun(xs, ..., log = TRUE),
    log(dfun(xs, ...)),
    tolerance = 1e-10,
    label = "log=TRUE matches log(density)"
  )

  # 2. q(p(x)) round-trip (skipped when qfun = NULL) --------------------------
  if (!is.null(qfun)) {
    expect_equal(
      qfun(pfun(xs, ...), ...),
      xs,
      tolerance = 1e-6,
      label = "q(p(x)) round-trip"
    )
  }

  # 3. lower.tail complement --------------------------------------------------
  expect_equal(
    pfun(xs, ..., lower.tail = FALSE),
    1 - pfun(xs, ...),
    tolerance = 1e-10,
    label = "lower.tail=FALSE is complement of lower.tail=TRUE"
  )

  # 4. density integrates to 1 ------------------------------------------------
  # integrate() forwards ... to dfun; log and eps use their defaults.
  result <- tryCatch(
    integrate(dfun, lower = lower, upper = upper, ...),
    error = function(e) list(value = NA_real_)
  )
  expect_equal(
    result$value, 1,
    tolerance = 1e-4,
    label = "density integrates to 1"
  )
}


#' Run tests 1, 3, and 4 for a zero-inflated continuous distribution
#'
#' Zero-inflated distributions have a point mass at 0 and a continuous density
#' for x > 0, so the standard normalisation test (integrate == 1) does not apply.
#' Instead we check: d(0) + integral_of_continuous_part == 1.
#'
#' @param dfun  density function, e.g. dzigamma
#' @param pfun  CDF function,     e.g. pzigamma
#' @param xs    strictly positive test points
#' @param upper upper integration bound for the continuous part (default Inf)
#' @param ...   named distribution parameters including zeroprob
check_zeroinfl_dist <- function(dfun, pfun, xs, upper = Inf, ...) {

  # 1. log=TRUE consistency at x > 0 -----------------------------------------
  expect_equal(
    dfun(xs, ..., log = TRUE),
    log(dfun(xs, ...)),
    tolerance = 1e-10,
    label = "log=TRUE matches log(density) for x > 0"
  )

  # 3. lower.tail complement --------------------------------------------------
  expect_equal(
    pfun(xs, ..., lower.tail = FALSE),
    1 - pfun(xs, ...),
    tolerance = 1e-10,
    label = "lower.tail=FALSE is complement of lower.tail=TRUE"
  )

  # 4. point mass + continuous integral = 1 -----------------------------------
  # integrate() assigns measure zero to the single point x=0, so the integral
  # gives the (1 - zeroprob) continuous mass; adding d(0) recovers 1.
  result <- tryCatch(
    integrate(dfun, lower = 0, upper = upper, ...),
    error = function(e) list(value = NA_real_)
  )
  expect_equal(
    dfun(0, ...) + result$value,
    1,
    tolerance = 1e-4,
    label = "point mass d(0) + continuous integral = 1"
  )
}


#' Run tests 1, 3, and 4 for a mixed (inflated) continuous distribution
#'
#' Handles distributions with point masses at one or more boundary values
#' (e.g. zero-inflated, one-inflated, zero-one-inflated beta).
#' The normalisation check verifies:
#'   sum(d(point_masses)) + integral of continuous part = 1
#'
#' @param dfun         density function
#' @param pfun         CDF function
#' @param xs           test points strictly inside the continuous support
#' @param lower        lower integration bound for the continuous part (default 0)
#' @param upper        upper integration bound for the continuous part (default 1)
#' @param point_masses numeric vector of inflated-mass locations (e.g. 0, 1, c(0,1))
#' @param ...          named distribution parameters
check_inflated_dist <- function(dfun, pfun, xs,
                                lower = 0, upper = 1,
                                point_masses = 0,
                                ...) {

  # 1. log=TRUE consistency at interior points --------------------------------
  expect_equal(
    dfun(xs, ..., log = TRUE),
    log(dfun(xs, ...)),
    tolerance = 1e-10,
    label = "log=TRUE matches log(density) at interior points"
  )

  # 3. lower.tail complement --------------------------------------------------
  expect_equal(
    pfun(xs, ..., lower.tail = FALSE),
    1 - pfun(xs, ...),
    tolerance = 1e-10,
    label = "lower.tail=FALSE is complement of lower.tail=TRUE"
  )

  # 4. point masses + continuous integral = 1 ---------------------------------
  pm_total <- sum(dfun(point_masses, ...))
  result <- tryCatch(
    integrate(dfun, lower = lower, upper = upper, ...),
    error = function(e) list(value = NA_real_)
  )
  expect_equal(
    pm_total + result$value,
    1,
    tolerance = 1e-4,
    label = "sum of point masses + continuous integral = 1"
  )
}


#' Run tests 1, 3, and 4 for a discrete distribution
#'
#' The q(p(x)) round-trip (test 2) is skipped because most discrete
#' distributions in the package do not expose a quantile function.
#' Pass pfun = NULL to also skip the lower.tail complement test (test 3),
#' e.g. for distributions that have no CDF function at all.
#'
#' @param dfun         PMF function, e.g. dzipois
#' @param pfun         CDF function, e.g. pzipois; pass NULL to skip test 3
#' @param xs_int       integer vector of test points
#' @param sum_support  integer vector spanning the full support (default 0:500)
#' @param ...          named distribution parameters
check_discrete_dist <- function(dfun, pfun = NULL,
                                xs_int,
                                sum_support = 0:500,
                                ...) {

  # 1. log = TRUE consistency -------------------------------------------------
  expect_equal(
    dfun(xs_int, ..., log = TRUE),
    log(dfun(xs_int, ...)),
    tolerance = 1e-10,
    label = "log=TRUE matches log(pmf)"
  )

  # 3. lower.tail complement (skipped when pfun = NULL) -----------------------
  if (!is.null(pfun)) {
    expect_equal(
      pfun(xs_int, ..., lower.tail = FALSE),
      1 - pfun(xs_int, ...),
      tolerance = 1e-10,
      label = "lower.tail=FALSE is complement of lower.tail=TRUE"
    )
  }

  # 4. PMF sums to 1 ----------------------------------------------------------
  expect_equal(
    sum(dfun(sum_support, ...)),
    1,
    tolerance = 1e-5,
    label = "PMF sums to 1 over support"
  )
}
