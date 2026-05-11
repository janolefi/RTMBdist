# Reusable helpers for standard distribution property tests.
# testthat loads every helper-*.R file automatically before running tests,
# so these functions are available in all test files without source().

# Tests covered:
#   1. log = TRUE  — d(x, log=TRUE) == log(d(x))
#   2. round-trip  — q(p(x)) == x            (continuous only)
#   3. lower.tail  — p(x, lower.tail=FALSE) == 1 - p(x)
#   4. normalised  — integral of pdf == 1  /  sum of pmf == 1


#' Run tests 1-4 for a continuous distribution
#'
#' @param dfun  density function, e.g. dbeta2
#' @param pfun  CDF function,     e.g. pbeta2
#' @param qfun  quantile function, e.g. qbeta2
#' @param xs    numeric vector of test points inside the support
#' @param lower lower integration bound (default -Inf)
#' @param upper upper integration bound (default  Inf)
#' @param ...   named distribution parameters, e.g. mu = 0.4, phi = 5
#'
#' @note Do NOT include log, log.p, or lower.tail in ... — they are handled
#'   internally by each check.
check_continuous_dist <- function(dfun, pfun, qfun,
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

  # 2. q(p(x)) round-trip -----------------------------------------------------
  expect_equal(
    qfun(pfun(xs, ...), ...),
    xs,
    tolerance = 1e-6,
    label = "q(p(x)) round-trip"
  )

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


#' Run tests 1, 3, and 4 for a discrete distribution
#'
#' The q(p(x)) round-trip (test 2) is skipped because most discrete
#' distributions in the package do not expose a quantile function.
#'
#' @param dfun         PMF function, e.g. dzipois
#' @param pfun         CDF function, e.g. pzipois
#' @param xs_int       integer vector of test points
#' @param sum_support  integer vector spanning the full support (default 0:500)
#' @param ...          named distribution parameters
check_discrete_dist <- function(dfun, pfun,
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

  # 3. lower.tail complement --------------------------------------------------
  expect_equal(
    pfun(xs_int, ..., lower.tail = FALSE),
    1 - pfun(xs_int, ...),
    tolerance = 1e-10,
    label = "lower.tail=FALSE is complement of lower.tail=TRUE"
  )

  # 4. PMF sums to 1 ----------------------------------------------------------
  expect_equal(
    sum(dfun(sum_support, ...)),
    1,
    tolerance = 1e-5,
    label = "PMF sums to 1 over support"
  )
}
