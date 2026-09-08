# Automatic-differentiation dispatch for functions that also exist in stats.
#
# RTMBdist needs AD-compatible versions of a few stats functions, but exporting
# them would mask the stats originals for anyone who loads the package. Instead
# the stats functions are promoted to S4 generics here, and an AD method is
# attached to each. The generics are deliberately NOT exported, so a user's
# stats::pt() and stats::plnorm() are untouched; inside this package a bare call
# to pt() or plnorm() dispatches on the arguments.
#
# The class scheme is RTMB's: "num" covers ordinary numeric input and "ad" covers
# advectors as well, so the more specific "num" method wins whenever nothing is
# taped, and the "ad" method is reached only when some argument is an advector.
# A trailing dot ("num.", "ad.") additionally allows the argument to be missing.
#
# The AD implementations themselves are exported under the names pt.ad() and
# plnorm.ad(), so they can also be called directly.

#' @importFrom methods setGeneric setMethod signature
#' @importFrom stats dgeom pgeom
NULL

setGeneric("pt", signature = c("q", "df"))

setMethod("pt", signature(q = "num", df = "num."), stats::pt)

setMethod("pt", signature(q = "ad", df = "ad."),
          function(q, df, ncp, lower.tail = TRUE, log.p = FALSE) {
            if (!missing(ncp)) {
              stop("non-centrality is not supported under automatic differentiation")
            }
            p <- pt.ad(q, df)
            if (!lower.tail) p <- 1 - p
            if (log.p) p <- log(p)
            p
          })

setGeneric("plnorm", signature = c("q", "meanlog", "sdlog"))

setMethod("plnorm", signature(q = "num", meanlog = "num.", sdlog = "num."), stats::plnorm)

# the AD implementations are wrapped rather than passed by name, so that they
# are looked up at call time and the files may be collated in any order
setMethod("plnorm", signature(q = "ad", meanlog = "ad.", sdlog = "ad."),
          function(q, meanlog = 0, sdlog = 1, lower.tail = TRUE, log.p = FALSE) {
            plnorm.ad(q, meanlog = meanlog, sdlog = sdlog,
                      lower.tail = lower.tail, log.p = log.p)
          })

setGeneric("dgeom", signature = c("x", "prob"))

setMethod("dgeom", signature(x = "num", prob = "num."), stats::dgeom)

setMethod("dgeom", signature(x = "ad", prob = "ad."),
          function(x, prob, log = FALSE) dgeom.ad(x, prob, log = log))

setGeneric("pgeom", signature = c("q", "prob"))

setMethod("pgeom", signature(q = "num", prob = "num."), stats::pgeom)

setMethod("pgeom", signature(q = "ad", prob = "ad."),
          function(q, prob, lower.tail = TRUE, log.p = FALSE) {
            pgeom.ad(q, prob, lower.tail = lower.tail, log.p = log.p)
          })
