# Tests for the skew normal distribution

test_that("skewnorm passes standard distribution checks (xi=0, omega=1, alpha=0)", {
  check_continuous_dist(
    dfun  = dskewnorm,
    pfun  = pskewnorm,
    qfun  = qskewnorm,
    xs    = c(-2, -1, 0, 1, 2),
    xi = 0, omega = 1, alpha = 0
  )
})

test_that("skewnorm passes standard distribution checks (xi=1, omega=2, alpha=3)", {
  check_continuous_dist(
    dfun  = dskewnorm,
    pfun  = pskewnorm,
    qfun  = qskewnorm,
    xs    = c(0.7, 1.6, 2.3, 3.3, 4.9),
    xi = 1, omega = 2, alpha = 3
  )
})
