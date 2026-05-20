# Tests for the Gumbel distribution

test_that("gumbel passes standard distribution checks (location=0, scale=1)", {
  check_continuous_dist(
    dfun  = dgumbel,
    pfun  = pgumbel,
    qfun  = qgumbel,
    xs    = c(-1, 0, 1, 2, 4),
    location = 0, scale = 1
  )
})

test_that("gumbel passes standard distribution checks (location=2, scale=3)", {
  check_continuous_dist(
    dfun  = dgumbel,
    pfun  = pgumbel,
    qfun  = qgumbel,
    xs    = c(-2, 0, 2, 5, 8),
    location = 2, scale = 3
  )
})
