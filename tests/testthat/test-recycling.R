# Argument recycling for the distributions ported from gamlss.dist.
#
# gamlss.dist branches on length(nu) and uses ifelse(), whose result takes the
# length of the *test*. A single x evaluated with vectorised parameters then
# collapsed to length 1, and a long x with a shorter parameter vector was
# silently truncated. These checks pin the expected behaviour: the result has
# the length of the longest argument, and each element matches the scalar call.

# fn, first argument values, and named parameter vectors
recycling_specs <- list(
  list(f = "qpareto",    x = c(.1, .3, .5, .7, .9), pars = list(mu = c(1, 2, 3, 4, 5))),
  list(f = "qbccg",      x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(1, 2, 3, 4, 5), sigma = c(.1, .2, .3, .4, .5), nu = c(-2, -1, 0, 1, 2))),
  list(f = "qbct",       x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(1, 2, 3, 4, 5), sigma = c(.1, .2, .3, .4, .5),
                   nu = c(-2, -1, 0, 1, 2), tau = c(1, 2, 3, 4, 5))),
  list(f = "qbcpe",      x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(1, 2, 3, 4, 5), sigma = c(.1, .2, .3, .4, .5),
                   nu = c(-2, -1, 0, 1, 2), tau = c(1, 2, 3, 4, 5))),
  list(f = "qpowerexp",  x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(0, 1, 2, 3, 4), sigma = c(1, 2, 3, 4, 5), nu = c(1, 2, 3, 4, 5))),
  list(f = "qpowerexp2", x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(0, 1, 2, 3, 4), sigma = c(1, 2, 3, 4, 5), nu = c(1, 2, 3, 4, 5))),
  list(f = "qjsu",       x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(0, 1, 2, 3, 4), sigma = c(1, 2, 3, 4, 5),
                   nu = c(-2, -1, 0, 1, 2), tau = c(1, 2, 3, 4, 5))),
  list(f = "qjsu2",      x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(0, 1, 2, 3, 4), sigma = c(1, 2, 3, 4, 5),
                   nu = c(-2, -1, 0, 1, 2), tau = c(1, 2, 3, 4, 5))),
  list(f = "qexgauss",   x = c(.1, .3, .5, .7, .9),
       pars = list(mu = c(0, 1, 2, 3, 4), sigma = c(1, 2, 3, 4, 5), lambda = c(1, 2, 3, 4, 5))),
  list(f = "qgenpois",   x = c(.1, .3, .5, .7, .9),
       pars = list(lambda = c(1, 2, 3, 4, 5), phi = c(.1, .2, .3, .4, .5))),
  list(f = "pgenpois",   x = c(0, 1, 2, 3, 4),
       pars = list(lambda = c(1, 2, 3, 4, 5), phi = c(.1, .2, .3, .4, .5)))
)

# elementwise reference: call the function once per recycled argument tuple
elementwise <- function(fun, args, n) {
  vapply(seq_len(n), function(i) {
    one <- lapply(args, function(v) v[[((i - 1L) %% length(v)) + 1L]])
    as.numeric(do.call(fun, one))
  }, numeric(1))
}

for (spec in recycling_specs) {
  local({
    sp <- spec
    fun <- get(sp$f)

    test_that(paste(sp$f, "recycles a scalar first argument against vector parameters"), {
      for (pn in names(sp$pars)) {
        args <- c(list(sp$x[1]), lapply(sp$pars, function(v) v[1]))
        names(args) <- c("", names(sp$pars))
        args[[pn]] <- sp$pars[[pn]]
        n <- max(lengths(args))
        out <- do.call(fun, args)
        expect_length(out, n)
        expect_equal(as.numeric(out), elementwise(fun, args, n), tolerance = 1e-9)
      }
    })

    test_that(paste(sp$f, "does not truncate when a parameter is shorter than x"), {
      for (pn in names(sp$pars)) {
        args <- c(list(rep(sp$x, 2)), lapply(sp$pars, function(v) v[1]))
        names(args) <- c("", names(sp$pars))
        args[[pn]] <- sp$pars[[pn]]          # length 5 against x of length 10
        n <- max(lengths(args))
        out <- do.call(fun, args)
        expect_length(out, n)
        expect_equal(as.numeric(out), elementwise(fun, args, n), tolerance = 1e-9)
      }
    })
  })
}
