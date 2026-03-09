# Wishart distribution

Density and random generation for the wishart distribution

## Usage

``` r
dwishart(x, nu, Sigma, log = FALSE)

rwishart(n, nu, Sigma)
```

## Arguments

- x:

  positive definite \\p \times p\\ matrix or array of such matrices of
  dimension \\p \times p \times n\\ (for \\n\\ density evaluations)

- nu:

  degrees of freedom, needs to be greater than `p - 1`

- Sigma:

  scale matrix, needs to be positive definite and match the dimension of
  `x`.

- log:

  logical; if `TRUE`, densities \\p\\ are returned as \\\log(p)\\.

- n:

  number of random deviates to return

## Value

`dwishart` gives the density, `rwishart` generates random deviates
(matrix for `n = 1`, array with `n` slices for `n > 1`)

## Examples

``` r
# single input: matrix-valued
x <- rwishart(1, nu = 5, Sigma = diag(3))
d <- dwishart(x, nu = 5, Sigma = diag(3))

# multiple inputs: array of matrices
x <- rwishart(4, nu = 5, Sigma = diag(3))
d <- dwishart(x, nu = 5, Sigma = diag(3))

# multiple inputs for x, nu and Sigma
nu <- c(7,5,8,9)
Sigma <- array(dim = c(3,3,4))
for(i in 1:4) Sigma[,,i] <- (i + 3) * diag(3)
x <- rwishart(4, nu, Sigma)
d <- dwishart(x, nu, Sigma)
```
