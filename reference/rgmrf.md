# Sample from a multivariate Gaussian with a sparse precision matrix

Draw samples from a multivariate Gaussian distribution specified by a
sparse precision matrix. This is numerically efficient for
high-dimensional but sparse systems.

## Usage

``` r
rgmrf(n, mean = 0, Q)
```

## Arguments

- n:

  Number of samples to draw.

- mean:

  Mean vector (or scalar, which will be recycled to match the dimension
  of `Q`).

- Q:

  Sparse precision matrix (\\\Sigma^{-1}\\).

## Value

A matrix of samples with rows corresponding to samples and columns to
dimensions.

## Examples

``` r
rgmrf(3, mean = c(1, 2, 3), Q = Matrix::Diagonal(3))
#>          [,1]     [,2]     [,3]
#> [1,] 2.239496 1.891034 2.882758
#> [2,] 1.183083 3.280555 1.272729
#> [3,] 2.690184 2.503812 5.528337
```
