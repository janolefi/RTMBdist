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
#>           [,1]     [,2]     [,3]
#> [1,] 2.5587083 2.070508 3.129288
#> [2,] 2.7150650 2.460916 1.734939
#> [3,] 0.3131471 1.554338 4.224082
```
