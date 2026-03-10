# Skellam distribution

Probability mass function, distribution function, and random generation
for the Skellam distribution.

## Usage

``` r
dskellam(x, mu1, mu2, log = FALSE)

rskellam(n, mu1, mu2)
```

## Arguments

- x:

  integer vector of counts

- mu1, mu2:

  Poisson means

- log:

  logical; return log-density if TRUE

- n:

  number of random values to return.

## Value

`dskellam` gives the probability mass function and `rskellam` generates
random deviates.

## Details

The Skellam distribution is the distribution of the difference of two
Poisson random variables. Specifically, if \\X_1 \sim
\text{Pois}(\mu_1)\\ and \\X_2 \sim \text{Pois}(\mu_2)\\, then \\X_1 -
X_2 \sim \text{Skellam}(\mu_1, \mu_2)\\.

This implementation of `dskellam` allows for automatic differentiation
with `RTMB`.

## Examples

``` r
x <- rskellam(1, 2, 3)
d <- dskellam(x, 2, 3)
```
