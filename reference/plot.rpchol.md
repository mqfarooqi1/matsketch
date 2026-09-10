# Plot the error of a randomly pivoted Cholesky approximation

Draws the relative trace error after each pivot on a logarithmic scale,
which shows how quickly the approximation improves with its rank.

## Usage

``` r
# S3 method for class 'rpchol'
plot(x, ...)
```

## Arguments

- x:

  An object from
  [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md).

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
plot(rpchol(kernel_matrix(X), k = 60))
```
