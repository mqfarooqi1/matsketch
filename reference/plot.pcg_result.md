# Plot the convergence of a conjugate gradient solve

Draws the relative residual after each iteration on a logarithmic scale,
one line per system solved.

## Usage

``` r
# S3 method for class 'pcg_result'
plot(x, tol = NULL, ...)
```

## Arguments

- x:

  An object from
  [`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md).

- tol:

  Optional tolerance to mark with a horizontal line.

- ...:

  Passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
K <- as.matrix(kernel_matrix(X))
plot(pcg(K, rnorm(500), mu = 1e-3), tol = 1e-8)
```
