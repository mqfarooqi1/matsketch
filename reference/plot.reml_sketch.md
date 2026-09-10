# Plot the path of a REML fit

Draws the genetic and residual variance estimates at each iteration,
which shows whether the fit settled or was still moving when it stopped.

## Usage

``` r
# S3 method for class 'reml_sketch'
plot(x, ...)
```

## Arguments

- x:

  An object from
  [`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)
  or
  [`reml_exact()`](https://mqfarooqi1.github.io/matsketch/reference/reml_exact.md).

- ...:

  Passed to
  [`graphics::matplot()`](https://rdrr.io/r/graphics/matplot.html).

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
dat <- sim_genomic(n = 200, p = 400, h2 = 0.5)
plot(reml_exact(dat$y, dat$G))
```
