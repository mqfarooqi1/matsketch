# Variance components by exact dense REML

Fits the same model as
[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)
with the same average-information updates, but forms and factorizes
\\V\\ and computes every trace exactly. It is meant for checking
[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)
on problems small enough to factorize.

## Usage

``` r
reml_exact(y, G, X = NULL, start = NULL, tol = 1e-08, maxit = 100L)
```

## Arguments

- y:

  Numeric response.

- G:

  A positive-semidefinite relationship matrix, or a lazy matrix from
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md)
  or
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md),
  which is formed in full.

- X:

  Fixed-effect design matrix. Defaults to an intercept.

- start:

  Starting values `c(genetic, residual)`. Defaults to half the residual
  variance of a least-squares fit for each.

- tol:

  Stop when no variance component changes by more than this fraction
  between iterations.

- maxit:

  Maximum number of REML iterations.

## Value

An object of class `reml_sketch` whose `approx` is `"exact"`.

## References

Gilmour, A. R., Thompson, R. & Cullis, B. R. (1995) Average information
REML: an efficient algorithm for variance parameter estimation in linear
mixed models. Biometrics 51, 1440-1450.
[doi:10.2307/2533274](https://doi.org/10.2307/2533274)

## See also

[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)

## Examples

``` r
set.seed(1)
dat <- sim_genomic(n = 300, p = 600, h2 = 0.5)
reml_exact(dat$y, dat$G)
#> <reml_sketch> converged in 8 iterations (exact dense fit)
#>          estimate std.error
#> genetic    0.3934    0.1231
#> residual   0.6216    0.1098
#> h2         0.3876    0.1088
```
