# Stochastic diagonal estimation

Estimates the diagonal of a matrix that is available only through
products \\A X\\, using the XDiag estimator, which applies the same
low-rank-plus-correction and leave-one-out ideas as
[`trace_est()`](https://mqfarooqi1.github.io/matsketch/reference/trace_est.md).

## Usage

``` r
diag_est(A, m, n = NULL, adjoint = NULL)
```

## Arguments

- A:

  A square matrix, a function computing \\A X\\ for a matrix \\X\\, or a
  lazy kernel from
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md).

- m:

  Number of matrix-vector products to spend.

- n:

  Dimension of \\A\\. Required only when `A` is a function.

- adjoint:

  For a non-symmetric `A` given as a function, a function computing
  \\A^\top X\\. When omitted, `A` is assumed symmetric.

## Value

A numeric vector of length \\n\\.

## References

Epperly, E. N., Tropp, J. A. & Webber, R. J. (2024) XTrace: making the
most of every sample in stochastic trace estimation. SIAM Journal on
Matrix Analysis and Applications 45, 1-23.
[doi:10.1137/23m1548323](https://doi.org/10.1137/23m1548323)

## See also

[`trace_est()`](https://mqfarooqi1.github.io/matsketch/reference/trace_est.md)

## Examples

``` r
set.seed(1)
U <- qr.Q(qr(matrix(rnorm(200 * 200), 200)))
A <- U %*% diag((1:200)^-1.5) %*% t(U)
est <- diag_est(A, m = 60)
cor(est, diag(A))
#> [1] 0.996457
```
