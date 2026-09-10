# Stochastic trace estimation

Estimates \\\mathrm{tr}(A)\\ for a matrix that is available only through
products \\A X\\, spending about `m` such products.

## Usage

``` r
trace_est(
  A,
  m,
  method = c("xtrace", "xnystrace", "hutchpp", "hutchinson"),
  n = NULL
)
```

## Arguments

- A:

  A square matrix, a function computing \\A X\\ for a matrix \\X\\, or a
  lazy kernel from
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md).

- m:

  Number of matrix-vector products to spend.

- method:

  Estimator; see Details.

- n:

  Dimension of \\A\\. Required only when `A` is a function.

## Value

An object of class `trace_est` holding the `estimate`, its `std_error`
(not available for Hutch++), the `method`, and the number of products
actually used, `matvecs`.

## Details

- `"xtrace"` (default) combines a low-rank approximation of \\A\\ with a
  correction for what it misses, and uses every product twice through a
  leave-one-out construction. It is typically far more accurate than the
  older estimators at the same cost, and it reports its own standard
  error. It works for any square matrix.

- `"xnystrace"` is the counterpart for positive-semidefinite matrices.
  It uses a Nystrom approximation, which lets it spend all `m` products
  on a single sketch.

- `"hutchpp"` is Hutch++, the estimator XTrace improves on.

- `"hutchinson"` is the classical Girard-Hutchinson estimator, the
  average of \\\omega^\top A \omega\\ over random sign vectors. It is
  included as the baseline the others are measured against.

## References

Epperly, E. N., Tropp, J. A. & Webber, R. J. (2024) XTrace: making the
most of every sample in stochastic trace estimation. SIAM Journal on
Matrix Analysis and Applications 45, 1-23.
[doi:10.1137/23m1548323](https://doi.org/10.1137/23m1548323)

Meyer, R. A., Musco, C., Musco, C. & Woodruff, D. P. (2021) Hutch++:
optimal stochastic trace estimation. Symposium on Simplicity in
Algorithms, 142-155.
[doi:10.1137/1.9781611976496.16](https://doi.org/10.1137/1.9781611976496.16)

Hutchinson, M. F. (1989) A stochastic estimator of the trace of the
influence matrix for Laplacian smoothing splines. Communications in
Statistics - Simulation and Computation 18, 1059-1076.
[doi:10.1080/03610918908812806](https://doi.org/10.1080/03610918908812806)

## See also

[`diag_est()`](https://mqfarooqi1.github.io/matsketch/reference/diag_est.md)

## Examples

``` r
set.seed(1)
U <- qr.Q(qr(matrix(rnorm(300 * 300), 300)))
A <- U %*% diag((1:300)^-2) %*% t(U)
sum(diag(A))
#> [1] 1.641606
trace_est(A, m = 40)
#> <trace_est> xtrace, 40 matrix-vector products
#>   estimate       : 1.63371
#>   standard error : 0.002271
trace_est(A, m = 40, method = "hutchinson")
#> <trace_est> hutchinson, 40 matrix-vector products
#>   estimate       : 1.62383
#>   standard error : 0.208
```
