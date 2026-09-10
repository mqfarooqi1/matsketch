# Randomized Nystrom approximation

Computes a rank-\\\ell\\ approximation \\A \approx U \hat\Lambda
U^\top\\ of a positive-semidefinite matrix from \\\ell\\ products with a
random test matrix. A tiny shift keeps the Cholesky step stable and is
removed from the eigenvalues afterwards.

## Usage

``` r
nystrom(A, l, n = NULL)
```

## Arguments

- A:

  A positive-semidefinite matrix, a function computing \\A X\\, or a
  lazy matrix from
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md)
  or
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md).

- l:

  Rank of the approximation, which is also the number of products used.

- n:

  Dimension of \\A\\. Required only when `A` is a function.

## Value

An object of class `nystrom` with the orthonormal eigenvectors `U` and
eigenvalues `values`, in decreasing order.

## References

Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized Nystrom
preconditioning. SIAM Journal on Matrix Analysis and Applications 44,
718-752. [doi:10.1137/21m1466244](https://doi.org/10.1137/21m1466244)

## See also

[`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md),
[`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(600), ncol = 3)
nys <- nystrom(kernel_matrix(X), l = 30)
head(nys$values)
#> [1] 122.601162  21.361927  17.685484  15.729190   3.594044   3.193150
```
