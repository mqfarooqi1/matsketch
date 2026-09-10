# Preconditioned conjugate gradients

Solves \\(A + \mu I) x = b\\ for a positive-semidefinite \\A\\ using
only products with \\A\\. Without a preconditioner this is the ordinary
conjugate gradient method, whose iteration count grows with the
condition number; with a preconditioner from
[`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md)
the count stays small and nearly independent of the size of the problem.

## Usage

``` r
pcg(A, b, mu = 0, precond = NULL, tol = 1e-08, maxit = 1000L, x0 = NULL)
```

## Arguments

- A:

  A positive-semidefinite matrix, a function computing \\A X\\, or a
  lazy matrix from
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md)
  or
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md).
  With `mu = 0`, \\A\\ must be positive definite.

- b:

  Right-hand side: a vector, or a matrix with one system per column.

- mu:

  Non-negative regularization parameter.

- precond:

  `NULL` for no preconditioning, an object from
  [`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md),
  or a function applying the inverse preconditioner to a matrix of
  residuals, column by column.

- tol:

  Stop each system when its residual norm falls below `tol` times the
  norm of its right-hand side.

- maxit:

  Maximum number of iterations.

- x0:

  Optional starting value, the same shape as `b`.

## Value

An object of class `pcg_result` with the solution `x` (the same shape as
`b`), the `iterations` and whether each system `converged`, and the
relative `residuals` after each iteration.

## Details

When `b` is a matrix, each column is solved as a separate system, but
all of them advance together, so every step multiplies \\A\\ by a block
of vectors. That is much faster in R than solving the columns one at a
time.

## References

Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized Nystrom
preconditioning. SIAM Journal on Matrix Analysis and Applications 44,
718-752. [doi:10.1137/21m1466244](https://doi.org/10.1137/21m1466244)

## See also

[`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
K <- as.matrix(kernel_matrix(X))
b <- rnorm(500)
plain <- pcg(K, b, mu = 1e-3)
pre <- nystrom_precond(rpchol(K, k = 60), mu = 1e-3)
fast <- pcg(K, b, mu = 1e-3, precond = pre)
c(plain = plain$iterations, preconditioned = fast$iterations)
#>          plain preconditioned 
#>            143              3 
```
