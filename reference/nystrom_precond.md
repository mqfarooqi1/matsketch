# Nystrom preconditioner

Builds the preconditioner for the regularized system \\(A + \mu I) x =
b\\ from a low-rank approximation \\A \approx U \hat\Lambda U^\top\\:
\$\$P^{-1} = (\hat\lambda\_\ell + \mu) U (\hat\Lambda + \mu I)^{-1}
U^\top + (I - U U^\top),\$\$ where \\\hat\lambda\_\ell\\ is the smallest
retained eigenvalue.

## Usage

``` r
nystrom_precond(approx, mu)
```

## Arguments

- approx:

  An object from
  [`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md)
  or
  [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md).

- mu:

  Positive regularization parameter.

## Value

An object of class `nystrom_precond`, to pass to
[`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md).

## Details

The preconditioned system has a small condition number once the rank
reaches about the effective dimension \\d\_{\mathrm{eff}}(\mu) =
\mathrm{tr}(A (A + \mu I)^{-1})\\: Frangella, Tropp and Udell show that
a rank of \\2 \lceil 1.5\\ d\_{\mathrm{eff}} (\mu) \rceil + 1\\ keeps
the expected condition number below 28, whatever the size of the matrix.
[`effective_dim()`](https://mqfarooqi1.github.io/matsketch/reference/effective_dim.md)
estimates that rank.

The approximation can come from
[`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md)
or from
[`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md).
Because only \\\mu\\ enters the formula after the approximation is
built, one approximation serves any number of values of \\\mu\\, which
is what makes it cheap to re-use inside an iterative fit.

## References

Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized Nystrom
preconditioning. SIAM Journal on Matrix Analysis and Applications 44,
718-752. [doi:10.1137/21m1466244](https://doi.org/10.1137/21m1466244)

## See also

[`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md),
[`effective_dim()`](https://mqfarooqi1.github.io/matsketch/reference/effective_dim.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
K <- kernel_matrix(X)
pre <- nystrom_precond(nystrom(K, l = 40), mu = 1e-3)
pre
#> <nystrom_precond> rank 40, mu = 0.001
#>   smallest retained eigenvalue : 6.628e-05
```
