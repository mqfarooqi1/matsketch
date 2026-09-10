# Effective dimension and recommended preconditioner rank

Estimates \\d\_{\mathrm{eff}}(\mu) = \sum_j \lambda_j / (\lambda_j +
\mu)\\ from the eigenvalues of a low-rank approximation, and the rank
\\2 \lceil 1.5\\ d\_{\mathrm{eff}} \rceil + 1\\ that Frangella, Tropp
and Udell show is enough for a well-conditioned preconditioner.

## Usage

``` r
effective_dim(approx, mu)
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

A list with `d_eff`, `recommended_rank`, and `sufficient`, which is
`TRUE` when `approx` already has at least the recommended rank.

## Details

The estimate uses only the retained eigenvalues, so it can only
understate the true effective dimension. If the recommended rank exceeds
the rank of `approx`, build a larger approximation and ask again.

## See also

[`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
effective_dim(nystrom(kernel_matrix(X), l = 60), mu = 1e-2)
#> $d_eff
#> [1] 25.95871
#> 
#> $recommended_rank
#> [1] 79
#> 
#> $sufficient
#> [1] FALSE
#> 
```
