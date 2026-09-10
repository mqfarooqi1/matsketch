# Randomly pivoted Cholesky

Builds a rank-\\k\\ approximation \\A \approx F F^\top\\ of a
positive-semidefinite matrix by choosing \\k\\ pivot columns at random,
each with probability proportional to the diagonal of the part of \\A\\
not yet explained.

## Usage

``` r
rpchol(
  A,
  k,
  method = c("accelerated", "simple", "greedy", "uniform"),
  block = NULL,
  tol = 0
)
```

## Arguments

- A:

  A symmetric positive-semidefinite matrix, or a lazy matrix from
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md)
  or
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md).

- k:

  Target rank.

- method:

  Pivoting rule; see Details.

- block:

  Proposals per round for the accelerated method. Defaults to
  `max(10, ceiling(k / 10))`.

- tol:

  Stop early once the unexplained trace falls below `tol` times the
  trace of `A`. A floor of `1e-13` always applies: below it the residual
  is rounding error, and further pivots would be chosen from noise.

## Value

An object of class `rpchol` containing the factor `F`, whose number of
columns is the rank achieved; the chosen `pivots`; the relative trace
error `trace_error`, and `trace_path`, its value after each pivot; and
`entries`, the number of matrix entries read.

## Details

Sampling in proportion to the residual diagonal is what separates the
method from its predecessors. Greedy pivoting always takes the largest
residual and can fixate on outliers; uniform sampling ignores where the
matrix actually has mass. Randomly pivoted Cholesky balances the two,
and reaches near-optimal approximations while reading only about
\\(k + 1) n\\ entries of the matrix, so it never needs \\A\\ in full.

Four pivoting rules are available:

- `"accelerated"` (default) proposes a block of pivots at once and
  accepts each by rejection sampling. Its output has the same
  distribution as `"simple"`, but most of its work is done in block
  operations, which is much faster when \\k\\ is large.

- `"simple"` draws one pivot at a time in proportion to the residual
  diagonal.

- `"greedy"` always takes the largest residual diagonal, which is the
  classical pivoted partial Cholesky decomposition.

- `"uniform"` takes pivots uniformly at random, which gives the
  column-sampling Nystrom approximation.

The last two are included as baselines, so that the gain from random
pivoting can be measured on the problem at hand.

The relative trace error reported is exact, not estimated: the residual
\\A - F F^\top\\ is positive semidefinite, so its trace norm is the sum
of the residual diagonal the algorithm maintains anyway.

## References

Chen, Y., Epperly, E. N., Tropp, J. A. & Webber, R. J. (2025) Randomly
pivoted Cholesky: practical approximation of a kernel matrix with few
entry evaluations. Communications on Pure and Applied Mathematics 78,
995-1041. [doi:10.1002/cpa.22234](https://doi.org/10.1002/cpa.22234)

Epperly, E. N., Tropp, J. A. & Webber, R. J. (2025) Embrace rejection:
kernel matrix approximation by accelerated randomly pivoted Cholesky.
SIAM Journal on Matrix Analysis and Applications 46, 2527-2557.
[doi:10.1137/24m1699048](https://doi.org/10.1137/24m1699048)

## See also

[`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md),
[`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md)

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(1000), ncol = 2)
K <- kernel_matrix(X)
fit <- rpchol(K, k = 40)
fit
#> <rpchol> rank-40 approximation of a 500 x 500 matrix (accelerated)
#>   relative trace error : 8.069e-06
#>   entries read         : 21,300 of 250,000 (8.52%)

# the same budget spent on greedy or uniform pivots
rpchol(K, k = 40, method = "greedy")$trace_error
#> [1] 5.541327e-06
rpchol(K, k = 40, method = "uniform")$trace_error
#> [1] 0.0003005479
```
