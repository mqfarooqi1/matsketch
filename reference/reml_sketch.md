# Variance components by sketched REML

Fits \\y = X\beta + g + e\\, with \\g \sim N(0, \sigma^2_g G)\\ and \\e
\sim N(0, \sigma^2_e I)\\, by average-information REML, without forming
or factorizing the covariance matrix \\V = \sigma^2_g G + \sigma^2_e
I\\.

## Usage

``` r
reml_sketch(
  y,
  G,
  X = NULL,
  rank = 100L,
  m = 40L,
  approx = c("auto", "rpchol", "nystrom"),
  estimator = c("xtrace", "hutchinson"),
  start = NULL,
  tol = 1e-04,
  maxit = 50L,
  cg_tol = 1e-06
)
```

## Arguments

- y:

  Numeric response.

- G:

  Relationship or kernel matrix: a positive-semidefinite matrix, a lazy
  matrix from
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md)
  or
  [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md),
  or a function computing \\G X\\ for a matrix \\X\\.

- X:

  Fixed-effect design matrix. Defaults to an intercept.

- rank:

  Rank of the approximation of `G` used as a preconditioner.

- m:

  Matrix-vector products per trace estimate. Each costs one linear solve
  at every iteration.

- approx:

  How to approximate `G`. `"auto"` uses
  [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md)
  when entries of `G` can be read and
  [`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md)
  when `G` is a function.

- estimator:

  Trace estimator, `"xtrace"` or `"hutchinson"`. XTrace is much more
  accurate when \\PG\\ has a few dominant eigenvalues, as under
  population structure; when the spectrum is flat the two are close.

- start:

  Starting values `c(genetic, residual)`. Defaults to half the residual
  variance of a least-squares fit for each.

- tol:

  Stop when no variance component changes by more than this fraction
  between iterations.

- maxit:

  Maximum number of REML iterations.

- cg_tol:

  Relative residual tolerance for each linear solve.

## Value

An object of class `reml_sketch` with the variance components `sigma2`,
the heritability `h2`, their standard errors `se`, the fixed effects
`beta`, the iteration `history`, the number of linear systems solved
`solves`, and the total conjugate gradient iterations `cg_iterations`
they took.

## Details

Three ideas from the package do the work.

- **Solves.** A system in \\V\\ is the system \\(G + \mu I) x = b /
  \sigma^2_g\\ with \\\mu = \sigma^2_e / \sigma^2_g\\. It is solved by
  conjugate gradients, preconditioned with a low-rank approximation of
  \\G\\ that is built once, by
  [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md)
  or
  [`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md),
  and reused for every value of \\\mu\\ the fit visits. The right-hand
  sides needed at each step are solved together, each starting from its
  solution at the previous step.

- **One stochastic trace.** The score needs \\\mathrm{tr}(PG)\\, which
  is estimated by XTrace or by Hutchinson's estimator. The test matrix
  is drawn once and kept for the whole fit, so successive iterations see
  the same sketch and converge smoothly rather than jittering.

- **One exact trace.** \\PV\\ is idempotent with rank \\n -
  \mathrm{rank}(X)\\, so \\\mathrm{tr}(P) = (n - \mathrm{rank}(X) -
  \sigma^2_g \mathrm{tr}(PG)) / \sigma^2_e\\ holds exactly, and the
  second trace costs nothing.

Here \\P = V^{-1} - V^{-1}X(X^\top V^{-1}X)^{-1}X^\top V^{-1}\\. The
average-information matrix needs only quadratic forms in \\P\\, and
these are computed from solves rather than estimated.

Because the test matrix is fixed, the fit converges to the exact
solution of a slightly perturbed set of REML equations. The `trace_se`
column of `history` shows the size of the perturbation, which shrinks as
`m` grows.
[`reml_exact()`](https://mqfarooqi1.github.io/matsketch/reference/reml_exact.md)
fits the same model by dense linear algebra, for checking results on
problems small enough to factorize.

## References

Gilmour, A. R., Thompson, R. & Cullis, B. R. (1995) Average information
REML: an efficient algorithm for variance parameter estimation in linear
mixed models. Biometrics 51, 1440-1450.
[doi:10.2307/2533274](https://doi.org/10.2307/2533274)

Bermann, M., Legarra, A., Aguilar, I., Alvarez-Munera, A., Misztal, I. &
Lourenco, D. (2025) Estimation of (co)variance components for very large
datasets and complex single-step genomic models. Genetics Selection
Evolution 57.
[doi:10.1186/s12711-025-01006-9](https://doi.org/10.1186/s12711-025-01006-9)

## See also

[`reml_exact()`](https://mqfarooqi1.github.io/matsketch/reference/reml_exact.md),
[`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md),
[`sim_genomic()`](https://mqfarooqi1.github.io/matsketch/reference/sim_genomic.md)

## Examples

``` r
set.seed(1)
dat <- sim_genomic(n = 300, p = 600, h2 = 0.5)
fit <- reml_sketch(dat$y, dat$G, rank = 60, m = 30)
fit
#> <reml_sketch> converged in 5 iterations (rpchol rank-60 preconditioner, XTrace with 30 products)
#>          estimate std.error
#> genetic    0.4100    0.1242
#> residual   0.6081    0.1083
#> h2         0.4027    0.1082
#>   linear systems solved : 174 (mean 6.2 CG iterations)

# the relationship matrix need never be formed
reml_sketch(dat$y, grm_matrix(dat$M), rank = 60, m = 30)$h2
#> [1] 0.3878727
```
