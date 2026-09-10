# matsketch: randomized matrix computations from few entries and products

Tools for answering questions about a large positive-semidefinite matrix
without forming or factorizing it: a low-rank approximation built from a
few of its rows, its trace and diagonal from a handful of matrix-vector
products, and fast solutions of regularized linear systems.

## Main functions

- [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md):
  low-rank approximation by randomly pivoted Cholesky.

- [`trace_est()`](https://mqfarooqi1.github.io/matsketch/reference/trace_est.md)
  and
  [`diag_est()`](https://mqfarooqi1.github.io/matsketch/reference/diag_est.md):
  trace and diagonal estimation by XTrace, XNysTrace, XDiag and Hutch++.

- [`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md),
  [`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md),
  [`effective_dim()`](https://mqfarooqi1.github.io/matsketch/reference/effective_dim.md)
  and
  [`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md):
  preconditioned conjugate gradients.

- [`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md):
  variance components on a relationship matrix, with
  [`reml_exact()`](https://mqfarooqi1.github.io/matsketch/reference/reml_exact.md)
  as a dense reference.

- [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md)
  and
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md):
  kernel and genomic relationship matrices that are never formed, only
  evaluated where needed.

## See also

Useful links:

- <https://github.com/mqfarooqi1/matsketch>

- Report bugs at <https://github.com/mqfarooqi1/matsketch/issues>

## Author

**Maintainer**: Muhammad Farooqi <mqfarooqi@gmail.com>
([ORCID](https://orcid.org/0000-0003-4918-9791))

Authors:

- Muhammad Farooqi <mqfarooqi@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-4918-9791))
