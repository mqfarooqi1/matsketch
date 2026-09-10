# Changelog

## matsketch 0.1.0

- First release.
- [`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md):
  randomly pivoted Cholesky, in its simple and accelerated
  (rejection-sampling) forms, with greedy and uniform pivoting as
  baselines.
- [`trace_est()`](https://mqfarooqi1.github.io/matsketch/reference/trace_est.md)
  and
  [`diag_est()`](https://mqfarooqi1.github.io/matsketch/reference/diag_est.md):
  the XTrace, XNysTrace and XDiag estimators, with Hutch++ and
  Girard-Hutchinson for comparison.
- [`nystrom()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom.md),
  [`nystrom_precond()`](https://mqfarooqi1.github.io/matsketch/reference/nystrom_precond.md),
  [`effective_dim()`](https://mqfarooqi1.github.io/matsketch/reference/effective_dim.md)
  and
  [`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md):
  randomized Nystrom preconditioning for conjugate gradients, solving
  many right-hand sides at once.
- [`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)
  and
  [`reml_exact()`](https://mqfarooqi1.github.io/matsketch/reference/reml_exact.md):
  variance components by average-information REML, with the covariance
  matrix never formed in the sketched fit.
- [`kernel_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/kernel_matrix.md)
  and
  [`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md):
  kernel and genomic relationship matrices evaluated only where needed.
