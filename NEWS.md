# matsketch 0.1.0

* First release.
* `rpchol()`: randomly pivoted Cholesky, in its simple and accelerated
  (rejection-sampling) forms, with greedy and uniform pivoting as baselines.
* `trace_est()` and `diag_est()`: the XTrace, XNysTrace and XDiag estimators,
  with Hutch++ and Girard-Hutchinson for comparison.
* `nystrom()`, `nystrom_precond()`, `effective_dim()` and `pcg()`: randomized
  Nystrom preconditioning for conjugate gradients, solving many right-hand
  sides at once.
* `reml_sketch()` and `reml_exact()`: variance components by
  average-information REML, with the covariance matrix never formed in the
  sketched fit.
* `kernel_matrix()` and `grm_matrix()`: kernel and genomic relationship
  matrices evaluated only where needed.
