#' matsketch: randomized matrix computations from few entries and products
#'
#' Tools for answering questions about a large positive-semidefinite matrix
#' without forming or factorizing it: a low-rank approximation built from a
#' few of its rows, its trace and diagonal from a handful of matrix-vector
#' products, and fast solutions of regularized linear systems.
#'
#' @section Main functions:
#' * [rpchol()]: low-rank approximation by randomly pivoted Cholesky.
#' * [trace_est()] and [diag_est()]: trace and diagonal estimation by XTrace,
#'   XNysTrace, XDiag and Hutch++.
#' * [nystrom()], [nystrom_precond()], [effective_dim()] and [pcg()]:
#'   preconditioned conjugate gradients.
#' * [reml_sketch()]: variance components on a relationship matrix, with
#'   [reml_exact()] as a dense reference.
#' * [kernel_matrix()] and [grm_matrix()]: kernel and genomic relationship
#'   matrices that are never formed, only evaluated where needed.
#'
#' @keywords internal
"_PACKAGE"
