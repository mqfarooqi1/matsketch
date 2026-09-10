#' Randomly pivoted Cholesky
#'
#' Builds a rank-\eqn{k} approximation \eqn{A \approx F F^\top} of a
#' positive-semidefinite matrix by choosing \eqn{k} pivot columns at random,
#' each with probability proportional to the diagonal of the part of \eqn{A}
#' not yet explained.
#'
#' Sampling in proportion to the residual diagonal is what separates the
#' method from its predecessors. Greedy pivoting always takes the largest
#' residual and can fixate on outliers; uniform sampling ignores where the
#' matrix actually has mass. Randomly pivoted Cholesky balances the two, and
#' reaches near-optimal approximations while reading only about
#' \eqn{(k + 1) n} entries of the matrix, so it never needs \eqn{A} in full.
#'
#' Four pivoting rules are available:
#'
#' * `"accelerated"` (default) proposes a block of pivots at once and accepts
#'   each by rejection sampling. Its output has the same distribution as
#'   `"simple"`, but most of its work is done in block operations, which is
#'   much faster when \eqn{k} is large.
#' * `"simple"` draws one pivot at a time in proportion to the residual
#'   diagonal.
#' * `"greedy"` always takes the largest residual diagonal, which is the
#'   classical pivoted partial Cholesky decomposition.
#' * `"uniform"` takes pivots uniformly at random, which gives the
#'   column-sampling Nystrom approximation.
#'
#' The last two are included as baselines, so that the gain from random
#' pivoting can be measured on the problem at hand.
#'
#' The relative trace error reported is exact, not estimated: the residual
#' \eqn{A - F F^\top} is positive semidefinite, so its trace norm is the sum of
#' the residual diagonal the algorithm maintains anyway.
#'
#' @param A A symmetric positive-semidefinite matrix, or a lazy matrix from
#'   [kernel_matrix()] or [grm_matrix()].
#' @param k Target rank.
#' @param method Pivoting rule; see Details.
#' @param block Proposals per round for the accelerated method. Defaults to
#'   `max(10, ceiling(k / 10))`.
#' @param tol Stop early once the unexplained trace falls below `tol` times
#'   the trace of `A`. A floor of `1e-13` always applies: below it the residual
#'   is rounding error, and further pivots would be chosen from noise.
#' @return An object of class `rpchol` containing the factor `F`, whose number
#'   of columns is the rank achieved; the chosen `pivots`; the relative trace
#'   error `trace_error`, and `trace_path`, its value after each pivot; and
#'   `entries`, the number of matrix entries read.
#' @references
#' Chen, Y., Epperly, E. N., Tropp, J. A. & Webber, R. J. (2025) Randomly
#' pivoted Cholesky: practical approximation of a kernel matrix with few entry
#' evaluations. Communications on Pure and Applied Mathematics 78, 995-1041.
#' \doi{10.1002/cpa.22234}
#'
#' Epperly, E. N., Tropp, J. A. & Webber, R. J. (2025) Embrace rejection:
#' kernel matrix approximation by accelerated randomly pivoted Cholesky. SIAM
#' Journal on Matrix Analysis and Applications 46, 2527-2557.
#' \doi{10.1137/24m1699048}
#' @seealso [kernel_matrix()], [nystrom_precond()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' K <- kernel_matrix(X)
#' fit <- rpchol(K, k = 40)
#' fit
#'
#' # the same budget spent on greedy or uniform pivots
#' rpchol(K, k = 40, method = "greedy")$trace_error
#' rpchol(K, k = 40, method = "uniform")$trace_error
#' @export
rpchol <- function(A, k, method = c("accelerated", "simple", "greedy",
                                    "uniform"),
                   block = NULL, tol = 0) {
    method <- match.arg(method)
    ent <- .as_entries(A)
    n <- ent$n
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k < 1L) {
        stop("`k` must be a positive integer.", call. = FALSE)
    }
    if (!is.numeric(tol) || length(tol) != 1L || tol < 0) {
        stop("`tol` must be a single non-negative number.", call. = FALSE)
    }
    k <- min(k, n)
    d <- ent$diag()
    if (any(d < 0)) {
        stop("`A` has a negative diagonal entry, so it is not positive ",
             "semidefinite.", call. = FALSE)
    }
    tr0 <- sum(d)
    if (tr0 <= 0) stop("`A` has zero trace.", call. = FALSE)
    stop_at <- max(tol, 1e-13) * tr0

    out <- if (method == "accelerated") {
        .rpchol_accelerated(ent, n, k, d, tr0, stop_at, block)
    } else {
        .rpchol_seq(ent, n, k, d, tr0, stop_at, method)
    }
    structure(c(out, list(n = n, k = k, method = method, trace = tr0)),
              class = "rpchol")
}

## One pivot at a time. A pivot whose residual has already fallen to rounding
## level adds nothing, so it is skipped rather than divided by.
#' @keywords internal
#' @noRd
.rpchol_seq <- function(ent, n, k, d, tr0, stop_at, method) {
    F <- matrix(0, n, k)
    piv <- integer(k)
    path <- numeric(k)
    d0 <- d
    entries <- n
    got <- 0L
    ord <- if (method == "uniform") sample.int(n) else NULL
    pos <- 0L
    while (got < k && sum(d) > stop_at) {
        if (method == "uniform") {
            pos <- pos + 1L
            if (pos > n) break
            s <- ord[pos]
        } else if (method == "greedy") {
            s <- which.max(d)
        } else {
            s <- sample.int(n, 1L, prob = d)
        }
        if (d[s] <= 1e-12 * d0[s]) next
        g <- ent$rows(s)[1L, ]
        entries <- entries + n
        if (got > 0L) {
            j <- seq_len(got)
            g <- g - drop(F[, j, drop = FALSE] %*% F[s, j])
        }
        got <- got + 1L
        F[, got] <- g / sqrt(d[s])
        d <- d - F[, got]^2
        d[d < 0] <- 0
        piv[got] <- s
        path[got] <- sum(d) / tr0
    }
    list(F = F[, seq_len(got), drop = FALSE], pivots = piv[seq_len(got)],
         trace_error = sum(d) / tr0, trace_path = path[seq_len(got)],
         entries = entries)
}

## Accept or reject each proposed pivot in turn. A proposal is accepted with
## probability equal to its current Schur-complement diagonal over its
## original one, which makes the accepted pivots distributed exactly as if they
## had been drawn one at a time by the simple method.
#' @keywords internal
#' @noRd
.rejection_chol <- function(H) {
    b <- nrow(H)
    u <- diag(H)
    L <- matrix(0, b, b)
    acc <- integer(0)
    for (j in seq_len(b)) {
        hjj <- H[j, j]
        ## a repeated proposal has a residual of zero up to rounding; treat
        ## that as a rejection rather than divide by a tiny number
        if (hjj > u[j] * 1e-12 && stats::runif(1L) * u[j] < hjj) {
            acc <- c(acc, j)
            L[j:b, j] <- H[j:b, j] / sqrt(hjj)
            if (j < b) {
                r <- (j + 1L):b
                H[r, r] <- H[r, r] - tcrossprod(L[r, j])
            }
        }
    }
    list(L = L[acc, acc, drop = FALSE], idx = acc)
}

#' @keywords internal
#' @noRd
.rpchol_accelerated <- function(ent, n, k, d, tr0, stop_at, block) {
    b <- if (is.null(block)) max(10L, ceiling(k / 10)) else as.integer(block)
    if (is.na(b) || b < 1L) {
        stop("`block` must be a positive integer.", call. = FALSE)
    }
    F <- matrix(0, n, k)
    piv <- integer(k)
    path <- numeric(k)
    entries <- n
    counter <- 0L
    empty <- 0L
    while (counter < k) {
        if (sum(d) <= stop_at) break
        idx <- sample.int(n, b, replace = TRUE, prob = d)
        H <- ent$block(idx, idx)
        entries <- entries + b * b
        if (counter > 0L) {
            j <- seq_len(counter)
            H <- H - tcrossprod(F[idx, j, drop = FALSE])
        }
        rj <- .rejection_chol(H)
        acc <- rj$idx
        if (!length(acc)) {
            ## every proposal rejected: normal occasionally, but a long run
            ## means what is left is rounding error
            empty <- empty + 1L
            if (empty >= 50L) break
            next
        }
        empty <- 0L
        if (length(acc) > k - counter) {
            keep <- seq_len(k - counter)
            acc <- acc[keep]
            rj$L <- rj$L[keep, keep, drop = FALSE]
        }
        sel <- idx[acc]
        Rw <- ent$rows(sel)
        entries <- entries + length(sel) * n
        if (counter > 0L) {
            j <- seq_len(counter)
            Rw <- Rw - F[sel, j, drop = FALSE] %*% t(F[, j, drop = FALSE])
        }
        Gb <- forwardsolve(rj$L, Rw)
        cols <- counter + seq_along(sel)
        F[, cols] <- t(Gb)
        ## the trace falls by the squared norm of each new column in turn
        path[cols] <- pmax(sum(d) - cumsum(rowSums(Gb^2)), 0) / tr0
        d <- d - colSums(Gb^2)
        d[d < 0] <- 0
        path[cols[length(cols)]] <- sum(d) / tr0
        piv[cols] <- sel
        counter <- counter + length(sel)
    }
    list(F = F[, seq_len(counter), drop = FALSE],
         pivots = piv[seq_len(counter)],
         trace_error = sum(d) / tr0, trace_path = path[seq_len(counter)],
         entries = entries)
}

#' @export
print.rpchol <- function(x, ...) {
    cat(sprintf("<rpchol> rank-%d approximation of a %d x %d matrix (%s)\n",
                ncol(x$F), x$n, x$n, x$method))
    cat(sprintf("  relative trace error : %.3e\n", x$trace_error))
    cat(sprintf("  entries read         : %s of %s (%.2f%%)\n",
                format(x$entries, big.mark = ",", scientific = FALSE),
                format(x$n^2, big.mark = ",", scientific = FALSE),
                100 * x$entries / x$n^2))
    invisible(x)
}

#' @export
as.matrix.rpchol <- function(x, ...) tcrossprod(x$F)

#' @keywords internal
#' @noRd
.lowrank_eigen <- function(F) {
    s <- svd(F, nv = 0L)
    list(U = s$u, values = s$d^2)
}
