#' Randomized Nystrom approximation
#'
#' Computes a rank-\eqn{\ell} approximation \eqn{A \approx U \hat\Lambda
#' U^\top} of a positive-semidefinite matrix from \eqn{\ell} products with a
#' random test matrix. A tiny shift keeps the Cholesky step stable and is
#' removed from the eigenvalues afterwards.
#'
#' @param A A positive-semidefinite matrix, a function computing \eqn{A X},
#'   or a lazy matrix from [kernel_matrix()] or [grm_matrix()].
#' @param l Rank of the approximation, which is also the number of products
#'   used.
#' @param n Dimension of \eqn{A}. Required only when `A` is a function.
#' @return An object of class `nystrom` with the orthonormal eigenvectors `U`
#'   and eigenvalues `values`, in decreasing order.
#' @references Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized
#'   Nystrom preconditioning. SIAM Journal on Matrix Analysis and Applications
#'   44, 718-752. \doi{10.1137/21m1466244}
#' @seealso [nystrom_precond()], [rpchol()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(600), ncol = 3)
#' nys <- nystrom(kernel_matrix(X), l = 30)
#' head(nys$values)
#' @export
nystrom <- function(A, l, n = NULL) {
    op <- .as_operator(A, n = n)
    l <- as.integer(l)
    if (length(l) != 1L || is.na(l) || l < 1L || l > op$n) {
        stop("`l` must be an integer between 1 and the dimension of `A`.",
             call. = FALSE)
    }
    Om <- qr.Q(qr(matrix(stats::rnorm(op$n * l), op$n, l)))
    Y <- op$matvec(Om)
    nu <- .Machine$double.eps * sqrt(sum(Y^2))
    Yv <- Y + nu * Om
    M <- crossprod(Om, Yv)
    C <- chol((M + t(M)) / 2)
    B <- Yv %*% backsolve(C, diag(l))
    s <- svd(B, nv = 0L)
    structure(list(U = s$u, values = pmax(0, s$d^2 - nu), l = l, n = op$n),
              class = "nystrom")
}

#' @export
print.nystrom <- function(x, ...) {
    cat(sprintf("<nystrom> rank-%d approximation of a %d x %d matrix\n",
                x$l, x$n, x$n))
    cat(sprintf("  largest eigenvalue  : %.4g\n", x$values[1L]))
    cat(sprintf("  smallest retained   : %.4g\n", x$values[x$l]))
    invisible(x)
}

#' @keywords internal
#' @noRd
.eigen_of <- function(approx) {
    if (inherits(approx, "nystrom")) {
        return(list(U = approx$U, values = approx$values))
    }
    if (inherits(approx, "rpchol")) return(.lowrank_eigen(approx$F))
    stop("`approx` must come from nystrom() or rpchol().", call. = FALSE)
}

#' Nystrom preconditioner
#'
#' Builds the preconditioner for the regularized system
#' \eqn{(A + \mu I) x = b} from a low-rank approximation
#' \eqn{A \approx U \hat\Lambda U^\top}:
#' \deqn{P^{-1} = (\hat\lambda_\ell + \mu) U (\hat\Lambda + \mu I)^{-1}
#'   U^\top + (I - U U^\top),}
#' where \eqn{\hat\lambda_\ell} is the smallest retained eigenvalue.
#'
#' The preconditioned system has a small condition number once the rank
#' reaches about the effective dimension
#' \eqn{d_{\mathrm{eff}}(\mu) = \mathrm{tr}(A (A + \mu I)^{-1})}: Frangella,
#' Tropp and Udell show that a rank of \eqn{2 \lceil 1.5\, d_{\mathrm{eff}}
#' (\mu) \rceil + 1} keeps the expected condition number below 28, whatever
#' the size of the matrix. [effective_dim()] estimates that rank.
#'
#' The approximation can come from [nystrom()] or from [rpchol()]. Because
#' only \eqn{\mu} enters the formula after the approximation is built, one
#' approximation serves any number of values of \eqn{\mu}, which is what
#' makes it cheap to re-use inside an iterative fit.
#'
#' @param approx An object from [nystrom()] or [rpchol()].
#' @param mu Positive regularization parameter.
#' @return An object of class `nystrom_precond`, to pass to [pcg()].
#' @references Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized
#'   Nystrom preconditioning. SIAM Journal on Matrix Analysis and Applications
#'   44, 718-752. \doi{10.1137/21m1466244}
#' @seealso [pcg()], [effective_dim()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' K <- kernel_matrix(X)
#' pre <- nystrom_precond(nystrom(K, l = 40), mu = 1e-3)
#' pre
#' @export
nystrom_precond <- function(approx, mu) {
    e <- .eigen_of(approx)
    if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu <= 0) {
        stop("`mu` must be a single positive number.", call. = FALSE)
    }
    U <- e$U
    lam <- e$values
    lam_l <- min(lam)
    apply_inv <- function(r) {
        r <- as.matrix(r)
        Ur <- crossprod(U, r)
        (lam_l + mu) * (U %*% (Ur / (lam + mu))) + r - U %*% Ur
    }
    structure(list(U = U, values = lam, mu = mu, lambda_l = lam_l,
                   apply = apply_inv),
              class = "nystrom_precond")
}

#' @export
print.nystrom_precond <- function(x, ...) {
    cat(sprintf("<nystrom_precond> rank %d, mu = %.4g\n", ncol(x$U), x$mu))
    cat(sprintf("  smallest retained eigenvalue : %.4g\n", x$lambda_l))
    invisible(x)
}

#' Effective dimension and recommended preconditioner rank
#'
#' Estimates \eqn{d_{\mathrm{eff}}(\mu) = \sum_j \lambda_j / (\lambda_j +
#' \mu)} from the eigenvalues of a low-rank approximation, and the rank
#' \eqn{2 \lceil 1.5\, d_{\mathrm{eff}} \rceil + 1} that Frangella, Tropp and
#' Udell show is enough for a well-conditioned preconditioner.
#'
#' The estimate uses only the retained eigenvalues, so it can only
#' understate the true effective dimension. If the recommended rank exceeds
#' the rank of `approx`, build a larger approximation and ask again.
#'
#' @param approx An object from [nystrom()] or [rpchol()].
#' @param mu Positive regularization parameter.
#' @return A list with `d_eff`, `recommended_rank`, and `sufficient`, which
#'   is `TRUE` when `approx` already has at least the recommended rank.
#' @seealso [nystrom_precond()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' effective_dim(nystrom(kernel_matrix(X), l = 60), mu = 1e-2)
#' @export
effective_dim <- function(approx, mu) {
    e <- .eigen_of(approx)
    if (!is.numeric(mu) || length(mu) != 1L || mu <= 0) {
        stop("`mu` must be a single positive number.", call. = FALSE)
    }
    d <- sum(e$values / (e$values + mu))
    rec <- 2L * as.integer(ceiling(1.5 * d)) + 1L
    list(d_eff = d, recommended_rank = rec,
         sufficient = length(e$values) >= rec)
}

## Conjugate gradients run on every column of B at once, so that each step
## multiplies by a block of vectors rather than by one, and columns drop out
## as they converge. `Amul` and `Pinv` take and return matrices. Shared by
## pcg() and reml_sketch().
#' @keywords internal
#' @noRd
.pcg_core <- function(Amul, B, Pinv, tol, maxit, X0 = NULL) {
    B <- as.matrix(B)
    k <- ncol(B)
    bn <- sqrt(colSums(B^2))
    if (is.null(X0)) {
        X <- matrix(0, nrow(B), k)
        R <- B
    } else {
        X <- as.matrix(X0)
        X[, bn == 0] <- 0
        R <- B - Amul(X)
    }
    res <- matrix(NA_real_, maxit, k)
    its <- integer(k)
    act <- which(bn > 0)
    act <- act[sqrt(colSums(R[, act, drop = FALSE]^2)) > tol * bn[act]]
    P <- matrix(0, nrow(B), k)
    rz <- numeric(k)
    if (length(act)) {
        Z <- as.matrix(Pinv(R[, act, drop = FALSE]))
        P[, act] <- Z
        rz[act] <- colSums(R[, act, drop = FALSE] * Z)
    }
    it <- 0L
    while (length(act) && it < maxit) {
        it <- it + 1L
        Pa <- P[, act, drop = FALSE]
        V <- Amul(Pa)
        alpha <- rz[act] / colSums(Pa * V)
        X[, act] <- X[, act, drop = FALSE] + sweep(Pa, 2L, alpha, "*")
        Ra <- R[, act, drop = FALSE] - sweep(V, 2L, alpha, "*")
        R[, act] <- Ra
        rr <- sqrt(colSums(Ra^2)) / bn[act]
        res[it, act] <- rr
        its[act] <- it
        keep <- rr > tol
        act <- act[keep]
        if (!length(act)) break
        Ra <- Ra[, keep, drop = FALSE]
        Z <- as.matrix(Pinv(Ra))
        rzn <- colSums(Ra * Z)
        P[, act] <- Z + sweep(P[, act, drop = FALSE], 2L, rzn / rz[act], "*")
        rz[act] <- rzn
    }
    list(x = X, iterations = its, converged = !(seq_len(k) %in% act),
         residuals = res[seq_len(it), , drop = FALSE])
}

#' Preconditioned conjugate gradients
#'
#' Solves \eqn{(A + \mu I) x = b} for a positive-semidefinite \eqn{A} using
#' only products with \eqn{A}. Without a preconditioner this is the ordinary
#' conjugate gradient method, whose iteration count grows with the condition
#' number; with a preconditioner from [nystrom_precond()] the count stays
#' small and nearly independent of the size of the problem.
#'
#' When `b` is a matrix, each column is solved as a separate system, but all
#' of them advance together, so every step multiplies \eqn{A} by a block of
#' vectors. That is much faster in R than solving the columns one at a time.
#'
#' @param A A positive-semidefinite matrix, a function computing \eqn{A X},
#'   or a lazy matrix from [kernel_matrix()] or [grm_matrix()]. With
#'   `mu = 0`, \eqn{A} must be positive definite.
#' @param b Right-hand side: a vector, or a matrix with one system per
#'   column.
#' @param mu Non-negative regularization parameter.
#' @param precond `NULL` for no preconditioning, an object from
#'   [nystrom_precond()], or a function applying the inverse preconditioner
#'   to a matrix of residuals, column by column.
#' @param tol Stop each system when its residual norm falls below `tol` times
#'   the norm of its right-hand side.
#' @param maxit Maximum number of iterations.
#' @param x0 Optional starting value, the same shape as `b`.
#' @return An object of class `pcg_result` with the solution `x` (the same
#'   shape as `b`), the `iterations` and whether each system `converged`,
#'   and the relative `residuals` after each iteration.
#' @references Frangella, Z., Tropp, J. A. & Udell, M. (2023) Randomized
#'   Nystrom preconditioning. SIAM Journal on Matrix Analysis and Applications
#'   44, 718-752. \doi{10.1137/21m1466244}
#' @seealso [nystrom_precond()]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' K <- as.matrix(kernel_matrix(X))
#' b <- rnorm(500)
#' plain <- pcg(K, b, mu = 1e-3)
#' pre <- nystrom_precond(rpchol(K, k = 60), mu = 1e-3)
#' fast <- pcg(K, b, mu = 1e-3, precond = pre)
#' c(plain = plain$iterations, preconditioned = fast$iterations)
#' @export
pcg <- function(A, b, mu = 0, precond = NULL, tol = 1e-8, maxit = 1000L,
                x0 = NULL) {
    single <- is.null(dim(b))
    B <- as.matrix(b)
    storage.mode(B) <- "double"
    op <- .as_operator(A, n = nrow(B))
    if (op$n != nrow(B)) {
        stop("`b` must have one row per row of `A`.", call. = FALSE)
    }
    if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < 0) {
        stop("`mu` must be a single non-negative number.", call. = FALSE)
    }
    maxit <- as.integer(maxit)
    if (length(maxit) != 1L || is.na(maxit) || maxit < 1L) {
        stop("`maxit` must be a positive integer.", call. = FALSE)
    }
    Pinv <- if (is.null(precond)) {
        function(R) R
    } else if (inherits(precond, "nystrom_precond")) {
        if (!isTRUE(all.equal(precond$mu, mu))) {
            warning("the preconditioner was built for mu = ", precond$mu,
                    " but the system uses mu = ", mu, "; it will still be ",
                    "correct but may be less effective.", call. = FALSE)
        }
        precond$apply
    } else if (is.function(precond)) {
        function(R) as.matrix(precond(R))
    } else {
        stop("`precond` must be NULL, a nystrom_precond object or a function.",
             call. = FALSE)
    }
    X0 <- NULL
    if (!is.null(x0)) {
        X0 <- as.matrix(x0)
        if (!identical(dim(X0), dim(B))) {
            stop("`x0` must have the same shape as `b`.", call. = FALSE)
        }
    }
    Amul <- function(V) op$matvec(V) + mu * V
    out <- .pcg_core(Amul, B, Pinv, tol, maxit, X0)
    if (single) {
        out$x <- drop(out$x)
        out$residuals <- out$residuals[, 1L]
    }
    structure(c(out, list(mu = mu, preconditioned = !is.null(precond))),
              class = "pcg_result")
}

#' @export
print.pcg_result <- function(x, ...) {
    k <- length(x$iterations)
    its <- if (k == 1L) {
        sprintf("%d iterations", x$iterations)
    } else {
        sprintf("%d systems, %d to %d iterations", k, min(x$iterations),
                max(x$iterations))
    }
    status <- if (all(x$converged)) "converged" else
        sprintf("%d NOT converged", sum(!x$converged))
    cat(sprintf("<pcg_result> %s CG, %s, %s\n",
                if (x$preconditioned) "preconditioned" else "plain", its,
                status))
    r <- as.matrix(x$residuals)
    last <- apply(r, 2L, function(v) {
        v <- v[!is.na(v)]
        if (length(v)) v[length(v)] else 0
    })
    cat(sprintf("  largest final relative residual : %.3e\n",
                max(c(0, last))))
    invisible(x)
}
