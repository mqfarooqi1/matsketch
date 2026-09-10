## Implementations follow the reference code published with Epperly, Tropp and
## Webber (2024), translated line for line so that the estimators behave as
## the paper describes. Test matrices are drawn on the sphere of radius
## sqrt(n), which is the "improved" normalisation the authors recommend.

#' @keywords internal
#' @noRd
.cnormc <- function(M) sweep(M, 2L, sqrt(colSums(M^2)), "/")

## diag(t(A) %*% B), computed without forming the product
#' @keywords internal
#' @noRd
.dp <- function(A, B) colSums(A * B)

#' @keywords internal
#' @noRd
.sphere <- function(n, m) sqrt(n) * .cnormc(matrix(stats::rnorm(n * m), n, m))

#' @keywords internal
#' @noRd
.signs <- function(n, m) {
    matrix(sample(c(-1, 1), n * m, replace = TRUE), n, m)
}

## Thin QR decomposition. R's qr() may pivot columns, so the pivot is returned
## too; callers permute their test matrix to match, which leaves every
## estimator unchanged because each averages over the columns.
#' @keywords internal
#' @noRd
.qr_thin <- function(Y) {
    q <- qr(Y, tol = 1e-14)
    list(Q = qr.Q(q), R = qr.R(q), pivot = q$pivot)
}

#' @keywords internal
#' @noRd
.xtrace_core <- function(op, Om) {
    n <- op$n
    Y <- op$matvec(Om)
    qd <- .qr_thin(Y)
    Q <- qd$Q
    R <- qd$R
    Om <- Om[, qd$pivot, drop = FALSE]
    m <- ncol(Om)
    Z <- op$matvec(Q)

    W <- crossprod(Q, Om)
    S <- .cnormc(t(backsolve(R, diag(m))))
    scale <- (n - m + 1) /
        (n - colSums(W^2) + (abs(.dp(S, W) * sqrt(colSums(S^2))))^2)

    H <- crossprod(Q, Z)
    HW <- H %*% W
    T <- crossprod(Z, Om)
    dSW <- .dp(S, W)
    dSHS <- .dp(S, H %*% S)
    dTW <- .dp(T, W)
    dWHW <- .dp(W, HW)
    dSRmHW <- .dp(S, R - HW)
    dTmHRS <- .dp(T - crossprod(H, W), S)

    ests <- sum(diag(H)) - dSHS +
        (-dTW + dWHW + dSW * dSRmHW + dSW^2 * dSHS + dTmHRS * dSW) * scale
    list(estimate = mean(ests), std_error = stats::sd(ests) / sqrt(m),
         matvecs = 2L * m)
}

#' @keywords internal
#' @noRd
.xnystrace_core <- function(op, Om) {
    n <- op$n
    Y <- op$matvec(Om)
    nu <- .Machine$double.eps * sqrt(sum(Y^2)) / sqrt(n)
    Y <- Y + nu * Om
    qd <- .qr_thin(Y)
    Q <- qd$Q
    R <- qd$R
    Om <- Om[, qd$pivot, drop = FALSE]
    Y <- Y[, qd$pivot, drop = FALSE]
    m <- ncol(Om)

    H <- crossprod(Om, Y)
    C <- chol((H + t(H)) / 2)
    Cinv <- backsolve(C, diag(m))
    B <- R %*% Cinv

    qo <- .qr_thin(Om)
    WW <- crossprod(qo$Q, Om[, qo$pivot, drop = FALSE])
    SS <- .cnormc(t(backsolve(qo$R, diag(m))))
    scale <- (n - m + 1) /
        (n - colSums(WW^2) + (abs(.dp(SS, WW) * sqrt(colSums(SS^2))))^2)
    scale <- scale[order(qo$pivot)]

    W <- crossprod(Q, Om)
    Hinv_diag <- diag(chol2inv(C))
    S <- sweep(B %*% t(Cinv), 2L, Hinv_diag^(-0.5), "*")
    dSW <- .dp(S, W)
    ests <- sum(B^2) - colSums(S^2) + dSW^2 * scale - nu * n
    list(estimate = mean(ests), std_error = stats::sd(ests) / sqrt(m),
         matvecs = m)
}

#' @keywords internal
#' @noRd
.hutchpp_core <- function(op, m) {
    n <- op$n
    s <- ceiling(m / 3)
    g <- floor(m / 3)
    Q <- .qr_thin(op$matvec(.signs(n, s)))$Q
    G <- .signs(n, g)
    G <- G - Q %*% crossprod(Q, G)
    est <- sum(diag(crossprod(Q, op$matvec(Q)))) +
        sum(diag(crossprod(G, op$matvec(G)))) / g
    list(estimate = est, std_error = NA_real_, matvecs = as.integer(2 * s + g))
}

#' @keywords internal
#' @noRd
.hutch_core <- function(op, m) {
    Om <- .signs(op$n, m)
    vals <- colSums(Om * op$matvec(Om))
    list(estimate = mean(vals), std_error = stats::sd(vals) / sqrt(m),
         matvecs = as.integer(m))
}

#' Stochastic trace estimation
#'
#' Estimates \eqn{\mathrm{tr}(A)} for a matrix that is available only through
#' products \eqn{A X}, spending about `m` such products.
#'
#' * `"xtrace"` (default) combines a low-rank approximation of \eqn{A} with a
#'   correction for what it misses, and uses every product twice through a
#'   leave-one-out construction. It is typically far more accurate than the
#'   older estimators at the same cost, and it reports its own standard
#'   error. It works for any square matrix.
#' * `"xnystrace"` is the counterpart for positive-semidefinite matrices. It
#'   uses a Nystrom approximation, which lets it spend all `m` products on a
#'   single sketch.
#' * `"hutchpp"` is Hutch++, the estimator XTrace improves on.
#' * `"hutchinson"` is the classical Girard-Hutchinson estimator, the average
#'   of \eqn{\omega^\top A \omega} over random sign vectors. It is included as
#'   the baseline the others are measured against.
#'
#' @param A A square matrix, a function computing \eqn{A X} for a matrix
#'   \eqn{X}, or a lazy kernel from [kernel_matrix()].
#' @param m Number of matrix-vector products to spend.
#' @param method Estimator; see Details.
#' @param n Dimension of \eqn{A}. Required only when `A` is a function.
#' @return An object of class `trace_est` holding the `estimate`, its
#'   `std_error` (not available for Hutch++), the `method`, and the number of
#'   products actually used, `matvecs`.
#' @references
#' Epperly, E. N., Tropp, J. A. & Webber, R. J. (2024) XTrace: making the most
#' of every sample in stochastic trace estimation. SIAM Journal on Matrix
#' Analysis and Applications 45, 1-23. \doi{10.1137/23m1548323}
#'
#' Meyer, R. A., Musco, C., Musco, C. & Woodruff, D. P. (2021) Hutch++:
#' optimal stochastic trace estimation. Symposium on Simplicity in Algorithms,
#' 142-155. \doi{10.1137/1.9781611976496.16}
#'
#' Hutchinson, M. F. (1989) A stochastic estimator of the trace of the
#' influence matrix for Laplacian smoothing splines. Communications in
#' Statistics - Simulation and Computation 18, 1059-1076.
#' \doi{10.1080/03610918908812806}
#' @seealso [diag_est()]
#' @examples
#' set.seed(1)
#' U <- qr.Q(qr(matrix(rnorm(300 * 300), 300)))
#' A <- U %*% diag((1:300)^-2) %*% t(U)
#' sum(diag(A))
#' trace_est(A, m = 40)
#' trace_est(A, m = 40, method = "hutchinson")
#' @export
trace_est <- function(A, m, method = c("xtrace", "xnystrace", "hutchpp",
                                       "hutchinson"),
                      n = NULL) {
    method <- match.arg(method)
    op <- .as_operator(A, n = n)
    m <- as.integer(m)
    need <- c(xtrace = 4L, xnystrace = 2L, hutchpp = 3L, hutchinson = 2L)
    if (length(m) != 1L || is.na(m) || m < need[[method]]) {
        stop(sprintf("`m` must be at least %d for method '%s'.",
                     need[[method]], method), call. = FALSE)
    }
    if (m > op$n) {
        stop("`m` cannot exceed the dimension of `A`.", call. = FALSE)
    }
    res <- switch(method,
                  xtrace = .xtrace_core(op, .sphere(op$n, m %/% 2L)),
                  xnystrace = .xnystrace_core(op, .sphere(op$n, m)),
                  hutchpp = .hutchpp_core(op, m),
                  hutchinson = .hutch_core(op, m))
    structure(c(res, list(method = method, n = op$n)), class = "trace_est")
}

#' @export
print.trace_est <- function(x, ...) {
    se <- if (is.na(x$std_error)) "not available" else
        sprintf("%.4g", x$std_error)
    cat(sprintf("<trace_est> %s, %d matrix-vector products\n", x$method,
                x$matvecs))
    cat(sprintf("  estimate       : %.6g\n", x$estimate))
    cat(sprintf("  standard error : %s\n", se))
    invisible(x)
}

#' @keywords internal
#' @noRd
.xdiag_core <- function(op, Om) {
    Y <- op$matvec(Om)
    qd <- .qr_thin(Y)
    Q <- qd$Q
    R <- qd$R
    Om <- Om[, qd$pivot, drop = FALSE]
    Y <- Y[, qd$pivot, drop = FALSE]
    m <- ncol(Om)

    Z <- op$adjvec(Q)
    T <- crossprod(Z, Om)
    S <- .cnormc(t(backsolve(R, diag(m))))
    QS <- Q %*% S
    dQZ <- rowSums(Q * Z)
    dQSSZ <- rowSums(QS * (Z %*% S))
    dOmQT <- rowSums(Om * (Q %*% T))
    dOmY <- rowSums(Om * Y)
    dOmQSST <- rowSums(Om * sweep(QS, 2L, .dp(S, T), "*"))
    dQZ + (-dQSSZ + dOmY - dOmQT + dOmQSST) / m
}

#' Stochastic diagonal estimation
#'
#' Estimates the diagonal of a matrix that is available only through
#' products \eqn{A X}, using the XDiag estimator, which applies the same
#' low-rank-plus-correction and leave-one-out ideas as [trace_est()].
#'
#' @inheritParams trace_est
#' @param adjoint For a non-symmetric `A` given as a function, a function
#'   computing \eqn{A^\top X}. When omitted, `A` is assumed symmetric.
#' @return A numeric vector of length \eqn{n}.
#' @references Epperly, E. N., Tropp, J. A. & Webber, R. J. (2024) XTrace:
#'   making the most of every sample in stochastic trace estimation. SIAM
#'   Journal on Matrix Analysis and Applications 45, 1-23.
#'   \doi{10.1137/23m1548323}
#' @seealso [trace_est()]
#' @examples
#' set.seed(1)
#' U <- qr.Q(qr(matrix(rnorm(200 * 200), 200)))
#' A <- U %*% diag((1:200)^-1.5) %*% t(U)
#' est <- diag_est(A, m = 60)
#' cor(est, diag(A))
#' @export
diag_est <- function(A, m, n = NULL, adjoint = NULL) {
    op <- .as_operator(A, n = n, adjoint = adjoint)
    m <- as.integer(m)
    if (length(m) != 1L || is.na(m) || m < 4L) {
        stop("`m` must be at least 4.", call. = FALSE)
    }
    if (m > op$n) stop("`m` cannot exceed the dimension of `A`.", call. = FALSE)
    .xdiag_core(op, .signs(op$n, m %/% 2L))
}
