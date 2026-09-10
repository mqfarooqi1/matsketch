psd_matrix <- function(n, decay = 1, rank = n) {
    U <- qr.Q(qr(matrix(stats::rnorm(n * n), n)))
    lam <- c(seq_len(rank)^(-decay), rep(0, n - rank))
    A <- U %*% (lam * t(U))
    (A + t(A)) / 2
}

## The estimators written out from their definitions, one leave-one-out
## low-rank approximation at a time: slow, but unambiguous.
xtrace_brute <- function(A, Om) {
    n <- nrow(A)
    m <- ncol(Om)
    Y <- A %*% Om
    e <- vapply(seq_len(m), function(i) {
        Q <- qr.Q(qr(Y[, -i, drop = FALSE]))
        w <- Om[, i]
        r <- w - Q %*% crossprod(Q, w)
        sum(diag(crossprod(Q, A %*% Q))) +
            (n - m + 1) / sum(r^2) * drop(crossprod(r, A %*% r))
    }, 0)
    c(estimate = mean(e), std_error = stats::sd(e) / sqrt(m))
}

xnystrace_brute <- function(A, Om) {
    n <- nrow(A)
    m <- ncol(Om)
    Y <- A %*% Om
    e <- vapply(seq_len(m), function(i) {
        O <- Om[, -i, drop = FALSE]
        Yi <- Y[, -i, drop = FALSE]
        Ahat <- Yi %*% solve(crossprod(O, Yi), t(Yi))
        w <- Om[, i]
        Q <- qr.Q(qr(O))
        r <- w - Q %*% crossprod(Q, w)
        sum(diag(Ahat)) +
            (n - m + 1) / sum(r^2) * drop(crossprod(w, (A - Ahat) %*% w))
    }, 0)
    c(estimate = mean(e), std_error = stats::sd(e) / sqrt(m))
}

xdiag_brute <- function(A, Om) {
    Y <- A %*% Om
    rowMeans(vapply(seq_len(ncol(Om)), function(i) {
        Q <- qr.Q(qr(Y[, -i, drop = FALSE]))
        w <- Om[, i]
        Aw <- A %*% w
        rowSums(Q * crossprod(A, Q)) + w * drop(Aw - Q %*% crossprod(Q, Aw))
    }, numeric(nrow(A))))
}

with_null_device <- function(code) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    force(code)
}
