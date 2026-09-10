test_that("XTrace equals its leave-one-out definition", {
    set.seed(10)
    A <- psd_matrix(80)
    N <- A + matrix(rnorm(80 * 80), 80) / 9
    for (M in list(A, N)) {
        Om <- matsketch:::.sphere(80, 8)
        got <- matsketch:::.xtrace_core(matsketch:::.as_operator(M), Om)
        ref <- xtrace_brute(M, Om)
        expect_equal(got$estimate, ref[["estimate"]], tolerance = 1e-10)
        expect_equal(got$std_error, ref[["std_error"]], tolerance = 1e-8)
    }
})

test_that("XNysTrace equals its leave-one-out definition", {
    set.seed(11)
    A <- psd_matrix(80)
    for (m in c(4, 12)) {
        Om <- matsketch:::.sphere(80, m)
        got <- matsketch:::.xnystrace_core(matsketch:::.as_operator(A), Om)
        ref <- xnystrace_brute(A, Om)
        expect_equal(got$estimate, ref[["estimate"]], tolerance = 1e-10)
        expect_equal(got$std_error, ref[["std_error"]], tolerance = 1e-8)
    }
})

test_that("XDiag equals its leave-one-out definition", {
    set.seed(12)
    A <- psd_matrix(80)
    N <- A + matrix(rnorm(80 * 80), 80) / 9
    for (M in list(A, N)) {
        Om <- matsketch:::.signs(80, 8)
        got <- matsketch:::.xdiag_core(matsketch:::.as_operator(M), Om)
        expect_equal(got, xdiag_brute(M, Om), tolerance = 1e-10)
    }
})

test_that("Hutch++ follows its definition", {
    set.seed(13)
    A <- psd_matrix(60)
    set.seed(99)
    got <- trace_est(A, 12, "hutchpp")
    set.seed(99)
    S <- matrix(sample(c(-1, 1), 60 * 4, TRUE), 60)
    G <- matrix(sample(c(-1, 1), 60 * 4, TRUE), 60)
    Q <- qr.Q(qr(A %*% S))
    G <- G - Q %*% crossprod(Q, G)
    ref <- sum(diag(crossprod(Q, A %*% Q))) + sum(diag(crossprod(G, A %*% G))) / 4
    expect_equal(got$estimate, ref, tolerance = 1e-12)
    expect_equal(got$matvecs, 12L)
})

test_that("every estimator is exact on a low-rank matrix", {
    set.seed(14)
    A <- psd_matrix(100, rank = 6)
    tr <- sum(diag(A))
    expect_equal(trace_est(A, 20, "xtrace")$estimate, tr, tolerance = 1e-8)
    expect_equal(trace_est(A, 10, "xnystrace")$estimate, tr, tolerance = 1e-8)
    expect_equal(trace_est(A, 24, "hutchpp")$estimate, tr, tolerance = 1e-8)
    expect_equal(diag_est(A, 20), diag(A), tolerance = 1e-8)
})

test_that("XTrace and XNysTrace beat Hutchinson on a decaying spectrum", {
    set.seed(15)
    A <- psd_matrix(200, decay = 2)
    tr <- sum(diag(A))
    err <- function(method) {
        median(replicate(15, abs(trace_est(A, 40, method)$estimate - tr)))
    }
    hutch <- err("hutchinson")
    expect_lt(err("xtrace"), hutch / 10)
    expect_lt(err("xnystrace"), hutch / 10)
})

test_that("A may be given as a function", {
    set.seed(16)
    A <- psd_matrix(100, decay = 2)
    f <- function(X) A %*% X
    expect_error(trace_est(f, 10), "`n` must be supplied")
    set.seed(1)
    a <- trace_est(f, 20, n = 100)
    set.seed(1)
    b <- trace_est(A, 20)
    expect_equal(a$estimate, b$estimate)
})

test_that("invalid budgets are rejected", {
    expect_error(trace_est(diag(10), 3), "at least 4")
    expect_error(trace_est(diag(10), 11), "cannot exceed")
    expect_error(diag_est(diag(10), 2), "at least 4")
    expect_error(trace_est(matrix(1, 2, 3), 2, "hutchinson"), "square")
})

test_that("the print method reports the estimate", {
    set.seed(17)
    expect_output(print(trace_est(psd_matrix(30), 10)), "estimate")
    expect_output(print(trace_est(psd_matrix(30), 9, "hutchpp")),
                  "not available")
})
