test_that("every pivoting rule recovers a low-rank matrix exactly", {
    set.seed(1)
    A <- psd_matrix(120, rank = 8)
    for (method in c("accelerated", "simple", "greedy", "uniform")) {
        fit <- rpchol(A, k = 20, method = method)
        expect_lte(ncol(fit$F), 20)
        expect_equal(as.matrix(fit), A, tolerance = 1e-8)
        expect_lt(fit$trace_error, 1e-10)
    }
})

test_that("the reported trace error is the true relative trace error", {
    set.seed(2)
    A <- psd_matrix(150, decay = 1.5)
    for (method in c("accelerated", "simple")) {
        fit <- rpchol(A, k = 25, method = method)
        truth <- sum(diag(A - as.matrix(fit))) / sum(diag(A))
        expect_equal(fit$trace_error, truth, tolerance = 1e-10)
        expect_length(fit$trace_path, ncol(fit$F))
        expect_equal(fit$trace_path[ncol(fit$F)], fit$trace_error)
        expect_true(all(diff(fit$trace_path) <= 1e-12))
    }
})

test_that("the accelerated and simple methods are equally accurate", {
    set.seed(3)
    A <- psd_matrix(150, decay = 1.5)
    simple <- replicate(20, rpchol(A, 25, "simple")$trace_error)
    accel <- replicate(20, rpchol(A, 25, "accelerated")$trace_error)
    expect_lt(abs(mean(accel) / mean(simple) - 1), 0.2)
})

test_that("a lazy kernel is read only where needed", {
    set.seed(4)
    X <- matrix(rnorm(800), ncol = 2)
    K <- kernel_matrix(X, bandwidth = 1)
    Kd <- exp(-as.matrix(dist(X))^2 / 2)
    fit <- rpchol(K, 30)
    expect_lt(fit$entries, 0.2 * 400^2)
    expect_lt(norm(Kd - as.matrix(fit), "F") / norm(Kd, "F"), 0.05)
})

test_that("greedy pivoting is deterministic and takes outliers first", {
    set.seed(5)
    X <- rbind(matrix(rnorm(200, sd = 0.2), ncol = 2), c(50, 50), c(-50, 50))
    K <- kernel_matrix(X, bandwidth = 1)
    a <- rpchol(K, 5, "greedy")
    expect_identical(a$pivots, rpchol(K, 5, "greedy")$pivots)
    expect_true(all(c(101L, 102L) %in% a$pivots[1:3]))
})

test_that("invalid input is rejected", {
    expect_error(rpchol(matrix(1, 2, 3), 1), "square")
    expect_error(rpchol(diag(c(1, -1)), 1), "negative diagonal")
    expect_error(rpchol(diag(3), 0), "positive integer")
    expect_error(rpchol(function(x) x, 2), "not a function")
})

test_that("print and plot methods work", {
    set.seed(6)
    fit <- rpchol(psd_matrix(50), 10)
    expect_output(print(fit), "rank-10")
    expect_invisible(with_null_device(plot(fit)))
})
