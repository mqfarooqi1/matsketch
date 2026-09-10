test_that("nystrom() matches the Nystrom formula", {
    set.seed(20)
    A <- psd_matrix(90)
    set.seed(7)
    ny <- nystrom(A, 15)
    set.seed(7)
    Om <- qr.Q(qr(matrix(rnorm(90 * 15), 90)))
    Y <- A %*% Om
    Ahat <- Y %*% solve(crossprod(Om, Y), t(Y))
    expect_equal(ny$U %*% (ny$values * t(ny$U)), Ahat, tolerance = 1e-10)
    expect_true(all(diff(ny$values) <= 0))
})

test_that("the preconditioner applies the inverse of its matrix", {
    set.seed(21)
    A <- psd_matrix(60, rank = 10)
    mu <- 0.1
    pre <- nystrom_precond(nystrom(A, 10), mu)
    U <- pre$U
    lam <- pre$values
    P <- U %*% ((lam + mu) * t(U)) / (min(lam) + mu) + diag(60) - tcrossprod(U)
    expect_equal(pre$apply(P), diag(60), tolerance = 1e-10)
})

test_that("pcg() solves regularized systems, alone and in blocks", {
    set.seed(22)
    A <- psd_matrix(150, decay = 2)
    B <- matrix(rnorm(150 * 3), 150)
    mu <- 1e-4
    exact <- solve(A + mu * diag(150), B)
    pre <- nystrom_precond(rpchol(A, 40), mu)
    plain <- pcg(A, B, mu = mu, tol = 1e-10, maxit = 3000)
    fast <- pcg(A, B, mu = mu, precond = pre, tol = 1e-10)
    expect_true(all(plain$converged))
    expect_true(all(fast$converged))
    expect_equal(plain$x, exact, tolerance = 1e-5)
    expect_equal(fast$x, exact, tolerance = 1e-5)
    expect_lt(max(fast$iterations), max(plain$iterations) / 3)

    one <- pcg(A, B[, 1], mu = mu, precond = pre, tol = 1e-10)
    expect_null(dim(one$x))
    expect_equal(one$x, exact[, 1], tolerance = 1e-5)

    warm <- pcg(A, B[, 1], mu = mu, precond = pre, tol = 1e-10, x0 = one$x)
    expect_lte(warm$iterations, 1L)
})

test_that("a preconditioner may be a function", {
    set.seed(23)
    A <- psd_matrix(80, decay = 1.5)
    b <- rnorm(80)
    pre <- nystrom_precond(nystrom(A, 20), 1e-2)
    a <- pcg(A, b, mu = 1e-2, precond = pre)
    f <- pcg(A, b, mu = 1e-2, precond = function(R) pre$apply(R))
    expect_equal(a$x, f$x)
})

test_that("pcg() checks its input", {
    set.seed(24)
    A <- psd_matrix(40)
    pre <- nystrom_precond(rpchol(A, 5), 1e-3)
    expect_warning(pcg(A, rnorm(40), mu = 1, precond = pre), "built for mu")
    expect_error(pcg(A, rnorm(39)), "one row per row")
    expect_error(pcg(A, rnorm(40), mu = -1), "non-negative")
    expect_error(pcg(A, rnorm(40), precond = "x"), "precond")
    expect_error(pcg(A, rnorm(40), x0 = rnorm(3)), "same shape")
})

test_that("effective_dim() follows its formula", {
    set.seed(25)
    A <- psd_matrix(60, rank = 10)
    ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values[1:10]
    ed <- effective_dim(nystrom(A, 10), mu = 0.05)
    expect_equal(ed$d_eff, sum(ev / (ev + 0.05)), tolerance = 1e-8)
    expect_equal(ed$recommended_rank,
                 2L * as.integer(ceiling(1.5 * ed$d_eff)) + 1L)
    expect_false(ed$sufficient)
})

test_that("print and plot methods work", {
    set.seed(26)
    A <- psd_matrix(50)
    ny <- nystrom(A, 10)
    expect_output(print(ny), "rank-10")
    expect_output(print(nystrom_precond(ny, 0.1)), "rank 10")
    res <- pcg(A, matrix(rnorm(100), 50), mu = 0.1)
    expect_output(print(res), "2 systems")
    expect_invisible(with_null_device(plot(res, tol = 1e-8)))
})
