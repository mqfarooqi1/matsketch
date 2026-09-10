test_that("grm_matrix() equals the VanRaden relationship matrix", {
    set.seed(30)
    M <- matrix(rbinom(50 * 200, 2, 0.3), 50)
    f <- colMeans(M) / 2
    Z <- sweep(M, 2, 2 * f)
    G <- tcrossprod(Z) / (2 * sum(f * (1 - f)))
    Gl <- grm_matrix(M)
    expect_equal(as.matrix(Gl), G)
    expect_equal(Gl$diag(), diag(G))
    expect_equal(Gl$block(3:5, 7:9), G[3:5, 7:9])
    V <- matrix(rnorm(100), 50)
    expect_equal(Gl$matvec(V), G %*% V)
    expect_output(print(Gl), "50 individuals from 200 markers")
})

test_that("grm_matrix() rejects bad genotypes", {
    expect_error(grm_matrix(matrix(c(0, 1, NA, 2), 2)), "missing")
    expect_error(grm_matrix(matrix(c(0, 1, 3, 2), 2)), "between 0 and 2")
    expect_error(grm_matrix(matrix(0, 3, 2)), "monomorphic")
    expect_error(grm_matrix(matrix(1, 3, 2), freq = 0.5), "one frequency")
})

test_that("kernel families match their formulas", {
    set.seed(31)
    X <- matrix(rnorm(40), ncol = 2)
    D <- unname(as.matrix(dist(X)))
    expect_equal(as.matrix(kernel_matrix(X, "gaussian", 2)),
                 exp(-D^2 / 8))
    expect_equal(as.matrix(kernel_matrix(X, "laplace", 2)), exp(-D / 2))
    r <- sqrt(3) * D / 2
    expect_equal(as.matrix(kernel_matrix(X, "matern32", 2)),
                 (1 + r) * exp(-r))
    r <- sqrt(5) * D / 2
    expect_equal(as.matrix(kernel_matrix(X, "matern52", 2)),
                 (1 + r + r^2 / 3) * exp(-r))
})

test_that("kernel rows, blocks and products agree", {
    set.seed(32)
    X <- matrix(rnorm(60), ncol = 3)
    K <- kernel_matrix(X, "laplace", block = 7)
    Kd <- as.matrix(K)
    expect_equal(diag(Kd), K$diag())
    expect_equal(K$block(c(2, 5, 5), c(5, 9)), Kd[c(2, 5, 5), c(5, 9)])
    expect_equal(K$matvec(diag(20)), Kd)
    expect_output(print(K), "laplace kernel on 20 points")
})

test_that("operators reject what they cannot use", {
    expect_error(kernel_matrix(matrix(c(1, NA), 1)), "missing")
    expect_error(kernel_matrix(matrix(1:4, 2), bandwidth = -1), "positive")
    expect_error(trace_est(list(1), 4), "square matrix")
})

test_that("sim_genomic() returns consistent pieces", {
    set.seed(33)
    dat <- sim_genomic(n = 60, p = 100, h2 = 0.3, pops = 3)
    expect_equal(dim(dat$M), c(60L, 100L))
    expect_equal(dim(dat$G), c(60L, 60L))
    expect_equal(nlevels(dat$pop), 3L)
    expect_equal(as.matrix(grm_matrix(dat$M)), dat$G)
    expect_null(sim_genomic(n = 20, p = 30, form_G = FALSE)$G)
    expect_error(sim_genomic(h2 = 2), "between 0 and 1")
    expect_error(sim_genomic(n = 5, pops = 6), "between 1 and `n`")
})
