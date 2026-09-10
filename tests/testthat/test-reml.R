test_that("exact AI-REML reaches the maximum of the REML likelihood", {
    set.seed(40)
    dat <- sim_genomic(n = 300, p = 600, h2 = 0.5)
    fit <- reml_exact(dat$y, dat$G)
    expect_true(fit$converged)
    nll <- function(par) {
        -matsketch:::.reml_loglik(dat$y, dat$G, dat$X, exp(par))
    }
    op <- stats::optim(log(fit$sigma2) + 0.3, nll,
                       control = list(reltol = 1e-12, maxit = 2000))
    expect_equal(unname(exp(op$par)), unname(fit$sigma2), tolerance = 1e-3)
    expect_equal(fit$h2, unname(fit$sigma2[1] / sum(fit$sigma2)))
})

test_that("sketched REML agrees with exact REML however G is supplied", {
    set.seed(41)
    dat <- sim_genomic(n = 400, p = 800, h2 = 0.6, pops = 3)
    exact <- reml_exact(dat$y, dat$G)
    fits <- list(
        dense = reml_sketch(dat$y, dat$G, rank = 60, m = 40),
        lazy = reml_sketch(dat$y, grm_matrix(dat$M), rank = 60, m = 40,
                           estimator = "hutchinson"),
        fun = reml_sketch(dat$y, function(V) dat$G %*% V, rank = 60, m = 40)
    )
    for (f in fits) {
        expect_true(f$converged)
        expect_lt(abs(f$h2 - exact$h2), 0.5 * exact$se[["h2"]])
        expect_equal(f$beta, exact$beta, tolerance = 0.05)
    }
    expect_equal(fits$dense$approx, "rpchol")
    expect_equal(fits$fun$approx, "nystrom")
    expect_gt(fits$dense$solves, 0)
})

test_that("REML input is checked", {
    expect_error(reml_sketch(c(1, NA, 3), diag(3)), "missing")
    expect_error(reml_exact(1:5, diag(4)), "one row and column")
    expect_error(reml_sketch(rnorm(10), diag(10), m = 2), "between 4")
    expect_error(reml_exact(rnorm(4), diag(4), X = cbind(1, rep(1, 4))),
                 "full column rank")
    expect_error(reml_exact(rnorm(4), diag(4), start = c(1, -1)),
                 "two positive")
})

test_that("print and plot methods work", {
    set.seed(42)
    dat <- sim_genomic(n = 150, p = 300, h2 = 0.5)
    fit <- reml_sketch(dat$y, dat$G, rank = 30, m = 20)
    expect_output(print(fit), "linear systems solved")
    expect_output(print(reml_exact(dat$y, dat$G)), "exact dense fit")
    expect_invisible(with_null_device(plot(fit)))
})
