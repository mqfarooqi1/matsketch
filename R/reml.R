#' Variance components by sketched REML
#'
#' Fits \eqn{y = X\beta + g + e}, with \eqn{g \sim N(0, \sigma^2_g G)} and
#' \eqn{e \sim N(0, \sigma^2_e I)}, by average-information REML, without
#' forming or factorizing the covariance matrix
#' \eqn{V = \sigma^2_g G + \sigma^2_e I}.
#'
#' Three ideas from the package do the work.
#'
#' * **Solves.** A system in \eqn{V} is the system
#'   \eqn{(G + \mu I) x = b / \sigma^2_g} with \eqn{\mu = \sigma^2_e /
#'   \sigma^2_g}. It is solved by conjugate gradients, preconditioned with a
#'   low-rank approximation of \eqn{G} that is built once, by [rpchol()] or
#'   [nystrom()], and reused for every value of \eqn{\mu} the fit visits. The
#'   right-hand sides needed at each step are solved together, each starting
#'   from its solution at the previous step.
#' * **One stochastic trace.** The score needs \eqn{\mathrm{tr}(PG)}, which is
#'   estimated by XTrace or by Hutchinson's estimator. The test matrix is
#'   drawn once and kept for the whole fit, so successive iterations see the
#'   same sketch and converge smoothly rather than jittering.
#' * **One exact trace.** \eqn{PV} is idempotent with rank
#'   \eqn{n - \mathrm{rank}(X)}, so
#'   \eqn{\mathrm{tr}(P) = (n - \mathrm{rank}(X) - \sigma^2_g
#'   \mathrm{tr}(PG)) / \sigma^2_e} holds exactly, and the second trace costs
#'   nothing.
#'
#' Here \eqn{P = V^{-1} - V^{-1}X(X^\top V^{-1}X)^{-1}X^\top V^{-1}}. The
#' average-information matrix needs only quadratic forms in \eqn{P}, and these
#' are computed from solves rather than estimated.
#'
#' Because the test matrix is fixed, the fit converges to the exact solution
#' of a slightly perturbed set of REML equations. The `trace_se` column of
#' `history` shows the size of the perturbation, which shrinks as `m` grows.
#' [reml_exact()] fits the same model by dense linear algebra, for checking
#' results on problems small enough to factorize.
#'
#' @param y Numeric response.
#' @param G Relationship or kernel matrix: a positive-semidefinite matrix, a
#'   lazy matrix from [grm_matrix()] or [kernel_matrix()], or a function
#'   computing \eqn{G X} for a matrix \eqn{X}.
#' @param X Fixed-effect design matrix. Defaults to an intercept.
#' @param rank Rank of the approximation of `G` used as a preconditioner.
#' @param m Matrix-vector products per trace estimate. Each costs one linear
#'   solve at every iteration.
#' @param approx How to approximate `G`. `"auto"` uses [rpchol()] when entries
#'   of `G` can be read and [nystrom()] when `G` is a function.
#' @param estimator Trace estimator, `"xtrace"` or `"hutchinson"`. XTrace is
#'   much more accurate when \eqn{PG} has a few dominant eigenvalues, as under
#'   population structure; when the spectrum is flat the two are close.
#' @param start Starting values `c(genetic, residual)`. Defaults to half the
#'   residual variance of a least-squares fit for each.
#' @param tol Stop when no variance component changes by more than this
#'   fraction between iterations.
#' @param maxit Maximum number of REML iterations.
#' @param cg_tol Relative residual tolerance for each linear solve.
#' @return An object of class `reml_sketch` with the variance components
#'   `sigma2`, the heritability `h2`, their standard errors `se`, the
#'   fixed effects `beta`, the iteration `history`, the number of linear
#'   systems solved `solves`, and the total conjugate gradient iterations
#'   `cg_iterations` they took.
#' @references
#' Gilmour, A. R., Thompson, R. & Cullis, B. R. (1995) Average information
#' REML: an efficient algorithm for variance parameter estimation in linear
#' mixed models. Biometrics 51, 1440-1450. \doi{10.2307/2533274}
#'
#' Bermann, M., Legarra, A., Aguilar, I., Alvarez-Munera, A., Misztal, I. &
#' Lourenco, D. (2025) Estimation of (co)variance components for very large
#' datasets and complex single-step genomic models. Genetics Selection
#' Evolution 57. \doi{10.1186/s12711-025-01006-9}
#' @seealso [reml_exact()], [grm_matrix()], [sim_genomic()]
#' @examples
#' set.seed(1)
#' dat <- sim_genomic(n = 300, p = 600, h2 = 0.5)
#' fit <- reml_sketch(dat$y, dat$G, rank = 60, m = 30)
#' fit
#'
#' # the relationship matrix need never be formed
#' reml_sketch(dat$y, grm_matrix(dat$M), rank = 60, m = 30)$h2
#' @export
reml_sketch <- function(y, G, X = NULL, rank = 100L, m = 40L,
                        approx = c("auto", "rpchol", "nystrom"),
                        estimator = c("xtrace", "hutchinson"),
                        start = NULL, tol = 1e-4, maxit = 50L,
                        cg_tol = 1e-6) {
    approx <- match.arg(approx)
    estimator <- match.arg(estimator)
    mod <- .reml_setup(y, X, start)
    n <- mod$n
    opG <- .as_operator(G, n = n)
    if (opG$n != n) {
        stop("`G` must have one row and column per observation.",
             call. = FALSE)
    }
    m <- as.integer(m)
    if (length(m) != 1L || is.na(m) || m < 4L || m > n) {
        stop("`m` must be between 4 and the number of observations.",
             call. = FALSE)
    }
    if (approx == "auto") {
        approx <- if (is.function(G)) "nystrom" else "rpchol"
    }
    rank <- min(as.integer(rank), n)
    fitG <- if (approx == "rpchol") rpchol(G, rank) else
        nystrom(G, rank, n = n)
    Om <- if (estimator == "xtrace") .sphere(n, m %/% 2L) else .signs(n, m)

    Gmat <- function(M) opG$matvec(as.matrix(M))
    counts <- new.env()
    counts$solves <- 0
    counts$cg <- 0
    counts$failed <- 0
    warm <- new.env()

    make_P <- function(theta) {
        s2g <- theta[1L]
        mu <- theta[2L] / s2g
        pre <- nystrom_precond(fitG, mu)
        Amul <- function(V) Gmat(V) + mu * V
        Vinv <- function(B, slot) {
            B <- as.matrix(B)
            x0 <- warm[[slot]]
            if (!is.null(x0) && !identical(dim(x0), dim(B))) x0 <- NULL
            out <- .pcg_core(Amul, B / s2g, pre$apply, cg_tol, 1000L, x0)
            warm[[slot]] <- out$x
            counts$solves <- counts$solves + ncol(B)
            counts$cg <- counts$cg + sum(out$iterations)
            counts$failed <- counts$failed + sum(!out$converged)
            out$x
        }
        .reml_projection(Vinv, mod$X)
    }

    trace_PG <- function(Pm, theta) {
        calls <- 0L
        op <- list(n = n, matvec = function(M) {
            calls <<- calls + 1L
            Pm$Pmul(Gmat(M), paste0("trace", calls))
        })
        if (estimator == "xtrace") return(.xtrace_core(op, Om))
        v <- colSums(Om * op$matvec(Om))
        list(estimate = mean(v), std_error = stats::sd(v) / sqrt(length(v)))
    }

    fit <- .reml_iterate(mod, Gmat, make_P, tol, maxit, trace_PG)
    if (counts$failed > 0) {
        warning(sprintf(paste("%d linear solve(s) stopped at 1000 conjugate",
                              "gradient iterations before reaching `cg_tol`;",
                              "consider a larger `rank`."), counts$failed),
                call. = FALSE)
    }
    structure(c(fit, list(approx = approx, estimator = estimator,
                          rank = ncol(.eigen_of(fitG)$U),
                          m = if (estimator == "xtrace") 2L * ncol(Om) else m,
                          solves = counts$solves,
                          cg_iterations = counts$cg)),
              class = "reml_sketch")
}

#' Variance components by exact dense REML
#'
#' Fits the same model as [reml_sketch()] with the same average-information
#' updates, but forms and factorizes \eqn{V} and computes every trace
#' exactly. It is meant for checking [reml_sketch()] on problems small
#' enough to factorize.
#'
#' @inheritParams reml_sketch
#' @param G A positive-semidefinite relationship matrix, or a lazy matrix from
#'   [grm_matrix()] or [kernel_matrix()], which is formed in full.
#' @return An object of class `reml_sketch` whose `approx` is `"exact"`.
#' @references Gilmour, A. R., Thompson, R. & Cullis, B. R. (1995) Average
#'   information REML: an efficient algorithm for variance parameter
#'   estimation in linear mixed models. Biometrics 51, 1440-1450.
#'   \doi{10.2307/2533274}
#' @seealso [reml_sketch()]
#' @examples
#' set.seed(1)
#' dat <- sim_genomic(n = 300, p = 600, h2 = 0.5)
#' reml_exact(dat$y, dat$G)
#' @export
reml_exact <- function(y, G, X = NULL, start = NULL, tol = 1e-8,
                       maxit = 100L) {
    mod <- .reml_setup(y, X, start)
    G <- as.matrix(G)
    if (!identical(dim(G), c(mod$n, mod$n))) {
        stop("`G` must have one row and column per observation.",
             call. = FALSE)
    }
    Gmat <- function(M) G %*% as.matrix(M)
    make_P <- function(theta) {
        V <- theta[1L] * G
        diag(V) <- diag(V) + theta[2L]
        Vi <- chol2inv(chol(V))
        .reml_projection(function(B, slot) Vi %*% as.matrix(B), mod$X,
                         Vi = Vi)
    }
    fit <- .reml_iterate(mod, Gmat, make_P, tol, maxit, function(Pm, theta) {
        list(estimate = sum(Pm$P * G), std_error = 0)
    })
    structure(c(fit, list(approx = "exact", estimator = "exact",
                          rank = NA_integer_, m = NA_integer_,
                          solves = NA_real_, cg_iterations = NA_real_)),
              class = "reml_sketch")
}

#' @keywords internal
#' @noRd
.reml_setup <- function(y, X, start) {
    y <- as.numeric(y)
    if (anyNA(y)) stop("`y` contains missing values.", call. = FALSE)
    n <- length(y)
    X <- if (is.null(X)) matrix(1, n, 1L) else as.matrix(X)
    if (nrow(X) != n) {
        stop("`X` must have one row per observation.", call. = FALSE)
    }
    qX <- qr(X)
    if (qX$rank < ncol(X)) {
        stop("`X` does not have full column rank.", call. = FALSE)
    }
    if (is.null(start)) {
        v0 <- sum(qr.resid(qX, y)^2) / (n - qX$rank)
        start <- c(v0, v0) / 2
    }
    if (length(start) != 2L || any(!is.finite(start)) || any(start <= 0)) {
        stop("`start` must be two positive numbers.", call. = FALSE)
    }
    list(y = y, X = X, n = n, p = qX$rank, start = as.numeric(start))
}

## The REML projection P as a function, built from a solver for V that takes
## a matrix of right-hand sides and a label under which its last solution is
## kept. Passing the dense inverse as well lets the exact fit read off P.
#' @keywords internal
#' @noRd
.reml_projection <- function(Vinv, X, Vi = NULL) {
    VinvX <- Vinv(X, "X")
    XtVX <- crossprod(X, VinvX)
    XtVX <- (XtVX + t(XtVX)) / 2
    Pmul <- function(M, slot) {
        M <- as.matrix(M)
        Vinv(M, slot) - VinvX %*% solve(XtVX, crossprod(VinvX, M))
    }
    P <- if (is.null(Vi)) NULL else
        Vi - VinvX %*% solve(XtVX, t(VinvX))
    list(Pmul = Pmul, VinvX = VinvX, XtVX = XtVX, P = P)
}

## Quadratic forms and the average-information matrix at one parameter value.
#' @keywords internal
#' @noRd
.reml_ai <- function(Pm, y, Gmat) {
    Py <- drop(Pm$Pmul(y, "y"))
    GPy <- drop(Gmat(Py))
    PP <- Pm$Pmul(cbind(GPy, Py), "ai")
    PGPy <- PP[, 1L]
    PPy <- PP[, 2L]
    AI <- 0.5 * matrix(c(sum(GPy * PGPy), sum(GPy * PPy),
                         sum(GPy * PPy), sum(Py * PPy)), 2L, 2L)
    list(Py = Py, GPy = GPy, AI = AI)
}

#' @keywords internal
#' @noRd
.reml_iterate <- function(mod, Gmat, make_P, tol, maxit, trace_PG) {
    theta <- mod$start
    n <- mod$n
    hist <- vector("list", maxit)
    converged <- FALSE
    it <- 0L
    while (it < maxit) {
        it <- it + 1L
        Pm <- make_P(theta)
        q <- .reml_ai(Pm, mod$y, Gmat)
        tr <- trace_PG(Pm, theta)
        trPG <- tr$estimate
        ## tr(PV) = n - rank(X) exactly, which gives tr(P) for free
        trP <- (n - mod$p - theta[1L] * trPG) / theta[2L]
        score <- -0.5 * c(trPG - sum(q$Py * q$GPy), trP - sum(q$Py^2))
        step <- tryCatch(solve(q$AI, score), error = function(e) {
            stop("the average-information matrix is singular; the genetic ",
                 "variance may be at zero.", call. = FALSE)
        })
        new <- theta + step
        halvings <- 0L
        while (any(new <= 0) && halvings < 30L) {
            step <- step / 2
            new <- theta + step
            halvings <- halvings + 1L
        }
        new <- pmax(new, 1e-10 * sum(theta))
        change <- max(abs(new - theta) / theta)
        hist[[it]] <- data.frame(iteration = it, genetic = new[1L],
                                 residual = new[2L], trace_PG = trPG,
                                 trace_se = tr$std_error, change = change)
        theta <- new
        if (change < tol) {
            converged <- TRUE
            break
        }
    }

    Pm <- make_P(theta)
    q <- .reml_ai(Pm, mod$y, Gmat)
    cov <- solve(q$AI)
    tot <- sum(theta)
    h2 <- theta[1L] / tot
    grad <- c(theta[2L], -theta[1L]) / tot^2
    beta <- drop(solve(Pm$XtVX, crossprod(Pm$VinvX, mod$y)))
    names(beta) <- colnames(mod$X)
    list(sigma2 = c(genetic = theta[1L], residual = theta[2L]),
         h2 = h2,
         se = c(genetic = sqrt(cov[1L, 1L]), residual = sqrt(cov[2L, 2L]),
                h2 = sqrt(drop(crossprod(grad, cov %*% grad)))),
         beta = beta, converged = converged, iterations = it,
         history = do.call(rbind, hist[seq_len(it)]), n = n)
}

#' @export
print.reml_sketch <- function(x, ...) {
    how <- if (identical(x$approx, "exact")) "exact dense fit" else
        sprintf("%s rank-%d preconditioner, %s with %d products",
                x$approx, x$rank,
                if (x$estimator == "xtrace") "XTrace" else "Hutchinson", x$m)
    cat(sprintf("<reml_sketch> %s in %d iterations (%s)\n",
                if (x$converged) "converged" else "NOT converged",
                x$iterations, how))
    tab <- data.frame(estimate = c(x$sigma2, h2 = x$h2),
                      std.error = x$se[c("genetic", "residual", "h2")])
    print(round(tab, 4))
    if (!identical(x$approx, "exact")) {
        cat(sprintf("  linear systems solved : %s (mean %.1f CG iterations)\n",
                    format(x$solves, big.mark = ","),
                    x$cg_iterations / max(1, x$solves)))
    }
    invisible(x)
}

## REML log-likelihood up to a constant, used to check the fits.
#' @keywords internal
#' @noRd
.reml_loglik <- function(y, G, X, sigma2) {
    V <- sigma2[1L] * G
    diag(V) <- diag(V) + sigma2[2L]
    cV <- chol(V)
    Vi <- chol2inv(cV)
    XtVX <- crossprod(X, Vi %*% X)
    P <- Vi - Vi %*% X %*% solve(XtVX, crossprod(X, Vi))
    -0.5 * (2 * sum(log(diag(cV))) +
                as.numeric(determinant(XtVX)$modulus) +
                drop(crossprod(y, P %*% y)))
}
