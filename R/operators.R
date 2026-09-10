## A matrix reaches the algorithms in one of two ways. Trace and diagonal
## estimation and the conjugate gradient method only multiply by it, so they
## accept anything that computes A %*% X. Randomly pivoted Cholesky reads
## individual entries instead, so it needs the diagonal and chosen rows.

#' @keywords internal
#' @noRd
.as_operator <- function(A, n = NULL, adjoint = NULL) {
    if (inherits(A, "matsketch_kernel")) {
        return(list(matvec = A$matvec, adjvec = A$matvec, n = A$n))
    }
    if (is.function(A)) {
        if (is.null(n)) {
            stop("`n` must be supplied when `A` is a function.", call. = FALSE)
        }
        mv <- function(X) as.matrix(A(X))
        av <- if (is.null(adjoint)) mv else function(X) as.matrix(adjoint(X))
        return(list(matvec = mv, adjvec = av, n = as.integer(n)))
    }
    if (is.matrix(A) || inherits(A, "Matrix")) {
        if (nrow(A) != ncol(A)) stop("`A` must be square.", call. = FALSE)
        return(list(matvec = function(X) as.matrix(A %*% X),
                    adjvec = function(X) as.matrix(crossprod(A, X)),
                    n = nrow(A)))
    }
    stop("`A` must be a square matrix, a function computing A %*% X, or a ",
         "lazy matrix from kernel_matrix() or grm_matrix().", call. = FALSE)
}

#' @keywords internal
#' @noRd
.as_entries <- function(A) {
    if (inherits(A, "matsketch_kernel")) return(A)
    if (is.matrix(A) || inherits(A, "Matrix")) {
        if (nrow(A) != ncol(A)) stop("`A` must be square.", call. = FALSE)
        d <- as.numeric(diag(A))
        return(list(n = nrow(A),
                    diag = function() d,
                    rows = function(i) as.matrix(A[i, , drop = FALSE]),
                    block = function(i, j) as.matrix(A[i, j, drop = FALSE])))
    }
    stop("randomly pivoted Cholesky reads entries of `A`, so it needs a ",
         "matrix or a lazy matrix from kernel_matrix() or grm_matrix(), not ",
         "a function.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.sqdist <- function(A, B) {
    d <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B)
    d[d < 0] <- 0
    d
}

#' Lazy kernel matrix
#'
#' Describes the kernel matrix \eqn{K} with entries \eqn{k(x_i, x_j)} without
#' computing it. Randomly pivoted Cholesky then evaluates only the entries it
#' needs, roughly \eqn{(k + 1) n} of the \eqn{n^2}, which is what makes it
#' practical when the full kernel matrix would not fit in memory.
#'
#' Each product with the whole matrix recomputes the kernel, block by block.
#' For methods that multiply many times, such as [pcg()], it is faster to
#' form the matrix once with `as.matrix()` whenever it fits in memory.
#'
#' @param X Numeric matrix with one row per point.
#' @param kernel Kernel family: `"gaussian"` \eqn{\exp(-r^2 / 2h^2)},
#'   `"laplace"` \eqn{\exp(-r / h)}, or the Matern kernels `"matern32"` and
#'   `"matern52"`, where \eqn{r} is the Euclidean distance.
#' @param bandwidth Length scale \eqn{h}. Defaults to the median distance
#'   between points, computed on at most the first 500 rows.
#' @param block Number of rows computed at a time when multiplying by the
#'   whole matrix, which bounds memory use.
#' @return An object of class `matsketch_kernel`, accepted wherever
#'   `matsketch` expects a matrix. `as.matrix()` forms the full matrix.
#' @seealso [rpchol()], [grm_matrix()]
#' @examples
#' X <- matrix(rnorm(400), ncol = 2)
#' K <- kernel_matrix(X, "gaussian")
#' K
#' K$block(1:3, 1:3)
#' @export
kernel_matrix <- function(X, kernel = c("gaussian", "laplace", "matern32",
                                        "matern52"),
                          bandwidth = NULL, block = 1000L) {
    kernel <- match.arg(kernel)
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    if (anyNA(X)) stop("`X` contains missing values.", call. = FALSE)
    n <- nrow(X)
    if (is.null(bandwidth)) {
        s <- X[seq_len(min(n, 500L)), , drop = FALSE]
        d <- sqrt(.sqdist(s, s))
        bandwidth <- stats::median(d[upper.tri(d)])
        if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 1
    }
    if (!is.numeric(bandwidth) || length(bandwidth) != 1L || bandwidth <= 0) {
        stop("`bandwidth` must be a single positive number.", call. = FALSE)
    }
    h <- bandwidth
    kfun <- switch(
        kernel,
        gaussian = function(d2) exp(-d2 / (2 * h^2)),
        laplace = function(d2) exp(-sqrt(d2) / h),
        matern32 = function(d2) {
            r <- sqrt(3 * d2) / h
            (1 + r) * exp(-r)
        },
        matern52 = function(d2) {
            r <- sqrt(5 * d2) / h
            (1 + r + r^2 / 3) * exp(-r)
        })
    block <- max(1L, as.integer(block))

    ## a point's distance to itself is set to exactly zero, so that the rows
    ## agree with the diagonal rather than differ from it by rounding
    krows <- function(i) {
        d2 <- .sqdist(X[i, , drop = FALSE], X)
        d2[cbind(seq_along(i), i)] <- 0
        kfun(d2)
    }
    matvec <- function(V) {
        V <- as.matrix(V)
        out <- matrix(0, n, ncol(V))
        for (start in seq(1L, n, by = block)) {
            r <- start:min(n, start + block - 1L)
            out[r, ] <- krows(r) %*% V
        }
        out
    }
    structure(
        list(n = n, kernel = kernel, bandwidth = h,
             diag = function() rep(1, n),
             rows = krows,
             block = function(i, j) {
                 d2 <- .sqdist(X[i, , drop = FALSE], X[j, , drop = FALSE])
                 d2[outer(i, j, "==")] <- 0
                 kfun(d2)
             },
             matvec = matvec),
        class = "matsketch_kernel")
}

#' @export
print.matsketch_kernel <- function(x, ...) {
    cat(sprintf("<matsketch_kernel> %s kernel on %d points, bandwidth %.4g\n",
                x$kernel, x$n, x$bandwidth))
    cat(sprintf("  full matrix would hold %s entries; nothing is stored\n",
                format(x$n^2, big.mark = ",", scientific = FALSE)))
    invisible(x)
}

#' @export
as.matrix.matsketch_kernel <- function(x, ...) x$rows(seq_len(x$n))

#' Lazy genomic relationship matrix
#'
#' Describes the genomic relationship matrix of VanRaden (2008),
#' \deqn{G = \frac{Z Z^\top}{2 \sum_j f_j (1 - f_j)},}
#' where \eqn{Z} holds the marker genotypes centred by twice the allele
#' frequencies \eqn{f_j}, without forming it. A product with \eqn{G} costs
#' two products with \eqn{Z}, and rows are computed on demand, so the
#' functions in this package can work with \eqn{G} while only the
#' \eqn{n \times p} genotypes are held in memory.
#'
#' This matters once \eqn{n} is large: for 50,000 individuals the full matrix
#' takes 20 GB, while the genotypes on a 10,000-marker panel take 4 GB.
#'
#' @param M Genotype matrix with one row per individual and one column per
#'   marker, coded as allele counts between 0 and 2. Missing genotypes must be
#'   imputed first.
#' @param freq Allele frequencies used for centring. Defaults to the observed
#'   frequencies, `colMeans(M) / 2`.
#' @return An object of class `matsketch_grm`, accepted wherever `matsketch`
#'   expects a matrix. `as.matrix()` forms the full matrix.
#' @references VanRaden, P. M. (2008) Efficient methods to compute genomic
#'   predictions. Journal of Dairy Science 91, 4414-4423.
#'   \doi{10.3168/jds.2007-0980}
#' @seealso [reml_sketch()], [sim_genomic()]
#' @examples
#' set.seed(1)
#' M <- matrix(rbinom(200 * 500, 2, 0.3), 200)
#' G <- grm_matrix(M)
#' G
#' G$block(1:3, 1:3)
#' @export
grm_matrix <- function(M, freq = NULL) {
    M <- as.matrix(M)
    if (!is.numeric(M)) {
        stop("`M` must be a numeric matrix of allele counts.", call. = FALSE)
    }
    if (anyNA(M)) {
        stop("`M` contains missing genotypes; impute them first.",
             call. = FALSE)
    }
    if (any(M < 0 | M > 2)) {
        stop("`M` must hold allele counts between 0 and 2.", call. = FALSE)
    }
    if (is.null(freq)) freq <- colMeans(M) / 2
    if (length(freq) != ncol(M) || anyNA(freq) || any(freq < 0 | freq > 1)) {
        stop("`freq` must hold one frequency between 0 and 1 per marker.",
             call. = FALSE)
    }
    denom <- 2 * sum(freq * (1 - freq))
    if (denom <= 0) stop("every marker is monomorphic.", call. = FALSE)
    Z <- sweep(M, 2L, 2 * freq)
    storage.mode(Z) <- "double"
    dg <- rowSums(Z^2) / denom
    structure(
        list(n = nrow(Z), kernel = "genomic relationship", markers = ncol(Z),
             denom = denom,
             diag = function() dg,
             rows = function(i) tcrossprod(Z[i, , drop = FALSE], Z) / denom,
             block = function(i, j) {
                 tcrossprod(Z[i, , drop = FALSE], Z[j, , drop = FALSE]) / denom
             },
             matvec = function(V) Z %*% crossprod(Z, as.matrix(V)) / denom),
        class = c("matsketch_grm", "matsketch_kernel"))
}

#' @export
print.matsketch_grm <- function(x, ...) {
    cat(sprintf("<matsketch_grm> relationships among %d individuals from %d ",
                x$n, x$markers), "markers\n", sep = "")
    cat(sprintf("  full matrix would hold %s entries; only genotypes are ",
                format(x$n^2, big.mark = ",", scientific = FALSE)),
        "stored\n", sep = "")
    invisible(x)
}
