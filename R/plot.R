#' @importFrom graphics abline legend matplot plot
#' @importFrom grDevices hcl.colors
NULL

#' Plot the error of a randomly pivoted Cholesky approximation
#'
#' Draws the relative trace error after each pivot on a logarithmic scale,
#' which shows how quickly the approximation improves with its rank.
#'
#' @param x An object from [rpchol()].
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' plot(rpchol(kernel_matrix(X), k = 60))
#' @export
plot.rpchol <- function(x, ...) {
    r <- x$trace_path
    plot(seq_along(r), r, log = "y", type = "l", lwd = 2,
         col = hcl.colors(3L, "Dark 3")[1L], xlab = "rank",
         ylab = "relative trace error",
         main = sprintf("Randomly pivoted Cholesky (%s)", x$method), ...)
    invisible(x)
}

#' Plot the convergence of a conjugate gradient solve
#'
#' Draws the relative residual after each iteration on a logarithmic scale,
#' one line per system solved.
#'
#' @param x An object from [pcg()].
#' @param tol Optional tolerance to mark with a horizontal line.
#' @param ... Passed to [graphics::matplot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(1000), ncol = 2)
#' K <- as.matrix(kernel_matrix(X))
#' plot(pcg(K, rnorm(500), mu = 1e-3), tol = 1e-8)
#' @export
plot.pcg_result <- function(x, tol = NULL, ...) {
    r <- as.matrix(x$residuals)
    cols <- hcl.colors(max(3L, ncol(r)), "Dark 3")[seq_len(ncol(r))]
    matplot(seq_len(nrow(r)), r, log = "y", type = "l", lty = 1, lwd = 2,
            col = cols, xlab = "iteration", ylab = "relative residual",
            main = sprintf("%s conjugate gradients",
                           if (x$preconditioned) "Preconditioned" else
                               "Plain"),
            ...)
    if (!is.null(tol)) abline(h = tol, lty = 2, col = "grey40")
    invisible(x)
}

#' Plot the path of a REML fit
#'
#' Draws the genetic and residual variance estimates at each iteration, which
#' shows whether the fit settled or was still moving when it stopped.
#'
#' @param x An object from [reml_sketch()] or [reml_exact()].
#' @param ... Passed to [graphics::matplot()].
#' @return `x`, invisibly.
#' @examples
#' set.seed(1)
#' dat <- sim_genomic(n = 200, p = 400, h2 = 0.5)
#' plot(reml_exact(dat$y, dat$G))
#' @export
plot.reml_sketch <- function(x, ...) {
    h <- x$history
    cols <- hcl.colors(3L, "Dark 3")[1:2]
    matplot(h$iteration, cbind(h$genetic, h$residual), type = "b", pch = 19,
            lty = 1, lwd = 2, col = cols, xlab = "iteration",
            ylab = "variance component", main = "REML iterations", ...)
    legend("right", legend = c("genetic", "residual"), col = cols, lwd = 2,
           pch = 19, bty = "n")
    invisible(x)
}
