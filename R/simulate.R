#' Simulate a genomic data set
#'
#' Draws biallelic marker genotypes, forms the genomic relationship matrix of
#' VanRaden (2008), and simulates a trait with the requested heritability, so
#' that [reml_sketch()] can be checked against known variance components.
#'
#' With `pops > 1` the individuals come from that many subpopulations whose
#' allele frequencies have drifted apart under the Balding-Nichols model with
#' fixation index `fst`. Population structure gives the relationship matrix a
#' few large eigenvalues, as real breeding and human cohorts do, and those are
#' exactly what a low-rank preconditioner captures.
#'
#' Genetic values are sums of marker effects on the centred genotypes, with
#' variance chosen so that the genetic and residual variance components on the
#' scale of the relationship matrix are `h2` and `1 - h2`.
#'
#' The function draws random numbers but does not set the seed; call
#' [base::set.seed()] first for a reproducible data set.
#'
#' @param n Number of individuals.
#' @param p Number of markers.
#' @param h2 Heritability, between 0 and 1.
#' @param pops Number of subpopulations.
#' @param fst Fixation index between subpopulations, used when `pops > 1`.
#' @param form_G Return the relationship matrix itself. Set to `FALSE` for
#'   large `n` and pass `grm_matrix(M)` to the fitting functions instead.
#' @return A list with the phenotype `y`, an intercept design matrix `X`, the
#'   genotype matrix `M`, the relationship matrix `G` (when `form_G = TRUE`),
#'   the subpopulation of each individual `pop`, and the true `h2`.
#' @references
#' VanRaden, P. M. (2008) Efficient methods to compute genomic predictions.
#' Journal of Dairy Science 91, 4414-4423. \doi{10.3168/jds.2007-0980}
#'
#' Balding, D. J. & Nichols, R. A. (1995) A method for quantifying
#' differentiation between populations at multi-allelic loci and its
#' implications for investigating identity and paternity. Genetica 96, 3-12.
#' \doi{10.1007/bf01441146}
#' @seealso [grm_matrix()], [reml_sketch()]
#' @examples
#' set.seed(1)
#' dat <- sim_genomic(n = 200, p = 500, h2 = 0.4, pops = 3)
#' dim(dat$G)
#' table(dat$pop)
#' @export
sim_genomic <- function(n = 1000L, p = 2000L, h2 = 0.5, pops = 1L,
                        fst = 0.05, form_G = TRUE) {
    if (!is.numeric(h2) || length(h2) != 1L || h2 < 0 || h2 > 1) {
        stop("`h2` must be between 0 and 1.", call. = FALSE)
    }
    if (!is.numeric(fst) || length(fst) != 1L || fst <= 0 || fst >= 1) {
        stop("`fst` must be strictly between 0 and 1.", call. = FALSE)
    }
    n <- as.integer(n)
    p <- as.integer(p)
    pops <- as.integer(pops)
    if (is.na(pops) || pops < 1L || pops > n) {
        stop("`pops` must be between 1 and `n`.", call. = FALSE)
    }
    anc <- stats::runif(p, 0.05, 0.5)
    pop <- sort(rep_len(seq_len(pops), n))
    f <- if (pops == 1L) {
        matrix(anc, 1L, p)
    } else {
        a <- anc * (1 - fst) / fst
        b <- (1 - anc) * (1 - fst) / fst
        matrix(stats::rbeta(pops * p, rep(a, each = pops),
                            rep(b, each = pops)), pops, p)
    }
    M <- matrix(stats::rbinom(n * p, 2L, f[pop, , drop = FALSE]), n, p)
    freq <- colMeans(M) / 2
    Z <- sweep(M, 2L, 2 * freq)
    denom <- 2 * sum(freq * (1 - freq))
    u <- stats::rnorm(p, sd = sqrt(h2 / denom))
    y <- 10 + drop(Z %*% u) + stats::rnorm(n, sd = sqrt(1 - h2))
    list(y = y, X = matrix(1, n, 1L), M = M,
         G = if (form_G) tcrossprod(Z) / denom else NULL,
         pop = factor(pop), h2 = h2)
}
