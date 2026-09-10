## The timing and accuracy study behind vignette("genomic-reml").
##
## Run from the package root. With R's reference BLAS on one core it takes
## one to two hours, most of it in the exact fits, and it writes two files to
## inst/extdata. The scaling file is rewritten after every fit, so a partial
## run still leaves usable results.

library(matsketch)

out <- file.path("inst", "extdata")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

timed <- function(expr) {
    gc()
    t <- system.time(value <- expr)[["elapsed"]]
    list(seconds = t, value = value)
}

## Accuracy: the same data sets fitted exactly and by the sketch, with the
## package defaults (rank 100, 40 products per trace estimate).
set.seed(20260911)
accuracy <- do.call(rbind, lapply(seq_len(40), function(r) {
    h2 <- c(0.3, 0.6)[1 + r %% 2]
    dat <- sim_genomic(n = 2000, p = 4000, h2 = h2, pops = 4, fst = 0.05)
    ex <- reml_exact(dat$y, dat$G)
    sk <- reml_sketch(dat$y, dat$G)
    message(sprintf("accuracy %2d: exact %.4f, sketch %.4f", r, ex$h2, sk$h2))
    data.frame(replicate = r, true_h2 = h2,
               exact_h2 = ex$h2, exact_se = ex$se[["h2"]],
               sketch_h2 = sk$h2, sketch_se = sk$se[["h2"]],
               sketch_trace_se = sk$history$trace_se[nrow(sk$history)])
}))
write.csv(accuracy, file.path(out, "reml-accuracy.csv"), row.names = FALSE)

## Scaling: time and memory as the number of individuals grows. Memory is the
## size of the n x n (or n x p) matrices each method holds, not a measurement.
set.seed(20260912)
p <- 5000
rows <- list()
gb <- function(doubles) doubles * 8 / 2^30
record <- function(n, method, res, memory_gb) {
    fit <- res$value
    is_fit <- inherits(fit, "reml_sketch")
    rows[[length(rows) + 1L]] <<- data.frame(
        n = n, method = method, seconds = res$seconds, memory_gb = memory_gb,
        h2 = if (is_fit) fit$h2 else NA_real_,
        iterations = if (is_fit) fit$iterations else NA_integer_,
        solves = if (is_fit) fit$solves else NA_real_,
        cg_iterations = if (is_fit) fit$cg_iterations else NA_real_)
    message(sprintf("n = %5d  %-12s %8.1f s", n, method, res$seconds))
    write.csv(do.call(rbind, rows), file.path(out, "reml-scaling.csv"),
              row.names = FALSE)
}

for (n in c(1000, 2000, 4000, 8000, 16000)) {
    dat <- sim_genomic(n = n, p = p, h2 = 0.5, pops = 4, fst = 0.05,
                       form_G = FALSE)
    Gl <- grm_matrix(dat$M)
    record(n, "sketch_lazy", timed(reml_sketch(dat$y, Gl)), gb(n * p))
    if (n <= 8000) {
        fg <- timed(as.matrix(Gl))
        G <- fg$value
        record(n, "form_G", fg, gb(n^2))
        record(n, "sketch_dense", timed(reml_sketch(dat$y, G)), gb(n^2))
        record(n, "exact", timed(reml_exact(dat$y, G)), gb(4 * n^2))
        if (n <= 4000) {
            record(n, "eigen", timed(eigen(G, symmetric = TRUE)), gb(2 * n^2))
        }
        rm(G, fg)
    }
    rm(dat, Gl)
}
writeLines(capture.output(sessionInfo()), file.path("data-raw",
                                                     "reml-benchmark-session.txt"))
