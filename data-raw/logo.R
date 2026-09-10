## Draws man/figures/logo.png: a kernel matrix with the few columns that
## randomly pivoted Cholesky reads picked out in orange.
out <- file.path("man", "figures", "logo.png")
png(out, width = 518, height = 600, bg = "transparent", res = 144)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(-0.87, 0.87), ylim = c(-1, 1), asp = 1,
            xaxs = "i", yaxs = "i")
th <- pi / 2 + (0:5) * pi / 3
polygon(0.97 * cos(th), 0.97 * sin(th), col = "#1F2A44", border = "#D55E00",
        lwd = 7)

x <- c(-1.1, -0.8, -0.55, -0.35, -0.1, 0.15, 2.7, 2.9, 3.15)
K <- exp(-outer(x, x, "-")^2 / 0.9)
n <- length(x)
cell <- 0.1
gap <- 0.014
x0 <- -n * cell / 2
y0 <- 0.6
piv <- c(2, 6, 8)
for (i in seq_len(n)) for (j in seq_len(n)) {
    col <- if (j %in% piv) {
        grDevices::adjustcolor("#F28C28", alpha.f = 0.3 + 0.7 * K[i, j])
    } else {
        grDevices::adjustcolor("white", alpha.f = 0.07 + 0.5 * K[i, j])
    }
    rect(x0 + (j - 1) * cell, y0 - i * cell, x0 + j * cell - gap,
         y0 - (i - 1) * cell - gap, col = col, border = NA)
}
text(0, -0.56, "matsketch", col = "white", font = 2, cex = 1.75)
dev.off()
