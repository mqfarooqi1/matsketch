# Lazy kernel matrix

Describes the kernel matrix \\K\\ with entries \\k(x_i, x_j)\\ without
computing it. Randomly pivoted Cholesky then evaluates only the entries
it needs, roughly \\(k + 1) n\\ of the \\n^2\\, which is what makes it
practical when the full kernel matrix would not fit in memory.

## Usage

``` r
kernel_matrix(
  X,
  kernel = c("gaussian", "laplace", "matern32", "matern52"),
  bandwidth = NULL,
  block = 1000L
)
```

## Arguments

- X:

  Numeric matrix with one row per point.

- kernel:

  Kernel family: `"gaussian"` \\\exp(-r^2 / 2h^2)\\, `"laplace"`
  \\\exp(-r / h)\\, or the Matern kernels `"matern32"` and `"matern52"`,
  where \\r\\ is the Euclidean distance.

- bandwidth:

  Length scale \\h\\. Defaults to the median distance between points,
  computed on at most the first 500 rows.

- block:

  Number of rows computed at a time when multiplying by the whole
  matrix, which bounds memory use.

## Value

An object of class `matsketch_kernel`, accepted wherever `matsketch`
expects a matrix. [`as.matrix()`](https://rdrr.io/r/base/matrix.html)
forms the full matrix.

## Details

Each product with the whole matrix recomputes the kernel, block by
block. For methods that multiply many times, such as
[`pcg()`](https://mqfarooqi1.github.io/matsketch/reference/pcg.md), it
is faster to form the matrix once with
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) whenever it fits in
memory.

## See also

[`rpchol()`](https://mqfarooqi1.github.io/matsketch/reference/rpchol.md),
[`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md)

## Examples

``` r
X <- matrix(rnorm(400), ncol = 2)
K <- kernel_matrix(X, "gaussian")
K
#> <matsketch_kernel> gaussian kernel on 200 points, bandwidth 1.632
#>   full matrix would hold 40,000 entries; nothing is stored
K$block(1:3, 1:3)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.7960616 0.8923518
#> [2,] 0.7960616 1.0000000 0.5203043
#> [3,] 0.8923518 0.5203043 1.0000000
```
