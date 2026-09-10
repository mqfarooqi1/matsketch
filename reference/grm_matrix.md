# Lazy genomic relationship matrix

Describes the genomic relationship matrix of VanRaden (2008), \$\$G =
\frac{Z Z^\top}{2 \sum_j f_j (1 - f_j)},\$\$ where \\Z\\ holds the
marker genotypes centred by twice the allele frequencies \\f_j\\,
without forming it. A product with \\G\\ costs two products with \\Z\\,
and rows are computed on demand, so the functions in this package can
work with \\G\\ while only the \\n \times p\\ genotypes are held in
memory.

## Usage

``` r
grm_matrix(M, freq = NULL)
```

## Arguments

- M:

  Genotype matrix with one row per individual and one column per marker,
  coded as allele counts between 0 and 2. Missing genotypes must be
  imputed first.

- freq:

  Allele frequencies used for centring. Defaults to the observed
  frequencies, `colMeans(M) / 2`.

## Value

An object of class `matsketch_grm`, accepted wherever `matsketch`
expects a matrix. [`as.matrix()`](https://rdrr.io/r/base/matrix.html)
forms the full matrix.

## Details

This matters once \\n\\ is large: for 50,000 individuals the full matrix
takes 20 GB, while the genotypes on a 10,000-marker panel take 4 GB.

## References

VanRaden, P. M. (2008) Efficient methods to compute genomic predictions.
Journal of Dairy Science 91, 4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

## See also

[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md),
[`sim_genomic()`](https://mqfarooqi1.github.io/matsketch/reference/sim_genomic.md)

## Examples

``` r
set.seed(1)
M <- matrix(rbinom(200 * 500, 2, 0.3), 200)
G <- grm_matrix(M)
G
#> <matsketch_grm> relationships among 200 individuals from 500 markers
#>   full matrix would hold 40,000 entries; only genotypes are stored
G$block(1:3, 1:3)
#>             [,1]        [,2]         [,3]
#> [1,]  0.98263803 0.090872547 -0.040163453
#> [2,]  0.09087255 1.003208604  0.008303878
#> [3,] -0.04016345 0.008303878  1.051365706
```
