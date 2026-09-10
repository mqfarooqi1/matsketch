# Simulate a genomic data set

Draws biallelic marker genotypes, forms the genomic relationship matrix
of VanRaden (2008), and simulates a trait with the requested
heritability, so that
[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)
can be checked against known variance components.

## Usage

``` r
sim_genomic(
  n = 1000L,
  p = 2000L,
  h2 = 0.5,
  pops = 1L,
  fst = 0.05,
  form_G = TRUE
)
```

## Arguments

- n:

  Number of individuals.

- p:

  Number of markers.

- h2:

  Heritability, between 0 and 1.

- pops:

  Number of subpopulations.

- fst:

  Fixation index between subpopulations, used when `pops > 1`.

- form_G:

  Return the relationship matrix itself. Set to `FALSE` for large `n`
  and pass `grm_matrix(M)` to the fitting functions instead.

## Value

A list with the phenotype `y`, an intercept design matrix `X`, the
genotype matrix `M`, the relationship matrix `G` (when `form_G = TRUE`),
the subpopulation of each individual `pop`, and the true `h2`.

## Details

With `pops > 1` the individuals come from that many subpopulations whose
allele frequencies have drifted apart under the Balding-Nichols model
with fixation index `fst`. Population structure gives the relationship
matrix a few large eigenvalues, as real breeding and human cohorts do,
and those are exactly what a low-rank preconditioner captures.

Genetic values are sums of marker effects on the centred genotypes, with
variance chosen so that the genetic and residual variance components on
the scale of the relationship matrix are `h2` and `1 - h2`.

The function draws random numbers but does not set the seed; call
[`base::set.seed()`](https://rdrr.io/r/base/Random.html) first for a
reproducible data set.

## References

VanRaden, P. M. (2008) Efficient methods to compute genomic predictions.
Journal of Dairy Science 91, 4414-4423.
[doi:10.3168/jds.2007-0980](https://doi.org/10.3168/jds.2007-0980)

Balding, D. J. & Nichols, R. A. (1995) A method for quantifying
differentiation between populations at multi-allelic loci and its
implications for investigating identity and paternity. Genetica 96,
3-12. [doi:10.1007/bf01441146](https://doi.org/10.1007/bf01441146)

## See also

[`grm_matrix()`](https://mqfarooqi1.github.io/matsketch/reference/grm_matrix.md),
[`reml_sketch()`](https://mqfarooqi1.github.io/matsketch/reference/reml_sketch.md)

## Examples

``` r
set.seed(1)
dat <- sim_genomic(n = 200, p = 500, h2 = 0.4, pops = 3)
dim(dat$G)
#> [1] 200 200
table(dat$pop)
#> 
#>  1  2  3 
#> 67 67 66 
```
