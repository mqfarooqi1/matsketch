## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* Local: Windows 11, R 4.6.1
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release and oldrel-1)

## Notes for the reviewer

* The implementations follow the reference code published with the cited
  papers and are tested against brute-force versions of each estimator's
  definition.
* The timings shown in the "genomic-reml" vignette are read from
  `inst/extdata`, so building the vignettes takes well under a minute.
