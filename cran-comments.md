## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Test environments

* Local: Windows 11, R 4.6.1
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R devel,
  release and oldrel-1)

## Notes for the reviewer

* Words that may be flagged as misspelled in DESCRIPTION (Epperly,
  Frangella, Musco, Nystrom, Tropp, Udell, XDiag, XNysTrace, XTrace) are the
  names of authors and of the methods the package implements.
* The implementations follow the reference code published with the cited
  papers and are tested against brute-force versions of each estimator's
  definition.
* The timings shown in the "genomic-reml" vignette are read from
  `inst/extdata`, so building both vignettes takes about a minute.
