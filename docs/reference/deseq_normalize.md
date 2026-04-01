# DESeq2 Size Factor Normalization

This normalizes a count matrix using the same approach as
mStat_normalize_data from the MicrobiomeStat package. This can provide a
useful alternative to purely compositional transformations when
constructing the stacked bar plots. The only reason we don't use
MicrobiomeStat directly is that the current CRAN version does not have
this function (only the development GitHub version does).

## Usage

``` r
deseq_normalize(otu)
```

## Arguments

- otu:

  A numeric matrix with taxa as rows and samples as columns. Zero-sum
  rows are dropped before normalization.

## Value

A numeric matrix with taxa as rows and samples as columns containing the
normalized counts. Any NaN or Inf entries (which can arise when a sample
has all-zero counts) are replaced with 0.

## Examples

``` r
otu <- matrix(c(10L, 0L, 5L, 20L, 3L, 0L), nrow = 3,
              dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
deseq_normalize(otu)
#>          s1       s2
#> t1 4.714045 6.735623
#> t2 0.000000 1.010343
#> t3 2.357023 0.000000
```
