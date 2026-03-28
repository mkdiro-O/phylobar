# Bray-Curtis dissimilarity matrix

Pairwise Bray-Curtis dissimilarity between matrix rows. \\\sum \|x_i -
y_i\| / (\sum x_i + \sum y_i)\\.

## Usage

``` r
bray_curtis_dist(x)
```

## Arguments

- x:

  A matrix with samples as rows.

## Value

A `dist` object of pairwise Bray-Curtis dissimilarities.
