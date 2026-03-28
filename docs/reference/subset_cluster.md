# Subsample rows of a matrix using clustering

This function takes a matrix `x` and subsamples its rows by clustering
them. One representative row is selected from each cluster, which helps
in visualizing a large dataset in a semi-representative way.

## Usage

``` r
subset_cluster(x, k = 100, method = c("hclust", "medoid"))
```

## Arguments

- x:

  A matrix whose rows will be clustered and subsampled.

- k:

  The number of clusters to form (default is 100).

- method:

  Clustering method to use. `"hclust"` (the default) uses Euclidean
  distances and hierarchical clustering, selecting an arbitrary member
  from each cluster. `"medoid"` uses Bray-Curtis dissimilarity and
  K-medoids (PAM) clustering, selecting the medoid of each cluster.

## Value

A matrix containing one representative row from each cluster.

## Examples

``` r
mat <- matrix(rnorm(1000), nrow = 100)
rownames(mat) <- seq_len(100)
result <- subset_cluster(mat, k = 10)
dim(result) # only 10 representatives
#> [1] 10 10

# Using Bray-Curtis + K-medoids (better for compositional data)
counts <- matrix(rpois(1000, 5), nrow = 100)
rownames(counts) <- seq_len(100)
result2 <- subset_cluster(counts, k = 10, method = "medoid")
```
