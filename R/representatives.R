#' Subsample rows of a matrix using hierarchical clustering
#'
#' This function takes a matrix `x` and subsamples its rows by clustering
#' them using hierarchical clustering. One representative row is selected
#' from each cluster, which helps in visualizing a large dataset in a
#' semi-representative way.
#'
#' @param x A matrix whose rows will be clustered and subsampled.
#' @param k The number of clusters to form (default is 100).
#' @return A matrix containing one representative row from each cluster.
#' @export
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100)
#' result <- subset_cluster(mat, k = 10)
#' dim(result) # only 10 representatives
subset_cluster <- function(x, k = 100) {
    fit <- hclust(dist(x))
    clusters <- cutree(fit, k = k)

    # one representative sample per cluster
    representatives <- tapply(rownames(x), clusters, \(k) k[1])
    x[representatives, ]
}