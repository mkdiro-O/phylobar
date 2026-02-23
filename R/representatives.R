#' Subsample rows of a matrix using clustering
#'
#' This function takes a matrix `x` and subsamples its rows by clustering
#' them. One representative row is selected from each cluster, which helps in
#' visualizing a large dataset in a semi-representative way.
#'
#' @param x A matrix whose rows will be clustered and subsampled.
#' @param k The number of clusters to form (default is 100).
#' @param method Clustering method to use. `"hclust"` (the default) uses
#'   Euclidean distances and hierarchical clustering, selecting an arbitrary
#'   member from each cluster. `"medoid"` uses Bray-Curtis dissimilarity and
#'   K-medoids (PAM) clustering, selecting the medoid of each cluster.
#' @return A matrix containing one representative row from each cluster.
#' @importFrom cluster pam
#' @importFrom stats dist hclust cutree as.dist
#' @export
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100)
#' rownames(mat) <- seq_len(100)
#' result <- subset_cluster(mat, k = 10)
#' dim(result) # only 10 representatives
#'
#' # Using Bray-Curtis + K-medoids (better for compositional data)
#' counts <- matrix(rpois(1000, 5), nrow = 100)
#' rownames(counts) <- seq_len(100)
#' result2 <- subset_cluster(counts, k = 10, method = "medoid")
subset_cluster <- function(x, k = 100, method = c("hclust", "medoid")) {
    method <- match.arg(method)

    if (method == "hclust") {
        fit <- hclust(dist(x))
        clusters <- cutree(fit, k = k)
        representatives <- tapply(rownames(x), clusters, \(k) k[1])
        return(x[representatives, ])
    }

    # Bray-Curtis dissimilarity + K-medoids
    d <- bray_curtis_dist(x)
    fit <- pam(d, k = k, diss = TRUE)
    x[fit$medoids, ]
}

#' Bray-Curtis dissimilarity matrix
#'
#' Pairwise Bray-Curtis dissimilarity between matrix rows.
#' \eqn{\sum |x_i - y_i| / (\sum x_i + \sum y_i)}.
#'
#' @param x A matrix with samples as rows.
#' @return A `dist` object of pairwise Bray-Curtis dissimilarities.
bray_curtis_dist <- function(x) {
    l1 <- as.matrix(dist(x, method = "manhattan"))
    rs <- rowSums(x)
    denom <- outer(rs, rs, "+")
    denom[denom == 0] <- 1
    as.dist(l1 / denom)
}
