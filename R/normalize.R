#' DESeq2 Size Factor Normalization
#'
#' This normalizes a count matrix using the same approach as
#' mStat_normalize_data from the MicrobiomeStat package. This can provide a
#' useful alternative to purely compositional transformations when constructing
#' the stacked bar plots. The only reason we don't use MicrobiomeStat directly
#' is that the current CRAN version does not have this function (only the
#' development GitHub version does).
#'
#' @param otu A numeric matrix with taxa as rows and samples as columns.
#'   Zero-sum rows are dropped before normalization.
#' @return A numeric matrix with taxa as rows and samples as columns containing
#'   the normalized counts.  Any NaN or Inf entries (which can arise when a
#'   sample has all-zero counts) are replaced with 0.
#' @examples
#' otu <- matrix(c(10L, 0L, 5L, 20L, 3L, 0L), nrow = 3,
#'               dimnames = list(c("t1", "t2", "t3"), c("s1", "s2")))
#' deseq_normalize(otu)
#' @export
deseq_normalize <- function(otu) {
    # Drop taxa with zero total counts (mirrors mStat_convert_phyloseq_to_data_obj)
    otu <- otu[rowSums(otu) > 0, , drop = FALSE]

    # Per-sample scale factor: sum(x) / geometric_mean(positive counts)
    scale_factors <- apply(otu, 2, function(x) {
        pos <- x[x > 0]
        if (length(pos) == 0L) return(NA_real_)
        sum(x) / exp(mean(log(pos)))
    })

    normalized <- sweep(otu, 2, scale_factors, "/")
    normalized[is.nan(normalized) | is.infinite(normalized)] <- 0
    normalized
}
