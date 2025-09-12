#' Convert an edge list to an ape phylo object.
#'
#' Takes a two column matrix (parents -> descendants) and constructs an ape
#' phylo object from it.
#'
#' @param edgelist A two-column matrix where each row represents a
#' parent-descendant relationship.
#' @return An object of class 'phylo' representing the tree structure.
#' @export
edgelist_to_phylo <- function(edgelist) {
    # Identify tips and internal nodes
    tips <- setdiff(edgelist[,2], edgelist[,1])
    internals <- setdiff(unique(edgelist[,1]), tips)

    # Assign indices: tips 1:n, internals n+1:n+Nnode
    tip_ids <- setNames(seq_along(tips), tips)
    internal_ids <- setNames(length(tips) + seq_along(internals), internals)
    node_ids <- c(tip_ids, internal_ids)

    # construct the tree object
    edge <- cbind(node_ids[edgelist[, 1]], node_ids[edgelist[, 2]])
    phylo <- list(
        edge = edge,
        tip.label = tips,
        node.label = internals,
        Nnode = length(internals)
    )

    class(phylo) <- "phylo"
    phylo
}

#' Convert a taxonomic table to an ape phylo object.
#'
#' Creates a phylo from a taxonomic tree, skipping over any NA assignments.
#' Assumes that the columns are sorted from coarsest to finest taxonomic
#' resolution.
#'
#' @param taxa A data.frame or matrix with columns representing taxonomic
#' ranks (sorted coarsest to finest) and rows representing taxa. NA or empty
#' values are skipped.
#' @return An object of class 'phylo' representing the taxonomic tree.
#' @export
taxonomy_to_tree <- function(taxa) {
    edges <- list()
    ranks <- colnames(taxa)

    for (i in seq_len(nrow(taxa))) {
        path <- taxa[i, ranks]
        path <- as.character(path)
        path <- path[!is.na(path) & path != ""]
        if (length(path) > 1) {
            for (j in seq_len(length(path) - 1)) {
                edges[[length(edges) + 1]] <- c(path[j], path[j + 1])
            }
        }
    }

    # merge into a phylogeny
    do.call(rbind, edges) |>
        unique() |>
        edgelist_to_phylo()
}

#' Add rank prefixes to taxonomic matrix values.
#'
#' For each value in the taxonomic matrix (except the last column), adds a
#' prefix based on the first character of the column (rank) name. For example,
#' if the column is "Genus", the prefix will be "G_".
#'
#' @param taxa A character matrix with columns representing taxonomic ranks and
#' rows representing taxa.
#' @return A character matrix with prefixes added to each value (except the last
#' column).
#' @export
add_prefix <- function(taxa) {
    for (j in seq_len(ncol(taxa) - 1)) {
        prefix <- substr(colnames(taxa)[j], 1, 1)
        taxa[, j] <- paste0(prefix, "_", taxa[, j])
    }
    taxa
}