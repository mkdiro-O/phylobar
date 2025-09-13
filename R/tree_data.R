#' Compute Node Totals
#'
#' This loops over all internal nodes in the tree and takes the sum over all
#' descendant taxa, for each sample.
#' @param tree An object of class phylo, representing the tree structure. Must
#'   have tip labels matching the columns of x_mat.
#' @param x_mat A numeric matrix of abundances, with samples in rows and
#'   features (tips) in columns. Column names should correspond to tree tip
#    labels.
#' @return A named list where each element corresponds to an internal node (by
#'   node label) and contains a vector of totals for each sample, computed by
#'   summing abundances over all descendant tips.
#' @importFrom phangorn Descendants
#' @importFrom ape Ntip
#' @examples
#' library(ape)
#' set.seed(1)
#' tree <- rtree(5)
#' x_mat <- matrix(runif(15), ncol = 5)
#' colnames(x_mat) <- tree$tip.label
#' node_totals(tree, x_mat)
#' @export
node_totals <- function(tree, x_mat) {
    node_ids <- tree$node.label
    totals <- list()
    for (node in node_ids) {
        desc_tips <- Descendants(tree, node, type = "tips")[[1]]
        tip_labels <- tree$tip.label[desc_tips]
        totals[[as.character(node)]] <- rowSums(
            x_mat[, tip_labels, drop = FALSE]
        )
    }
    totals
}

#' Construct Initial Hierarchy
#'
#' This reshapes the list output from node_totals into the hierarchical format
#' needed for the d3 tree visualization.
#' @param tree An object of class phylo, representing the tree structure.
#' @param totals A named list of node totals, as returned by node_totals.
#' @param node Name of the node from which to start a recursion. Defaults to the
#'   root node.
#' @return A nested list representing the hierarchy, with each node containing '
#'   its name, value, summary, and children (if any).
#' @importFrom phangorn Children
#' @importFrom ape Ntip
#' @importFrom purrr map
#' @examples
#' library(ape)
#' set.seed(1)
#' tree <- rtree(5)
#' x_mat <- matrix(runif(15), ncol = 5)
#' colnames(x_mat) <- tree$tip.label
#'
#' totals <- node_totals(tree, x_mat)
#' node_hierarchy(tree, totals)
#' @export
node_hierarchy <- function(tree, totals, node = NULL) {
    if (is.null(node)) node <- Ntip(tree) + 1
    if (node <= Ntip(tree)) {
        tip_name <- tree$tip.label[node]
        list(
            name = tip_name,
            value = matrix(totals[[tip_name]]),
            summary = sum(totals[[tip_name]]),
            children = NA
        )
    } else {
        node_idx <- node - Ntip(tree)
        node_name <- tree$node.label[node_idx]
        children_nodes <- Children(tree, node)
        names(children_nodes) <- NULL

        children <- map(children_nodes, ~ node_hierarchy(tree, totals, .x))
        list(
            name = node_name,
            value = matrix(totals[[node_name]]),
            summary = sum(totals[[node_name]]),
            children = children
        )
    }
}

#' Prepare tree data for phylobar visualization
#'
#' @param tree A n object of class phylo, representing the tree structure.
#' @param x A matrix of abundances. Samples along rows, features along columns.
#' @param hclust_order Logical; if TRUE, reorder rows/columns by hierarchical
#'      clustering.
#' @return A list with tree_data and labels.
#' @examples
#' library(ape)
#' set.seed(1)
#' tree <- rtree(5)
#' x <- matrix(runif(15), nrow = 3)
#' colnames(x) <- tree$tip.label
#' rownames(x) <- paste0("sample", 1:3)
#' phylobar_data(x, tree)
#' @export
phylobar_data <- function(x, tree, hclust_order = TRUE) {
    if (tree$node.label[1] == "") {
        tree$node.label[1] <- "root"
    }

    totals <- node_totals(tree, x)
    for (j in seq_len(ncol(x))) {
        tax_id <- colnames(x)[j]
        totals[[tax_id]] <- x[, j]
    }
    if (hclust_order) {
        x <- x[hclust(dist(x))$order, ]
        x <- x[, hclust(dist(t(x)))$order]
    }

    tree_data <- node_hierarchy(tree, totals)
    labels <- rownames(x)
    list(tree_data = tree_data, labels = labels)
}