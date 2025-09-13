
setup_tree_data <- function() {
    library(ape)
    set.seed(1)
    tree <- rtree(5)
    tree$node.label <- paste("node", seq_len(tree$Nnode))

    x_mat <- matrix(runif(15), ncol = 5)
    colnames(x_mat) <- tree$tip.label
    list(tree = tree, x_mat = x_mat)
}

test_that("node_totals output names match internal node labels", {
    td <- setup_tree_data()
    result <- node_totals(td$tree, td$x_mat)
    expect_equal(names(result), td$tree$node.label)
})

test_that("root node total equals sum of x_mat", {
    td <- setup_tree_data()
    result <- node_totals(td$tree, td$x_mat)
    root_label <- td$tree$node.label[1]
    root_total <- sum(result[[1]])
    expect_equal(root_total, sum(td$x_mat))
})

test_that("internal node with two tip children equals sum of tips", {
    # first use the function
    td <- setup_tree_data()
    tree <- td$tree
    x_mat <- td$x_mat
    result <- node_totals(tree, x_mat)

    # we'll compare with a manual calculation for two sibling tips.
    # we need to search over tips until we find such a sibling.
    for (i in seq_along(tree$node.label)) {
        node_num <- length(tree$tip.label) + i
        desc <- phangorn::Descendants(tree, node_num, type = "tips")[[1]]
        if (length(desc) == 2) {
            tip_labels <- tree$tip.label[desc]
            node_label <- tree$node.label[i]

            tip_totals <- rowSums(x_mat[, tip_labels])
            expect_equal(result[[node_label]], tip_totals)
            break
        }
    }
})