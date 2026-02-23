library(testthat)
library(phylobar)

test_that("edgelist_to_phylo works with basic example", {
    # Example from documentation
    edgelist <- matrix(
        c("A", "B",
          "A", "C",
          "B", "D",
          "B", "E"),
        ncol = 2, byrow = TRUE,
        dimnames = list(NULL, c("parent", "child"))
    )

    tree <- edgelist_to_phylo(edgelist)

    # Test structure
    expect_s3_class(tree, "phylo")
    expect_named(tree, c("edge", "tip.label", "node.label", "Nnode"))

    # Test tip labels
    expect_equal(sort(tree$tip.label), c("C", "D", "E"))

    # Test node labels
    expect_equal(sort(tree$node.label), c("A", "B"))

    # Test number of internal nodes
    expect_equal(tree$Nnode, 2)

    # Test edge matrix dimensions
    expect_equal(nrow(tree$edge), 4)
    expect_equal(ncol(tree$edge), 2)
})

test_that("taxonomy_to_tree works with basic example", {
    # Basic example from documentation
    taxa <- matrix(
        c("Firmicutes", "Bacilli", "Lactobacillales",
          "Proteobacteria", "Gammaproteobacteria", "Enterobacterales"),
        ncol = 3, byrow = TRUE,
        dimnames = list(NULL, c("Phylum", "Class", "Order"))
    )

    tree <- taxonomy_to_tree(taxa)

    # Test structure
    expect_s3_class(tree, "phylo")
    expect_named(tree, c("edge", "tip.label", "node.label", "Nnode"))

    # Test that all taxonomic names are preserved
    all_names <- c("Firmicutes", "Bacilli", "Lactobacillales",
                   "Proteobacteria", "Gammaproteobacteria", "Enterobacterales")
    tree_names <- c(tree$tip.label, tree$node.label)
    expect_true(all(all_names %in% tree_names))
})

test_that("taxonomy_to_tree works with missing values example", {
    # More complex example with missing values
    taxa <- matrix(
        c("Firmicutes", "Bacilli", "Lactobacillales", "ASV1",
          "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "ASV2",
          "Firmicutes", "Bacilli", NA, "ASV3"),
        ncol = 4, byrow = TRUE,
        dimnames = list(NULL, c("Phylum", "Class", "Order", "ASV"))
    )

    taxmat <- add_prefix(taxa)
    tree <- taxonomy_to_tree(taxmat)

    # Test structure
    expect_s3_class(tree, "phylo")

    # Test that ASVs are in tip labels
    expect_true("ASV1" %in% tree$tip.label)
    expect_true("ASV2" %in% tree$tip.label)
    expect_true("ASV3" %in% tree$tip.label)

    # Test that prefixed taxonomic names are present
    expect_true("P_Firmicutes" %in% c(tree$tip.label, tree$node.label))
    expect_true("C_Bacilli" %in% c(tree$tip.label, tree$node.label))
})

test_that("add_prefix works correctly", {
    # Example from documentation
    taxa <- matrix(
        c("Firmicutes", "Bacilli", "Lactobacillales",
          "Proteobacteria", "Gammaproteobacteria", "Enterobacterales"),
        ncol = 3, byrow = TRUE,
        dimnames = list(NULL, c("Phylum", "Class", "Order"))
    )

    result <- add_prefix(taxa)

    # Test that prefixes are added correctly (using as.character to remove names)
    expect_equal(as.character(result[1, 1]), "P_Firmicutes")
    expect_equal(as.character(result[1, 2]), "C_Bacilli")
    expect_equal(as.character(result[2, 1]), "P_Proteobacteria")
    expect_equal(as.character(result[2, 2]), "C_Gammaproteobacteria")

    # Test that last column is unchanged
    expect_equal(as.character(result[1, 3]), "Lactobacillales")
    expect_equal(as.character(result[2, 3]), "Enterobacterales")

    # Test dimensions are preserved
    expect_equal(dim(result), dim(taxa))
    expect_equal(colnames(result), colnames(taxa))
})

test_that("add_prefix handles NA values correctly", {
    taxa <- matrix(
        c("Firmicutes", "Bacilli", NA,
          "Proteobacteria", NA, "Enterobacterales"),
        ncol = 3, byrow = TRUE,
        dimnames = list(NULL, c("Phylum", "Class", "Order"))
    )

    result <- add_prefix(taxa)

    # Test that prefixes are added to non-NA values (using as.character to remove names)
    expect_equal(as.character(result[1, 1]), "P_Firmicutes")
    expect_equal(as.character(result[1, 2]), "C_Bacilli")
    expect_equal(as.character(result[2, 1]), "P_Proteobacteria")

    # Test that NA values remain NA
    expect_true(is.na(result[1, 3]))
    expect_true(is.na(result[2, 2]))

    # Test that last column is unchanged
    expect_equal(as.character(result[2, 3]), "Enterobacterales")
})

test_that("edgelist_to_phylo handles edge cases", {
    # Test with minimal tree (single edge)
    edgelist <- matrix(c("root", "tip"), ncol = 2)
    tree <- edgelist_to_phylo(edgelist)

    expect_s3_class(tree, "phylo")
    expect_equal(tree$tip.label, "tip")
    expect_equal(tree$node.label, "root")
    expect_equal(tree$Nnode, 1)
})

test_that("taxonomy_to_tree handles empty or single column input", {
    # Test with single taxonomic rank
    taxa <- matrix(c("ASV1", "ASV2"), ncol = 1,
                   dimnames = list(NULL, "ASV"))

    # This should handle gracefully (though may produce empty tree)
    expect_no_error(taxonomy_to_tree(taxa))
})