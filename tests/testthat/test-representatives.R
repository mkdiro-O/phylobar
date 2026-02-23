library(testthat)
library(phylobar)

test_that("subset_cluster works with basic example", {
    # Example from documentation
    set.seed(123) # For reproducible results
    mat <- matrix(rnorm(1000), nrow = 100)
    rownames(mat) <- seq_len(100)
    result <- subset_cluster(mat, k = 10)

    # Test dimensions - should have 10 rows (representatives)
    expect_equal(nrow(result), 10)
    expect_equal(ncol(result), ncol(mat))

    # Test that result is a matrix
    expect_true(is.matrix(result))

    # Test that all row names are from original matrix
    expect_true(all(rownames(result) %in% rownames(mat)))

    # Test that row names are unique (no duplicates)
    expect_equal(length(rownames(result)), length(unique(rownames(result))))
})

test_that("subset_cluster handles different k values", {
    set.seed(123)
    mat <- matrix(rnorm(500), nrow = 50)
    rownames(mat) <- paste0("row_", seq_len(50))

    # Test with k = 5
    result_5 <- subset_cluster(mat, k = 5)
    expect_equal(nrow(result_5), 5)

    # Test with k equal to number of rows (each row is its own cluster)
    result_50 <- subset_cluster(mat, k = 50)
    expect_equal(nrow(result_50), 50)
})

test_that("subset_cluster preserves column structure", {
    set.seed(123)
    mat <- matrix(rnorm(200), nrow = 20)
    colnames(mat) <- paste0("col_", seq_len(ncol(mat)))
    rownames(mat) <- paste0("row_", seq_len(20))

    result <- subset_cluster(mat, k = 5)

    # Test that column names are preserved
    expect_equal(colnames(result), colnames(mat))

    # Test that all values in result are from original matrix
    expect_true(all(result %in% mat))
})

test_that("subset_cluster uses default k value", {
    set.seed(123)
    # Create matrix with more than 100 rows to test default k = 100
    mat <- matrix(rnorm(1500), nrow = 150)
    rownames(mat) <- seq_len(150)

    # Test with default k (should be 100)
    result_default <- subset_cluster(mat)
    expect_equal(nrow(result_default), 100)

    # Test that explicit k = 100 gives same result
    result_explicit <- subset_cluster(mat, k = 100)
    expect_equal(nrow(result_explicit), 100)
})

test_that("subset_cluster handles edge cases", {
    set.seed(123)

    # Test with very small matrix
    small_mat <- matrix(rnorm(6), nrow = 3)
    rownames(small_mat) <- c("a", "b", "c")
    result_small <- subset_cluster(small_mat, k = 2)
    expect_equal(nrow(result_small), 2)
})

test_that("subset_cluster clustering is deterministic with same input", {
    set.seed(123)
    mat <- matrix(rnorm(200), nrow = 20)
    rownames(mat) <- paste0("row_", seq_len(20))

    # Run twice with same seed
    set.seed(456)
    result1 <- subset_cluster(mat, k = 5)

    set.seed(456)
    result2 <- subset_cluster(mat, k = 5)

    # Results should be identical
    expect_identical(result1, result2)
})

test_that("subset_cluster selects one representative per cluster", {
    set.seed(123)
    mat <- matrix(rnorm(300), nrow = 30)
    rownames(mat) <- paste0("row_", seq_len(30))

    k_val <- 6
    result <- subset_cluster(mat, k = k_val)

    # Should have exactly k representatives
    expect_equal(nrow(result), k_val)

    # Each representative should be from original matrix
    for (i in seq_len(nrow(result))) {
        row_name <- rownames(result)[i]
        original_row <- mat[row_name, ]
        result_row <- result[i, ]
        expect_equal(as.numeric(result_row), as.numeric(original_row))
    }
})