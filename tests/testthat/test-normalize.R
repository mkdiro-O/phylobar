library(testthat)
library(phylobar)

test_that("deseq_normalize returns taxa x samples matrix with correct dimnames", {
    otu <- matrix(
        c(10L, 0L, 5L, 20L, 3L, 0L), nrow = 3,
        dimnames = list(c("t1", "t2", "t3"), c("s1", "s2"))
    )
    result <- deseq_normalize(otu)

    expect_true(is.matrix(result))
    expect_equal(rownames(result), c("t1", "t2", "t3"))
    expect_equal(colnames(result), c("s1", "s2"))
})

test_that("deseq_normalize scale factors match the manual formula", {
    otu <- matrix(
        c(10L, 2L, 0L,  4L, 5L,  6L), nrow = 3, byrow = TRUE,
        dimnames = list(c("t1", "t2", "t3"), c("s1", "s2"))
    )
    result <- deseq_normalize(otu)

    # s1: positive counts are 10, 5; geometric mean = exp(mean(log(c(10, 5))))
    sf_s1 <- sum(c(10, 0, 5)) / exp(mean(log(c(10, 5))))
    expect_equal(
        result[, "s1"],
        c(t1 = 10 / sf_s1, t2 = 0 / sf_s1, t3 = 5 / sf_s1)
    )

    # s2: all counts are positive
    sf_s2 <- sum(c(2, 4, 6)) / exp(mean(log(c(2, 4, 6))))
    expect_equal(
        result[, "s2"],
        c(t1 = 2 / sf_s2, t2 = 4 / sf_s2, t3 = 6 / sf_s2)
    )
})
