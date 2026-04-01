# Prepare tree data for phylobar visualization

Prepare tree data for phylobar visualization

## Usage

``` r
phylobar_data(x, tree, hclust_order = TRUE)
```

## Arguments

- x:

  A matrix of abundances. Samples along rows, features along columns.

- tree:

  A n object of class phylo, representing the tree structure.

- hclust_order:

  Logical; if TRUE, reorder rows/columns by hierarchical clustering.

## Value

A list with tree_data and labels.

## Examples

``` r
library(ape)
tree <- rtree(5)
tree$node.label <- paste0("node", seq_len(4))
x <- matrix(runif(15), nrow = 3)
colnames(x) <- tree$tip.label
rownames(x) <- paste0("sample", seq_len(3))
phylobar_data(x, tree)
#> $tree_data
#> $tree_data$name
#> [1] "node1"
#> 
#> $tree_data$value
#>          [,1]
#> [1,] 4.422935
#> [2,] 2.530391
#> [3,] 2.382729
#> 
#> $tree_data$summary
#> [1] 9.336055
#> 
#> $tree_data$children
#> $tree_data$children[[1]]
#> $tree_data$children[[1]]$name
#> [1] "node2"
#> 
#> $tree_data$children[[1]]$value
#>          [,1]
#> [1,] 2.645466
#> [2,] 1.096146
#> [3,] 1.881193
#> 
#> $tree_data$children[[1]]$summary
#> [1] 5.622805
#> 
#> $tree_data$children[[1]]$children
#> $tree_data$children[[1]]$children[[1]]
#> $tree_data$children[[1]]$children[[1]]$name
#> [1] "node3"
#> 
#> $tree_data$children[[1]]$children[[1]]$value
#>          [,1]
#> [1,] 1.870553
#> [2,] 0.511671
#> [3,] 1.247217
#> 
#> $tree_data$children[[1]]$children[[1]]$summary
#> [1] 3.62944
#> 
#> $tree_data$children[[1]]$children[[1]]$children
#> $tree_data$children[[1]]$children[[1]]$children[[1]]
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$name
#> [1] "t1"
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$value
#>            [,1]
#> [1,] 0.87420594
#> [2,] 0.01147954
#> [3,] 0.88824957
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$summary
#> [1] 1.773935
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$name
#> [1] "t2"
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 0.9963469
#> [2,] 0.5001915
#> [3,] 0.3589670
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$summary
#> [1] 1.855505
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $tree_data$children[[1]]$children[[2]]
#> $tree_data$children[[1]]$children[[2]]$name
#> [1] "t4"
#> 
#> $tree_data$children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 0.7749130
#> [2,] 0.5844753
#> [3,] 0.6339764
#> 
#> $tree_data$children[[1]]$children[[2]]$summary
#> [1] 1.993365
#> 
#> $tree_data$children[[1]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $tree_data$children[[2]]
#> $tree_data$children[[2]]$name
#> [1] "node4"
#> 
#> $tree_data$children[[2]]$value
#>           [,1]
#> [1,] 1.7774694
#> [2,] 1.4342445
#> [3,] 0.5015357
#> 
#> $tree_data$children[[2]]$summary
#> [1] 3.71325
#> 
#> $tree_data$children[[2]]$children
#> $tree_data$children[[2]]$children[[1]]
#> $tree_data$children[[2]]$children[[1]]$name
#> [1] "t3"
#> 
#> $tree_data$children[[2]]$children[[1]]$value
#>           [,1]
#> [1,] 0.8586662
#> [2,] 0.5668943
#> [3,] 0.2529970
#> 
#> $tree_data$children[[2]]$children[[1]]$summary
#> [1] 1.678558
#> 
#> $tree_data$children[[2]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $tree_data$children[[2]]$children[[2]]
#> $tree_data$children[[2]]$children[[2]]$name
#> [1] "t5"
#> 
#> $tree_data$children[[2]]$children[[2]]$value
#>           [,1]
#> [1,] 0.9188032
#> [2,] 0.8673502
#> [3,] 0.2485387
#> 
#> $tree_data$children[[2]]$children[[2]]$summary
#> [1] 2.034692
#> 
#> $tree_data$children[[2]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> 
#> 
#> $labels
#> [1] "sample3" "sample1" "sample2"
#> 
```
