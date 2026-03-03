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
#> [1,] 3.294990
#> [2,] 1.836605
#> [3,] 2.950020
#> 
#> $tree_data$summary
#> [1] 8.081616
#> 
#> $tree_data$children
#> $tree_data$children[[1]]
#> $tree_data$children[[1]]$name
#> [1] "node2"
#> 
#> $tree_data$children[[1]]$value
#>          [,1]
#> [1,] 2.368444
#> [2,] 1.284228
#> [3,] 2.372954
#> 
#> $tree_data$children[[1]]$summary
#> [1] 6.025626
#> 
#> $tree_data$children[[1]]$children
#> $tree_data$children[[1]]$children[[1]]
#> $tree_data$children[[1]]$children[[1]]$name
#> [1] "node3"
#> 
#> $tree_data$children[[1]]$children[[1]]$value
#>           [,1]
#> [1,] 1.4342445
#> [2,] 0.5015357
#> [3,] 1.3216844
#> 
#> $tree_data$children[[1]]$children[[1]]$summary
#> [1] 3.257465
#> 
#> $tree_data$children[[1]]$children[[1]]$children
#> $tree_data$children[[1]]$children[[1]]$children[[1]]
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$name
#> [1] "t5"
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$value
#>           [,1]
#> [1,] 0.5668943
#> [2,] 0.2529970
#> [3,] 0.9188032
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$summary
#> [1] 1.738695
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$name
#> [1] "t3"
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 0.8673502
#> [2,] 0.2485387
#> [3,] 0.4028812
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$summary
#> [1] 1.51877
#> 
#> $tree_data$children[[1]]$children[[1]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $tree_data$children[[1]]$children[[2]]
#> $tree_data$children[[1]]$children[[2]]$name
#> [1] "node4"
#> 
#> $tree_data$children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 0.9341994
#> [2,] 0.7826920
#> [3,] 1.0512700
#> 
#> $tree_data$children[[1]]$children[[2]]$summary
#> [1] 2.768161
#> 
#> $tree_data$children[[1]]$children[[2]]$children
#> $tree_data$children[[1]]$children[[2]]$children[[1]]
#> $tree_data$children[[1]]$children[[2]]$children[[1]]$name
#> [1] "t4"
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[1]]$value
#>           [,1]
#> [1,] 0.7696302
#> [2,] 0.1194854
#> [3,] 0.1946950
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[1]]$summary
#> [1] 1.083811
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[2]]
#> $tree_data$children[[1]]$children[[2]]$children[[2]]$name
#> [1] "t2"
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[2]]$value
#>           [,1]
#> [1,] 0.1645692
#> [2,] 0.6632066
#> [3,] 0.8565750
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[2]]$summary
#> [1] 1.684351
#> 
#> $tree_data$children[[1]]$children[[2]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> 
#> 
#> $tree_data$children[[2]]
#> $tree_data$children[[2]]$name
#> [1] "t1"
#> 
#> $tree_data$children[[2]]$value
#>           [,1]
#> [1,] 0.9265464
#> [2,] 0.5523776
#> [3,] 0.5770657
#> 
#> $tree_data$children[[2]]$summary
#> [1] 2.05599
#> 
#> $tree_data$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $labels
#> [1] "sample1" "sample2" "sample3"
#> 
```
