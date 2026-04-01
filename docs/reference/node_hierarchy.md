# Construct Initial Hierarchy

This reshapes the list output from node_totals into the hierarchical
format needed for the d3 tree visualization.

## Usage

``` r
node_hierarchy(tree, totals, node = NULL)
```

## Arguments

- tree:

  An object of class phylo, representing the tree structure.

- totals:

  A named list of node totals, as returned by node_totals.

- node:

  Name of the node from which to start a recursion. Defaults to the root
  node.

## Value

A nested list representing the hierarchy, with each node containing '
its name, value, summary, and children (if any).

## Examples

``` r
library(ape)
tree <- rtree(5)
x_mat <- matrix(runif(15), ncol = 5)
colnames(x_mat) <- tree$tip.label

tree$node.label <- as.character(seq_len(tree$Nnode))
totals <- c(
  node_totals(tree, x_mat),
  as.list(data.frame(x_mat))
)
node_hierarchy(tree, totals)
#> $name
#> [1] "1"
#> 
#> $value
#>          [,1]
#> [1,] 1.289513
#> [2,] 2.609064
#> [3,] 3.373775
#> 
#> $summary
#> [1] 7.272352
#> 
#> $children
#> $children[[1]]
#> $children[[1]]$name
#> [1] "t5"
#> 
#> $children[[1]]$value
#>            [,1]
#> [1,] 0.06366146
#> [2,] 0.38870131
#> [3,] 0.97554784
#> 
#> $children[[1]]$summary
#> [1] 1.427911
#> 
#> $children[[1]]$children
#> [1] NA
#> 
#> 
#> $children[[2]]
#> $children[[2]]$name
#> [1] "2"
#> 
#> $children[[2]]$value
#>          [,1]
#> [1,] 1.225851
#> [2,] 2.220363
#> [3,] 2.398228
#> 
#> $children[[2]]$summary
#> [1] 5.844442
#> 
#> $children[[2]]$children
#> $children[[2]]$children[[1]]
#> $children[[2]]$children[[1]]$name
#> [1] "3"
#> 
#> $children[[2]]$children[[1]]$value
#>           [,1]
#> [1,] 0.5372953
#> [2,] 2.1891326
#> [3,] 2.1726650
#> 
#> $children[[2]]$children[[1]]$summary
#> [1] 4.899093
#> 
#> $children[[2]]$children[[1]]$children
#> $children[[2]]$children[[1]]$children[[1]]
#> $children[[2]]$children[[1]]$children[[1]]$name
#> [1] "4"
#> 
#> $children[[2]]$children[[1]]$children[[1]]$value
#>          [,1]
#> [1,] 0.485849
#> [2,] 1.658920
#> [3,] 1.476841
#> 
#> $children[[2]]$children[[1]]$children[[1]]$summary
#> [1] 3.62161
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children
#> $children[[2]]$children[[1]]$children[[1]]$children[[1]]
#> $children[[2]]$children[[1]]$children[[1]]$children[[1]]$name
#> [1] "t2"
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[1]]$value
#>           [,1]
#> [1,] 0.2898923
#> [2,] 0.6783804
#> [3,] 0.7353196
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[1]]$summary
#> [1] 1.703592
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[2]]
#> $children[[2]]$children[[1]]$children[[1]]$children[[2]]$name
#> [1] "t3"
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 0.1959567
#> [2,] 0.9805397
#> [3,] 0.7415215
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[2]]$summary
#> [1] 1.918018
#> 
#> $children[[2]]$children[[1]]$children[[1]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $children[[2]]$children[[1]]$children[[2]]
#> $children[[2]]$children[[1]]$children[[2]]$name
#> [1] "t1"
#> 
#> $children[[2]]$children[[1]]$children[[2]]$value
#>            [,1]
#> [1,] 0.05144628
#> [2,] 0.53021246
#> [3,] 0.69582388
#> 
#> $children[[2]]$children[[1]]$children[[2]]$summary
#> [1] 1.277483
#> 
#> $children[[2]]$children[[1]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> $children[[2]]$children[[2]]
#> $children[[2]]$children[[2]]$name
#> [1] "t4"
#> 
#> $children[[2]]$children[[2]]$value
#>            [,1]
#> [1,] 0.68855600
#> [2,] 0.03123033
#> [3,] 0.22556253
#> 
#> $children[[2]]$children[[2]]$summary
#> [1] 0.9453489
#> 
#> $children[[2]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> 
```
