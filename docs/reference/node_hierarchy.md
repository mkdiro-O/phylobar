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
#> [1,] 2.877252
#> [2,] 1.668961
#> [3,] 2.884882
#> 
#> $summary
#> [1] 7.431095
#> 
#> $children
#> $children[[1]]
#> $children[[1]]$name
#> [1] "2"
#> 
#> $children[[1]]$value
#>          [,1]
#> [1,] 2.172665
#> [2,] 0.935959
#> [3,] 1.541982
#> 
#> $children[[1]]$summary
#> [1] 4.650606
#> 
#> $children[[1]]$children
#> $children[[1]]$children[[1]]
#> $children[[1]]$children[[1]]$name
#> [1] "t4"
#> 
#> $children[[1]]$children[[1]]$value
#>           [,1]
#> [1,] 0.7353196
#> [2,] 0.1959567
#> [3,] 0.9805397
#> 
#> $children[[1]]$children[[1]]$summary
#> [1] 1.911816
#> 
#> $children[[1]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $children[[1]]$children[[2]]
#> $children[[1]]$children[[2]]$name
#> [1] "3"
#> 
#> $children[[1]]$children[[2]]$value
#>           [,1]
#> [1,] 1.4373454
#> [2,] 0.7400023
#> [3,] 0.5614428
#> 
#> $children[[1]]$children[[2]]$summary
#> [1] 2.73879
#> 
#> $children[[1]]$children[[2]]$children
#> $children[[1]]$children[[2]]$children[[1]]
#> $children[[1]]$children[[2]]$children[[1]]$name
#> [1] "t2"
#> 
#> $children[[1]]$children[[2]]$children[[1]]$value
#>            [,1]
#> [1,] 0.74152153
#> [2,] 0.05144628
#> [3,] 0.53021246
#> 
#> $children[[1]]$children[[2]]$children[[1]]$summary
#> [1] 1.32318
#> 
#> $children[[1]]$children[[2]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $children[[1]]$children[[2]]$children[[2]]
#> $children[[1]]$children[[2]]$children[[2]]$name
#> [1] "t1"
#> 
#> $children[[1]]$children[[2]]$children[[2]]$value
#>            [,1]
#> [1,] 0.69582388
#> [2,] 0.68855600
#> [3,] 0.03123033
#> 
#> $children[[1]]$children[[2]]$children[[2]]$summary
#> [1] 1.41561
#> 
#> $children[[1]]$children[[2]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> 
#> 
#> $children[[2]]
#> $children[[2]]$name
#> [1] "4"
#> 
#> $children[[2]]$value
#>           [,1]
#> [1,] 0.7045871
#> [2,] 0.7330021
#> [3,] 1.3428995
#> 
#> $children[[2]]$summary
#> [1] 2.780489
#> 
#> $children[[2]]$children
#> $children[[2]]$children[[1]]
#> $children[[2]]$children[[1]]$name
#> [1] "t3"
#> 
#> $children[[2]]$children[[1]]$value
#>           [,1]
#> [1,] 0.2255625
#> [2,] 0.3008308
#> [3,] 0.6364656
#> 
#> $children[[2]]$children[[1]]$summary
#> [1] 1.162859
#> 
#> $children[[2]]$children[[1]]$children
#> [1] NA
#> 
#> 
#> $children[[2]]$children[[2]]
#> $children[[2]]$children[[2]]$name
#> [1] "t5"
#> 
#> $children[[2]]$children[[2]]$value
#>           [,1]
#> [1,] 0.4790245
#> [2,] 0.4321713
#> [3,] 0.7064338
#> 
#> $children[[2]]$children[[2]]$summary
#> [1] 1.61763
#> 
#> $children[[2]]$children[[2]]$children
#> [1] NA
#> 
#> 
#> 
#> 
#> 
```
