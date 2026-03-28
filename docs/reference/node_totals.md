# Compute Node Totals

This loops over all internal nodes in the tree and takes the sum over
all descendant taxa, for each sample.

## Usage

``` r
node_totals(tree, x_mat)
```

## Arguments

- tree:

  An object of class phylo, representing the tree structure. Must have
  tip labels matching the columns of x_mat.

- x_mat:

  A numeric matrix of abundances, with samples in rows and features
  (tips) in columns. Column names should correspond to tree tip

## Value

A named list where each element corresponds to an internal node (by node
label) and contains a vector of totals for each sample, computed by
summing abundances over all descendant tips.

## Examples

``` r
library(ape)
tree <- rtree(5)
x_mat <- matrix(runif(15), ncol = 5)
colnames(x_mat) <- tree$tip.label
tree$node.label <- as.character(seq_len(tree$Nnode))
node_totals(tree, x_mat)
#> $`1`
#> [1] 3.315717 2.257086 1.864552
#> 
#> $`2`
#> [1] 2.032234 1.118343 1.128430
#> 
#> $`3`
#> [1] 1.2070343 0.8445247 0.5583851
#> 
#> $`4`
#> [1] 1.2834830 1.1387432 0.7361214
#> 
```
