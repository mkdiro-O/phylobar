# Convert an edge list to an ape phylo object.

Takes a two column matrix (parents -\> descendants) and constructs an
ape phylo object from it.

## Usage

``` r
edgelist_to_phylo(edgelist)
```

## Arguments

- edgelist:

  A two-column matrix where each row represents a parent-descendant
  relationship.

## Value

An object of class 'phylo' representing the tree structure.

## Examples

``` r
# Example edge list: parent -> child
edgelist <- matrix(
  c("A", "B",
    "A", "C",
    "B", "D",
    "B", "E"),
  ncol = 2, byrow = TRUE,
  dimnames = list(NULL, c("parent", "child"))
)
tree <- edgelist_to_phylo(edgelist)
str(tree)
#> List of 4
#>  $ edge      : int [1:4, 1:2] 4 4 5 5 5 1 2 3
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "A" "A" "B" "B"
#>   .. ..$ : NULL
#>  $ tip.label : chr [1:3] "C" "D" "E"
#>  $ node.label: chr [1:2] "A" "B"
#>  $ Nnode     : int 2
#>  - attr(*, "class")= chr "phylo"
```
