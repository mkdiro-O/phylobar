# Convert a taxonomic table to an ape phylo object.

Creates a phylo from a taxonomic tree, skipping over any NA assignments.
Assumes that the columns are sorted from coarsest to finest taxonomic
resolution.

## Usage

``` r
taxonomy_to_tree(taxa)
```

## Arguments

- taxa:

  A data.frame or matrix with columns representing taxonomic ranks
  (sorted coarsest to finest) and rows representing taxa. NA or empty
  values are skipped.

## Value

An object of class 'phylo' representing the taxonomic tree.

## Examples

``` r
taxa <- matrix(
  c("Firmicutes", "Bacilli", "Lactobacillales",
    "Proteobacteria", "Gammaproteobacteria", "Enterobacterales"),
  ncol = 3, byrow = TRUE,
  dimnames = list(NULL, c("Phylum", "Class", "Order"))
)
tree <- taxonomy_to_tree(taxa)
str(tree)
#> List of 4
#>  $ edge      : int [1:4, 1:2] 3 4 5 6 4 1 6 2
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:4] "Firmicutes" "Bacilli" "Proteobacteria" "Gammaproteobacteria"
#>   .. ..$ : NULL
#>  $ tip.label : chr [1:2] "Lactobacillales" "Enterobacterales"
#>  $ node.label: chr [1:4] "Firmicutes" "Bacilli" "Proteobacteria" "Gammaproteobacteria"
#>  $ Nnode     : int 4
#>  - attr(*, "class")= chr "phylo"

# A more involved example with missing values
taxa <- matrix(
  c("Firmicutes", "Bacilli", "Lactobacillales", "ASV1",
    "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "ASV2",
    "Firmicutes", "Bacilli", NA, "ASV3"),
  ncol = 4, byrow = TRUE,
  dimnames = list(NULL, c("Phylum", "Class", "Order", "ASV"))
)
taxmat <- add_prefix(taxa)
tree <- taxonomy_to_tree(taxmat)
str(tree)
#> List of 4
#>  $ edge      : int [1:7, 1:2] 4 5 6 7 8 9 5 5 6 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:7] "P_Firmicutes" "C_Bacilli" "O_Lactobacillales" "P_Proteobacteria" ...
#>   .. ..$ : NULL
#>  $ tip.label : chr [1:3] "ASV1" "ASV2" "ASV3"
#>  $ node.label: chr [1:6] "P_Firmicutes" "C_Bacilli" "O_Lactobacillales" "P_Proteobacteria" ...
#>  $ Nnode     : int 6
#>  - attr(*, "class")= chr "phylo"
```
