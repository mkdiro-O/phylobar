# Add rank prefixes to taxonomic matrix values.

For each value in the taxonomic matrix (except the last column), adds a
prefix based on the first character of the column (rank) name. For
example, if the column is "Genus", the prefix will be "G\_".

## Usage

``` r
add_prefix(taxa)
```

## Arguments

- taxa:

  A character matrix with columns representing taxonomic ranks and rows
  representing taxa.

## Value

A character matrix with prefixes added to each value (except the last
column).

## Examples

``` r
taxa <- matrix(
  c("Firmicutes", "Bacilli", "Lactobacillales",
    "Proteobacteria", "Gammaproteobacteria", "Enterobacterales"),
  ncol = 3, byrow = TRUE,
  dimnames = list(NULL, c("Phylum", "Class", "Order"))
)
add_prefix(taxa)
#>      Phylum             Class                   Order             
#> [1,] "P_Firmicutes"     "C_Bacilli"             "Lactobacillales" 
#> [2,] "P_Proteobacteria" "C_Gammaproteobacteria" "Enterobacterales"
```
