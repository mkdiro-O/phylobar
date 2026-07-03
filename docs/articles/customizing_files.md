# Customizing Style in phylobar

## Overview

Phylobar draws a phylogenetic tree alongside a stacked bar chart. This
vignette shows how to customize visual style: text sizes, color
palettes, layout ratios, legend placement, and more.

We will use a small random dataset:

``` r
library(ape)
library(phylobar)

set.seed(1)
tree <- rtree(20)
samples <- matrix(rpois(100 * 20, 1), nrow = 100, ncol = 20)

phylobar(samples, tree)
```

The function signature and key styling parameters we will explore:

- `palette`: colors used for painted subtrees (stacked bars)
- `width`, `height`: widget size in pixels
- `sample_font_size`, `sample_label_margin`, `sample_label_space`
- `sample_magnify`, `sample_show_all`
- `rel_width`, `rel_height`, `rel_space`
- `legend_mode`, `legend_x_start`, `legend_spacing`
- `hclust_order`: optional hierarchical reordering of rows/columns

Internally,
[`phylobar()`](http://mkdiro-o.github.io/phylobar/reference/phylobar.md)
ensures that the tree has node labels and that the abundance matrix has
row/column names. If they’re missing, helper `check_inputs()` creates
sensible defaults (e.g., `sample_1`, `sample_2`).

## Color Palette

If a palette is not provided, phylobar uses a default set of six colors.
It is possible to supply own vector of hex colors or R color names.

``` r
# Custom qualitative palette (as many colors as you want; repeats if needed)
my_palette <- c(
    "#4E79A7", "#F28E2B", "#E15759",
    "#76B7B2", "#59A14F", "#EDC948"
)
phylobar(samples, tree, palette = my_palette)
```

## Widget Size

By default the widget adapts to the container. Fix its size by:

``` r
phylobar(samples, tree,
    width = 800,
    height = 500
)
```

## Sample Label Styling

Control sample label font size, spacing, and hover magnification.

``` r
phylobar(samples, tree,
    width = 800, height = 500,
    sample_font_size = 10,
    sample_label_margin = 10, # space between labels and bars
    sample_label_space = 100, # reserved margin for labels
    sample_magnify = 1.3, # how much to enlarge labels on hover
    sample_show_all = TRUE
)
```

## Tree-bar Layout Ratio

Change how much horizontal and vertical space the tree occupies:

- `rel_width`: fraction of total width reserved for the tree panel
  (default 0.4)
- `rel_height`: fraction of total height for the tree panel (default
  0.85)
- `rel_space`: pixels between the tree and bar panels (default 10)

``` r
phylobar(
    samples, tree,
    width = 800, height = 500,
    sample_label_space = 100,
    sample_magnify = 1.3,
    rel_width = 0.2, # narrower tree
    rel_height = 0.70, # shorter tree (more space for legend)
    rel_space = 14 # larger gap between panels
)
```

## Legend Placement

Choose whether painted subtree labels appear in a separate legend
(`legend_mode = TRUE`, default) or placed inside the tree
(`legend_mode = FALSE`). `legend_x_start` and `legend_spacing` can also
be adjusted to control the legen position.

``` r
# Legend below the tree (default)
phylobar(
    samples, tree,
    width = 800, height = 500,
    sample_label_space = 100,
    sample_magnify = 1.3,
    legend_mode = TRUE,
    legend_x_start = 20, # horizontal start in pixels
    legend_spacing = 20 # vertical spacing between legend items in pixels
)
```

``` r
# Put labels inside the tree instead of a separate legend
phylobar(
    samples, tree,
    width = 800, height = 500,
    sample_label_space = 100,
    sample_magnify = 1.3,
    legend_mode = FALSE,
    legend_x_start = 20, # horizontal start in pixels
    legend_spacing = 20 # vertical spacing between legend items in pixels
)
```

## Optional Hierarchical Reordering

`hclust_order = TRUE` (default) reorders samples/features by
hierarchical clustering, which can make patterns more visible. Turn it
off to preserve original order:

``` r
# Keep the input order as-is
phylobar(
    samples, tree,
    width = 800, height = 500,
    sample_label_space = 100,
    sample_magnify = 1.3,
    hclust_order = FALSE
)
```

## Session Info

``` r
sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] phylobar_0.99.12 ape_5.8-1        BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5        jsonlite_2.0.0      compiler_4.6.1     
#>  [4] BiocManager_1.30.27 Rcpp_1.1.1-1.1      parallel_4.6.1     
#>  [7] cluster_2.1.8.2     jquerylib_0.1.4     systemfonts_1.3.2  
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] lattice_0.22-9      R6_2.6.1            generics_0.1.4     
#> [16] igraph_2.3.3        knitr_1.51          htmlwidgets_1.6.4  
#> [19] bookdown_0.47       desc_1.4.3          bslib_0.11.0       
#> [22] rlang_1.2.0         fastmatch_1.1-8     cachem_1.1.0       
#> [25] xfun_0.59           quadprog_1.5-8      fs_2.1.0           
#> [28] sass_0.4.10         otel_0.2.0          cli_3.6.6          
#> [31] pkgdown_2.2.0       magrittr_2.0.5      phangorn_2.12.1    
#> [34] digest_0.6.39       grid_4.6.1          lifecycle_1.0.5    
#> [37] nlme_3.1-169        vctrs_0.7.3         evaluate_1.0.5     
#> [40] codetools_0.2-20    ragg_1.5.2          purrr_1.2.2        
#> [43] rmarkdown_2.31      tools_4.6.1         pkgconfig_2.0.3    
#> [46] htmltools_0.5.9
```
