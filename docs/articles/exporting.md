# Exporting to Vector Graphics Format

## Introduction

Once you identify an interesting view, how can you share it with others?
First, you could press the Escape key to freeze the view, so that it
doesn’t respond to any other mouse interactions. From there, you could
always take a screenshot, but if you want to preserve more details (or
add some notes within the figure itself), you could instead export it to
a Scalable Graphics (SVG), as discussed in this vignette. For this to
work, you will need to use Google Chrome and install the bookmarklet
[SVG Crowbar 2](https://nytimes.github.io/svg-crowbar/). We’ll
illustrate the process using the simple randomly generatred tree below.

``` r

library(ape)
library(phylobar)

tree <- rtree(20)
samples <- matrix(rpois(100 * 20, 1), 100, 20)
phylobar(samples, tree, width=800)
```

## Step 1: Exporting from Browser

Once you have the bookmarklet installed, you can press it to download
any of the phylobar figures in a Quarto or Rmarkdown notebook. In the
recording below, we pressed Escape after painting the blue subtree. Then
we could move to our bookmarks bar without changing the view. Clicking
the bookmarklet brings up a “SVG \# (Download)” button that you can
click, and the associated graphic will appear in your Downloads folder.
Note that you can scroll down to other phylobar visualizations and you
would see other buttons for downloading those plots.

# An error occurred.

Unable to execute JavaScript.

## Step 2: Open in Editor

The downloaded file is an ordinary SVG, which means it can be edited in
any image editor designed for vector graphics. For example, in the
recording below, we open the downloaded figure using Inkscape, a free
and open source vector graphics editor. Notice that the components of
the phylobar plot are already organized into layers.. For example, there
is a layer for the barplot, sublayers for each sample, and a final
sublayer for the taxa within the samples. This layer structure reflects
the original DOM structure created by the D3 visualization code. The
recording below shows how we can navigate this layer structure within
Inkscape.

# An error occurred.

Unable to execute JavaScript.

## Step 3: Modify and Save

From here, we could edit if we want (e.g., to add a comment about an
interesting pattern). We are then able to export into other formats that
might be more accessible to others, e.g., one that we can embed in
slides or add to a word document. In the recording below, we use the
layer structure shown in the previous snippet to modify the green
subtree so that it uses dashed (rather than solid) lines to link the
nodes. Once we have our final version, we export the graphic to PNG
format.

# An error occurred.

Unable to execute JavaScript.

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Sequoia 15.7.4
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] phylobar_0.99.10 ape_5.8-1        BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5        jsonlite_2.0.0      compiler_4.6.0     
#>  [4] BiocManager_1.30.27 Rcpp_1.1.1-1.1      parallel_4.6.0     
#>  [7] cluster_2.1.8.2     jquerylib_0.1.4     systemfonts_1.3.2  
#> [10] textshaping_1.0.5   yaml_2.3.12         fastmap_1.2.0      
#> [13] lattice_0.22-9      R6_2.6.1            generics_0.1.4     
#> [16] igraph_2.3.1        knitr_1.51          htmlwidgets_1.6.4  
#> [19] bookdown_0.46       desc_1.4.3          bslib_0.10.0       
#> [22] rlang_1.2.0         fastmatch_1.1-8     cachem_1.1.0       
#> [25] xfun_0.57           quadprog_1.5-8      fs_2.1.0           
#> [28] sass_0.4.10         otel_0.2.0          cli_3.6.6          
#> [31] pkgdown_2.2.0       magrittr_2.0.5      phangorn_2.12.1    
#> [34] digest_0.6.39       grid_4.6.0          lifecycle_1.0.5    
#> [37] nlme_3.1-169        vctrs_0.7.3         evaluate_1.0.5     
#> [40] codetools_0.2-20    ragg_1.5.2          purrr_1.2.2        
#> [43] rmarkdown_2.31      tools_4.6.0         pkgconfig_2.0.3    
#> [46] htmltools_0.5.9
```
