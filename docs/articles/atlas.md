# Tipping Points in the Atlas 1006 Dataset

The Human Intestinal Tract (HITChip) Atlas dataset is a widely used
reference for studying variation in the gut microbiome across diverse
individuals. The dataset is available in the microbiome package as
atlas1006, and includes microbial abundances along with metadata such as
extraction method, subject information, and taxonomic annotations.

In this vignette, we demonstrate how to use phylobar to create
phylogeny-aware visualizations of the atlas data. We will walk through
tree construction, highlight potential pitfalls in taxonomic labeling,
and generate visualizations that summarize variation across samples.

``` r

library(ape)
library(phylobar)
library(microbiome)
#> Warning: package 'ggplot2' was built under R version 4.5.2
```

## Setup

We begin by loading the atlas data, which is available through the
microbiome package as a
[phyloseq](https://bioconductor.org/packages/phyloseq/) object (a
Bioconductor container for microbiome data).

``` r

data(atlas1006)
pseq3 <- atlas1006 |>
    subset_samples(DNA_extraction_method == "r") |>
    transform(transform = "compositional")
```

Next, we load the atlas1006 dataset. For simplicity, we focus on samples
with the r DNA extraction method and transform abundances into
compositional form (so each sample sums to 1).

``` r

tree <- taxonomy_to_tree(tax_table(pseq3))
checkValidPhylo(tree)
#> Starting checking the validity of tree...
#> Found number of tips: n = 128 
#> Found number of nodes: m = 23 
#>   FATAL: nodes and tips should appear only once in the 2nd column of 'edge'
#>   FATAL: the root node should not appear in the 2nd column of 'edge'
#> Done.
```

A naive attempt to build a taxonomy-based tree can fail, because
duplicate names across different levels of the taxonomy may result in
loops. To address this, we add level-specific prefixes (e.g., p\_ for
phylum, f\_ for family), ensuring that ancestor and descendant nodes
have unique labels.

``` r

taxa <- tax_table(pseq3)
taxa <- cbind(Kingdom = "Bacteria", taxa)
taxa <- phylobar::add_prefix(taxa)

tree <- taxonomy_to_tree(taxa)
checkValidPhylo(tree)
#> Starting checking the validity of tree...
#> Found number of tips: n = 130 
#> Found number of nodes: m = 31 
#> Done.
```

With a valid tree and abundance matrix, we can generate a phylobar plot.
Here, samples are arranged along the x-axis, with stacked bars
representing taxonomic abundances, aligned with the corresponding
branches of the tree.

``` r

x <- t(otu_table(pseq3))
phylobar(x, tree, sample_show_all = FALSE, rel_width = 0.2)
```

Because the full dataset is complex, we can simplify the visualization
by subsampling representative samples. This preserves the main
conclusions while reducing visual clutter.

``` r

x_sub <- subset_cluster(x)
x_sub <- x_sub[, colSums(x_sub) > 0]

leaves_to_keep <- intersect(tree$tip.label, colnames(x_sub))
filtered_tree <- drop.tip(tree, setdiff(tree$tip.label, leaves_to_keep))
phylobar(x_sub, filtered_tree)
```

By default, the [`subset_cluster()`](../reference/subset_cluster.md)
function subsets to samples that are close to a cut hierarchical
clustering tree’s centroids. Alternatively, setting `method = "medoid"`
will choose representatives using $`K`$-medoids according to Bray-Curtis
distance. The result is shown below, though our interpretation will
focus on the hierarchical clustering version.

``` r

x_sub_bc <- subset_cluster(x, method = "medoid")
x_sub_bc <- x_sub[, colSums(x_sub) > 0]

leaves_to_keep <- intersect(tree$tip.label, colnames(x_sub_bc))
filtered_tree <- drop.tip(tree, setdiff(tree$tip.label, leaves_to_keep))
phylobar(x_sub_bc, filtered_tree)
```

## Interpretation

There are many interesting patterns in the atlas dataset revealing
distinct patterns in the human intestinal. For example, Enterobacter
aerogenes appears across a broad range of samples, but generally at low
relative abundance. Only one sample shows notable representation,
creating a narrow band in the stacked barplot. This patchy distribution
is a characteristic of opportunistic or conditionally rare taxa, which
may be present in many individuals but rarely dominate the community
structure. Such patterns have been previously been described as examples
of bimodal or “tipping element” behavior, where a taxon is abundant in
some hosts but nearly absent in others.
(<https://pmc.ncbi.nlm.nih.gov/articles/PMC4102116/pdf/ncomms5344.pdf>)

![](https://i.imgur.com/uj1avu1.png)

Next, this screenshot highlights the dominance of Firmicutes (yellow)
and Bacteroidetes (purple) in the human gut microbiome. While these two
phyla are prevalent across nearly all samples, their relative
proportions vary widely. Some samples are Firmicutes-dominant, while
others are Bacteroidetes-dominant. This variability reflects the complex
interplay of host genetics, diet, and environmental factors that shape
the gut microbiome.

![](https://i.imgur.com/qUZoXZh.png)

## Session Info

``` r

sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Tahoe 26.3
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] microbiome_1.30.0 ggplot2_4.0.2     phyloseq_1.52.0   phylobar_0.99.0  
#> [5] ape_5.8-1         BiocStyle_2.36.0 
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-23             tidyselect_1.2.1        dplyr_1.2.0            
#>  [4] farver_2.1.2            Biostrings_2.76.0       S7_0.2.1               
#>  [7] fastmap_1.2.0           digest_0.6.39           lifecycle_1.0.5        
#> [10] cluster_2.1.8.2         survival_3.8-6          magrittr_2.0.4         
#> [13] compiler_4.5.1          rlang_1.1.7             sass_0.4.10            
#> [16] tools_4.5.1             igraph_2.2.1            yaml_2.3.12            
#> [19] data.table_1.18.2.1     knitr_1.51              phangorn_2.12.1        
#> [22] htmlwidgets_1.6.4       plyr_1.8.9              RColorBrewer_1.1-3     
#> [25] Rtsne_0.17              withr_3.0.2             purrr_1.2.1            
#> [28] BiocGenerics_0.54.1     desc_1.4.3              grid_4.5.1             
#> [31] stats4_4.5.1            multtest_2.64.0         biomformat_1.36.0      
#> [34] Rhdf5lib_1.30.0         scales_1.4.0            iterators_1.0.14       
#> [37] MASS_7.3-65             cli_3.6.5               vegan_2.7-2            
#> [40] rmarkdown_2.30          crayon_1.5.3            ragg_1.5.0             
#> [43] generics_0.1.4          otel_0.2.0              httr_1.4.7             
#> [46] reshape2_1.4.5          cachem_1.1.0            rhdf5_2.52.1           
#> [49] stringr_1.6.0           splines_4.5.1           parallel_4.5.1         
#> [52] BiocManager_1.30.27     XVector_0.48.0          vctrs_0.7.1            
#> [55] Matrix_1.7-4            jsonlite_2.0.0          bookdown_0.46          
#> [58] IRanges_2.42.0          S4Vectors_0.48.0        systemfonts_1.3.1      
#> [61] foreach_1.5.2           tidyr_1.3.2             jquerylib_0.1.4        
#> [64] glue_1.8.0              pkgdown_2.2.0           codetools_0.2-20       
#> [67] stringi_1.8.7           gtable_0.3.6            GenomeInfoDb_1.44.3    
#> [70] UCSC.utils_1.4.0        quadprog_1.5-8          tibble_3.3.1           
#> [73] pillar_1.11.1           htmltools_0.5.9         rhdf5filters_1.20.0    
#> [76] GenomeInfoDbData_1.2.14 R6_2.6.1                textshaping_1.0.4      
#> [79] evaluate_1.0.5          lattice_0.22-9          Biobase_2.68.0         
#> [82] bslib_0.10.0            Rcpp_1.1.1              fastmatch_1.1-8        
#> [85] permute_0.9-10          nlme_3.1-168            mgcv_1.9-4             
#> [88] xfun_0.56               fs_1.6.6                pkgconfig_2.0.3
```
