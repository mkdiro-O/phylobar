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

By default, the
[`subset_cluster()`](http://mkdiro-o.github.io/phylobar/reference/subset_cluster.md)
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
#> [1] microbiome_1.34.0 ggplot2_4.0.3     phyloseq_1.56.0   phylobar_0.99.12 
#> [5] ape_5.8-1         BiocStyle_2.40.0 
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-24                 tidyselect_1.2.1           
#>  [3] dplyr_1.2.1                 farver_2.1.2               
#>  [5] Biostrings_2.80.0           S7_0.2.2                   
#>  [7] fastmap_1.2.0               digest_0.6.39              
#>  [9] lifecycle_1.0.5             cluster_2.1.8.2            
#> [11] survival_3.8-6              magrittr_2.0.5             
#> [13] compiler_4.6.0              rlang_1.2.0                
#> [15] sass_0.4.10                 tools_4.6.0                
#> [17] igraph_2.3.1                yaml_2.3.12                
#> [19] data.table_1.18.2.1         knitr_1.51                 
#> [21] phangorn_2.12.1             S4Arrays_1.12.0            
#> [23] htmlwidgets_1.6.4           DelayedArray_0.38.1        
#> [25] plyr_1.8.9                  RColorBrewer_1.1-3         
#> [27] abind_1.4-8                 Rtsne_0.17                 
#> [29] withr_3.0.2                 purrr_1.2.2                
#> [31] BiocGenerics_0.58.0         desc_1.4.3                 
#> [33] grid_4.6.0                  stats4_4.6.0               
#> [35] multtest_2.68.0             biomformat_1.40.0          
#> [37] scales_1.4.0                iterators_1.0.14           
#> [39] MASS_7.3-65                 SummarizedExperiment_1.42.0
#> [41] cli_3.6.6                   vegan_2.7-3                
#> [43] rmarkdown_2.31              crayon_1.5.3               
#> [45] ragg_1.5.2                  generics_0.1.4             
#> [47] otel_0.2.0                  reshape2_1.4.5             
#> [49] cachem_1.1.0                stringr_1.6.0              
#> [51] splines_4.6.0               parallel_4.6.0             
#> [53] BiocManager_1.30.27         XVector_0.52.0             
#> [55] matrixStats_1.5.0           vctrs_0.7.3                
#> [57] Matrix_1.7-5                jsonlite_2.0.0             
#> [59] bookdown_0.46               IRanges_2.46.0             
#> [61] S4Vectors_0.50.0            systemfonts_1.3.2          
#> [63] foreach_1.5.2               tidyr_1.3.2                
#> [65] jquerylib_0.1.4             glue_1.8.1                 
#> [67] pkgdown_2.2.0               codetools_0.2-20           
#> [69] stringi_1.8.7               gtable_0.3.6               
#> [71] GenomicRanges_1.64.0        quadprog_1.5-8             
#> [73] tibble_3.3.1                pillar_1.11.1              
#> [75] htmltools_0.5.9             Seqinfo_1.2.0              
#> [77] R6_2.6.1                    textshaping_1.0.5          
#> [79] evaluate_1.0.5              lattice_0.22-9             
#> [81] Biobase_2.72.0              bslib_0.10.0               
#> [83] Rcpp_1.1.1-1.1              fastmatch_1.1-8            
#> [85] permute_0.9-10              SparseArray_1.12.2         
#> [87] nlme_3.1-169                mgcv_1.9-4                 
#> [89] xfun_0.57                   fs_2.1.0                   
#> [91] MatrixGenerics_1.24.0       pkgconfig_2.0.3
```
