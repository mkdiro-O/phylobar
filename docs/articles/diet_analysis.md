# Differences in the Community Composition due to Diet

Diet is a key factor shaping the gut microbiome. In a recent mouse
study, researchers profiled gut microbial communities under different
dietary conditions to examine how nutrient intake influences community
structure. The dataset, available as a
[phyloseq](https://bioconductor.org/packages/phyloseq/) object (a
Bioconductor container for microbiome data), includes microbial feature
tables, taxonomic annotations, and sample metadata.

In this vignette, we use these data to demonstrate how phylobar can
produce phylogenetically informed visualizations of microbiome
composition. We begin by normalizing feature abundances and
reconstructing the taxonomic tree. We then generate a phylobar plot to
highlight community shifts across diet groups, with a focus on clades
where dietary effects are most pronounced.

``` r

library(phylobar)
library(phyloseq)
library(DESeq2)
```

## Setup

We begin by loading the mouse diet study data, which is stored remotely
and can be loaded into R as a phyloseq object.

``` r

data("dietswap", package = "microbiome")
diet_temp <- subset_samples(dietswap, timepoint == 1)
diet <- subset_taxa(diet_temp, taxa_sums(diet_temp) > 0)
```

To prepare the data for visualization, we first normalize the phyloseq
count matrix using
[`deseq_normalize()`](http://mkdiro-o.github.io/phylobar/reference/deseq_normalize.md),
which corrects for sequencing depth differences across samples using the
DESeq size-factor method. We then transpose the normalized feature
matrix so that rows correspond to taxa and columns to samples, as
expected by
[`phylobar()`](http://mkdiro-o.github.io/phylobar/reference/phylobar.md).

Unlike relative-abundance normalization, DESeq2 does not constrain
sample totals to equal one. This distinction highlights the flexibility
of our functions, since phylobar can accommodate both compositional and
non-compositional inputs. After normalization, we transpose the feature
matrix so that rows correspond to taxa and columns to samples, as
required by
[`phylobar()`](http://mkdiro-o.github.io/phylobar/reference/phylobar.md).

``` r

otu <- as(otu_table(diet), "matrix")
x <- t(deseq_normalize(otu))
```

Next, we construct a taxonomy-based tree from the available annotations.
To ensure that ancestor and descendant nodes remain uniquely
identifiable, we add prefixes to each taxonomic level (e.g., `p_` for
phylum, `f_` for family). Without this step, different levels can share
the same names, which would prevent the taxonomy from being converted
into a valid tree. The prefixed taxonomy table is then passed to
[`taxonomy_to_tree()`](http://mkdiro-o.github.io/phylobar/reference/taxonomy_to_tree.md),
producing a hierarchical tree that will serve as the backbone for
visualization.

``` r

taxa <- tax_table(diet) |>
    phylobar::add_prefix()
taxa <- cbind(Kingdom = "k_Bacteria", taxa)
tree <- taxonomy_to_tree(taxa)
```

At this point, we have two aligned objects ready for visualization: a
normalized abundance matrix (x) and a taxonomy-derived tree (tree). We
now pass these objects to
[`phylobar()`](http://mkdiro-o.github.io/phylobar/reference/phylobar.md).
The resulting plot shows stacked bar charts of taxa abundances, aligned
with the hierarchical tree structure.

``` r

phylobar(x, tree, width = 800)
```

## Interpretation

We showcase several notable patterns using screenshots from the
interactive plots: most mice display a co-dominance of Firmicutes
(purple) and Bacteroidetes (teal), but sample_52 shows a clear
Firmicutes dominance, underscoring how dietary interventions can shift
the balance between these two phyla.

![](https://i.imgur.com/6ZYEBEA.png)

Next, we observe the abundance of the Clostridium cluster IV family
acros samples. The variability in its abundance suggests that while
Clostridium cluster IV is widespread, its contribution to community
structure differs across individuals or conditions.

![](https://i.imgur.com/RVot4qv.png)

## Session Info

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.7.4
#> 
#> Matrix products: default
#> BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: America/Chicago
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] DESeq2_1.50.2               SummarizedExperiment_1.40.0
#>  [3] Biobase_2.70.0              MatrixGenerics_1.22.0      
#>  [5] matrixStats_1.5.0           GenomicRanges_1.62.1       
#>  [7] Seqinfo_1.0.0               IRanges_2.44.0             
#>  [9] S4Vectors_0.48.0            BiocGenerics_0.56.0        
#> [11] generics_0.1.4              phyloseq_1.54.2            
#> [13] phylobar_0.99.10            BiocStyle_2.38.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-24         tidyselect_1.2.1    dplyr_1.2.0        
#>  [4] farver_2.1.2        Biostrings_2.78.0   S7_0.2.1           
#>  [7] fastmap_1.2.0       digest_0.6.39       lifecycle_1.0.5    
#> [10] cluster_2.1.8.2     survival_3.8-6      magrittr_2.0.5     
#> [13] compiler_4.5.2      rlang_1.2.0         sass_0.4.10        
#> [16] tools_4.5.2         igraph_2.3.0        yaml_2.3.12        
#> [19] data.table_1.18.2.1 knitr_1.51          phangorn_2.12.1    
#> [22] S4Arrays_1.10.1     htmlwidgets_1.6.4   DelayedArray_0.36.0
#> [25] plyr_1.8.9          RColorBrewer_1.1-3  BiocParallel_1.44.0
#> [28] abind_1.4-8         purrr_1.2.2         desc_1.4.3         
#> [31] grid_4.5.2          multtest_2.66.0     biomformat_1.38.3  
#> [34] Rhdf5lib_1.32.0     ggplot2_4.0.2       scales_1.4.0       
#> [37] iterators_1.0.14    MASS_7.3-65         dichromat_2.0-0.1  
#> [40] cli_3.6.6           rmarkdown_2.31      vegan_2.7-3        
#> [43] crayon_1.5.3        ragg_1.5.2          otel_0.2.0         
#> [46] reshape2_1.4.5      ape_5.8-1           cachem_1.1.0       
#> [49] rhdf5_2.54.1        stringr_1.6.0       splines_4.5.2      
#> [52] parallel_4.5.2      BiocManager_1.30.27 XVector_0.50.0     
#> [55] vctrs_0.7.3         Matrix_1.7-5        jsonlite_2.0.0     
#> [58] bookdown_0.46       systemfonts_1.3.2   locfit_1.5-9.12    
#> [61] foreach_1.5.2       jquerylib_0.1.4     glue_1.8.1         
#> [64] pkgdown_2.2.0       codetools_0.2-20    stringi_1.8.7      
#> [67] gtable_0.3.6        quadprog_1.5-8      tibble_3.3.1       
#> [70] pillar_1.11.1       rhdf5filters_1.22.0 htmltools_0.5.9    
#> [73] R6_2.6.1            textshaping_1.0.5   evaluate_1.0.5     
#> [76] lattice_0.22-9      bslib_0.10.0        Rcpp_1.1.1-1.1     
#> [79] fastmatch_1.1-8     SparseArray_1.10.9  nlme_3.1-168       
#> [82] permute_0.9-10      mgcv_1.9-4          xfun_0.57          
#> [85] fs_2.1.0            pkgconfig_2.0.3
```
