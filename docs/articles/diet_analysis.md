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
library(MicrobiomeStat)
```

## Setup

We begin by loading the mouse diet study data, which is stored remotely
and can be loaded into R as a phyloseq object.

``` r

data("dietswap", package = "microbiome")
diet_temp <- subset_samples(dietswap, timepoint == 1)
diet <- subset_taxa(diet_temp, taxa_sums(diet_temp) > 0)
```

To prepare the data for visualization, we first convert the phyloseq
object into a MicrobiomeStat object using
[`mStat_convert_phyloseq_to_data_obj()`](https://cafferychen777.github.io/MicrobiomeStat/reference/mStat_convert_phyloseq_to_data_obj.html).
This enables us to apply DESeq2 normalization, which corrects for
sequencing depth differences across samples. We then transpose the
normalized feature matrix so that rows correspond to taxa and columns to
samples, as expected by [`phylobar()`](../reference/phylobar.md).

Unlike relative-abundance normalization, DESeq2 does not constrain
sample totals to equal one. This distinction highlights the flexibility
of our functions, since phylobar can accommodate both compositional and
non-compositional inputs. After normalization, we transpose the feature
matrix so that rows correspond to taxa and columns to samples, as
required by [`phylobar()`](../reference/phylobar.md).

``` r

diet_mstat <- mStat_convert_phyloseq_to_data_obj(diet)
diet_mstat <- mStat_normalize_data(diet_mstat, "DESeq")
x <- t(diet_mstat$data.obj.norm$feature.tab)
```

Next, we construct a taxonomy-based tree from the available annotations.
To ensure that ancestor and descendant nodes remain uniquely
identifiable, we add prefixes to each taxonomic level (e.g., `p_` for
phylum, `f_` for family). Without this step, different levels can share
the same names, which would prevent the taxonomy from being converted
into a valid tree. The prefixed taxonomy table is then passed to
[`taxonomy_to_tree()`](../reference/taxonomy_to_tree.md), producing a
hierarchical tree that will serve as the backbone for visualization.

``` r

taxa <- tax_table(diet) |>
    phylobar::add_prefix()
taxa <- cbind(Kingdom = "k_Bacteria", taxa)
tree <- taxonomy_to_tree(taxa)
```

At this point, we have two aligned objects ready for visualization: a
normalized abundance matrix (x) and a taxonomy-derived tree (tree). We
now pass these objects to [`phylobar()`](../reference/phylobar.md). The
resulting plot shows stacked bar charts of taxa abundances, aligned with
the hierarchical tree structure.

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
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] MicrobiomeStat_1.5.0        tibble_3.3.1               
#>  [3] rlang_1.1.7                 DESeq2_1.48.2              
#>  [5] SummarizedExperiment_1.38.1 Biobase_2.68.0             
#>  [7] MatrixGenerics_1.20.0       matrixStats_1.5.0          
#>  [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.3        
#> [11] IRanges_2.42.0              S4Vectors_0.48.0           
#> [13] BiocGenerics_0.54.1         generics_0.1.4             
#> [15] phyloseq_1.52.0             phylobar_0.99.0            
#> [17] BiocStyle_2.36.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6            phangorn_2.12.1         permute_0.9-10         
#>   [4] magrittr_2.0.4          clue_0.3-66             ade4_1.7-23            
#>   [7] otel_0.2.0              compiler_4.5.1          mgcv_1.9-4             
#>  [10] systemfonts_1.3.1       vctrs_0.7.1             reshape2_1.4.5         
#>  [13] rmutil_1.1.10           quadprog_1.5-8          stringr_1.6.0          
#>  [16] pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
#>  [19] XVector_0.48.0          pander_0.6.6            modeest_2.4.0          
#>  [22] rmarkdown_2.30          nloptr_2.2.1            UCSC.utils_1.4.0       
#>  [25] ragg_1.5.0              purrr_1.2.1             xfun_0.56              
#>  [28] cachem_1.1.0            jsonlite_2.0.0          biomformat_1.36.0      
#>  [31] rhdf5filters_1.20.0     DelayedArray_0.34.1     Rhdf5lib_1.30.0        
#>  [34] BiocParallel_1.42.2     parallel_4.5.1          cluster_2.1.8.2        
#>  [37] R6_2.6.1                bslib_0.10.0            stringi_1.8.7          
#>  [40] RColorBrewer_1.1-3      rpart_4.1.24            boot_1.3-32            
#>  [43] numDeriv_2016.8-1.1     jquerylib_0.1.4         Rcpp_1.1.1             
#>  [46] bookdown_0.46           iterators_1.0.14        knitr_1.51             
#>  [49] Matrix_1.7-4            splines_4.5.1           igraph_2.2.1           
#>  [52] tidyselect_1.2.1        abind_1.4-8             yaml_2.3.12            
#>  [55] timeDate_4052.112       vegan_2.7-2             codetools_0.2-20       
#>  [58] lattice_0.22-9          lmerTest_3.2-0          plyr_1.8.9             
#>  [61] withr_3.0.2             S7_0.2.1                stable_1.1.7           
#>  [64] evaluate_1.0.5          desc_1.4.3              survival_3.8-6         
#>  [67] Biostrings_2.76.0       pillar_1.11.1           BiocManager_1.30.27    
#>  [70] foreach_1.5.2           reformulas_0.4.4        ggplot2_4.0.2          
#>  [73] tidytree_0.4.7          scales_1.4.0            timeSeries_4052.112    
#>  [76] minqa_1.2.8             glue_1.8.0              statip_0.2.3           
#>  [79] lazyeval_0.2.2          tools_4.5.1             data.table_1.18.2.1    
#>  [82] spatial_7.3-18          lme4_1.1-38             fBasics_4052.98        
#>  [85] locfit_1.5-9.12         fs_1.6.6                fastmatch_1.1-8        
#>  [88] rhdf5_2.52.1            grid_4.5.1              tidyr_1.3.2            
#>  [91] ape_5.8-1               rbibutils_2.4.1         nlme_3.1-168           
#>  [94] GenomeInfoDbData_1.2.14 cli_3.6.5               rappdirs_0.3.4         
#>  [97] textshaping_1.0.4       S4Arrays_1.8.1          dplyr_1.2.0            
#> [100] gtable_0.3.6            yulab.utils_0.2.4       stabledist_0.7-2       
#> [103] sass_0.4.10             digest_0.6.39           SparseArray_1.8.1      
#> [106] htmlwidgets_1.6.4       farver_2.1.2            htmltools_0.5.9        
#> [109] pkgdown_2.2.0           multtest_2.64.0         lifecycle_1.0.5        
#> [112] httr_1.4.7              MASS_7.3-65
```
