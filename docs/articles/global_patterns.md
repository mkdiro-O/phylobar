# Working with phyloseq and TreeSummarizedExperiment versions of the Global Patterns Dataset

This vignette demonstrates how to create phylobar plots using either
[phyloseq](https://bioconductor.org/packages/phyloseq/) or
[TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)
objects, two popular Bioconductor containers for microbiome data. We’ll
work with the same dataset (GlobalPatterns) stored in both formats.

## Setup

``` r

library(ape)
library(phylobar)
library(phyloseq)
library(miaViz)
library(scater)
```

## phyloseq Inputs

We will first study how to use phylobar on phyloseq objects. The block
below sets the global pattern status and subsets to the Chlamydia
subtree. We can extract the relevant sample composition using the
otu_table accessor function.

``` r

data(GlobalPatterns, package = "phyloseq")
chlamydiae <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
x_phylo <- t(otu_table(chlamydiae))
```

### Phylogenetic Tree Tree

Our first approach will use the phylogenetic tree to support
hierarchical interaction. The tree can be extracted from the phyloseq
object using `phy_tree`. We first subsetted to only the chlamydia
subtree. There are now some samples without any observed counts.
Therefore, we will remove the samples from the view.

``` r

tree_phylo <- phy_tree(chlamydiae)
x_phylo <- x_phylo[rowSums(x_phylo) > 0, ]
x_phylo <- x_phylo / rowSums(x_phylo)
```

We can now create our phylobar plot.

``` r

phylobar(x_phylo, tree_phylo, width = 900)
```

### Using Taxonomic Hierarchy

This creates the same type of visualization but using a taxonomic
hierarchy rather than a phylogeny. We need some problem-specific
preprocessing to deal with the fact that some ancestor and descendant
names look identical in this dataset. The `add_prefix` helper ensures
that the names are distinct between parents and children by adding in
the taxonomic rank to each name.

``` r

taxa <- tax_table(chlamydiae) |>
    as.data.frame()
taxa$ASV <- rownames(taxa)
taxa <- phylobar::add_prefix(taxa)
```

Now we can create a phylo tree object associated with the original
taxonomy. Note that this is not a binary tree, since several taxonomic
categories descended from a single node.

``` r

tax_tree <- taxonomy_to_tree(taxa)
plot(tax_tree)
```

![](global_patterns_files/figure-html/construct_taxonomy-1.png)

Given our new taxonomy and the earlier sample information, we can now
create a phylobar plot.

``` r

phylobar(x_phylo, tax_tree, width = 900)
```

## TreeSummarizedExperiment Inputs

This example comes from the [miaViz
documentation](https://microbiome.github.io/miaViz/articles/miaViz.html)
and shows how to create a stacked bar tree from a
TreeSummarizedExperiment. We’ll filter down to only those species that
are present in at least one percent of all samples. In this example
we’ll use a relative abundance transformation, so that we left with a
composition of our plot rather than a more general stacked bar plot.

``` r

data(GlobalPatterns, package = "mia")
prev_species <- getPrevalent(GlobalPatterns, rank = "Species", detection = 0.01)
GlobalPatterns_tse <- GlobalPatterns[
    rowData(GlobalPatterns)$Species %in% prev_species,
]
GlobalPatterns_tse <- transformAssay(
    GlobalPatterns_tse,
    method = "relabundance"
)
```

Let’s keep only those taxes that are observed in at least one sample. We
also need to filter the tree to reflect this reduced subset of taxa.

``` r

x_tse <- t(assay(GlobalPatterns_tse, "relabundance"))
x_tse <- x_tse[, colSums(x_tse) > 0]
tree_tse <- rowTree(GlobalPatterns_tse) |>
    keep.tip(colnames(x_tse))
```

With these inputs, we can generate the phylobar visualization. This view
didn’t subset to the chlamydia class, which is why we have a larger tree
here. But hopefully it’s clear that TreeSummarizedExperiment and
phyloseq can essentially be used interchangeably when constructing the
necessary inputs for these visualizations.

``` r

phylobar(x_tse, tree_tse, width = 900)
```

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
#>  [1] scater_1.38.1                   scuttle_1.20.0                 
#>  [3] miaViz_1.18.1                   mia_1.18.0                     
#>  [5] TreeSummarizedExperiment_2.18.0 Biostrings_2.78.0              
#>  [7] XVector_0.50.0                  SingleCellExperiment_1.32.0    
#>  [9] MultiAssayExperiment_1.36.1     SummarizedExperiment_1.40.0    
#> [11] Biobase_2.70.0                  GenomicRanges_1.62.1           
#> [13] Seqinfo_1.0.0                   IRanges_2.44.0                 
#> [15] S4Vectors_0.48.0                BiocGenerics_0.56.0            
#> [17] generics_0.1.4                  MatrixGenerics_1.22.0          
#> [19] matrixStats_1.5.0               ggraph_2.2.2                   
#> [21] ggplot2_4.0.2                   phyloseq_1.54.2                
#> [23] phylobar_0.99.10                ape_5.8-1                      
#> [25] BiocStyle_2.38.0               
#> 
#> loaded via a namespace (and not attached):
#>   [1] splines_4.5.2               ggplotify_0.1.3            
#>   [3] cellranger_1.1.0            tibble_3.3.1               
#>   [5] polyclip_1.10-7             DirichletMultinomial_1.52.0
#>   [7] lifecycle_1.0.5             lattice_0.22-9             
#>   [9] MASS_7.3-65                 SnowballC_0.7.1            
#>  [11] magrittr_2.0.5              sass_0.4.10                
#>  [13] rmarkdown_2.31              jquerylib_0.1.4            
#>  [15] yaml_2.3.12                 otel_0.2.0                 
#>  [17] DBI_1.3.0                   RColorBrewer_1.1-3         
#>  [19] ade4_1.7-24                 multcomp_1.4-30            
#>  [21] abind_1.4-8                 quadprog_1.5-8             
#>  [23] purrr_1.2.2                 fillpattern_1.0.3          
#>  [25] yulab.utils_0.2.4           TH.data_1.1-5              
#>  [27] tweenr_2.0.3                rappdirs_0.3.4             
#>  [29] sandwich_3.1-1              gdtools_0.5.0              
#>  [31] ggrepel_0.9.8               tokenizers_0.3.0           
#>  [33] irlba_2.3.7                 tidytree_0.4.7             
#>  [35] vegan_2.7-3                 rbiom_2.2.1                
#>  [37] parallelly_1.46.1           pkgdown_2.2.0              
#>  [39] permute_0.9-10              DelayedMatrixStats_1.32.0  
#>  [41] codetools_0.2-20            DelayedArray_0.36.0        
#>  [43] ggtext_0.1.2                xml2_1.5.2                 
#>  [45] ggforce_0.5.0               tidyselect_1.2.1           
#>  [47] aplot_0.2.9                 farver_2.1.2               
#>  [49] ScaledMatrix_1.18.0         viridis_0.6.5              
#>  [51] jsonlite_2.0.0              BiocNeighbors_2.4.0        
#>  [53] multtest_2.66.0             decontam_1.30.0            
#>  [55] tidygraph_1.3.1             survival_3.8-6             
#>  [57] iterators_1.0.14            emmeans_2.0.2              
#>  [59] systemfonts_1.3.2           foreach_1.5.2              
#>  [61] tools_4.5.2                 ggnewscale_0.5.2           
#>  [63] treeio_1.34.0               ragg_1.5.2                 
#>  [65] Rcpp_1.1.1-1.1              glue_1.8.1                 
#>  [67] gridExtra_2.3               SparseArray_1.10.9         
#>  [69] BiocBaseUtils_1.12.0        xfun_0.57                  
#>  [71] mgcv_1.9-4                  dplyr_1.2.0                
#>  [73] withr_3.0.2                 BiocManager_1.30.27        
#>  [75] fastmap_1.2.0               rhdf5filters_1.22.0        
#>  [77] bluster_1.20.0              digest_0.6.39              
#>  [79] rsvd_1.0.5                  gridGraphics_0.5-1         
#>  [81] R6_2.6.1                    estimability_1.5.1         
#>  [83] textshaping_1.0.5           dichromat_2.0-0.1          
#>  [85] tidyr_1.3.2                 fontLiberation_0.1.0       
#>  [87] data.table_1.18.2.1         DECIPHER_3.6.0             
#>  [89] graphlayouts_1.2.3          htmlwidgets_1.6.4          
#>  [91] S4Arrays_1.10.1             pkgconfig_2.0.3            
#>  [93] gtable_0.3.6                S7_0.2.1                   
#>  [95] janeaustenr_1.0.0           htmltools_0.5.9            
#>  [97] fontBitstreamVera_0.1.1     bookdown_0.46              
#>  [99] biomformat_1.38.3           scales_1.4.0               
#> [101] ggfun_0.2.0                 knitr_1.51                 
#> [103] tzdb_0.5.0                  reshape2_1.4.5             
#> [105] coda_0.19-4.1               nlme_3.1-168               
#> [107] cachem_1.1.0                zoo_1.8-15                 
#> [109] rhdf5_2.54.1                stringr_1.6.0              
#> [111] parallel_4.5.2              vipor_0.4.7                
#> [113] desc_1.4.3                  pillar_1.11.1              
#> [115] grid_4.5.2                  vctrs_0.7.3                
#> [117] slam_0.1-55                 BiocSingular_1.26.1        
#> [119] beachmat_2.26.0             xtable_1.8-8               
#> [121] cluster_2.1.8.2             beeswarm_0.4.0             
#> [123] evaluate_1.0.5              readr_2.2.0                
#> [125] mvtnorm_1.3-6               cli_3.6.6                  
#> [127] compiler_4.5.2              rlang_1.2.0                
#> [129] crayon_1.5.3                tidytext_0.4.3             
#> [131] plyr_1.8.9                  fs_2.1.0                   
#> [133] ggbeeswarm_0.7.3            ggiraph_0.9.6              
#> [135] stringi_1.8.7               viridisLite_0.4.3          
#> [137] BiocParallel_1.44.0         lazyeval_0.2.2             
#> [139] fontquiver_0.2.1            Matrix_1.7-5               
#> [141] hms_1.1.4                   patchwork_1.3.2            
#> [143] sparseMatrixStats_1.22.0    Rhdf5lib_1.32.0            
#> [145] gridtext_0.1.6              igraph_2.3.0               
#> [147] memoise_2.0.1               bslib_0.10.0               
#> [149] ggtree_4.0.5                phangorn_2.12.1            
#> [151] fastmatch_1.1-8             readxl_1.4.5
```
