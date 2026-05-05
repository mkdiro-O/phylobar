# Cell-type Hierarchy Informs COVID-19 Immunology

In December 2020, while COVID-19 was still raging across the globe, [Su
et al. ( 2020)](https://doi.org/10.1016/j.cell.2020.10.037) published
their pioneering study examining the immune and metabolic landscapes of
mild to severe COVID patients.

See the [DIVAS
project](https://byronsyun.github.io/DIVAS_COVID19_CaseStudy/) for a
more comprehensive introduction to their study design and a more
in-depth analysis of their multiomic data. Here we just took the immune
cell type composition information to illustrate how phylobar enables us
to visualise immune cell type composition changes from milder to severer
COVID-19.

``` r

library(ape)
library(dplyr)
library(phylobar)
library(readr)
library(tibble)
library(tidyr)
library(seriation)
library(zen4R)
```

First load the cell and patient metadata.

``` r

download_zenodo("10.5281/zenodo.17477876", tempdir())
#> [zen4R][INFO] ZenodoRecord - Download in sequential mode
#> [zen4R][INFO] ZenodoRecord - Will download 1 file from record '17477876' (doi: '10.5281/zenodo.17477876') - total size: 12.8 MiB
#> [zen4R][INFO] Downloading file 'su-2020.rds' - size: 12.8 MiB
#> [zen4R][INFO] File downloaded at '/private/var/folders/mt/w0x_hxms30qch78n14v9zz0h0000gn/T/RtmpXZMniw'.
#> [zen4R][INFO] ZenodoRecord - Verifying file integrity...
#> [zen4R][INFO] File 'su-2020.rds': integrity verified (md5sum: d35b5955552582e34ee8d722e96d70f5)
#> [zen4R][INFO] ZenodoRecord - End of download
all_data <- readRDS(file.path(tempdir(), "/su-2020.rds"))
cell_counts <- all_data$cell_counts
metadata <- all_data$metadata
```

Calculate cell type compositions.

``` r

x <- cell_counts |>
    dplyr::count(patient_id, majority_voting) |>
    pivot_wider(
        names_from = majority_voting,
        values_from = n,
        values_fill = 0
    ) |>
    column_to_rownames("patient_id")
x <- x / rowSums(x)
```

Based on biological knowledge, we manually define a phylo object for the
hierarchy of immune cell types.

``` r

edges_text <- "source,target
Hematopoietic stem cell,Common myeloid progenitor
Hematopoietic stem cell,Common lymphoid progenitor
Common myeloid progenitor,Megakaryocyte
Megakaryocyte,Platelet
Common myeloid progenitor,RBC
Common myeloid progenitor,Myeloblast
Myeloblast,Monocyte
Monocyte,CD14 Monocyte
Monocyte,CD16 Monocyte
Myeloblast,SC & Eosinophil
Myeloblast,DC
Myeloblast,pDC
Common lymphoid progenitor,NK
Common lymphoid progenitor,Small lymphocyte
Small lymphocyte,T lymphocyte
T lymphocyte,CD4 T lymphocyte
T lymphocyte,CD8 T lymphocyte
CD4 T lymphocyte,CD4 T
CD4 T lymphocyte,CD4m T
CD4 T lymphocyte,CD4n T
CD8 T lymphocyte,CD8m T
CD8 T lymphocyte,CD8eff T
T lymphocyte,gd T
Small lymphocyte,B lymphocyte
B lymphocyte,B
B lymphocyte,IgA PB
B lymphocyte,IgG PB"

edge <- read_csv(edges_text)
tips <- setdiff(edge$target, edge$source)

# preparation for phylo construction
internal_nodes <- setdiff(unlist(edge), tips)
node_order <- c(tips, internal_nodes)
ix <- setNames(seq_along(node_order), node_order)

# create phylo
tree <- list(
    edge = matrix(c(ix[edge$source], ix[edge$target]), ncol = 2),
    Nnode = length(internal_nodes),
    tip.label = node_order[seq_along(tips)],
    node.label = node_order[seq(length(tips) + 1, length(ix))]
)
class(tree) <- "phylo"
```

Let’s order patients by COVID severity (mild, moderate and severe).
Then, within each severity group, use a hierarchical clustering leaf
order to place samples with similar cell type compositions closer.

``` r

md <- metadata |>
    distinct(clinical_id, severity) |>
    mutate(severity = factor(
        severity, levels = c("mild", "moderate", "severe")
    )) |>
    arrange(severity)

patient_order <- md %>%
    split(.$severity) |>
    lapply(\(df) {
        ids <- df$clinical_id
        m <- x[ids, , drop = FALSE]
        ids[hclust(dist(m))$order]
    }) |>
    unlist(use.names = FALSE)
```

Now we can make a phylobar visualization. We can see that T cells and
monocytes made up the majority of cellular populations in patient
samples. We can observe an overall increase of monocytes and decrease of
NK and T cells, especially CD8 T cells, in moderate and severe COVID
patients. These changes were found statistically significant in Su et
al. (2020), and where highlighted in their Fig 1.

Another feature that was not characterised in Su et al. (2020) was the
presence of gd (γδ, i.e. gamma delta) T cells in milder COVID patients
exclusively. This echos some more recent findings suggesting that gd T
cells may play an important role in the immune response to the
pathogenic SARS-Cov-2 virus ([Massow et al.,
2021](https://doi.org/10.3389/fimmu.2021.741218); [Terzoli et al.,
2024](https://www.nature.com/articles/s41541-024-00853-9)).

``` r

phylobar(x[patient_order, ], tree)
```

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
#> [1] zen4R_0.10.5     seriation_1.5.8  tidyr_1.3.2      tibble_3.3.1    
#> [5] readr_2.2.0      phylobar_0.99.11 dplyr_1.2.1      ape_5.8-1       
#> [9] BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] fastmatch_1.1-8     xfun_0.57           bslib_0.10.0       
#>  [4] htmlwidgets_1.6.4   lattice_0.22-9      tzdb_0.5.0         
#>  [7] quadprog_1.5-8      vctrs_0.7.3         tools_4.6.0        
#> [10] generics_0.1.4      curl_7.1.0          parallel_4.6.0     
#> [13] ca_0.71.1           cluster_2.1.8.2     pkgconfig_2.0.3    
#> [16] Matrix_1.7-5        desc_1.4.3          lifecycle_1.0.5    
#> [19] compiler_4.6.0      textshaping_1.0.5   keyring_1.4.1      
#> [22] codetools_0.2-20    htmltools_0.5.9     sass_0.4.10        
#> [25] yaml_2.3.12         crayon_1.5.3        pillar_1.11.1      
#> [28] pkgdown_2.2.0       jquerylib_0.1.4     cachem_1.1.0       
#> [31] iterators_1.0.14    TSP_1.2.7           foreach_1.5.2      
#> [34] nlme_3.1-169        phangorn_2.12.1     tidyselect_1.2.1   
#> [37] digest_0.6.39       purrr_1.2.2         bookdown_0.46      
#> [40] fastmap_1.2.0       grid_4.6.0          cli_3.6.6          
#> [43] magrittr_2.0.5      XML_3.99-0.23       utf8_1.2.6         
#> [46] withr_3.0.2         bit64_4.8.0         registry_0.5-1     
#> [49] rmarkdown_2.31      httr_1.4.8          bit_4.6.0          
#> [52] igraph_2.3.1        otel_0.2.0          ragg_1.5.2         
#> [55] hms_1.1.4           evaluate_1.0.5      knitr_1.51         
#> [58] rlang_1.2.0         Rcpp_1.1.1-1.1      glue_1.8.1         
#> [61] BiocManager_1.30.27 xml2_1.5.2          vroom_1.7.1        
#> [64] jsonlite_2.0.0      R6_2.6.1            plyr_1.8.9         
#> [67] systemfonts_1.3.2   fs_2.1.0
```
