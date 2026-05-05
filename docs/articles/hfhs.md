# The effect of a high-fat high-sugar (HFHS) diet on the mouse microbiome

## Overview

High-fat, high-sugar (HFHS) diets are known to profoundly alter the gut
microbiome, influencing both microbial composition and functional
capacity. As described by Susin et al. (2020), the effects of diet on
the gut microbiome were evaluated in 47 C57BL/6 female mice. The mice
were fed either an HFHS diet or a normal diet, and fecal samples were
collected on days 0, 1, 4, and 7. Illumina MiSeq was used to generate
16S rRNA sequencing data, which were subsequently processed with QIIME
1.9.0. In this vignette, we use the preprocessed dataset from Kodikara
et al. (2025). The dataset can be accessed through the GitHub repository
of the [`LUPINE` R package](https://github.com/SarithaKodikara/LUPINE)

Here, we demonstrate how to apply the phylobar package to create
phylogeny-aware visualizations of the HFHS diet study data.
Specifically, we walk through tree construction and generate
visualizations that reproduce key findings from previous studies,
including decreased relative abundance of but an increase in the
relative abundance of (Jo et al. 2021, Yang et al. 2024).

## Why phylogeny-aware plots?

Static stacked bar plots often hide changes that occur at specific
taxonomic levels. Phylogeny-aware plots maintain hierarchical structure,
so zooming to a level (for example, class Bacteroidia or family S24-7)
lets you see whether shifts are local to a clade or spread across
multiple branches. This supports exploratory analysis and helps decide
which levels warrant downstream testing.

## Setup

The following code block loads the packages required for this analysis.

``` r

library(phylobar)   # phylogeny-aware bar plots
library(dplyr)      # data wrangling
library(ape)        # phylogenetic tree manipulation
library(stringr)    # string operations
library(zen4R)      # download from zenodo link
```

Next, we load the HFHS dataset, which is stored in the LUPINE package as
`HFHSdata`. This dataset is stored as a list and contains OTU tables for
mice fed either a normal diet or an HFHS diet, along with the
corresponding taxonomic annotations.

``` r

download_zenodo("10.5281/zenodo.18791960", tempdir())
#> [zen4R][INFO] ZenodoRecord - Download in sequential mode
#> [zen4R][INFO] ZenodoRecord - Will download 1 file from record '18791960' (doi: '10.5281/zenodo.18791960') - total size: 76.2 KiB
#> [zen4R][INFO] Downloading file 'HFHS-data.rds' - size: 76.2 KiB
#> [zen4R][INFO] File downloaded at '/private/var/folders/mt/w0x_hxms30qch78n14v9zz0h0000gn/T/Rtmpy1jRY8'.
#> [zen4R][INFO] ZenodoRecord - Verifying file integrity...
#> [zen4R][INFO] File 'HFHS-data.rds': integrity verified (md5sum: 3266a55a3d0e01db0f0c99c7cb7a8e06)
#> [zen4R][INFO] ZenodoRecord - End of download
HFHSdata <- readRDS(str_c(tempdir(), "/HFHS-data.rds"))
normal <- HFHSdata$OTUdata_Normal
hfhs <- HFHSdata$OTUdata_HFHS
taxa <- HFHSdata$filtered_taxonomy
```

## Taxonomy-based tree construction

We begin by constructing a taxonomy-based tree from the taxonomic
annotations provided in the dataset. The `taxonomy_to_tree` function
from the phylobar package is used to convert the taxonomic data into a
phylogenetic tree structure.

``` r

tree <- taxa |>
    select(-X1, X1) |>
    mutate(
        across(
            everything(),
            ~if_else(str_ends(., "_"), NA, .)
        )
    ) |>
    taxonomy_to_tree()

checkValidPhylo(tree)
#> Starting checking the validity of tree...
#> Found number of tips: n = 212 
#> Found number of nodes: m = 80 
#> Done.
```

## Diet contrasts on day 7

To visualize the effect of diet on the mouse microbiome, we will focus
on samples collected on day 7. These samples are expected to show the
most pronounced differences between the normal diet and HFHS diet group
(Kodikara et al. 2025).

We begin by extracting the relevant slices from the 3D arrays
representing the OTU tables for both diet groups. We also rename the
samples to indicate their diet group (N for normal diet and H for HFHS
diet).

``` r

# Extract day 7 samples and rename
normal_day7 <- normal[,, 4]
rownames(normal_day7) <- str_replace(rownames(normal_day7), "M_", "N")
hfhs_day7 <- hfhs[,, 4]
rownames(hfhs_day7) <- str_replace(rownames(hfhs_day7), "M_", "H")
# Combine the two datasets
all_day7 <- rbind(normal_day7, hfhs_day7)

# Order samples and species using hierarchical clustering
comp_norm <- normal_day7
comp_norm <- comp_norm / rowSums(comp_norm)
sample_order1 <- hclust(dist(comp_norm))$order
species_order <- hclust(dist(t(comp_norm)))$order

# Now order the HFHS samples
comp_hfhs <- hfhs_day7
comp_hfhs <- comp_hfhs / rowSums(comp_hfhs)
sample_order2 <- hclust(dist(comp_hfhs))$order

# Combine the sample orders
sample_order <- c(rownames(normal_day7)[sample_order1], rownames(hfhs_day7)[sample_order2])

# Reorder the combined data
x <- all_day7[sample_order, species_order]
```

Finally we convert to compositions and visualize.

``` r

comp <- x / rowSums(x)
colnames(comp) <- colnames(all_day7)[species_order]
phylobar(comp, tree, hclust_order = FALSE, sample_font_size = 9)
```

Lets select Bacteroidia class from the dropdown and see the effect of
diet on this class.

### Interpretation

![Phylobar plot of the mouse microbiome on day 7, comparing normal diet
(N) and HFHS diet (H). The plot highlights the relative abundances of
Bacteroidia class, which is known to decrease in response to an HFHS
diet.](Figures/fig1.png)

Phylobar plot of the mouse microbiome on day 7, comparing normal diet
(N) and HFHS diet (H). The plot highlights the relative abundances of
Bacteroidia class, which is known to decrease in response to an HFHS
diet.

As shown in Figure @ref(fig:fig1), the relative abundance of Bacteroidia
decreases in the HFHS diet group compared to the normal diet group,
consistent with previous findings (Jo et al. 2021, Yang et al. 2024).
Similarly, we can visualize the Firmicutes phylum to observe its
increase in relative abundance under the HFHS diet (Figure
@ref(fig:fig2)).

![Phylobar plot of the mouse microbiome on day 7, comparing normal diet
(N) and HFHS diet (H). The plot highlights the relative abundances of
Firmicutes phylum, which is known to increase in response to an HFHS
diet.](Figures/fig2.png)

Phylobar plot of the mouse microbiome on day 7, comparing normal diet
(N) and HFHS diet (H). The plot highlights the relative abundances of
Firmicutes phylum, which is known to increase in response to an HFHS
diet.

## Longitudinal trajectories within HFHS

To explore longitudinal patterns in the HFHS diet group, we visualized
samples collected on days 0, 1, 4, and 7. This allows us to track how
the microbiome composition changes over time in response to the HFHS
diet. Samples were deliberately sorted by mouse and time point, enabling
us to observe within-mouse changes across the study period. This is
important because the initial relative abundances can differ between
mice, and aligning samples by individual helps highlight true
longitudinal shifts.

``` r

y <- hfhs[sample_order2, species_order, ] |>
    aperm(c(1, 3, 2)) |>
    matrix(nrow = nrow(hfhs)*dim(hfhs)[3] , ncol = ncol(hfhs))

rownames(y) <- apply(
    expand.grid(c(0, 1, 4, 7), str_replace(rownames(hfhs),"M_","")[sample_order2]),
    1, paste, collapse = "-"
)
```

We now convert to compositions and visualize the longitudinal data for
the HFHS diet group.

``` r

comp <- y / rowSums(y)
colnames(comp) <- colnames(hfhs)[species_order]
phylobar(comp, tree, hclust_order = FALSE, sample_font_size = 9)
```

``` r


# p <- phylobar(
#     comp, tree, hclust_order = FALSE,
#     width = 1200, height = 800, sample_font_size = 8
# )
# htmlwidgets::saveWidget(
#     p, "interactive_scatter.html", selfcontained = TRUE
# )
```

### Interpretations

Even though the variability between mice is high, we can observe a
general decrease in the relative abundance of the S24-7 family (within
the Bacteroidetes phylum) over time in response to the HFHS diet (Figure
@ref(fig:fig3)). This pattern aligns with previous findings indicating
that members of the S24-7 family are sensitive to dietary changes and
tend to decrease in abundance under high-fat, high-sugar dietary
conditions (Kodikara et al. 2025).

![Phylobar plot of the mouse microbiome for HFHS diet. The plot
highlights the relative abundances of S24-7 family that belong to the
Bacteroidetes phylum, which is known to decrease over time in response
to an HFHS diet.](Figures/fig3.png)

Phylobar plot of the mouse microbiome for HFHS diet. The plot highlights
the relative abundances of S24-7 family that belong to the Bacteroidetes
phylum, which is known to decrease over time in response to an HFHS
diet.

### Alternative Sorting

We can create the same phylobar visualization first sorting by time and
then by mouse. The variation in abundance within continguous blocks
gives a sense of mouse-to-mouse variability, while changes in area from
the left to the right side of the stacked bars gives a sense of the
change over time within these HFHS diet mice.

``` r

time <- as.numeric(sub("-.*", "", rownames(comp)))
mouse <- factor(
    sub(".*-", "", rownames(comp)),
    levels = unique(sub(".*-", "", rownames(comp)))
)
comp_sort <- comp[order(time, mouse), ]
phylobar(comp_sort, tree, hclust_order = FALSE, sample_font_size = 9)
```

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
#> [1] zen4R_0.10.5     stringr_1.6.0    ape_5.8-1        dplyr_1.2.1     
#> [5] phylobar_0.99.12 BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] utf8_1.2.6          sass_0.4.10         generics_0.1.4     
#>  [4] xml2_1.5.2          stringi_1.8.7       lattice_0.22-9     
#>  [7] digest_0.6.39       magrittr_2.0.5      evaluate_1.0.5     
#> [10] grid_4.6.0          bookdown_0.46       fastmap_1.2.0      
#> [13] plyr_1.8.9          jsonlite_2.0.0      Matrix_1.7-5       
#> [16] BiocManager_1.30.27 httr_1.4.8          purrr_1.2.2        
#> [19] XML_3.99-0.23       codetools_0.2-20    textshaping_1.0.5  
#> [22] jquerylib_0.1.4     cli_3.6.6           rlang_1.2.0        
#> [25] withr_3.0.2         cachem_1.1.0        yaml_2.3.12        
#> [28] otel_0.2.0          tools_4.6.0         parallel_4.6.0     
#> [31] fastmatch_1.1-8     curl_7.1.0          png_0.1-9          
#> [34] vctrs_0.7.3         R6_2.6.1            lifecycle_1.0.5    
#> [37] fs_2.1.0            htmlwidgets_1.6.4   ragg_1.5.2         
#> [40] cluster_2.1.8.2     pkgconfig_2.0.3     desc_1.4.3         
#> [43] pkgdown_2.2.0       bslib_0.10.0        pillar_1.11.1      
#> [46] glue_1.8.1          phangorn_2.12.1     Rcpp_1.1.1-1.1     
#> [49] systemfonts_1.3.2   xfun_0.57           tibble_3.3.1       
#> [52] tidyselect_1.2.1    keyring_1.4.1       knitr_1.51         
#> [55] htmltools_0.5.9     nlme_3.1-169        igraph_2.3.1       
#> [58] rmarkdown_2.31      compiler_4.6.0      quadprog_1.5-8
```
