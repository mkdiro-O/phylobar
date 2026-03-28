# Getting Started

## Overview

phylobar is a visualization package for constructing stacked barplots by
interactively “painting” an associated tree. This is an alternative to
defining a color palette at a fixed taxonomic resolution. It also helps
avoid the issue of grouping all rare taxa into a color for “other” — you
can selectively highlight specific rare taxa while leaving others at
their background color.

The main function, [`phylobar()`](../reference/phylobar.md), takes two
inputs: a table of (potentially normalized) counts and a tree structure
which can come from either a phylogeny or a table of taxonomic
assignments. phylobar is designed for submission to
[Bioconductor](https://bioconductor.org/) and is compatible with popular
Bioconductor containers such as
[TreeSummarizedExperiment](https://bioconductor.org/packages/TreeSummarizedExperiment/)
and [phyloseq](https://bioconductor.org/packages/phyloseq/), so the
visualization works seamlessly with a wide range of existing datasets.
Several vignettes — including [Global Patterns](global_patterns.md),
[Diet Analysis](diet_analysis.md), [Taxonomies](taxonomies.md), and
[Atlas](atlas.md) — demonstrate working with data from these containers.
You can read more about the package in our preprint:

> Kuo, M., Lê Cao, K.-A., Kodikara, S., Mao, J., & Sankaran, K. (2025).
> phylobar: an R package for multiresolution compositional barplots in
> omics studies. <doi:10.1101/2025.11.05.686662>

## Installation

You can install the development version of phylobar from GitHub:

``` r

remotes::install_github("mkdiro-o/phylobar")
```

## Quick Start

The block below applies phylobar to a small random dataset. Click a node
to collapse a subtree, press the Control key to add a new color, and
press the Escape key to freeze the view.

``` r

library(ape)
library(phylobar)

tree <- rtree(20)
samples <- matrix(rpois(100 * 20, 1), 100, 20)
phylobar(samples, tree)
```

## Vignettes

The vignettes below walk through common workflows and fully developed
case studies. They are grouped into **workflow** vignettes (how to
prepare, style, and export phylobar plots) and **case-study** vignettes
(biological analyses that illustrate phylobar in practice).

### Workflow Vignettes

#### [Building Trees from Taxonomies](taxonomies.md)

Rather than interacting with abstract phylogenetic trees, it can be
helpful to choose stacked bar colors using taxonomic assignments. This
vignette gives a comprehensive guide to converting taxonomy tables into
valid `phylo` objects via
[`taxonomy_to_tree()`](../reference/taxonomy_to_tree.md). It covers
three common pitfalls — missing root nodes, missing assignments encoded
as strings instead of `NA`, and duplicated names across taxonomic levels
— and shows how to resolve each with
[`add_prefix()`](../reference/add_prefix.md) and
[`checkValidPhylo()`](https://rdrr.io/pkg/ape/man/checkValidPhylo.html).

#### [Customizing Style](customizing_files.md)

Phylobar supports several styling options, including custom color
palettes, widget dimensions, sample label font and margin settings,
tree-bar layout ratios (`rel_width`, `rel_height`, `rel_space`), legend
placement, and hierarchical reordering with `hclust_order`. This
vignette walks through all currently available parameters.

#### [Exporting Views](exporting.md)

While interactivity is useful for exploration, we often need to export a
specific static view to discuss with collaborators. This vignette shows
how to export phylobar visualizations to SVG using the SVG Crowbar
bookmarklet, edit the resulting layered SVG in Inkscape, and export to
PNG — preserving image quality for publication.

#### [Runtime Evaluation](runtime_evaluation.md)

This vignette benchmarks [`phylobar()`](../reference/phylobar.md) across
a grid of sample counts and taxa counts, with and without hierarchical
clustering. It shows that disabling `hclust_order` improves performance
at large sample sizes, and recommends using
[`subset_cluster()`](../reference/subset_cluster.md) to keep the number
of samples browser-friendly.

### Case-Study Vignettes

#### [HFHS Diet and the Mouse Microbiome](hfhs.md)

Analyzes 16S rRNA data from 47 mice on either normal or high-fat,
high-sugar (HFHS) diets across four timepoints. The vignette reproduces
known findings of decreased Bacteroidetes and increased Firmicutes under
HFHS and adds a longitudinal view tracking individual mice over time.

#### [Atlas 1006 — Tipping Points in the Human Gut](atlas.md)

Demonstrates phylobar on the HITChip Atlas 1006 human gut microbiome
dataset, highlighting “tipping element” taxa and the
Firmicutes/Bacteroidetes balance. It also shows how
[`subset_cluster()`](../reference/subset_cluster.md) can reduce visual
clutter by selecting representative samples.

#### [Diet Analysis with DESeq2 Normalization](diet_analysis.md)

Examines how diet shapes gut microbiome composition using the `dietswap`
dataset. It applies DESeq2 normalization via MicrobiomeStat, builds a
taxonomy tree, and visualizes the resulting abundances with phylobar.

#### [Global Patterns — phyloseq and TreeSummarizedExperiment](global_patterns.md)

Shows how to create phylobar plots from both `phyloseq` and
`TreeSummarizedExperiment` input objects using the Global Patterns
dataset, demonstrating that both container formats work interchangeably.

#### [COVID-19 Immunology — Beyond Microbiome Data](covid_immunology.md)

Applies phylobar to immune cell-type compositions from COVID-19 patient
data, using a manually defined cell-type hierarchy in place of a
microbial taxonomy. This vignette illustrates phylobar’s applicability
beyond microbiome studies.

## Function Reference

The full function reference is available at
[https://mkdiro-o.github.io/phylobar/reference/](https://mkdiro-o.github.io/phylobar/reference/index.html).

## Contact

You can reach us by creating an
[Issue](https://github.com/mkdiro-O/phylobar/issues) or emailing
<ksankaran@wisc.edu>. We appreciate your interest and will respond
promptly.

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
#> [1] BiocStyle_2.36.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39       desc_1.4.3          R6_2.6.1           
#>  [4] bookdown_0.46       fastmap_1.2.0       xfun_0.56          
#>  [7] cachem_1.1.0        knitr_1.51          htmltools_0.5.9    
#> [10] rmarkdown_2.30      lifecycle_1.0.5     cli_3.6.5          
#> [13] sass_0.4.10         pkgdown_2.2.0       textshaping_1.0.4  
#> [16] jquerylib_0.1.4     systemfonts_1.3.1   compiler_4.5.1     
#> [19] tools_4.5.1         ragg_1.5.0          bslib_0.10.0       
#> [22] evaluate_1.0.5      yaml_2.3.12         BiocManager_1.30.27
#> [25] otel_0.2.0          jsonlite_2.0.0      rlang_1.1.7        
#> [28] fs_1.6.6            htmlwidgets_1.6.4
```
