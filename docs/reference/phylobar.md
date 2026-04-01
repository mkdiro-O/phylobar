# Tree + Stacked Bars

`phylobar` is a visualization package that makes it possible to
construct a stacked barplot by interactively "painting" an associated
tree. This is an alternative to defining a color palette using a fixed
taxonomic resolution. It also helps avoid the issue of grouping all rare
taxa into a color for "other" since species can be chosen selectively,
we can paint a few rare taxa but not the rest.

## Usage

``` r
phylobar(
  x,
  tree,
  hclust_order = TRUE,
  palette = NULL,
  width = NULL,
  height = NULL,
  sample_font_size = 8,
  sample_label_margin = 10,
  sample_label_space = 50,
  sample_magnify = 1.5,
  sample_show_all = TRUE,
  element_id = NULL,
  rel_width = 0.4,
  rel_height = 0.85,
  rel_space = 10,
  legend_mode = TRUE,
  legend_x_start = 5,
  legend_spacing = 16
)
```

## Arguments

- x:

  A matrix of abundances. Samples along rows, features along columns.

- tree:

  An object of class phylo, representing the tree structure.

- hclust_order:

  Logical; if TRUE, reorder rows/columns by hierarchical clustering.

- palette:

  Character vector of colors for stacked bars. If NULL, uses default
  palette: c("#9c7bbaff", "#6eb8acff", "#ce7b7bff", "#7b9cc4ff",
  "#c47ba0ff", "#e1d07eff").

- width:

  Width of the widget in pixels. If NULL, uses window default.

- height:

  Height of the widget in pixels. If NULL, uses window default.

- sample_font_size:

  Font size for sample labels (integer).

- sample_label_margin:

  Margin between sample labels and bars in pixels.

- sample_label_space:

  Space allocated for sample labels in pixels.

- sample_magnify:

  Magnification factor for hovered sample labels.

- sample_show_all:

  Logical; if TRUE, show all sample labels.

- element_id:

  Optional HTML element ID to attach the widget to.

- rel_width:

  Width of the tree panel relative to the overall visualization.
  Defaults to 0.4.

- rel_height:

  Relative height of the tree in the overall visualization. Defaults to
  0.85. Adjust this if you need more/less space for the legend.

- rel_space:

  Space between tree and barplot panels in pixels.

- legend_mode:

  Logical; if TRUE (default), display labels for the painted subtrees in
  a legend near the bottom of the tree. If FALSE, include the labels
  within the tree itself.

- legend_x_start:

  Horizontal starting position (in pixels) for the legend. Defaults to
  4.

- legend_spacing:

  Vertical spacing (in pixels) between legend entries.

## Value

An htmlwidget visualization attached to the element element_id on the
output HTML page.

## Examples

``` r
library(ape)
tree <- rtree(5)
x <- matrix(rpois(15, 1), ncol = 5)
phylobar(x, tree)

{"x":{"tree_data":{"name":"node1","value":[[3],[9],[5]],"summary":17,"children":[{"name":"node2","value":[[1],[2],[3]],"summary":6,"children":[{"name":"t1","value":[[1],[1],[0]],"summary":2,"children":null},{"name":"t4","value":[[0],[1],[3]],"summary":4,"children":null}]},{"name":"node3","value":[[2],[7],[2]],"summary":11,"children":[{"name":"t3","value":[[1],[0],[1]],"summary":2,"children":null},{"name":"node4","value":[[1],[7],[1]],"summary":9,"children":[{"name":"t2","value":[[0],[5],[0]],"summary":5,"children":null},{"name":"t5","value":[[1],[2],[1]],"summary":4,"children":null}]}]}]},"labels":["sample_2","sample_1","sample_3"],"palette":["#9c7bbaff","#6eb8acff","#ce7b7bff","#7b9cc4ff","#c47ba0ff","#e1d07eff"],"opts":{"legend_mode":true,"legend_spacing":16,"legend_x_start":5,"rel_height":0.85,"rel_space":10,"rel_width":0.4,"sample_font_size":8,"sample_label_margin":10,"sample_label_space":50,"sample_magnify":1.5,"sample_show_all":true}},"evals":[],"jsHooks":[]}
```
