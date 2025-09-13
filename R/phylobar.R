#' Check that the input has the necessary labels
#'
#' The internal data structures in phylobar expect the tree to have node labels
#' and for the input sample composition matrix to have row and column names.
#' This helper provides defaults in case they are not present in the original
#' input.
#' @noRd
check_inputs <- function(x, tree) {
    if (is.null(rownames(x))) {
        rownames(x) <- paste0("sample_", seq_len(nrow(x)))
    }
    if (is.null(colnames(x))) {
        colnames(x) <- tree$tip.label
    }
    if (is.null(tree$node.label)) {
        tree$node.label <- paste0("node", seq_len(tree$Nnode))
    }
    list(x = x, tree = tree)
}

#' Tree + Stacked Bars
#'
#' `phylobar` is visualization package that makes it possible to construct a
#' stacked barplot by interactively "painting" an associated tree. This is an
#' alternative to defining a color palette using a fixed taxonomic resolution.
#' It also helps avoid the issue of grouping all rare taxa into a color for
#' "other" since species can be chosen selectively, we can paint a few rare taxa
#' but not the rest.
#'
#' @param tree A n object of class phylo, representing the tree structure.
#' @param x A matrix of abundances. Samples along rows, features along columns.
#' @param hclust_order Logical; if TRUE, reorder rows/columns by hierarchical
#'      clustering.
#' @param palette Character vector of colors for stacked bars. If NULL, uses
#'   default palette: c( "#9c7bbaff", "#6eb8acff", "#ce7b7bff", "#7b9cc4ff",
#'    "#c47ba0ff", "#e1d07eff").
#' @param width Width of the widget in pixels. If NULL, uses window default.
#' @param height Height of the widget in pixels. If NULL, uses window default.
#' @param sample_font_size Font size for sample labels (integer).
#' @param sample_label_margin Margin between sample labels and bars in pixels.
#' @param sample_label_space Space allocated for sample labels in pixels.
#' @param sample_magnify Magnification factor for hovered sample labels.
#' @param sample_show_all Logical; if TRUE, show all sample labels.
#' @param element_id Optional HTML element ID to attach the widget to.
#' @param rel_width Width of the tree panel relative to the overall
#'    visualization. Defaults to 0.4.
#' @param rel_space Space between tree and barplot panels in pixels.
#' @return An htmlwidget visualization attached to the element element_id on
#     the output HTML page.
#' @importFrom htmlwidgets createWidget
#' @examples
#' library(ape)
#' set.seed(1)
#' tree <- rtree(5)
#' x <- matrix(rpois(15), ncol = 5)
#' phylobar(x, tree)
#' @export
phylobar <- function(
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
    rel_space = 10
) {
    if (is.null(palette)) {
        palette <- c(
            "#9c7bbaff", "#6eb8acff", "#ce7b7bff",
            "#7b9cc4ff", "#c47ba0ff", "#e1d07eff"
        )
    }
    checked <- check_inputs(x, tree)

    inputs <- phylobar_data(checked$x, checked$tree, hclust_order)
    opts <- list(
        rel_width = rel_width,
        rel_space = rel_space,
        sample_font_size  = sample_font_size,
        sample_label_margin = sample_label_margin,
        sample_label_space = sample_label_space,
        sample_magnify = sample_magnify,
        sample_show_all  = sample_show_all
    )

    createWidget(
        name = "phylobar",
        x = list(
            tree_data = inputs$tree_data,
            labels = inputs$labels,
            palette = palette,
            opts = opts
        ),
        width = width,
        height = height,
        package = "phylobar",
        elementId = element_id
    )
}