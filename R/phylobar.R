#' Tree + Stacked Bars
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
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
    elementId = NULL,
    rel_width = 0.4,
    rel_space = 10
) {
    if (is.null(palette)) {
        palette <- c(
            "#9c7bbaff", "#6eb8acff", "#ce7b7bff",
            "#7b9cc4ff", "#c47ba0ff", "#e1d07eff"
        )
    }

    inputs <- phylobar_data(x, tree, hclust_order)
    opts <- list(
        rel_width = rel_width,
        rel_space = rel_space,
        sample_font_size  = sample_font_size,
        sample_label_margin = sample_label_margin,
        sample_label_space = sample_label_space,
        sample_magnify = sample_magnify,
        sample_show_all  = sample_show_all
    )

    htmlwidgets::createWidget(
        name = 'phylobar',
        x = list(
            tree_data = inputs$tree_data,
            labels = inputs$labels,
            palette = palette,
            opts = opts
        ),
        width = width,
        height = height,
        package = 'phylobar',
        elementId = elementId
    )
}