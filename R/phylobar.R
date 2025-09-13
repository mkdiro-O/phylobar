#" check that the input has the necessary labels
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

#" Tree + Stacked Bars
#"
#" <Add Description>
#"
#" @import htmlwidgets
#"
#" @export
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

    htmlwidgets::createWidget(
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