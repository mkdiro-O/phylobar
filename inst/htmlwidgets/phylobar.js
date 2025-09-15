HTMLWidgets.widget({

  name: 'phylobar',

  type: 'output',

  factory: function(el, width, height) {
    let tree_data, neighborhoods, feature_map, b_scale, rscale, stacks, labels,
      palette, opts;
    // initialize the HTML skeleton
    let svg = d3.select(el).append("svg")
      .attr("width", width)
      .attr("height", height)
    el.frozen = false

    svg.append("g").attr("id", "barplot");
    svg.append("g").attr("id", "tree");
    svg.append("g").attr("id", "tree_labels");
    svg.select("#tree").append("g").attr("id", "links");
    svg.select("#tree").append("g").attr("id", "nodes");
    svg.select("#barplot").append("g").attr("id", "sample_labels");
    svg.select("#barplot").append("g").attr("id", "bars");

    // some variables used throughout
    el.color_ix = 0;
    let color_sets = new Map();
    const link_gen = d3.linkVertical()
      .x(d => d.x)
      .y(d => d.y);

    // focus on the current widget element
    el.setAttribute("tabindex", "-1");
    el.focus();

    return {
      renderValue: function(x) {
        labels = x.labels;
        palette = x.palette;
        opts = x.opts;
        tree_data = x.tree_data;

        x_max = opts.rel_width * width
        let tree = phylobar.make_tree(tree_data, x_max, height);
        tree.each(d => { d._children = null});

        feature_map = phylobar.create_feature_map(tree);
        b_scale = phylobar.stack_scales(
          labels, phylobar.stack_data(tree.leaves(), labels),
          x_extent=[opts.rel_space + x_max, width],
          y_extent=[height - opts.sample_label_space, 0]
        );

        rscale = phylobar.radius_scale(tree);
        neighborhoods = d3.Delaunay.from(tree.descendants().map(d => [d.x, d.y + 10]));
        stacks = phylobar.stack_data(tree.leaves(), labels);

        phylobar.update_tree(el, tree, rscale, link_gen,  palette, color_sets);
        phylobar.update_stack(el, stacks, b_scale, labels, feature_map, palette, color_sets);
        phylobar.update_event_listeners(el, tree, neighborhoods, feature_map, palette, color_sets, true, x_max);
        phylobar.update_tree_labels(el, color_sets, feature_map, tree.descendants());
        phylobar.update_sample_labels(el, b_scale, labels, opts, x_max);
        phylobar.update_resolution(
          el, neighborhoods, tree, stacks, labels, palette, color_sets, 
          feature_map, opts, b_scale, rscale, link_gen, width, height, x_max
        )
      },

      resize: function(width, height) {
        if (!tree_data || !feature_map || !b_scale) return;
        let tree = phylobar.make_tree(tree_data, opts.rel_width * width, height);
        feature_map = phylobar.create_feature_map(tree);
        b_scale = phylobar.stack_scales(
          labels, phylobar.stack_data(tree.leaves(), labels),
          x_extent=[opts.rel_space + opts.rel_width * width, width],
          y_extent=[height, 0]
        );

        neighborhoods = d3.Delaunay.from(tree.descendants().map(d => [d.x, d.y + 10]));
        phylobar.update_tree(el, tree, rscale, link_gen, palette, color_sets);
        phylobar.update_stack(el, stacks, b_scale, labels, feature_map, palette, color_sets);
        phylobar.update_event_listeners(el, tree, neighborhoods, feature_map, palette, color_sets, false, x_max);
      }
    };
  }
});
