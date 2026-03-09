# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/multi_omics_modules.rda")
# merge_result <- multi_omics_modules
# module_id <- "Functional_module_92"
# {
#   merge_result = multi_omics_modules
#   module_id = "Functional_module_92"
#   node_colors = c(
#     "gene" = "#4E79A7",
#     "metabolite" = "#F28E2B",
#     "pathway" = "#59A14F"
#   )
#   node_shapes = c(
#     "gene" = 21, # circle (fill-able)
#     "metabolite" = 24, # triangle up
#     "pathway" = 22 # square
#   )
#   edge_colors = c(
#     "TF-target" = "#E15759",
#     "PPI" = "#76B7B2",
#     "Reaction" = "#B07AA1",
#     "molecule_pathway" = "#F1CE63",
#     "diffusion_similarity" = "grey"
#   )
#   node_size = 5
#   label_size = 3
#   show_rwr_edge = FALSE
#   show_labels = TRUE
#   title = NULL
# }
#
# p1 <- plot_multi_omics_module_info(merge_result = multi_omics_modules,
#                                    module_id = c("Functional_module_89"),
#                                    node_colors = node_colors,
#                                    node_shapes = node_shapes,
#                                    edge_colors = edge_colors,
#                                    node_size = node_size,
#                                    label_size = label_size,
#                                    show_rwr_edge = show_rwr_edge,
#                                    show_labels = show_labels,
#                                    title = title
#                                    )
# p1

#' Plot Multi-Omics Functional Module Information
#'
#' Visualises a single functional module as a network graph. Nodes represent
#' genes, metabolites, or enriched pathways. Edges are drawn in two categories:
#'
#' * **Knowledge edges** (solid lines, coloured by type):
#'   TF-target, PPI, Reaction (metabolite_reaction + enzyme_metabolite),
#'   molecule_pathway (pathway–molecule links).
#'
#' * **Computational edges** (dashed grey lines):
#'   Diffusion-based similarity edges with no matching knowledge-layer
#'   connection (edge_type == "diffusion_similarity").
#'
#' All edge information — including `edge_type`, `weight` (knowledge layer
#' weight), and `diff_weight` (cosine similarity from the diffusion matrix) —
#' is read directly from the `graph_data` tidygraph object produced by
#' [merge_multi_omics_nodes()]. No `mnet_obj` is required.
#'
#' @param merge_result A `list`. Direct output of [merge_multi_omics_nodes()],
#'   containing elements `graph_data`, `functional_module_result`, and
#'   `result_with_module`. The `graph_data` edge table must include `edge_type`,
#'   `weight`, and `diff_weight` columns (produced by the current version of
#'   [merge_multi_omics_nodes()]).
#' @param module_id `character(1)`. Module name to plot, e.g. `"Functional_module_1"`.
#' @param node_colors Named character vector of fill colours for each node type
#'   (`"gene"`, `"metabolite"`, `"pathway"`).
#' @param node_shapes Named integer vector of point shapes for each node type.
#'   Use fill-able shapes (21 = circle, 22 = square, 24 = triangle up).
#' @param edge_colors Named character vector of colours for each edge type.
#' @param node_size `numeric(1)`. Base node size. Default `5`.
#' @param label_size `numeric(1)`. Text label size (pt). Default `3`.
#' @param show_rwr_edge `logical(1)`. Whether to display diffusion-similarity
#'   edges in addition to knowledge edges. Default `FALSE`.
#' @param show_labels `logical(1)`. Whether to display node labels. Default `TRUE`.
#' @param title `character(1)` or `NULL`. Plot title. When `NULL`
#'   (default) the module name is used.
#'
#' @return A `ggplot` / `ggraph` object.
#'
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_colour_manual theme_graph
#' @importFrom ggplot2 aes scale_colour_manual labs theme element_text
#'   guides guide_legend
#' @importFrom dplyr filter select mutate left_join bind_rows distinct case_when
#' @importFrom tidygraph tbl_graph activate as_tibble
#' @importFrom purrr map_chr
#'
#' @export
plot_multi_omics_module_info <- function(
    merge_result,
    module_id,
    node_colors = c(
      "gene" = "#4E79A7",
      "metabolite" = "#F28E2B",
      "pathway" = "#59A14F"
    ),
    node_shapes = c(
      "gene" = 21, # circle (fill-able)
      "metabolite" = 24, # triangle up
      "pathway" = 22 # square
    ),
    edge_colors = c(
      "TF-target" = "#E15759",
      "PPI" = "#76B7B2",
      "Reaction" = "#B07AA1",
      "molecule_pathway" = "#F1CE63",
      "diffusion_similarity" = "grey"
    ),
    node_size = 5,
    label_size = 3,
    show_rwr_edge = FALSE,
    show_labels = TRUE,
    title = NULL
) {

  # Resolve module ID and extract nodes
  result_with_module <- merge_result$result_with_module

  plot_nodes <- result_with_module |>
    dplyr::filter(module %in% module_id) |>
    dplyr::mutate(
      label = dplyr::if_else(
        node_type == "metabolite",
        purrr::map_chr(node_info, ~ .x[["cpd_name"]] %||% NA_character_),
        node_id
      )
    )

  if (nrow(plot_nodes) == 0) {
    stop(sprintf("Module '%s' not found.", module_id))
  }

  node_ids <- unique(plot_nodes$node_id)

  gd_nodes <- merge_result$graph_data |>
    tidygraph::activate("nodes") |>
    tibble::as_tibble()

  gd_edges_raw <- merge_result$graph_data |>
    tidygraph::activate("edges") |>
    tibble::as_tibble()

  if (nrow(gd_edges_raw) > 0 && is.numeric(gd_edges_raw$from)) {
    gd_edges_raw <- gd_edges_raw |>
      dplyr::mutate(
        from = gd_nodes$node_id[from],
        to = gd_nodes$node_id[to]
      )
  }

  all_edges <- gd_edges_raw |>
    dplyr::filter(from %in% node_ids, to %in% node_ids) |>
    dplyr::mutate(
      edge_category = dplyr::if_else(
        edge_type == "diffusion_similarity", "computational", "knowledge"
      ),
      # Use knowledge weight when available, fall back to diffusion similarity
      display_weight = dplyr::if_else(!is.na(weight), weight, diff_weight)
    ) |>
    dplyr::select(from, to, diff_weight, edge_type, weight, display_weight, edge_category) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  # Optionally drop diffusion-only edges
  if (!show_rwr_edge) {
    all_edges <- all_edges |>
      dplyr::filter(edge_category == "knowledge")
  }

  # Build tidygraph object
  if (nrow(all_edges) == 0) {
    message(sprintf("No edges found for module '%s'. Plotting isolated nodes.", module_id))
    all_edges <- tibble::tibble(
      from = character(0),
      to = character(0),
      diff_weight = numeric(0),
      edge_type = character(0),
      weight = numeric(0),
      display_weight = numeric(0),
      edge_category = character(0)
    )
  }

  g <- tidygraph::tbl_graph(
    nodes = plot_nodes |> dplyr::distinct(node_id, .keep_all = TRUE),
    edges = all_edges,
    directed = FALSE,
    node_key = "node_id"
  )

  edge_lty <- c(
    "TF-target" = "solid",
    "PPI" = "solid",
    "Reaction" = "solid",
    "molecule_pathway" = "solid",
    "diffusion_similarity" = "dashed"
  )

  # Count stats for subtitle
  n_genes <- sum(plot_nodes$node_type == "gene", na.rm = TRUE)
  n_mets  <- sum(plot_nodes$node_type == "metabolite", na.rm = TRUE)
  n_paths <- sum(plot_nodes$node_type == "pathway", na.rm = TRUE)
  n_ke <- all_edges |> dplyr::filter(edge_category == "knowledge") |> nrow()
  n_de <- all_edges |> dplyr::filter(edge_category == "computational") |> nrow()

  plot_title <- if (!is.null(title)) title else module_id
  plot_sub <- sprintf(
    "%d nodes  (%d genes  \u00b7  %d metabolites  \u00b7  %d pathways)   |   %d knowledge  \u00b7  %d diffusion edges",
    nrow(plot_nodes), n_genes, n_mets, n_paths, n_ke, n_de
  )

  # Draw
  set.seed(42)

  p <- ggraph::ggraph(g, layout = "fr") +

    ggraph::geom_edge_link(
      ggplot2::aes(
        colour = edge_type,
        linetype = edge_type,
        edge_width = display_weight
      ),
      alpha = 0.6
    ) +
    ggraph::scale_edge_width(range = c(0.2, 1)) +
    ggraph::scale_edge_colour_manual(
      name = "Edge type",
      values = edge_colors
    ) +
    ggraph::scale_edge_linetype_manual(
      name = "Edge type",
      values = edge_lty
    ) +

    # Nodes
    ggraph::geom_node_point(
      ggplot2::aes(fill = node_type, shape = node_type),
      size = node_size,
      colour = "black",
      stroke = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name = "Node type",
      values = node_colors,
    ) +
    ggplot2::scale_shape_manual(
      name  = "Node type",
      values = node_shapes
    ) +

    # Labels
    {
      if (show_labels)
        ggraph::geom_node_text(
          ggplot2::aes(label = label),
          size = label_size,
          repel = TRUE,
          colour = "grey15",
          bg.colour = "white",
          bg.r = 0.12,
          max.overlaps = 20
        )
    } +

    # Theme
    ggraph::theme_graph() +
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_sub
    ) +
    ggplot2::theme(
      # plot.title    = ggplot2::element_text(face = "bold", size = 14),
      # plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"),
      # plot.caption  = ggplot2::element_text(size = 8, colour = "grey55",
      #                                       face = "italic"),
      # legend.title  = ggplot2::element_text(size = 9, face = "bold"),
      # legend.text   = ggplot2::element_text(size = 8),
      legend.position = "right"
    )

  p
}
