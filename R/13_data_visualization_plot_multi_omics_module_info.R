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
#'   molecule_pathway (pathway–molecule links), pathway_similarity (pathway–pathway
#'   biotext-embedding cosine similarity).
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
#' @param metabolite_colors `character(2)`. A length-2 vector of colours for
#'   metabolite nodes. When `diff_metric` values are present, the two colours
#'   define the diverging gradient (first = most negative, second = most
#'   positive; midpoint `"#F2F2F2"`). When all `diff_metric` values are `NA`
#'   (e.g. no differential metric was supplied), the first colour is used as a
#'   solid fill for all metabolite nodes. Default `c("#71b7ed", "#f57c6e")`.
#' @param llm_text `logical(1)`. When `TRUE` and a `llm_module_name` column is
#'   present in `merge_result$functional_module_result`, the plot title is
#'   formatted as `"<module_id>: <llm_module_name>"`. Ignored when `title` is
#'   supplied. Default `FALSE`.
#'
#' @return A `ggplot` / `ggraph` object.
#'
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   scale_edge_colour_manual theme_graph
#' @importFrom ggplot2 aes scale_colour_manual labs theme element_text
#'   guides guide_legend scale_fill_gradient2
#' @importFrom ggnewscale new_scale_fill
#' @importFrom dplyr filter select mutate left_join bind_rows distinct case_when
#' @importFrom tidygraph tbl_graph activate as_tibble
#' @importFrom purrr map_chr map_dbl
#'
#' @export
plot_multi_omics_module_info <- function(
    merge_result,
    module_id,
    node_colors = c(
      "gene_Transcriptome" = "#fae69e",
      "gene_Proteome"      = "#f2b56f",
      "gene_T_and_P"       = "#b8aeeb",
      "pathway"            = "#bcd59b"
    ),
    node_shapes = c(
      "gene_Transcriptome" = 21, # circle (fill-able)
      "gene_Proteome"      = 21,
      "gene_T_and_P"       = 21,
      "metabolite"         = 24, # triangle up
      "pathway"            = 22  # square
    ),
    edge_colors = c(
      "TF-target" = "#E15759",
      "PPI" = "#76B7B2",
      "Reaction" = "#B07AA1",
      "molecule_pathway" = "#F1CE63",
      "pathway_similarity" = "#439222",
      "diffusion_similarity" = "grey"
    ),
    node_size = 5,
    label_size = 3,
    show_rwr_edge = FALSE,
    show_labels = TRUE,
    title = NULL,
    metabolite_colors = c("#71b7ed", "#f57c6e"),
    llm_text = FALSE
) {

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

  # Derive per-node type detail for gene nodes based on dt_src in node_info
  plot_nodes <- plot_nodes |>
    dplyr::mutate(
      .dt_src = purrr::map_chr(node_info, ~ {
        v <- .x[["dt_src"]]
        if (is.null(v) || length(v) == 0 || all(is.na(v))) NA_character_
        else as.character(v[[1]])
      }),
      node_type_detail = dplyr::case_when(
        node_type == "metabolite"                                    ~ "metabolite",
        node_type == "pathway"                                       ~ "pathway",
        node_type == "gene" & !is.na(.dt_src) &
          grepl("T", .dt_src) & grepl("P", .dt_src)                ~ "gene_T_and_P",
        node_type == "gene" & !is.na(.dt_src) & .dt_src == "P"     ~ "gene_Proteome",
        node_type == "gene"                                         ~ "gene_Transcriptome"
      ),
      diff_metric = purrr::map_dbl(node_info, ~ {
        v <- .x[["diff_metric"]]
        if (is.null(v) || length(v) == 0 || all(is.na(v))) NA_real_
        else as.numeric(v[[1]])
      })
    )

  # Compute symmetric gradient limits so sign of diff_metric is preserved even
  # when there is only one metabolite (degenerate range would otherwise map the
  # single value to the midpoint colour regardless of sign).
  met_diff_vals <- plot_nodes$diff_metric[plot_nodes$node_type == "metabolite"]
  met_diff_vals <- met_diff_vals[!is.na(met_diff_vals)]
  has_diff_metric <- length(met_diff_vals) > 0
  met_limits <- if (has_diff_metric) {
    max_abs <- max(abs(met_diff_vals), na.rm = TRUE)
    if (is.finite(max_abs) && max_abs > 0) c(-max_abs, max_abs) else NULL
  } else {
    NULL
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
        to   = gd_nodes$node_id[to]
      )
  }

  all_edges <- gd_edges_raw |>
    dplyr::filter(from %in% node_ids, to %in% node_ids) |>
    dplyr::mutate(
      edge_category = dplyr::if_else(
        edge_type == "diffusion_similarity", "computational", "knowledge"
      ),
      display_weight = dplyr::if_else(!is.na(weight), weight, diff_weight)
    ) |>
    dplyr::select(from, to, diff_weight, edge_type, weight, display_weight, edge_category) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  # Build knowledge-only graph and compute FR layout
  knowledge_edges <- all_edges |>
    dplyr::filter(edge_category == "knowledge")

  layout_edges <- if (nrow(knowledge_edges) == 0) {
    tibble::tibble(
      from = character(0), to = character(0),
      diff_weight = numeric(0), edge_type = character(0),
      weight = numeric(0), display_weight = numeric(0),
      edge_category = character(0)
    )
  } else {
    knowledge_edges
  }

  g_layout <- tidygraph::tbl_graph(
    nodes = plot_nodes |> dplyr::distinct(node_id, .keep_all = TRUE),
    edges = layout_edges,
    directed = FALSE,
    node_key = "node_id"
  )

  # Compute FR layout positions on the knowledge-only graph
  layout_coords <- ggraph::create_layout(g_layout, layout = "fr")
  fixed_xy <- layout_coords[, c("x", "y")]  # save positions keyed by row order

  # Build final graph (knowledge + optional RWR edges)
  plot_edges <- if (show_rwr_edge) all_edges else knowledge_edges

  if (nrow(plot_edges) == 0) {
    message(sprintf("No edges found for module '%s'. Plotting isolated nodes.", module_id))
    plot_edges <- tibble::tibble(
      from = character(0), to = character(0),
      diff_weight = numeric(0), edge_type = character(0),
      weight = numeric(0), display_weight = numeric(0),
      edge_category = character(0)
    )
  }

  g_final <- tidygraph::tbl_graph(
    nodes = plot_nodes |> dplyr::distinct(node_id, .keep_all = TRUE),
    edges = plot_edges,
    directed = FALSE,
    node_key = "node_id"
  )

  # Inject fixed coordinates as a manual layout
  final_layout <- ggraph::create_layout(g_final, layout = "manual",
                                        x = fixed_xy$x, y = fixed_xy$y)

  n_genes <- sum(plot_nodes$node_type == "gene",       na.rm = TRUE)
  n_mets  <- sum(plot_nodes$node_type == "metabolite", na.rm = TRUE)
  n_paths <- sum(plot_nodes$node_type == "pathway",    na.rm = TRUE)
  n_ke    <- dplyr::filter(plot_edges, edge_category == "knowledge")      |> nrow()
  n_de    <- dplyr::filter(plot_edges, edge_category == "computational")  |> nrow()

  if (!is.null(title)) {
    plot_title <- title
  } else if (llm_text) {
    fmr      <- merge_result$functional_module_result
    llm_name <- if (!is.null(fmr) && "llm_module_name" %in% colnames(fmr)) {
      fmr$llm_module_name[fmr$module == module_id][1]
    } else {
      NA_character_
    }
    plot_title <- if (!is.na(llm_name) && nzchar(llm_name)) {
      paste0(module_id, ": ", llm_name)
    } else {
      module_id
    }
  } else {
    plot_title <- module_id
  }

  plot_sub <- sprintf(
    "%d nodes  (%d genes  \u00b7  %d metabolites  \u00b7  %d pathways)   |   %d knowledge  \u00b7  %d diffusion edges",
    nrow(plot_nodes), n_genes, n_mets, n_paths, n_ke, n_de
  )

  edge_lty <- c(
    "TF-target"           = "solid",
    "PPI"                 = "solid",
    "Reaction"            = "solid",
    "molecule_pathway"    = "solid",
    "pathway_similarity"  = "solid",
    "diffusion_similarity" = "dashed"
  )

  # Internal node_colors with a neutral placeholder for metabolite so the first
  # layer can render all nodes; metabolites are overdrawn by the second layer.
  node_colors_all <- c(node_colors, "metabolite" = "#F2F2F2")

  # Fill values for the shape-legend override, in node_shapes order.
  used_types <- names(node_shapes)[
    names(node_shapes) %in% unique(plot_nodes$node_type_detail)
  ]
  fill_for_legend <- vapply(used_types, function(nt) {
    if (nt == "metabolite") "#F2F2F2" else node_colors[[nt]] %||% "grey50"
  }, character(1))

  p <- ggraph::ggraph(final_layout) +

    ggraph::geom_edge_link(
      ggplot2::aes(
        colour    = edge_type,
        linetype  = edge_type,
        edge_width = display_weight
      ),
      alpha = 0.6
    ) +
    ggraph::scale_edge_width(range = c(0.4, 0.8)) +
    ggraph::scale_edge_colour_manual(name = "Edge type", values = edge_colors) +
    ggraph::scale_edge_linetype_manual(name = "Edge type", values = edge_lty) +

    # All nodes: placeholder fill (metabolites overdrawn below); drives shape legend.
    ggraph::geom_node_point(
      ggplot2::aes(fill = node_type_detail, shape = node_type_detail),
      size   = node_size,
      colour = "black",
      stroke = 0.5
    ) +
    ggplot2::scale_fill_manual(
      name   = "Node type",
      values = node_colors_all,
      guide  = "none"          # legend driven by shape scale instead
    ) +
    ggplot2::scale_shape_manual(
      name  = "Node type",
      values = node_shapes,
      guide = ggplot2::guide_legend(
        override.aes = list(fill = fill_for_legend, colour = "black", stroke = 0.5)
      )
    ) +

    # Metabolite nodes: gradient if diff_metric exists, solid blue otherwise.
    {if (has_diff_metric) {
        list(
          ggnewscale::new_scale_fill(),
          ggraph::geom_node_point(
            data   = \(x) dplyr::filter(x, .data$node_type == "metabolite"),
            ggplot2::aes(fill = .data$diff_metric),
            shape  = 24,
            size   = node_size,
            colour = "black",
            stroke = 0.5
          ),
          ggplot2::scale_fill_gradient2(
            name     = "diff_metric",
            low      = metabolite_colors[1],
            mid      = "#F2F2F2",
            high     = metabolite_colors[2],
            midpoint = 0,
            limits   = met_limits,
            na.value = "grey80"
          )
        )
      } else {
        ggraph::geom_node_point(
          data   = \(x) dplyr::filter(x, .data$node_type == "metabolite"),
          fill   = "#71b7ed",
          shape  = 24,
          size   = node_size,
          colour = "black",
          stroke = 0.5
        )
      }
    } +

    {
      if (show_labels)
        ggraph::geom_node_text(
          ggplot2::aes(label = label),
          size        = label_size,
          repel       = TRUE,
          colour      = "grey15",
          bg.colour   = "white",
          bg.r        = 0.12,
          max.overlaps = 20
        )
    } +

    ggraph::theme_graph() +
    ggplot2::labs(title = plot_title, subtitle = plot_sub) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"),
      legend.title  = ggplot2::element_text(size = 9, face = "bold"),
      legend.position = "right"
    )

  p
}
