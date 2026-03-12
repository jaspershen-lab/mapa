# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/multi_omics_modules.rda")
# load("demo_data/demo_multi-omics/llm_interpreted_object.rda")
# module_id <- c("Functional_module_79", "Functional_module_89")

# Internal helper: grouped FR layout ─────────────────────────────────────────
# Places N module centroids evenly on a large circle, runs an independent FR
# sub-layout for each module's induced subgraph, normalises local coords to
# [-3, 3], then shifts each cluster to its centroid.
# Returns list(x, y) in the original node order of `graph`.
.make_grouped_layout <- function(graph, modules, module_circle_radius = 10, scale_factor = 1.2, seed = 42) {

  set.seed(seed)

  node_tbl <- graph |> tidygraph::activate("nodes") |> tibble::as_tibble()
  edge_tbl <- graph |> tidygraph::activate("edges") |> tibble::as_tibble()
  # edge_tbl$from / $to are 1-based row indices into node_tbl

  n_mod  <- length(modules)
  radius <- if (n_mod > 1) module_circle_radius else 0
  angles <- seq(0, 2 * pi, length.out = n_mod + 1)[seq_len(n_mod)]
  cx <- stats::setNames(radius * cos(angles), modules)
  cy <- stats::setNames(radius * sin(angles), modules)

  local_list <- vector("list", n_mod)

  for (i in seq_along(modules)) {
    mod <- modules[i]
    node_idx <- which(node_tbl$module == mod)
    n <- length(node_idx)

    if (n == 1L) {
      local_list[[i]] <- data.frame(row_idx = node_idx, lx = 0, ly = 0)
      next
    }

    # Induced edges: both endpoints inside this module
    keep <- edge_tbl$from %in% node_idx & edge_tbl$to %in% node_idx
    local_edges <- edge_tbl[keep, ]

    if (nrow(local_edges) == 0L) {
      # No internal edges: scatter on a small circle
      ang <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
      local_list[[i]] <- data.frame(
        row_idx = node_idx,
        lx = 1.5 * cos(ang),
        ly = 1.5 * sin(ang)
      )
      next
    }

    # Re-index from/to to local positions 1:n
    from_local <- match(local_edges$from, node_idx)
    to_local <- match(local_edges$to,   node_idx)

    ig <- igraph::graph_from_data_frame(
      d = data.frame(from = from_local, to = to_local),
      directed = FALSE,
      vertices = seq_len(n)
    )
    pos <- igraph::layout_with_fr(ig)

    local_list[[i]] <- data.frame(row_idx = node_idx, lx = pos[, 1L], ly = pos[, 2L])
  }

  # Normalise each module to [-3, 3] then shift to its centroid
  result_list <- vector("list", n_mod)

  for (i in seq_along(modules)) {
    mod <- modules[i]
    df <- local_list[[i]]
    lx <- df$lx
    ly <- df$ly

    n <- nrow(df)
    local_scale <- max(1, log2(n) * scale_factor)

    rx <- diff(range(lx))
    ry <- diff(range(ly))
    lx <- if (rx > 0) (lx - mean(range(lx))) / (rx / 2) * local_scale else rep(0, length(lx))
    ly <- if (ry > 0) (ly - mean(range(ly))) / (ry / 2) * local_scale else rep(0, length(ly))

    result_list[[i]] <- data.frame(
      row_idx = df$row_idx,
      x = lx + cx[mod],
      y = ly + cy[mod]
    )
  }

  result_df <- do.call(rbind, result_list)
  result_df <- result_df[order(result_df$row_idx), ]

  list(x = result_df$x, y = result_df$y)
}


#' Plot Multi-Omics Functional Module Similarity Network
#'
#' Creates a similarity network visualisation of multi-omics functional modules.
#'
#' @param object `list`. Direct output of [merge_multi_omics_nodes()] or
#'   [llm_interpret_module()], containing `graph_data`, `functional_module_result`, and `result_with_module`.
#' @param module_size_cutoff `numeric(1)` or `NULL`. Retain only modules whose
#'   `module_content_number` is **>=** this value. Mutually exclusive with
#'   `module_id`. When both are `NULL` all modules are shown. Default `NULL`.
#' @param module_id `character` vector or `NULL`. Retain only the listed module
#'   IDs (e.g. `c("Functional_module_1", "Functional_module_3")`). Mutually
#'   exclusive with `module_size_cutoff`. Default `NULL`.
#' @param node_shapes Named integer vector of point shapes for each node type.
#'   Use fill-able shapes (21 = circle, 22 = square, 24 = triangle up).
#' @param edge_colors Named character vector of colours for each edge type.
#' @param node_size `numeric(1)`. Base point size. Default `4`.
#' @param label_size `numeric(1)`. Module label text size (pt). Default `3`.
#' @param show_rwr_edge `logical(1)`. Whether to display diffusion-similarity
#'   edges in addition to knowledge edges. Default `TRUE`.
#' @param show_labels `logical(1)`. Whether to display module-level labels
#'   at the centroid of each module. Default `TRUE`.
#' @param use_llm_label `logical(1)`. When `TRUE` (default), module labels use
#'   the LLM-generated name (`llm_module_name`) when available, falling back to
#'   the module ID. When `FALSE`, always use the module ID as the label.
#' @param show_node_labels `logical(1)`. Whether to display individual node ID
#'   labels on each node. Default `FALSE`.
#' @param node_label_size `numeric(1)`. Text size for node labels (pt).
#'   Default `2.5`.
#' @param module_circle_radius `numeric(1)`. Radius of the large circle on which
#'   module centroids are placed. Increase to spread modules further apart;
#'   ignored when only one module is shown. Default `10`.
#' @param module_scale_factor `numeric(1)`. Controls how spread-out nodes are
#'   within each module. Larger values expand each module's internal layout;
#'   the spread also scales with `log2(node count)`. Default `1.2`.
#' @param outline_expand `numeric(1)`. Extra space (mm) added around each
#'   convex-hull outline. Default `5`.
#' @param seed `integer(1)`. Random seed for the Fruchterman-Reingold layout.
#'   Default `42`.
#'
#' @return A `ggplot` / `ggraph` object.
#'
#' @importFrom ggraph ggraph geom_edge_link geom_node_point
#'   scale_edge_colour_manual scale_edge_linetype_manual scale_edge_width
#'   theme_graph create_layout
#' @importFrom ggplot2 aes scale_fill_manual scale_shape_manual
#'   scale_colour_manual scale_size_identity theme unit ggplot geom_blank guides
#'   guide_legend
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter mutate if_else pull left_join group_by summarise
#'   first any_of select
#' @importFrom tidygraph activate tbl_graph
#' @importFrom tibble as_tibble
#' @importFrom stringr str_replace str_wrap
#'
#' @export
plot_multi_omics_similarity_network <- function(
    object,
    module_size_cutoff = NULL,
    module_id = NULL,
    node_shapes = c(
      "gene" = 21,
      "metabolite" = 24,
      "pathway" = 22
    ),
    edge_colors = c(
      "TF-target" = "#E15759",
      "PPI" = "#76B7B2",
      "Reaction" = "#B07AA1",
      "molecule_pathway" = "#F1CE63",
      "pathway_similarity" = "#439222",
      "diffusion_similarity" = "grey80"
    ),
    node_size = 4,
    label_size = 3,
    show_rwr_edge = TRUE,
    show_labels = TRUE,
    use_llm_label = TRUE,
    show_node_labels = FALSE,
    node_label_size = 2.5,
    module_circle_radius = 10,
    module_scale_factor = 1.2,
    outline_expand = 5,
    seed = 42
) {

  required_slots <- c("graph_data", "functional_module_result", "result_with_module")
  if (!is.list(object) || !all(required_slots %in% names(object))) {
    stop(
      "'object' must be the output of merge_multi_omics_nodes() or ",
      "llm_interpret_module(), containing: ",
      paste(required_slots, collapse = ", "), "."
    )
  }
  if (!is.null(module_size_cutoff) && !is.null(module_id)) {
    stop("Provide either 'module_size_cutoff' or 'module_id', not both.")
  }

  functional_module_result <- object$functional_module_result

  # Determine selected modules
  if (!is.null(module_size_cutoff)) {
    selected_modules <- functional_module_result |>
      dplyr::filter(module_content_number >= module_size_cutoff) |>
      dplyr::pull(module)

  } else if (!is.null(module_id)) {
    invalid <- setdiff(module_id, functional_module_result$module)
    if (length(invalid) > 0) {
      warning("Module IDs not found and will be skipped: ",
              paste(invalid, collapse = ", "))
    }
    selected_modules <- intersect(module_id, functional_module_result$module)
  } else {
    selected_modules <- unique(functional_module_result$module)
  }

  if (length(selected_modules) == 0) {
    warning("No modules remain after filtering.")
    return(ggplot2::ggplot() + ggplot2::geom_blank())
  }

  # Extract subnetwork from graph_data
  subgraph <- object$graph_data |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(module = stringr::str_replace(module, "^Module_", "Functional_module_")) |>
    dplyr::filter(module %in% selected_modules)

  if (igraph::gorder(subgraph) == 0) {
    warning("No nodes found for the selected modules in graph_data.")
    return(ggplot2::ggplot() + ggplot2::geom_blank())
  }

  # Optionally remove diffusion-similarity edges
  if (!show_rwr_edge) {
    subgraph <- subgraph |>
      tidygraph::activate("edges") |>
      dplyr::filter(edge_type != "diffusion_similarity")
  }

  # Add edge display weight (knowledge weight when available, else diff_weight)
  subgraph <- subgraph |>
    tidygraph::activate("edges") |>
    dplyr::mutate(display_weight = dplyr::if_else(!is.na(weight), weight, diff_weight))

  # Module label lookup: LLM name with module ID fallback
  has_llm <- "llm_module_name" %in% colnames(functional_module_result)

  module_label_df <- functional_module_result |>
    dplyr::filter(module %in% selected_modules) |>
    dplyr::select(module, dplyr::any_of("llm_module_name")) |>
    dplyr::mutate(
      module_label = if (use_llm_label && has_llm) {
        dplyr::if_else(
          !is.na(llm_module_name) & nchar(llm_module_name) > 0,
          llm_module_name,
          module
        )
      } else {
        module
      }
    )

  # Layout — grouped FR: each module is laid out independently then placed on a
  # large circle so clusters are spatially separated.
  # create_layout() is called first just to obtain the correct ggraph layout
  # object structure/class; its x/y are then overwritten with the grouped coords.
  set.seed(seed)
  lay <- ggraph::create_layout(subgraph, layout = "fr")
  grouped_coords <- .make_grouped_layout(subgraph, selected_modules,
                                         module_circle_radius = module_circle_radius,
                                         scale_factor = module_scale_factor,
                                         seed = seed)
  lay$x <- grouped_coords$x
  lay$y <- grouped_coords$y

  # Module centroids for labels
  centroids <- lay |>
    dplyr::left_join(module_label_df, by = "module") |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      cx = mean(x),
      cy = mean(y),
      module_label = dplyr::first(module_label),
      .groups = "drop"
    )

  # Combined label dataframe for mutual repulsion between node and module labels
  label_df_parts <- list()
  if (show_node_labels) {
    label_df_parts[["nodes"]] <- data.frame(
      lx       = lay$x,
      ly       = lay$y,
      lbl      = lay$node_id,
      txt_size = node_label_size,
      stringsAsFactors = FALSE
    )
  }
  if (show_labels) {
    label_df_parts[["modules"]] <- data.frame(
      lx       = centroids$cx,
      ly       = centroids$cy,
      lbl      = stringr::str_wrap(centroids$module_label, width = 25),
      txt_size = label_size,
      stringsAsFactors = FALSE
    )
  }
  all_label_df <- if (length(label_df_parts) > 0) do.call(rbind, label_df_parts) else NULL

  # Per-module node fill palette
  n_mod <- length(selected_modules)
  node_fill_palette <- colorRampPalette(c(
    "#0ca9ce", "#ff6f81", "#59A14F", "#B07AA1", "#F28E2B",
    "#76B7B2", "#E15759", "#4E79A7", "#F1CE63", "#d386bf",
    "#00b1a5", "#ffa68f", "#acd295", "#448c99", "#db888e"
  ))(n_mod)
  names(node_fill_palette) <- sort(selected_modules)

  # Edge linetype lookup
  edge_lty <- c(
    "TF-target" = "solid",
    "PPI" = "solid",
    "Reaction" = "solid",
    "molecule_pathway" = "solid",
    "pathway_similarity" = "solid",
    "diffusion_similarity" = "dashed"
  )

  # Assemble plot
  p <- ggraph::ggraph(lay) +

    # Edges
    ggraph::geom_edge_link(
      ggplot2::aes(
        colour = edge_type,
        linetype = edge_type,
        edge_width = display_weight
      ),
      alpha = 0.5
    ) +
    ggraph::scale_edge_width(range = c(0.2, 1), guide = "none") +
    ggraph::scale_edge_colour_manual(name = "Edge type", values = edge_colors) +
    ggraph::scale_edge_linetype_manual(name = "Edge type", values = edge_lty) +

    # Nodes
    ggraph::geom_node_point(
      ggplot2::aes(fill = module, shape = node_type),
      size = node_size
    ) +
    ggplot2::scale_fill_manual(name = "Module", values = node_fill_palette) +
    ggplot2::scale_shape_manual(name = "Node type", values = node_shapes) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21))
    ) +

    # Node and module labels — single repel layer so both sets push each other
    {
      if (!is.null(all_label_df))
        ggrepel::geom_text_repel(
          data        = all_label_df,
          mapping     = ggplot2::aes(x = lx, y = ly, label = lbl, size = txt_size),
          colour      = "grey10",
          max.overlaps = 30,
          seed        = seed,
          inherit.aes = FALSE,
          show.legend = FALSE
        )
    } +
    ggplot2::scale_size_identity(guide = "none") +

    ggraph::theme_graph() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::coord_cartesian(clip = "off")

  p
}
