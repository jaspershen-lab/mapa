# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/multiplex_het.rda")
# diffusion_profile <- compute_diffusion_profiles(multiplex_het)
# save(diffusion_profile, file = "demo_data/demo_multi-omics/diffusion_profile.rda")
# load("demo_data/demo_multi-omics/diffusion_profile.rda")
# diffusion_profile |> dim()
# colnames(diffusion_profile) |> unique() |> length()
# hist(diffusion_profile)
# hist(diffusion_profile[diffusion_profile > 0.5])
# range(diffusion_profile)
# cosine_sim <- .calculate_sim(diffusion_profile)
# hist(cosine_sim@x[cosine_sim@x > 0.6 & cosine_sim@x < 0.9])
# node_meta <- rbind(object@mol_nodes, object@path_nodes)
# save(graph_data, file = "demo_data/demo_multi-omics/graph_data.rda")
#
# load("demo_data/demo_multi-omics/mnet_obj.rda")
# load("demo_data/demo_multi-omics/sim_matrix_multi_omics.rda")
# set.seed(123)
# multi_omics_modules <- merge_multi_omics_nodes(
#   object = mnet_obj,
#   sim_matrix = sim_matrix,
#   sim_cutoff = 0.45,
#   cluster_method = "louvain"
# )
# save(multi_omics_modules, file = "demo_data/demo_multi-omics/multi_omics_modules.rda")

#' Merge Multi-Omics Nodes into Functional Modules
#'
#' Clusters molecules (genes, metabolites) and enriched pathways from a
#' \code{multi_omics_functional_module} object into functional modules using a
#' pre-computed cosine-similarity matrix.
#'
#' @param object A \code{multi_omics_functional_module} S4 object produced by
#'   [build_MNetwork()].
#' @param sim_matrix A square cosine-similarity matrix (rows/columns = node IDs),
#'   as returned by [get_multi_omics_sim()].
#' @param sim_cutoff Numeric. Minimum similarity to retain an edge. Default \code{0.55}.
#' @param cluster_method Character. Clustering algorithm. One of \code{"louvain"},
#'   \code{"walktrap"}, \code{"h_ward.D"}, \code{"binary_cut"}, etc.
#'   Default \code{"louvain"}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{graph_data}}{A \code{tbl_graph} with node and edge annotations.}
#'   \item{\code{functional_module_result}}{A data frame summarising each module.}
#'   \item{\code{result_with_module}}{A data frame of all nodes with module assignments.}
#' }
#'
#' @export
merge_multi_omics_nodes <- function(
    object,
    sim_matrix,
    sim_cutoff = 0.55,
    cluster_method = "louvain",
    verbose = TRUE
) {
  # similarity matrix → clustered tidygraph
  if (verbose) message("\n Clustering multi-omics molecules and enriched pathways ...")
  node_meta <- rbind(object@mol_nodes, object@path_nodes)

  graph_data <- cluster_nodes(
    object = object,
    sim_matrix = sim_matrix,
    sim_cutoff = sim_cutoff,
    node_meta = node_meta,
    cluster_method = cluster_method
  )

  if (verbose) message("\n Merge multi-omics nodes complete!")

  # result_with_module
  result_with_module <- graph_data |>
    tidygraph::activate(what = "nodes") |>
    tibble::as_tibble() |>
    dplyr::mutate(module = stringr::str_replace(module, "^Module_", "Functional_module_"))

  # functional_module_result
  x <- result_with_module |> dplyr::mutate(node_type = tolower(node_type))

  module_sizes <- x |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      module_content_number = max(module_size, na.rm = TRUE),
      .groups = "drop"
    )

  module_wide <- x |>
    dplyr::group_by(module, node_type) |>
    dplyr::summarise(
      ids = paste(unique(node_id), collapse = "/"),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(names_from = node_type, values_from = ids)

  needed_cols <- c("gene", "metabolite", "pathway")
  missing_cols <- setdiff(needed_cols, colnames(module_wide))
  if (length(missing_cols) > 0) module_wide[missing_cols] <- NA_character_

  functional_module_result <- module_sizes |>
    dplyr::left_join(module_wide, by = "module") |>
    dplyr::transmute(
      module,
      module_content_number,
      include_genes = !is.na(gene) & gene != "",
      genes = ifelse(include_genes, gene, NA_character_),
      include_metabolites = !is.na(metabolite) & metabolite != "",
      metabolites = ifelse(include_metabolites, metabolite, NA_character_),
      include_pathways = !is.na(pathway) & pathway != "",
      pathways = ifelse(include_pathways, pathway, NA_character_)
    ) |>
    dplyr::mutate(multi_omics_num = include_genes + include_metabolites + include_pathways) |>
    dplyr::select(module, module_content_number, multi_omics_num, everything())

  list(
    graph_data = graph_data,
    functional_module_result = functional_module_result,
    result_with_module = result_with_module
  )
}

#' Cluster Nodes by Similarity
#'
#' @param object A \code{multi_omics_functional_module} S4 object.
#' @param sim_matrix Square similarity matrix with node IDs as row/column names.
#' @param sim_cutoff Numeric. Minimum similarity to retain an edge. Default \code{0.55}.
#' @param node_meta Data frame of node metadata (must contain \code{node_id} and \code{node_type}).
#' @param cluster_method Character. Clustering algorithm. Default \code{"louvain"}.
#'
#' @return A \code{tbl_graph} with nodes annotated by module membership.
#'
#' @noRd
cluster_nodes <- function(
    object,
    sim_matrix,
    sim_cutoff = 0.55,
    node_meta,
    cluster_method  = "louvain"
) {

  available_methods <- c(
    "h_ward.D", "h_ward.D2", "h_single", "h_complete",
    "h_average", "h_mcquitty", "h_median", "h_centroid",
    "binary_cut",
    "louvain", "walktrap", "infomap", "edge_betweenness",
    "fast_greedy", "label_prop", "leading_eigen", "optimal"
  )

  if (!(cluster_method %in% available_methods)) {
    stop(paste0(
      "Invalid cluster_method: '", cluster_method, "'.\n",
      "Available: ", paste(available_methods, collapse = ", ")
    ))
  }

  if (inherits(sim_matrix, "Matrix")) sim_matrix <- as.matrix(sim_matrix)
  stopifnot(is.matrix(sim_matrix), nrow(sim_matrix) == ncol(sim_matrix))

  all_nodes <- rownames(sim_matrix)
  if (is.null(all_nodes)) stop("sim_matrix must have row/col names = node IDs.")

  message("Building edge table...")
  edge_data <- as.data.frame(as.table(sim_matrix), responseName = "sim")
  colnames(edge_data) <- c("from", "to", "sim")
  edge_data$from <- as.character(edge_data$from)
  edge_data$to <- as.character(edge_data$to)
  edge_data <- edge_data[edge_data$from < edge_data$to, ]

  if (nrow(edge_data) == 0) stop("Edge table empty – need >= 2 nodes.")

  # Align node_meta to matrix row order
  node_data <- merge(
    data.frame(node_id = all_nodes, stringsAsFactors = FALSE),
    node_meta, by = "node_id", all.x = TRUE
  )
  node_data <- node_data[match(all_nodes, node_data$node_id), ]
  node_data$node <- node_data$node_id
  rownames(node_data) <- NULL

  # Handle hierarchical prefix
  hclust_method <- NULL
  if (grepl("^h_", cluster_method)) {
    hclust_method <- sub("^h_", "", cluster_method)
    cluster_method <- "hierarchical"
  }

  .assign_alone <- function() {
    data.frame(node = node_data$node,
               module = paste0("Module_", seq_len(nrow(node_data))),
               stringsAsFactors = FALSE)
  }

  .make_graph <- function(fe) {
    igraph::graph_from_data_frame(
      fe, directed = FALSE,
      vertices = node_data[, c("node", setdiff(colnames(node_data), "node"))]
    )
  }

  .membership_df <- function(comm) {
    data.frame(node = node_data$node,
               module = paste0("Module_", as.character(igraph::membership(comm))),
               stringsAsFactors = FALSE)
  }

  # Clustering
  message(sprintf("Clustering: method = '%s', sim_cutoff = %.2f",
                  cluster_method, sim_cutoff))

  cluster_result <- switch(
    cluster_method,

    "hierarchical" = {
      dist_mat <- as.dist(1 - sim_matrix)
      hc <- hclust(dist_mat, method = hclust_method)
      groups <- cutree(hc, h = 1 - sim_cutoff)
      data.frame(node = names(groups),
                 module = paste0("Module_", groups),
                 stringsAsFactors = FALSE)
    },

    "binary_cut" = {
      if (!requireNamespace("simplifyEnrichment", quietly = TRUE))
        stop("Install 'simplifyEnrichment': BiocManager::install('simplifyEnrichment')")
      km <- simplifyEnrichment::binary_cut(sim_matrix, cutoff = sim_cutoff)
      data.frame(node = rownames(sim_matrix),
                 module = paste0("Module_", km),
                 stringsAsFactors = FALSE)
    },

    {
      fe <- edge_data[edge_data$sim >= sim_cutoff, ]
      if (nrow(fe) == 0) return(.assign_alone())
      g  <- .make_graph(fe)

      comm <- switch(
        cluster_method,
        "louvain" = igraph::cluster_louvain(g, weights = igraph::E(g)$sim),
        "walktrap" = igraph::cluster_walktrap(g, weights = igraph::E(g)$sim),
        "infomap" = igraph::cluster_infomap(g, e.weights = igraph::E(g)$sim),
        "edge_betweenness" = {
          igraph::E(g)$dist <- 1 - igraph::E(g)$sim   # distance weight
          igraph::cluster_edge_betweenness(g, weights = igraph::E(g)$dist)
        },
        "fast_greedy" = igraph::cluster_fast_greedy(g, weights = igraph::E(g)$sim),
        "label_prop" = igraph::cluster_label_prop(g, weights = igraph::E(g)$sim),
        "leading_eigen" = igraph::cluster_leading_eigen(g, weights = igraph::E(g)$sim),
        "optimal" = igraph::cluster_optimal(g, weights = igraph::E(g)$sim),
        stop("Unknown method: ", cluster_method)
      )
      .membership_df(comm)
    }
  )

  message("Annotating edges with knowledge layer information...")

  .ke <- function(df, etype) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df |>
      dplyr::select(from, to, weight) |>
      dplyr::mutate(
        from = as.character(from),
        to = as.character(to)
      ) |>
      dplyr::mutate(
        edge_type = etype,
        from_c = pmin(from, to),
        to_c = pmax(from, to)
      ) |>
      dplyr::select(from = from_c, to = to_c, weight, edge_type)
  }

  knowledge_edges <- dplyr::bind_rows(
    .ke(object@mol_layers_edge_weight$tf_target, "TF-target"),
    .ke(object@mol_layers_edge_weight$ppi, "PPI"),
    .ke(object@mol_layers_edge_weight$metabolite_reaction, "Reaction"),
    .ke(object@mol_layers_edge_weight$enzyme_metabolite, "Reaction"),
    .ke(object@pathway_mol_edge_weight |>
          dplyr::rename(from = pathway_id, to = mol_id) |>
          dplyr::mutate(weight = 1), "molecule_pathway"),
    .ke(object@path_layer_edge_weight |>
          dplyr::filter(weight >= sim_cutoff), "pathway_similarity")
  ) |>
    dplyr::distinct(from, to, edge_type, .keep_all = TRUE)

  # Annotate diffusion edges
  # Filtered diffusion edges (sim > sim_cutoff) are the backbone of the graph.
  # Each diffusion edge is left-joined to knowledge_edges on (from, to):
  #   - Matching edges expand into one row per knowledge edge_type, carrying
  #     both the knowledge `weight` and the diffusion similarity as `diff_weight`.
  #   - Non-matching edges receive edge_type = "diffusion_similarity" and
  #     weight = NA_real_, retaining only the `diff_weight`.
  filtered_edges <- edge_data[edge_data$sim > sim_cutoff, ] |>
    dplyr::rename(diff_weight = sim)

  annotated_edges <- filtered_edges |>
    dplyr::left_join(knowledge_edges, by = c("from", "to")) |>
    dplyr::mutate(
      edge_type = dplyr::if_else(is.na(edge_type), "diffusion_similarity", edge_type),
      weight = dplyr::if_else(edge_type == "diffusion_similarity", NA_real_, weight)
    ) |>
    dplyr::select(from, to, diff_weight, edge_type, weight)

  # Build tidygraph
  message("Building tidygraph object...")

  graph_data <-
    tidygraph::tbl_graph(
      nodes = node_data,
      edges = annotated_edges,
      directed = FALSE,
      node_key = "node"
    ) |>
    dplyr::mutate(degree = tidygraph::centrality_degree()) |>
    dplyr::left_join(cluster_result, by = "node")

  # result_with_module
  result_with_module <-
    igraph::vertex_attr(graph_data) |>
    do.call(what = cbind) |>
    as.data.frame() |>
    dplyr::mutate(module = as.character(unlist(module)))

  mod_count <- dplyr::count(result_with_module, module, name = "module_size")

  result_with_module <- dplyr::left_join(result_with_module, mod_count, by = "module") |>
    dplyr::arrange(module)

  graph_data <- graph_data |>
    tidygraph::activate("nodes") |>
    dplyr::left_join(mod_count, by = "module")

  n_modules <- result_with_module |>
    dplyr::filter(module_size > 3) |>
    dplyr::pull(module) |>
    unique() |>
    length()
  message("Done. Multi-omics modules with size > 3: ", n_modules)

  graph_data
}
