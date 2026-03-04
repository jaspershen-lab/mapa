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
# node_meta <- rbind(mnet_obj@mol_nodes, mnet_obj@path_nodes)
# save(graph_data, file = "demo_data/demo_multi-omics/graph_data.rda")
#
# load("demo_data/demo_multi-omics/mnet_obj.rda")
# load("demo_data/demo_multi-omics/sim_matrix_multi_omics.rda")
# multi_omics_modules <- get_functional_modules(
#   object = mnet_obj,
#   sim_matrix = sim_matrix,
#   sim_cutoff = 0.45,
#   cluster_method = "louvain"
# )
# save(multi_omics_modules, file = "demo_data/demo_multi-omics/multi_omics_modules.rda")

merge_multi_omics_nodes <- function(
    mnet_obj,
    sim_matrix,
    sim_cutoff = 0.55,
    cluster_method = "louvain",
    verbose = TRUE
) {
  # similarity matrix → clustered tidygraph
  if (verbose) message("\n Clustering multi-omics molecules and enriched pathways ...")
  node_meta <- rbind(mnet_obj@mol_nodes, mnet_obj@path_nodes)

  graph_data <- cluster_nodes(
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
    dplyr::mutate(
      module = stringr::str_replace(module, "^Module_", "Functional_module_")
    )

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
      genes = if_else(include_genes, gene, NA_character_),
      include_metabolites = !is.na(metabolite) & metabolite != "",
      metabolites = if_else(include_metabolites, metabolite, NA_character_),
      include_pathways = !is.na(pathway) & pathway != "",
      pathways = if_else(include_pathways, pathway, NA_character_)
    ) |>
    dplyr::mutate(multi_omics_num = include_genes + include_metabolites + include_pathways) |>
    dplyr::select(module, module_content_number, multi_omics_num, everything())

  list(
    graph_data = graph_data,
    functional_module_result = functional_module_result,
    result_with_module = result_with_module
  )
}

cluster_nodes <- function(
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

  # Build edge table (upper triangle only, no self-loops)
  message("Building edge table...")
  edge_data <- as.data.frame(as.table(sim_matrix), responseName = "sim")
  colnames(edge_data) <- c("from", "to", "sim")
  edge_data$from <- as.character(edge_data$from)
  edge_data$to <- as.character(edge_data$to)
  edge_data <- edge_data[edge_data$from < edge_data$to, ]   # upper triangle

  if (nrow(edge_data) == 0) stop("Edge table empty – need >= 2 nodes.")

  # Align node_meta to matrix row order
  node_data <- merge(
    data.frame(node_id = all_nodes, stringsAsFactors = FALSE),
    node_meta, by = "node_id", all.x = TRUE
  )
  node_data <- node_data[match(all_nodes, node_data$node_id), ]
  node_data$node <- node_data$node_id   # alias used by tidygraph key
  rownames(node_data) <- NULL

  message(sprintf("Nodes: %d  |  Edges before filter: %d",
                  nrow(node_data), nrow(edge_data)))

  # Handle hierarchical prefix
  hclust_method <- NULL
  if (grepl("^h_", cluster_method)) {
    hclust_method  <- sub("^h_", "", cluster_method)
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
    data.frame(node   = node_data$node,
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

  # Build tidygraph
  message("Building tidygraph object...")

  graph_data <-
    tidygraph::tbl_graph(
      nodes    = node_data,
      edges    = edge_data[edge_data$sim > sim_cutoff, ],
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

  message("Done. Modules found: ", length(unique(cluster_result$module)))

  graph_data
}

