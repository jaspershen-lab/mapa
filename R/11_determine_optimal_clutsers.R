# setwd(r4projects::get_project_wd())
# library(mapa)
# source("R/6_utils.R")
# load("demo_data/pregnancy_data/results/results_biotext/openai_sim_matrix_met.rda")
# load("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/ora_openai_semantic_sim_matrix.rda")
# {
#   object <- openai_sim_matrix_met
#   hclust_method = "complete"
# }
# results <- determine_optimal_clusters(object = openai_sim_matrix_met)
# results$evaluation_plot

#' Determine Optimal Clusters for Pathway Enrichment Analysis Results
#'
#' This function evaluates different clustering methods and similarity cutoffs to
#' determine the optimal clustering parameters for pathway enrichment results.
#' It tests hierarchical clustering, binary cut, and Girvan-Newman methods across
#' a range of similarity cutoffs, evaluating each combination using modularity
#' and silhouette metrics. The function also provides informative messages to track
#' its progress.
#'
#' @param object An object containing enrichment results with a similarity matrix
#'   and enriched pathway data. The object should have components:
#'   \itemize{
#'     \item \code{sim_matrix}: A similarity matrix for pathways.
#'     \item \code{enriched_pathway}: Enrichment results containing process
#'       information and results from GO, KEGG, Reactome, or HMDB databases.
#'   }
#' @param cutoff_increment Numeric value specifying the step size for the
#'   similarity cutoff sequence, which ranges from 0.2 to 0.9. Default is 0.05.
#' @param hclust_method Character string specifying the agglomeration method for
#'   hierarchical clustering. Default is "complete". Options include "ward.D",
#'   "ward.D2", "single", "average", "mcquitty", "median", or "centroid".
#'
#' @return A list containing:
#'   \describe{
#'     \item{cluster_result}{A data frame with clustering evaluation results,
#'       including method, cutoff, metric type, and metric values.}
#'     \item{evaluation_plot}{A ggplot object showing validation metrics vs.
#'       similarity cutoffs for the different clustering methods.}
#'     \item{best_combination}{A data frame showing the best method-cutoff
#'       combination(s) for each evaluation metric. Note: This may include
#'       multiple rows per metric if there is a tie for the highest score.}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming 'enrichment_obj' is your enrichment analysis result
#' cluster_eval <- determine_optimal_clusters(enrichment_obj)
#'
#' # View the evaluation plot
#' print(cluster_eval$evaluation_plot)
#'
#' # Check the best combination(s)
#' print(cluster_eval$best_combination)
#'
#' # Run with a finer cutoff increment for a more detailed search
#' cluster_eval_fine <- determine_optimal_clusters(enrichment_obj,
#'                                                 cutoff_increment = 0.01)
#'
#' # Use a different hierarchical clustering method
#' cluster_eval_ward <- determine_optimal_clusters(enrichment_obj,
#'                                                 hclust_method = "ward.D2")
#' }
#'
#' @importFrom dplyr filter rename mutate across select bind_rows group_by slice_max ungroup
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_shape_manual scale_color_brewer scale_x_continuous labs theme_minimal theme element_text guides guide_legend
#' @importFrom igraph graph_from_adjacency_matrix
#'
#' @export

determine_optimal_clusters <- function(object,
                                       cutoff_increment = 0.05,
                                       hclust_method = "complete") {

  message("Starting optimal cluster determination...")
  message("Extracting similarity matrix and pathway data...")
  # 1. Data preparation ====
  ## 1.1 Collect sim_matrix
  sim_matrix <- object$sim_matrix

  ## 1.2 Collect node and edge data for Girvan-Newman clustering
  parameters <- object$enriched_pathway@process_info$merge_pathways@parameter
  query_type <- parameters$query_type
  ### Collect edge data
  edge_data <-
    as.data.frame.table(object$sim_matrix, responseName = "sim") |>
    dplyr::filter(Var1 != Var2) |>                 # Remove selfs-edges
    dplyr::rename(from = Var1, to = Var2) |>
    dplyr::mutate(dplyr::across(c(from, to), as.character)) |>
    dplyr::filter(from < to)

  ### Collect node data
  result <- data.frame()

  if (query_type == "gene") {
    if (!is.null(object$enriched_pathway@enrichment_go_result)) {
      result <-
        object$enriched_pathway@enrichment_go_result@result |>
        dplyr::select(-ONTOLOGY) |>
        dplyr::filter(p_adjust < parameters$p.adjust.cutoff.go) |>
        dplyr::filter(Count > parameters$count.cutoff.go) |>
        dplyr::mutate(database = "GO") |>
        rbind(result)
    }

    if (!is.null(object$enriched_pathway@enrichment_kegg_result)) {
      result <-
        object$enriched_pathway@enrichment_kegg_result@result |>
        (\(x) if ("enrich_pathway" %in% names(object$enriched_pathway@process_info))
          dplyr::select(x, -c(category, subcategory)) else x)() |>
        dplyr::filter(.data$p_adjust < parameters$p.adjust.cutoff.kegg,
                      .data$Count > parameters$count.cutoff.kegg) |>
        dplyr::mutate(database = "KEGG") |>
        rbind(result)
    }

    if (!is.null(object$enriched_pathway@enrichment_reactome_result)) {
      result <-
        object$enriched_pathway@enrichment_reactome_result@result |>
        dplyr::filter(p_adjust < parameters$p.adjust.cutoff.reactome) |>
        dplyr::filter(Count > parameters$count.cutoff.reactome) |>
        dplyr::mutate(database = "Reactome") |>
        rbind(result)
    }

    node_data <- result |> dplyr::rename(node = ID)

  } else if (query_type == "metabolite") {
    if (!is.null(object$enriched_pathway@enrichment_hmdb_result)) {
      result <-
        object$enriched_pathway@enrichment_hmdb_result@result |>
        dplyr::filter(p_adjust < parameters$p.adjust.cutoff.hmdb) |>
        dplyr::filter(mapped_number > parameters$count.cutoff.hmdb) |>
        dplyr::mutate(database = "HMDB") |>
        rbind(result)
    }

    if (!is.null(object$enriched_pathway@enrichment_metkegg_result)) {
      result <-
        object$enriched_pathway@enrichment_metkegg_result@result |>
        dplyr::filter(p_adjust < parameters$p.adjust.cutoff.metkegg) |>
        dplyr::filter(mapped_number > parameters$count.cutoff.metkegg) |>
        dplyr::mutate(database = "KEGG") |>
        rbind(result)
    }

    node_data <- result |> dplyr::rename(node = pathway_id)
  }

  # 2. Generate clusters and evaluate clusters ====
  message("Starting clustering evaluation across methods and cutoffs...")

  results <- list()

  methods <- c("hierarchical", "binary_cut", "girvan_newman")
  cutoff <- seq(0.2, 0.9, by = cutoff_increment)
  for (method in methods) {
    message(sprintf("Processing method: %s ...", method))
    method_results <- data.frame(
      cutoff = cutoff,
      modularity = numeric(length(cutoff)),
      silhouette = numeric(length(cutoff))
    )

    for (i in seq_along(cutoff)) {
      cf <- cutoff[i]

      # Generate clustering
      clusters <- generate_clustering(
        sim_matrix = sim_matrix,
        node_data = node_data,
        edge_data = edge_data,
        method = method,
        cutoff = cf,
        hclust_method = hclust_method
      )

      # Calculate metrics
      dist_matrix <- 1 - sim_matrix
      dist_obj <- as.dist(dist_matrix)

      # Create graph object for modularity
      graph_obj <- igraph::graph_from_adjacency_matrix(
        sim_matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )

      method_results$modularity[i] <- calculate_modularity(graph_obj, clusters)
      method_results$silhouette[i] <- calculate_silhouette(dist_obj, clusters)
    }

    results[[method]] <- method_results
  }

  # 3. Generate evaluation plot and find the best combination ====
  # Combine results for plotting
  message("Clustering evaluation completed! Generating plots and finding optimal parameters...")

  plot_data <- dplyr::bind_rows(results, .id = "method") |>
    tidyr::pivot_longer(cols = c(modularity, silhouette),
                        names_to = "metric", values_to = "value") |>
    # Ensure no NA/NaN/Inf values interfere with finding the maximum
    dplyr::filter(is.finite(value))

  # Create the plot
  eva_plot <- create_validation_plots(plot_data = plot_data)

  # Find the best combination for each metric
  best_combination <- plot_data |>
    dplyr::group_by(metric) |>
    dplyr::slice_max(order_by = value, n = 1, with_ties = TRUE) |>
    dplyr::ungroup()

  # 4. Consolidate results ====
  all_res <- list(
    cluster_result = plot_data,
    evaluation_plot = eva_plot,
    best_combination = best_combination
  )

  message("Analysis complete!")

  return(all_res)
}


#' Generate Clustering Results Using Different Methods
#'
#' Internal function that generates clustering results using one of three methods:
#' hierarchical clustering, binary cut, or Girvan-Newman algorithm.
#'
#' @param sim_matrix Similarity matrix for clustering
#' @param node_data Data frame containing node information
#' @param edge_data Data frame containing edge information with columns
#'   'from', 'to', and 'sim'
#' @param method Character string specifying clustering method. Options are:
#'   "hierarchical", "binary_cut", or "girvan_newman"
#' @param cutoff Numeric similarity cutoff threshold
#' @param hclust_method Character string specifying hierarchical clustering method.
#'   Only used when method = "hierarchical"
#'
#' @return Integer vector of cluster assignments
#'
#' @noRd

generate_clustering <- function(sim_matrix,
                                node_data,
                                edge_data,
                                method,
                                cutoff,
                                hclust_method = NULL) {

  switch(method,
         "hierarchical" = {
           cluster_res <- merge_by_hierarchical(sim_matrix = sim_matrix,
                                                hclust.method = hclust_method,
                                                sim.cutoff = cutoff)$module
           as.integer(sub("Functional_module_", "", cluster_res))
         },
         "binary_cut" = {
           cluster_res <- merge_by_binary_cut(sim_matrix = sim_matrix,
                                              sim.cutoff = cutoff)$module
           as.integer(sub("Functional_module_", "", cluster_res))
         },
         "girvan_newman" = {
           cluster_res <- merge_by_Girvan_Newman(edge_data = edge_data,
                                                 node_data = node_data,
                                                 sim.cutoff = cutoff)$module
           as.integer(sub("Functional_module_", "", cluster_res))
         }
  )
}

#' Create Validation Plots for Clustering Evaluation
#'
#' Internal function that creates ggplot2 visualization showing clustering
#' validation metrics (modularity and silhouette) across different similarity
#' cutoffs and clustering methods.
#'
#' @param plot_data Data frame containing clustering evaluation results with
#'   columns: method, cutoff, metric, value
#'
#' @return A ggplot object showing validation metrics vs similarity cutoffs
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_shape_manual
#'   scale_color_brewer scale_x_continuous labs theme_minimal theme element_text
#'   guides guide_legend
#' @importFrom dplyr filter
#'
#' @noRd

create_validation_plots <- function(plot_data) {
  # Modularity plot
  p <- plot_data |>
    # Remove NA values to avoid warnings in plotting
    filter(!is.na(value)) |>
    ggplot(aes(x = cutoff, y = value, color = method, shape = metric)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(aes(group = interaction(method, metric)), linewidth = 1, alpha = 0.7) +
    scale_shape_manual(
      name = "Metric",
      values = c("modularity" = 16, "silhouette" = 17),  # circle for modularity, triangle for silhouette
      labels = c("modularity" = "Modularity", "silhouette" = "Silhouette")
    ) +
    scale_color_brewer(
      name = "Method",
      type = "qual",
      palette = "Set1"
    ) +
    scale_x_continuous(breaks = seq(0, 1, 0.05)) +
    labs(
      title = "Clustering Validation Metrics vs Similarity Cutoff",
      x = "Similarity Cutoff",
      y = "Metric Value",
      caption = "Points represent different metrics; colors represent clustering methods"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    guides(
      color = guide_legend(title = "Clustering Method", override.aes = list(size = 3)),
      shape = guide_legend(title = "Metric", override.aes = list(size = 3))
    )

  return(p)
}
