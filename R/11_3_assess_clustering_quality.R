# setwd(r4projects::get_project_wd())
#
# load("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/biotext_functional_modules.rda")
# assess_result <- assess_clustering_quality(object = biotext_functional_modules)
# save(assess_result, file = "demo_data/updated_object_results_for_genes_ora/biotext_sim_result/assess_result.rda")
# assess_result$size_plot
# assess_result$evaluation_plot
# assess_result$quality_metrics

# head(assess_result$quality_metrics)
# assess_result$evaluation_plot

#' Assess Clustering Quality and Visualize Results
#'
#' @description
#' This function evaluates the quality of functional module clustering. It computes
#' several standard clustering validation metrics, including the average silhouette
#' score, Calinski-Harabasz index, and Davies-Bouldin index.
#'
#' @details
#' The function operates on a `functional_module` object. It first identifies
#' and filters out singleton modules (modules with only one node).
#'
#' Using the `sim` edge attribute from the object's graph data, it calculates a
#' distance matrix (`distance = 1 - similarity`). This matrix and the cluster
#' assignments are then used to compute the quality metrics.
#'
#' The final output is a list containing the calculated metrics in a data frame
#' and multiple plots for visual assessment.
#'
#' @param object A `functional_module` class object, which should be the result
#'   of running the function `get_functional_modules()`.
#' @param save_plot A logical value indicating whether to save the generated plot
#'   to a file. Defaults to `FALSE`. The plot saved is the `evaluation_plot`.
#' @param plot_path A character string specifying the file path where the plot
#'   will be saved. Defaults to `"clustering_quality_plot.png"`.
#' @param width A numeric value for the width of the saved plot in inches.
#'   Defaults to `12`.
#' @param height A numeric value for the height of the saved plot in inches.
#'   Defaults to `8`.
#' @param ignore_singleton A logical value indicating whether to exclude singleton
#'   modules (modules with only one node) from the quality assessment. Defaults to
#'   `FALSE`.
#'
#' @return
#' A list containing three elements:
#' \item{quality_metrics}{
#'   A data frame providing detailed metrics for each node.
#'   The columns include:
#'   \itemize{
#'     \item \code{module}: The identifier for the functional module.
#'     \item \code{size}: The total number of nodes in the module.
#'     \item \code{node}: The identifier for the individual node (e.g., gene or term).
#'     \item \code{cluster}: The numeric cluster ID assigned to the node.
#'     \item \code{neighbor}: The neighboring cluster used for silhouette width calculation.
#'     \item \code{sil_width}: The silhouette width of the node, indicating how well it fits within its cluster.
#'   }
#' }
#' \item{size_plot}{
#'   A `ggplot` object displaying a bar chart of the module size distribution,
#'   showing the number of nodes in each functional module.
#' }
#' \item{evaluation_plot}{
#'   A `ggplot` object that visualizes the silhouette scores for all nodes in
#'   the clusters. The plot title summarizes the overall clustering quality with:
#'   \itemize{
#'     \item Average silhouette score.
#'     \item The proportion of modules that are non-singletons.
#'     \item Calinski-Harabasz (CH) Index.
#'     \item Davies-Bouldin (DB) Index.
#'   }
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming 'enriched_functional_module' is a pre-existing object
#' # from get_functional_modules()
#'
#' # Run the quality assessment
#' quality_results <- assess_clustering_quality(
#'   object = enriched_functional_module
#' )
#'
#' # Print the evaluation plot to the console
#' print(quality_results$evaluation_plot)
#'
#' # Print the module size distribution plot
#' print(quality_results$size_plot)
#'
#' # Inspect the first few rows of the quality metrics table
#' head(quality_results$quality_metrics)
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom cluster silhouette
#' @importFrom igraph as_adjacency_matrix vertex_attr
#' @importFrom factoextra fviz_silhouette
#' @importFrom clusterSim index.G1 index.DB
#'
#' @export

assess_clustering_quality <- function(object,
                                      save_plot = FALSE,
                                      plot_path = "clustering_quality_plot.png",
                                      width = 12,
                                      height = 8,
                                      ignore_singleton = FALSE
                                      ) {

  graph_data <- object@merged_module$graph_data
  result_with_module <- object@merged_module$result_with_module

  module_sizes <- result_with_module |>
    dplyr::count(module, name = "size") |>
    dplyr::arrange(dplyr::desc(size))

  total_modules <- nrow(module_sizes)
  non_singleton_modules <- module_sizes |> dplyr::filter(size > 1)
  non_singleton_count <- nrow(non_singleton_modules)
  prop_non_singleton <- non_singleton_count / total_modules

  if (!ignore_singleton) {
    non_singleton_modules <- module_sizes |>
      dplyr::filter(size > 0)
  }

  result_filtered <- result_with_module |>
    dplyr::filter(module %in% non_singleton_modules$module)

  node_data <- igraph::vertex_attr(graph_data)
  node_data_df <- data.frame(
    node = node_data$node,
    module = node_data$module,
    stringsAsFactors = FALSE
  )

  filtered_node_data <- node_data_df |>
    dplyr::filter(module %in% non_singleton_modules$module)

  numeric_clusters <- as.numeric(sub("Functional_module_", "", filtered_node_data$module))
  names(numeric_clusters) <- filtered_node_data$node

  # Calculate silhouette scores for non-singleton clusters only
  adj_matrix <- igraph::as_adjacency_matrix(graph_data, attr = "sim", sparse = FALSE)
  sim_matrix <- as.data.frame(adj_matrix)
  diag(sim_matrix) <- 1
  colnames(sim_matrix) <- node_data$node
  rownames(sim_matrix) <- node_data$node

  filtered_sim_matrix <- sim_matrix[filtered_node_data$node, filtered_node_data$node]

  dist_matrix <- 1 - filtered_sim_matrix
  dist_obj <- as.dist(dist_matrix)

  # Calculate silhouette scores
  sil_scores <- cluster::silhouette(numeric_clusters, dist_obj)

  # Calinski-Harabasz Index
  ch_index <- clusterSim::index.G1(x = filtered_sim_matrix,
                                   cl = numeric_clusters)

  # Davies-Bouldin Index
  db_index <- clusterSim::index.DB(x = filtered_sim_matrix,
                                   cl = numeric_clusters)$DB

  if (any(is.na(sil_scores))) {
    new_avg_score <- NA
    sil_df <- data.frame(
      node = filtered_node_data$node,
      cluster = NA,
      neighbor = NA,
      sil_width = NA,
      stringsAsFactors = FALSE
    )
  } else {
    new_avg_score <- mean(sil_scores[, "sil_width"])

    sil_df <- data.frame(
      node = filtered_node_data$node,
      cluster = sil_scores[, "cluster"],
      neighbor = sil_scores[, "neighbor"],
      sil_width = sil_scores[, "sil_width"],
      stringsAsFactors = FALSE
    )
  }

  sil_score_df <-
    result_filtered |>
    dplyr::select(node, module) |>
    dplyr::right_join(sil_df, by = "node")

  quality_metrics <- non_singleton_modules |>
    dplyr::left_join(sil_score_df, by = "module") |>
    dplyr::arrange(dplyr::desc(size))

  # Create module size distribution plot
  size_plot <-
    ggplot(module_sizes, aes(x = reorder(module, size), y = size)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = size),
              hjust = -0.1,
              size = 3,
              color = "black") +
    coord_flip() +
    labs(
      title = "Module Size Distribution",
      x = "Functional Module",
      y = "Number of Nodes"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    )

  # Generate colors for clusters
  n_clusters <- length(unique(numeric_clusters))
  colors <- colorRampPalette(c('#0ca9ce', '#78cfe5', '#c6ecf1', '#ff6f81', '#ff9c8f', '#ffc2c0','#d386bf',
                               '#cdb1d2', '#fae6f0', '#eb6fa6', '#ff88b5', '#00b1a5',"#ffa68f","#ffca75","#97bc83","#acd295",
                               "#00ada1","#009f93","#ace2da","#448c99","#00b3bc","#b8d8c9","#db888e","#e397a4","#ead0c7",
                               "#8f9898","#bfcfcb"))(n_clusters)

  names(colors) <- sort(unique(numeric_clusters))

  # Reorder silhouette scores for visualization
  # sil_scores_reordered <- sil_scores
  # sil_scores_reordered[, 1] <- forcats::fct_reorder(
  #   factor(sil_scores[, 1]),
  #   sil_scores[, 3],
  #   .desc = TRUE
  # )

  if (!is.na(sil_scores)) {
    # Create the evaluation plot
    evaluation_plot <-
      factoextra::fviz_silhouette(sil_scores, label = FALSE, print.summary = FALSE) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      scale_y_continuous(limits = c(-1, 1)) +
      theme_bw() +
      labs(
        title = paste(
          "Avg Silhouette =", round(new_avg_score, 3),
          "| Non-singleton Prop =", round(prop_non_singleton, 3),
          "| CH Index =", round(ch_index, 3),
          "| DB Index =", round(db_index, 3)
        )
      ) +
      theme(
        panel.grid = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(size = 10, hjust = 0.5)
      )
  } else {
    evaluation_plot <- ggplot() +
      theme_void() +
      labs(
        title = paste(
          "Silhouette plot unavailable |",
          "Non-singleton Prop =", round(prop_non_singleton, 3),
          "| CH Index =", round(ch_index, 3),
          "| DB Index =", round(db_index, 3)
        )
      ) +
      theme(plot.title = element_text(size = 10, hjust = 0.5))
  }

  if (save_plot) {
    ggsave(plot_path, evaluation_plot, width = width, height = height)
    message(paste("Plot saved to:", plot_path))
  }

  message("=== Clustering Quality Assessment Summary ===")
  message(paste("Total modules:", total_modules))
  message(paste("Non-singleton modules:", non_singleton_count))
  message(paste("Proportion of non-singleton modules:", round(prop_non_singleton, 3)))
  message(paste("Average silhouette score:", round(new_avg_score, 3)))
  message(paste("Calinski-Harabasz Index:", round(ch_index, 3)))
  message(paste("Davies-Bouldin Index:", round(db_index, 3)))

  # Return results
  return(list(
    quality_metrics = quality_metrics,
    size_plot = size_plot,
    evaluation_plot = evaluation_plot
  ))
}
