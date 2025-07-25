# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# devtools::load_all()
# object is list from biotext embedding
# setwd("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/")
# load("openai_semantic_sim_matrix.rda")
# results <- determine_optimal_clusters(object = openai_semantic_sim_matrix)
# results$evaluation_plot
# results$best_combination
# results$cluster_result
# save(results, file = "results.rda")
####
# object is functional module object from overlap
# load("demo_data/updated_object_results_for_genes_ora/gene_overlap_result/ora_enriched_modules.rda")
# results <- determine_optimal_clusters(object = enriched_modules,
#                                       cutoff_increment = 0.1,
#                                       cutoff_range = c(0.1, 0.8))
# results$evaluation_plot
# results$best_combination

#' Determine Optimal Clustering Parameters for Pathway Enrichment Results
#'
#' This function systematically evaluates various clustering strategies to identify the
#' optimal parameters for grouping pathway enrichment results. It iterates through
#' different clustering methods, including network-based community detection,
#' distance-based hierarchical clustering, and a recursive divisive algorithm.
#' Each parameter combination is assessed using modularity and silhouette scores
#' to quantify cluster quality.
#'
#' @details
#' The function handles two types of input for the `object` parameter:
#' 1.  A `list` that must contain a pre-computed `sim_matrix` (a similarity matrix)
#'     and an `enriched_pathway` object. This is typically used for results from
#'     biotext embedding.
#' 2.  An S4 `functional_module` object from the `mapa` package, which has been
#'     processed by the `mapa::merge_pathways()` function. In this case,
#'     the function will first calculate a Jaccard similarity matrix internally
#'     based on the overlapping genes/metabolites between pathways before
#'     proceeding with the evaluation.
#'
#' Available clustering methods:
#' - Network-based: "louvain", "walktrap", "infomap", "edge_betweenness", "fast_greedy",
#'   "label_prop", "leading_eigen", "optimal"
#' - Distance-based: "hierarchical" (supports various linkage methods)
#' - Divisive: "binary_cut"
#'
#' @param object An object containing pathway enrichment results. This can be a `list`
#'   with `sim_matrix` and `enriched_pathway` elements, or a `functional_module` S4 object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A list containing three elements:
#'   \item{cluster_result}{A data frame in long format with the complete
#'     evaluation results, including method, cutoff, metric type (modularity or
#'     silhouette), and the calculated metric values.}
#'   \item{evaluation_plot}{A comprehensive `ggplot` object with a heatmap-like
#'     visualization showing performance across methods and cutoffs. This plot
#'     makes optimal parameter selection visually intuitive, with color representing
#'     modularity and point size representing the silhouette score.}
#'   \item{best_combination}{A data frame identifying the best-performing
#'     method and cutoff combination for each metric (modularity and silhouette).
#'     Ties for the top score will result in multiple rows per metric.}
#'
#' @examples
#' \dontrun{
#' # Assuming 'enrichment_obj' is your S4 enrichment analysis result from mapa
#' # that has been processed by merge_pathways().
#' cluster_eval <- determine_optimal_clusters(enrichment_obj)
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom igraph graph_from_data_frame cluster_louvain cluster_walktrap cluster_infomap cluster_edge_betweenness cluster_fast_greedy cluster_label_prop cluster_leading_eigen cluster_optimal membership V E
#' @importFrom cluster silhouette
#'
#' @export

determine_optimal_clusters <- function(object, ...) {
  UseMethod("determine_optimal_clusters")
}


#' Determine Optimal Clustering for Functional Modules
#'
#' Identify the best clustering configuration for a **`functional_module`** object.
#' The function builds a Jaccard similarity matrix for all pathways and evaluates
#' multiple clustering algorithms across a range of similarity cutoffs, selecting
#' the combination that maximizes a composite of modularity and silhouette scores.
#'
#' @param object An **S4** *functional_module* object (typically the output of
#'   `mapa::get_functional_modules()`).
#' @param cutoff_range Numeric vector of length 2 giving the lower and upper
#'   bounds (inclusive) of the similarity threshold to test.
#'   Must be in `[0, 1]`; default is `c(0.2, 0.9)`.
#' @param cutoff_increment Numeric step size used to generate the sequence of
#'   cutoffs inside `cutoff_range`. Default is `0.05`.
#' @param methods Character vector of clustering approaches to evaluate.
#'   The following values are recognised (any mix allowed):
#'   \itemize{
#'     \item **Hierarchical** ‒ supply `"h_<agglom.method>"`, where
#'       `<agglom.method>` is one of
#'       `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
#'       `"mcquitty"`, `"median"`, `"centroid"`.
#'     \item `"binary_cut"`
#'     \item Graph-based: `"louvain"`, `"walktrap"`, `"infomap"`,
#'       `"edge_betweenness"`, `"fast_greedy"`, `"label_prop"`,
#'       `"leading_eigen"`, `"optimal"`
#'   }
#'   Default: `c("h_ward.D", "binary_cut", "louvain")`.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A **list** with three elements:
#' * **`cluster_result`** – a `data.frame` summarising every
#'   method/cut-off combination with its modularity and silhouette scores;
#' * **`evaluation_plot`** – a `ggplot2` object visualising the score landscape;
#' * **`best_combination`** – the row of `cluster_result` corresponding to the
#'   globally optimal settings.
#'
#'
#' @examples
#' \dontrun{
#' # Assume `enrichment_obj` is a processed functional_module object
#' # Explore graph-based algorithms over a narrower similarity range
#' res <- determine_optimal_clusters(
#'   enrichment_obj,
#'   cutoff_range = c(0.1, 0.7),
#'   methods = c("louvain", "walktrap", "fast_greedy", "optimal")
#' )
#'
#' # Retrieve the winning settings
#' res$best_combination
#' }
#'
#' @method determine_optimal_clusters functional_module
#' @export

determine_optimal_clusters.functional_module <-
  function(object,
           cutoff_range = c(0.2, 0.9),
           cutoff_increment = 0.05,
           methods = c("h_ward.D", "binary_cut", "louvain"),
           ...) {
    message("Starting optimal cluster determination for functional_module object...")

    # Validate inputs
    if (!is.numeric(cutoff_range) ||
        length(cutoff_range) != 2 ||
        cutoff_range[1] >= cutoff_range[2]) {
      stop(
        "`cutoff_range` must be a numeric vector of length 2 with the first element smaller than the second."
      )
    }

    available_methods <- c(
      "h_ward.D", "h_ward.D2", "h_single", "h_complete",
      "h_average", "h_mcquitty", "h_median", "h_centroid",
      "binary_cut", "louvain", "walktrap", "infomap",
      "edge_betweenness", "fast_greedy", "label_prop", "leading_eigen",
      "optimal"
    )

    if (!all(methods %in% available_methods)) {
      invalid_methods <- methods[!methods %in% available_methods]
      stop(paste(
        "Invalid methods:",
        paste(invalid_methods, collapse = ", "),
        "\nAvailable methods:",
        paste(available_methods, collapse = ", ")
      ))
    }

    message("Extracting similarity matrix and pathway data...")

    # Data preparation for S4 object
    variable_info <- object@variable_info

    ## Check if it has been processed by merge_pathways
    if (all(names(object@process_info) != "merge_pathways")) {
      stop("Use merge_pathways() function to process first")
    }

    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <-
        object@process_info$enrich_pathway@parameter$query_type
    } else{
      analysis_type <- "do_gsea"
      query_type <- object@process_info$do_gsea@parameter$query_type
    }

    ## Calculate the similarity (jaccard index) between all the pathways
    if (length(object@merged_pathway_go) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_go <-
          object@merged_pathway_go$module_result |>
          # dplyr::filter(ONTOLOGY != "CC") |>
          dplyr::arrange(p_adjust) |>
          dplyr::mutate(database = "GO") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            BgRatio,
            RichFactor,
            FoldEnrichment,
            zScore,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            database,
            geneID,
            pathway_id = node
          ) |>
          dplyr::mutate(Count = as.numeric(Count)) |>
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_go <-
          object@merged_pathway_go$module_result |>
          # dplyr::filter(ONTOLOGY != "CC") |>
          dplyr::arrange(dplyr::desc(abs(NES))) |>
          dplyr::mutate(database = "GO") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            setSize,
            enrichmentScore,
            NES,
            rank,
            leading_edge,
            core_enrichment,
            database,
            pathway_id = node
          ) |>
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_go <- NULL
    }

    if (length(object@merged_pathway_kegg) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_kegg <-
          object@merged_pathway_kegg$module_result |>
          dplyr::arrange(p_adjust) |>
          dplyr::mutate(database = "KEGG") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            BgRatio,
            RichFactor,
            FoldEnrichment,
            zScore,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            database,
            geneID,
            pathway_id = node
          ) |>
          dplyr::mutate(Count = as.numeric(Count)) |>
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_kegg <-
          object@merged_pathway_kegg$module_result |>
          dplyr::arrange(dplyr::desc(abs(NES))) |>
          dplyr::mutate(database = "KEGG") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            setSize,
            enrichmentScore,
            NES,
            rank,
            leading_edge,
            core_enrichment,
            database,
            pathway_id = node
          ) |>
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_kegg <- NULL
    }

    if (length(object@merged_pathway_reactome) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_reactome <-
          object@merged_pathway_reactome$module_result |>
          dplyr::arrange(p_adjust) |>
          dplyr::mutate(database = "Reactome") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            BgRatio,
            RichFactor,
            FoldEnrichment,
            zScore,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            database,
            geneID,
            pathway_id = node
          ) |>
          dplyr::mutate(Count = as.numeric(Count)) |>
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_reactome <-
          object@merged_pathway_reactome$module_result |>
          dplyr::arrange(dplyr::desc(abs(NES))) |>
          dplyr::mutate(database = "Reactome") |>
          dplyr::select(
            module_annotation,
            module,
            Description,
            pvalue,
            qvalue,
            p_adjust,
            Count,
            setSize,
            enrichmentScore,
            NES,
            rank,
            leading_edge,
            core_enrichment,
            database,
            pathway_id = node
          ) |>
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_reactome <- NULL
    }

    if (length(object@merged_pathway_hmdb) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_hmdb <-
          object@merged_pathway_hmdb$module_result |>
          dplyr::arrange(p_adjust) |>
          dplyr::mutate(database = "HMDB") |>
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_hmdb <-
          object@merged_pathway_hmdb$module_result |>
          dplyr::arrange(dplyr::desc(abs(NES))) |>
          dplyr::mutate(database = "HMDB") |>
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_hmdb <- NULL
    }

    if (length(object@merged_pathway_metkegg) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_metkegg <-
          object@merged_pathway_metkegg$module_result |>
          dplyr::arrange(p_adjust) |>
          dplyr::mutate(database = "KEGG") |>
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_metkegg <-
          object@merged_pathway_metkegg$module_result |>
          dplyr::arrange(dplyr::desc(abs(NES))) |>
          dplyr::mutate(database = "KEGG") |>
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_metkegg <- NULL
    }

    message("Calculating the similarity matrix...")
    jaccard_index <-
      get_jaccard_index_for_diff_databases(
        query_type = query_type,
        variable_info = object@variable_info,
        analysis_type = analysis_type,
        module_result_go = module_result_go,
        module_result_kegg = module_result_kegg,
        module_result_reactome = module_result_reactome,
        module_result_hmdb = module_result_hmdb,
        module_result_metkegg = module_result_metkegg
      )

    edge_data <-
      jaccard_index |>
      dplyr::rename(from = name1, to = name2, sim = value)

    rownames(edge_data) <- NULL

    # All unique module names from both columns
    all_names <- unique(c(jaccard_index$name1, jaccard_index$name2))

    # Initialize a square matrix with 0s on diagonal
    sim_matrix <- matrix(
      0,
      nrow = length(all_names),
      ncol = length(all_names),
      dimnames = list(all_names, all_names)
    )

    # Set diagonal to 1 (self-similarity)
    diag(sim_matrix) <- 1

    # Fill in the values
    for (i in seq_len(nrow(jaccard_index))) {
      row <- jaccard_index[i, ]
      sim_matrix[row$name1, row$name2] <- row$value
      sim_matrix[row$name2, row$name1] <- row$value
    }

    # Call the shared clustering evaluation function
    result <- perform_clustering_evaluation(
      sim_matrix = sim_matrix,
      edge_data = edge_data,
      methods = methods,
      cutoff_range = cutoff_range,
      cutoff_increment = cutoff_increment
    )

    message("Analysis complete!")
    return(result)
  }

#' Determine Optimal Clustering for Biotext Embedding Similarity
#'
#' Select the best clustering configuration when you already have a *semantic*
#' similarity matrix (derived from BioText embeddings).
#'
#' @param object A **list** with two named elements:
#'   * **`enriched_pathway`** – a *functional_module* object (output of
#'     `mapa::get_functional_modules()`);
#'   * **`sim_matrix`** – a square numeric matrix whose row/column names are the
#'     pathway identifiers to be clustered.
#' @param cutoff_range Numeric vector of length 2 giving the lower and upper
#'   bounds of similarity thresholds to test (inclusive), constrained to `[0, 1]`.
#'   Default: `c(0.2, 0.9)`.
#' @param cutoff_increment Numeric step size used to generate the sequence inside
#'   `cutoff_range`. Default: `0.05`.
#' @param methods Character vector of clustering algorithms to evaluate.
#'   Recognised identifiers (any mix allowed):
#'   \itemize{
#'     \item **Hierarchical** – prefix the agglomeration algorithm with `"h_"`:
#'       `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
#'       `"mcquitty"`, `"median"`, `"centroid"`.
#'     \item `"binary_cut"`
#'     \item Graph-based: `"louvain"`, `"walktrap"`, `"infomap"`,
#'       `"edge_betweenness"`, `"fast_greedy"`, `"label_prop"`,
#'       `"leading_eigen"`, `"optimal"`
#'   }
#'   Default: `c("h_ward.D", "binary_cut", "louvain")`.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return A **list** with three components:
#' * **`cluster_result`** – a data frame summarising every tested
#'   method/cut-off pair with their modularity and silhouette scores;
#' * **`evaluation_plot`** – a `ggplot2` object visualising the score landscape;
#' * **`best_combination`** – the row of `cluster_result` that achieved the
#'   overall optimum.
#'
#' @examples
#' \dontrun{
#' biotext_results <- list(
#'   enriched_pathway = my_pathway_results,  # functional_module object
#'   sim_matrix = embedding_similarity_matrix
#' )
#'
#' # Evaluate several graph-based algorithms on a finer cut-off grid
#' res <- determine_optimal_clusters(
#'   biotext_results,
#'   cutoff_increment = 0.01,
#'   cutoff_range = c(0.4, 0.95),
#'   methods = c("louvain", "walktrap", "fast_greedy", "optimal")
#' )
#'
#' res$best_combination
#' }
#'
#' @method determine_optimal_clusters list
#' @export

determine_optimal_clusters.list <-
  function(object,
           cutoff_range = c(0.2, 0.9),
           cutoff_increment = 0.05,
           methods = c("h_ward.D", "binary_cut", "louvain"),
           ...) {
    # Check if this is a bioembedding similarity object
    if (!all(c("enriched_pathway", "sim_matrix") %in% names(object))) {
      stop("List object must contain 'enriched_pathway' and 'sim_matrix' components.")
    }

    message("Starting optimal cluster determination for biotext embedding object...")

    # Validate inputs
    if (!is.numeric(cutoff_range) ||
        length(cutoff_range) != 2 ||
        cutoff_range[1] >= cutoff_range[2]) {
      stop(
        "`cutoff_range` must be a numeric vector of length 2 with the first element smaller than the second."
      )
    }

    available_methods <- c(
      "h_ward.D", "h_ward.D2", "h_single", "h_complete",
      "h_average", "h_mcquitty", "h_median", "h_centroid",
      "binary_cut", "louvain", "walktrap", "infomap",
      "edge_betweenness", "fast_greedy", "label_prop", "leading_eigen",
      "optimal"
    )

    if (!all(methods %in% available_methods)) {
      invalid_methods <- methods[!methods %in% available_methods]
      stop(paste(
        "Invalid methods:",
        paste(invalid_methods, collapse = ", "),
        "\nAvailable methods:",
        paste(available_methods, collapse = ", ")
      ))
    }

    message("Extracting similarity matrix and pathway data...")

    ## Data preparation for list object
    sim_matrix <- object$sim_matrix
    edge_data <-
      as.data.frame.table(object$sim_matrix, responseName = "sim") |>
      dplyr::filter(as.character(Var1) < as.character(Var2)) |> # Remove self-edges and duplicates
      dplyr::rename(from = Var1, to = Var2) |>
      dplyr::mutate(dplyr::across(c(from, to), as.character))

    # Call the shared clustering evaluation function
    result <- perform_clustering_evaluation(
      sim_matrix = sim_matrix,
      edge_data = edge_data,
      methods = methods,
      cutoff_range = cutoff_range,
      cutoff_increment = cutoff_increment
    )

    message("Analysis complete!")
    return(result)
  }


#' Perform Clustering Evaluation (Internal)
#'
#' @description
#' This internal helper function performs the core clustering evaluation across
#' various methods and cutoff thresholds.
#'
#' @param sim_matrix A numeric similarity matrix used for clustering.
#' @param edge_data A data frame of network edges with columns 'from', 'to', and 'sim'.
#' @param methods A character vector of clustering methods to evaluate.
#' @param cutoff_range A numeric vector of length 2 defining the start and end of the cutoff sequence.
#' @param cutoff_increment A numeric step size for the cutoff sequence.
#'
#' @return A list containing the full `cluster_result` data frame, the
#'   `evaluation_plot`, and the `best_combination` data frame.
#'
#' @noRd
perform_clustering_evaluation <- function(sim_matrix,
                                          edge_data,
                                          methods,
                                          cutoff_range,
                                          cutoff_increment) {
  message("Starting clustering evaluation across methods and cutoffs...")

  if (any(grepl("^h_", methods))) {
    hclust_methods <- gsub("^h_", "", methods[grepl("^h_", methods)])
    methods <- c("hierarchical", methods[!grepl("^h_", methods)])
  }

  # Validate hierarchical clustering methods
  valid_hclust_methods <-
    c("ward.D",
      "ward.D2",
      "single",
      "complete",
      "average",
      "mcquitty",
      "median",
      "centroid")

  if (!all(hclust_methods %in% valid_hclust_methods)) {
    invalid_hclust <-
      hclust_methods[!hclust_methods %in% valid_hclust_methods]
    stop(
      paste(
        "Invalid hierarchical clustering methods:",
        paste(invalid_hclust, collapse = ", "),
        "\nValid methods:",
        paste(valid_hclust_methods, collapse = ", ")
      )
    )
  }

  # Expand hierarchical methods if "hierarchical" is in methods
  expanded_methods <- c()
  for (method in methods) {
    if (method == "hierarchical") {
      # Add each hierarchical method as a separate method
      for (hclust_method in hclust_methods) {
        expanded_methods <-
          c(expanded_methods, paste("hierarchical", hclust_method, sep = "_"))
      }
    } else {
      expanded_methods <- c(expanded_methods, method)
    }
  }

  results <- list()
  cutoff <- seq(cutoff_range[1], cutoff_range[2], by = cutoff_increment)

  for (method in expanded_methods) {
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
        edge_data = edge_data,
        method = method,
        cutoff = cf
      )

      # Skip if clustering failed, resulted in a single cluster, or had issues
      if (is.null(clusters) ||
          length(unique(clusters)) <= 1 ||
          any(is.na(clusters)) || length(clusters) == 0) {
        method_results$modularity[i] <- NA
        method_results$silhouette[i] <- NA
        next
      }

      # Calculate metrics
      # Modularity: order of membership must match vertex order in the graph
      method_results$modularity[i] <-
        calculate_modularity(
          sim_matrix = sim_matrix,
          edge_data = edge_data,
          sim.cutoff = cf,
          clusters = clusters
        )
      # Silhouette: order of elements must match objects in the distance matrix
      method_results$silhouette[i] <-
        calculate_silhouette(sim_matrix = sim_matrix,
                             clusters = clusters)
    }

    results[[method]] <- method_results
  }

  # 3. Generate evaluation plots and find the best combination
  message("Clustering evaluation completed! Generating plots and finding optimal parameters...")

  # Combine results for plotting
  plot_data <- dplyr::bind_rows(results, .id = "method") |>
    tidyr::pivot_longer(
      cols = c(modularity, silhouette),
      names_to = "metric",
      values_to = "value"
    )

  # Create single comprehensive plot for easy parameter selection
  evaluation_plot <- create_evaluation_plot(plot_data = plot_data)

  # Find the best combination for each metric
  best_combination <- plot_data |>
    dplyr::group_by(metric) |>
    dplyr::slice_max(order_by = value,
                     n = 1,
                     with_ties = TRUE) |>
    dplyr::ungroup()

  # 4. Consolidate results
  all_res <- list(
    cluster_result = plot_data,
    evaluation_plot = evaluation_plot,
    best_combination = best_combination
  )

  return(all_res)
}

#' Generate Clusters Using a Specific Method
#'
#' @description
#' This internal function dispatches to various clustering algorithms based on the
#' specified method. It handles network-based, hierarchical, and binary cut methods.
#'
#' @param sim_matrix A numeric similarity matrix.
#' @param edge_data A data frame of network edges with 'from', 'to', and 'sim' columns.
#' @param method A character string specifying the clustering method (e.g., "louvain",
#'   "hierarchical_ward.D2").
#' @param cutoff A numeric similarity cutoff threshold.
#'
#' @return An integer vector of cluster assignments for each item in the
#'   similarity matrix. Returns `NULL` if the method fails.
#'
#' @importFrom igraph graph_from_data_frame cluster_louvain cluster_walktrap cluster_infomap cluster_edge_betweenness cluster_fast_greedy cluster_label_prop cluster_leading_eigen cluster_optimal membership E
#'
#' @noRd
generate_clustering <- function(sim_matrix,
                                edge_data,
                                method,
                                cutoff) {

  tryCatch({
    # Handle hierarchical methods with linkage specification
    if (grepl("hierarchical_", method)) {
      # Extract the linkage method from the method name
      hclust_method <- sub("hierarchical_", "", method)
      cluster_res <- merge_by_hierarchical(sim_matrix = sim_matrix,
                                           hclust.method = hclust_method,
                                           sim.cutoff = cutoff)$module
      return(as.integer(sub("Functional_module_", "", cluster_res)))
    }

    switch(method,
           "hierarchical" = {
             # Default hierarchical (shouldn't be called with new structure, but kept for safety)
             cluster_res <- merge_by_hierarchical(sim_matrix = sim_matrix,
                                                  hclust.method = "ward.D2",
                                                  sim.cutoff = cutoff)$module
             as.integer(sub("Functional_module_", "", cluster_res))
           },
           "binary_cut" = {
             # Use 1-cutoff for distance-based clustering
             cluster_res <- merge_by_binary_cut(sim_matrix = sim_matrix,
                                                sim.cutoff = cutoff)$module
             as.integer(sub("Functional_module_", "", cluster_res))
           },
           # Network-based methods
           "louvain" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_louvain(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "walktrap" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_walktrap(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "infomap" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_infomap(graph_obj, e.weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "edge_betweenness" = {
             # Use distance weights (1 - similarity) for edge betweenness
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ] |> dplyr::mutate(sim = 1 - sim)
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_edge_betweenness(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "fast_greedy" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_fast_greedy(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "label_prop" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_label_prop(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "leading_eigen" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_leading_eigen(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           },
           "optimal" = {
             filtered_edges <- edge_data[edge_data$sim >= cutoff, ]
             if (nrow(filtered_edges) == 0) return(rep(1, nrow(sim_matrix)))

             graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                        vertices = rownames(sim_matrix))
             comm <- igraph::cluster_optimal(graph_obj, weights = igraph::E(graph_obj)$sim)
             igraph::membership(comm)[rownames(sim_matrix)]
           }
    )
  }, error = function(e) {
    warning(paste("Clustering method", method, "failed at cutoff", cutoff, ":", e$message))
    return(NULL)
  })
}

#' Create Evaluation Plot
#'
#' @description
#' Generates a ggplot object to visualize the clustering evaluation results.
#' The plot uses a heatmap-style layout where the x-axis is the cutoff, the y-axis
#' is the clustering method, the fill color represents modularity, and the point
#' size represents the silhouette score.
#'
#' @param plot_data A data frame in long format containing the columns 'method',
#'   'cutoff', 'metric', and 'value'.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_wider
#'
#' @noRd
create_evaluation_plot <- function(plot_data) {
  plot_data <- plot_data |>
    tidyr::pivot_wider(names_from = metric, values_from = value) |>
    dplyr::arrange(method, cutoff)

  n_methods <- length(unique(plot_data$method))
  n_cutoffs <- length(unique(plot_data$cutoff))
  ratio <- n_methods / n_cutoffs

  heatmap_plot <- ggplot2::ggplot(data = plot_data,
                                  aes(x = cutoff, y = method)) +
    # Add a tile layer for the background grid
    ggplot2::geom_tile(color = "black",
                       fill = "white",
                       linewidth = 0.2) +

    # Add points (circles) for valid results
    ggplot2::geom_point(
      data = . %>% dplyr::filter(!is.na(modularity) & !is.na(silhouette)),
      aes(fill = modularity, size = silhouette),
      shape = 21,
      color = "black"
    ) +

    # Add "NA" text for cells with missing data
    ggplot2::geom_text(
      data = . %>% dplyr::filter(is.na(modularity) | is.na(silhouette)),
      aes(label = "NA"),
      size = 3,
      color = "gray50"
    ) +

    # --- Customize Scales and Colors ---
    # Color gradient for fill (modularity)
    ggplot2::scale_fill_gradient2(
      low = "#4877b5",
      mid = "#fbf6bb",
      high = "#d73226",
      midpoint = 0.5,
      name = "Modularity",
      na.value = "transparent"
    ) +

    # Scale for circle size (silhouette score)
    ggplot2::scale_size_continuous(range = c(4, 10),
                                   name = "Silhouette\nScore") +
    ggplot2::scale_x_continuous(breaks = unique(plot_data$cutoff)) +

    # --- Theming and Labels ---
    labs(x = "Cutoff",
         y = "Clustering Method") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      aspect.ratio = ratio
    )

  return(heatmap_plot)
}
