# ####pathway enrichment
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
#
# load("demo_data/updated_object_results_for_genes_ora/gene_overlap_result/ora_enriched_modules.rda")
# enriched_functional_module <-
#   merge_modules(
#     object = enriched_modules,
#     sim.cutoff = 0.5,
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_functional_module, file = "result/enriched_functional_module")
#
#
# ##gsea results
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/covid_data/")
#
# load("result/enriched_modules")
#
# enriched_functional_module <-
#   merge_modules(
#     object = gsea_enriched_modules,
#     sim.cutoff = 0.5,
#     measure_method = c("jaccard"),
#     cluster_method = "girvan newman",
#     path = "result",
#     save_to_local = FALSE
# )
#
# save(enriched_functional_module, file = "result/enriched_functional_module")

# For metabolite
# object <- merged_pathways
# enriched_functional_module_met <- merge_modules(
#   object = merged_pathways,
#   sim.cutoff = 0,
#   measure_method = "jaccard",
#   cluster_method = "girvan newman"
# )

#' Identify Functional Modules Across Different Databases
#'
#' This function identify functional modules from various databases, such as GO, KEGG, and Reactome.
#' It calculates the similarity matrix between all the pathways and clusters them into functional modules.
#'
#' @param object An S4 object, expected to be processed by `merge_pathways()`.
#' @param sim.cutoff A numerical value for similarity cutoff, default is 0.5. This is the similarity cutoff for pathways.
#'  Only edges with similarity above the cutoff will be stored in the graph data.
#' @param measure_method A character vector indicating the method for measuring similarity between modules, default is "jaccard".
#' @param cluster_method Character, clustering method options:  "louvain" (default), "hierarchical",
#'   "binary_cut", "walktrap", "infomap", "edge_betweenness", "fast_greedy", "label_prop",
#'   "leading_eigen", "optimal".
#' @param hclust.method Character, agglomeration method for hierarchical clustering including:
#'   "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#'   Only used when `cluster_method = "hierarchical"`.
#' @param path A character string for the directory where to save results, default is "result".
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#'
#' @examples
#' \dontrun{
#' merge_modules(variable_info = my_var_info, object = my_object)
#' }
#'
#' @importFrom dplyr filter arrange mutate select left_join count rename distinct
#' @importFrom tidygraph tbl_graph centrality_degree activate
#' @importFrom igraph cluster_edge_betweenness edge_attr membership vertex_attr
#' @importFrom simplifyEnrichment binary_cut
#' @importFrom purrr map
#' @importFrom stringr str_split str_replace
#' @importFrom plyr dlply
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export

merge_modules <-
  function(object,
           sim.cutoff = 0.5,
           measure_method = c("jaccard"),
           cluster_method = "louvain",
           hclust.method = NULL,
           path = "result",
           save_to_local = FALSE) {

    variable_info <- object@variable_info

    available_methods <- c(
      "hierarchical",
      "binary_cut",
      "louvain",
      "walktrap",
      "infomap",
      "edge_betweenness",
      "fast_greedy",
      "label_prop",
      "leading_eigen",
      "optimal"
    )

    if (!(cluster_method %in% available_methods)) {
      stop(paste(
        "Invalid methods:",
        paste(cluster_method, collapse = ", "),
        "\nAvailable methods:",
        paste(available_methods, collapse = ", ")
      ))
    }

    if (save_to_local) {
      path <- file.path(path, "functional_modules")
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
      dir.create(
        file.path(path, "intermediate_data"),
        showWarnings = FALSE,
        recursive = TRUE
      )
    }

    ###check if it is been processed by merge_pathways
    if (all(names(object@process_info) != "merge_pathways")) {
      stop("Use merge_pathways() function to process first")
    }

    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <- object@process_info$enrich_pathway@parameter$query_type
    } else{
      analysis_type <- "do_gsea"
      query_type <- object@process_info$do_gsea@parameter$query_type
    }

    ######calculate the similarity (jaccard index) between all the pathways
    if (length(object@merged_pathway_go) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_go <-
          object@merged_pathway_go$module_result %>%
          # dplyr::filter(ONTOLOGY != "CC") %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::mutate(database = "GO") %>%
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
          ) %>%
          dplyr::mutate(Count = as.numeric(Count)) %>%
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_go <-
          object@merged_pathway_go$module_result %>%
          # dplyr::filter(ONTOLOGY != "CC") %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::mutate(database = "GO") %>%
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
          ) %>%
          dplyr::filter(!is.na(module_annotation))
      }

    } else{
      module_result_go <- NULL
    }

    if (length(object@merged_pathway_kegg) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_kegg <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::mutate(database = "KEGG") %>%
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
          ) %>%
          dplyr::mutate(Count = as.numeric(Count)) %>%
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_kegg <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::mutate(database = "KEGG") %>%
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
          ) %>%
          dplyr::filter(!is.na(module_annotation))
      }

    } else{
      module_result_kegg <- NULL
    }

    if (length(object@merged_pathway_reactome) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_reactome <-
          object@merged_pathway_reactome$module_result %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::mutate(database = "Reactome") %>%
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
          ) %>%
          dplyr::mutate(Count = as.numeric(Count)) %>%
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_reactome <-
          object@merged_pathway_reactome$module_result %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::mutate(database = "Reactome") %>%
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
          ) %>%
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_reactome <- NULL
    }

    if (length(object@merged_pathway_hmdb) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_hmdb <-
          object@merged_pathway_hmdb$module_result %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::mutate(database = "HMDB") %>%
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_kegg <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::mutate(database = "KEGG") %>%
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
          ) %>%
          dplyr::filter(!is.na(module_annotation))
      }
    } else{
      module_result_hmdb <- NULL
    }

    if (length(object@merged_pathway_metkegg) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_metkegg <-
          object@merged_pathway_metkegg$module_result %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::mutate(database = "KEGG") %>%
          dplyr::filter(!is.na(module_annotation))
      } else{
        module_result_kegg <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::mutate(database = "KEGG") %>%
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
          ) %>%
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

    ####module detection
    message("Identifying funcitonal modules...")

    node_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome,
            module_result_hmdb,
            module_result_metkegg)

    if (is.null(node_data)) {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "merge_modules()",
        parameter = list(
          sim.cutoff = sim.cutoff,
          measure_method = measure_method,
          cluster_method = cluster_method,
          path = path
        ),
        time = Sys.time()
      )

      process_info <-
        slot(object, "process_info")

      process_info$merge_modules <-
        parameter

      slot(object, "process_info") <-
        process_info

      message("Done")
      return(object)
    }

    functional_modules <-
      identify_functional_modules(
        query_type = query_type,
        sim_matrix = jaccard_index,
        variable_info = object@variable_info,
        node_data = node_data,
        sim.cutoff = sim.cutoff,
        cluster_method = cluster_method,
        hclust.method = hclust.method,
        save_to_local = save_to_local,
        path = file.path(path, "intermediate_data"),
        analysis_type = analysis_type
      )

    slot(object, "merged_module") <-
      functional_modules

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "merge_modules()",
      parameter = list(
        sim.cutoff = sim.cutoff,
        measure_method = measure_method,
        clustering_method = cluster_method,
        path = path
      ),
      time = Sys.time()
    )

    process_info <-
      slot(object, "process_info")

    process_info$merge_modules <-
      parameter

    slot(object, "process_info") <-
      process_info

    message("Done")
    object
}

#' Identify Functional Modules in a Similarity Matrix
#'
#' This function identifies functional modules from a similarity matrix and node data. It constructs a network graph, applies clustering algorithms, and generates module-level results, which can be saved locally. The function supports both pathway enrichment and gene set enrichment analysis (GSEA).
#'
#' @param query_type Character, the category of biological entity to query. Must be either "gene" or "metabolite".
#' @param variable_info A data frame containing mapping information between different ID types.
#'   For genes, should contain columns: "ensembl", "uniprot", and "entrezid".
#'   For metabolites, should contain columns: "hmdbid" and "keggid".
#' @param sim_matrix A data frame containing the similarity matrix with columns `name1`, `name2`, and `value`, representing the similarity between entities.
#' @param node_data A data frame containing node information, such as module assignments and other relevant attributes. This parameter cannot be `NULL`.
#' @param sim.cutoff A numerical value for similarity cutoff, default is 0.5. This is the similarity cutoff for pathways.
#'  Only edges with similarity above the cutoff will be stored in the graph data.
#' @param cluster_method Character, clustering method options:  "louvain" (default), "hierarchical",
#'   "binary_cut", "walktrap", "infomap", "edge_betweenness", "fast_greedy", "label_prop",
#'   "leading_eigen", "optimal".
#' @param hclust.method Character, agglomeration method for hierarchical clustering including:
#'   "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#'   Only used when `cluster_method = "hierarchical"`.
#' @param save_to_local Logical. Whether to save the resulting data to local files. Default is `FALSE`.
#' @param path Character. The directory path where the results will be saved, if `save_to_local = TRUE`. Default is `"."` (the current working directory).
#' @param analysis_type Character. Type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`. Default is `"enrich_pathway"`.
#'
#' @return A list containing:
#' \item{graph_data}{A `tidygraph` object representing the constructed network, including nodes and edges.}
#' \item{functional_module_result}{A data frame with the identified functional modules and their associated attributes, such as p-values and pathway descriptions.}
#' \item{result_with_module}{A data frame with the node data enriched with module information.}
#'
#' @details
#' The function first constructs a graph using the similarity matrix and node data, then applies clustering algorithms to identify functional modules. The function supports three clustering methods:
#' 1. Girvan-Newman: Community detection based on edge betweenness
#' 2. Binary cut: Uses the simplifyEnrichment binary_cut algorithm
#' 3. Hierarchical clustering: Groups pathways based on distance thresholds
#' The module information is further processed to summarize results such as the number of nodes and enriched pathways for each module. If `save_to_local` is `TRUE`, the results are saved as files in the specified `path`.
#'
#' @noRd

identify_functional_modules <-
  function(query_type = c("gene", "metabolite"),
           variable_info = NULL,
           sim_matrix,
           node_data,
           sim.cutoff = 0.5,
           cluster_method = "louvain",
           hclust.method = NULL,
           save_to_local = FALSE,
           path = ".",
           analysis_type = c("enrich_pathway", "do_gsea")) {

    if (missing(query_type)) {
      stop("query_type is required")
    }

    query_type <- match.arg(query_type)
    cluster_method <- match.arg(cluster_method)
    analysis_type <- match.arg(analysis_type)

    edge_data <-
      sim_matrix %>%
      dplyr::filter(value > sim.cutoff) %>%
      dplyr::rename(from = name1, to = name2, sim = value)

    rownames(edge_data) <- NULL

    if (is.null(node_data)) {
      return(NULL)
    }

    if ("node" %in% colnames(node_data)) {
      node_data <-
        node_data %>%
        dplyr::rename(pathway_id = node) %>%
        dplyr::select(module, dplyr::everything()) %>%
        dplyr::rename(node = module)
    } else {
      node_data <-
        node_data %>%
        dplyr::select(module, dplyr::everything()) %>%
        dplyr::rename(node = module)
    }

    ## Get clustering results based on selected method ====
    # Create similarity matrix for clustering methods that need it
    if (cluster_method %in% c("binary cut", "hierarchical")) {
      # Convert sim_matrix from edge list to matrix format
      sim_matrix_mat <- matrix(0, nrow = length(unique(node_data$node)),
                               ncol = length(unique(node_data$node)))
      rownames(sim_matrix_mat) <- colnames(sim_matrix_mat) <- unique(node_data$node)

      # Fill the matrix with similarity values
      for (i in 1:nrow(sim_matrix)) {
        name1 <- sim_matrix$name1[i]
        name2 <- sim_matrix$name2[i]
        value <- sim_matrix$value[i]
        if (name1 %in% rownames(sim_matrix_mat) && name2 %in% colnames(sim_matrix_mat)) {
          sim_matrix_mat[name1, name2] <- value
          sim_matrix_mat[name2, name1] <- value
        }
      }
      # Set diagonal to 1
      diag(sim_matrix_mat) <- 1
    }

    # cluster_result <- switch(cluster_method,
    #                          "binary cut" = merge_by_binary_cut(
    #                            sim_matrix = sim_matrix_mat,
    #                            sim.cutoff = sim.cutoff
    #                          ),
    #                          "girvan newman" = merge_by_Girvan_Newman(
    #                            edge_data = edge_data,
    #                            node_data = node_data,
    #                            sim.cutoff = sim.cutoff
    #                          ),
    #                          "hierarchical" = merge_by_hierarchical(
    #                            sim_matrix = sim_matrix_mat,
    #                            hclust.method = hclust.method,
    #                            sim.cutoff = sim.cutoff
    #                          ))

    cluster_result <-
      switch(cluster_method,
             "hierarchical" = {
               merge_by_hierarchical(sim_matrix = sim_matrix_mat,
                                     hclust.method = hclust.method,
                                     sim.cutoff = sim.cutoff)
             },
             "binary_cut" = {
               merge_by_binary_cut(sim_matrix = sim_matrix_mat,
                                   sim.cutoff = sim.cutoff)
             },
             # Network-based methods
             "louvain" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_louvain(graph_obj, weights = igraph::E(graph_obj)$sim)

               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "walktrap" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_walktrap(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "infomap" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_infomap(graph_obj, e.weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "edge_betweenness" = {
               # Use distance weights (1 - similarity) for edge betweenness
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ] |> dplyr::mutate(sim = 1 - sim)
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_edge_betweenness(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "fast_greedy" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_fast_greedy(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "label_prop" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_label_prop(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "leading_eigen" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_leading_eigen(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             },
             "optimal" = {
               filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
               if (nrow(filtered_edges) == 0) return(data.frame(node = node_data$node,
                                                                module = paste("Functional_module", "1", sep = "_")))

               graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                                          vertices = node_data)
               comm <- igraph::cluster_optimal(graph_obj, weights = igraph::E(graph_obj)$sim)
               data.frame(node = node_data$node,
                          module = paste("Functional_module", as.character(igraph::membership(comm)), sep = "_"))
             }
      )

    ## Create and update graph data with clustering result ====
    edge_data <-
      edge_data %>%
      dplyr::filter(sim > sim.cutoff)
    ## Create and update tidygraph object
    graph_data <-
      tidygraph::tbl_graph(nodes = if (analysis_type == "do_gsea") node_data[, -c(12)] else node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree()) %>%
      dplyr::left_join(cluster_result, by = "node")

    ###clustered different GO terms
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>%
      dplyr::arrange(module, p_adjust)

    ##add module content number
    module_content_number <-
      result_with_module %>%
      dplyr::count(module) %>%
      dplyr::rename(module_content_number = n)

    result_with_module <-
      result_with_module %>%
      dplyr::left_join(module_content_number, by = "module")

    if (save_to_local) {
      save(result_with_module, file = file.path(path, "result_with_module.RData"))
    }

    graph_data <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::left_join(module_content_number, by = "module")

    if (save_to_local) {
      save(graph_data, file = file.path(path, "graph_data.RData"))
    }

    if (query_type == "gene") {
      if (analysis_type == "enrich_pathway") {
        functional_module_result <-
          result_with_module %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            if (nrow(x) == 1) {
              x$module_content <-
                paste(x$node, collapse = ";")

              x$geneID <-
                x$geneID %>%
                stringr::str_split(pattern = "/") %>%
                unlist() %>%
                unify_id_internal(variable_info = variable_info,
                                  query_type = "gene") %>%
                unique() %>%
                paste(collapse = '/')

              x <-
                x %>%
                dplyr::select(module, everything()) %>%
                dplyr::distinct(module, .keep_all = TRUE) %>%
                dplyr::select(-node)
              return(x)
            }

            x =
              x %>%
              dplyr::arrange(p_adjust)

            x$module_content <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$Description, collapse = ";")

            x$BgRatio <-
              paste(x$BgRatio, collapse = ";")

            x$pvalue <- min(as.numeric(x$pvalue))
            x$p_adjust <- min(as.numeric(x$p_adjust))
            x$qvalue <- min(as.numeric(x$qvalue))
            x$geneID =
              x$geneID %>%
              stringr::str_split(pattern = "/") %>%
              unlist() %>%
              unify_id_internal(variable_info = variable_info,
                                query_type = "gene") %>%
              unique() %>%
              paste(collapse = '/')

            x$Count <-
              length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])

            x$pathway_id <-
              paste(x$pathway_id, collapse = ";")

            x <-
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE) %>%
              dplyr::select(-node)
            x
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
            Description, ";"
          ), `[`, 1))) %>%
          dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::select(module_annotation, everything())

        functional_module_result$module_annotation <-
          stringr::str_split(functional_module_result$Description, ";") %>%
          purrr::map(function(x) {
            x[1]
          }) %>%
          unlist()

        functional_module_result$Count <-
          as.numeric(functional_module_result$Count)
        functional_module_result$pvalue <-
          as.numeric(functional_module_result$pvalue)
        functional_module_result$degree <-
          as.numeric(functional_module_result$degree)
        functional_module_result$qvalue <-
          as.numeric(functional_module_result$qvalue)
      } else{
        functional_module_result <-
          result_with_module %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            if (nrow(x) == 1) {
              x$module_content <-
                paste(x$node, collapse = ";")
              x <-
                x %>%
                dplyr::select(module, everything()) %>%
                dplyr::distinct(module, .keep_all = TRUE) %>%
                dplyr::select(-node)
              x$geneID =
                x$core_enrichment %>%
                stringr::str_replace(";", "/") %>%
                stringr::str_split(pattern = "/") %>%
                unlist() %>%
                unify_id_internal(variable_info = variable_info,
                                  query_type = "gene") %>%
                unique() %>%
                paste(collapse = '/')
              # x$Count <-
              #   length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
              x <-
                x %>%
                dplyr::select(module, everything()) %>%
                dplyr::distinct(module, .keep_all = TRUE)
              return(x)
            }

            x <-
              x %>%
              dplyr::arrange(p_adjust)

            x$module_content <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$Description, collapse = ";")

            x$pvalue <- as.numeric(x$pvalue)[1]
            x$p_adjust <- as.numeric(x$p_adjust)[1]
            x$qvalue <- as.numeric(x$qvalue)[1]


            x$setSize <- as.numeric(x$setSize)[1]
            x$enrichmentScore <- as.numeric(x$enrichmentScore)[1]
            x$NES <- as.numeric(x$NES)[1]
            x$rank <- as.numeric(x$rank)[1]
            x$leading_edge <- x$leading_edge[1]
            x$core_enrichment <-
              paste0(unique(unlist(
                stringr::str_split(x$core_enrichment, pattern = "/")
              )), collapse = "/")
            x$degree <- as.numeric(x$degree)[1]

            x$geneID =
              x$core_enrichment %>%
              stringr::str_replace(";", "/") %>%
              stringr::str_split(pattern = "/") %>%
              unlist() %>%
              unify_id_internal(variable_info = variable_info,
                                query_type = "gene") %>%
              unique() %>%
              paste(collapse = '/')

            x$pathway_id <-
              paste(x$pathway_id, collapse = ";")

            x <-
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE) %>%
              dplyr::select(-node)
            x
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
            Description, ";"
          ), `[`, 1))) %>%
          dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::select(module_annotation, everything())

        functional_module_result$module_annotation <-
          stringr::str_split(functional_module_result$Description, ";") %>%
          purrr::map(function(x) {
            x[1]
          }) %>%
          unlist()

        functional_module_result$Count <-
          purrr::map(functional_module_result$geneID, function(x) {
            length(stringr::str_split(x, pattern = "/")[[1]])
          }) %>%
          unlist()

        functional_module_result$Count <-
          as.numeric(functional_module_result$Count)
        functional_module_result$pvalue <-
          as.numeric(functional_module_result$pvalue)
        functional_module_result$degree <-
          as.numeric(functional_module_result$degree)
        functional_module_result$qvalue <-
          as.numeric(functional_module_result$qvalue)
        functional_module_result$NES <-
          as.numeric(functional_module_result$NES)
      }
    } else if (query_type == "metabolite") {
      if (analysis_type == "enrich_pathway") {
        result_with_module <-
          result_with_module %>%
          dplyr::select(-c(all_id, all_number, mapped_number, mapped_percentage))

        functional_module_result <-
          result_with_module %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            # cat(unique(x$module), " ")
            if (nrow(x) == 1) {
              x$module_content <-
                paste(x$node, collapse = ";")

              x$mapped_id <-
                x$mapped_id %>%
                stringr::str_split(pattern = "/") %>%
                unlist() %>%
                unify_id_internal(variable_info = variable_info,
                                  query_type = "metabolite") %>%
                unique() %>%
                paste(collapse = '/')

              x$Count <-
                length(stringr::str_split(x$mapped_id[1], pattern = "/")[[1]])

              x <-
                x %>%
                dplyr::select(module, everything()) %>%
                dplyr::distinct(module, .keep_all = TRUE) %>%
                dplyr::select(-node)
              return(x)
            }

            x =
              x %>%
              dplyr::arrange(p_adjust)

            x$module_content <-
              paste(x$node, collapse = ";")

            x$pathway_id <-
              paste(x$pathway_id, collapse = ";")

            x$Description <-
              paste(x$Description, collapse = ";")

            x$p_value <- min(as.numeric(x$p_value))
            x$p_adjust <- min(as.numeric(x$p_adjust))

            x$mapped_id <-
              x$mapped_id %>%
              stringr::str_split(pattern = "/") %>%
              unlist() %>%
              unify_id_internal(variable_info = variable_info,
                                query_type = "metabolite") %>%
              unique() %>%
              paste(collapse = '/')

            x$Count <-
              length(stringr::str_split(x$mapped_id[1], pattern = "/")[[1]])

            x =
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE) %>%
              dplyr::select(-node)

            x

          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
            Description, ";"
          ), `[`, 1))) %>%
          dplyr::mutate(
            p_adjust = as.numeric(p_adjust),
            Count = as.numeric(Count),
            p_value = as.numeric(p_value),
            degree = as.numeric(degree)
          ) %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::select(module_annotation, everything())
      }
    }

    if (save_to_local) {
      save(functional_module_result,
           file = file.path(path, "functional_module_result.RData"))
    }

    list(
      graph_data = graph_data,
      functional_module_result = functional_module_result %>% dplyr::select(-database),
      result_with_module = result_with_module
    )

}
