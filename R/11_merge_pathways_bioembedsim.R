# setwd(r4projects::get_project_wd())
# library(mapa)
# source("R/6_utils.R")

# For gene
# load("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/ora_openai_semantic_sim_matrix.rda")
# enriched_functional_module <-
#   merge_pathways_bioembedsim(
#     object = openai_semantic_sim_matrix,
#     # cluster_method = "edge_betweenness",
#     sim.cutoff = 0.5,
#     save_to_local = FALSE
# )

# For metabolite
# object <- openai_sim_matrix_met
# enriched_functional_module <-
#   merge_pathways_bioembedsim(
#     object = openai_sim_matrix_met,
#     sim.cutoff = 0.5
#   )

#' Merge Pathways based on biotext embedding similarity
#'
#' @description
#' Identifies functional modules by clustering pathways based on their biotext embedding similarity.
#' This function groups related pathways into functional modules using various clustering
#' methods including binary cut, graph-based clustering, or hierarchical clustering.
#'
#' @param object An object containing functional enrichment result and similarity matrix.
#'   The object must have `enriched_pathway` and `sim_matrix` components.
#' @param sim.cutoff A numerical value for similarity cutoff, default is 0.5. This is the similarity cutoff for pathways.
#'  Only edges with similarity above the cutoff will be stored in the graph data.
#' @param cluster_method Character, clustering method options:
#'        \itemize{
#'          \item **Hierarchical** ‒ supply `"h_<agglom.method>"`, where
#'            `<agglom.method>` is one of
#'            `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
#'            `"mcquitty"`, `"median"`, `"centroid"`.
#'          \item `"binary_cut"`
#'          \item Graph-based: `"louvain"`, `"walktrap"`, `"infomap"`,
#'            `"edge_betweenness"`, `"fast_greedy"`, `"label_prop"`,
#'            `"leading_eigen"`, `"optimal"`
#'        }
#' @param save_to_local Logical, whether to save intermediate results locally (default: FALSE).
#' @param path Character, directory to save results when save_to_local is TRUE (default: "result").
#'
#' @return Returns the updated enriched_pathway object with functional modules added.
#'   The object will contain a new slot "merged_module" with:
#'   \itemize{
#'     \item graph_data: A tidygraph object representing the network of related pathways
#'     \item functional_module_result: Data frame with module information and annotations
#'     \item result_with_module: Data frame with raw results assigned to modules
#'   }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming 'result' is an object with enriched pathways and similarity matrix
#' enriched_with_modules <- merge_pathways_bioembedsim(
#'   object = result,
#'   sim.cutoff = 0.6,
#'   cluster_method = "binary_cut"
#' )
#'
#' # Using hierarchical clustering with complete linkage
#' enriched_with_modules <- merge_pathways_bioembedsim(
#'   object = result,
#'   sim.cutoff = 0.5,
#'   cluster_method = "h_ward.D"
#' )
#' }
#'
#' @importFrom dplyr filter rename mutate select left_join arrange count distinct
#' @importFrom tidygraph tbl_graph activate
#' @importFrom igraph vertex_attr cluster_edge_betweenness edge_attr membership
#' @importFrom simplifyEnrichment binary_cut
#' @importFrom stats hclust cutree as.dist
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom plyr dlply
#'
#' @export

merge_pathways_bioembedsim <-
  function(object,
           sim.cutoff = 0.5,
           cluster_method = "louvain",
           save_to_local = FALSE,
           path = "result") {

    available_methods <- c(
      "h_ward.D", "h_ward.D2", "h_single", "h_complete",
      "h_average", "h_mcquitty", "h_median", "h_centroid",
      "binary_cut", "louvain", "walktrap", "infomap",
      "edge_betweenness", "fast_greedy", "label_prop", "leading_eigen",
      "optimal"
    )
    if (!(cluster_method %in% available_methods)) {
      invalid_methods <- cluster_method
      stop(paste(
        "Invalid methods:",
        paste(invalid_methods, collapse = ", "),
        "\nAvailable methods:",
        paste(available_methods, collapse = ", ")
      ))
    }
    if (grepl("^h_", cluster_method)) {
      hclust_method <- gsub("^h_", "", cluster_method[grepl("^h_", cluster_method)])
      cluster_method <- c("hierarchical")
    }

    message("Identifying funcitonal modules...")

    ## Collect edge data ====
    parameters <- object$enriched_pathway@process_info$merge_pathways@parameter
    query_type <- parameters$query_type

    edge_data <-
      as.data.frame.table(object$sim_matrix, responseName = "sim") %>%
      dplyr::filter(Var1 != Var2) %>%                 # Remove selfs-edges
      dplyr::rename(from = Var1, to = Var2) %>%
      dplyr::mutate(across(c(from, to), as.character)) %>%
      dplyr::filter(from < to)

    ## Collect node data
    result <- data.frame()

    if (query_type == "gene") {
      if (!is.null(object$enriched_pathway@enrichment_go_result)) {
        result <-
          object$enriched_pathway@enrichment_go_result@result %>%
          dplyr::select(-ONTOLOGY) %>%
          dplyr::filter(p_adjust < parameters$p.adjust.cutoff.go) %>%
          dplyr::filter(Count > parameters$count.cutoff.go) %>%
          dplyr::mutate(database = "GO") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_kegg_result)) {
        result <-
          object$enriched_pathway@enrichment_kegg_result@result %>%
          {if ("enrich_pathway" %in% names(object$enriched_pathway@process_info))
            dplyr::select(., -c(category, subcategory)) else .} %>%
          dplyr::filter(.data$p_adjust < parameters$p.adjust.cutoff.kegg,
                 .data$Count > parameters$count.cutoff.kegg) %>%
          dplyr::mutate(database = "KEGG") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_reactome_result)) {
        result <-
          object$enriched_pathway@enrichment_reactome_result@result %>%
          dplyr::filter(p_adjust < parameters$p.adjust.cutoff.reactome) %>%
          dplyr::filter(Count > parameters$count.cutoff.reactome) %>%
          dplyr::mutate(database = "Reactome") %>%
          rbind(result)
      }

      node_data <- result %>% dplyr::rename(node = ID)

    } else if (query_type == "metabolite") {
      if (!is.null(object$enriched_pathway@enrichment_hmdb_result)) {
        result <-
          object$enriched_pathway@enrichment_hmdb_result@result %>%
          dplyr::filter(p_adjust < parameters$p.adjust.cutoff.hmdb) %>%
          dplyr::filter(mapped_number > parameters$count.cutoff.hmdb) %>%
          dplyr::mutate(database = "HMDB") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_metkegg_result)) {
        result <-
          object$enriched_pathway@enrichment_metkegg_result@result %>%
          dplyr::filter(p_adjust < parameters$p.adjust.cutoff.metkegg) %>%
          dplyr::filter(mapped_number > parameters$count.cutoff.metkegg) %>%
          dplyr::mutate(database = "KEGG") %>%
          rbind(result)
      }

      node_data <- result %>% dplyr::rename(node = pathway_id)
    }

    ## Get clustering results ====
    # cluster_result <- switch(cluster_method,
    #                          "binary cut" = merge_by_binary_cut(
    #                            sim_matrix = object$sim_matrix,
    #                            sim.cutoff = sim.cutoff
    #                          ),
    #                          "girvan newman" = merge_by_Girvan_Newman(
    #                            edge_data = edge_data,
    #                            node_data = node_data,
    #                            sim.cutoff = sim.cutoff
    #                          ),
    #                          "hierarchical" = merge_by_hierarchical(
    #                            sim_matrix = object$sim_matrix,
    #                            hclust.method = hclust.method,
    #                            sim.cutoff = sim.cutoff
    #                          ))
    cluster_result <-
      switch(cluster_method,
             "hierarchical" = {
               merge_by_hierarchical(sim_matrix = object$sim_matrix,
                                     hclust.method = hclust_method,
                                     sim.cutoff = sim.cutoff)
             },
             "binary_cut" = {
               merge_by_binary_cut(sim_matrix = object$sim_matrix,
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
      edge_data |>
      dplyr::filter(sim > sim.cutoff)
    ## Create and update tidygraph object
    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE,
                           node_key = "node") |>
      dplyr::mutate(degree = tidygraph::centrality_degree()) |>
      dplyr::left_join(cluster_result)

    ## Get result with module ====
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      # # Check which column exists and standardize the name
      # dplyr::mutate(
      #   p.adjust = case_when(
      #     "p.adjust" %in% names(.) ~ as.numeric(`p.adjust`),
      #     "p_value_adjust" %in% names(.) ~ as.numeric(`p_value_adjust`)
      #   )
      # ) %>%
      # Process based on which column exists
      {
        df <- .
        if("p_adjust" %in% names(df)) {
          df$p_adjust <- as.numeric(as.character(df$p_adjust))
        }
        df
      } %>%
      dplyr::arrange(module, p_adjust)

    ### Add module content number
    module_content_number <-
      result_with_module %>%
      dplyr::count(module) %>%
      dplyr::rename(module_content_number = n)

    result_with_module <-
      result_with_module %>%
      dplyr::left_join(module_content_number, by = "module")

    graph_data <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::left_join(module_content_number, by = "module")

    ## Get functional module result ====
    if (query_type == "gene") {
      if ("enrich_pathway" %in% names(object$enriched_pathway@process_info)) {
        functional_module_result <-
          result_with_module %>%
          dplyr::select(-database) %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            # cat(unique(x$module), " ")
            if (nrow(x) == 1) {
              x$geneID <-
                x$geneID %>%
                stringr::str_split(pattern = "/") %>%
                unlist() %>%
                unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                  query_type = "gene") %>%
                unique() %>%
                paste(collapse = '/')

              x$Count <-
                length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])

              return(x)
            }

            x =
              x %>%
              dplyr::arrange(p_adjust)

            x$node <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$Description, collapse = ";")

            x$BgRatio <-
              paste(x$BgRatio, collapse = ";")

            x$pvalue <- min(as.numeric(x$pvalue))
            x$p_adjust <- min(as.numeric(x$p_adjust))
            x$qvalue <- min(as.numeric(x$qvalue))
            x$geneID <-
              x$geneID %>%
              stringr::str_split(pattern = "/") %>%
              unlist() %>%
              unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                query_type = "gene") %>%
              unique() %>%
              paste(collapse = '/')

            x$Count <-
              length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])

            x =
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE)

            x

          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
            Description, ";"
          ), `[`, 1))) %>%
          dplyr::mutate(
            pvalue = as.numeric(pvalue),
            Count = as.numeric(Count),
            p_adjust = as.numeric(p_adjust),
            qvalue = as.numeric(qvalue),
            degree = as.numeric(degree)
          ) %>%
          dplyr::arrange(p_adjust) %>%
          dplyr::select(module_annotation, everything())
      } else{
        functional_module_result <-
          result_with_module %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            # cat(unique(x$module), " ")
            if (nrow(x) == 1) {
              x$geneID <-
                x$core_enrichment %>%
                stringr::str_split(pattern = "/") %>%
                unlist() %>%
                unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                  query_type = "gene") %>%
                unique() %>%
                paste(collapse = '/')

              return(x)
            }

            x <-
              x %>%
              dplyr::mutate(NES = as.numeric(NES)) %>%
              dplyr::arrange(dplyr::desc(abs(NES)))

            x$node <-
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

            x$geneID <-
              x$core_enrichment %>%
              stringr::str_split(pattern = "/") %>%
              unlist() %>%
              unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                query_type = "gene") %>%
              unique() %>%
              paste(collapse = '/')

            x$degree <- as.numeric(x$degree)[1]

            x <-
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE)

            x

          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
            Description, ";"
          ), `[`, 1))) %>%
          dplyr::mutate(
            p_adjust = as.numeric(p_adjust),
            NES = as.numeric(NES),
            qvalue = as.numeric(qvalue),
            degree = as.numeric(degree)
          ) %>%
          dplyr::mutate() %>%
          dplyr::arrange(dplyr::desc(abs(NES))) %>%
          dplyr::select(module_annotation, everything())

        functional_module_result$Count <-
          purrr::map(functional_module_result$geneID, function(x) {
            length(stringr::str_split(x, pattern = "/")[[1]])
          }) %>%
          unlist()
      }
    } else if (query_type == "metabolite") {
      if ("enrich_pathway" %in% names(object$enriched_pathway@process_info)) {
        result_with_module <-
          result_with_module %>%
          dplyr::rename(Count = mapped_number) %>%
          dplyr::select(-c(all_id, all_number, mapped_percentage))

        functional_module_result <-
          result_with_module %>%
          dplyr::select(-database) %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            # cat(unique(x$module), " ")
            if (nrow(x) == 1) {
              x <- x %>%
                dplyr::mutate(Description = pathway_name)

              # x$mapped_id =
              #   x$mapped_id %>%
              #   stringr::str_split(pattern = ";") %>%
              #   unlist() %>%
              #   unify_id_internal(variable_info = object$enriched_pathway@variable_info,
              #                     query_type = "metabolite") %>%
              #   unique() %>%
              #   paste(collapse = '/')
              x$mapped_id =
                x$mapped_id %>%
                stringr::str_split(pattern = ";") %>%
                unlist() %>%
                {
                  if (length(object$enriched_pathway@merged_pathway_hmdb) > 0) {
                    unify_id_internal(., variable_info = object$enriched_pathway@variable_info,
                                      query_type = "metabolite")
                  } else {
                    .
                  }
                } %>%
                unique() %>%
                paste(collapse = '/')

              x$Count <-
                length(stringr::str_split(x$mapped_id, pattern = "/")[[1]])

              return(x)
            }

            x =
              x %>%
              dplyr::arrange(p_adjust)

            x$node <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$pathway_name, collapse = ";")

            x$p_value <- min(as.numeric(x$p_value))
            x$p_adjust <- min(as.numeric(x$p_adjust))

            # x$mapped_id =
            #   x$mapped_id %>%
            #   stringr::str_split(pattern = ";") %>%
            #   unlist() %>%
            #   unify_id_internal(variable_info = object$enriched_pathway@variable_info,
            #                     query_type = "metabolite") %>%
            #   unique() %>%
            #   paste(collapse = '/')
            x$mapped_id =
              x$mapped_id %>%
              stringr::str_split(pattern = ";") %>%
              unlist() %>%
              {
                if (length(object$enriched_pathway@merged_pathway_hmdb) > 0) {
                  unify_id_internal(., variable_info = object$enriched_pathway@variable_info,
                                    query_type = "metabolite")
                } else {
                  .
                }
              } %>%
              unique() %>%
              paste(collapse = '/')

            x$Count <-
              length(stringr::str_split(x$mapped_id[1], pattern = "/")[[1]])

            x =
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE)

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

    functional_module_result$module_annotation <-
      stringr::str_split(functional_module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()

    if (query_type == "metabolite") {
      functional_module_result <-
        functional_module_result %>%
        dplyr::select(-c(describtion))
    }

    ## Collect all results (graph_data, functional_module_result, result_with_module) and parameters ====
    functional_modules <- list(
      graph_data = graph_data,
      functional_module_result = functional_module_result,
      result_with_module = result_with_module
    )

    slot(object$enriched_pathway, "merged_module") <-
      functional_modules

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "get_bioembedsim() and merge_pathways_bioembedsim()",
      parameter = list(
        sim.cutoff = sim.cutoff,
        clustering_method = cluster_method
      ),
      time = Sys.time()
    )

    process_info <-
      slot(object$enriched_pathway, "process_info")

    process_info$merge_modules <-
      parameter

    slot(object$enriched_pathway, "process_info") <-
      process_info

    message("Done")

    ## Save clustering result
    if (save_to_local) {
      path <- file.path(path, "functional_modules")
      dir.create(path, showWarnings = FALSE, recursive = TRUE)

      save(graph_data,
           file = file.path(path, "graph_data.RData"))
      save(functional_module_result,
           file = file.path(path, "functional_module_result.RData"))
      save(result_with_module,
           file = file.path(path, "result_with_module.RData"))
    }

    return(object$enriched_pathway)
}
