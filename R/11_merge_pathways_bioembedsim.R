# object <- openai_semantic_sim_matrix
# sim_matrix <- object$sim_matrix
# sim.cutoff <- 0.8
# hclust.method <- "complete"

# enriched_functional_module <-
#   merge_pathways_bioembedsim(
#     object = openai_semantic_sim_matrix,
#     cluster_method = "girvan newman",
#     sim.cutoff = 0.5,
#     save_to_local = TRUE,
#     path = "~/Desktop/result"
#   )

# For metabolite
# object <- openai_sim_matrix_met
# enriched_functional_module <-
#   merge_pathways_bioembedsim(
#     object = object,
#     cluster_method = "girvan newman",
#     sim.cutoff = 0.5
#   )

#' Merge Pathways based on biotext embedding similarity
#'
#' @description
#' Identifies functional modules by clustering pathways based on their biotext embedding similarity.
#' This function groups related pathways into functional modules using various clustering
#' methods including binary cut, Girvan-Newman, or hierarchical clustering.
#'
#' @param object An object containing functional enrichment result and similarity matrix.
#'   The object must have `enriched_pathway` and `sim_matrix` components.
#' @param sim.cutoff Numeric, similarity cutoff for clustering (default: 0.5).
#'   For binary cut and hierarchical clustering, sim.cutoff is the cutoff for cutting the dendrogram.
#'   For Girvan-Newman, sim.cutoff is the similarity cutoff for pathways. Only edges with similarity above the cutoff will be stored in the graph data.

#' @param cluster_method Character, method for clustering pathways. Options are:
#'   \itemize{
#'     \item "binary cut": Uses simplifyEnrichment binary_cut algorithm
#'     \item "girvan newman": Uses Girvan-Newman community detection
#'     \item "hierarchical": Uses hierarchical clustering
#'   }
#' @param hclust.method Character, the agglomeration method for hierarchical clustering.
#'   Only used when cluster_method is "hierarchical". See \code{\link[stats]{hclust}} for options.
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
#' @details
#' The function takes biotext embedding similarity scores between pathways
#' and groups related pathways into functional modules.
#' It supports three clustering methods:
#'
#' 1. Binary cut: Uses the simplifyEnrichment binary_cut algorithm
#' 2. Girvan-Newman: Community detection based on edge betweenness
#' 3. Hierarchical clustering: Groups pathways based on distance thresholds
#'
#' The function handles GO terms, KEGG pathways, and Reactome pathways, combining
#' them into unified functional modules based on similarity.
#'
#' @examples
#' \dontrun{
#' # Assuming 'result' is an object with enriched pathways and similarity matrix
#' enriched_with_modules <- merge_pathways_bioembedsim(
#'   object = result,
#'   sim.cutoff = 0.6,
#'   cluster_method = "binary cut"
#' )
#'
#' # Using hierarchical clustering with complete linkage
#' enriched_with_modules <- merge_pathways_bioembedsim(
#'   object = result,
#'   sim.cutoff = 0.5,
#'   cluster_method = "hierarchical",
#'   hclust.method = "complete"
#' )
#'
#' # Save results to local directory
#' enriched_with_modules <- merge_pathways_bioembedsim(
#'   object = result,
#'   cluster_method = "girvan newman",
#'   save_to_local = TRUE,
#'   path = "project/results"
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
           cluster_method = c("binary cut", "girvan newman", "hierarchical"),
           hclust.method = NULL,
           save_to_local = FALSE,
           path = "result") {

    cluster_method <- match.arg(cluster_method)

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
          dplyr::filter(p.adjust < parameters$p.adjust.cutoff.go) %>%
          dplyr::filter(Count > parameters$count.cutoff.go) %>%
          dplyr::mutate(database = "GO") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_kegg_result)) {
        result <-
          object$enriched_pathway@enrichment_kegg_result@result %>%
          dplyr::select(-c(category, subcategory)) %>%
          dplyr::filter(p.adjust < parameters$p.adjust.cutoff.kegg) %>%
          dplyr::filter(Count > parameters$count.cutoff.kegg) %>%
          dplyr::mutate(database = "KEGG") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_reactome_result)) {
        result <-
          object$enriched_pathway@enrichment_reactome_result@result %>%
          dplyr::filter(p.adjust < parameters$p.adjust.cutoff.reactome) %>%
          dplyr::filter(Count > parameters$count.cutoff.reactome) %>%
          dplyr::mutate(database = "Reactome") %>%
          rbind(result)
      }

      node_data <- result %>% dplyr::rename(node = ID)

    } else if (query_type == "metabolite") {
      if (!is.null(object$enriched_pathway@enrichment_hmdb_result)) {
        result <-
          object$enriched_pathway@enrichment_hmdb_result@result %>%
          dplyr::filter(p_value_adjust < parameters$p.adjust.cutoff.hmdb) %>%
          dplyr::filter(mapped_number > parameters$count.cutoff.hmdb) %>%
          dplyr::mutate(database = "HMDB") %>%
          rbind(result)
      }

      if (!is.null(object$enriched_pathway@enrichment_metkegg_result)) {
        result <-
          object$enriched_pathway@enrichment_metkegg_result@result %>%
          dplyr::filter(p_value_adjust < parameters$p.adjust.cutoff.metkegg) %>%
          dplyr::filter(mapped_number > parameters$count.cutoff.metkegg) %>%
          dplyr::mutate(database = "KEGG") %>%
          rbind(result)
      }

      node_data <- result %>% dplyr::rename(node = pathway_id)
    }

    ## Get clustering results ====
    cluster_result <- switch(cluster_method,
                             "binary cut" = merge_by_binary_cut(
                               sim_matrix = object$sim_matrix,
                               sim.cutoff = sim.cutoff
                             ),
                             "girvan newman" = merge_by_Girvan_Newman(
                               edge_data = edge_data,
                               node_data = node_data,
                               sim.cutoff = sim.cutoff
                             ),
                             "hierarchical" = merge_by_hierarchical(
                               sim_matrix = object$sim_matrix,
                               hclust.method = hclust.method,
                               sim.cutoff = sim.cutoff
                             ))

    ## Create and update graph data with clustering result ====
    if (cluster_method == "girvan newman") {
      ## Filter graph data according to sim.cutoff
      edge_data <-
        edge_data %>%
        dplyr::filter(sim > sim.cutoff)
      ## Create and update tidygraph object
      graph_data <-
        tidygraph::tbl_graph(nodes = node_data,
                             edges = edge_data,
                             directed = FALSE,
                             node_key = "node") %>%
        dplyr::mutate(degree = tidygraph::centrality_degree()) %>%
        dplyr::left_join(cluster_result)
    } else {
      ## Create and update tidygraph object
      graph_data <-
        tidygraph::tbl_graph(nodes = node_data,
                             edges = edge_data,
                             directed = FALSE,
                             node_key = "node") %>%
        dplyr::left_join(cluster_result) %>%
        activate(edges) %>%
        filter(
          ## Get the module attribute for the from node
          tidygraph::.N()$module[from] == tidygraph::.N()$module[to]
        ) %>%
        activate(nodes) %>%
        dplyr::mutate(degree = tidygraph::centrality_degree())
    }

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
        if("p.adjust" %in% names(df)) {
          df$p.adjust <- as.numeric(as.character(df$p.adjust))
        } else if("p_value_adjust" %in% names(df)) {
          df$p.adjust <- as.numeric(as.character(df$p_value_adjust))
          df$p_value_adjust <- NULL  # Remove original column
        }
        df
      } %>%
      dplyr::arrange(module, p.adjust)

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
                x$core_enrichment %>%
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
              dplyr::arrange(p.adjust)

            x$node <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$Description, collapse = ";")

            x$BgRatio <-
              paste(x$BgRatio, collapse = ";")

            x$pvalue <- min(as.numeric(x$pvalue))
            x$p.adjust <- min(as.numeric(x$p.adjust))
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
            p.adjust = as.numeric(p.adjust),
            qvalue = as.numeric(qvalue),
            degree = as.numeric(degree)
          ) %>%
          dplyr::arrange(p.adjust) %>%
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
            x$p.adjust <- as.numeric(x$p.adjust)[1]
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
            p.adjust = as.numeric(p.adjust),
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

              x$mapped_id =
                x$mapped_id %>%
                stringr::str_split(pattern = ";") %>%
                unlist() %>%
                unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                  query_type = "metabolite") %>%
                unique() %>%
                paste(collapse = '/')

              x$Count <-
                length(stringr::str_split(x$mapped_id, pattern = "/")[[1]])

              return(x)
            }

            x =
              x %>%
              dplyr::arrange(p.adjust)

            x$node <-
              paste(x$node, collapse = ";")

            x$Description <-
              paste(x$pathway_name, collapse = ";")

            x$p_value <- min(as.numeric(x$p_value))
            x$p.adjust <- min(as.numeric(x$p.adjust))

            x$mapped_id =
              x$mapped_id %>%
              stringr::str_split(pattern = ";") %>%
              unlist() %>%
              unify_id_internal(variable_info = object$enriched_pathway@variable_info,
                                query_type = "metabolite") %>%
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
            p.adjust = as.numeric(p.adjust),
            Count = as.numeric(Count),
            p_value = as.numeric(p_value),
            degree = as.numeric(degree)
          ) %>%
          dplyr::arrange(p.adjust) %>%
          dplyr::select(module_annotation, everything())
      }
    }

    functional_module_result$module_annotation <-
      stringr::str_split(functional_module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()
    functional_module_result <-
      functional_module_result %>%
      dplyr::select(-c(describtion))

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
        measure_method = cluster_method
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

## Internal functions to get clustering results
merge_by_binary_cut <- function(sim_matrix,
                                sim.cutoff) {
  clusters <- simplifyEnrichment::binary_cut(mat = sim_matrix, cutoff = sim.cutoff)
  cluster_result <-
    data.frame(node = rownames(sim_matrix),
               module = paste("Functional_module", as.character(clusters), sep = "_"))

  return(cluster_result)
}

merge_by_Girvan_Newman <- function(edge_data,
                                   node_data,
                                   sim.cutoff) {
  ## Filter graph data according to sim.cutoff
  edge_data <-
    edge_data %>%
    dplyr::filter(sim > sim.cutoff)
  ## Create tidygraph object
  graph_data <-
    tidygraph::tbl_graph(nodes = node_data,
                         edges = edge_data,
                         directed = FALSE,
                         node_key = "node") %>%
    dplyr::mutate(degree = tidygraph::centrality_degree())

  ## Perform clustering
  subnetwork <-
    suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))
  ## Assign functional module label for pathways
  cluster_result <-
    data.frame(node = node_data$node,
               module = paste("Functional_module", as.character(igraph::membership(subnetwork)), sep = "_"))

  return(cluster_result)
}

merge_by_hierarchical <- function(sim_matrix,
                                  hclust.method,
                                  sim.cutoff) {
  cosine_dist <- 1 - sim_matrix
  ## Convert distance matrix to a 'dist' object
  cosine_dist_obj <- as.dist(cosine_dist)
  ## Perform hierarchical clustering
  hc <- hclust(cosine_dist_obj, method = hclust.method)

  clusters <- cutree(hc, h = sim.cutoff)
  cluster_result <-
    data.frame(node = hc$labels,
               module = paste("Functional_module", as.character(clusters), sep = "_"))

  return(cluster_result)
}
