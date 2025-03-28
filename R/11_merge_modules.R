# ####pathway enrichment
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")
#
# load("demo_data.rda")
# load("result/enriched_modules")
#
# enriched_functional_module <-
#   merge_modules(
#     object = enriched_modules,
#     sim.cutoff = 0.5,
#     measure_method = c("jaccard"),
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
# gsea_enriched_functional_module <-
#   merge_modules(
#     object = gsea_enriched_modules,
#     sim.cutoff = 0.5,
#     measure_method = c("jaccard"),
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_functional_module, file = "result/enriched_functional_module")

#' Identify Functional Modules Across Different Databases
#'
#' This function identify functional modules from various databases, such as GO, KEGG, and Reactome.
#' It calculates the similarity matrix between all the pathways and clusters them into functional modules.
#'
#' @param object An S4 object, expected to be processed by `merge_pathways()`.
#' @param sim.cutoff A numerical value for similarity cutoff, default is 0.5.
#' @param measure_method A character vector indicating the method for measuring similarity, default is "jaccard".
#' @param path A character string for the directory where to save results, default is "result".
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#'
#' @examples
#' \dontrun{
#' merge_modules(variable_info = my_var_info, object = my_object)
#' }
#'
#' @note Please ensure that `object` has been processed by `merge_pathways()`
#' before using this function.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @export
merge_modules <-
  function(object,
           sim.cutoff = 0.5,
           measure_method = c("jaccard"),
           path = "result",
           save_to_local = TRUE) {
    variable_info <- object@variable_info
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
    } else{
      analysis_type <- "do_gsea"
    }

    ######calculate the similarity (jaccard index) between all the pathways
    if (length(object@merged_pathway_go) != 0) {
      if (analysis_type == "enrich_pathway") {
        module_result_go <-
          object@merged_pathway_go$module_result %>%
          # dplyr::filter(ONTOLOGY != "CC") %>%
          dplyr::arrange(p.adjust) %>%
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
            p.adjust,
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
            p.adjust,
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
          dplyr::arrange(p.adjust) %>%
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
            p.adjust,
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
            p.adjust,
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
          dplyr::arrange(p.adjust) %>%
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
            p.adjust,
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
            p.adjust,
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

    message("Calculating the similarity matrix...")
    jaccard_index <-
      get_jaccard_index_for_three_databases(
        module_result_go = module_result_go,
        module_result_kegg = module_result_kegg,
        module_result_reactome = module_result_reactome,
        variable_info = variable_info,
        analysis_type = analysis_type
      )

    ####module detection
    message("Identifying funcitonal modules...")

    node_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome)

    if (is.null(node_data)) {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "merge_modules()",
        parameter = list(
          sim.cutoff = sim.cutoff,
          measure_method = measure_method,
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
        sim_matrix = jaccard_index,
        node_data = node_data,
        sim.cutoff = sim.cutoff,
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
#' @param sim_matrix A data frame containing the similarity matrix with columns `name1`, `name2`, and `value`, representing the similarity between entities.
#' @param node_data A data frame containing node information, such as module assignments and other relevant attributes. This parameter cannot be `NULL`.
#' @param sim.cutoff Numeric. A similarity cutoff value used to filter edges in the similarity matrix. Only edges with similarity values greater than `sim.cutoff` will be included in the analysis. Default is `0.5`.
#' @param save_to_local Logical. Whether to save the resulting data to local files. Default is `TRUE`.
#' @param path Character. The directory path where the results will be saved, if `save_to_local = TRUE`. Default is `"."` (the current working directory).
#' @param analysis_type Character. Type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`. Default is `"enrich_pathway"`.
#'
#' @return A list containing:
#' \item{graph_data}{A `tidygraph` object representing the constructed network, including nodes and edges.}
#' \item{functional_module_result}{A data frame with the identified functional modules and their associated attributes, such as p-values and pathway descriptions.}
#' \item{result_with_module}{A data frame with the node data enriched with module information.}
#'
#' @details
#' The function first constructs a graph using the similarity matrix and node data, then applies the edge betweenness clustering algorithm to identify functional modules. The module information is further processed to summarize results such as the number of nodes and enriched pathways for each module. If `save_to_local` is `TRUE`, the results are saved as files in the specified `path`.
#'
#' @examples
#' \dontrun{
#' sim_matrix <-
#'   data.frame(
#'     name1 = c("A", "B", "C"),
#'     name2 = c("B", "C", "A"),
#'     value = c(0.6, 0.7, 0.8)
#'   )
#' node_data <-
#'   data.frame(module = c("M1", "M2", "M3"),
#'              p.adjust = c(0.01, 0.05, 0.03))
#' result <-
#'   identify_functional_modules(sim_matrix,
#'                               node_data,
#'                               sim.cutoff = 0.5,
#'                               save_to_local = FALSE)
#' }
#'
#' @importFrom dplyr filter rename select mutate count arrange everything left_join distinct case_when
#' @importFrom tidygraph tbl_graph centrality_degree activate
#' @importFrom igraph cluster_edge_betweenness edge_attr membership upgrade_graph vertex_attr
#' @importFrom purrr map
#' @importFrom stringr str_split str_replace
#' @importFrom plyr dlply
#' @export



identify_functional_modules <-
  function(sim_matrix,
           node_data,
           sim.cutoff = 0.5,
           save_to_local = TRUE,
           path = ".",
           analysis_type = c("enrich_pathway", "do_gsea")) {
    analysis_type <- match.arg(analysis_type)

    edge_data <-
      sim_matrix %>%
      dplyr::filter(value > sim.cutoff) %>%
      dplyr::rename(from = name1, to = name2, sim = value)

    rownames(edge_data) <- NULL

    if (is.null(node_data)) {
      return(NULL)
    }

    node_data <-
      node_data %>%
      dplyr::select(module, dplyr::everything()) %>%
      dplyr::rename(node = module)

    graph_data <-
      tidygraph::tbl_graph(nodes = if (analysis_type == "do_gsea") node_data[, -c(12)] else node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())

    subnetwork <-
      suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))

    cluster <-
      paste("Functional_module",
            as.character(igraph::membership(subnetwork)),
            sep = "_")

    graph_data <-
      graph_data %>%
      igraph::upgrade_graph() %>%
      tidygraph::activate(what = "nodes") %>%
      tidygraph::mutate(module = cluster)

    ###clustered different GO terms
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(module, p.adjust)

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

    if (analysis_type == "enrich_pathway") {
      functional_module_result <-
        result_with_module %>%
        plyr::dlply(.variables = .(module)) %>%
        purrr::map(function(x) {
          # cat(unique(x$module), " ")
          if (nrow(x) == 1) {
            x$module_content <-
              paste(x$node, collapse = ";")
            x <-
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE) %>%
              dplyr::select(-node)
            return(x)
          }

          x =
            x %>%
            dplyr::arrange(p.adjust)

          x$module_content <-
            paste(x$node, collapse = ";")

          x$Description <-
            paste(x$Description, collapse = ";")

          x$BgRatio <-
            paste(x$BgRatio, collapse = ";")

          x$pvalue <- min(as.numeric(x$pvalue))
          x$p.adjust <- min(as.numeric(x$p.adjust))
          x$qvalue <- min(as.numeric(x$qvalue))
          x$geneID =
            x$geneID %>%
            stringr::str_split(pattern = "/") %>%
            unlist() %>%
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
        dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
        dplyr::arrange(p.adjust) %>%
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
              unique() %>%
              paste(collapse = '/')
            x$Count <-
              length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])
            x <-
              x %>%
              dplyr::select(module, everything()) %>%
              dplyr::distinct(module, .keep_all = TRUE)
            return(x)
          }

          x <-
            x %>%
            dplyr::arrange(p.adjust)

          x$module_content <-
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
          x$degree <- as.numeric(x$degree)[1]

          x$geneID =
            x$core_enrichment %>%
            stringr::str_replace(";", "/") %>%
            stringr::str_split(pattern = "/") %>%
            unlist() %>%
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
        dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::select(module_annotation, everything())

      functional_module_result$module_annotation <-
        stringr::str_split(functional_module_result$Description, ";") %>%
        purrr::map(function(x) {
          x[1]
        }) %>%
        unlist()

      functional_module_result$Count <-
        purrr::map(functional_module_result$core_enrichment, function(x) {
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
