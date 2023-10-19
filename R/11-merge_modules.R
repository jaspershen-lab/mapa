# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("demo_data.rda")
# load("result/modules")
#
# variable_info <-
#   extract_variable_info(demo_data)
#
# functional_module <-
#   merge_modules(
#     variable_info = variable_info,
#     object = modules,
#     sim.cutoff = 0.5,
#     measure_method = c("jaccard"),
#     path = "result",
#     save_to_local = TRUE
#   )
#
# save(functional_module, file = "result/functional_module")

#' Identify Functional Modules Across Different Databases
#'
#' This function identify functional modules from various databases, such as GO, KEGG, and Reactome.
#' It calculates the similarity matrix between all the pathways and clusters them into functional modules.
#'
#' @param variable_info A data frame or tibble containing variable information.
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
  function(variable_info,
           object,
           sim.cutoff = 0.5,
           measure_method = c("jaccard"),
           path = "result",
           save_to_local = FALSE) {
    if(save_to_local){
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

    ######calculate the similarity (jaccard index) between all the pathways
    if (length(object@merged_pathway_go) != 0) {
      module_result_go <-
        object@merged_pathway_go$module_result %>%
        dplyr::filter(ONTOLOGY != "CC") %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "GO") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
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
      module_result_go <- NULL
    }

    if (length(object@merged_pathway_kegg) != 0) {
      module_result_kegg <-
        object@merged_pathway_kegg$module_result %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "KEGG") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
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
      module_result_kegg <- NULL
    }

    if (length(object@merged_pathway_reactome) != 0) {
      module_result_reactome <-
        object@merged_pathway_reactome$module_result %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "Reactome") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
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
      module_result_reactome <- NULL
    }

    message("Calculating the similarity matrix...")
    jaccard_index <-
      get_jaccard_index_for_three_databases(
        module_result_go = module_result_go,
        module_result_kegg = module_result_kegg,
        module_result_reactome = module_result_reactome,
        variable_info = variable_info
      )

    edge_data =
      jaccard_index %>%
      dplyr::filter(value > sim.cutoff) %>%
      dplyr::rename(from = name1, to = name2, sim = value)

    node_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome) %>%
      dplyr::select(module, dplyr::everything()) %>%
      dplyr::rename(node = module)

    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())

    subnetwork <-
      igraph::cluster_edge_betweenness(graph = graph_data,
                                       weights = abs(edge_attr(graph_data,
                                                               "sim")))
    cluster <-
      paste("Functional_module",
            as.character(igraph::membership(subnetwork)),
            sep = "_")

    graph_data <-
      graph_data %>%
      mutate(module = cluster)

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

    if(save_to_local){
      save(result_with_module,
           file = file.path(path, "intermediate_data/result_with_module"))
    }

    graph_data <-
      graph_data %>%
      activate(what = "nodes") %>%
      dplyr::left_join(module_content_number, by = "module")

    if(save_to_local){
      save(graph_data, file = file.path(path, "intermediate_data/graph_data"))
    }

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

        x <-
          x %>%
          dplyr::select(module, everything()) %>%
          dplyr::distinct(module, .keep_all = TRUE) %>%
          dplyr::select(-node)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(
        module_annotation = case_when(
          module == "Other" ~ Description,
          module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
        )
      ) %>%
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

    if(save_to_local){
      save(
        functional_module_result,
        file = file.path(path, "intermediate_data/functional_module_result")
      )
    }

    slot(object, "merged_module") <-
      list(
        graph_data = graph_data,
        functional_module_result = functional_module_result,
        result_with_module = result_with_module
      )

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

    message("\nDone")
    object
  }
