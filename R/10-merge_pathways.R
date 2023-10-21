# setwd(r4projects::get_project_wd())
# source("R/6-utils.R")
# source("R/8-functional_module_class.R")
# setwd("demo_data/")
# library(ggraph)
# library(igraph)
# library(tidygraph)
# library(tidyverse)
# library(extrafont)
# library(simplifyEnrichment)
# library(GOSim)
# library(plyr)
#
# load("result/enriched_pathways")
#
# object = enriched_pathways
# p.adjust.cutoff.go = 0.05
# p.adjust.cutoff.kegg = 0.05
# p.adjust.cutoff.reactome = 0.05
# count.cutoff.go = 5
# count.cutoff.kegg = 5
# count.cutoff.reactome = 5
# sim.cutoff.go = 0.5
# sim.cutoff.kegg = 0.5
# sim.cutoff.reactome = 0.5
# measure.method.go = "Wang"
# measure.method.kegg = "jaccard"
# measure.method.reactome = "jaccard"
# path = "result"
# save_to_local = FALSE
#
#
# enriched_modules <-
#   merge_pathways(
#     object = enriched_pathways,
#     p.adjust.cutoff.go = 0.05,
#     p.adjust.cutoff.kegg = 0.05,
#     p.adjust.cutoff.reactome = 0.05,
#     count.cutoff.go = 5,
#     count.cutoff.kegg = 5,
#     count.cutoff.reactome = 5,
#     sim.cutoff.go = 0.5,
#     sim.cutoff.kegg = 0.5,
#     sim.cutoff.reactome = 0.5,
#     measure.method.go = "Wang",
#     measure.method.kegg = "jaccard",
#     measure.method.reactome = "jaccard",
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_modules, file = "result/enriched_modules")


#' Merge Pathways from Multiple Databases
#'
#' This function merges pathway enrichment results from multiple databases (GO, KEGG, and Reactome)
#' into a single object. The function takes an object of class "functional_module" and various
#' parameters to filter and compute similarity measures.
#'
#' @param object An object of class "functional_module", typically a result from enrich_pathway function.
#' @param p.adjust.cutoff.go Adjusted p-value cutoff for GO database. Default is 0.05.
#' @param p.adjust.cutoff.kegg Adjusted p-value cutoff for KEGG database. Default is 0.05.
#' @param p.adjust.cutoff.reactome Adjusted p-value cutoff for Reactome database. Default is 0.05.
#' @param count.cutoff.go Count cutoff for GO database. Default is 5.
#' @param count.cutoff.kegg Count cutoff for KEGG database. Default is 5.
#' @param count.cutoff.reactome Count cutoff for Reactome database. Default is 5.
#' @param sim.cutoff.go Similarity cutoff for GO database. Default is 0.5.
#' @param sim.cutoff.kegg Similarity cutoff for KEGG database. Default is 0.5.
#' @param sim.cutoff.reactome Similarity cutoff for Reactome database. Default is 0.5.
#' @param measure.method.go A character vector specifying the similarity measure method for GO. Choices are "Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS", "jaccard". Default is "Wang".
#' @param measure.method.kegg A character vector specifying the similarity measure method for KEGG. Default is "jaccard".
#' @param measure.method.reactome A character vector specifying the similarity measure method for Reactome. Default is "jaccard".
#' @param path Directory path to save the results. Default is "result".
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#'
#' @return An object of class "functional_module" with slots for merged pathways from each database.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @export
#'

merge_pathways <-
  function(object,
           p.adjust.cutoff.go = 0.05,
           p.adjust.cutoff.kegg = 0.05,
           p.adjust.cutoff.reactome = 0.05,
           count.cutoff.go = 5,
           count.cutoff.kegg = 5,
           count.cutoff.reactome = 5,
           sim.cutoff.go = 0.5,
           sim.cutoff.kegg = 0.5,
           sim.cutoff.reactome = 0.5,
           measure.method.go = c("Wang", "Resnik",
                                 "Rel", "Jiang",
                                 "Lin", "TCSS",
                                 "jaccard"),
           measure.method.kegg = c("jaccard"),
           measure.method.reactome = c("jaccard"),
           path = "result",
           save_to_local = FALSE) {
    measure.method.go <-
      match.arg(measure.method.go)
    measure.method.kegg <-
      match.arg(measure.method.kegg)
    measure.method.reactome <-
      match.arg(measure.method.reactome)

    if (missing(object)) {
      stop("object is required")
    }

    if (!is(object, "functional_module")) {
      stop("object must be result from enrich_pathway function")
    }

    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    ###GO database
    message(rep("-", 20))
    message("GO database...")
    merged_pathway_go <-
      merge_pathways_internal(
        pathway_result = object@enrichment_go_result,
        p.adjust.cutoff = p.adjust.cutoff.go,
        count.cutoff = count.cutoff.go,
        database = "go",
        sim.cutoff = sim.cutoff.go,
        measure.method = measure.method.go,
        path = path,
        save_to_local = save_to_local
      )
    ###KEGG database
    message(rep("-", 20))
    message("KEGG database...")
    merged_pathway_kegg <-
      merge_pathways_internal(
        pathway_result = object@enrichment_kegg_result,
        p.adjust.cutoff = p.adjust.cutoff.kegg,
        count.cutoff = count.cutoff.kegg,
        database = "kegg",
        sim.cutoff = sim.cutoff.kegg,
        measure.method = measure.method.kegg,
        path = path,
        save_to_local = save_to_local
      )

    ###Reactome database
    message(rep("-", 20))
    message("Reactome database...")
    merged_pathway_reactome <-
      merge_pathways_internal(
        pathway_result = object@enrichment_reactome_result,
        p.adjust.cutoff = p.adjust.cutoff.reactome,
        count.cutoff = count.cutoff.reactome,
        database = "reactome",
        sim.cutoff = sim.cutoff.reactome,
        measure.method = measure.method.reactome,
        path = path,
        save_to_local = save_to_local
      )

    slot(object, "merged_pathway_go") <-
      merged_pathway_go
    slot(object, "merged_pathway_kegg") <-
      merged_pathway_kegg
    slot(object, "merged_pathway_reactome") <-
      merged_pathway_reactome
    slot(object, "merged_module") <-
      list()

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "merge_pathways()",
      parameter = list(
        p.adjust.cutoff.go = p.adjust.cutoff.go,
        p.adjust.cutoff.kegg = p.adjust.cutoff.kegg,
        p.adjust.cutoff.reactome = p.adjust.cutoff.reactome,
        count.cutoff.go = count.cutoff.go,
        count.cutoff.kegg = count.cutoff.kegg,
        count.cutoff.reactome = count.cutoff.reactome,
        sim.cutoff.go = sim.cutoff.go,
        sim.cutoff.kegg = sim.cutoff.kegg,
        sim.cutoff.reactome = sim.cutoff.reactome,
        measure.method.go = measure.method.go,
        measure.method.kegg = measure.method.kegg,
        measure.method.reactome = measure.method.reactome,
        path = path
      ),
      time = Sys.time()
    )

    process_info <-
      slot(object, "process_info")

    process_info$merge_pathways <-
      parameter

    slot(object, "process_info") <-
      process_info

    message("Done")

    object
  }




#' Merge Pathway Enrichment Results Internally
#'
#' This function merges pathway enrichment results obtained through various databases (GO, KEGG, Reactome).
#' It applies similarity measures to find closely related pathways and categorizes them into modules.
#'
#' @param pathway_result A required object containing results from the `enrich_pathway` function.
#' @param p.adjust.cutoff Numeric, p-adjusted value cutoff for filtering enriched pathways.
#' @param count.cutoff Numeric, count cutoff for filtering enriched pathways.
#' @param database Character vector, the database from which the enrichment results were obtained ('go', 'kegg', 'reactome').
#' @param sim.cutoff Numeric, similarity cutoff for clustering pathways.
#' @param measure.method Character, method for calculating term similarity.
#' @param path Character, directory to save intermediate and final results.
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#'
#' @return A list containing `graph_data`, `module_result`, and `result_with_module`.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Load pathway results obtained through `enrich_pathway` function
#' pathway_results <- load_pathway_results("path/to/results")
#'
#' # Merge pathways and find modules
#' merged_results <- merge_pathways_internal(pathway_result = pathway_results)
#'}

merge_pathways_internal <-
  function(pathway_result,
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database = c("go", "kegg", "reactome"),
           sim.cutoff = 0.5,
           measure.method = c("Wang", "Resnik",
                              "Rel", "Jiang",
                              "Lin", "TCSS",
                              "jaccard"),
           path = "result",
           save_to_local = FALSE) {
    measure.method <-
      match.arg(measure.method)
    database <- match.arg(database)
    path <- file.path(path, database)

    if (missing(pathway_result)) {
      stop("pathway_result is required")
    }

    if(is.null(pathway_result)){
      return(list())
    }

    if (!is(pathway_result, "enrichResult")) {
      stop("pathway_result must be result from enrich_pathway function")
    }

    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      dir.create(
        file.path(path, "intermediate_data"),
        showWarnings = FALSE,
        recursive = TRUE
      )
    }

    result <-
      pathway_result@result %>%
      dplyr::filter(p.adjust < p.adjust.cutoff) %>%
      dplyr::filter(Count > count.cutoff) %>%
      dplyr::arrange(p.adjust)

    if (nrow(result) == 0) {
      return(NULL)
    }

    ##get the similartiy matrix
    message("Calculating similartiy matrix, it may take a while...")
    if (database == "go") {
      sim_matrix <-
        get_go_result_sim(
          result = result,
          sim.cutoff = sim.cutoff,
          measure.method = measure.method
        )
    }

    if (database == "kegg") {
      sim_matrix <-
        tryCatch(
          sim_matrix <-
            simplifyEnrichment::term_similarity_from_KEGG(term_id = c(result$ID),
                                                          method = "jaccard") %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "name1") %>%
            tidyr::pivot_longer(
              cols = -name1,
              names_to = "name2",
              values_to = "sim"
            ) %>%
            dplyr::filter(name1 != name2),
          error = function(x) {
            data.frame(name1 = character(),
                       name2 = character(),
                       sim = numeric())
          }
        )
    }

    if (database == "reactome") {
      sim_matrix <-
        tryCatch(
          sim_matrix <-
            simplifyEnrichment::term_similarity_from_Reactome(term_id = c(result$ID),
                                                              method = "jaccard") %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "name1") %>%
            tidyr::pivot_longer(
              cols = -name1,
              names_to = "name2",
              values_to = "sim"
            ) %>%
            dplyr::filter(name1 != name2),
          error = function(x) {
            data.frame(name1 = character(),
                       name2 = character(),
                       sim = numeric())
          }
        )
    }

    if (save_to_local) {
      save(sim_matrix,
           file = file.path(path, "intermediate_data/sim_matrix"))
    }

    ####module detection
    message("Identifying modules...")

    edge_data <-
      rbind(sim_matrix) %>%
      dplyr::rename(from = name1, to = name2) %>%
      dplyr::filter(sim > sim.cutoff)

    node_data <-
      rbind(result) %>%
      as.data.frame() %>%
      dplyr::select(ID, everything()) %>%
      dplyr::rename(node = ID)

    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())

    subnetwork <-
      suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data,
                                                        weights = abs(edge_attr(graph_data,
                                                                                "sim"))))

    # save(subnetwork, file = file.path(path, "subnetwork"))
    cluster <-
      paste(database,
            "Module",
            as.character(igraph::membership(subnetwork)),
            sep = "_")

    graph_data <-
      graph_data %>%
      igraph::upgrade_graph() %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::mutate(module = cluster)

    ###clustered different GO terms
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(module, p.adjust)

    if (database == "go") {
      result_with_module <-
        result_with_module %>%
        dplyr::arrange(ONTOLOGY, module, p.adjust)
    }

    ###add module content number
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

    if (save_to_local) {
      save(result_with_module,
           file = file.path(path, "intermediate_data/result_with_module"))
      save(graph_data,
           file = file.path(path, "intermediate_data/graph_data"))
    }

    module_result <-
      result_with_module %>%
      plyr::dlply(.variables = .(module)) %>%
      purrr::map(function(x) {
        # cat(unique(x$module), " ")
        if (nrow(x) == 1) {
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
        x$geneID =
          x$geneID %>%
          stringr::str_split(pattern = "/") %>%
          unlist() %>%
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
      dplyr::mutate(
        module_annotation = case_when(
          module == "Other" ~ Description,
          module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
        )
      ) %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(module_annotation, everything())

    module_result$module_annotation <-
      stringr::str_split(module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()

    module_result$Count <-
      as.numeric(module_result$Count)

    module_result$p.adjust <-
      as.numeric(module_result$p.adjust)

    module_result$qvalue <-
      as.numeric(module_result$qvalue)

    module_result$degree <-
      as.numeric(module_result$degree)

    if (save_to_local) {
      save(module_result,
           file = file.path(path, "intermediate_data/module_result"))
    }

    message("Done")

    list(
      graph_data = graph_data,
      module_result = module_result,
      result_with_module = result_with_module
    )
  }
