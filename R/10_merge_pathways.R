# ##pathway enrichment analysis
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")
# library(ggraph)
# library(igraph)
# library(tidygraph)
# library(tidyverse)
# library(extrafont)
# library(simplifyEnrichment)
# # library(GOSim)
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
#     measure.method.go = "Sim_Wang_2007",
#     measure.method.kegg = "jaccard",
#     measure.method.reactome = "jaccard",
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_modules, file = "result/enriched_modules")
#
#
#

# ##GSEA analysis results
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/covid_data/")
# library(ggraph)
# library(igraph)
# library(tidygraph)
# library(tidyverse)
# library(extrafont)
# library(simplifyEnrichment)
# # library(GOSim)
# library(plyr)
#
# load("result/gsea_pathways")
#
# object = gsea_pathways
# p.adjust.cutoff.go = 0.05
# p.adjust.cutoff.kegg = 0.05
# p.adjust.cutoff.reactome = 0.05
# count.cutoff.go = 5
# count.cutoff.kegg = 5
# count.cutoff.reactome = 5
# sim.cutoff.go = 0.5
# sim.cutoff.kegg = 0.5
# sim.cutoff.reactome = 0.5
# measure.method.go = "Sim_Wang_2007"
# measure.method.kegg = "jaccard"
# measure.method.reactome = "jaccard"
# path = "result"
# save_to_local = FALSE
#
# gsea_enriched_modules <-
#   merge_pathways(
#     object = gsea_pathways,
#     p.adjust.cutoff.go = 0.05,
#     p.adjust.cutoff.kegg = 0.05,
#     p.adjust.cutoff.reactome = 0.05,
#     count.cutoff.go = 5,
#     count.cutoff.kegg = 5,
#     count.cutoff.reactome = 5,
#     sim.cutoff.go = 0.5,
#     sim.cutoff.kegg = 0.5,
#     sim.cutoff.reactome = 0.5,
#     measure.method.go = "Sim_Wang_2007",
#     measure.method.kegg = "jaccard",
#     measure.method.reactome = "jaccard",
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_modules, file = "result/enriched_modules")
# #


####GSEA analysis

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
#' @param measure.method.go A character vector specifying the similarity measure method for GO. Choices are "Sim_Wang_2007", "Sim_Lin_1998", "Sim_Resnik_1999", "Sim_FaITH_2010", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_XGraSM_2013", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Zhang_2006", "Sim_universal", "Sim_GOGO_2018", "Sim_Rada_1989", "Sim_Resnik_edge_2005", "Sim_Leocock_1998", "Sim_WP_1994", "Sim_Slimani_2006", "Sim_Shenoy_2012", "Sim_Pekar_2002", "Sim_Stojanovic_2001", "Sim_Wang_edge_2012", "Sim_Zhong_2002", "Sim_AlMubaid_2006", "Sim_Li_2003", "Sim_RSS_2013", "Sim_HRSS_2013", "Sim_Shen_2010", "Sim_SSDD_2013", "Sim_Jiang_1997", "Sim_Kappa", "Sim_Jaccard", "Sim_Dice",  "Sim_Overlap", "Sim_Ancestor".
#' @param measure.method.kegg A character vector specifying the similarity measure method for KEGG. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#' @param measure.method.reactome A character vector specifying the similarity measure method for Reactome. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
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
           measure.method.go = c("Sim_Wang_2007", "Sim_Lin_1998", "Sim_Resnik_1999", "Sim_FaITH_2010", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_XGraSM_2013", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Zhang_2006", "Sim_universal", "Sim_GOGO_2018", "Sim_Rada_1989", "Sim_Resnik_edge_2005", "Sim_Leocock_1998", "Sim_WP_1994", "Sim_Slimani_2006", "Sim_Shenoy_2012", "Sim_Pekar_2002", "Sim_Stojanovic_2001", "Sim_Wang_edge_2012", "Sim_Zhong_2002", "Sim_AlMubaid_2006", "Sim_Li_2003", "Sim_RSS_2013", "Sim_HRSS_2013", "Sim_Shen_2010", "Sim_SSDD_2013", "Sim_Jiang_1997", "Sim_Kappa", "Sim_Jaccard", "Sim_Dice",  "Sim_Overlap", "Sim_Ancestor"),
           measure.method.kegg = c("jaccard", "dice", "overlap", "kappa"),
           measure.method.reactome = c("jaccard", "dice", "overlap", "kappa"),
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

    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
    } else{
      analysis_type <- "do_gsea"
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
        analysis_type = analysis_type,
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
        analysis_type = analysis_type,
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
        analysis_type = analysis_type,
        p.adjust.cutoff = p.adjust.cutoff.reactome,
        count.cutoff = count.cutoff.reactome,
        database = "reactome",
        sim.cutoff = sim.cutoff.reactome,
        measure.method = measure.method.reactome,
        path = path,
        save_to_local = save_to_local
      )

    if (is.null(merged_pathway_go)) {
      merged_pathway_go <- list()
    }

    if (is.null(merged_pathway_kegg)) {
      merged_pathway_kegg <- list()
    }

    if (is.null(merged_pathway_reactome)) {
      merged_pathway_reactome <- list()
    }

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
#' @param analysis_type Character, type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`.
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
           analysis_type = c("enrich_pathway", "do_gsea"),
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database = c("go", "kegg", "reactome"),
           sim.cutoff = 0.5,
           measure.method = c("Sim_Wang_2007", "Sim_Lin_1998", "Sim_Resnik_1999", "Sim_FaITH_2010", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_XGraSM_2013", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Zhang_2006", "Sim_universal", "Sim_GOGO_2018", "Sim_Rada_1989", "Sim_Resnik_edge_2005", "Sim_Leocock_1998", "Sim_WP_1994", "Sim_Slimani_2006", "Sim_Shenoy_2012", "Sim_Pekar_2002", "Sim_Stojanovic_2001", "Sim_Wang_edge_2012", "Sim_Zhong_2002", "Sim_AlMubaid_2006", "Sim_Li_2003", "Sim_RSS_2013", "Sim_HRSS_2013", "Sim_Shen_2010", "Sim_SSDD_2013", "Sim_Jiang_1997", "Sim_Kappa", "Sim_Jaccard", "Sim_Dice",  "Sim_Overlap", "Sim_Ancestor", "jaccard", "dice", "overlap", "kappa"),
           path = "result",
           save_to_local = FALSE) {

    analysis_type <- match.arg(analysis_type)
    measure.method <-
      match.arg(measure.method)
    database <- match.arg(database)
    path <- file.path(path, database)

    if (missing(pathway_result)) {
      stop("pathway_result is required")
    }

    if (is.null(pathway_result)) {
      return(list())
    }

    if (!is(pathway_result, "enrichResult") &
        !is(pathway_result, "gseaResult")) {
      stop("pathway_result must be result from enrich_pathway or do_gsea function")
    }

    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      dir.create(
        file.path(path, "intermediate_data"),
        showWarnings = FALSE,
        recursive = TRUE
      )
    }

    if (is(pathway_result, "enrichResult")) {
      result <-
        pathway_result@result %>%
        dplyr::filter(p.adjust < p.adjust.cutoff) %>%
        dplyr::filter(Count > count.cutoff) %>%
        dplyr::arrange(p.adjust)
    } else{
      result <-
        pathway_result@result %>%
        dplyr::filter(p.adjust < p.adjust.cutoff) %>%
        dplyr::arrange(p.adjust)
    }


    if (database == "go") {
      result <-
        dplyr::filter(result, ONTOLOGY != "CC")
    }

    if (nrow(result) == 0) {
      return(NULL)
    }

    ##get the similartiy matrix
    message("Calculating similartiy matrix, it may take a while...")
    if (database == "go") {
      sim_matrix <-
        tryCatch(
          get_go_result_sim(
            result = dplyr::filter(result, ONTOLOGY != "CC"),
            sim.cutoff = sim.cutoff,
            measure.method = measure.method
          ),
          error = function(x) {
            data.frame(name1 = character(),
                       name2 = character(),
                       sim = numeric())
          }
        )
    }

    if (database == "kegg") {
      sim_matrix <-
        tryCatch(
          sim_matrix <-
            term_similarity_KEGG(term_id = c(result$ID),
                                 measure.method = measure.method) %>%
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
            term_similarity_Reactome(term_id = c(result$ID),
                                     measure.method = measure.method) %>%
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
      save(sim_matrix, file = file.path(path, "intermediate_data/sim_matrix"))
    }

    ####module detection
    message("Identifying modules...")

    identify_modules(
      sim_matrix = sim_matrix,
      analysis_type = analysis_type,
      result = result,
      database = database,
      sim.cutoff = sim.cutoff,
      save_to_local = save_to_local,
      path = path
    )
  }




#' Identify Modules in Similarity Matrix
#'
#' This function identifies modules (clusters) in a given similarity matrix based on pathway enrichment or gene set enrichment analysis (GSEA). It constructs a network graph using similarity and result data, applies clustering algorithms, and optionally saves the results to a specified path.
#'
#' @param sim_matrix A data frame containing the similarity matrix with columns `name1`, `name2`, and `sim` representing the similarity between entities.
#' @param analysis_type Character. Type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`. Default is `"enrich_pathway"`.
#' @param result A data frame containing the enrichment analysis results, including columns like `ID`, `Description`, `p.adjust`, and other relevant data.
#' @param database Character. The database from which the enrichment results were obtained (`go`, `kegg`, `reactome`).
#' @param sim.cutoff Numeric. The similarity cutoff value used to filter the edges in the similarity matrix. Default is `0.5`.
#' @param save_to_local Logical. Whether to save the resulting data to local files. Default is `TRUE`.
#' @param path Character. The directory path where intermediate results will be saved, if `save_to_local = TRUE`. Default is an empty string (current working directory).
#'
#' @return A list containing:
#' \item{graph_data}{A tidygraph object representing the network with nodes and edges.}
#' \item{module_result}{A data frame with the identified modules and their associated information, including pathway descriptions and p-values.}
#' \item{result_with_module}{A data frame with the original result data enriched with module information.}
#'
#' @details
#' The function first constructs a graph from the similarity matrix and result data, then applies clustering to identify modules. For pathway enrichment, the function organizes and processes the result data for easier interpretation. If `save_to_local` is `TRUE`, it saves intermediate results in the specified `path`.
#'
#' @examples
#' \dontrun{
#' sim_matrix <- data.frame(
#'   name1 = c("A", "B", "C"),
#'   name2 = c("B", "C", "A"),
#'   sim = c(0.6, 0.7, 0.8)
#' )
#' result <- data.frame(
#'   ID = c("P1", "P2", "P3"),
#'   p.adjust = c(0.01, 0.05, 0.03),
#'   Description = c("Pathway 1", "Pathway 2", "Pathway 3")
#' )
#' modules <- identify_modules(
#'   sim_matrix,
#'   "enrich_pathway",
#'   result,
#'   sim.cutoff = 0.5,
#'   save_to_local = FALSE
#' )
#' }
#'
#' @importFrom dplyr rename filter select mutate count arrange everything left_join case_when distinct
#' @importFrom tidygraph tbl_graph centrality_degree activate
#' @importFrom igraph cluster_edge_betweenness edge_attr membership upgrade_graph vertex_attr
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @importFrom plyr dlply
#' @export


identify_modules <-
  function(sim_matrix,
           analysis_type = c("enrich_pathway", "do_gsea"),
           result,
           database = c("go", "kegg", "reactome"),
           sim.cutoff = 0.5,
           save_to_local = TRUE,
           path = "") {
    analysis_type <- match.arg(analysis_type)
    database <- match.arg(database)

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
      suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))

    # save(subnetwork, file = file.path(path, "subnetwork"))
    cluster <-
      paste(database, "Module", as.character(igraph::membership(subnetwork)), sep = "_")

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
      save(graph_data, file = file.path(path, "intermediate_data/graph_data"))
    }

    if (analysis_type == "enrich_pathway") {
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
        dplyr::mutate(module_annotation = ifelse(module == "Other", Description, sapply(strsplit(
          Description, ";"
        ), `[`, 1))) %>%
        dplyr::mutate(
          p.adjust = as.numeric(p.adjust),
          Count = as.numeric(Count),
          p.adjust = as.numeric(p.adjust),
          qvalue = as.numeric(qvalue),
          degree = as.numeric(degree)
        ) %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::select(module_annotation, everything())
    } else{
      module_result <-
        result_with_module %>%
        plyr::dlply(.variables = .(module)) %>%
        purrr::map(function(x) {
          # cat(unique(x$module), " ")
          if (nrow(x) == 1) {
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

      module_result$Count <-
        purrr::map(module_result$core_enrichment, function(x) {
          length(stringr::str_split(x, pattern = "/")[[1]])
        }) %>%
        unlist()
    }

    module_result$module_annotation <-
      stringr::str_split(module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()

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
