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
# library(massdataset)
#
# load("result/enriched_pathways")
#
# merged_pathways <- merge_pathways(
#   object = enrich_pathway_res,
#   database = c("go", "kegg"),
#   sim.cutoff.go = 0.5,
#   sim.cutoff.kegg = 0.5,
#   measure.method.go = "Sim_Lin_1998",
#   go.orgdb = mf.orgdb,
#   measure.method.kegg = "jaccard"
# )
#
# enriched_modules <-
#   merge_pathways(
#    object = enriched_pathways,
#    p.adjust.cutoff.go = 0.05,
#    go.orgdb = org.Hs.eg.db,
#    p.adjust.cutoff.kegg = 0.05,
#    p.adjust.cutoff.reactome = 0.05,
#    count.cutoff.go = 5,
#    count.cutoff.kegg = 5,
#    count.cutoff.reactome = 5,
#    sim.cutoff.go = 0.5,
#    sim.cutoff.kegg = 0.5,
#    sim.cutoff.reactome = 0.5,
#    measure.method.go = "Sim_Wang_2007",
#    control.method.go = list(contribution_factor = c("is_a" = 0.8,
#                                           "part_of" = 0.6,
#                                           "regulates" = 0.7,
#                                           "negatively_regulates" = 0.75,
#                                           "positively_regulates" = 0.75)),
#    measure.method.kegg = "jaccard",
#    measure.method.reactome = "jaccard",
#    path = "result",
#    save_to_local = FALSE
#   )
#
# save(enriched_modules, file = "result/enriched_modules")

# ##GSEA analysis results
# load("result/gsea_pathways")
#
# object = gsea_pathways
# database = c("go", "kegg", "reactome")
# p.adjust.cutoff.go = 0.05
# p.adjust.cutoff.kegg = 0.05
# p.adjust.cutoff.reactome = 0.05
# count.cutoff.go = 5
# count.cutoff.kegg = 5
# count.cutoff.reactome = 5
# sim.cutoff.go = 0.5
# sim.cutoff.kegg = 0.5
# sim.cutoff.reactome = 0.5
# measure.method.go = "Sim_XGraSM_2013"
# measure.method.kegg = "jaccard"
# measure.method.reactome = "jaccard"
# path = "result"
# save_to_local = FALSE

# gsea_enriched_modules <-
#   merge_pathways(
#     object = gsea_pathways,
#     database = c("go", "kegg", "reactome"),
#     p.adjust.cutoff.go = 0.05,
#     p.adjust.cutoff.kegg = 0.05,
#     p.adjust.cutoff.reactome = 0.05,
#     count.cutoff.go = 5,
#     count.cutoff.kegg = 5,
#     count.cutoff.reactome = 5,
#     sim.cutoff.go = 0.5,
#     sim.cutoff.kegg = 0.5,
#     sim.cutoff.reactome = 0.5,
#     measure.method.go = "Sim_XGraSM_2013",
#     measure.method.kegg = "jaccard",
#     measure.method.reactome = "jaccard",
#     path = "result",
#     save_to_local = FALSE
#   )
#
# save(enriched_modules, file = "result/enriched_modules")

# For metabolite
# object <- enriched_pathways
# merged_pathways <-
#   merge_pathways(
#     object = enriched_pathways,
#     database = c("hmdb", "kegg"),
#     sim.cutoff.hmdb = 0,
#     sim.cutoff.metkegg = 0
#   )

#' This function merges enrichment analysis results from different databases (GO, KEGG, and Reactome) into a single object respectively. The function takes an object of class "functional_module" and set up various parameters to compute and filter term similarity scores.
#'
#' @param object An object of class "functional_module", typically a result from enrich_pathway function.
#' @param database Character vector, specify which database(s) to use for merging pathways.
#'   - For genes: 'go', 'kegg', 'reactome'
#'   - For metabolites: 'hmdb', 'metkegg'
#' @param p.adjust.cutoff.go Adjusted p-value cutoff for GO database. Default is 0.05.
#' @param p.adjust.cutoff.kegg Adjusted p-value cutoff for KEGG database. Default is 0.05.
#' @param p.adjust.cutoff.reactome Adjusted p-value cutoff for Reactome database. Default is 0.05.
#' @param p.adjust.cutoff.hmdb Adjusted p-value cutoff for HMDB pathways from metabolite enrichment analysis result. Default is 0.05.
#' @param p.adjust.cutoff.metkegg Adjusted p-value cutoff for KEGG pathways from metabolite enrichment analysis result. Default is 0.05.
#' @param count.cutoff.go Gene count cutoff for GO terms from gene enrichment analysis. Default is 5.
#' @param count.cutoff.kegg Gene count cutoff for KEGG pathways from gene enrichment analysis. Default is 5.
#' @param count.cutoff.reactome Count cutoff for Reactome pathways from gene enrichment analysis. Default is 5.
#' @param count.cutoff.hmdb Metabolite count cutoff for HMDB pathways from metabolite enrichment analysis result. Default is 5.
#' @param count.cutoff.metkegg Metabolite count cutoff for KEGG pathways from metabolite enrichment analysis result. Default is 5.
#' @param sim.cutoff.go Similarity cutoff for GO database. Default is 0.5.
#' @param sim.cutoff.kegg Similarity cutoff for KEGG database. Default is 0.5.
#' @param sim.cutoff.reactome Similarity cutoff for Reactome database. Default is 0.5.
#' @param sim.cutoff.hmdb Similarity cutoff for HMDB database when interpreting metabolite enrichment result. Default is 0.5.
#' @param sim.cutoff.metkegg Similarity cutoff for KEGG database when interpreting metabolite enrichment result. Default is 0.5.
#' @param measure.method.go A character vector specifying the term semantic similarity measure method for GO terms. Default is `"Sim_XGraSM_2013"`. See `simona::term_sim()` for available measures.
#' @param control.method.go a list of parameters passing to specified measure method for GO term semantic similarity. For details about how to set this parameter, please go to https://jokergoo.github.io/simona/articles/v05_term_similarity.html.
#' @param go.orgdb OrgDb object or character string naming the *OrgDb* annotation package used to derive gene–GO mappings, required when database includes "go".
#' @param measure.method.kegg A character vector specifying the similarity measure method for KEGG. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#' @param measure.method.reactome A character vector specifying the similarity measure method for Reactome. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#' @param measure.method.hmdb A character vector specifying the similarity measure method for HMDB when interpreting metabolite enrichment result. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#' @param measure.method.metkegg A character vector specifying the similarity measure method for KEGG when interpreting metabolite enrichment result. Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#' @param path Directory path to save the results. Default is "result".
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#'
#' @return An object of class "functional_module" with slots for merged pathways from each database.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export
#'

merge_pathways <-
  function(object,
           database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
           p.adjust.cutoff.go = 0.05,
           p.adjust.cutoff.kegg = 0.05,
           p.adjust.cutoff.reactome = 0.05,
           p.adjust.cutoff.hmdb = 0.05,
           p.adjust.cutoff.metkegg = 0.05,
           count.cutoff.go = 5,
           count.cutoff.kegg = 5,
           count.cutoff.reactome = 5,
           count.cutoff.hmdb = 5,
           count.cutoff.metkegg = 5,
           sim.cutoff.go = 0.5,
           sim.cutoff.kegg = 0.5,
           sim.cutoff.reactome = 0.5,
           sim.cutoff.hmdb = 0.5,
           sim.cutoff.metkegg = 0.5,
           measure.method.go = c("Sim_XGraSM_2013", "Sim_Wang_2007", "Sim_Lin_1998", "Sim_Resnik_1999", "Sim_FaITH_2010", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Zhang_2006", "Sim_universal", "Sim_GOGO_2018", "Sim_Rada_1989", "Sim_Resnik_edge_2005", "Sim_Leocock_1998", "Sim_WP_1994", "Sim_Slimani_2006", "Sim_Shenoy_2012", "Sim_Pekar_2002", "Sim_Stojanovic_2001", "Sim_Wang_edge_2012", "Sim_Zhong_2002", "Sim_AlMubaid_2006", "Sim_Li_2003", "Sim_RSS_2013", "Sim_HRSS_2013", "Sim_Shen_2010", "Sim_SSDD_2013", "Sim_Jiang_1997", "Sim_Kappa", "Sim_Jaccard", "Sim_Dice",  "Sim_Overlap", "Sim_Ancestor"),
           control.method.go = list(),
           go.orgdb = NULL,
           measure.method.kegg = c("jaccard", "dice", "overlap", "kappa"),
           measure.method.reactome = c("jaccard", "dice", "overlap", "kappa"),
           measure.method.hmdb = c("jaccard", "dice", "overlap", "kappa"),
           measure.method.metkegg = c("jaccard", "dice", "overlap", "kappa"),
           path = "result",
           save_to_local = FALSE) {

    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <- object@process_info$enrich_pathway@parameter$query_type
    } else{
      analysis_type <- "do_gsea"
      query_type <- object@process_info$do_gsea@parameter$query_type
    }

    if (query_type == "gene" && missing(database)) {
      stop("Please specify databases to merge pathways.")
      database <- match.arg(database, choices = c("go", "kegg", "reactome"), several.ok = TRUE)
    } else if (query_type == "metabolite") {
      database <- c("hmdb", "metkegg")
    }

    ## Check input parameter for different query type
    if (query_type == "gene") {
      if ("go" %in% database) {
        measure.method.go <-
          match.arg(measure.method.go)
      }

      if ("kegg" %in% database) {
        measure.method.kegg <-
          match.arg(measure.method.kegg)
      }

      if ("reactome" %in% database) {
        measure.method.reactome <-
          match.arg(measure.method.reactome)
      }
    } else if (query_type == "metabolite") {
      if ("hmdb" %in% database) {
        measure.method.hmdb <-
          match.arg(measure.method.hmdb)
      }

      if ("metkegg" %in% database) {
        measure.method.metkegg <-
          match.arg(measure.method.metkegg)
      }
    }

    if (missing(object)) {
      stop("object is required")
    }

    if (!is(object, "functional_module")) {
      stop("object must be result from enrich_pathway function")
    }

    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }

    # Merge pathway
    ###GO database
    if (query_type == "gene") {
      merged_pathway_hmdb <- NULL
      merged_pathway_metkegg <- NULL

      if ("go" %in% database) {
        message(rep("-", 20))
        message("GO database...")

        merged_pathway_go <-
          merge_pathways_internal(
            query_type = query_type,
            pathway_result = object@enrichment_go_result,
            analysis_type = analysis_type,
            go.orgdb = go.orgdb,
            p.adjust.cutoff = p.adjust.cutoff.go,
            count.cutoff = count.cutoff.go,
            database = "go",
            sim.cutoff = sim.cutoff.go,
            measure.method = measure.method.go,
            control.method.go = control.method.go,
            path = path,
            save_to_local = save_to_local
          )
      } else {
        merged_pathway_go <- NULL
      }

      if ("kegg" %in% database) {
        message(rep("-", 20))
        message("KEGG database...")

        merged_pathway_kegg <-
          merge_pathways_internal(
            query_type = query_type,
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
      } else {
        merged_pathway_kegg <- NULL
      }

      if ("reactome" %in% database) {
        message(rep("-", 20))
        message("Reactome database...")

        merged_pathway_reactome <-
          merge_pathways_internal(
            query_type = query_type,
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
      } else {
        merged_pathway_reactome <- NULL
      }
    } else if (query_type == "metabolite") {
      merged_pathway_go <- NULL
      merged_pathway_kegg <- NULL
      merged_pathway_reactome <- NULL

      if ("hmdb" %in% database) {
        message(rep("-", 20))
        message("HMDB database...")

        merged_pathway_hmdb <-
          merge_pathways_internal(
            query_type = query_type,
            pathway_result = object@enrichment_hmdb_result,
            analysis_type = analysis_type,
            p.adjust.cutoff = p.adjust.cutoff.hmdb,
            count.cutoff = count.cutoff.hmdb,
            database = "hmdb",
            sim.cutoff = sim.cutoff.hmdb,
            measure.method = measure.method.hmdb,
            path = path,
            save_to_local = save_to_local
          )
      } else {
        merged_pathway_hmdb <- NULL
      }

      if ("metkegg" %in% database) {
        message(rep("-", 20))
        message("KEGG database...")

        merged_pathway_metkegg <-
          merge_pathways_internal(
            query_type = query_type,
            pathway_result = object@enrichment_metkegg_result,
            analysis_type = analysis_type,
            p.adjust.cutoff = p.adjust.cutoff.metkegg,
            count.cutoff = count.cutoff.metkegg,
            database = "metkegg",
            sim.cutoff = sim.cutoff.metkegg,
            measure.method = measure.method.metkegg,
            path = path,
            save_to_local = save_to_local
          )
      } else {
        merged_pathway_metkegg <- NULL
      }
    }

    if (is.null(merged_pathway_go)) {
      merged_pathway_go <- list()
    }

    if (is.null(merged_pathway_kegg)) {
      merged_pathway_kegg <- list()
    }

    if (is.null(merged_pathway_reactome)) {
      merged_pathway_reactome <- list()
    }

    if (is.null(merged_pathway_hmdb)) {
      merged_pathway_hmdb <- list()
    }

    if (is.null(merged_pathway_metkegg)) {
      merged_pathway_metkegg <- list()
    }

    slot(object, "merged_pathway_go") <-
      merged_pathway_go
    slot(object, "merged_pathway_kegg") <-
      merged_pathway_kegg
    slot(object, "merged_pathway_reactome") <-
      merged_pathway_reactome
    slot(object, "merged_pathway_hmdb") <-
      merged_pathway_hmdb
    slot(object, "merged_pathway_metkegg") <-
      merged_pathway_metkegg
    slot(object, "merged_module") <-
      list()

    if (query_type == "gene") {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "merge_pathways()",
        parameter = list(
          query_type = query_type,
          database = database,
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
          measure.method.reactome = measure.method.reactome
        ),
        time = Sys.time()
      )
    } else if (query_type == "metabolite") {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "merge_pathways()",
        parameter = list(
          query_type = query_type,
          database = database,
          p.adjust.cutoff.hmdb = p.adjust.cutoff.hmdb,
          p.adjust.cutoff.metkegg = p.adjust.cutoff.metkegg,
          count.cutoff.hmdb = count.cutoff.hmdb,
          count.cutoff.metkegg = count.cutoff.metkegg,
          sim.cutoff.hmdb = sim.cutoff.hmdb,
          sim.cutoff.metkegg = sim.cutoff.metkegg,
          measure.method.hmdb = measure.method.hmdb,
          measure.method.metkegg = measure.method.metkegg
        ),
        time = Sys.time()
      )
    }

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
#' @description
#' Internal function that merges pathway enrichment results obtained through various databases
#' (GO, KEGG, Reactome, HMDB). It applies similarity measures to find closely related pathways
#' and categorizes them into modules.
#'
#' @param query_type Character, the category of biological entity to query: either "gene" or "metabolite".
#' @param pathway_result A required object containing results from the `enrich_pathway` or `do_gsea` function.
#' @param analysis_type Character, type of analysis performed: either "enrich_pathway" or "do_gsea".
#' @param go.orgdb An optional organism-specific database for GO annotations, required when database="go".
#' @param p.adjust.cutoff Numeric, p-adjusted value cutoff for filtering enriched pathways (default: 0.05).
#' @param count.cutoff Numeric, count cutoff for filtering enriched pathways (default: 5).
#' @param database Character, the database from which the enrichment results were obtained:
#'                one of 'go', 'kegg', 'reactome', or 'hmdb'.
#' @param sim.cutoff Numeric, similarity cutoff for clustering pathways (default: 0.5).
#' @param measure.method Character, method for calculating term similarity.
#'                       See `simona::all_term_sim_methods()` for available measures.
#' @param control.method.go A list of parameters passing to the specified measure method for GO term semantic
#'                         similarity. For details, see https://jokergoo.github.io/simona/articles/v05_term_similarity.html.
#' @param path Character, directory to save intermediate and final results (default: "result").
#' @param save_to_local Logical, if TRUE the results will be saved to local disk (default: FALSE).
#'
#' @return A list containing `graph_data`, `module_result`, and `result_with_module`, or NULL if no
#'         pathways meet the filtering criteria.
#'
#' @importFrom dplyr filter arrange rename
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @keywords internal
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}

merge_pathways_internal <-
  function(query_type = c("gene", "metabolite"),
           pathway_result,
           analysis_type = c("enrich_pathway", "do_gsea"),
           go.orgdb = NULL,
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
           sim.cutoff = 0.5,
           measure.method,
           control.method.go = list(),
           path = "result",
           save_to_local = FALSE) {

    query_type <- match.arg(query_type)
    analysis_type <- match.arg(analysis_type)
    database <- match.arg(database)
    if (database == "go" && (missing(go.orgdb) || is.null(go.orgdb))) {
      stop("Please provide the OrgDb object to merge GO terms.")
    }

    if (save_to_local) {
      path <- file.path(path, database)
    }

    if (missing(pathway_result)) {
      stop("pathway_result is required")
    }

    if (is.null(pathway_result)) {
      return(list())
    }

    if (!is(pathway_result, "enrichResult") &
        !is(pathway_result, "gseaResult") &
        !is(pathway_result, "enrich_result")) {
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

    ## Collect enrichment result ====
    if (query_type == "gene") {
      if (is(pathway_result, "enrichResult")) {
        result <-
          pathway_result@result %>%
          dplyr::filter(p_adjust < p.adjust.cutoff) %>%
          dplyr::filter(Count > count.cutoff) %>%
          dplyr::arrange(p_adjust)
      } else{
        result <-
          pathway_result@result %>%
          dplyr::filter(p_adjust < p.adjust.cutoff) %>%
          dplyr::arrange(p_adjust)
        ## TO DO: Whether to add filter for the count of core_enrichment based on the parameter count.cutoff?
        # result$Count <- purrr::map_int(pathway_result@result$core_enrichment,
        #                                function(x) {
        #                                  length(stringr::str_split(x, pattern = "/")[[1]])
        #                                  })
      }
    } else if (query_type == "metabolite") {
      result <-
        pathway_result@result %>%
        dplyr::filter(p_adjust < p.adjust.cutoff) %>%
        dplyr::filter(mapped_number > count.cutoff) %>%
        dplyr::arrange(p_adjust)
    }


    # if (database == "go") {
    #   result <-
    #     dplyr::filter(result, ONTOLOGY != "CC")
    # }

    if (nrow(result) == 0) {
      return(NULL)
    }

    ## Get the similarity matrix ====
    message("Calculating similartiy matrix, it may take a while...")
    if (query_type == "gene") {
      if (database == "go") {
        sim_matrix <-
          tryCatch(
            expr = {
              result_sim <- get_go_result_sim(
                # result = dplyr::filter(result, ONTOLOGY != "CC"),
                result = result,
                go.orgdb = go.orgdb,
                sim.cutoff = sim.cutoff,
                measure.method = measure.method,
                control.method = control.method.go
              )
              message("Completed GO term similarity calculation successfully!")
              result_sim
            },
            error = function(e) {
              message(paste("Error in GO term similarity calculation:", e$message))
              message("Setting GO term sim_matrix as empty data frame and continuing...")
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }

      if (database == "kegg") {
        sim_matrix <-
          tryCatch(
            {
              result_sim <-
                term_similarity_KEGG(term_id = c(result$ID),
                                     measure.method = measure.method) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "name1") %>%
                tidyr::pivot_longer(
                  cols = -name1,
                  names_to = "name2",
                  values_to = "sim"
                ) %>%
                dplyr::filter(name1 < name2)
              message("Completed KEGG pathway similarity calculation successfully!")
              result_sim
            },
            error = function(e) {
              message(paste("Error in KEGG pathway similarity calculation:", e$message))
              message("Setting KEGG pathway sim_matrix as empty data frame and continuing...")
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }

      if (database == "reactome") {
        sim_matrix <-
          tryCatch(
            {
              result_sim <-
                term_similarity_Reactome(term_id = c(result$ID),
                                         measure.method = measure.method) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "name1") %>%
                tidyr::pivot_longer(
                  cols = -name1,
                  names_to = "name2",
                  values_to = "sim"
                ) %>%
                dplyr::filter(name1 < name2)
              message("Completed Reactome pathway similarity calculation successfully!")
              result_sim
            },
            error = function(e) {
              message(paste("Error in Reactome pathway similarity calculation:", e$message))
              message("Setting Reactome pathway sim_matrix as empty data frame and continuing...")
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }
    } else if (query_type == "metabolite") {
      if (database == "hmdb") {
        sim_matrix <-
          tryCatch(
            {
              result_sim <-
                term_similarity_metabolite(
                  enrichment_result = result,
                  measure.method = measure.method) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "name1") %>%
                tidyr::pivot_longer(
                  cols = -name1,
                  names_to = "name2",
                  values_to = "sim"
                ) %>%
                dplyr::filter(name1 != name2)

              message("Completed SMPDB pathway similarity calculation successfully!")
              result_sim
            },
            error = function(e) {
              message(paste("Error in SMPDB pathway similarity calculation:", e$message))
              message("Setting SMPDB pathway sim_matrix as empty data frame and continuing...")
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }

      if (database == "metkegg") {
        sim_matrix <-
          tryCatch(
            {
              result_sim <-
                term_similarity_metabolite(
                  enrichment_result = result,
                  measure.method = measure.method) %>%
                as.data.frame() %>%
                tibble::rownames_to_column(var = "name1") %>%
                tidyr::pivot_longer(
                  cols = -name1,
                  names_to = "name2",
                  values_to = "sim"
                ) %>%
                dplyr::filter(name1 != name2)

              message("Completed KEGG pathway similarity calculation successfully!")
              result_sim
            },
            error = function(e) {
              message(paste("Error in KEGG pathway similarity calculation:", e$message))
              message("Setting KEGG pathway sim_matrix as empty data frame and continuing...")
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }
    }

    if (save_to_local) {
      save(sim_matrix, file = file.path(path, "intermediate_data/sim_matrix.RData"))
    }

    ####module detection
    message("Identifying modules...")

    identify_modules(
      sim_matrix = sim_matrix,
      query_type = query_type,
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
#' @param query_type Character, the category of biological entity to query ("gene", "metabolite") for merging pathway enrichment result.
#' @param analysis_type Character. Type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`. Default is `"enrich_pathway"`.
#' @param result A data frame containing the enrichment analysis results, including columns like `ID`, `Description`, `p_adjust`, and other relevant data.
#' @param database Character. The database from which the enrichment results were obtained (`go`, `kegg`, `reactome`, `hmdb`, `metkegg`).
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
#'   p_adjust = c(0.01, 0.05, 0.03),
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
           query_type = c("gene", "metabolite"),
           analysis_type = c("enrich_pathway", "do_gsea"),
           result,
           database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
           sim.cutoff = 0.5,
           save_to_local = TRUE,
           path = "") {
    analysis_type <- match.arg(analysis_type)
    database <- match.arg(database)
    query_type <- match.arg(query_type)

    edge_data <-
      rbind(sim_matrix) %>%
      dplyr::rename(from = name1, to = name2) %>%
      dplyr::filter(sim > sim.cutoff)

    if (database == "go") {
      node_data <-
        rbind(result) %>%
        as.data.frame() %>%
        dplyr::filter(!(ID %in% attr(sim_matrix, "obsolete_terms"))) %>%
        dplyr::select(ID, everything()) %>%
        dplyr::rename(node = ID)
    } else {
      if (query_type == "gene") {
        node_data <-
          rbind(result) %>%
          as.data.frame() %>%
          dplyr::select(ID, everything()) %>%
          dplyr::rename(node = ID)
      } else if (query_type == "metabolite") {
        node_data <-
          result %>%
          dplyr::select(pathway_id, everything()) %>%
          dplyr::rename(node = pathway_id)
      }
        }

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

    ### Collect pathway merging result for gene enrichment result ====
    if (query_type == "gene") {
      result_with_module <-
        igraph::vertex_attr(graph_data) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>%
        dplyr::arrange(module, p_adjust)

      if (database == "go") {
        result_with_module <-
          result_with_module %>%
          dplyr::arrange(ONTOLOGY, module, p_adjust)
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
            pvalue = as.numeric(pvalue),
            Count = as.numeric(Count),
            p_adjust = as.numeric(p_adjust),
            qvalue = as.numeric(qvalue),
            degree = as.numeric(degree)
          ) %>%
          dplyr::arrange(p_adjust) %>%
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
    } else if (query_type == "metabolite") {
      ### Collect pathway merging result for metabolite enrichment result ====
      result_with_module <-
        igraph::vertex_attr(graph_data) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        dplyr::mutate(p_adjust = as.numeric(p_adjust)) %>%
        dplyr::arrange(module, p_adjust)

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

      if (analysis_type == "enrich_pathway") {
        module_result <-
          result_with_module %>%
          plyr::dlply(.variables = .(module)) %>%
          purrr::map(function(x) {
            if (nrow(x) == 1) {
              x <- x %>%
                dplyr::mutate(Description = pathway_name,
                              Count = 1)

              x$mapped_id =
                x$mapped_id %>%
                stringr::str_split(pattern = ";") %>%
                unlist() %>%
                unique() %>%
                paste(collapse = '/')

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

            x$mapped_id =
              x$mapped_id %>%
              stringr::str_split(pattern = ";") %>%
              unlist() %>%
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

      module_result <-
        module_result %>%
        dplyr::select(-c(describtion))
    }

    if (save_to_local) {
      save(result_with_module,
           file = file.path(path, "intermediate_data/result_with_module.RData"))
      save(graph_data, file = file.path(path, "intermediate_data/graph_data.RData"))
      save(module_result,
           file = file.path(path, "intermediate_data/module_result.RData"))
    }

    message("Done")

    list(
      graph_data = graph_data,
      module_result = module_result,
      result_with_module = result_with_module
    )
  }
