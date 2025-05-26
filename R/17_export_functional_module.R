# setwd(r4projects::get_project_wd())
# source("R/8_functional_module_class.R")
# source("R/6_utils.R")
# setwd("demo_data/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
# library(igraph)
#
# load("result/enriched_functional_module")
#
# object <-
#   enriched_functional_module
#
# export_functional_module(object = object)

#' Export functional-module results to CSV
#'
#' Writes every enrichment/table contained in a
#' `functional_module` object to **individual UTF-8 CSV files**.
#'
#' @param object  A `functional_module` object returned by the MAPA
#'   workflow, containing the enrichment and module results.
#' @param path    Directory in which to create the `.csv` files.
#'   Created if it does not already exist.
#'   **Default:** `"result"`.
#'
#' @details
#' * Files are written with \code{\link[readr]{write_csv}}, overwriting any
#'   existing files of the same name.
#' * Slots that are `NULL` are written as empty data-frames so the list
#'   of output files is always predictable.
#'
#' @return
#' Invisibly returns `NULL`; the function is called for its side-effect
#' of writing files.
#'
#' @author Xiaotao Shen (\email{shenxt1990@outlook.com})
#' @author Yifei Ge (\email{yifeii.ge@outlook.com})
#'
#'
#' @export

export_functional_module <-
  function(object,
           path = "result") {
    ## Check object and variable_info ====
    if (!is(object, "functional_module")) {
      stop("object should be functional_module class")
    }

    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    ## Save analysis results to Excel file ====
    ## pathway result
    if (!is.null(object@enrichment_go_result)) {
      enriched_pathway_go <-
        object@enrichment_go_result@result
    } else {
      enriched_pathway_go <-
        data.frame()
    }

    if (!is.null(object@enrichment_kegg_result)) {
      enriched_pathway_kegg <-
        object@enrichment_kegg_result@result
    } else {
      enriched_pathway_kegg <-
        data.frame()
    }


    if (!is.null(object@enrichment_reactome_result)) {
      enriched_pathway_reactome <-
        object@enrichment_reactome_result@result
    } else {
      enriched_pathway_reactome <-
        data.frame()
    }

    if (!is.null(object@enrichment_hmdb_result)) {
      enriched_pathway_hmdb <-
        object@enrichment_hmdb_result@result
    } else {
      enriched_pathway_hmdb <-
        data.frame()
    }

    if (!is.null(object@enrichment_metkegg_result)) {
      enriched_pathway_metkegg <-
        object@enrichment_metkegg_result@result
    } else {
      enriched_pathway_metkegg <-
        data.frame()
    }


    ### module result
    enriched_module_go <-
      object@merged_pathway_go$module_result
    if (is.null(enriched_module_go)) {
      enriched_module_go <-
        data.frame()
    }

    enriched_module_kegg <-
      object@merged_pathway_kegg$module_result
    if (is.null(enriched_module_kegg)) {
      enriched_module_kegg <-
        data.frame()
    }

    enriched_module_reactome <-
      object@merged_pathway_reactome$module_result
    if (is.null(enriched_module_reactome)) {
      enriched_module_reactome <-
        data.frame()
    }

    enriched_module_hmdb <-
      object@merged_pathway_hmdb$module_result
    if (is.null(enriched_module_hmdb)) {
      enriched_module_hmdb <-
        data.frame()
    }

    enriched_module_metkegg <-
      object@merged_pathway_metkegg$module_result
    if (is.null(enriched_module_metkegg)) {
      enriched_module_metkegg <-
        data.frame()
    }

    ### functional module result
    enriched_functional_module <-
      object@merged_module$functional_module_result
    if (is.null(enriched_functional_module)) {
      enriched_functional_module <-
        data.frame()
    }

    ### llm interpretation result
    if ("llm_interpret_module" %in% names(object@process_info)) {
      llm_module_interpretation <- extract_llm_module_data(object@llm_module_interpretation)
    } else {
      llm_module_interpretation <-
        data.frame()
    }

    query_type <- object@process_info$merge_pathways@parameter$query_type

    if (query_type == "gene") {
      sheet_names <- c(
        "enriched_pathway_go",
        "enriched_pathway_kegg",
        "enriched_pathway_reactome",
        "enriched_module_go",
        "enriched_module_kegg",
        "enriched_module_reactome",
        "enriched_functional_module",
        "llm_module_interpretation"
      )

      data_list <- list(
        enriched_pathway_go,
        enriched_pathway_kegg,
        enriched_pathway_reactome,
        enriched_module_go,
        enriched_module_kegg,
        enriched_module_reactome,
        enriched_functional_module,
        llm_module_interpretation
      )
    } else { # query_type == "metabolite"
      sheet_names <- c(
        "enriched_pathway_hmdb",
        "enriched_pathway_metkegg",
        "enriched_module_hmdb",
        "enriched_module_metkegg",
        "enriched_functional_module",
        "llm_module_interpretation"
      )

      data_list <- list(
        enriched_pathway_hmdb,
        enriched_pathway_metkegg,
        enriched_module_hmdb,
        enriched_module_metkegg,
        enriched_functional_module,
        llm_module_interpretation
      )
    }

    ## write each element to CSV using the same names
    purrr::walk2(
      data_list,
      sheet_names,
      ~ readr::write_csv(.x, file.path(path, paste0(.y, ".csv")))
    )

    # wb <-
    #   openxlsx::createWorkbook()
    # openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
    #
    # query_type <- object@process_info$merge_pathways@parameter$query_type
    # if (query_type == "gene") {
    #   for (idx in seq_len(8)) {
    #     openxlsx::addWorksheet(
    #       wb,
    #       sheetName = c(
    #         "enriched_pathway_go",
    #         "enriched_pathway_kegg",
    #         "enriched_pathway_reactome",
    #         "enriched_module_go",
    #         "enriched_module_kegg",
    #         "enriched_module_reactome",
    #         "enriched_functional_module",
    #         "llm_module_interpretation"
    #       )[idx],
    #       gridLines = TRUE
    #     )
    #     openxlsx::freezePane(wb,
    #                          sheet = idx,
    #                          firstRow = TRUE,
    #                          firstCol = TRUE)
    #     openxlsx::writeData(
    #       wb,
    #       sheet = idx,
    #       x = list(
    #         enriched_pathway_go,
    #         enriched_pathway_kegg,
    #         enriched_pathway_reactome,
    #         enriched_module_go,
    #         enriched_module_kegg,
    #         enriched_module_reactome,
    #         enriched_functional_module,
    #         llm_module_interpretation
    #       )[[idx]],
    #       colNames = TRUE,
    #       rowNames = FALSE
    #     )
    #   }
    # } else if (query_type == "metabolite") {
    #   for (idx in seq_len(6)) {
    #     openxlsx::addWorksheet(
    #       wb,
    #       sheetName = c(
    #         "enriched_pathway_hmdb",
    #         "enriched_pathway_metkegg",
    #         "enriched_module_hmdb",
    #         "enriched_module_metkegg",
    #         "enriched_functional_module",
    #         "llm_module_interpretation"
    #       )[idx],
    #       gridLines = TRUE
    #     )
    #     openxlsx::freezePane(wb,
    #                          sheet = idx,
    #                          firstRow = TRUE,
    #                          firstCol = TRUE)
    #     openxlsx::writeData(
    #       wb,
    #       sheet = idx,
    #       x = list(
    #         enriched_pathway_hmdb,
    #         enriched_pathway_metkegg,
    #         enriched_module_hmdb,
    #         enriched_module_metkegg,
    #         enriched_functional_module,
    #         llm_module_interpretation
    #       )[[idx]],
    #       colNames = TRUE,
    #       rowNames = FALSE
    #     )
    #   }
    # }
    #
    # openxlsx::saveWorkbook(
    #   wb,
    #   file = file.path(path, "Functional_module_analysis_results.xlsx"),
    #   overwrite = TRUE
    # )
}
