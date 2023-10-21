# setwd(r4projects::get_project_wd())
# source("R/8-functional_module_class.R")
# source("R/6-utils.R")
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


#' Export Functional Module Analysis Results to Excel
#'
#' This function exports the results of functional module analysis to an Excel file.
#' The function takes an object of class `functional_module` and writes multiple
#' sheets corresponding to various types of enrichment results.
#'
#' @param object An object of class `functional_module`. This object should contain
#'        enrichment results and other relevant data.
#' @param path The directory path where the resulting Excel file will be saved.
#'        Defaults to "result".
#'
#' @return Saves an Excel file to the specified path.
#'         The Excel file contains sheets for GO, KEGG, and Reactome enrichment results
#'         for pathways and modules, as well as enriched functional modules.
#'
#' @author Xiaotao Shen, \email{shenxt1990@outloo.com}
#'
#' @export

export_functional_module <-
  function(object,
           path = "result") {
    ###check object and variable_info
    if (!is(object, "functional_module")) {
      stop("object should be functional_module class")
    }

    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    enriched_pathway_go <-
      object@enrichment_go_result@result

    enriched_pathway_kegg <-
      object@enrichment_kegg_result@result

    enriched_pathway_reactome <-
      object@enrichment_reactome_result@result

    enriched_module_go <-
      object@merged_pathway_go$module_result

    enriched_module_kegg <-
      object@merged_pathway_kegg$module_result

    enriched_module_reactome <-
      object@merged_pathway_reactome$module_result

    enriched_functional_module <-
      object@merged_module$functional_module_result

    wb <-
      openxlsx::createWorkbook()
    openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")

    for (idx in seq_len(7)) {
      openxlsx::addWorksheet(
        wb,
        sheetName = c(
          "enriched_pathway_go",
          "enriched_pathway_kegg",
          "enriched_pathway_reactome",
          "enriched_module_go",
          "enriched_module_kegg",
          "enriched_module_reactome",
          "enriched_functional_module"
        )[idx],
        gridLines = TRUE
      )
      openxlsx::freezePane(wb,
                           sheet = idx,
                           firstRow = TRUE,
                           firstCol = TRUE)
      openxlsx::writeDataTable(
        wb,
        sheet = idx,
        x = list(
          enriched_pathway_go,
          enriched_pathway_kegg,
          enriched_pathway_reactome,
          enriched_module_go,
          enriched_module_kegg,
          enriched_module_reactome,
          enriched_functional_module
        )[[idx]],
        colNames = TRUE,
        rowNames = FALSE
      )
    }
    openxlsx::saveWorkbook(
      wb,
      file = file.path(path, "Functional_module_analysis_results.xlsx"),
      overwrite = TRUE
    )
  }
