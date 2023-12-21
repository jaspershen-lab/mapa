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
# intepretation_result <-
#   inteprete_pathways(
#     object = object,
#     p.adjust.cutoff = 0.05,
#     disease = "pregnancy",
#     count.cutoff = 5,
#     top_n = 3
#   )


#' Interpret Enriched Pathways in Relation to a Specific Disease
#'
#' This function interprets the enriched pathways from an object of class
#' 'functional_module' in relation to a specific disease. It filters pathways
#' based on p-value and count thresholds and returns their functions and
#' relationships with a given disease.
#'
#' @param object An object of class 'functional_module'. The function will
#'               interpret pathways based on this object.
#' @param disease A character string specifying the disease of interest.
#'                Default is "pregnancy".
#' @param p.adjust.cutoff A numeric value for the adjusted p-value cutoff
#'                        for pathway enrichment. Default is 0.05.
#' @param count.cutoff A numeric value specifying the minimum count cutoff
#'                     for pathway enrichment. Default is 5.
#' @param top_n An integer indicating the number of top pathways to consider
#'              based on the p.adjust.cutoff and count.cutoff. Default is 3.
#'
#' @return A list containing two elements: 'functions' and
#'         'relationship_with_disease'. Each element is a list of character
#'         strings, with each string providing information for a specific
#'         pathway in markdown format.
#'
#' @details The function first checks if the input object is of the correct
#'          class. It then filters the enriched pathways based on the provided
#'          p-value and count thresholds. For each of the top pathways, the
#'          function retrieves information about their functions and their
#'          relationships with the specified disease. This information is
#'          gathered using an external service ('request_chatgpt_response').
#'
#'
#' @importFrom dplyr filter arrange
#' @importFrom stringr str_split
#'
#' @export

inteprete_pathways <-
  function(object,
           disease = "pregnancy",
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           top_n = 3) {
    ###check object and variable_info
    if (!is(object, "functional_module")) {
      stop("object should be functional_module class")
    }

    enriched_functional_module <-
      object@merged_module$functional_module_result

    if (is.null(enriched_functional_module)) {
      enriched_functional_module <-
        data.frame()
    }

    enriched_functional_module <-
      enriched_functional_module %>%
      dplyr::filter(p.adjust < p.adjust.cutoff,
                    Count >= count.cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      head(top_n)

    functions_list <-
      relationship_with_disease_list <-
      vector(mode = "list",
             length = nrow(enriched_functional_module))

    names(functions_list) <-
      names(relationship_with_disease_list) <-
      enriched_functional_module$module

    for (idx in seq_len(nrow(enriched_functional_module))) {
      message("Module: ", enriched_functional_module$module[idx], "\n")
      pathway_names <-
        enriched_functional_module$Description[idx] %>%
        stringr::str_split(";") %>%
        `[[`(1)

      pathway_ids <-
        enriched_functional_module$pathway_id[idx] %>%
        stringr::str_split(";") %>%
        `[[`(1)

      pathway_info <-
        paste0(pathway_names,
               " (ID is: ",
               pathway_ids,
               ")") %>%
        paste0(collapse = "; ")

      functions <-
        request_chatgpt_response(
          prompt = paste0(
            "What is the function of these pathways: ",
            pathway_info,
            "?\n The output result should be listed, and using the markdown format.",
            "\n Please give the links of the references."
          )
        )

      relationship_with_disease <-
        request_chatgpt_response(
          prompt = paste0(
            "What is the relationship between ",
            disease,
            " and these pathways: ",
            pathway_info,
            "?\n The output result should be listed, and using the markdown format.",
            "\n Please give the links of the references."
          )
        )

      relationship_with_disease_list[[idx]] <-
        relationship_with_disease

      functions_list[[idx]] <-
        functions
    }

    return(
      list(
        functions = functions_list,
        relationship_with_disease = relationship_with_disease_list
      )
    )

    message("Done.")

  }
