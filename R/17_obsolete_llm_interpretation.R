#' # setwd(r4projects::get_project_wd())
#' # source("R/8_functional_module_class.R")
#' # source("R/6_utils.R")
#' # setwd("demo_data/")
#' # library(clusterProfiler)
#' # library(org.Hs.eg.db)
#' # library(ReactomePA)
#' # library(igraph)
#' #
#' # load("result/enriched_functional_module")
#' #
#' # object <-
#' #   enriched_functional_module
#' #
#' # interpretation_result <-
#' #   interpret_pathways(
#' #     object = object,
#' #     p.adjust.cutoff = 0.05,
#' #     disease = "pregnancy",
#' #     count.cutoff = 5,
#' #     top_n = 3
#' #   )
#' #
#' # save(interpretation_result,
#' #      file = "interpretation_result.rda")
#'
#' #' Interpret Pathways in Functional Modules
#' #'
#' #' This function interprets the pathways in functional modules, particularly focused on a specific disease.
#' #' It filters and arranges functional modules based on p-value adjustment and count cutoffs,
#' #' then retrieves and formats information about the top pathways related to the disease.
#' #'
#' #' @param object An object of class "functional_module". This object contains the functional module data that needs to be interpreted.
#' #' @param disease A character string specifying the disease to focus on. Default is "pregnancy".
#' #' @param p.adjust.cutoff A numeric value representing the p-value adjustment cutoff for filtering functional modules. Default is 0.05.
#' #' @param count.cutoff A numeric value representing the count cutoff for filtering functional modules. Default is 5.
#' #' @param top_n An integer indicating the number of top pathways to return. Default is 3.
#' #'
#' #' @return A character string or list containing the interpreted pathway information. If no enriched functional modules are found,
#' #' a message is returned indicating that no enriched functional modules are found and suggesting different cutoffs.
#' #'
#' #' @importFrom dplyr filter arrange
#' #' @importFrom stringr str_split
#' #' @export
#'
#' interpret_pathways <-
#'   function(object,
#'            disease = "pregnancy",
#'            p.adjust.cutoff = 0.05,
#'            count.cutoff = 5,
#'            top_n = 3) {
#'     ###check object and variable_info
#'     if (!is(object, "functional_module")) {
#'       stop("object should be functional_module class")
#'     }
#'
#'     enriched_functional_module <-
#'       object@merged_module$functional_module_result
#'
#'     if (is.null(enriched_functional_module)) {
#'       enriched_functional_module <-
#'         data.frame()
#'     }
#'
#'     enriched_functional_module <-
#'       enriched_functional_module %>%
#'       dplyr::filter(p.adjust < p.adjust.cutoff,
#'                     Count >= count.cutoff) %>%
#'       dplyr::arrange(p.adjust) %>%
#'       head(top_n)
#'
#'     if (nrow(enriched_functional_module) == 0) {
#'       return(
#'         "No enriched functional modules are found. Please try different p.adjust.cutoff and count.cutoff."
#'       )
#'     }
#'
#'     functions_list <-
#'       vector(mode = "list",
#'              length = nrow(enriched_functional_module))
#'
#'     names(functions_list) <-
#'       enriched_functional_module$module
#'
#'     for (idx in seq_len(nrow(enriched_functional_module))) {
#'       message("Module: ", enriched_functional_module$module[idx], "\n")
#'       pathway_names <-
#'         enriched_functional_module$Description[idx] %>%
#'         stringr::str_split(";") %>%
#'         `[[`(1)
#'
#'       pathway_ids <-
#'         enriched_functional_module$pathway_id[idx] %>%
#'         stringr::str_split(";") %>%
#'         `[[`(1)
#'
#'       pathway_info <-
#'         paste0(pathway_names,
#'                " (ID is: ",
#'                pathway_ids,
#'                ")") %>%
#'         paste0(collapse = "; ")
#'
#'       functions <-
#'         request_chatgpt_response(
#'           prompt = paste0(
#'             "I am interested in understanding the role of these pathways: ",
#'             pathway_info,
#'             ", in relation to ",
#'             disease,
#'             ".\nHere are my requirements for the response:",
#'             "\n1. Present the findings in markdown format.",
#'             "\n2. Organize the information in a list format, but please avoid using markdown headings.",
#'             "\n3. Give the response directly, don't say anything else.",
#'             "\n4. The response should be simple, conise, and easy to understand, no more than 300 words, and at least 100 words.",
#'             "\n5. Include references for the information provided.
#'             The references should be formatted similarly to those in Nature Journal and should be listed at the end of the response.",
#'             "\n6. In case there is a lack of published papers directly linking this pathway to ",
#'             disease,
#'             " I would appreciate an informed hypothesis based on the known functions of the pathway."
#'           )
#'         )
#'
#'       functions_list[[idx]] <-
#'         functions
#'     }
#'
#'     message("Done.")
#'
#'     for (idx in seq_along(functions_list)) {
#'       functions_list[[idx]] <-
#'         paste0("## ", names(functions_list)[idx], "\n\n", functions_list[[idx]])
#'     }
#'
#'     functions_list <-
#'       paste(functions_list, collapse = "\n\n")
#'
#'     return(functions_list)
#'   }
