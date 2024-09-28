#' Filter Functional Modules
#'
#' This function filters functional modules, modules, pathways, or molecules in the given object based on specified IDs.
#'
#' @param object An object containing data on functional modules, modules, pathways, or molecules.
#' @param level A character string specifying the level at which filtering should occur. Valid levels are "functional_module", "module", "pathway", and "molecule".
#' @param remain_id A vector of IDs to retain in the filtering process.
#'
#' @return Returns the modified object with filtered data based on the given level and IDs.
#'
#' @importFrom dplyr filter
#' @importFrom stringr str_split
#' @export
#'
#' @note This function modifies the input object in-place,
#' depending on the specified level and IDs to remain.
#'

filter_functional_module <-
  function(object,
           level = c("functional_module",
                     "module",
                     "pathway",
                     "molecule"),
           remain_id) {
    if (missing(remain_id)) {
      return(object)
    }

    if (level == "functional_module") {
      object@merged_module$functional_module_result <-
        dplyr::filter(object@merged_module$functional_module_result,
                      module %in% remain_id)

      object@merged_module$result_with_module <-
        dplyr::filter(object@merged_module$result_with_module,
                      module %in% remain_id)

    }

    if (level == "module") {
      ###remove functional modules
      remain_idx <-
        which(unlist(
          lapply(object@merged_module$functional_module_result$module_content,
                 function(x) {
                   any(stringr::str_split(x, ";")[[1]] == remain_id)
                 })
        ))

      object@merged_module$functional_module_result <-
        object@merged_module$functional_module_result[remain_idx, , drop = FALSE]

      object@merged_module$result_with_module <-
        dplyr::filter(object@merged_module$result_with_module,
                      node %in% remain_id)

      ####remove GO modules
      if (length(object@merged_pathway_go) > 0) {
        object@merged_pathway_go$module_result <-
          dplyr::filter(object@merged_pathway_go$module_result,
                        module %in% remain_id)

        object@merged_pathway_go$result_with_module <-
          dplyr::filter(object@merged_pathway_go$result_with_module,
                        module %in% remain_id)
      }

      ####remove KEGG modules
      if (length(object@merged_pathway_kegg) > 0) {
        object@merged_pathway_kegg$module_result <-
          dplyr::filter(object@merged_pathway_kegg$module_result,
                        module %in% remain_id)

        object@merged_pathway_kegg$result_with_module <-
          dplyr::filter(object@merged_pathway_kegg$result_with_module,
                        module %in% remain_id)
      }

      ####remove Reactome modules
      if (length(object@merged_pathway_reactome) > 0) {
        object@merged_pathway_reactome$module_result <-
          dplyr::filter(object@merged_pathway_reactome$module_result,
                        module %in% remain_id)

        object@merged_pathway_reactome$result_with_module <-
          dplyr::filter(object@merged_pathway_reactome$result_with_module,
                        module %in% remain_id)
      }
    }

    return(object)

  }
