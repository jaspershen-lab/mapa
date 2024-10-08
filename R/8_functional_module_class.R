#'An S4 class to represent enriched pathways
#' @docType class
#' @slot variable_info variable information contains marker information.
#' @slot enrichment_go_result ANY Object storing GO enrichment results.
#' @slot enrichment_kegg_result ANY Object storing KEGG enrichment results.
#' @slot enrichment_reactome_result ANY Object storing Reactome enrichment results.
#' @slot merged_pathway_go list List containing merged GO pathway information.
#' @slot merged_pathway_kegg list List containing merged KEGG pathway information.
#' @slot merged_pathway_reactome list List containing merged Reactome pathway information.
#' @slot merged_module list List containing merged modules.
#' @slot process_info list List containing information about the processes.
#'
#' @exportClass functional_module

setClass(
  Class = "functional_module",
  representation(
    variable_info = "data.frame",
    enrichment_go_result = "ANY",
    enrichment_kegg_result = "ANY",
    enrichment_reactome_result = "ANY",
    merged_pathway_go = "list",
    merged_pathway_kegg = "list",
    merged_pathway_reactome = "list",
    merged_module = "list",
    process_info = "list"
  )
)

# a validity method to enforce that `enrichment_go_result` is either a enrichResult, gseaResult or a NULL
setValidity("functional_module", function(object) {
  if (!(
    is(object@enrichment_go_result, "enrichResult") ||
    is(object@enrichment_go_result, "gseaResult") ||
    is.null(object@enrichment_go_result)
  )) {
    return(
      "The 'enrichment_go_result' slot must be either a enrichResult, gseaResult or a NULL."
    )
  }
  if (!(
    is(object@enrichment_kegg_result, "enrichResult") ||
    is(object@enrichment_kegg_result, "gseaResult") ||
    is.null(object@enrichment_kegg_result)
  )) {
    return(
      "The 'enrichment_kegg_result' slot must be either a enrichResult, gseaResult or a NULL."
    )
  }
  if (!(
    is(object@enrichment_reactome_result, "enrichResult") ||
    is(object@enrichment_reactome_result, "gseaResult") ||
    is.null(object@enrichment_reactome_result)
  )) {
    return(
      "The 'enrichment_reactome_result' slot must be either a enrichResult, gseaResult or a NULL."
    )
  }
  TRUE
})

setMethod(
  f = "show",
  signature = "functional_module",
  definition = function(object) {
    enrichment_go_result <-
      try(object@enrichment_go_result, silent = TRUE)
    enrichment_kegg_result <-
      try(object@enrichment_kegg_result, silent = TRUE)
    enrichment_reactome_result <-
      try(object@enrichment_reactome_result, silent = TRUE)
    parameter <-
      try(object@parameter, silent = TRUE)

    ###pathway enrichment or gsea
    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
    } else{
      analysis_type <- "do_gsea"
    }

    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("Analysis method:", analysis_type, "\n"))
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))

    message(crayon::green("-----------Variable information------------"))
    message(crayon::green(nrow(object@variable_info)),
            crayon::green(" features/markers in total"))


    p.adjust.cutoff <-
      tryCatch(
        object@process_info$merge_pathways@parameter$p.adjust.cutoff.go,
        error = function(e) {
          tryCatch(
            object@process_info$enrich_pathway@parameter$pvalueCutoff,
            error = function(e) {
              object@process_info$do_gse@parameter$pvalueCutoff
            }
          )
        }
      )

    count.cutoff <-
      tryCatch(
        object@process_info$merge_pathways@parameter$count.cutoff.go,
        error = function(e) {
          tryCatch(
            object@process_info$enrich_pathway@parameter$pvalueCutoff,
            error = function(e) {
              object@process_info$do_gse@parameter$pvalueCutoff
            }
          )
        }
      )

    message(crayon::green("-----------GO------------"))
    if (is.null(enrichment_go_result)) {
      message(crayon::green('No GO results'))
    } else{
      if (analysis_type == "enrich_pathway") {
        message(
          crayon::green(
            nrow(
              enrichment_go_result@result %>%
                dplyr::filter(
                  p.adjust < p.adjust.cutoff &
                    Count > count.cutoff &
                    ONTOLOGY != "CC"
                )
            ),
            "GO terms (BP and MF) with p.adjust <",
            p.adjust.cutoff,
            "and Count >",
            count.cutoff
          )
        )
      } else{
        message(crayon::green(
          nrow(
            enrichment_go_result@result %>%
              dplyr::filter(p.adjust < p.adjust.cutoff &
                              ONTOLOGY != "CC")
          ),
          "GO terms (BP and MF) with p.adjust <",
          p.adjust.cutoff
        ))
      }
    }

    if (length(object@merged_pathway_go) == 0) {
      message(crayon::green('No GO modules'))
    } else{
      message(crayon::green(
        length(object@merged_pathway_go$module_result$module),
        "GO modules"
      ))
    }

    message(crayon::green("-----------KEGG------------"))
    if (is.null(enrichment_kegg_result)) {
      message(crayon::green('No KEGG results'))
    } else{
      if (analysis_type == "enrich_pathway") {
        message(
          crayon::green(
            nrow(
              enrichment_kegg_result@result %>%
                dplyr::filter(p.adjust < p.adjust.cutoff &
                                Count > count.cutoff)
            ),
            "KEGG pathways with p.adjust <",
            p.adjust.cutoff,
            "and Count >",
            count.cutoff
          )
        )
      } else{
        message(crayon::green(
          nrow(
            enrichment_kegg_result@result %>%
              dplyr::filter(p.adjust < p.adjust.cutoff)
          ),
          "KEGG pathways with p.adjust <",
          p.adjust.cutoff
        ))
      }
    }

    if (length(object@merged_pathway_kegg) == 0) {
      message(crayon::green('No KEGG modules'))
    } else{
      message(crayon::green(
        length(object@merged_pathway_kegg$module_result$module),
        "KEGG modules"
      ))
    }

    message(crayon::green("-----------Reactome------------"))
    if (is.null(enrichment_reactome_result)) {
      message(crayon::green('No Reactome results'))
    } else{
      if (analysis_type == "enrich_pathway") {
        message(
          crayon::green(
            nrow(
              enrichment_reactome_result@result %>%
                dplyr::filter(p.adjust < p.adjust.cutoff &
                                Count > count.cutoff)
            ),
            "Reactome pathways with p.adjust <",
            p.adjust.cutoff,
            "and Count >",
            count.cutoff
          )
        )
      } else{
        message(crayon::green(
          nrow(
            enrichment_reactome_result@result %>%
              dplyr::filter(p.adjust < p.adjust.cutoff)
          ),
          "Reactome pathways with p.adjust <",
          p.adjust.cutoff
        ))
      }
    }

    if (length(object@merged_pathway_reactome) == 0) {
      message(crayon::green('No Reactome modules'))
    } else{
      message(crayon::green(
        length(object@merged_pathway_reactome$module_result$module),
        "Reactome modules"
      ))
    }

    message(crayon::green("-----------Functional modules------------"))
    if (length(object@merged_module) == 0) {
      message(crayon::green('No Functional modules'))
    } else{
      message(crayon::green(
        length(object@merged_module$functional_module_result$module),
        "Functional modules"
      ))
    }

    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("Processing information\n"))
    if (.hasSlot(object = object, name = "process_info") &
        length(object@process_info) != 0) {
      process_info <- object@process_info
      cat(crayon::green(length(process_info), "processings in total\n"))
      if (length(process_info) > 5) {
        cat(crayon::green("Latest 3 processings show\n"))
        process_info <- tail(process_info, 3)
      }
      for (idx in seq_along(process_info)) {
        cat(crayon::green(names(process_info)[idx], paste(rep("-", 10), collapse = ""), "\n"))
        if (length(process_info[[idx]]) == 1) {
          data.frame(
            "Package" = process_info[[idx]]@pacakge_name,
            "Function used" = process_info[[idx]]@function_name,
            "Time" = process_info[[idx]]@time
          ) %>%
            print()
        } else{
          data.frame(
            "Package" = process_info[[idx]] %>% lapply(function(x)
              x@pacakge_name) %>% unlist(),
            "Function used" = process_info[[idx]] %>% lapply(function(x)
              x@function_name) %>% unlist(),
            "Time" = process_info[[idx]] %>% lapply(function(x)
              as.character(x@time)) %>% unlist()
          ) %>%
            print()
        }
      }
    } else{
      cat(crayon::red("There are no processing for your data.\n"))
    }
  }
)
