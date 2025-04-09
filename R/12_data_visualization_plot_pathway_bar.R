#######pathway enrichment
# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# load("enriched_functional_module.rda")
#
# object <-
#   enriched_functional_module
# object <-
#   enriched_pathways
#
# library(showtext)
# showtext_auto(enable = TRUE)
#
# parameters
# object = enriched_functional_module
# top_n = 10
# x_axis_name = "RichFactor" #"qscore", "RichFactor", or "FoldEnrichment"
# y_label_width = 30
# level = "functional_module"
# line_type = "straight"
# database = c("go", "kegg", "reactome", "hmdb")
# database_color =
#   c(
#     GO = "#eeca40",
#     KEGG = "#fd7541",
#     Reactome = "#23b9c7",
#     HMDB = "#7998ad"
#   )
# p.adjust.cutoff = 0.05
# count.cutoff = 5

# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   x_axis_name = "RichFactor", #"qscore", "RichFactor", or "FoldEnrichment"
#   y_label_width = 30,
#   level = "functional_module",
#   line_type = "straight"
# )
#
#
######GSEA analysis
# setwd(r4projects::get_project_wd())
# setwd("demo_data/covid_data/")
# load("result/enriched_functional_module")
#
# object <-
#   enriched_functional_module
#
# object <- gsea_enriched_functional_module
#
#
# library(showtext)
# showtext_auto(enable = TRUE)
#
# plot_pathway_bar(
#   object = gsea_pathways,
#   top_n = 10,
#   y_label_width = 30,
#   level = "pathway",
#   translation = FALSE,
#   line_type = "straight"
# )
#
# plot_pathway_bar(
#   object = object,
#   top_n = 10,
#   y_label_width = 30,
#   level = "module",
#   line_type = "meteor"
# )
#
# plot_pathway_bar(
#   object = object,
#   top_n = 10,
#   y_label_width = 30,
#   level = "functional_module"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_label_width = 30,
#   level = "functional_module",
#   line_type = "meteor"
# )
#
# Metabolite enrichment result
# plot_pathway_bar(object =  enriched_functional_modules,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  level = "pathway",
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 5,
#                  database = c("kegg", "hmdb"))
#

# plot_pathway_bar(object = enriched_functional_modules,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  level = "module",
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 5,
#                  database = c("hmdb", "kegg"))

# plot_pathway_bar(object = enriched_functional_modules,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  level = "functional_module",
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 5,
#                  database = c("hmdb", "kegg"))

#' Plot Enriched Pathway/Module Bar Chart
#'
#' This function generates a bar chart to visualize enriched pathways, modules, or functional modules.
#' The bars are colored according to the database of origin (GO, KEGG, Reactome, or HMDB).
#' Different metrics are displayed on the x-axis depending on the type of enrichment analysis
#' (ORA or GSEA) and the query type (gene or metabolite).
#'
#' @param object An object of class "functional_module" containing the enrichment results and other relevant data.
#' @param top_n An integer specifying the top N pathways to display. Default is 10.
#' @param x_axis_name A character string indicating the metric to use for the x-axis.
#'   For gene ORA, options include \code{"qscore"} (i.e. \eqn{-\log10} of the FDR-adjusted p-values),
#'   \code{"RichFactor"}, or \code{"FoldEnrichment"}.
#'   For metabolite enrichment, this parameter should be \code{"qscore"}.
#'   For GSEA results, this parameter should be \code{"NES"}.
#' @param y_label_width An integer specifying the width of the Y-axis labels for text wrapping. Default is 50.
#' @param translation Logical indicating whether to use translated descriptions. If TRUE,
#'   the function will use the translated descriptions from the `translate_language` function. Default is FALSE.
#' @param level A character string specifying the level of analysis.
#'   One of "pathway", "module", or "functional_module". Default is "pathway".
#'   For biotext embedding results, only "pathway" or "functional_module" can be used.
#' @param line_type A character string specifying the type of line to use for the bar chart.
#'   One of "straight" or "meteor". Default is "straight".
#' @param p.adjust.cutoff A numeric value for the FDR adjusted P-value cutoff. Default is 0.05.
#' @param count.cutoff A numeric value for the minimum count of genes/metabolites to include. Default is 5.
#' @param database_color A named vector containing the colors for different databases.
#'   Default colors are provided for GO, KEGG, Reactome, and HMDB.
#' @param database A character vector indicating which databases to include in the plot.
#'   For gene enrichment analysis, valid options are \code{"go"}, \code{"kegg"}, and/or \code{"reactome"}.
#'   For metabolite enrichment analysis, valid options are \code{"hmdb"} and/or \code{"kegg"}.
#'   Default includes all databases.
#'
#' @return A ggplot object representing the enrichment bar chart.
#'
#' @import ggplot2
#' @importFrom ggforce geom_link
#' @importFrom dplyr mutate filter select rename full_join arrange
#' @importFrom stringr str_wrap str_split str_detect
#' @importFrom purrr map map_chr
#' @importFrom methods is
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @examples
#' \dontrun{
#' data("enriched_functional_module", package = "mapa")
#' # Basic pathway plot
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_label_width = 30,
#'   level = "pathway"
#' )
#'
#' # Using meteor line type
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_label_width = 30,
#'   level = "pathway",
#'   line_type = "meteor"
#' )
#'
#' # Module level plot
#' plot_pathway_bar(
#'   enriched_functional_module,
#'   top_n = 10,
#'   level = "module"
#' )
#'
#' # Functional module plot with custom p-value cutoff
#' plot_pathway_bar(
#'   enriched_functional_module,
#'   top_n = 10,
#'   level = "functional_module",
#'   p.adjust.cutoff = 0.05
#' )
#'
#' # Specifying x-axis metric for gene enrichment
#' plot_pathway_bar(
#'   enriched_functional_module,
#'   top_n = 10,
#'   x_axis_name = "RichFactor",
#'   level = "pathway"
#' )
#' }
#'
#' @export

plot_pathway_bar <-
  function(object,
           top_n = 10,
           x_axis_name = NULL,
           y_label_width = 50,
           translation = FALSE,
           level = c("pathway", "module", "functional_module"),
           line_type = c("straight", "meteor"),
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database_color =
             c(
               GO = "#eeca40",
               KEGG = "#fd7541",
               Reactome = "#23b9c7",
               HMDB = "#7998ad"
             ),
           database = c("go", "kegg", "reactome", "hmdb")) {

    level <- match.arg(level)
    sim_method <- object@process_info$merge_pathways@function_name
    if (sim_method == "get_bioembedsim()" & level == "module") {
      stop("level `module` can not be applied to results generated by get_bioembedsim() and merge_pathways_bioembedsim() since these functions directly generate functional modules from enriched pathways. Please use either `patwhay` or `functional_module` level.")
    }
    line_type <- match.arg(line_type)

    if ("enrich_pathway" %in% names(object@process_info)) {
      query_type <- object@process_info$enrich_pathway@parameter$query_type
      analysis_type <- "enrich_pathway"
      if (query_type == "gene") {
        if (is.null(x_axis_name)){
          stop("Please provide a variable for x axis.")
        } else {
          x_axis_name <- match.arg(x_axis_name, choices = c("qscore", "RichFactor", "FoldEnrichment"))
        }
      } else if (query_type == "metabolite") {
        if (x_axis_name != "qscore") {
          message("For metabolite enrichment, please use qscore(adjusted p value) as x axis.")
        }
      }
    } else {
      query_type <- "gene"
      analysis_type <- "do_gsea"
      x_axis_name <- "NES"
    }

    if (translation) {
    #   if (all(names(object@process_info) != "translate_language")) {
    #     stop("Please use the 'translate_language' function to translate first.")
    #   } else{
    #     if (length(object@enrichment_go_result) > 0) {
    #       object@enrichment_go_result@result <-
    #         object@enrichment_go_result@result %>%
    #         dplyr::select(-Description) %>%
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #
    #     if (length(object@enrichment_kegg_result) > 0) {
    #       object@enrichment_kegg_result@result <-
    #         object@enrichment_kegg_result@result %>%
    #         dplyr::select(-Description) %>%
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #     if (length(object@enrichment_reactome_result) > 0) {
    #       object@enrichment_reactome_result@result <-
    #         object@enrichment_reactome_result@result %>%
    #         dplyr::select(-Description) %>%
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #
    #     if (length(object@merged_pathway_go) > 0) {
    #       object@merged_pathway_go$module_result <-
    #         object@merged_pathway_go$module_result %>%
    #         dplyr::select(-c(Description, module_annotation)) %>%
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_pathway_kegg) > 0) {
    #       object@merged_pathway_kegg$module_result <-
    #         object@merged_pathway_kegg$module_result %>%
    #         dplyr::select(-c(Description, module_annotation)) %>%
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_pathway_reactome) > 0) {
    #       object@merged_pathway_reactome$module_result <-
    #         object@merged_pathway_reactome$module_result %>%
    #         dplyr::select(-c(Description, module_annotation)) %>%
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_module) > 0) {
    #       object@merged_module$functional_module_result <-
    #         object@merged_module$functional_module_result %>%
    #         dplyr::select(-c(Description, module_annotation)) %>%
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #   }
    }

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    if (level == "pathway") {
      ## For gene enrichment analysis
      if (query_type == "gene") {
        enrichment_go_result <-
          tryCatch(
            object@enrichment_go_result@result %>%
              dplyr::mutate(class = "GO"),
            error = function(e) {
              NULL
            }
          )
        enrichment_kegg_result <-
          tryCatch(
            object@enrichment_kegg_result@result %>%
              dplyr::mutate(class = "KEGG"),
            error = function(e) {
              NULL
            }
          )
        enrichment_reactome_result <-
          tryCatch(
            object@enrichment_reactome_result@result %>%
              dplyr::mutate(class = "Reactome"),
            error = function(e) {
              NULL
            }
          )

        if (is.null(enrichment_go_result) &
            is.null(enrichment_kegg_result) &
            is.null(enrichment_reactome_result)) {
          stop("No enriched pathways for all the datasets")
        }

        if (!is.null(enrichment_reactome_result)) {
          if (!is.null(enrichment_kegg_result)) {
            if (ncol(enrichment_reactome_result) != ncol(enrichment_kegg_result)) {
              enrichment_reactome_result <-
                data.frame(
                  enrichment_reactome_result,
                  category = NA,
                  subcategory = NA
                ) %>%
                dplyr::select(category, subcategory, dplyr::everything())

              enrichment_kegg_result <-
                enrichment_kegg_result %>%
                dplyr::select(category, subcategory, dplyr::everything())
            }
          }
        }

        temp_data <-
          rbind(enrichment_reactome_result, enrichment_kegg_result)

        if (is.null(temp_data)) {
          temp_data <-
            enrichment_go_result
        } else{
          if (!is.null(enrichment_go_result)) {
            temp_data <-
              enrichment_go_result %>%
              dplyr::full_join(temp_data, by = intersect(colnames(.), colnames(temp_data)))

            if (analysis_type == "do_gsea") {
              temp_data$Count <-
                unlist(lapply(
                  stringr::str_split(temp_data$core_enrichment, "/"),
                  length
                ))
            }

            temp_data <-
              temp_data %>%
              dplyr::filter(p.adjust < p.adjust.cutoff &
                              Count > count.cutoff)
          }
        }
      } else if (query_type == "metabolite") {
        enrichment_hmdb_result <-
          tryCatch(
            object@enrichment_hmdb_result@result %>%
              dplyr::mutate(class = "HMDB"),
            error = function(e) {
              NULL
            }
          )

        enrichment_metkegg_result <-
          tryCatch(
            object@enrichment_metkegg_result@result %>%
              dplyr::mutate(class = "KEGG"),
            error = function(e) {
              NULL
            }
          )

        if (is.null(enrichment_hmdb_result) &
            is.null(enrichment_metkegg_result)) {
          stop("No enriched pathways for all the datasets.")
        }

        temp_data <-
          rbind(enrichment_hmdb_result, enrichment_metkegg_result)

        colnames(temp_data)[6] <- "p.adjust"
        colnames(temp_data)[10] <- "Count"

        temp_data <-
          temp_data %>%
          dplyr::filter(p.adjust < p.adjust.cutoff &
                          Count > count.cutoff)
      }
    }

    if (level == "module") {
      if (query_type == "gene") {
        if (length(object@merged_pathway_go) == 0 &
            length(object@merged_pathway_kegg) == 0 &
            length(object@merged_pathway_reactome) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else{
          module_result_go <-
            tryCatch(
              object@merged_pathway_go$module_result %>%
                dplyr::mutate(class = "GO"),
              error = function(e) {
                NULL
              }
            )
          module_result_kegg <-
            tryCatch(
              object@merged_pathway_kegg$module_result %>%
                dplyr::mutate(class = "KEGG"),
              error = function(e) {
                NULL
              }
            )

          module_result_reactome <-
            tryCatch(
              object@merged_pathway_reactome$module_result %>%
                dplyr::mutate(class = "Reactome"),
              error = function(e) {
                NULL
              }
            )

          if (!is.null(module_result_reactome)) {
            if (!is.null(module_result_kegg)) {
              if (ncol(module_result_reactome) != ncol(module_result_kegg)) {
                module_result_reactome <-
                  data.frame(module_result_reactome,
                             category = NA,
                             subcategory = NA) %>%
                  dplyr::select(category, subcategory, dplyr::everything())

                module_result_kegg <-
                  module_result_kegg %>%
                  dplyr::select(category, subcategory, dplyr::everything())
              }
            }
          }

          temp_data <-
            rbind(module_result_kegg, module_result_reactome)

          if (is.null(temp_data)) {
            temp_data <-
              module_result_go
          } else{
            if (!is.null(module_result_go)) {
              temp_data <-
                temp_data %>%
                dplyr::full_join(module_result_go, by = intersect(colnames(.), colnames(module_result_go)))
            }
          }

          temp_data <-
            temp_data %>%
            dplyr::filter(p.adjust < p.adjust.cutoff &
                            Count > count.cutoff) %>%
            dplyr::select(-Description) %>%
            dplyr::rename(Description = module_annotation)
        }
      } else if (query_type == "metabolite") {
        if (length(object@merged_pathway_hmdb) == 0 &
            length(object@merged_pathway_metkegg) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else {
          module_result_hmdb <-
            tryCatch(
              object@merged_pathway_hmdb$module_result%>%
                dplyr::mutate(class = "HMDB"),
              error = function(e) {
                NULL
              }
            )
          module_result_metkegg <-
            tryCatch(
              object@merged_pathway_metkegg$module_result %>%
                dplyr::mutate(class = "KEGG"),
              error = function(e) {
                NULL
              }
            )

          temp_data <-
            rbind(module_result_hmdb, module_result_metkegg)

          temp_data <-
            temp_data %>%
            dplyr::filter(p.adjust < p.adjust.cutoff &
                            Count > count.cutoff) %>%
            dplyr::select(-Description) %>%
            dplyr::rename(Description = module_annotation)
        }
      }
    }

    if (level == "functional_module") {
      if (length(object@merged_module) == 0 &
          all(names(object@process_info) != "merge_modules")) {
        stop("Please use the merge_modules() function to process first")
      }

      if (length(object@merged_module) == 0 &
          any(names(object@process_info) == "merge_modules")) {
        stop('No enriched functional modules')
      }

      functional_module_result <-
        object@merged_module$functional_module_result

      if (any(colnames(functional_module_result) == "module_content")) {
        temp_data <-
          functional_module_result %>%
          dplyr::filter(p.adjust < p.adjust.cutoff &
                          Count > count.cutoff) %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = module_annotation) %>%
          dplyr::mutate(
            class = purrr::map(module_content, \(x) {
              modules <- stringr::str_split(x, ";")[[1]]
              dbs <- unique(dplyr::case_when(
                stringr::str_detect(modules, "^go_Module") ~ "GO",
                stringr::str_detect(modules, "^kegg_Module") ~ "KEGG",
                stringr::str_detect(modules, "^reactome_Module") ~ "Reactome",
                stringr::str_detect(modules, "^hmdb_Module") ~ "HMDB"
              ))
              paste(dbs, collapse = "/")
            })
          )
      } else {
        ## For biotext emebdding result, the ndoe is the module_content
        temp_data <-
          functional_module_result %>%
          dplyr::filter(p.adjust < p.adjust.cutoff &
                          Count > count.cutoff) %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = module_annotation) %>%
          dplyr::mutate(
            class = purrr::map(node, \(x) {
              modules <- stringr::str_split(x, ";")[[1]]
              dbs <- unique(dplyr::case_when(
                stringr::str_detect(modules, "GO") ~ "GO",
                stringr::str_detect(modules, "hsa") ~ "KEGG",
                stringr::str_detect(modules, "R-HSA") ~ "Reactome",
                stringr::str_detect(modules, "SMP") ~ "HMDB"
              ))
              paste(dbs, collapse = "/")
            })
          )
      }

    }

    database2 <-
      database

    database2[database2 == "go"] <- "GO"
    database2[database2 == "kegg"] <- "KEGG"
    database2[database2 == "reactome"] <- "Reactome"
    database2[database2 == "hmdb"] <- "HMDB"

    ### Select the representative pathway (min adjusted p value) in the module to plot the barplot
    if (level == "functional_module") {
      temp_data <-
        temp_data %>%
        dplyr::filter(sapply(stringr::str_split(class, "/"), function(x) {
          all(x %in% database2)
        })) %>%
        dplyr::mutate(
          class = purrr::map(class, \(x) {
            dbs <- stringr::str_split(x, "/")[[1]]
            dbs[1]
          })
        )

    } else {
      temp_data <-
        temp_data %>%
        dplyr::filter(class %in% database2)
    }


    if (nrow(temp_data) == 0) {
      warning("No suitable pathways or modules are available")
      return(
        temp_data %>%
          dplyr::arrange(p.adjust) %>%
          head(top_n) %>%
          dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
          dplyr::arrange(log.p) %>%
          dplyr::mutate(Description = factor(Description, levels = Description)) %>%
          ggplot(aes(log.p, Description)) +
          scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
          geom_segment(
            aes(
              x = 0,
              y = Description,
              xend = log.p,
              yend = Description,
              color = class
            )
          ) +
          geom_point(
            aes(size = Count, fill = class),
            shape = 21,
            alpha = 1
          ) +
          scale_size_continuous(range = c(3, 7)) +
          scale_fill_manual(values = database_color) +
          scale_color_manual(values = database_color) +
          theme_bw() +
          labs(y = "", x = "-log10(FDR adjusted P-values)") +
          geom_vline(xintercept = 0) +
          theme(panel.grid.minor = element_blank())
      )
    }

    if (analysis_type == "enrich_pathway") {
      if (query_type == "gene") {
        temp_data <-
          temp_data %>%
          dplyr::mutate(qscore = -log(p.adjust, 10)) %>%
          dplyr::arrange(.data[[x_axis_name]]) %>%
          tail(top_n)

        temp_data[[x_axis_name]] <- as.numeric(temp_data[[x_axis_name]])
      } else if (query_type == "metabolite") {
        temp_data <-
          temp_data %>%
          dplyr::mutate(qscore = -log(p.adjust, 10)) %>%
          dplyr::arrange(dplyr::desc(qscore)) %>%
          head(top_n)
      }
    } else{
      temp_data <-
        temp_data %>%
        dplyr::arrange(dplyr::desc(abs(NES))) %>%
        head(top_n)
      temp_data$NES <- as.numeric(temp_data$NES)
    }

    plot <-
      plot4pathway_enrichment(
        temp_data = temp_data,
        line_type = line_type,
        level = level,
        analysis_type = analysis_type,
        query_type = query_type,
        x_axis_name = x_axis_name,
        y_label_width = y_label_width,
        database_color = database_color
      )

    plot
  }


plot4pathway_enrichment <-
  function(temp_data = temp_data,
           line_type = c("straight", "meteor"),
           level = level,
           analysis_type = c("enrich_pathway", "do_gsea"),
           query_type = NULL,
           x_axis_name = c("NES", "qscore", "RichFactor", "FoldEnrichment"),
           y_label_width = 50,
           database_color =
             c(GO = "#eeca40",
               KEGG = "#fd7541",
               Reactome = "#23b9c7",
               HMDB = "#7998ad")
             ) {
    line_type <-
      match.arg(line_type)
    analysis_type <-
      match.arg(analysis_type)
    x_axis_name <-
      match.arg(x_axis_name)

    temp_data$class <- as.character(temp_data$class)

    ###line_type is straight

    if (line_type == "straight") {
      if (analysis_type == "enrich_pathway") {
        if (query_type == "gene") {
          plot <-
            ggplot(data = temp_data,
                   aes(x = .data[[x_axis_name]],
                       y = reorder(Description, .data[[x_axis_name]], decreasing = FALSE))) +
            scale_y_discrete(
              labels = function(x)
                stringr::str_wrap(x, width = y_label_width)
            ) +
            #scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
            coord_cartesian(xlim = c(0, max(as.numeric(temp_data[[x_axis_name]]))*1.1),
                            ylim = c(0.5, nrow(temp_data) + 0.5),
                            expand = FALSE) +
            geom_segment(
              aes(
                x = 0,
                y = Description,
                xend = .data[[x_axis_name]],
                yend = Description,
                color = class
              ),
              show.legend = FALSE
            ) +
            geom_point(aes(size = Count,
                           fill = class),
                       shape = 21) +
            #alpha = 1) +
            scale_size_continuous(range = c(3, 7)) +
            scale_fill_manual(values = database_color) +
            scale_color_manual(values = database_color) +
            theme_bw() +
            labs(
              title = paste0("Level: ", level),
              y = "",
              x = x_axis_name,
              size = "Gene number",
              fill = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(fill = guide_legend(override.aes = list(size = 5)))
        } else if (query_type == "metabolite") {
          plot <-
            ggplot(data = temp_data,
                   aes(x = temp_data$qscore,
                       y = reorder(temp_data$pathway_name, temp_data$qscore, decreasing = FALSE))) +
            scale_y_discrete(
              labels = function(x)
                stringr::str_wrap(x, width = y_label_width)
            ) +
            #scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
            coord_cartesian(xlim = c(0, max(as.numeric(temp_data$qscore))*1.1),
                            ylim = c(0.5, nrow(temp_data) + 0.5),
                            expand = FALSE) +
            geom_segment(
              aes(
                x = 0,
                y = temp_data$pathway_name,
                xend = temp_data$qscore,
                yend = temp_data$pathway_name,
                color = temp_data$class
              ),
              show.legend = FALSE
            ) +
            geom_point(aes(size = temp_data$Count,
                           fill = temp_data$class),
                       shape = 21) +
            #alpha = 1) +
            scale_size_continuous(range = c(3, 7)) +
            scale_fill_manual(values = database_color) +
            scale_color_manual(values = database_color) +
            theme_bw() +
            labs(
              title = paste0("Level: ", level),
              y = "",
              x = "qscore",
              size = "Gene number",
              fill = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(fill = guide_legend(override.aes = list(size = 5)))
        }
      }

      if (analysis_type == "do_gsea") {
        plot <-
          temp_data %>%
          dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
          dplyr::arrange(NES) %>%
          dplyr::mutate(Description = factor(Description, levels = Description)) %>%
          ggplot(aes(NES, Description)) +
          scale_y_discrete(
            labels = function(x)
              stringr::str_wrap(x, width = y_label_width)
          ) +
          scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
          geom_segment(
            aes(
              x = 0,
              y = Description,
              xend = NES,
              yend = Description,
              color = class
            ),
            show.legend = FALSE
          ) +
          geom_point(aes(size = Count, fill = class),
                     shape = 21) +
                     #alpha = 1) +
          scale_size_continuous(range = c(3, 7),
                                breaks = round(seq(min(temp_data$Count), max(temp_data$Count), length.out = 4))) +
          scale_fill_manual(values = database_color) +
          scale_color_manual(values = database_color) +
          theme_bw() +
          labs(
            title = paste0("Level: ", level),
            y = "",
            x = "NES",
            size = "Gene number",
            fill = "Database"
          ) +
          geom_vline(xintercept = 0) +
          theme(panel.grid.minor = element_blank(),
                axis.ticks.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
          guides(fill = guide_legend(override.aes = list(size = 5)))
      }
    }



    if (line_type == "meteor") {

      if (analysis_type == "enrich_pathway") {
        if (query_type == "gene") {
          plot <-
            ggplot(data = temp_data,
                   aes(x = .data[[x_axis_name]],
                       y = reorder(Description, .data[[x_axis_name]], decreasing = FALSE))) +
            scale_y_discrete(
              labels = function(x)
                stringr::str_wrap(x, width = y_label_width)
            ) +
            #scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
            coord_cartesian(xlim = c(0, max(as.numeric(temp_data[[x_axis_name]]))*1.1),
                            ylim = c(0.5, nrow(temp_data) + 0.5),
                            expand = FALSE) +
            # geom_segment(aes(
            #   x = 0,
            #   y = Description,
            #   xend = log.p,
            #   yend = Description,
            #   color = class
            # ),
            # show.legend = FALSE) +
            ggforce::geom_link(
              aes(
                x = 0,
                y = Description,
                xend = .data[[x_axis_name]],
                yend = Description,
                alpha = after_stat(index),
                size = after_stat(index),
                color = class
              ),
              show.legend = FALSE
            ) +
            geom_point(
              aes(size = Count, color = class),
              shape = 21,
              size = 6,
              fill = "white",
              alpha = 1
            ) +
            geom_text(
              aes(
                x = .data[[x_axis_name]],
                y = Description,
                label = paste("Gene number:", Count)
              ),
              size = 2.5,
              color = "white",
              hjust = 1.2,
              nudge_x = 0.05
            ) +
            scale_size_continuous(range = c(3, 7)) +
            scale_fill_manual(values = database_color) +
            scale_color_manual(values = database_color) +
            theme_bw() +
            labs(
              title = paste0("Level: ", level),
              y = "",
              x = x_axis_name,
              size = "Gene number",
              fill = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(fill = guide_legend(override.aes = list(size = 5)))

        } else if (query_type == "metabolite") {
          plot <-
            ggplot(data = temp_data,
                   aes(x = temp_data$qscore,
                       y = reorder(temp_data$pathway_name, temp_data$qscore, decreasing = FALSE))) +
            scale_y_discrete(
              labels = function(x)
                stringr::str_wrap(x, width = y_label_width)
            ) +
            #scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
            coord_cartesian(xlim = c(0, max(as.numeric(temp_data$qscore))*1.1),
                            ylim = c(0.5, nrow(temp_data) + 0.5),
                            expand = FALSE) +
            # geom_segment(aes(
            #   x = 0,
            #   y = Description,
            #   xend = log.p,
            #   yend = Description,
            #   color = class
            # ),
            # show.legend = FALSE) +
            ggforce::geom_link(
              aes(
                x = 0,
                y = temp_data$pathway_name,
                xend = temp_data$qscore,
                yend = temp_data$pathway_name,
                alpha = after_stat(index),
                size = after_stat(index),
                color = class
              ),
              show.legend = FALSE
            ) +
            geom_point(
              aes(size = temp_data$Count, color = temp_data$class),
              shape = 21,
              size = 6,
              fill = "white",
              alpha = 1
            ) +
            geom_text(
              aes(
                x = temp_data$qscore,
                y = temp_data$pathway_name,
                label = paste("Gene number:", temp_data$Count)
              ),
              size = 2.5,
              color = "white",
              hjust = 1.2,
              nudge_x = 0.05
            ) +
            scale_size_continuous(range = c(3, 7)) +
            scale_fill_manual(values = database_color) +
            scale_color_manual(values = database_color) +
            theme_bw() +
            labs(
              title = paste0("Level: ", level),
              y = "",
              x = "qscore",
              size = "Gene number",
              fill = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(fill = guide_legend(override.aes = list(size = 5)))
        }
      }

      if (analysis_type == "do_gsea") {
        plot <-
          temp_data %>%
          dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
          dplyr::arrange(NES) %>%
          dplyr::mutate(Description = factor(Description, levels = Description)) %>%
          ggplot(aes(NES, Description)) +
          scale_y_discrete(
            labels = function(x)
              stringr::str_wrap(x, width = y_label_width)
          ) +
          scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
          ggforce::geom_link(
            aes(
              x = 0,
              y = Description,
              xend = NES,
              yend = Description,
              alpha = after_stat(index),
              size = after_stat(index),
              color = class
            ),
            show.legend = FALSE
          ) +
          geom_point(
            aes(size = Count, color = class),
            shape = 21,
            size = 6,
            fill = "white",
            alpha = 1
          ) +
          geom_text(
            aes(
              x = NES * 0.65, # Position text halfway along the link
              y = Description,
              label = paste("Gene number:", Count)
            ),
            size = 2.5,
            color = "black"
          ) +
          scale_size_continuous(range = c(3, 7)) +
          scale_fill_manual(values = database_color) +
          scale_color_manual(values = database_color) +
          theme_bw() +
          labs(
            title = paste0("Level: ", level),
            y = "",
            x = "NES",
            size = "Gene number",
            fill = "Database"
          ) +
          geom_vline(xintercept = 0) +
          theme(panel.grid.minor = element_blank(),
                axis.ticks.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
          guides(fill = guide_legend(override.aes = list(size = 5)))
      }

    }

    plot

  }
