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
#   x_axis_name = "FoldEnrichment", #"qscore", "RichFactor", or "FoldEnrichment"
#   y_label_width = 30,
#   level = "module", # "pathway", "module", "functional_module"
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
#   object = object,
#   top_n = 10,
#   y_label_width = 30,
#   level = "pathway",
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
#   object = functional_module_annotation,
#   top_n = 10,
#   y_label_width = 30,
#   level = "functional_module",
#   # level = "module",
#   line_type = "straight",
#   llm_text = TRUE
# )
#
# Metabolite enrichment result
# plot_pathway_bar(object =  enriched_functional_module_met,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  # level = "pathway",
#                  level = "functional_module",
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 0,
#                  database = c("kegg", "hmdb"))

# plot_pathway_bar(object = llm_interpreted_enriched_functional_module,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  level = "pathway",
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 5,
#                  database = c("metkegg"))

# plot_pathway_bar(object = llm_interpreted_enriched_functional_module,
#                  top_n = 5,
#                  x_axis_name = "qscore",
#                  level = "functional_module",
#                  llm_text = TRUE,
#                  line_type = "straight",
#                  p.adjust.cutoff = 0.05,
#                  count.cutoff = 5,
#                  database = c("hmdb", "kegg"))

#' Plot enrichment results as a horizontal bar chart
#'
#' `plot_pathway_bar()` visualises the *top N* enriched items at any chosen
#' level—**pathway**, **module**, or **functional module**—as a horizontal bar
#' chart.  Bars are coloured by the originating database (GO, KEGG, Reactome,
#' HMDB) and scaled by gene/metabolite count.  The x-axis metric adapts to the
#' analysis type:
#'
#' * **Over-representation analysis (ORA)**
#'   *Genes*: choose one of `"qscore"` (−log<sub>10</sub> FDR), `"RichFactor"`,
#'   or `"FoldEnrichment"`;
#'   *Metabolites*: `"qscore"` only.
#' * **GSEA** results always use `"NES"`.
#'
#'
#' @param object  A **functional_module** object containing enrichment and merge
#'   information.
#' @param top_n  Number of top items to display (default `10`).
#' @param x_axis_name  Metric for the x-axis (see Description, default `qscore`).  Ignored for
#'   GSEA (always `"NES"`).
#' @param y_label_width  Maximum character width before y-axis labels are
#'   wrapped (default `50`).
#' @param level  One of `"pathway"`, `"module"`, or `"functional_module"`.
#' @param llm_text  Use GPT-generated labels (`llm_module_name`) instead of
#'   `module_annotation` at *functional_module* level?  (Logical)
#' @param line_type  `"straight"` or `"meteor"`.
#' @param p.adjust.cutoff  FDR threshold used to filter results (default `0.05`).
#' @param count.cutoff  Minimum gene/metabolite count to keep (default `5`).
#' @param database_color  Named vector mapping databases to fill colours.
#' @param database  Character vector of databases to include; values are
#'   case-insensitive versions of `"go"`, `"kegg"`, `"reactome"`, `"hmdb"`, `"metkegg"`.
#'
#' @return A **ggplot** object.
#'
#' @author Xiaotao Shen, \email{shenxt1990@outlook.com}
#' @author Yifei Ge, \email{yifeii.ge@outlook.com}
#'
#' @examples
#' \dontrun{
#' data("enriched_functional_module", package = "mapa")
#' ## Basic pathway-level bar chart
#' plot_pathway_bar(enriched_functional_module, level = "pathway", top_n = 15)
#'
#' ## GSEA results, meteor style
#' plot_pathway_bar(enriched_functional_module,
#'                  level      = "pathway",
#'                  line_type  = "meteor")
#' }
#'
#' @import ggplot2
#' @import ggforce
#' @import dplyr
#' @importFrom stringr str_wrap str_split str_detect
#' @importFrom purrr map map_chr
#' @importFrom methods is
#' @export

plot_pathway_bar <-
  function(object,
           top_n = 10,
           x_axis_name = "qscore",
           y_label_width = 50,
           level = c("pathway", "module", "functional_module"),
           llm_text = FALSE,
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
           database = c("go", "kegg", "reactome", "hmdb", "metkegg")) {

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

    # if (translation) {
    #   if (all(names(object@process_info) != "translate_language")) {
    #     stop("Please use the 'translate_language' function to translate first.")
    #   } else{
    #     if (length(object@enrichment_go_result) > 0) {
    #       object@enrichment_go_result@result <-
    #         object@enrichment_go_result@result |>
    #         dplyr::select(-Description) |>
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #
    #     if (length(object@enrichment_kegg_result) > 0) {
    #       object@enrichment_kegg_result@result <-
    #         object@enrichment_kegg_result@result |>
    #         dplyr::select(-Description) |>
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #     if (length(object@enrichment_reactome_result) > 0) {
    #       object@enrichment_reactome_result@result <-
    #         object@enrichment_reactome_result@result |>
    #         dplyr::select(-Description) |>
    #         dplyr::rename(Description = Description_trans)
    #     }
    #
    #
    #     if (length(object@merged_pathway_go) > 0) {
    #       object@merged_pathway_go$module_result <-
    #         object@merged_pathway_go$module_result |>
    #         dplyr::select(-c(Description, module_annotation)) |>
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_pathway_kegg) > 0) {
    #       object@merged_pathway_kegg$module_result <-
    #         object@merged_pathway_kegg$module_result |>
    #         dplyr::select(-c(Description, module_annotation)) |>
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_pathway_reactome) > 0) {
    #       object@merged_pathway_reactome$module_result <-
    #         object@merged_pathway_reactome$module_result |>
    #         dplyr::select(-c(Description, module_annotation)) |>
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #
    #     if (length(object@merged_module) > 0) {
    #       object@merged_module$functional_module_result <-
    #         object@merged_module$functional_module_result |>
    #         dplyr::select(-c(Description, module_annotation)) |>
    #         dplyr::rename(Description = Description_trans,
    #                       module_annotation = module_annotation_trans)
    #     }
    #   }
    # }

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    if (level == "pathway") {
      ## For gene enrichment analysis
      if (query_type == "gene") {
        enrichment_go_result <-
          tryCatch(
            object@enrichment_go_result@result |>
              dplyr::mutate(class = "GO"),
            error = function(e) {
              NULL
            }
          )
        enrichment_kegg_result <-
          tryCatch(
            object@enrichment_kegg_result@result |>
              dplyr::mutate(class = "KEGG"),
            error = function(e) {
              NULL
            }
          )
        enrichment_reactome_result <-
          tryCatch(
            object@enrichment_reactome_result@result |>
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
                ) |>
                dplyr::select(category, subcategory, dplyr::everything())

              enrichment_kegg_result <-
                enrichment_kegg_result |>
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
              enrichment_go_result |>
              dplyr::full_join(temp_data, by = intersect(colnames(.), colnames(temp_data)))

            if (analysis_type == "do_gsea") {
              temp_data$Count <-
                unlist(lapply(
                  stringr::str_split(temp_data$core_enrichment, "/"),
                  length
                ))
            }

            temp_data <-
              temp_data |>
              dplyr::filter(p_adjust < p.adjust.cutoff &
                              Count > count.cutoff)
          }
        }
      } else if (query_type == "metabolite") {
        enrichment_hmdb_result <-
          tryCatch(
            object@enrichment_hmdb_result@result |>
              dplyr::mutate(class = "HMDB"),
            error = function(e) {
              NULL
            }
          )

        enrichment_metkegg_result <-
          tryCatch(
            object@enrichment_metkegg_result@result |>
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

        colnames(temp_data)[10] <- "Count"

        temp_data <-
          temp_data |>
          dplyr::filter(p_adjust < p.adjust.cutoff &
                          Count > count.cutoff) |>
          dplyr::rename(Description = pathway_name)
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
              object@merged_pathway_go$module_result |>
                dplyr::mutate(class = "GO"),
              error = function(e) {
                NULL
              }
            )
          module_result_kegg <-
            tryCatch(
              object@merged_pathway_kegg$module_result |>
                dplyr::mutate(class = "KEGG"),
              error = function(e) {
                NULL
              }
            )

          module_result_reactome <-
            tryCatch(
              object@merged_pathway_reactome$module_result |>
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
                             subcategory = NA) |>
                  dplyr::select(category, subcategory, dplyr::everything())

                module_result_kegg <-
                  module_result_kegg |>
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
                temp_data |>
                dplyr::full_join(module_result_go, by = intersect(colnames(.), colnames(module_result_go)))
            }
          }

          temp_data <-
            temp_data |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
            dplyr::rename(Description = module_annotation)
        }
      } else if (query_type == "metabolite") {
        if (length(object@merged_pathway_hmdb) == 0 &
            length(object@merged_pathway_metkegg) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else {
          module_result_hmdb <-
            tryCatch(
              object@merged_pathway_hmdb$module_result|>
                dplyr::mutate(class = "HMDB"),
              error = function(e) {
                NULL
              }
            )
          module_result_metkegg <-
            tryCatch(
              object@merged_pathway_metkegg$module_result |>
                dplyr::mutate(class = "KEGG"),
              error = function(e) {
                NULL
              }
            )

          temp_data <-
            rbind(module_result_hmdb, module_result_metkegg)

          temp_data <-
            temp_data |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
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
        if (llm_text) {
          temp_data <-
            functional_module_result |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
            dplyr::rename(Description = llm_module_name) |>
            dplyr::mutate(
              class = purrr::map(module_content, \(x) {
                modules <- stringr::str_split(x, ";")[[1]]
                dbs <- unique(dplyr::case_when(
                  stringr::str_detect(modules, "^go_Module") ~ "GO",
                  stringr::str_detect(modules, "^kegg_Module") ~ "KEGG",
                  stringr::str_detect(modules, "^reactome_Module") ~ "Reactome",
                  stringr::str_detect(modules, "^hmdb_Module") ~ "HMDB",
                  stringr::str_detect(modules, "^metkegg_Module") ~ "KEGG"
                ))
                paste(dbs, collapse = "/")
              })
            )
        } else {
          temp_data <-
            functional_module_result |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
            dplyr::rename(Description = module_annotation) |>
            dplyr::mutate(
              class = purrr::map(module_content, \(x) {
                modules <- stringr::str_split(x, ";")[[1]]
                dbs <- unique(dplyr::case_when(
                  stringr::str_detect(modules, "^go_Module") ~ "GO",
                  stringr::str_detect(modules, "^kegg_Module") ~ "KEGG",
                  stringr::str_detect(modules, "^reactome_Module") ~ "Reactome",
                  stringr::str_detect(modules, "^hmdb_Module") ~ "HMDB",
                  stringr::str_detect(modules, "^metkegg_Module") ~ "KEGG"
                ))
                paste(dbs, collapse = "/")
              })
            )
        }

      } else {
        ## For biotext emebdding result, the node is the module_content
        if (llm_text) {
          temp_data <-
            functional_module_result |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
            dplyr::rename(Description = llm_module_name) |>
            dplyr::mutate(
              class = purrr::map(node, \(x) {
                modules <- stringr::str_split(x, ";")[[1]]
                dbs <- unique(dplyr::case_when(
                  stringr::str_detect(modules, "GO") ~ "GO",
                  stringr::str_detect(modules, "hsa") ~ "KEGG",
                  stringr::str_detect(modules, "R-HSA") ~ "Reactome",
                  stringr::str_detect(modules, "SMP") ~ "HMDB",
                  TRUE ~ "KEGG"
                ))
                paste(dbs, collapse = "/")
              })
            )
        } else {
          temp_data <-
            functional_module_result |>
            dplyr::filter(p_adjust < p.adjust.cutoff &
                            Count > count.cutoff) |>
            dplyr::select(-Description) |>
            dplyr::rename(Description = module_annotation) |>
            dplyr::mutate(
              class = purrr::map(node, \(x) {
                modules <- stringr::str_split(x, ";")[[1]]
                dbs <- unique(dplyr::case_when(
                  stringr::str_detect(modules, "GO") ~ "GO",
                  stringr::str_detect(modules, "hsa") ~ "KEGG",
                  stringr::str_detect(modules, "R-HSA") ~ "Reactome",
                  stringr::str_detect(modules, "SMP") ~ "HMDB",
                  TRUE ~ "KEGG"
                ))
                paste(dbs, collapse = "/")
              })
            )
        }
      }

    }

    database2 <-
      database

    database2[database2 == "go"] <- "GO"
    database2[database2 == "kegg"] <- "KEGG"
    database2[database2 == "reactome"] <- "Reactome"
    database2[database2 == "hmdb"] <- "HMDB"
    database2[database2 == "metkegg"] <- "KEGG"

    ### Select the database of the representative pathway (min adjusted p value) in the module to color the module
    if (level == "functional_module") {
      temp_data <-
        temp_data |>
        dplyr::filter(sapply(stringr::str_split(class, "/"), function(x) {
          all(x %in% database2)
        })) |>
        dplyr::mutate(
          class = purrr::map(class, \(x) {
            dbs <- stringr::str_split(x, "/")[[1]]
            dbs[1]
          })
        )

    } else {
      temp_data <-
        temp_data |>
        dplyr::filter(class %in% database2)
    }


    if (nrow(temp_data) == 0) {
      warning("No suitable pathways or modules are available")
      return(
        temp_data |>
          dplyr::arrange(p_adjust) |>
          head(top_n) |>
          dplyr::mutate(log.p = -log(p_adjust, 10)) |>
          dplyr::arrange(log.p) |>
          dplyr::mutate(Description = factor(Description, levels = Description)) |>
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
          temp_data |>
          dplyr::mutate(qscore = -log(p_adjust, 10)) |>
          dplyr::arrange(.data[[x_axis_name]]) |>
          tail(top_n)

        temp_data[[x_axis_name]] <- as.numeric(temp_data[[x_axis_name]])
      } else if (query_type == "metabolite") {
        temp_data <-
          temp_data |>
          dplyr::mutate(qscore = -log(p_adjust, 10)) |>
          dplyr::arrange(dplyr::desc(qscore)) |>
          head(top_n)
      }
    } else{
      temp_data <-
        temp_data |>
        dplyr::arrange(dplyr::desc(abs(NES))) |>
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
                   aes(x = qscore,
                       y = reorder(Description, qscore, decreasing = FALSE))) +
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
                y = Description,
                xend = qscore,
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
          temp_data |>
          dplyr::mutate(log.p = -log(p_adjust, 10)) |>
          dplyr::arrange(NES) |>
          dplyr::mutate(Description = factor(Description, levels = Description)) |>
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
              color = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(color = guide_legend(override.aes = list(size = 5)))

        } else if (query_type == "metabolite") {
          plot <-
            ggplot(data = temp_data,
                   aes(x = temp_data$qscore,
                       y = reorder(temp_data$Description, temp_data$qscore, decreasing = FALSE))) +
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
                y = temp_data$Description,
                xend = temp_data$qscore,
                yend = temp_data$Description,
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
                y = temp_data$Description,
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
              color = "Database"
            ) +
            #geom_vline(xintercept = 0) +
            theme(panel.grid.minor = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.title = element_text(hjust = 0.5)) +
            guides(color = guide_legend(override.aes = list(size = 5)))
        }
      }

      if (analysis_type == "do_gsea") {
        plot <-
          temp_data |>
          dplyr::mutate(log.p = -log(p_adjust, 10)) |>
          dplyr::arrange(NES) |>
          dplyr::mutate(Description = factor(Description, levels = Description)) |>
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
            color = "Database"
          ) +
          geom_vline(xintercept = 0) +
          theme(panel.grid.minor = element_blank(),
                axis.ticks.y = element_blank(),
                plot.title = element_text(hjust = 0.5)) +
          guides(color = guide_legend(override.aes = list(size = 5)))
      }

    }

    plot

  }
