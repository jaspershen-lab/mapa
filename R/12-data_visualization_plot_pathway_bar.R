# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# load("enriched_functional_module.rda")
#
# object <-
#   enriched_functional_module
#
# library(showtext)
# showtext_auto(enable = TRUE)
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "pathway",
#   translation = TRUE
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "pathway",
#   translation = TRUE
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "module"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "functional_module"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "functional_module",
#   line_type = "meteor"
# )

#' Plot Pathway Bar Chart
#'
#' This function plots a bar chart representing enriched pathways, modules, or functional modules.
#' The bars are colored according to the database of origin.
#'
#' @param object An object containing the enrichment results and other relevant data.
#' @param top_n An integer specifying the top N pathways to display.
#' @param y_lable_width An integer specifying the width of the Y-axis labels.
#' @param translation translation or not.
#' @param level A character string specifying the level of analysis.
#' One of "pathway", "module", or "functional_module".
#' @param line_type A character string specifying the type of line to use for the bar chart.
#' One of "straight" or "meteor".
#' @param p.adjust.cutoff A numeric value for the FDR adjusted P-value cutoff.
#' @param count.cutoff A numeric value for the count cutoff.
#' @param database_color A named vector containing the colors for different databases.
#' @param database A character vector specifying which databases to include in the plot.
#' Should be "go", "kegg", or/and "reactome".
#'
#' @return A ggplot object representing the bar chart.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#'
#' @examples
#' data("enriched_functional_module", package = "mapa")
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_lable_width = 30,
#'   level = "pathway"
#' )
#'
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_lable_width = 30,
#'   level = "pathway",
#'   line_type = "meteor"
#' )
#'
#' plot_pathway_bar(enriched_functional_module,
#'                  top_n = 10,
#'                  level = "module")
#'
#' plot_pathway_bar(
#'   enriched_functional_module,
#'   top_n = 10,
#'   level = "functional_module",
#'   p.adjust.cutoff = 0.05
#' )
#'@export

plot_pathway_bar <-
  function(object,
           top_n = 10,
           y_lable_width = 50,
           translation = FALSE,
           level = c("pathway",
                     "module",
                     "functional_module"),
           line_type = c("straight", "meteor"),
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database_color =
             c(
               GO = "#1F77B4FF",
               KEGG = "#FF7F0EFF",
               Reactome = "#2CA02CFF"
             ),
           database = c("go", "kegg", "reactome")) {
    level <-
      match.arg(level)
    line_type <-
      match.arg(line_type)

    if (translation) {
      if(all(names(object@process_info) != "translate_language")){
        stop("Please use the 'translate_language' function to translate first.")
      }else{
        object@enrichment_go_result@result <-
          object@enrichment_go_result@result %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = Description_trans)

        object@enrichment_kegg_result@result <-
          object@enrichment_kegg_result@result %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = Description_trans)

        object@enrichment_reactome_result@result <-
          object@enrichment_reactome_result@result %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = Description_trans)

        object@merged_pathway_go$module_result <-
          object@merged_pathway_go$module_result %>%
          dplyr::select(-c(Description, module_annotation)) %>%
          dplyr::rename(Description = Description_trans,
                        module_annotation = module_annotation_trans)

        object@merged_pathway_kegg$module_result <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::select(-c(Description, module_annotation)) %>%
          dplyr::rename(Description = Description_trans,
                        module_annotation = module_annotation_trans)

        object@merged_pathway_reactome$module_result <-
          object@merged_pathway_reactome$module_result %>%
          dplyr::select(-c(Description, module_annotation)) %>%
          dplyr::rename(Description = Description_trans,
                        module_annotation = module_annotation_trans)

        object@merged_module$functional_module_result <-
          object@merged_module$functional_module_result %>%
          dplyr::select(-c(Description, module_annotation)) %>%
          dplyr::rename(Description = Description_trans,
                        module_annotation = module_annotation_trans)
      }
    }

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    if (level == "pathway") {
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

      temp_data <-
        rbind(enrichment_reactome_result,
              enrichment_kegg_result)

      if (is.null(temp_data)) {
        temp_data <-
          enrichment_go_result
      } else{
        if (!is.null(enrichment_go_result)) {
          temp_data <-
            enrichment_go_result %>%
            dplyr::full_join(temp_data,
                             by = intersect(colnames(.),
                                            colnames(temp_data))) %>%
            dplyr::filter(p.adjust < p.adjust.cutoff &
                            Count > count.cutoff)
        }
      }
    }

    if (level == "module") {
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

        temp_data <-
          rbind(module_result_kegg,
                module_result_reactome)

        if (is.null(temp_data)) {
          temp_data <-
            module_result_go
        } else{
          if (!is.null(module_result_go)) {
            temp_data <-
              temp_data %>%
              dplyr::full_join(module_result_go,
                               by = intersect(colnames(.),
                                              colnames(module_result_go)))
          }
        }

        temp_data <-
          temp_data %>%
          dplyr::filter(p.adjust < p.adjust.cutoff &
                          Count > count.cutoff) %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = module_annotation)

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

      temp_data <-
        functional_module_result %>%
        dplyr::filter(p.adjust < p.adjust.cutoff &
                        Count > count.cutoff) %>%
        dplyr::select(-Description) %>%
        dplyr::rename(Description = module_annotation,
                      class = database)
    }

    database2 <-
      database

    database2[database2 == "go"] <- "GO"
    database2[database2 == "kegg"] <- "KEGG"
    database2[database2 == "reactome"] <- "Reactome"

    temp_data <-
      temp_data %>%
      dplyr::filter(class %in% database2)

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
            aes(size = Count,
                fill = class),
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

    temp_data <-
      temp_data %>%
      dplyr::arrange(p.adjust) %>%
      head(top_n)

    if (line_type == "straight") {
      plot <-
        temp_data %>%
        dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
        dplyr::arrange(log.p) %>%
        dplyr::mutate(Description = factor(Description, levels = Description)) %>%
        ggplot(aes(log.p, Description)) +
        scale_y_discrete(
          labels = function(x)
            stringr::str_wrap(x, width = y_lable_width)
        ) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
        geom_segment(
          aes(
            x = 0,
            y = Description,
            xend = log.p,
            yend = Description,
            color = class
          ),
          show.legend = FALSE
        ) +
        geom_point(aes(size = Count,
                       fill = class),
                   shape = 21,
                   alpha = 1) +
        scale_size_continuous(range = c(3, 7)) +
        scale_fill_manual(values = database_color) +
        scale_color_manual(values = database_color) +
        theme_bw() +
        labs(
          y = "",
          x = "-log10(FDR adjusted P-values)",
          size = "Gene number",
          fill = "Database"
        ) +
        geom_vline(xintercept = 0) +
        theme(panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size = 5)))
    } else{
      plot <-
        temp_data %>%
        dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
        dplyr::arrange(log.p) %>%
        dplyr::mutate(Description = factor(Description, levels = Description)) %>%
        ggplot(aes(log.p, Description)) +
        scale_y_discrete(
          labels = function(x)
            str_wrap(x, width = y_lable_width)
        ) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
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
            xend = log.p,
            yend = Description,
            alpha = stat(index),
            size = after_stat(index),
            color = class
          ),
          show.legend = FALSE
        ) +
        geom_point(
          aes(size = Count,
              color = class),
          shape = 21,
          size = 6,
          fill = "white",
          alpha = 1
        ) +
        geom_text(
          aes(
            x = log.p,
            y = Description,
            label = paste("Gene number:", Count)
          ),
          size = 2.5,
          color = "grey20",
          hjust = 1.2,
          nudge_x = 0.05
        ) +
        scale_size_continuous(range = c(3, 7)) +
        scale_fill_manual(values = database_color) +
        scale_color_manual(values = database_color) +
        theme_bw() +
        labs(
          y = "",
          x = "-log10(FDR adjusted P-values)",
          size = "Gene number",
          fill = "Database"
        ) +
        geom_vline(xintercept = 0) +
        theme(panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank()) +
        guides(fill = guide_legend(override.aes = list(size = 5)))
    }

    plot
  }
