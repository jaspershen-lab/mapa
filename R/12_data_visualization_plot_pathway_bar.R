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
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_label_width = 30,
#   level = "pathway",
#   translation = FALSE
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_label_width = 30,
#   level = "pathway",
#   translation = TRUE
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_label_width = 30,
#   level = "module"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
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
#
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
# object <- gsea_pathways
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
#   translation = FALSE,
#   line_type = "meteor"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 10,
#   y_label_width = 30,
#   level = "module"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
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



#' Plot Pathway Bar Chart
#'
#' This function plots a bar chart representing enriched pathways, modules, or functional modules.
#' The bars are colored according to the database of origin.
#'
#' @param object An object containing the enrichment results and other relevant data.
#' @param top_n An integer specifying the top N pathways to display.
#' @param x_axis_name A character vector specifying the variable for x-axis when plotting ORA reslut. One of "qscore" (i.e. -log10(FDR adjusted P-values)), "RichFactor", or "FoldEnrichment". For GSEA result, use the default value NULL (i.e. NES).
#' @param y_label_width An integer specifying the width of the Y-axis labels.
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
#' \dontrun{
#' data("enriched_functional_module", package = "mapa")
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_label_width = 30,
#'   level = "pathway"
#' )
#'
#' plot_pathway_bar(
#'   object = enriched_functional_module,
#'   top_n = 10,
#'   y_label_width = 30,
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
#' }
#'@export

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
               GO = "#1F77B4FF",
               KEGG = "#FF7F0EFF",
               Reactome = "#2CA02CFF"
             ),
           database = c("go", "kegg", "reactome")) {
    level <-
      match.arg(level)
    line_type <-
      match.arg(line_type)

    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      if (is.null(x_axis_name)){
        stop("Please provide a variable for x axis.")
      } else {
        x_axis_name <- match.arg(x_axis_name, choices = c("qscore", "RichFactor", "FoldEnrichment"))
      }
    } else {
      analysis_type <- "do_gsea"
      if (!is.null(x_axis_name)) {
        warning("x_axis_name is for ORA. For GSEA, this parameter should be NULL, and use NES as x axis.")
      }
    }

    if (translation) {
      if (all(names(object@process_info) != "translate_language")) {
        stop("Please use the 'translate_language' function to translate first.")
      } else{
        if (length(object@enrichment_go_result) > 0) {
          object@enrichment_go_result@result <-
            object@enrichment_go_result@result %>%
            dplyr::select(-Description) %>%
            dplyr::rename(Description = Description_trans)
        }


        if (length(object@enrichment_kegg_result) > 0) {
          object@enrichment_kegg_result@result <-
            object@enrichment_kegg_result@result %>%
            dplyr::select(-Description) %>%
            dplyr::rename(Description = Description_trans)
        }

        if (length(object@enrichment_reactome_result) > 0) {
          object@enrichment_reactome_result@result <-
            object@enrichment_reactome_result@result %>%
            dplyr::select(-Description) %>%
            dplyr::rename(Description = Description_trans)
        }


        if (length(object@merged_pathway_go) > 0) {
          object@merged_pathway_go$module_result <-
            object@merged_pathway_go$module_result %>%
            dplyr::select(-c(Description, module_annotation)) %>%
            dplyr::rename(Description = Description_trans,
                          module_annotation = module_annotation_trans)
        }

        if (length(object@merged_pathway_kegg) > 0) {
          object@merged_pathway_kegg$module_result <-
            object@merged_pathway_kegg$module_result %>%
            dplyr::select(-c(Description, module_annotation)) %>%
            dplyr::rename(Description = Description_trans,
                          module_annotation = module_annotation_trans)
        }

        if (length(object@merged_pathway_reactome) > 0) {
          object@merged_pathway_reactome$module_result <-
            object@merged_pathway_reactome$module_result %>%
            dplyr::select(-c(Description, module_annotation)) %>%
            dplyr::rename(Description = Description_trans,
                          module_annotation = module_annotation_trans)
        }

        if (length(object@merged_module) > 0) {
          object@merged_module$functional_module_result <-
            object@merged_module$functional_module_result %>%
            dplyr::select(-c(Description, module_annotation)) %>%
            dplyr::rename(Description = Description_trans,
                          module_annotation = module_annotation_trans)
        }
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
        dplyr::rename(Description = module_annotation) %>%
        dplyr::mutate(
          class = purrr::map(module_content, \(x) {
            modules <- stringr::str_split(x, ";")[[1]]
            dbs <- unique(dplyr::case_when(
              stringr::str_detect(modules, "^go_Module") ~ "GO",
              stringr::str_detect(modules, "^kegg_Module") ~ "KEGG",
              stringr::str_detect(modules, "^reactome_Module") ~ "Reactome"
            ))
            paste(sort(dbs), collapse = "/")
          })
      )
    }

    database2 <-
      database

    database2[database2 == "go"] <- "GO"
    database2[database2 == "kegg"] <- "KEGG"
    database2[database2 == "reactome"] <- "Reactome"

    if (level == "functional_module") {
      temp_data <-
        temp_data %>%
        dplyr::mutate(
          class = purrr::map_chr(class, \(x) {
            dbs <- stringr::str_split(x, "/")[[1]]
            dbs <- intersect(dbs, database2)
            if (length(dbs) == 0) return(NA_character_)
            if (length(dbs) > 1) {
              colors <- col2rgb(database_color[dbs])
              mixed_color <- rgb(
                mean(colors[1,]) / 255,
                mean(colors[2,]) / 255,
                mean(colors[3,]) / 255,
                alpha = 0.6
              )
              database_color[x] <- mixed_color
            }
            paste(sort(dbs), collapse = "/")
          })
        ) %>%
        dplyr::filter(!is.na(class))
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
      # temp_data <-
      #   temp_data %>%
      #   dplyr::arrange(p.adjust) %>%
      #   head(top_n)
      temp_data <-
        temp_data %>%
        dplyr::mutate(qscore = -log(p.adjust, 10)) %>%
        dplyr::arrange(.data[[x_axis_name]]) %>%
        tail(top_n)

    } else{
      temp_data <-
        temp_data %>%
        dplyr::arrange(dplyr::desc(abs(NES))) %>%
        head(top_n)

    }

    temp_data[[x_axis_name]] <- as.numeric(temp_data[[x_axis_name]])

    plot <-
      plot4pathway_enrichment(
        temp_data = temp_data,
        line_type = line_type,
        level = level,
        analysis_type = analysis_type,
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
           x_axis_name = c(NULL, "qscore", "RichFactor", "FoldEnrichment"),
           y_label_width = 50,
           database_color =
             c(
               GO = "#1F77B4FF",
               KEGG = "#FF7F0EFF",
               Reactome = "#2CA02CFF"
             )) {
    line_type <-
      match.arg(line_type)
    analysis_type <-
      match.arg(analysis_type)
    x_axis_name <-
      match.arg(x_axis_name)

    ###line_type is straight

    if (line_type == "straight") {
      if (analysis_type == "enrich_pathway") {

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
              x = NES,
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
