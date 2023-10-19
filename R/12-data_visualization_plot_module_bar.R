# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("demo_data.rda")
# load("result/modules")
# load("result/enriched_pathways")
# load("result/functional_module")
#
# object <-
#   functional_module
#
# plot_pathway_bar(
#   object = enriched_pathways,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "pathway"
# )
#
# plot_pathway_bar(
#   object = modules,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "pathway"
# )
#
# plot_pathway_bar(
#   object = modules,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "module"
# )
#
#
# plot_pathway_bar(
#   object = functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "pathway"
# )
#
# plot_pathway_bar(
#   object = functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "module"
# )
#
# plot_pathway_bar(
#   object = functional_module,
#   top_n = 10,
#   y_lable_width = 30,
#   level = "functional_module"
# )



#' Plot Pathway Bar Chart
#'
#' This function plots a bar chart representing enriched pathways, modules, or functional modules.
#' The bars are colored according to the database of origin.
#'
#' @param object An object containing the enrichment results and other relevant data.
#' @param top_n An integer specifying the top N pathways to display.
#' @param y_lable_width An integer specifying the width of the Y-axis labels.
#' @param level A character string specifying the level of analysis. One of "pathway", "module", or "functional_module".
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
#' # Assume `obj` is an object containing the relevant data
#' plot_pathway_bar(obj, top_n = 5, level = "pathway")
#' plot_pathway_bar(obj, top_n = 10, level = "module", p.adjust.cutoff = 0.05)
#' }
#'
#' @export

plot_pathway_bar <-
  function(object,
           top_n = 10,
           y_lable_width = 50,
           level = c("pathway",
                     "module",
                     "functional_module"),
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

    if (level == "pathway") {
      enrichment_go_result <-
        object@enrichment_go_result@result %>%
        dplyr::mutate(class = "GO")
      enrichment_kegg_result <-
        object@enrichment_kegg_result@result %>%
        dplyr::mutate(class = "KEGG")
      enrichment_reactome_result <-
        object@enrichment_reactome_result@result %>%
        dplyr::mutate(class = "Reactome")

      temp_data <-
        enrichment_go_result %>%
        dplyr::full_join(enrichment_kegg_result,
                         by = intersect(colnames(.),
                                        colnames(enrichment_kegg_result))) %>%
        dplyr::full_join(enrichment_reactome_result,
                         by = colnames(enrichment_reactome_result)) %>%
        dplyr::filter(p.adjust < p.adjust.cutoff &
                        Count > count.cutoff)
    }

    if (level == "module") {
      if (length(object@merged_pathway_go) == 0 &
          length(object@merged_pathway_go) == 0 &
          length(object@merged_pathway_go) == 0) {
        stop("Please use the merge_pathways() function to process first")
      } else{
        module_result_go <-
          object@merged_pathway_go$module_result %>%
          dplyr::mutate(class = "GO")
        module_result_kegg <-
          object@merged_pathway_kegg$module_result %>%
          dplyr::mutate(class = "KEGG")
        module_result_reactome <-
          object@merged_pathway_reactome$module_result %>%
          dplyr::mutate(class = "Reactome")

        temp_data <-
          module_result_go %>%
          dplyr::full_join(module_result_kegg,
                           by = intersect(colnames(.),
                                          colnames(module_result_kegg))) %>%
          dplyr::full_join(module_result_reactome,
                           by = colnames(module_result_reactome)) %>%
          dplyr::filter(p.adjust < p.adjust.cutoff &
                          Count > count.cutoff) %>%
          dplyr::select(-Description) %>%
          dplyr::rename(Description = module_annotation)

      }
    }


    if (level == "function_module") {
      if (length(object@merged_module) == 0) {
        stop("Please use the merge_modules() function to process first")
      } else{
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
    }

    database2 <-
      switch (database,
              go = "GO",
              kegg = "KEGG",
              reactome = "Reactome")

    plot <-
      temp_data %>%
      dplyr::filter(class %in% database2) %>%
      dplyr::arrange(p.adjust) %>%
      head(top_n) %>%
      dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
      dplyr::arrange(log.p) %>%
      dplyr::mutate(Description = factor(Description, levels = Description)) %>%
      ggplot(aes(log.p, Description)) +
      scale_y_discrete(
        labels = function(x)
          str_wrap(x, width = y_lable_width)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      geom_segment(aes(
        x = 0,
        y = Description,
        xend = log.p,
        yend = Description,
        color = class
      )) +
      geom_point(aes(size = Count,
                     fill = class),
                 shape = 21,
                 alpha = 1) +
      scale_size_continuous(range = c(3, 7)) +
      scale_fill_manual(values = database_color) +
      scale_color_manual(values = database_color) +
      theme_bw() +
      labs(y = "", x = "-log10(FDR adjusted P-values)") +
      geom_vline(xintercept = 0) +
      theme(panel.grid.minor = element_blank())

    plot
  }
