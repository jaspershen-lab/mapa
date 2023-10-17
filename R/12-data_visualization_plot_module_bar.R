# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("result/function_modules/intermediate_data/module_result")
#
# plot_module_bar(head(module_result, 5))

#' plot_module_bar Function
#'
#' This function generates a bar plot for the top_n modules based on their significance from the module_result data.
#' It visualizes the -log10 transformed FDR adjusted p-values of the modules alongside their count size.
#'
#' @param module_result A data frame containing module results.
#' @param top_n An integer specifying the number of top modules to be plotted. Default is 5.
#'
#' @return A ggplot object.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @export

plot_module_bar <-
  function(module_result,
           top_n = 5) {
    plot <-
      module_result %>%
      dplyr::arrange(p.adjust) %>%
      head(top_n) %>%
      dplyr::mutate(log.p = -log(p.adjust, 10)) %>%
      dplyr::mutate(Count = as.numeric(Count)) %>%
      dplyr::arrange(log.p) %>%
      dplyr::mutate(module_annotation = factor(module_annotation, levels = module_annotation)) %>%
      ggplot(aes(log.p, module_annotation)) +
      scale_y_discrete(
        labels = function(x)
          str_wrap(x, width = 50)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      geom_segment(aes(
        x = 0,
        y = module_annotation,
        xend = log.p,
        yend = module_annotation
      )) +
      geom_point(
        aes(size = Count),
        fill = "black",
        shape = 21,
        alpha = 1
      ) +
      scale_size_continuous(range = c(3, 7)) +
      theme_bw() +
      labs(y = "", x = "-log10(FDR adjusted P-values)") +
      geom_vline(xintercept = 0) +
      theme(panel.grid.minor = element_blank())

    plot
  }
