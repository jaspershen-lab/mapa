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
# functional_module@merged_pathway_go$module_result$module[1]
#
# plot <-
#   plot_module_info(
#     object = functional_module,
#     level = "module",
#     database = "go",
#     module_id = "go_Module_3"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(functional_module@merged_pathway_kegg$result_with_module$module)
#
# plot <-
#   plot_module_info(
#     object = functional_module,
#     level = "module",
#     database = "kegg",
#     module_id = "kegg_Module_15"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(functional_module@merged_pathway_reactome$result_with_module$module)
#
# plot <-
#   plot_module_info(
#     object = functional_module,
#     level = "module",
#     database = "reactome",
#     module_id = "reactome_Module_7"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(functional_module@merged_module$result_with_module$module)
#
# plot <-
#   plot_module_info(object = functional_module,
#                    level = "functional_module",
#                    module_id = "Functional_module_27")
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
#
# export_module_info_plot(object = object,
#                         path = "result")

#' Plot Module Information for Enrichment Analysis
#'
#' This function generates various visualizations including a network plot, bar plot,
#' and word cloud to provide detailed insights into the enriched functional or biological modules.
#'
#' @param object An object of a specific class containing enrichment results.
#' @param level A character string specifying the level of biological organization, either "module" or "functional_module".
#' @param database A character string specifying the source database, either "go", "kegg" or "reactome".
#' @param module_id A single identifier specifying the module of interest.
#'
#' @return A list containing three ggplot objects: network plot (`network`), bar plot (`barplot`), and word cloud (`wordcloud`).
#'
#' @export
#' @author Xiaotao Shen <shenxt1990@outlook.com>
#' @examples
#' \dontrun{
#' result <- plot_module_info(object = myObject,
#' level = "module",
#' database = "go",
#' module_id = "M123")
#' }

plot_module_info <-
  function(object,
           level = c("module",
                     "functional_module"),
           database = c("go", "kegg", "reactome"),
           module_id) {
    level <-
      match.arg(level)
    database <-
      match.arg(database)

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    if (missing(module_id)) {
      stop("module_id is required")
    } else{
      if (length(module_id) != 1) {
        stop("module_id' length should be one")
      }
    }

    if (level == "module") {
      ###GO
      if (database == "go") {
        if (length(object@merged_pathway_go) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else{
          graph_data <-
            object@merged_pathway_go$graph_data
          result_with_module <-
            object@merged_pathway_go$result_with_module

        }
      }

      ###KEGG
      if (database == "kegg") {
        if (length(object@merged_pathway_kegg) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else{
          graph_data <-
            object@merged_pathway_kegg$graph_data
          result_with_module <-
            object@merged_pathway_kegg$result_with_module
        }
      }

      ###Reactome
      if (database == "reactome") {
        if (length(object@merged_pathway_reactome) == 0) {
          stop("Please use the merge_pathways() function to process first")
        } else{
          graph_data <-
            object@merged_pathway_reactome$graph_data
          result_with_module <-
            object@merged_pathway_reactome$result_with_module
        }
      }
    }

    if (level == "functional_module") {
      if (length(object@merged_module) == 0) {
        stop("Please use the merge_modules() function to process first")
      } else{
        graph_data <-
          object@merged_module$graph_data
        result_with_module <-
          object@merged_module$result_with_module
      }
    }

    if (!missing(module_id)) {
      if (length(module_id) != 1) {
        stop("module_id should be one name")
      }
      graph_data <-
        graph_data %>%
        tidygraph::activate(what = "nodes") %>%
        dplyr::filter(module %in% module_id)
    }

    ###network
    plot1 <-
      graph_data %>%
      ggraph(layout = 'fr',
             circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = -log(p.adjust, 10),
            size = Count),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(3, 10)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "left",
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    if (level == "module") {
      plot1 <-
        plot1 +
        geom_node_text(aes(x = x,
                           y = y,
                           label = Description),
                       check_overlap = TRUE)
    }

    if (level == "functional_module") {
      plot1 <-
        plot1 +
        geom_node_text(aes(x = x,
                           y = y,
                           label = module_annotation),
                       check_overlap = TRUE)
    }



    ###barplot
    if (level == "module") {
      temp_data <-
        result_with_module %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(module == module_id) %>%
        dplyr::mutate(log.p = -log(as.numeric(p.adjust, 10))) %>%
        dplyr::arrange(log.p) %>%
        dplyr::mutate(Description = factor(Description, levels = Description))
    } else{
      temp_data <-
        result_with_module %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(module == module_id) %>%
        dplyr::mutate(log.p = -log(as.numeric(p.adjust, 10))) %>%
        dplyr::arrange(log.p) %>%
        dplyr::select(-Description) %>%
        dplyr::rename(Description = module_annotation) %>%
        dplyr::mutate(Description = factor(Description, levels = Description))
    }

    plot2 <-
      temp_data %>%
      ggplot(aes(log.p, Description)) +
      scale_y_discrete(
        labels = function(x)
          str_wrap(x, width = 40)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
      geom_segment(aes(
        x = 0,
        y = Description,
        xend = log.p,
        yend = Description
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


    ###wordcloud
    temp_data <-
      result_with_module %>%
      dplyr::filter(module == module_id) %>%
      dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
      dplyr::select(Description, p.adjust) %>%
      dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
      plyr::dlply(.variables = .(Description)) %>%
      purrr::map(function(x) {
        data.frame(word = stringr::str_split(x$Description, " ")[[1]],
                   p.adjust = x$p.adjust)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      plyr::dlply(.variables = .(word)) %>%
      purrr::map(function(x) {
        x$p.adjust <- sum(x$p.adjust)
        x %>%
          dplyr::distinct(word, .keep_all = TRUE)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::filter(!word %in% remove_words)

    plot3 <-
      temp_data %>%
      ggplot(aes(label = word, size = p.adjust)) +
      ggwordcloud::geom_text_wordcloud() +
      scale_radius(range = c(5, 15), limits = c(0, NA)) +
      theme_minimal()

    list(network = plot1,
         barplot = plot2,
         wordcloud = plot3)

  }

#' Export Module Information Plot
#'
#' This function exports plots for functional modules, saving them as PDF files in the specified directory.
#' Plots for modules from the GO, KEGG, Reactome, and custom functional modules are supported.
#'
#' @param object A functional_module object.
#' @param path A character string specifying the directory path to save the plots.
#'
#' @return Prints messages indicating the progress and saves plots in the specified directory.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#'
#' @examples
#' \dontrun{
#' obj <- functional_module$new(...)
#' export_module_info_plot(obj, path = "result")
#' }
#'

export_module_info_plot <-
  function(object,
           path = "result") {
    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(
      file.path(path, "go/module_plot"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(path, "kegg/module_plot"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(path, "reactome/module_plot"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    dir.create(
      file.path(path, "functional_modules/module_plot"),
      recursive = TRUE,
      showWarnings = FALSE
    )

    ##GO
    message('GO...')
    all_module_id <-
      object@merged_pathway_go$result_with_module$module

    all_module_id <-
      data.frame(all_module_id = all_module_id) %>%
      dplyr::count(all_module_id) %>%
      dplyr::filter(n > 1) %>%
      pull(all_module_id)

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = all_module_id, function(module_id) {
        message(module_id, " ")
        plot <-
          plot_module_info(
            object = object,
            level = "module",
            database = "go",
            module_id = module_id
          )

        plot <-
          plot[[1]] + plot[[2]] + plot[[3]] + patchwork::plot_layout(nrow = 1)

        suppressMessages(extrafont::loadfonts())
        ggsave(
          plot,
          filename = file.path(
            path,
            "go/module_plot",
            paste(module_id, "plot.pdf", sep = "_")
          ),
          width = 21,
          height = 7
        )
      })
    }





    ##KEGG
    message('KEGG...')
    all_module_id <-
      object@merged_pathway_kegg$result_with_module$module

    all_module_id <-
      data.frame(all_module_id = all_module_id) %>%
      dplyr::count(all_module_id) %>%
      dplyr::filter(n > 1) %>%
      pull(all_module_id)

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = all_module_id, function(module_id) {
        message(module_id, " ")
        plot <-
          plot_module_info(
            object = object,
            level = "module",
            database = "kegg",
            module_id = module_id
          )

        plot <-
          plot[[1]] + plot[[2]] + plot[[3]] + patchwork::plot_layout(nrow = 1)

        ggsave(
          plot,
          filename = file.path(
            path,
            "kegg/module_plot",
            paste(module_id, "plot.pdf", sep = "_")
          ),
          width = 21,
          height = 7
        )
      })
    }


    ##Reactome
    message('Reactome...')
    all_module_id <-
      object@merged_pathway_reactome$result_with_module$module

    all_module_id <-
      data.frame(all_module_id = all_module_id) %>%
      dplyr::count(all_module_id) %>%
      dplyr::filter(n > 1) %>%
      pull(all_module_id)

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = all_module_id, function(module_id) {
        message(module_id, " ")
        plot <-
          plot_module_info(
            object = object,
            level = "module",
            database = "reactome",
            module_id = module_id
          )

        plot <-
          plot[[1]] + plot[[2]] + plot[[3]] + patchwork::plot_layout(nrow = 1)

        ggsave(
          plot,
          filename = file.path(
            path,
            "reactome/module_plot",
            paste(module_id, "plot.pdf", sep = "_")
          ),
          width = 21,
          height = 7
        )
      })
    }

    ##Functional module
    message('Functional module...')
    all_module_id <-
      object@merged_module$result_with_module$module

    all_module_id <-
      data.frame(all_module_id = all_module_id) %>%
      dplyr::count(all_module_id) %>%
      dplyr::filter(n > 1) %>%
      pull(all_module_id)

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = all_module_id, function(module_id) {
        message(module_id, " ")
        plot <-
          plot_module_info(
            object = object,
            level = "functional_module",
            database = "reactome",
            module_id = module_id
          )

        plot <-
          plot[[1]] + plot[[2]] + plot[[3]] + patchwork::plot_layout(nrow = 1)

        ggsave(
          plot,
          filename = file.path(
            path,
            "functional_modules/module_plot",
            paste(module_id, "plot.pdf", sep = "_")
          ),
          width = 21,
          height = 7
        )
      })
    }

    message("Done")

  }
