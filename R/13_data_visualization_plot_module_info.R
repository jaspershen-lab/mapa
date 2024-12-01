# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# load("enriched_functional_module.rda")
#
# object <-
#   enriched_functional_module
#
# enriched_functional_module@merged_pathway_go$module_result$module[1]
# library(showtext)
# showtext_auto(enable = TRUE)
# plot <-
#   plot_module_info(
#     object = enriched_functional_module,
#     level = "module",
#     database = "go",
#     module_id = "go_Module_4",
#     translation = FALSE
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(enriched_functional_module@merged_pathway_kegg$result_with_module$module)
#
# plot <-
#   plot_module_info(
#     object = enriched_functional_module,
#     level = "module",
#     database = "kegg",
#     module_id = "kegg_Module_15"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(enriched_functional_module@merged_pathway_reactome$result_with_module$module)
#
# plot <-
#   plot_module_info(
#     object = enriched_functional_module,
#     level = "module",
#     database = "reactome",
#     module_id = "reactome_Module_7"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# library(showtext)
# showtext_auto(enable = TRUE)
#
# table(enriched_functional_module@merged_module$result_with_module$module)
#
# plot <-
#   plot_module_info(object = enriched_functional_module,
#                    level = "functional_module",
#                    module_id = "Functional_module_47")
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
#
# export_module_info_plot(object = object,
#                         path = "result")
#
##GSEA
# plot <-
#   plot_module_info(
#     object = gsea_enriched_functional_module,
#     level = "module",
#     database = "go",
#     module_id = "go_Module_1"
#   )
#
# plot[[1]]
# plot[[2]]
# plot[[3]]
#
# table(enriched_functional_module@merged_module$result_with_module$module)
#
# plot <-
#   plot_module_info(object = gsea_enriched_functional_module,
#                    level = "functional_module",
#                    module_id = "Functional_module_2")
#
# plot[[1]]
# plot[[2]]
# plot[[3]]

#' Plot Module Information for Enrichment Analysis
#'
#' This function generates various visualizations including a network plot, bar plot,
#' and word cloud to provide detailed insights into the enriched functional or biological modules.
#'
#' @param object An object of a specific class containing enrichment results.
#' @param level A character string specifying the level of biological organization, either "module" or "functional_module".
#' @param database A character string specifying the source database, either "go", "kegg" or "reactome".
#' @param module_id A single identifier specifying the module of interest.
#' @param translation translation or not.
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
           module_id,
           translation = FALSE) {
    level <-
      match.arg(level)
    database <-
      match.arg(database)

    # Determine analysis type based on process_info
    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
    } else{
      analysis_type <- "do_gsea"
    }

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

    if(translation){
      if(all(names(object@process_info) != "translate_language")){
        stop("Please use the 'translate_language' function to translate first.")
      }
    }

    if (level == "module") {
      if (all(names(object@process_info) != "merge_pathways")) {
        stop("Please use the merge_pathways() function to process first")
      }
      ###GO
      if (database == "go") {
        if (length(object@merged_pathway_go) == 0) {
          warning("No enriched GO modules")
          return(ggplot() +
                   geom_blank())
        } else{
          graph_data <-
            object@merged_pathway_go$graph_data
          result_with_module <-
            object@merged_pathway_go$result_with_module

          if(translation){
            result_with_module <-
              result_with_module %>%
              dplyr::select(-Description) %>%
              dplyr::rename(Description = Description_trans)
            graph_data <-
              graph_data %>%
              tidygraph::activate(what = "nodes") %>%
              dplyr::select(-Description) %>%
              dplyr::left_join(result_with_module[,c("node", "Description")],
                               by = "node")
          }

        }
      }

      ###KEGG
      if (database == "kegg") {
        if (length(object@merged_pathway_kegg) == 0) {
          warning("No enriched KEGG modules")
          return(ggplot() +
                   geom_blank())
        } else{
          graph_data <-
            object@merged_pathway_kegg$graph_data
          result_with_module <-
            object@merged_pathway_kegg$result_with_module

          if(translation){
            result_with_module <-
              result_with_module %>%
              dplyr::select(-Description) %>%
              dplyr::rename(Description = Description_trans)
            graph_data <-
              graph_data %>%
              tidygraph::activate(what = "nodes") %>%
              dplyr::select(-Description) %>%
              dplyr::left_join(result_with_module[,c("node", "Description")],
                               by = "node")
          }

        }
      }

      ###Reactome
      if (database == "reactome") {
        if (length(object@merged_pathway_reactome) == 0) {
          warning("No enriched Reactome modules")
          return(ggplot() +
                   geom_blank())
        } else{
          graph_data <-
            object@merged_pathway_reactome$graph_data
          result_with_module <-
            object@merged_pathway_reactome$result_with_module

          if(translation){
            result_with_module <-
              result_with_module %>%
              dplyr::select(-Description) %>%
              dplyr::rename(Description = Description_trans)
            graph_data <-
              graph_data %>%
              tidygraph::activate(what = "nodes") %>%
              dplyr::select(-Description) %>%
              dplyr::left_join(result_with_module[,c("node", "Description")],
                               by = "node")
          }

        }
      }
    }

    if (level == "functional_module") {
      if (all(names(object@process_info) != "merge_modules")) {
        stop("Please use the merge_modules() function to process first")
      }

      if (length(object@merged_module) == 0) {
        warning("No enriched functional modules")
        return(ggplot() +
                 geom_blank())
      } else{
        graph_data <-
          object@merged_module$graph_data
        result_with_module <-
          object@merged_module$result_with_module

        if(translation){
          result_with_module <-
            result_with_module %>%
            dplyr::select(-module_annotation) %>%
            dplyr::rename(module_annotation = module_annotation_trans)

          graph_data <-
            graph_data %>%
            tidygraph::activate(what = "nodes") %>%
            dplyr::select(-module_annotation)

          module_annotation <-
            result_with_module$module_annotation[match(igraph::vertex_attr(graph_data)$node,
                                                       result_with_module$node)]

          graph_data <-
            graph_data %>%
            tidygraph::activate(what = "nodes") %>%
            dplyr::mutate(module_annotation = module_annotation)
        }
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

    if (igraph::gorder(graph_data) == 0) {
      warning(module_id, " is not in the graph data")
      return(ggplot() +
               geom_blank())
    }

    # Set title based on level and database
    if (level == "module") {
      database2 <-
        database
      database2[database2 == "go"] <- "GO"
      database2[database2 == "kegg"] <- "KEGG"
      database2[database2 == "reactome"] <- "Reactome"
      plot_title <- paste0("Terms from ", database2)
    } else {
      # For functional module, determine databases from module content
      dbs <- object@merged_module$result_with_module %>%
        dplyr::filter(module == module_id) %>%
        dplyr::pull(database) %>%
        unique()
      plot_title <- paste0("Modules from ", paste(sort(dbs), collapse = " and "))
    }

    ###network
    plot1 <-
      graph_data %>%
      ggraph::ggraph(layout = 'fr',
                     circular = FALSE) +
      ggraph::geom_edge_link(
        aes(width = sim),
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      ggraph::geom_node_point(
        aes(fill = if(analysis_type == "enrich_pathway") -log(p.adjust, 10) else NES,
            size = Count),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(3, 10)) +
      labs(fill = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "NES") +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "left",
        legend.background = element_rect(fill = "transparent", color = NA),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(plot_title)

    if (level == "module") {
      plot1 <-
        plot1 +
        ggraph::geom_node_text(aes(x = x,
                                   y = y,
                                   label = stringr::str_wrap(Description, width = 30)),
                               check_overlap = TRUE,
                               size = 3,
                               repel = TRUE)  # Use repel to avoid overlapping
    }

    if (level == "functional_module") {
      plot1 <-
        plot1 +
        ggraph::geom_node_text(aes(x = x,
                                   y = y,
                                   label = stringr::str_wrap(module_annotation, width = 30)),
                               check_overlap = TRUE,
                               size = 3,
                               repel = TRUE)  # Use repel to avoid overlapping
    }

    ###barplot
    if (level == "module") {
      temp_data <-
        result_with_module %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(module == module_id)

      if(analysis_type == "enrich_pathway") {
        temp_data <- temp_data %>%
          dplyr::mutate(log.p = -log(as.numeric(p.adjust, 10))) %>%
          dplyr::arrange(log.p)
      } else {
        temp_data <- temp_data %>%
          dplyr::mutate(NES = as.numeric(NES)) %>%
          dplyr::arrange(NES)
      }

      temp_data <- temp_data %>%
        dplyr::mutate(Description = factor(Description, levels = Description))

    } else{
      temp_data <-
        result_with_module %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(module == module_id)

      if(analysis_type == "enrich_pathway") {
        temp_data <- temp_data %>%
          dplyr::mutate(log.p = -log(as.numeric(p.adjust, 10))) %>%
          dplyr::arrange(log.p)
      } else {
        temp_data <- temp_data %>%
          dplyr::mutate(NES = as.numeric(NES)) %>%
          dplyr::arrange(NES)
      }

      temp_data <- temp_data %>%
        dplyr::select(-Description) %>%
        dplyr::rename(Description = module_annotation) %>%
        dplyr::mutate(Description = factor(Description, levels = Description))
    }

    plot2 <-
      temp_data %>%
      ggplot(aes(if(analysis_type == "enrich_pathway") log.p else NES, Description)) +
      scale_y_discrete(
        labels = function(x)
          stringr::str_wrap(x, width = 40)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      geom_segment(aes(
        x = 0,
        y = Description,
        xend = if(analysis_type == "enrich_pathway") log.p else NES,
        yend = Description
      )) +
      geom_point(
        aes(size = Count,
            fill = if(level == "functional_module") database else "black"),
        shape = 21,
        alpha = 1
      ) +
      scale_size_continuous(range = c(3, 7)) +
      scale_fill_manual(values = c(GO = "#1F77B4FF",
                                   KEGG = "#FF7F0EFF",
                                   Reactome = "#2CA02CFF")) +
      theme_bw() +
      labs(y = "",
           x = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "NES",
           fill = "Database") +
      geom_vline(xintercept = 0) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      guides(fill = guide_legend(override.aes = list(size = 5))) +
      ggtitle(plot_title)

    ###wordcloud
    temp_data <-
      result_with_module %>%
      dplyr::filter(module == module_id) %>%
      dplyr::mutate(value = if(analysis_type == "enrich_pathway")
                            -log(as.numeric(p.adjust, 10))
                           else abs(as.numeric(NES))) %>%
      dplyr::select(Description, value) %>%
      dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
      plyr::dlply(.variables = .(Description)) %>%
      purrr::map(function(x) {
        data.frame(word = stringr::str_split(x$Description, " ")[[1]],
                   value = x$value)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      plyr::dlply(.variables = .(word)) %>%
      purrr::map(function(x) {
        x$value <- sum(x$value)
        x %>%
          dplyr::distinct(word, .keep_all = TRUE)
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::filter(!word %in% remove_words)

    plot3 <-
      suppressWarnings(
        ggplot(data = temp_data,
               aes(label = word, size = value)) +
          ggwordcloud::geom_text_wordcloud() +
          scale_radius(range = c(5, 15), limits = c(0, NA)) +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5)) +
          ggtitle(plot_title)
      )

    list(network = plot1,
         barplot = plot2,
         wordcloud = plot3)
  }

#' Export Module Information Plots
#'
#' This function exports plots of module information for different databases (GO, KEGG, Reactome, etc.)
#' from a given 'functional_module' object. It supports multiple image formats.
#'
#' @param object An object of class 'functional_module' containing the module information.
#' @param path The file path where the plots will be saved.
#' Default is "result".
#' @param image_type The format of the output image, can be 'pdf', 'png', or 'jpeg'.
#' Default is 'pdf'.
#' @param width Width of the output image in inches. Default is 21.
#' @param height Height of the output image in inches. Default is 7.
#'
#' @return This function does not return a value but outputs image files to the specified path.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'my_module' is a valid functional_module object
#'   export_module_info_plot(object = my_module, path = "my_results", image_type = "png")
#' }
#'
#' @export
#' @importFrom ggplot2 ggsave
#' @importFrom purrr walk
#' @importFrom dplyr count filter
#' @importFrom stringr str_sort
#' @importFrom patchwork plot_layout

export_module_info_plot <-
  function(object,
           path = "result",
           image_type = c("pdf", "png", "jpeg"),
           width = 21,
           height = 7) {
    image_type <-
      match.arg(image_type)
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
    message(rep("-", 20))
    message('GO...')
    all_module_id <-
      object@merged_pathway_go$result_with_module$module

    if (is.null(all_module_id)) {
      all_module_id <- character()
    } else{
      all_module_id <-
        data.frame(all_module_id = all_module_id) %>%
        dplyr::count(all_module_id) %>%
        dplyr::filter(n > 1) %>%
        pull(all_module_id)
    }

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = stringr::str_sort(all_module_id, numeric = TRUE),
                  function(module_id) {
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
                        paste(module_id, paste0("plot.", image_type), sep = "_")
                      ),
                      width = width,
                      height = height
                    )
                  })
    }

    ##KEGG
    message(rep("-", 20))
    message('KEGG...')
    all_module_id <-
      object@merged_pathway_kegg$result_with_module$module

    if (is.null(all_module_id)) {
      all_module_id <- character()
    } else{
      all_module_id <-
        data.frame(all_module_id = all_module_id) %>%
        dplyr::count(all_module_id) %>%
        dplyr::filter(n > 1) %>%
        pull(all_module_id)
    }

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = stringr::str_sort(all_module_id, numeric = TRUE),
                  function(module_id) {
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
                        paste(module_id, paste0("plot.", image_type), sep = "_")
                      ),
                      width = width,
                      height = height
                    )
                  })
    }


    ##Reactome
    message(rep("-", 20))
    message('Reactome...')
    all_module_id <-
      object@merged_pathway_reactome$result_with_module$module

    if (is.null(all_module_id)) {
      all_module_id <- character()
    } else{
      all_module_id <-
        data.frame(all_module_id = all_module_id) %>%
        dplyr::count(all_module_id) %>%
        dplyr::filter(n > 1) %>%
        pull(all_module_id)
    }

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = stringr::str_sort(all_module_id, numeric = TRUE),
                  function(module_id) {
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
                        paste(module_id, paste0("plot.", image_type), sep = "_")
                      ),
                      width = weight,
                      height = height
                    )
                  })
    }

    ##Functional module
    message(rep("-", 20))
    message('Functional module...')
    all_module_id <-
      object@merged_module$result_with_module$module

    if (is.null(all_module_id)) {
      all_module_id <- character()
    } else{
      all_module_id <-
        data.frame(all_module_id = all_module_id) %>%
        dplyr::count(all_module_id) %>%
        dplyr::filter(n > 1) %>%
        pull(all_module_id)
    }

    if (length(all_module_id) > 0) {
      message("Output module plot...")
      purrr::walk(.x = stringr::str_sort(all_module_id, numeric = TRUE),
                  function(module_id) {
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
                        paste(module_id, paste0("plot.", image_type), sep = "_")
                      ),
                      width = width,
                      height = height
                    )
                  })
    }

    message("Done")

  }
