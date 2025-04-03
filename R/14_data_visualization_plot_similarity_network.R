# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")
# library(tidyverse)
# load("enriched_functional_module.rda")
#
# object <-
#   enriched_functional_module
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   database = "go",
#   degree_cutoff = 10
# )
#
# library(showtext)
# showtext_auto(enable = TRUE)
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   database = "go",
#   degree_cutoff = 10,
#   translation = TRUE
# )
#
#
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "functional_module",
#   database = "go",
#   degree_cutoff = 1,
#   module_id = "Functional_module_53",
#   text_all = TRUE
# )
#
# GSEA result
# plot_similarity_network(
#   object = gsea_enriched_functional_module,
#   level = "functional_module",
#   database = "go",
#   degree_cutoff = 1,
#   module_id = "Functional_module_2"
# )
#
# plot_similarity_network(
#   object = gsea_enriched_functional_module,
#   level = "module",
#   degree_cutoff = 2,
#   database = "go",
# )
#
# metabolite
# plot_similarity_network(
#   object = enriched_functional_modules,
#   level = "module",
#   degree_cutoff = 0,
#   database = c("kegg")
# )


#' Plot Similarity Network
#'
#' This function plots a similarity network based on the given parameters.
#' It can operate on different levels and databases, and provides options for
#' text annotations.
#'
#' @param object An object of class "functional_module" containing
#' the necessary data for plotting.
#' @param level A character string indicating the level of analysis,
#' either "module" or "functional_module".
#' @param database A character string indicating the database to be used,
#' one of "go", "kegg", "reactome", or "hmdb".
#' @param degree_cutoff A numeric value indicating the degree cutoff for
#' filtering nodes in the network.
#' @param module_id Optional, a character or numeric vector of module IDs
#' to include in the plot.
#' @param text Logical, whether to include text annotations on the plot.
#' @param text_all Logical, whether to include text annotations for all nodes.
#' @param translation translation or not.
#' @return A ggplot object representing the similarity network, or an empty
#' ggplot if no modules pass the filtering criteria.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @author Yifei Ge <yifeii.ge@outlook.com>
#'
#' @examples
#' \dontrun{
#' # Assume `obj` is a prepared object containing relevant data
#' plot_similarity_network(obj, level = "module", database = "go")
#' plot_similarity_network(obj, level = "functional_module", degree_cutoff = 2)
#' }
#'
#' @export

plot_similarity_network <-
  function(object,
           level = c("module",
                     "functional_module"),
           database = c("go", "kegg", "reactome", "hmdb"),
           degree_cutoff = 0,
           module_id,
           text = TRUE,
           text_all = FALSE,
           translation = FALSE) {
    level <-
      match.arg(level)
    database <-
      match.arg(database)

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    if(translation){
      if(all(names(object@process_info) != "translate_language")){
        stop("Please use the 'translate_language' function to translate first.")
      }
    }

    # Determine analysis type
    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <- object@process_info$enrich_pathway@parameter$query_type
    } else {
      analysis_type <- "do_gsea"
      query_type <- "gene"
    }

    if (level == "module") {
      if (all(names(object@process_info) != "merge_pathways")) {
        stop("Please use the merge_pathways() function to process first")
      }
      if (query_type == "gene") {
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

      } else if (query_type == "metabolite") {
        ###HMDB
        if (database == "hmdb") {
          if (length(object@merged_pathway_hmdb) == 0) {
            warning("No enriched GO modules")
            return(ggplot() +
                     geom_blank())
          } else{
            graph_data <-
              object@merged_pathway_hmdb$graph_data
            result_with_module <-
              object@merged_pathway_hmdb$result_with_module

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
          if (length(object@merged_pathway_metkegg) == 0) {
            warning("No enriched KEGG modules")
            return(ggplot() +
                     geom_blank())
          } else{
            graph_data <-
              object@merged_pathway_metkegg$graph_data
            result_with_module <-
              object@merged_pathway_metkegg$result_with_module

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
        if (query_type == "metabolite") {
          graph_data <-
            graph_data %>%
            tidygraph::activate(what = "nodes") %>%
            dplyr::rename(module_content_number = module_content_number.y)
        }

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

    graph_data <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::filter(module_content_number > degree_cutoff) %>%
      dplyr::mutate(label = if (query_type == "gene") Description else pathway_name)

    if (igraph::gorder(graph_data) == 0) {
      warning("No functional modules have degree > ", degree_cutoff)
      return(ggplot() +
               geom_blank())
    }


    if (!missing(module_id)) {
      graph_data <-
        graph_data %>%
        tidygraph::activate(what = "nodes") %>%
        dplyr::filter(module %in% module_id)

      if (igraph::gorder(graph_data) == 0) {
        warning(module_id, " is not in the graph data")
        return(ggplot() +
                 geom_blank())
      }

    }

    ###plot to show the clusters of GO terms
    cluster_label_module <-
      tryCatch(
        expr = {
          df <- igraph::as_data_frame(graph_data, what = "vertices")

          if (analysis_type == "enrich_pathway") {
            df %>%
              dplyr::group_by(module) %>%
              dplyr::arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
              dplyr::slice_head(n = 1) %>%
              dplyr::pull(label)
          } else {
            df %>%
              dplyr::group_by(module) %>%
              dplyr::arrange(desc(abs(NES)), desc(Count), .by_group = TRUE) %>%
              dplyr::slice_head(n = 1) %>%
              dplyr::pull(label)
        }
          },
        error = function(e) {

        },
        warning = function(w) {
          # Handle warning here if needed
        }
      )

    if (is.null(cluster_label_module)) {
      cluster_label_module <- ""
    }

    cluster_label_all <-
      graph_data %>%
      tidygraph::activate(what = "nodes") %>%
      dplyr::pull(label)

    plot <-
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
        aes(fill = module,
            size = if(analysis_type == "enrich_pathway") -log(p.adjust, 10) else abs(NES)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      guides(fill = guide_legend(ncol = 1)) +
      ggraph::scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      labs(size = if(analysis_type == "enrich_pathway") "-log10(FDR adjusted P-values)" else "abs(NES)") +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    if (text) {
      if (text_all) {
        plot <-
          plot +
          ggraph::geom_node_text(aes(x = x,
                                     y = y,
                                     label = label),
                                 size = 3,
                                 repel = TRUE)
      } else{
        plot <-
          plot +
          ggraph::geom_node_text(aes(
            x = x,
            y = y,
            label = ifelse(label %in% cluster_label_module, label, NA)
          ),
          size = 3,
          repel = TRUE)
      }
    }
    plot
  }


####output some results
# dir.create(
#   file.path(path, "Similarity_plot"),
#   recursive = TRUE,
#   showWarnings = FALSE
# )
#


# ##matrix tow show the cluster GO terms
# result_with_module %>%
#   dplyr::group_by(ONTOLOGY) %>%
#   dplyr::summarise(n = n()) %>%
#   dplyr::mutate(n = n * 10 / max(n) + 2)
#output the correlation matrix

# if(database == "go"){
#   message("Output correlation matrix plot...")
#   for(ont in c('MF', "BP", "CC")) {
#     cat(ont, " ")
#     show_matrix_cluster(
#       result = result_with_module %>% dplyr::mutate(Direction = "UP"),
#       ont = ont,
#       measure = "Wang",
#       remove_words = remove_words,
#       margin = 15,
#       width = 14,
#       height = 8,
#       path = path,
#       top = 15
#     )
#   }
# }

###output the cluster annotation for each cluster
# if (database == "go") {
#   dir.create(file.path(path, "GO_module_graph"), showWarnings = FALSE)
#
#   unique(module_result$module) %>%
#     purrr::map(
#       .f = function(x) {
#         cat(x, " ")
#         number <- module_result %>%
#           dplyr::filter(module == x) %>%
#           pull(module_content_number) %>%
#           as.numeric()
#         if (number == 1) {
#           return(NULL)
#         }
#
#         temp_id <-
#           module_result %>%
#           dplyr::filter(module == x) %>%
#           dplyr::pull(node) %>%
#           stringr::str_split(";") %>%
#           `[[`(1) %>%
#           pRoloc::goIdToTerm(keepNA = FALSE) %>%
#           data.frame(id = ., class = "YES") %>%
#           tibble::rownames_to_column(var = "name")
#
#         temp_plot =
#           GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
#           igraph::igraph.from.graphNEL() %>%
#           tidygraph::as_tbl_graph() %>%
#           left_join(temp_id, by = "name") %>%
#           dplyr::mutate(class = dplyr::case_when(is.na(class) ~ "NO",
#                                           TRUE ~ class))
#
#         plot =
#           temp_plot %>%
#           ggraph(layout = 'kk',
#                  circular = FALSE) +
#           geom_edge_link(
#             color = "#3B4992FF",
#             alpha = 1,
#             arrow = grid::arrow(
#               angle = 10,
#               length = unit(0.2, "inches"),
#               type = "closed"
#             ),
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class),
#             shape = 21,
#             alpha = 1,
#             size = 6,
#             show.legend = FALSE
#           ) +
#           geom_node_text(aes(
#             x = x,
#             y = y,
#             label = ifelse(class == "YES", id, NA)
#           ),
#           size = 3,
#           repel = TRUE) +
#           scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "left",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#         # plot
#         ggsave(
#           plot,
#           filename = file.path(
#             path,
#             "GO_module_graph",
#             paste(x, "_GO graph.pdf", sep = "")
#           ),
#           width = 7,
#           height = 7
#         )
#       }
#     )
# }
