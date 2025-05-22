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
#   database = "reactome",
#   degree_cutoff = 0
# )

# library(showtext)
# showtext_auto(enable = TRUE)
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "functional_module",
#   # database = c("go", "kegg", "reactome"),
#   degree_cutoff = 5
# )
#
# plot_similarity_network(
#   object = llm_interpreted_enriched_functional_module,
#   level = "functional_module",
#   degree_cutoff = 0,
#   llm_text = TRUE,
#   text_all = FALSE
# )
#
# GSEA result
# plot_similarity_network(
#   object = functional_module_annotation,
#   llm_text = TRUE,
#   level = "functional_module",
#   module_id = "Functional_module_4"
#   # degree_cutoff = 5
# )
#
# plot_similarity_network(
#   object = gsea_enriched_functional_module,
#   level = "module",
#   degree_cutoff = 2,
#   database = "reactome",
# )
#
# metabolite
# plot_similarity_network(
#   object = merged_pathways,
#   level = "module",
#   degree_cutoff = 0,
#   database = c("metkegg")
#   # module_id = "metkegg_Module_1"
# )
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "functional_module",
#   degree_cutoff = 0
# )

#' Plot Similarity Network
#'
#' Creates a visualization of similarity networks between biological pathways or functional modules.
#' The function supports various biological databases and offers flexible visualization options.
#'
#' @param object A "functional_module" class object containing network data and analysis results.
#' @param level Character string specifying the analysis level: "module" (database-specific modules) or
#'   "functional_module" (merged modules across databases). For results from get_bioembedsim() or
#'   merge_pathways_bioembedsim(), use "functional_module".
#' @param database Character string specifying the database to visualize: "go", "kegg", "reactome", "hmdb", or "metkegg".
#'   Only required when level = "module".
#' @param degree_cutoff Numeric value for filtering nodes; only nodes with module_content_number > degree_cutoff
#'   are displayed.
#' @param module_id Optional character or numeric vector of specific module IDs to include in the plot.
#' @param text Logical indicating whether to display module names (one per module), i.e., the name of
#'   the pathway with the minimum adjusted p value.
#' @param llm_text Logical indicating whether to display LLM-generated module names (one per module) when available.
#' @param text_all Logical indicating whether to display text annotations for all nodes (overrides the default
#'   behavior of showing only one annotation per module).
#'
#' @return A ggplot object representing the similarity network, or an empty ggplot if no modules meet
#'   the filtering criteria.
#'
#' @details
#' This function visualizes similarity networks at either the pathway module level (database-specific) or
#' functional module level (merged across databases). Node size represents statistical significance
#' (-log10 of adjusted p-value for enrichment or absolute NES for GSEA), while edge width represents
#' similarity between nodes. Text labels can be customized to show pathway name or LLM-generated names.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Basic usage for gene enrichment with GO database
#' plot_similarity_network(enrichment_result,
#'                         level = "module",
#'                         database = "go")
#'
#' # Visualize functional modules with degree filter
#' plot_similarity_network(enrichment_result,
#'                         level = "functional_module",
#'                         degree_cutoff = 2)
#'
#' # Use LLM-generated module names and show all text labels
#' plot_similarity_network(enrichment_result,
#'                         level = "functional_module",
#'                         text = TRUE,
#'                         llm_text = TRUE,
#'                         text_all = TRUE)
#'
#' # Show only specific modules
#' plot_similarity_network(enrichment_result,
#'                         level = "functional_module",
#'                         module_id = c("FM1", "FM2", "FM3"))
#' }
#'
#' @importFrom dplyr filter mutate select rename pull group_by arrange slice_head left_join summarise
#' @importFrom tidygraph activate as_tibble
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph create_layout
#' @importFrom igraph as_data_frame vertex_attr gorder
#' @importFrom ggplot2 ggplot geom_blank element_rect labs scale_size_continuous theme guides guide_legend
#' @importFrom stringr str_wrap
#'
#' @export

plot_similarity_network <-
  function(object,
           level = c("module",
                     "functional_module"),
           database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
           degree_cutoff = 0,
           module_id,
           llm_text = FALSE,
           text = TRUE,
           text_all = FALSE) {
    level <-
      match.arg(level)
    sim_method <- object@process_info$merge_pathways@function_name
    if (sim_method == "get_bioembedsim()" & level == "module") {
      stop("level `module` can not be applied to results generated by get_bioembedsim() and merge_pathways_bioembedsim() since these functions directly generate functional modules from enriched pathways. Please use `functional_module` level.")
    }

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    # if(translation){
    #   if(all(names(object@process_info) != "translate_language")){
    #     stop("Please use the 'translate_language' function to translate first.")
    #   }
    # }

    # Determine analysis type
    if ("enrich_pathway" %in% names(object@process_info)) {
      analysis_type <- "enrich_pathway"
      query_type <- object@process_info$enrich_pathway@parameter$query_type
    } else {
      analysis_type <- "do_gsea"
      query_type <- "gene"
    }

    if (level == "module") {
      database <-
        match.arg(database)
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
              object@merged_pathway_go$graph_data|>
              tidygraph::activate(what = "nodes")
            result_with_module <-
              object@merged_pathway_go$result_with_module

            # if(translation){
            #   result_with_module <-
            #     result_with_module |>
            #     dplyr::select(-Description) |>
            #     dplyr::rename(Description = Description_trans)
            #   graph_data <-
            #     graph_data |>
            #     tidygraph::activate(what = "nodes") |>
            #     dplyr::select(-Description) |>
            #     dplyr::left_join(result_with_module[,c("node", "Description")],
            #                      by = "node")
            # }

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
              object@merged_pathway_kegg$graph_data|>
              tidygraph::activate(what = "nodes")
            result_with_module <-
              object@merged_pathway_kegg$result_with_module

            # if(translation){
            #   result_with_module <-
            #     result_with_module |>
            #     dplyr::select(-Description) |>
            #     dplyr::rename(Description = Description_trans)
            #   graph_data <-
            #     graph_data |>
            #     tidygraph::activate(what = "nodes") |>
            #     dplyr::select(-Description) |>
            #     dplyr::left_join(result_with_module[,c("node", "Description")],
            #                      by = "node")
            # }
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
              object@merged_pathway_reactome$graph_data|>
              tidygraph::activate(what = "nodes")
            result_with_module <-
              object@merged_pathway_reactome$result_with_module

            # if(translation){
            #   result_with_module <-
            #     result_with_module |>
            #     dplyr::select(-Description) |>
            #     dplyr::rename(Description = Description_trans)
            #   graph_data <-
            #     graph_data |>
            #     tidygraph::activate(what = "nodes") |>
            #     dplyr::select(-Description) |>
            #     dplyr::left_join(result_with_module[,c("node", "Description")],
            #                      by = "node")
            # }
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
              object@merged_pathway_hmdb$graph_data |>
              tidygraph::activate(what = "nodes")

            result_with_module <-
              object@merged_pathway_hmdb$result_with_module

            # if(translation){
            #   result_with_module <-
            #     result_with_module |>
            #     dplyr::select(-Description) |>
            #     dplyr::rename(Description = Description_trans)
            #   graph_data <-
            #     graph_data |>
            #     tidygraph::activate(what = "nodes") |>
            #     dplyr::select(-Description) |>
            #     dplyr::left_join(result_with_module[,c("node", "Description")],
            #                      by = "node")
            # }

          }
        }

        ###KEGG
        if (database == "metkegg") {
          if (length(object@merged_pathway_metkegg) == 0) {
            warning("No enriched KEGG modules")
            return(ggplot() +
                     geom_blank())
          } else{
            graph_data <-
              object@merged_pathway_metkegg$graph_data |>
              tidygraph::activate(what = "nodes")

            result_with_module <-
              object@merged_pathway_metkegg$result_with_module

            # if(translation){
            #   result_with_module <-
            #     result_with_module |>
            #     dplyr::select(-Description) |>
            #     dplyr::rename(Description = Description_trans)
            #   graph_data <-
            #     graph_data |>
            #     tidygraph::activate(what = "nodes") |>
            #     dplyr::select(-Description) |>
            #     dplyr::left_join(result_with_module[,c("node", "Description")],
            #                      by = "node")
            # }
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
          object@merged_module$graph_data|>
          tidygraph::activate(what = "nodes")

        if (query_type == "metabolite") {
          col_name_node <-
            graph_data |>
            tidygraph::activate(what = "nodes") |>
            as_tibble() |>
            colnames()

          if ("Count" %in% col_name_node) {
            graph_data <-
              graph_data |>
              tidygraph::activate(what = "nodes") |>
              dplyr::select(-Count) |>
              dplyr::mutate(Count = as.numeric(mapped_number))
          } else {
            graph_data <-
              graph_data |>
              tidygraph::activate(what = "nodes") |>
              dplyr::mutate(Count = as.numeric(mapped_number))
          }

          if ("module_content_number.y" %in% col_name_node) {
            graph_data <-
              graph_data |>
              tidygraph::activate(what = "nodes") |>
              dplyr::rename(module_content_number = module_content_number.y)
          }
        }

        result_with_module <-
          object@merged_module$result_with_module

        # if(translation){
        #   result_with_module <-
        #     result_with_module |>
        #     dplyr::select(-module_annotation) |>
        #     dplyr::rename(module_annotation = module_annotation_trans)
        #
        #   graph_data <-
        #     graph_data |>
        #     tidygraph::activate(what = "nodes") |>
        #     dplyr::select(-module_annotation)
        #
        #   module_annotation <-
        #     result_with_module$module_annotation[match(igraph::vertex_attr(graph_data)$node,
        #                                                result_with_module$node)]
        #
        #   graph_data <-
        #     graph_data |>
        #     tidygraph::activate(what = "nodes") |>
        #     dplyr::mutate(module_annotation = module_annotation)
        # }
      }
    }

    graph_data <-
      graph_data |>
      tidygraph::activate(what = "nodes") |>
      dplyr::filter(module_content_number > degree_cutoff) |>
      dplyr::mutate(label = if (query_type == "gene") Description else pathway_name)

    if (igraph::gorder(graph_data) == 0) {
      warning("No functional modules have degree > ", degree_cutoff)
      return(ggplot() +
               geom_blank())
    }


    if (!missing(module_id)) {
      graph_data <-
        graph_data |>
        tidygraph::activate(what = "nodes") |>
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
            df |>
              dplyr::group_by(module) |>
              dplyr::arrange(p_adjust, desc(Count), .by_group = TRUE) |>
              dplyr::slice_head(n = 1) |>
              dplyr::pull(label)
          } else {
            df |>
              dplyr::group_by(module) |>
              dplyr::arrange(desc(abs(NES)), desc(Count), .by_group = TRUE) |>
              dplyr::slice_head(n = 1) |>
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
      graph_data |>
      tidygraph::activate(what = "nodes") |>
      dplyr::pull(label)

    lay <- ggraph::create_layout(graph_data, layout = "fr")

    plot <-
      ggraph::ggraph(lay) +
      ggraph::geom_edge_link(
        aes(width = sim),
        color = "grey",
        alpha = 1,
        show.legend = TRUE
      ) +
      ggraph::geom_node_point(
        aes(fill = module,
            size = if(analysis_type == "enrich_pathway") -log(p_adjust, 10) else abs(NES)),
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

    if (text_all) {
      plot <-
        plot +
        ggraph::geom_node_text(aes(x = x,
                                   y = y,
                                   label = label),
                               size = 3,
                               repel = TRUE)
    } else if (llm_text && level == "functional_module") {
      node_tbl  <- as_tibble(lay)
      centroids <- node_tbl |>
        group_by(module) |>
        summarise(cx = mean(x), cy = mean(y), .groups = "drop")
      module_name_df <- object@merged_module$functional_module_result |> dplyr::select(module, llm_module_name)
      centroids <- centroids |> dplyr::left_join(module_name_df, by = "module")

      plot <-
        plot +
        ggraph::geom_node_text(
          data = centroids,
          aes(x = cx,
              y = cy,
              label = stringr::str_wrap(llm_module_name, 30)),
          check_overlap = TRUE,
          size = 3,
          repel = TRUE
        )
    } else if (text) {
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
# result_with_module |>
#   dplyr::group_by(ONTOLOGY) |>
#   dplyr::summarise(n = n()) |>
#   dplyr::mutate(n = n * 10 / max(n) + 2)
#output the correlation matrix

# if(database == "go"){
#   message("Output correlation matrix plot...")
#   for(ont in c('MF', "BP", "CC")) {
#     cat(ont, " ")
#     show_matrix_cluster(
#       result = result_with_module |> dplyr::mutate(Direction = "UP"),
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
#   unique(module_result$module) |>
#     purrr::map(
#       .f = function(x) {
#         cat(x, " ")
#         number <- module_result |>
#           dplyr::filter(module == x) |>
#           pull(module_content_number) |>
#           as.numeric()
#         if (number == 1) {
#           return(NULL)
#         }
#
#         temp_id <-
#           module_result |>
#           dplyr::filter(module == x) |>
#           dplyr::pull(node) |>
#           stringr::str_split(";") |>
#           `[[`(1) |>
#           pRoloc::goIdToTerm(keepNA = FALSE) |>
#           data.frame(id = ., class = "YES") |>
#           tibble::rownames_to_column(var = "name")
#
#         temp_plot =
#           GOSim::getGOGraph(term = temp_id$name, prune = Inf) |>
#           igraph::igraph.from.graphNEL() |>
#           tidygraph::as_tbl_graph() |>
#           left_join(temp_id, by = "name") |>
#           dplyr::mutate(class = dplyr::case_when(is.na(class) ~ "NO",
#                                           TRUE ~ class))
#
#         plot =
#           temp_plot |>
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
