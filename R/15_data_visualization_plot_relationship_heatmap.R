# setwd(r4projects::get_project_wd())
# rm(list = ls())
# library(mapa)
# source("R/6_utils.R")
# source("R/15_data_visualization_plot_relationship_network.R")
#
# load("demo_data/gene_ora_res/expression_data.rda")
# load("demo_data/gene_ora_res/llm_interpreted_functional_module.rda")
# object <- llm_interpreted_functional_module
# load("demo_data/gene_ora_res/enriched_functional_module.rda")
# object <- enriched_functional_module_res
#
# module_content_number_cutoff = 2
# module_ids = NULL
#
# module_ids = c("Functional_module_48", "Functional_module_61", "Functional_module_65", "Functional_module_72", "Functional_module_75")
# module_ids = c("Functional_module_48", "Functional_module_61", "Functional_module_65")
# module_ids = c("Functional_module_48", "Functional_module_61")
# module_ids = c("Functional_module_48")
# module_ids = c("Functional_module_2")
# module_content_number_cutoff = NULL
#
# {
#   object <- llm_interpreted_functional_module
#   # level = "pathway"
#   level = "molecule"
#   expression_data = expression_data
#   # module_content_number_cutoff = module_content_number_cutoff
#   module_ids = module_ids
#   scale_expression_data = TRUE
#   dist_method = "euclidean"
#   hclust_method = "complete"
#   cluster_rows = FALSE
#   show_cluster_tree = TRUE
#   wordcloud = TRUE
#   llm_text = FALSE
#   functional_module_position_limits = c(0.2, 0.8)
#   pathway_position_limits = c(0.1, 0.9)
#   molecule_position_limits = c(0, 1)
#   functional_module_color = "#DD4124FF"
#   pathway_color = "#43B284FF"
#   molecule_color = "#FAB255FF"
#   functional_module_text = TRUE
#   pathway_text = TRUE
#   molecule_text = FALSE
#   text_size = 3
#   text_col = c("#FE9B00FF", "#D8443CFF", "#9B3441FF", "#DE597CFF",
#                "#E87B89FF", "#E6A2A6FF", "#AA7AA1FF", "#9F5691FF",
#                "#633372FF", "#1F6E9CFF", "#2B9B81FF", "#92C051FF" )
#   num_text_show = 15
#   min_font_size = 8
#   max_font_size = 12
#   text_box_color = "grey"
#   text_box_fill = "white"
#   text_box_width = 4
#   heatmap_height_ratios = c(1,100)
#   network_height_ratios = NULL
# }

# plot_relationship_heatmap(
#   object = llm_interpreted_functional_module,
#   level = "pathway",
#   expression_data = expression_data,
#   module_ids = module_ids,
#   llm_text = TRUE
# )

#' Plot a Combined Relationship Network, Heatmap, and Word Cloud
#'
#' Create a composite visual that links a bipartite/tripartite network
#' of functional modules, pathways and/or molecules to an expression
#' heat‑map, with optional word‑cloud annotations that summarise each
#' functional module.
#'
#' @param object A **functional_module** object containing analysis results.
#'
#' @param level Character string selecting the visualisation level;
#'   must be `"pathway"` or `"molecule"`. Defaults to `"pathway"`.
#'
#' @param expression_data A data frame of expression values. It must
#'   contain an `id` column (ENSEMBL for genes/proteins, HMDB ID or
#'   KEGG ID for metabolites) followed by one column per sample/group.
#'
#' @param module_content_number_cutoff Optional numeric cutoff; keep
#'   only modules with more pathways than this value.
#'
#' @param module_ids Optional character vector of module IDs to plot.
#'   When supplied, `module_content_number_cutoff` is ignored.
#'
#' @param scale_expression_data Logical; if `TRUE` (default) the heat‑map
#'   matrix is scaled row‑wise.
#'
#' @param dist_method Distance metric for clustering heat‑map rows when
#'   `level = "molecule"` and `cluster_rows = TRUE`. Passed to
#'   [stats::dist()]. Default is `"euclidean"`.
#'
#' @param hclust_method Agglomeration method for hierarchical clustering
#'   when `level = "molecule"` and `cluster_rows = TRUE`. Passed to
#'   [stats::hclust()]. Default is `"complete"`.
#'
#' @param cluster_rows Logical; if `TRUE`, molecules (rows) are clustered
#'   in the heat‑map when `level = "molecule"`. Defaults to `TRUE`.
#'
#' @param show_cluster_tree Logical; when `level = "molecule"` and
#'   `cluster_rows = TRUE`, controls whether the row dendrogram is shown
#'   (`TRUE`, default) or hidden (`FALSE`).
#'
#' @param wordcloud Logical; if `TRUE` (default) a word‑cloud is drawn
#'   for each module and displayed as a row annotation (only when
#'   `level = "pathway"`).
#'
#' @param plot_widths Numeric length‑2 vector specifying the relative
#'   widths of the network plot and heatmap. Defaults to `c(1, 2)`.
#'
#' @param llm_text Logical; if `TRUE`, the word‑cloud is built from
#'   LLM‑generated summaries, otherwise from concatenated pathway
#'   descriptions. Defaults to `FALSE`.
#'
#' @param functional_module_position_limits,pathway_position_limits,molecule_position_limits
#'   Numeric length‑2 vectors setting the
#'   relative x‑axis span for the respective node types. Defaults are
#'   `c(0.2, 0.8)`, `c(0.1, 0.9)` and `c(0, 1)`.
#'
#' @param functional_module_color,pathway_color,molecule_color Colours
#'   for the three node classes. Defaults are
#'   `"#DD4124FF"`, `"#43B284FF"`, and `"#FAB255FF"`.
#'
#' @param functional_module_text,pathway_text,molecule_text Logical flags
#'   controlling whether text labels are drawn for each node class.
#'   Defaults are `TRUE`, `TRUE`, and `FALSE` respectively.
#'
#' @param text_size Numeric font size for node labels. Default is `3`.
#'
#' @param text_col A vector of colours to cycle through in the word‑cloud
#'   (only for `level = "pathway"`). Default is a palette from
#'   **MetBrewer::Signac**.
#'
#' @param num_text_show,max_font_size,min_font_size Word‑cloud layout
#'   controls. Show up to `num_text_show` words (default 15) with font
#'   sizes scaled between `min_font_size` (8) and `max_font_size` (12).
#'
#' @param text_box_color,text_box_fill,text_box_width Appearance of the
#'   annotation box that contains each word‑cloud. Defaults are
#'   `"grey"`, `"white"`, and `4` cm.
#'
#' @param heatmap_height_ratios Numeric length‑2 vector controlling the
#'   relative heights of the heatmap components (spacing vs heatmap)
#'   Defaults to `c(1, 100)`. Used for vertical alignment with the relationship plot.
#'
#' @param network_height_ratios Optional numeric length‑2 vector controlling
#'   the relative heights of the relationship network plot components (relationship plot vs spacing).
#'   If `NULL` (default), automatically calculated based on
#'   the number of rows in heatmap and modules. Used for vertical
#'   alignment with the heatmap plot.
#'
#' @return A **patchwork** object combining the network graph and
#'   heat‑map, ready for further composition or direct plotting.
#'
#' @author Yifei Ge <yifeii.ge@outlook.com>
#'
#' @details
#' **Pathway level** – edges connect functional modules to their
#' pathways. The heat‑map shows the mean expression of molecules within
#' each pathway and optional word‑clouds give functional summaries.
#'
#' **Molecule level** – edges connect modules → pathways → individual
#' molecules. The heat‑map displays each molecule’s expression profile.
#' Row clustering and dendrogram display can be toggled with
#' `cluster_rows` and `show_cluster_tree`.
#'
#' The network layout is manually spaced to improve readability.
#'
#' @examples
#' \dontrun{
#' # Load data
#' load("demo_data/input_data.rda")
#' load("demo_data/llm_interpreted_functional_module.rda")
#'
#' # Pathway‑level visual
#' plot_relationship_heatmap(
#'   object      = llm_interpreted_functional_module,
#'   level       = "pathway",
#'   expression_data = input_data,
#'   module_content_number_cutoff = 5
#' )
#'
#' # Molecule‑level visual with clustering
#' plot_relationship_heatmap(
#'   object      = llm_interpreted_functional_module,
#'   level       = "molecule",
#'   expression_data = input_data,
#'   module_ids  = c("Functional_module_10", "Functional_module_25"),
#'   cluster_rows      = TRUE,
#'   show_cluster_tree = TRUE
#' )
#' }
#'
#' @import dplyr ggplot2
#' @importFrom tidygraph tbl_graph
#' @importFrom ggraph ggraph geom_edge_diagonal geom_node_point
#'   geom_node_text theme_graph create_layout scale_edge_color_manual
#' @importFrom igraph V bipartite_mapping
#' @importFrom ComplexHeatmap pheatmap rowAnnotation anno_textbox
#' @importFrom grid gpar unit
#' @importFrom patchwork plot_layout
#' @importFrom ggplotify as.ggplot
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_wrap str_split str_remove
#' @importFrom textstem lemmatize_words
#'
#' @export

plot_relationship_heatmap <-
  function(object,
           level = c("pathway", "molecule"),
           expression_data,
           module_content_number_cutoff = NULL,
           module_ids = NULL,
           scale_expression_data = TRUE,
           dist_method = "euclidean",
           hclust_method = "complete",
           cluster_rows = FALSE,
           show_cluster_tree = TRUE,
           wordcloud = TRUE,
           llm_text = FALSE,
           plot_widths = c(1,2),
           functional_module_position_limits = c(0.2, 0.8),
           pathway_position_limits = c(0.1, 0.9),
           molecule_position_limits = c(0, 1),
           functional_module_color = "#DD4124FF",
           pathway_color = "#43B284FF",
           molecule_color = "#FAB255FF",
           functional_module_text = TRUE,
           pathway_text = TRUE,
           molecule_text = FALSE,
           text_size = 3,
           text_col = c("#FE9B00FF", "#D8443CFF", "#9B3441FF", "#DE597CFF",
                        "#E87B89FF", "#E6A2A6FF", "#AA7AA1FF", "#9F5691FF",
                        "#633372FF", "#1F6E9CFF", "#2B9B81FF", "#92C051FF" ),
           num_text_show = 15,
           min_font_size = 8,
           max_font_size = 12,
           text_box_color = "grey",
           text_box_fill = "white",
           text_box_width = 4,
           heatmap_height_ratios = c(1,100),
           network_height_ratios = NULL) {

    level <- match.arg(level)

    # 1. Data preparation ====
    if (!is.null(module_ids)) {
      object@merged_module$functional_module_result  <-
        object@merged_module$functional_module_result |>
        dplyr::filter(module %in% module_ids)
    } else if (!is.null(module_content_number_cutoff)) {
      object@merged_module$functional_module_result  <-
        object@merged_module$functional_module_result |>
        dplyr::filter(module_content_number > module_content_number_cutoff)
    }

    mods <- object@merged_module$functional_module_result$module

    # 2. Create network and heatmap matrix ====
    if (level == "molecule") {
      include_functional_modules <- TRUE
      include_pathways <- TRUE
      include_molecules <- TRUE
      graph_data <-
        create_relation_network(
          object = object,
          include_functional_modules = include_functional_modules,
          include_modules = FALSE,
          include_pathways = include_pathways,
          include_molecules = include_molecules,
          include_variables = FALSE
        )
      node_data <- graph_data$node_data
      molecules <- node_data |>
        dplyr::filter(class == "Molecule") |>
        dplyr::pull(node)
      expre_dt <- expression_data |> dplyr::filter(id %in% molecules)
      heatmap_matrix <- expre_dt |>
        dplyr::summarise(dplyr::across(where(is.numeric), \ (x)
                                       mean(x, na.rm = TRUE)),
                         .by = id) |>
        tibble::column_to_rownames("id") |>
        as.matrix()

      if (scale_expression_data) {
        heatmap_matrix <- t(scale(t(heatmap_matrix)))
      }

    } else if (level == "pathway") {
      include_functional_modules <- TRUE
      include_pathways <- TRUE
      include_molecules <- FALSE
      graph_data <-
        create_relation_network(
          object = object,
          include_functional_modules = include_functional_modules,
          include_modules = FALSE,
          include_pathways = include_pathways,
          include_molecules = include_molecules,
          include_variables = FALSE
        )
      node_data <- graph_data$node_data
      pathways <- node_data |>
        dplyr::filter(class == "Pathway") |>
        dplyr::pull(node)

      query_type <- object@process_info$merge_pathways@parameter$query_type
      if (query_type == "gene") {
        mapped_molecules <-
          object@merged_module$result_with_module |>
          dplyr::filter(node %in% pathways) |>
          dplyr::select(node, geneID) |>
          tidyr::separate_rows(geneID, sep = "/") |>
          dplyr::mutate(id = unify_id_internal(
            geneID,
            variable_info = object@variable_info,
            query_type = query_type
          ))
        heatmap_matrix <-
          mapped_molecules |>
          dplyr::left_join(expression_data, by = "id") |>
          dplyr::select(-geneID) |>
          dplyr::summarise(dplyr::across(where(is.numeric), \ (x)
                                         mean(x, na.rm = TRUE)),
                           .by = node) |>
          tibble::column_to_rownames("node") |>
          as.matrix()

      } else if (query_type == "metabolite") {
        mapped_molecules <-
          object@merged_module$result_with_module |>
          dplyr::filter(node %in% pathways) |>
          dplyr::select(node, mapped_id) |>
          tidyr::separate_rows(mapped_id, sep = "/") |>
          dplyr::rename(id = mapped_id)

        heatmap_matrix <-
          mapped_molecules |>
          dplyr::left_join(expression_data, by = "id") |>
          dplyr::select(-id) |>
          dplyr::summarise(dplyr::across(where(is.numeric), \ (x)
                                         mean(x, na.rm = TRUE)),
                           .by = node) |>
          tibble::column_to_rownames("node") |>
          as.matrix()
      }

      if (scale_expression_data) {
        heatmap_matrix <- t(scale(t(heatmap_matrix)))
      }
    }

    if (nrow(heatmap_matrix) == 0) {
      stop("No expression data available for your selected modules.")
    }

    # 3. Create graph object and layout ====
    if (level == "molecule") {
      if (cluster_rows) {
        dist_heatmap_matrix <- dist(heatmap_matrix, method = dist_method)
        hclust_result <- hclust(dist_heatmap_matrix, method = hclust_method)
        row_label <- hclust_result$labels[hclust_result$order]
        node_dt_molecule <- node_data |>
          dplyr::filter(class == "Molecule")
        node_dt_molecule <-
          node_dt_molecule[rev(match(row_label, node_dt_molecule$node)), ]
        node_data_remove_molecule <- node_data |>
          dplyr::filter(!class == "Molecule")
        updated_node_data <-
          rbind(node_data_remove_molecule, node_dt_molecule)

        g <- tidygraph::tbl_graph(nodes    = updated_node_data,
                                  edges    = graph_data$edge_data,
                                  directed = FALSE)
      } else {
        network_node_label <- node_data$node[node_data$class == "Molecule"]
        heatmap_matrix <-
          heatmap_matrix[rev(match(network_node_label, rownames(heatmap_matrix))), ]
        g <- graph_data$graph_data
      }

    } else if (level == "pathway") {
      if ("node" %in% colnames(object@merged_module$functional_module_result)) {
        pathway_in_module <-
          object@merged_module$functional_module_result |>
          dplyr::select(module, node) |>
          tidyr::separate_rows(node, sep = ";")
        fm_order <- unique(pathway_in_module$module)
      } else if ("pathway_id" %in% colnames(object@merged_module$functional_module_result)) {
        pathway_in_module <-
          object@merged_module$functional_module_result |>
          dplyr::select(module, pathway_id) |>
          tidyr::separate_rows(pathway_id, sep = ";") |>
          dplyr::rename(node = pathway_id)
        fm_order <- unique(pathway_in_module$module)
      }

      node_dt_fm <- node_data |> dplyr::filter(class == "Functional_module")
      node_dt_fm <- node_dt_fm[rev(match(fm_order, node_dt_fm$node)), ]

      node_dt_pathway <- node_data |> dplyr::filter(class == "Pathway")
      node_dt_pathway <-
        node_dt_pathway[rev(match(pathway_in_module$node, node_dt_pathway$node)), ]

      updated_node_data <- rbind(node_dt_fm, node_dt_pathway)

      g <- tidygraph::tbl_graph(nodes = updated_node_data,
                                edges = graph_data$edge_data,
                                directed = FALSE)

      heatmap_matrix <-
        heatmap_matrix[match(pathway_in_module$node, rownames(heatmap_matrix)), ]
    }

    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    coords <- ggraph::create_layout(g, layout = "bipartite")
    coords$index = 1:nrow(coords)
    coords$x <- coords$x + 1

    if (level == "molecule") {
      coords$y[coords$class == "Functional_module"] <- 2.5
      coords$y[coords$class == "Pathway"] <- 1.5
      coords$y[coords$class == "Molecule"] <- 1
      coords$y <- coords$y * 0.5

      x_span <- max(coords$x) - min(coords$x)
      coords$x[coords$class == "Functional_module"] <-
        seq(
          min(coords$x) + functional_module_position_limits[1] * x_span,
          min(coords$x) + functional_module_position_limits[2] * x_span,
          length.out = length(coords$x[coords$class == "Functional_module"])
        )

      coords$x[coords$class == "Pathway"] <-
        seq(
          min(coords$x) + pathway_position_limits[1] * x_span,
          min(coords$x) + pathway_position_limits[2] * x_span,
          length.out = length(coords$x[coords$class == "Pathway"])
        )

      coords$x[coords$class == "Molecule"] <-
        seq(
          min(coords$x) + molecule_position_limits[1] * x_span,
          min(coords$x) + molecule_position_limits[2] * x_span,
          length.out = length(coords$x[coords$class == "Molecule"])
        )

    } else if (level == "pathway") {
      coords$y[coords$class == "Functional_module"] <- 1.25
      coords$y[coords$class == "Pathway"] <- 1
      coords$y <- coords$y * 0.5

      x_span <- max(coords$x) - min(coords$x)
      coords$x[coords$class == "Functional_module"] <-
        seq(
          min(coords$x) + functional_module_position_limits[1] * x_span,
          min(coords$x) + functional_module_position_limits[2] * x_span,
          length.out = length(coords$x[coords$class == "Functional_module"])
        )

      coords$x[coords$class == "Pathway"] <-
        seq(
          min(coords$x) + pathway_position_limits[1] * x_span,
          min(coords$x) + pathway_position_limits[2] * x_span,
          length.out = length(coords$x[coords$class == "Pathway"])
        )
    }

    coords <- coords |> dplyr::arrange(index)

    my_graph <-
      ggraph::create_layout(
        graph = g,
        layout = "manual",
        x = coords$x,
        y = coords$y
      )

    # 4. Create network plot ====
    node_color <-
      c(
        "Functional_module" = functional_module_color,
        "Pathway" = pathway_color,
        "Molecule" = molecule_color
      )

    edge_color <-
      c(
        "Functional_module-Pathway" = unname(node_color["Functional_module"]),
        "Pathway-Molecule" = unname(node_color["Pathway"])
      )

    plot <-
      ggraph::ggraph(my_graph,
                     layout = 'bipartite') +
      ggraph::geom_edge_diagonal(
        strength = 1,
        aes(color = class),
        edge_width = 0.5,
        alpha = 1,
        show.legend = FALSE
      ) +
      ggraph::geom_node_point(
        aes(fill = class,
            size = Count),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      scale_fill_manual(values = node_color) +
      ggraph::scale_edge_color_manual(values = edge_color,
                                      guide = "none") +
      scale_size_continuous(range = c(2, 8)) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    text_parameter <-
      data.frame(
        class = c("Functional_module",
                  "Pathway",
                  "Molecule"),
        include = c(
          include_functional_modules,
          include_pathways,
          include_molecules
        ),
        text = c(
          functional_module_text,
          pathway_text,
          molecule_text
        )
      ) |>
      dplyr::filter(include)

    if (sum(text_parameter$text) > 0) {
      if (level == "pathway") {
        plot <-
          plot +
          ggraph::geom_node_text(
            aes(
              x = x * 1.02,
              y = y * 0.85,
              label = ifelse(
                class == "Functional_module",
                stringr::str_wrap(annotation, 30),
                node
              )
            ),
            hjust = 1,
            angle = 0,
            size = text_size,
            show.legend = FALSE
          )
      } else if (level == "molecule") {
        plot <-
          plot +
          ggraph::geom_node_text(
            aes(
              x = x * 1.02,
              y = y * 0.85,
              label = ifelse(
                class %in% text_parameter$class[text_parameter$text],
                stringr::str_wrap(annotation, 40),
                NA
              )
            ),
            hjust = 1,
            angle = 0,
            size = text_size,
            show.legend = FALSE
          )
      }
    }

    # 5. Create word cloud data (if applicable) ====
    split_words <- NULL
    if (wordcloud && level == "pathway") {
      if (llm_text) {
        if (!("llm_interpret_module" %in% names(object@process_info))) {
          stop(
            "Please perform llm_interpret_module() first to get LLM generated function summary."
          )
        }
        desc <-
          object@llm_module_interpretation[match(mods, names(object@llm_module_interpretation))] |>
          purrr::map(function(x) {
            fm_summary <- x$generated_name$summary
            gsub("\\.", "", fm_summary)
          }) |>
          unlist() |>
          unname()
      } else {
        desc <- object@merged_module$functional_module_result$Description
        desc <- gsub(";", " ", desc)
      }

      split_words <- setNames(strsplit(desc, " "), mods) |>
        purrr::map(function(x) {
          clean_words <- x |>
            stringr::str_split(pattern = " ") |>
            unlist() |>
            stringr::str_remove(",$")  |>
            stringr::str_remove("^\\(") |>
            stringr::str_remove("\\)$")
          clean_words <- textstem::lemmatize_words(tolower(clean_words))

          filtered_words <-
            clean_words[!(clean_words %in% remove_words)]

          word_freq <-
            as.data.frame(sort(table(filtered_words), decreasing = TRUE),
                          stringsAsFactors = FALSE)

          if (nrow(word_freq) == 0)
            return(NULL)

          box_text <- word_freq |>
            dplyr::rename(text = filtered_words) |>
            dplyr::mutate(
              fontsize = min_font_size + (max_font_size - min_font_size) * (Freq - min(Freq)) / (max(Freq) - min(Freq))
            ) |>
            dplyr::select(-Freq)

          current_text_col <- colorRampPalette(text_col)(num_text_show)
          if (nrow(box_text) < num_text_show) {
            box_text <-
              dplyr::mutate(box_text, col = current_text_col[1:nrow(box_text)])
          } else {
            box_text <- box_text |>
              head(num_text_show) |>
              dplyr::mutate(col = current_text_col)
          }
          return(box_text)
        })
    }

    # 6. Create heatmap plot ====
    if (level == "molecule") {
      if (show_cluster_tree) {
        heatmap_plot <-
          ComplexHeatmap::pheatmap(
            heatmap_matrix,
            cluster_rows = TRUE,
            row_dend_side = "right",
            cluster_cols = FALSE,
            display_numbers = FALSE,
            show_colnames = TRUE,
            show_rownames = FALSE,
            angle_col = "45",
            border_color = "white",
            border = TRUE,
            heatmap_legend_param = list(title = "Value",
                                        direction = "vertical")
          )
      } else {
        heatmap_plot <-
          ComplexHeatmap::pheatmap(
            heatmap_matrix,
            cluster_rows = FALSE,
            cluster_cols = FALSE,
            display_numbers = FALSE,
            show_colnames = TRUE,
            show_rownames = FALSE,
            angle_col = "45",
            border_color = "white",
            border = TRUE,
            heatmap_legend_param = list(title = "Value",
                                        direction = "vertical")
          )
      }
    } else if (level == "pathway") {
      r_split <-
        factor(pathway_in_module$module, levels = unique(pathway_in_module$module))

      row_ha <- NULL
      if (!is.null(split_words)) {
        row_ha <- ComplexHeatmap::rowAnnotation(textbox = ComplexHeatmap::anno_textbox(
          r_split,
          split_words,
          max_width = grid::unit(text_box_width, "cm"),
          background_gp = grid::gpar(fill = text_box_fill, col = text_box_color)
        ))
      }

      heatmap_plot <-
        ComplexHeatmap::pheatmap(
          heatmap_matrix,
          cluster_rows = FALSE,
          # row_split = r_split,
          # row_title = NULL,
          cluster_cols = FALSE,
          display_numbers = FALSE,
          show_colnames = TRUE,
          show_rownames = FALSE,
          right_annotation = row_ha,
          angle_col = "90",
          border_color = "white",
          border = TRUE,
          heatmap_legend_param = list(title = "Value",
                                      direction = "vertical")
        )
    }

    # 7. Combine heatmap and network ====
    heatmap_ggplot <- ggplotify::as.ggplot(heatmap_plot) +
      theme(plot.margin = margin(0, 0, 0, 0))

    reversed_plot <- plot +
      coord_flip() +
      scale_y_reverse() +
      theme(plot.margin = margin(0, 0, 0, 0))

    place_holder_plot <- ggplot() +
      annotate("text", x = 0, y = 0, label = "Plot Place holder",
               size = 0.5, fontface = "bold", color = "white") +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    if (is.null(network_height_ratios)) {
      add_place_holder_1 <- reversed_plot/place_holder_plot + plot_layout(heights = c(nrow(heatmap_matrix), length(module_ids)))
    } else {
      add_place_holder_1 <- reversed_plot/place_holder_plot + plot_layout(heights = network_height_ratios)
    }

    add_place_holder_2 <- place_holder_plot/heatmap_ggplot + plot_layout(heights = heatmap_height_ratios)

    combined_plot <-
      (add_place_holder_1 | add_place_holder_2) +
      patchwork::plot_layout(ncol = 2, nrow = 1,
                             widths = plot_widths,
                             guides = "collect")

    return(combined_plot)
  }


