#' Plot a Network of Modules, Pathways, and Molecules
#'
#' This function generates a network plot based on the input data which includes module results, pathway enrichments, and marker information.
#'
#' @param module_result A data frame of module results.
#' @param result_with_module A data frame of results with module information.
#' @param enrichment_go A data frame of GO enrichment results.
#' @param enrichment_kegg A data frame of KEGG enrichment results.
#' @param enrichment_reactome A data frame of Reactome enrichment results.
#' @param marker_info A data frame of marker information.
#' @param including_molecule Logical, whether to include molecules in the plot.
#' @param circle Logical, whether to use a circular layout for the plot.
#' @param node_color A vector of colors for different node types.
#' @param text_size Numeric, size of the text in the plot.
#' @param arrange_position Logical, whether to arrange position of nodes.
#' @param position_ratio Numeric, ratio for arranging position of nodes.
#'
#' @return A `ggplot` object representing the network plot.
#' @export
#'
#' @examples
#' \dontrun{
#' plot_network(module_result,
#' result_with_module,
#' enrichment_go,
#' enrichment_kegg,
#' enrichment_reactome,
#' marker_info)
#' }

plot_network <-
  function(module_result = module_result_all,
           result_with_module = result_with_module,
           enrichment_go = enrichment_go,
           enrichment_kegg = enrichment_kegg,
           enrichment_reactome = enrichment_reactome,
           marker_info = marker,
           including_molecule = TRUE,
           circle = FALSE,
           node_color = c(
             "Function_module" = "#F05C3BFF",
             "Module" = "#46732EFF",
             "Pathway" = "#197EC0FF",
             "Molecule" = "#3B4992FF"
           ),
           text_size = 3,
           arrange_position = TRUE,
           position_ratio = 0.95) {
    edge_color <-
      c(
        "Function_module-Module" = unname(node_color["Function_module"]),
        "Module-Pathway" = unname(node_color["Module"]),
        "Pathway-Molecule" = unname(node_color["Pathway"])
      )

    ####function_module vs module
    edge_data1 <-
      data.frame(from = module_result$module,
                 to = module_result$module_content) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], ";")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(class = "Function_module-Module")

    node_data1 <-
      module_result %>%
      dplyr::select(node = module,
                    annotation = module_annotation,
                    p.adjust,
                    Count) %>%
      dplyr::mutate(Count = as.numeric(Count),
                    class = "Function_module")

    ###module vs pathway
    edge_data2 <-
      result_with_module %>%
      dplyr::select(from = node,
                    to = pathway_id) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], ";")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(class = "Module-Pathway") %>%
      dplyr::filter(from %in% edge_data1$to)

    node_data2 <-
      result_with_module %>%
      dplyr::select(node,
                    annotation = module_annotation,
                    p.adjust,
                    Count,
                    database) %>%
      dplyr::mutate(Count = as.numeric(Count),
                    class = "Module") %>%
      dplyr::filter(node %in% edge_data2$from)


    #####pathway vs molecule
    edge_data3_go <-
      enrichment_go@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::filter(ONTOLOGY != "CC") %>%
      dplyr::select(from = ID,
                    to = geneID) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], "/")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    if (nrow(edge_data3_go) > 0) {
      edge_data3_go <-
        edge_data3_go %>%
        dplyr::filter(from %in% edge_data2$to)
    }

    node_data3_go <-
      enrichment_go@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::filter(ONTOLOGY != "CC") %>%
      dplyr::select(node = ID,
                    annotation = Description,
                    p.adjust,
                    Count) %>%
      dplyr::mutate(database = "GO")

    if (nrow(node_data3_go) > 0) {
      node_data3_go <-
        node_data3_go %>%
        dplyr::filter(node %in% edge_data3_go$from)
    }

    edge_data3_kegg <-
      enrichment_kegg@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::select(from = ID,
                    to = geneID) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], "/")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    if (nrow(edge_data3_kegg) > 0) {
      edge_data3_kegg <-
        edge_data3_kegg %>%
        dplyr::filter(from %in% edge_data2$to)
    }

    node_data3_kegg <-
      enrichment_kegg@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::select(node = ID,
                    annotation = Description,
                    p.adjust,
                    Count) %>%
      dplyr::mutate(database = "KEGG")

    if (nrow(node_data3_kegg) > 0) {
      node_data3_kegg <-
        node_data3_kegg %>%
        dplyr::filter(node %in% edge_data3_kegg$from)
    }

    edge_data3_reactome <-
      enrichment_reactome@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::select(from = ID,
                    to = geneID) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], "/")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame()

    if (nrow(edge_data3_reactome) > 0) {
      edge_data3_reactome <-
        edge_data3_reactome %>%
        dplyr::filter(from %in% edge_data2$to)
    }

    node_data3_reactome <-
      enrichment_reactome@result %>%
      dplyr::filter(p.adjust < 0.05 & Count > 5) %>%
      dplyr::select(node = ID,
                    annotation = Description,
                    p.adjust,
                    Count) %>%
      dplyr::mutate(database = "Reactome")

    if (nrow(node_data3_reactome) > 0) {
      node_data3_reactome <-
        node_data3_reactome %>%
        dplyr::filter(node %in% edge_data3_reactome$from)
    }

    node_data3_molecule <-
      marker_info %>%
      dplyr::select(
        node = ensembl,
        uniprot,
        symbol,
        entrezid,
        annotation = genename,
        p.adjust = fdr,
        SAM_score = score
      ) %>%
      dplyr::mutate(Count = 1,
                    class = "Molecule")

    edge_data3_kegg$to <-
      marker_info$ensembl[match(edge_data3_kegg$to, marker_info$uniprot)]

    edge_data3_reactome$to <-
      marker_info$ensembl[match(edge_data3_reactome$to, marker_info$entrezid)]

    edge_data3 <-
      rbind(edge_data3_go,
            edge_data3_kegg,
            edge_data3_reactome) %>%
      dplyr::filter(!is.na(to)) %>%
      dplyr::distinct(from, to, .keep_all = TRUE) %>%
      dplyr::mutate(class = "Pathway-Molecule")

    node_data3_pathway <-
      rbind(node_data3_go,
            node_data3_kegg,
            node_data3_reactome) %>%
      dplyr::mutate(class = "Pathway")

    node_data3_molecule <-
      node_data3_molecule %>%
      dplyr::filter(node %in% edge_data3$to)

    node_data3 <-
      node_data3_pathway %>%
      dplyr::full_join(node_data3_molecule,
                       by = intersect(
                         colnames(node_data3_pathway),
                         colnames(node_data3_molecule)
                       ))
    edge_data <-
      rbind(edge_data1,
            edge_data2,
            edge_data3)

    node_data <-
      node_data1 %>%
      dplyr::full_join(node_data2,
                       by = intersect(colnames(.),
                                      colnames(node_data2))) %>%
      dplyr::full_join(node_data3,
                       by = intersect(colnames(.),
                                      colnames(node_data3)))

    node_data <-
      node_data %>%
      dplyr::filter(node %in% c(edge_data$from, edge_data$to))

    edge_data <-
      edge_data %>%
      dplyr::filter(from %in% node_data$node &
                      to %in% node_data$node)

    # table(node_data$class)
    #
    # sum(!edge_data$from %in% node_data$node)
    # sum(!edge_data$to %in% node_data$node)
    # sum(!node_data$node %in% c(edge_data$from,
    #                            edge_data$to))

    if (!including_molecule) {
      temp_node_data <-
        node_data %>%
        dplyr::filter(class != "Molecule")
      temp_edge_data <-
        edge_data %>%
        dplyr::filter(class != "Pathway-Molecule")

      total_graph <-
        tidygraph::tbl_graph(nodes = temp_node_data,
                             edges = temp_edge_data,
                             directed = FALSE)

      g <- total_graph
      V(g)$type <- bipartite_mapping(g)$type
      coords <-
        ggraph::create_layout(g, layout = "bipartite")
      coords$index = 1:nrow(coords)

      coords$x <-
        coords$x + 1

      if (circle) {
        coords$y[coords$class == "Function_module"] <- 0
        coords$y[coords$class == "Module"] <- 1
        coords$y[coords$class == "Pathway"] <- 2
        coords$y[coords$class == "Molecule"] <- 3

        if (arrange_position) {
          coords <-
            arrange_coords(coords = coords,
                           ratio = position_ratio)
        }

        coords <-
          coords %>%
          dplyr::select(x, y) %>%
          dplyr::mutate(
            theta = x / (max(x) + 1) * 2 * pi,
            r = y + 1,
            x = r * cos(theta),
            y = r * sin(theta)
          )

        my_graph <-
          ggraph::create_layout(
            graph = g,
            layout = "manual",
            x = coords$x,
            y = coords$y
            # node.position = coords
          )

        plot <-
          ggraph(my_graph,
                 layout = 'bipartite') +
          geom_edge_diagonal(
            strength = 1,
            aes(color = class),
            edge_width = 0.5,
            alpha = 1,
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class,
                size = Count),
            shape = 21,
            alpha = 1,
            show.legend = TRUE
          ) +
          geom_node_text(
            aes(
              x = x * 1.03,
              y = y * 1.03,
              hjust = ifelse(class == "Pathway", "outward", 'inward'),
              angle = -((-node_angle(x, y) + 90) %% 180) + 90,
              label = annotation
            ),
            size = text_size,
            show.legend = FALSE
          ) +
          scale_fill_manual(values = node_color) +
          ggraph::scale_edge_color_manual(values = edge_color) +
          scale_size_continuous(range = c(1, 8)) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "right",
            legend.background = element_rect(fill = "transparent", color = NA)
          )

      } else{
        coords$y[coords$class == "Function_module"] <- 3
        coords$y[coords$class == "Module"] <- 2
        coords$y[coords$class == "Pathway"] <- 1
        coords$y[coords$class == "Molecule"] <- 0

        if (arrange_position) {
          coords <-
            arrange_coords(coords = coords,
                           ratio = position_ratio)
        }

        my_graph <-
          ggraph::create_layout(
            graph = g,
            layout = "manual",
            x = coords$x,
            y = coords$y
            # node.position = coords
          )

        plot <-
          ggraph(my_graph,
                 layout = 'bipartite') +
          geom_edge_diagonal(
            strength = 1,
            aes(color = class),
            edge_width = 0.5,
            alpha = 1,
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class,
                size = Count),
            shape = 21,
            alpha = 1,
            show.legend = TRUE
          ) +
          geom_node_text(
            aes(x = x,
                y = y,
                label = annotation),
            hjust = 1,
            angle = 90,
            size = text_size,
            show.legend = FALSE
          ) +
          scale_fill_manual(values = node_color) +
          ggraph::scale_edge_color_manual(values = edge_color) +
          scale_size_continuous(range = c(1, 8)) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "right",
            legend.background = element_rect(fill = "transparent", color = NA)
          )

      }
    } else{
      total_graph <-
        tidygraph::tbl_graph(nodes = node_data,
                             edges = edge_data,
                             directed = FALSE)

      g <- total_graph
      V(g)$type <- bipartite_mapping(g)$type
      coords <-
        ggraph::create_layout(g, layout = "bipartite")
      coords$index = 1:nrow(coords)

      if (circle) {
        coords$y[coords$class == "Function_module"] <- 0
        coords$y[coords$class == "Module"] <- 1
        coords$y[coords$class == "Pathway"] <- 2
        coords$y[coords$class == "Molecule"] <- 3

        if (arrange_position) {
          coords <-
            arrange_coords(coords = coords,
                           ratio = position_ratio)
        }

        coords <-
          coords %>%
          dplyr::select(x, y) %>%
          dplyr::mutate(
            theta = x / (max(x) + 1) * 2 * pi,
            r = y + 1,
            x = r * cos(theta),
            y = r * sin(theta)
          )

        my_graph <-
          ggraph::create_layout(
            graph = g,
            layout = "manual",
            x = coords$x,
            y = coords$y
            # node.position = coords
          )

        plot <-
          ggraph(my_graph,
                 layout = 'bipartite') +
          geom_edge_diagonal(
            strength = 1,
            aes(color = class),
            edge_width = 0.5,
            alpha = 1,
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class,
                size = Count),
            shape = 21,
            alpha = 1,
            show.legend = TRUE
          ) +
          geom_node_text(
            aes(
              x = x * 1.03,
              y = y * 1.03,
              hjust = ifelse(class == "Pathway", "outward", 'inward'),
              angle = -((-node_angle(x, y) + 90) %% 180) + 90,
              label = ifelse(class == "Molecule", NA, annotation)
            ),
            size = text_size,
            show.legend = FALSE
          ) +
          scale_fill_manual(values = node_color) +
          ggraph::scale_edge_color_manual(values = edge_color) +
          scale_size_continuous(range = c(1, 8)) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "right",
            legend.background = element_rect(fill = "transparent", color = NA)
          )

      } else{
        coords$y[coords$class == "Function_module"] <- 3
        coords$y[coords$class == "Module"] <- 2
        coords$y[coords$class == "Pathway"] <- 1
        coords$y[coords$class == "Molecule"] <- 0

        if (arrange_position) {
          coords <-
            arrange_coords(coords = coords,
                           ratio = position_ratio)
        }

        my_graph <-
          ggraph::create_layout(
            graph = g,
            layout = "manual",
            x = coords$x,
            y = coords$y
            # node.position = coords
          )

        plot <-
          ggraph(my_graph,
                 layout = 'bipartite') +
          geom_edge_diagonal(
            strength = 1,
            aes(color = class),
            edge_width = 0.5,
            alpha = 1,
            show.legend = FALSE
          ) +
          geom_node_point(
            aes(fill = class,
                size = Count),
            shape = 21,
            alpha = 1,
            show.legend = TRUE
          ) +
          geom_node_text(
            aes(
              x = x,
              y = y,
              label = ifelse(class == "Molecule", NA, annotation)
            ),
            hjust = 1,
            angle = 90,
            size = text_size,
            show.legend = FALSE
          ) +
          scale_fill_manual(values = node_color) +
          ggraph::scale_edge_color_manual(values = edge_color) +
          scale_size_continuous(range = c(1, 8)) +
          ggraph::theme_graph() +
          theme(
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.position = "right",
            legend.background = element_rect(fill = "transparent", color = NA)
          )

      }


    }
    plot
  }
