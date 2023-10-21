# setwd(r4projects::get_project_wd())
# source("R/8-functional_module_class.R")
# source("R/6-utils.R")
# setwd("demo_data/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
# library(igraph)
#
# load("result/enriched_functional_module")
#
# object <-
#   enriched_functional_module
#
# object@merged_module$functional_module_result <-
#   head(object@merged_module$functional_module_result, 2)
#
# {
#   include_functional_modules = TRUE
#   include_modules = FALSE
#   include_pathways = TRUE
#   include_molecules = TRUE
#   functional_module_color = "#F05C3BFF"
#   module_color = "#46732EFF"
#   pathway_color = "#197EC0FF"
#   molecule_color = "#3B4992FF"
#   functional_module_text = TRUE
#   module_text = FALSE
#   pathway_text = TRUE
#   molecule_text = TRUE
#   functional_module_text_size = 3
#   module_text_size = 3
#   pathway_text_size = 3
#   molecule_text_size = 3
#   circular_plot = FALSE
#   functional_module_arrange_position = TRUE
#   module_arrange_position = TRUE
#   pathway_arrange_position = TRUE
#   molecule_arrange_position = TRUE
#   functional_module_position_limits = c(0.3, 0.7)
#   module_position_limits = c(0.2, 0.8)
#   pathway_position_limits = c(0.1, 0.9)
#   molecule_position_limits = c(0, 1)
#   }
#
# plot_relationship_network(
#   object = object,
#   include_functional_modules = include_functional_modules,
#   include_modules = include_modules,
#   include_pathways = include_pathways,
#   include_molecules = include_molecules,
#   functional_module_text = functional_module_text,
#   module_text = module_text,
#   pathway_text = pathway_text,
#   molecule_text = molecule_text,
#   functional_module_text_size = functional_module_text_size,
#   module_text_size = module_text_size,
#   pathway_text_size = pathway_text_size,
#   molecule_text_size = molecule_text_size,
#   circular_plot = circular_plot,
#   functional_module_arrange_position = functional_module_arrange_position,
#   module_arrange_position = module_arrange_position,
#   pathway_arrange_position = pathway_arrange_position,
#   molecule_arrange_position = molecule_arrange_position,
#   functional_module_position_limits = functional_module_position_limits,
#   module_position_limits = module_position_limits,
#   pathway_position_limits = pathway_position_limits,
#   molecule_position_limits = molecule_position_limits
# )



#' Plot Relationship Network
#'
#' This function plots a relationship network that includes functional modules,
#' modules, pathways, and molecules.
#'
#' @param object An object of class "functional_module".
#' @param include_functional_modules Logical; include functional modules in the plot?
#' @param include_modules Logical; include modules in the plot?
#' @param include_pathways Logical; include pathways in the plot?
#' @param include_molecules Logical; include molecules in the plot?
#' @param functional_module_color Character; color for functional modules.
#' @param module_color Character; color for modules.
#' @param pathway_color Character; color for pathways.
#' @param molecule_color Character; color for molecules.
#' @param functional_module_text Logical; include functional module text labels?
#' @param module_text Logical; include module text labels?
#' @param pathway_text Logical; include pathway text labels?
#' @param molecule_text Logical; include molecule text labels?
#' @param functional_module_text_size Numeric; size of functional module text labels.
#' @param module_text_size Numeric; size of module text labels.
#' @param pathway_text_size Numeric; size of pathway text labels.
#' @param molecule_text_size Numeric; size of molecule text labels.
#' @param circular_plot Logical; make the plot circular?
#' @param functional_module_arrange_position Logical; should the position of functional modules be arranged?
#' @param module_arrange_position Logical; should the position of modules be arranged?
#' @param pathway_arrange_position Logical; should the position of pathways be arranged?
#' @param molecule_arrange_position Logical; should the position of molecules be arranged?
#' @param functional_module_position_limits Numeric vector; limits for functional module positions.
#' @param module_position_limits Numeric vector; limits for module positions.
#' @param pathway_position_limits Numeric vector; limits for pathway positions.
#' @param molecule_position_limits Numeric vector; limits for molecule positions.
#'
#' @return A ggplot object representing the relationship network.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#'
#'
#' @export

plot_relationship_network <-
  function(object,
           include_functional_modules = TRUE,
           include_modules = TRUE,
           include_pathways = TRUE,
           include_molecules = TRUE,
           functional_module_color = "#F05C3BFF",
           module_color = "#46732EFF",
           pathway_color = "#197EC0FF",
           molecule_color = "#3B4992FF",
           functional_module_text = TRUE,
           module_text = TRUE,
           pathway_text = TRUE,
           molecule_text = FALSE,
           functional_module_text_size = 3,
           module_text_size = 3,
           pathway_text_size = 3,
           molecule_text_size = 3,
           circular_plot = FALSE,
           functional_module_arrange_position = TRUE,
           module_arrange_position = TRUE,
           pathway_arrange_position = TRUE,
           molecule_arrange_position = TRUE,
           functional_module_position_limits = c(0, 1),
           module_position_limits = c(0, 1),
           pathway_position_limits = c(0, 1),
           molecule_position_limits = c(0, 1)) {
    ###at least two classes of nodes

    if (sum(
      c(
        include_functional_modules,
        include_modules,
        include_pathways,
        include_molecules
      )
    ) < 2) {
      stop("The network should includes at least two classes of nodes")
    }

    node_color <-
      c(
        "Functional_module" = functional_module_color,
        "Module" = module_color,
        "Pathway" = pathway_color,
        "Molecule" = molecule_color
      )

    edge_color <-
      c(
        "Functional_module-Module" = unname(node_color["Functional_module"]),
        "Functional_module-Pathway" = unname(node_color["Functional_module"]),
        "Functional_module-Molecule" = unname(node_color["Functional_module"]),
        "Module-Pathway" = unname(node_color["Module"]),
        "Module-Molecule" = unname(node_color["Module"]),
        "Pathway-Molecule" = unname(node_color["Pathway"])
      )

    ###check object and variable_info

    if (!is(object, "functional_module")) {
      stop("object should be functional_module class")
    }

    variable_info <-
      object@variable_info

    check_variable_info(variable_info = variable_info)

    if (all(names(object@process_info) != "merge_pathways")) {
      stop("Please use the merge_pathways() function to process first")
    }

    if (all(names(object@process_info) != "merge_modules")) {
      stop("Please use the merge_modules() function to process first")
    }

    ######create edge_data and node_data and network
    total_graph <-
      create_relation_network(
        object = object,
        include_functional_modules = include_functional_modules,
        include_modules = include_modules,
        include_pathways = include_pathways,
        include_molecules = include_molecules
      )

    g <- total_graph
    igraph::V(g)$type <- igraph::bipartite_mapping(g)$type
    coords <-
      ggraph::create_layout(g, layout = "bipartite")
    coords$index = 1:nrow(coords)

    coords$x <-
      coords$x + 1

    if (circular_plot) {
      coords$y[coords$class == "Functional_module"] <- 0
      coords$y[coords$class == "Module"] <- 1
      coords$y[coords$class == "Pathway"] <- 2
      coords$y[coords$class == "Molecule"] <- 3
    } else{
      coords$y[coords$class == "Functional_module"] <- 3
      coords$y[coords$class == "Module"] <- 2
      coords$y[coords$class == "Pathway"] <- 1
      coords$y[coords$class == "Molecule"] <- 0
    }

    # sort(coords$x[coords$class == "Functional_module"])
    # sort(coords$x[coords$class == "Molecule"])

    ###arrange position
    if (functional_module_arrange_position &
        include_functional_modules) {
      temp_coords <-
        coords %>%
        dplyr::filter(class == "Functional_module") %>%
        dplyr::arrange(x)
      temp_coords$x <-
        seq(
          min(coords$x) + functional_module_position_limits[1] * max(coords$x),
          functional_module_position_limits[2] * max(coords$x),
          length.out = nrow(temp_coords)
        )
      coords[coords$class == "Functional_module", ] <-
        temp_coords
    }

    if (module_arrange_position & include_modules) {
      temp_coords <-
        coords %>%
        dplyr::filter(class == "Module") %>%
        dplyr::arrange(x)
      temp_coords$x <-
        seq(
          min(coords$x) + module_position_limits[1] * max(coords$x),
          module_position_limits[2] * max(coords$x),
          length.out = nrow(temp_coords)
        )
      coords[coords$class == "Module", ] <-
        temp_coords
    }

    if (pathway_arrange_position & include_pathways) {
      temp_coords <-
        coords %>%
        dplyr::filter(class == "Pathway") %>%
        dplyr::arrange(x)
      temp_coords$x <-
        seq(
          min(coords$x) + pathway_position_limits[1] * max(coords$x),
          pathway_position_limits[2] * max(coords$x),
          length.out = nrow(temp_coords)
        )
      coords[coords$class == "Pathway", ] <-
        temp_coords
    }

    if (molecule_arrange_position & include_molecules) {
      temp_coords <-
        coords %>%
        dplyr::filter(class == "Molecule") %>%
        dplyr::arrange(x)
      temp_coords$x <-
        seq(
          min(coords$x) + molecule_position_limits[1] * max(coords$x),
          molecule_position_limits[2] * max(coords$x),
          length.out = nrow(temp_coords)
        )
      coords[coords$class == "Molecule", ] <-
        temp_coords
    }

    coords <-
      coords %>%
      dplyr::arrange(index)

    if (circular_plot) {
      coords <-
        coords %>%
        dplyr::select(x, y) %>%
        dplyr::mutate(
          theta = x / (max(x) + 1) * 2 * pi,
          r = y + 1,
          x = r * cos(theta),
          y = r * sin(theta)
        )
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
      ggraph::scale_edge_color_manual(values = edge_color) +
      scale_size_continuous(range = c(1, 8)) +
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
                  "Module",
                  "Pathway",
                  "Molecule"),
        include = c(
          include_functional_modules,
          include_modules,
          include_pathways,
          include_molecules
        ),
        text = c(
          functional_module_text,
          module_text,
          pathway_text,
          molecule_text
        ),
        size = c(
          functional_module_text_size,
          module_text_size,
          pathway_text_size,
          molecule_text
        )
      ) %>%
      dplyr::filter(include)

    ###which one is the last level?
    tail_position <-
      text_parameter %>%
      dplyr::arrange(dplyr::desc(dplyr::row_number())) %>%
      pull(include) %>%
      `==`(TRUE) %>%
      which() %>%
      min()

    ###add text
    if (sum(text_parameter$text) > 0) {
      tail_position <- nrow(text_parameter) + 1 - tail_position

      text_parameter <-
        text_parameter %>%
        dplyr::mutate(tail_position = FALSE)

      text_parameter$tail_position[tail_position] <- TRUE

      text_size <-
        c(
          Functional_module = functional_module_text_size,
          Module = module_text_size,
          Pathway = pathway_text_size,
          Molecule = molecule_text_size
        )

      if (circular_plot) {
        plot <-
          plot +
          ggraph::geom_node_text(
            aes(
              x = x * 1.03,
              y = y * 1.03,
              hjust = ifelse(class == text_parameter$class[text_parameter$tail_position],
                             "outward", 'inward'),
              angle = -((-ggraph::node_angle(x, y) + 90) %% 180) + 90,
              label = ifelse(class %in% text_parameter$class[text_parameter$text],
                             annotation, NA)
            ),
            size = functional_module_text_size,
            show.legend = FALSE
          )
      } else{
        plot <-
          plot +
          ggraph::geom_node_text(
            aes(
              x = x,
              y = y,
              label = ifelse(class %in% text_parameter$class[text_parameter$text],
                             annotation, NA)
            ),
            hjust = 1,
            angle = 90,
            size = functional_module_text_size,
            show.legend = FALSE
          )
      }
    }
    plot
  }



#' Create a Relationship Network
#'
#' This function constructs a network of relationships between different biological entities.
#' The network can include functional modules, modules, pathways, and molecules.
#'
#' @param object An object of class `functional_module` containing the functional module data.
#' @param include_functional_modules Logical; whether to include functional modules in the network. Defaults to TRUE.
#' @param include_modules Logical; whether to include modules in the network. Defaults to TRUE.
#' @param include_pathways Logical; whether to include pathways in the network. Defaults to TRUE.
#' @param include_molecules Logical; whether to include molecules in the network. Defaults to TRUE.
#'
#' @return A list containing `edge_data` and `node_data` data frames.
#'   - `edge_data` is a data frame containing the relationships between nodes, including their classes.
#'   - `node_data` is a data frame containing information about each node, such as annotation and class.
#'
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @export

create_relation_network <-
  function(object,
           include_functional_modules = TRUE,
           include_modules = TRUE,
           include_pathways = TRUE,
           include_molecules = TRUE) {
    ###at least two classes of nodes

    if (sum(
      c(
        include_functional_modules,
        include_modules,
        include_pathways,
        include_molecules
      )
    ) < 2) {
      stop("The network should includes at least two classes of nodes")
    }

    ###check object and variable_info

    if (!is(object, "functional_module")) {
      stop("object should be functional_module class")
    }

    variable_info <-
      object@variable_info

    check_variable_info(variable_info = variable_info)


    if (all(names(object@process_info) != "merge_pathways")) {
      stop("Please use the merge_pathways() function to process first")
    }

    if (all(names(object@process_info) != "merge_modules")) {
      stop("Please use the merge_modules() function to process first")
    }

    ######create edge_data and node_data
    ####1functional_module vs module
    edge_data1 <-
      tryCatch(
        expr = {
          data.frame(
            from = object@merged_module$functional_module_result$module,
            to = object@merged_module$functional_module_result$module_content
          ) %>%
            apply(1, function(x) {
              data.frame(from = as.character(x[1]),
                         to = as.character(stringr::str_split(x[2], ";")[[1]]))
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            dplyr::mutate(class = "Functional_module-Module")
        },
        error = function(e) {
          NULL
        }
      )

    node_data1 <-
      tryCatch(
        object@merged_module$functional_module_result %>%
          dplyr::select(node = module,
                        annotation = module_annotation,
                        p.adjust,
                        Count) %>%
          dplyr::mutate(Count = as.numeric(Count),
                        class = "Functional_module"),
        error = function(e) {
          NULL
        }
      )


    ###2 module vs pathway
    edge_data2 <-
      tryCatch(
        object@merged_module$result_with_module %>%
          dplyr::select(from = node,
                        to = pathway_id) %>%
          apply(1, function(x) {
            data.frame(from = as.character(x[1]),
                       to = as.character(stringr::str_split(x[2], ";")[[1]]))
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame() %>%
          dplyr::mutate(class = "Module-Pathway") %>%
          dplyr::filter(from %in% edge_data1$to),
        error = function(e) {
          NULL
        }
      )

    node_data2 <-
      tryCatch(
        object@merged_module$result_with_module %>%
          dplyr::select(node,
                        annotation = module_annotation,
                        p.adjust,
                        Count,
                        database) %>%
          dplyr::mutate(Count = as.numeric(Count),
                        class = "Module") %>%
          dplyr::filter(node %in% edge_data2$from),
        error = function(e) {
          NULL
        }
      )


    #####3 pathway vs molecule
    edge_data3_go <-
      tryCatch(
        object@enrichment_go_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.go &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.go
          ) %>%
          dplyr::filter(ONTOLOGY != "CC") %>%
          dplyr::select(from = ID,
                        to = geneID) %>%
          apply(1, function(x) {
            data.frame(from = as.character(x[1]),
                       to = as.character(stringr::str_split(x[2], "/")[[1]]))
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame(),
        error = function(e) {
          NULL
        }
      )

    if (!is.null(edge_data3_go)) {
      if (nrow(edge_data3_go) > 0) {
        edge_data3_go <-
          edge_data3_go %>%
          dplyr::filter(from %in% edge_data2$to)
      }
    }

    node_data3_go <-
      tryCatch(
        object@enrichment_go_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.go &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.go
          ) %>%
          dplyr::filter(ONTOLOGY != "CC") %>%
          dplyr::select(node = ID,
                        annotation = Description,
                        p.adjust,
                        Count) %>%
          dplyr::mutate(database = "GO"),
        error = function(e) {
          NULL
        }
      )

    if (!is.null(node_data3_go)) {
      if (nrow(node_data3_go) > 0) {
        node_data3_go <-
          node_data3_go %>%
          dplyr::filter(node %in% edge_data3_go$from)
      }
    }

    edge_data3_kegg <-
      tryCatch(
        object@enrichment_kegg_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.kegg &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.kegg
          ) %>%
          dplyr::select(from = ID,
                        to = geneID) %>%
          apply(1, function(x) {
            data.frame(from = as.character(x[1]),
                       to = as.character(stringr::str_split(x[2], "/")[[1]]))
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame(),
        error = function(e) {
          NULL
        }
      )

    if (!is.null(edge_data3_kegg)) {
      if (nrow(edge_data3_kegg) > 0) {
        edge_data3_kegg <-
          edge_data3_kegg %>%
          dplyr::filter(from %in% edge_data2$to)
      }
    }

    node_data3_kegg <-
      tryCatch(
        object@enrichment_kegg_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.kegg &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.kegg
          ) %>%
          dplyr::select(node = ID,
                        annotation = Description,
                        p.adjust,
                        Count) %>%
          dplyr::mutate(database = "KEGG"),
        error = function(e) {
          NULL
        }
      )

    if (!is.null(node_data3_kegg)) {
      if (nrow(node_data3_kegg) > 0) {
        node_data3_kegg <-
          node_data3_kegg %>%
          dplyr::filter(node %in% edge_data3_kegg$from)
      }
    }

    edge_data3_reactome <-
      tryCatch(
        object@enrichment_reactome_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.reactome &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.reactome
          ) %>%
          dplyr::select(from = ID,
                        to = geneID) %>%
          apply(1, function(x) {
            data.frame(from = as.character(x[1]),
                       to = as.character(stringr::str_split(x[2], "/")[[1]]))
          }) %>%
          do.call(rbind, .) %>%
          as.data.frame(),
        error = function(e) {
          NULL
        }
      )


    if (!is.null(edge_data3_reactome)) {
      if (nrow(edge_data3_reactome) > 0) {
        edge_data3_reactome <-
          edge_data3_reactome %>%
          dplyr::filter(from %in% edge_data2$to)
      }
    }

    node_data3_reactome <-
      tryCatch(
        object@enrichment_reactome_result@result %>%
          dplyr::filter(
            p.adjust < object@process_info$merge_pathways@parameter$p.adjust.cutoff.reactome &
              Count > object@process_info$merge_pathways@parameter$count.cutoff.reactome
          ) %>%
          dplyr::select(node = ID,
                        annotation = Description,
                        p.adjust,
                        Count) %>%
          dplyr::mutate(database = "Reactome"),
        error = function(e) {
          NULL
        }
      )

    if (!is.null(node_data3_reactome)) {
      if (nrow(node_data3_reactome) > 0) {
        node_data3_reactome <-
          node_data3_reactome %>%
          dplyr::filter(node %in% edge_data3_reactome$from)
      }
    }

    node_data3_molecule <-
      variable_info %>%
      dplyr::select(
        node = ensembl,
        uniprot,
        symbol,
        entrezid,
        p.adjust = fdr,
        SAM_score = score
      ) %>%
      dplyr::mutate(annotation = symbol,
                    Count = 1,
                    class = "Molecule")

    if (!is.null(edge_data3_kegg)) {
      edge_data3_kegg$to <-
        variable_info$ensembl[match(edge_data3_kegg$to, variable_info$uniprot)]

      edge_data3_reactome$to <-
        variable_info$ensembl[match(edge_data3_reactome$to, variable_info$entrezid)]
    }

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

    ####4. functional_module vs pathway
    edge_data4 <-
      data.frame(
        from = object@merged_module$functional_module_result$module,
        to = object@merged_module$functional_module_result$pathway_id
      ) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], ";")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(class = "Functional_module-Pathway")

    ####5. functional_module vs molecule
    edge_data5 <-
      data.frame(
        from = object@merged_module$functional_module_result$module,
        to = object@merged_module$functional_module_result$geneID
      ) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], "/")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(class = "Functional_module-Molecule")

    id2 <-
      variable_info$ensembl[match(edge_data5$to, variable_info$uniprot)]
    id3 <-
      variable_info$ensembl[match(edge_data5$to, variable_info$entrezid)]

    edge_data5$to <-
      data.frame(id1 = edge_data5$to,
                 id2,
                 id3) %>%
      apply(1, function(x) {
        x <- x[!is.na(x)]
        grep("ENSG", x, value = TRUE)
      })


    ####6. module vs molecule
    edge_data6 <-
      rbind(
        data.frame(
          from = object@merged_pathway_go$module_result$module,
          to = object@merged_pathway_go$module_result$geneID
        ),
        data.frame(
          from = object@merged_pathway_kegg$module_result$module,
          to = object@merged_pathway_kegg$module_result$geneID
        ),
        data.frame(
          from = object@merged_pathway_reactome$module_result$module,
          to = object@merged_pathway_reactome$module_result$geneID
        )
      ) %>%
      apply(1, function(x) {
        data.frame(from = as.character(x[1]),
                   to = as.character(stringr::str_split(x[2], "/")[[1]]))
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(class = "Module-Molecule")

    id2 <-
      variable_info$ensembl[match(edge_data6$to, variable_info$uniprot)]
    id3 <-
      variable_info$ensembl[match(edge_data6$to, variable_info$entrezid)]

    edge_data6$to <-
      data.frame(id1 = edge_data6$to,
                 id2,
                 id3) %>%
      apply(1, function(x) {
        x <- x[!is.na(x)]
        grep("ENSG", x, value = TRUE)
      })

    #####graph network
    edge_data <-
      rbind(edge_data1,
            edge_data2,
            edge_data3,
            edge_data4,
            edge_data5,
            edge_data6)

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
    # table(edge_data$class)
    #
    # sum(!edge_data$from %in% node_data$node)
    # sum(!edge_data$to %in% node_data$node)
    # sum(!node_data$node %in% c(edge_data$from,
    #                            edge_data$to))

    #####remove some nodes and edges according to parameters
    temp_node_data <-
      node_data

    temp_edge_data <-
      edge_data

    if (!include_functional_modules) {
      temp_node_data <-
        temp_node_data %>%
        dplyr::filter(class != "Functional_module")
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Module") %>%
        dplyr::filter(class != "Functional_module-Pathway") %>%
        dplyr::filter(class != "Functional_module-Molecule")
    }

    if (!include_modules) {
      temp_node_data <-
        temp_node_data %>%
        dplyr::filter(class != "Module")
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Module") %>%
        dplyr::filter(class != "Module-Pathway") %>%
        dplyr::filter(class != "Module-Molecule")
    }

    if (!include_pathways) {
      temp_node_data <-
        temp_node_data %>%
        dplyr::filter(class != "Pathway")
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Pathway") %>%
        dplyr::filter(class != "Module-Pathway") %>%
        dplyr::filter(class != "Pathway-Molecule")
    }

    if (!include_molecules) {
      temp_node_data <-
        temp_node_data %>%
        dplyr::filter(class != "Molecule")
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Molecule") %>%
        dplyr::filter(class != "Module-Molecule") %>%
        dplyr::filter(class != "Pathway-Molecule")
    }

    ###remove some edges
    if (any(temp_edge_data$class == "Functional_module-Module")) {
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Pathway") %>%
        dplyr::filter(class != "Functional_module-Molecule")
    }

    if (any(temp_edge_data$class == "Module-Pathway")) {
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Module-Molecule")
    }

    if (any(temp_edge_data$class == "Functional_module-Pathway")) {
      temp_edge_data <-
        temp_edge_data %>%
        dplyr::filter(class != "Functional_module-Molecule")
    }


    temp_node_data$class <-
      factor(
        temp_node_data$class,
        levels = c("Functional_module",
                   "Module",
                   "Pathway",
                   "Molecule")[c(
                     include_functional_modules,
                     include_modules,
                     include_pathways,
                     include_molecules
                   )]
      )

    #####network
    temp_node_data <-
      temp_node_data %>%
      dplyr::filter(node %in% unique(c(
        temp_edge_data$from,
        temp_edge_data$to
      )))

    total_graph <-
      tidygraph::tbl_graph(nodes = temp_node_data,
                           edges = temp_edge_data,
                           directed = FALSE)
    total_graph
  }
