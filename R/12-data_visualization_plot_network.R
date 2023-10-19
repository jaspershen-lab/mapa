# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("demo_data.rda")
# load("result/functional_modules/intermediate_data/functional_module")
# load("result/functional_modules/intermediate_data/result_with_module")
#
# load("result/go/enrichment_go_result")
# load("result/kegg/enrichment_kegg_result")
# load("result/reactome/enrichment_reactome_result")
#
# variable_info <-
#   demo_data@variable_info
#
# plot_relationship_network <-
#   function(variable_info,
#            functional_module,
#            result_with_module,
#            enrichment_go_result,
#            enrichment_kegg_result,
#            enrichment_reactome_result,
#            p.adjust.cutoff = 0.05,
#            count.cutoff = 5,
#            include_functional_modules = TRUE,
#            include_modules = TRUE,
#            include_pathways = TRUE,
#            include_molecules = TRUE,
#            functional_module_color = "#F05C3BFF",
#            module_color = "#46732EFF",
#            pathway_color = "#197EC0FF",
#            molecule_color = "#3B4992FF",
#            functional_module_text_size = 3,
#            module_text_size = 3,
#            pathway_text_size = 3,
#            molecule_text_size = 3,
#            circle = FALSE,
#            arrange_position = TRUE,
#            position_ratio = 0.95) {
#     node_color <-
#       c(
#         "Functional_module" = functional_module_color,
#         "Module" = module_color,
#         "Pathway" = pathway_color,
#         "Molecule" = molecule_color
#       )
#
#     edge_color <-
#       c(
#         "Functional_module-Module" = unname(node_color["Functional_module"]),
#         "Module-Pathway" = unname(node_color["Module"]),
#         "Pathway-Molecule" = unname(node_color["Pathway"])
#       )
#
#     ###check variable_info
#     check_variable_info(variable_info = variable_info)
#
#     ######create edge_data and node_data
#     ####1functional_module vs module
#     edge_data1 <-
#       data.frame(from = functional_module$module,
#                  to = functional_module$module_content) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], ";")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(class = "Functional_module-Module")
#
#     node_data1 <-
#       functional_module %>%
#       dplyr::select(node = module,
#                     annotation = module_annotation,
#                     p.adjust,
#                     Count) %>%
#       dplyr::mutate(Count = as.numeric(Count),
#                     class = "Functional_module")
#
#     ###2 module vs pathway
#     edge_data2 <-
#       result_with_module %>%
#       dplyr::select(from = node,
#                     to = pathway_id) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], ";")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(class = "Module-Pathway") %>%
#       dplyr::filter(from %in% edge_data1$to)
#
#     node_data2 <-
#       result_with_module %>%
#       dplyr::select(node,
#                     annotation = module_annotation,
#                     p.adjust,
#                     Count,
#                     database) %>%
#       dplyr::mutate(Count = as.numeric(Count),
#                     class = "Module") %>%
#       dplyr::filter(node %in% edge_data2$from)
#
#     #####3 pathway vs molecule
#     edge_data3_go <-
#       enrichment_go_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::filter(ONTOLOGY != "CC") %>%
#       dplyr::select(from = ID,
#                     to = geneID) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], "/")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#
#     if (nrow(edge_data3_go) > 0) {
#       edge_data3_go <-
#         edge_data3_go %>%
#         dplyr::filter(from %in% edge_data2$to)
#     }
#
#     node_data3_go <-
#       enrichment_go_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::filter(ONTOLOGY != "CC") %>%
#       dplyr::select(node = ID,
#                     annotation = Description,
#                     p.adjust,
#                     Count) %>%
#       dplyr::mutate(database = "GO")
#
#     if (nrow(node_data3_go) > 0) {
#       node_data3_go <-
#         node_data3_go %>%
#         dplyr::filter(node %in% edge_data3_go$from)
#     }
#
#     edge_data3_kegg <-
#       enrichment_kegg_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::select(from = ID,
#                     to = geneID) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], "/")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#
#     if (nrow(edge_data3_kegg) > 0) {
#       edge_data3_kegg <-
#         edge_data3_kegg %>%
#         dplyr::filter(from %in% edge_data2$to)
#     }
#
#     node_data3_kegg <-
#       enrichment_kegg_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::select(node = ID,
#                     annotation = Description,
#                     p.adjust,
#                     Count) %>%
#       dplyr::mutate(database = "KEGG")
#
#     if (nrow(node_data3_kegg) > 0) {
#       node_data3_kegg <-
#         node_data3_kegg %>%
#         dplyr::filter(node %in% edge_data3_kegg$from)
#     }
#
#     edge_data3_reactome <-
#       enrichment_reactome_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::select(from = ID,
#                     to = geneID) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], "/")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame()
#
#     if (nrow(edge_data3_reactome) > 0) {
#       edge_data3_reactome <-
#         edge_data3_reactome %>%
#         dplyr::filter(from %in% edge_data2$to)
#     }
#
#     node_data3_reactome <-
#       enrichment_reactome_result@result %>%
#       dplyr::filter(p.adjust < p.adjust.cutoff &
#                       Count > count.cutoff) %>%
#       dplyr::select(node = ID,
#                     annotation = Description,
#                     p.adjust,
#                     Count) %>%
#       dplyr::mutate(database = "Reactome")
#
#     if (nrow(node_data3_reactome) > 0) {
#       node_data3_reactome <-
#         node_data3_reactome %>%
#         dplyr::filter(node %in% edge_data3_reactome$from)
#     }
#
#     node_data3_molecule <-
#       variable_info %>%
#       dplyr::select(
#         node = ensembl,
#         uniprot,
#         symbol,
#         entrezid,
#         p.adjust = fdr,
#         SAM_score = score
#       ) %>%
#       dplyr::mutate(annotation = symbol,
#                     Count = 1,
#                     class = "Molecule")
#
#     edge_data3_kegg$to <-
#       variable_info$ensembl[match(edge_data3_kegg$to, variable_info$uniprot)]
#
#     edge_data3_reactome$to <-
#       variable_info$ensembl[match(edge_data3_reactome$to, variable_info$entrezid)]
#
#     edge_data3 <-
#       rbind(edge_data3_go,
#             edge_data3_kegg,
#             edge_data3_reactome) %>%
#       dplyr::filter(!is.na(to)) %>%
#       dplyr::distinct(from, to, .keep_all = TRUE) %>%
#       dplyr::mutate(class = "Pathway-Molecule")
#
#     node_data3_pathway <-
#       rbind(node_data3_go,
#             node_data3_kegg,
#             node_data3_reactome) %>%
#       dplyr::mutate(class = "Pathway")
#
#     node_data3_molecule <-
#       node_data3_molecule %>%
#       dplyr::filter(node %in% edge_data3$to)
#
#     node_data3 <-
#       node_data3_pathway %>%
#       dplyr::full_join(node_data3_molecule,
#                        by = intersect(
#                          colnames(node_data3_pathway),
#                          colnames(node_data3_molecule)
#                        ))
#
#     ####4. functional_module vs pathway
#     edge_data4 <-
#       data.frame(from = functional_module$module,
#                  to = functional_module$pathway_id) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], ";")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(class = "Functional_module-Pathway")
#
#     ####5. functional_module vs molecule
#     edge_data5 <-
#       data.frame(from = functional_module$module,
#                  to = functional_module$geneID) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], "/")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(class = "Functional_module-Molecule")
#
#     id2 <-
#       variable_info$ensembl[match(edge_data5$to, variable_info$uniprot)]
#     id3 <-
#       variable_info$ensembl[match(edge_data5$to, variable_info$entrezid)]
#
#     edge_data5$to <-
#       data.frame(id1 = edge_data5$to,
#                  id2,
#                  id3) %>%
#       apply(1, function(x) {
#         x <- x[!is.na(x)]
#         grep("ENSG", x, value = TRUE)
#       })
#
#
#     ####6. module vs molecule
#     edge_data6 <-
#       data.frame(from = functional_module$module,
#                  to = functional_module$geneID) %>%
#       apply(1, function(x) {
#         data.frame(from = as.character(x[1]),
#                    to = as.character(stringr::str_split(x[2], "/")[[1]]))
#       }) %>%
#       do.call(rbind, .) %>%
#       as.data.frame() %>%
#       dplyr::mutate(class = "Functional_module-Molecule")
#
#     id2 <-
#       variable_info$ensembl[match(edge_data5$to, variable_info$uniprot)]
#     id3 <-
#       variable_info$ensembl[match(edge_data5$to, variable_info$entrezid)]
#
#     edge_data5$to <-
#       data.frame(id1 = edge_data5$to,
#                  id2,
#                  id3) %>%
#       apply(1, function(x) {
#         x <- x[!is.na(x)]
#         grep("ENSG", x, value = TRUE)
#       })
#
#
#     edge_data <-
#       rbind(edge_data1,
#             edge_data2,
#             edge_data3)
#
#     node_data <-
#       node_data1 %>%
#       dplyr::full_join(node_data2,
#                        by = intersect(colnames(.),
#                                       colnames(node_data2))) %>%
#       dplyr::full_join(node_data3,
#                        by = intersect(colnames(.),
#                                       colnames(node_data3)))
#
#     node_data <-
#       node_data %>%
#       dplyr::filter(node %in% c(edge_data$from, edge_data$to))
#
#     edge_data <-
#       edge_data %>%
#       dplyr::filter(from %in% node_data$node &
#                       to %in% node_data$node)
#
#     # table(node_data$class)
#     #
#     # sum(!edge_data$from %in% node_data$node)
#     # sum(!edge_data$to %in% node_data$node)
#     # sum(!node_data$node %in% c(edge_data$from,
#     #                            edge_data$to))
#
#     #####remove some nodes and edges according to parameters
#     if (sum(
#       c(
#         include_functional_modules,
#         include_modules,
#         include_pathways,
#         include_molecules
#       )
#     ) < 2) {
#       stop("You need to include at least two type of molecules")
#     }
#
#     if (!include_functional_modules) {
#       temp_node_data <-
#         node_data %>%
#         dplyr::filter(class != "Functional_module")
#       temp_edge_data <-
#         edge_data %>%
#         dplyr::filter(class != "Functional_module-Module")
#     }
#
#     if (!include_modules) {
#       temp_node_data <-
#         node_data %>%
#         dplyr::filter(class != "Functional_module")
#       temp_edge_data <-
#         edge_data %>%
#         dplyr::filter(class != "Functional_module-Module")
#     }
#
#
#
#     if (!include_molecules) {
#       temp_node_data <-
#         node_data %>%
#         dplyr::filter(class != "Molecule")
#       temp_edge_data <-
#         edge_data %>%
#         dplyr::filter(class != "Pathway-Molecule")
#
#       total_graph <-
#         tidygraph::tbl_graph(nodes = temp_node_data,
#                              edges = temp_edge_data,
#                              directed = FALSE)
#
#       g <- total_graph
#       V(g)$type <- bipartite_mapping(g)$type
#       coords <-
#         ggraph::create_layout(g, layout = "bipartite")
#       coords$index = 1:nrow(coords)
#
#       coords$x <-
#         coords$x + 1
#
#       if (circle) {
#         coords$y[coords$class == "Functional_module"] <- 0
#         coords$y[coords$class == "Module"] <- 1
#         coords$y[coords$class == "Pathway"] <- 2
#         coords$y[coords$class == "Molecule"] <- 3
#
#         if (arrange_position) {
#           coords <-
#             arrange_coords(coords = coords,
#                            ratio = position_ratio)
#         }
#
#         coords <-
#           coords %>%
#           dplyr::select(x, y) %>%
#           dplyr::mutate(
#             theta = x / (max(x) + 1) * 2 * pi,
#             r = y + 1,
#             x = r * cos(theta),
#             y = r * sin(theta)
#           )
#
#         my_graph <-
#           ggraph::create_layout(
#             graph = g,
#             layout = "manual",
#             x = coords$x,
#             y = coords$y
#             # node.position = coords
#           )
#
#         plot <-
#           ggraph(my_graph,
#                  layout = 'bipartite') +
#           geom_edge_diagonal(
#             strength = 1,
#             aes(color = class),
#             edge_width = 0.5,
#             alpha = 1,
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class,
#                 size = Count),
#             shape = 21,
#             alpha = 1,
#             show.legend = TRUE
#           ) +
#           geom_node_text(
#             aes(
#               x = x * 1.03,
#               y = y * 1.03,
#               hjust = ifelse(class == "Pathway", "outward", 'inward'),
#               angle = -((-node_angle(x, y) + 90) %% 180) + 90,
#               label = annotation
#             ),
#             size = text_size,
#             show.legend = FALSE
#           ) +
#           scale_fill_manual(values = node_color) +
#           ggraph::scale_edge_color_manual(values = edge_color) +
#           scale_size_continuous(range = c(1, 8)) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "right",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#
#       } else{
#         coords$y[coords$class == "Functional_module"] <- 3
#         coords$y[coords$class == "Module"] <- 2
#         coords$y[coords$class == "Pathway"] <- 1
#         coords$y[coords$class == "Molecule"] <- 0
#
#         if (arrange_position) {
#           coords <-
#             arrange_coords(coords = coords,
#                            ratio = position_ratio)
#         }
#
#         my_graph <-
#           ggraph::create_layout(
#             graph = g,
#             layout = "manual",
#             x = coords$x,
#             y = coords$y
#             # node.position = coords
#           )
#
#         plot <-
#           ggraph(my_graph,
#                  layout = 'bipartite') +
#           geom_edge_diagonal(
#             strength = 1,
#             aes(color = class),
#             edge_width = 0.5,
#             alpha = 1,
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class,
#                 size = Count),
#             shape = 21,
#             alpha = 1,
#             show.legend = TRUE
#           ) +
#           geom_node_text(
#             aes(x = x,
#                 y = y,
#                 label = annotation),
#             hjust = 1,
#             angle = 90,
#             size = text_size,
#             show.legend = FALSE
#           ) +
#           scale_fill_manual(values = node_color) +
#           ggraph::scale_edge_color_manual(values = edge_color) +
#           scale_size_continuous(range = c(1, 8)) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "right",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#
#       }
#     } else{
#       total_graph <-
#         tidygraph::tbl_graph(nodes = node_data,
#                              edges = edge_data,
#                              directed = FALSE)
#
#       g <- total_graph
#       V(g)$type <- bipartite_mapping(g)$type
#       coords <-
#         ggraph::create_layout(g, layout = "bipartite")
#       coords$index = 1:nrow(coords)
#
#       if (circle) {
#         coords$y[coords$class == "Functional_module"] <- 0
#         coords$y[coords$class == "Module"] <- 1
#         coords$y[coords$class == "Pathway"] <- 2
#         coords$y[coords$class == "Molecule"] <- 3
#
#         if (arrange_position) {
#           coords <-
#             arrange_coords(coords = coords,
#                            ratio = position_ratio)
#         }
#
#         coords <-
#           coords %>%
#           dplyr::select(x, y) %>%
#           dplyr::mutate(
#             theta = x / (max(x) + 1) * 2 * pi,
#             r = y + 1,
#             x = r * cos(theta),
#             y = r * sin(theta)
#           )
#
#         my_graph <-
#           ggraph::create_layout(
#             graph = g,
#             layout = "manual",
#             x = coords$x,
#             y = coords$y
#             # node.position = coords
#           )
#
#         plot <-
#           ggraph(my_graph,
#                  layout = 'bipartite') +
#           geom_edge_diagonal(
#             strength = 1,
#             aes(color = class),
#             edge_width = 0.5,
#             alpha = 1,
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class,
#                 size = Count),
#             shape = 21,
#             alpha = 1,
#             show.legend = TRUE
#           ) +
#           geom_node_text(
#             aes(
#               x = x * 1.03,
#               y = y * 1.03,
#               hjust = ifelse(class == "Pathway", "outward", 'inward'),
#               angle = -((-node_angle(x, y) + 90) %% 180) + 90,
#               label = ifelse(class == "Molecule", NA, annotation)
#             ),
#             size = text_size,
#             show.legend = FALSE
#           ) +
#           scale_fill_manual(values = node_color) +
#           ggraph::scale_edge_color_manual(values = edge_color) +
#           scale_size_continuous(range = c(1, 8)) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "right",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#
#       } else{
#         coords$y[coords$class == "Functional_module"] <- 3
#         coords$y[coords$class == "Module"] <- 2
#         coords$y[coords$class == "Pathway"] <- 1
#         coords$y[coords$class == "Molecule"] <- 0
#
#         if (arrange_position) {
#           coords <-
#             arrange_coords(coords = coords,
#                            ratio = position_ratio)
#         }
#
#         my_graph <-
#           ggraph::create_layout(
#             graph = g,
#             layout = "manual",
#             x = coords$x,
#             y = coords$y
#             # node.position = coords
#           )
#
#         plot <-
#           ggraph(my_graph,
#                  layout = 'bipartite') +
#           geom_edge_diagonal(
#             strength = 1,
#             aes(color = class),
#             edge_width = 0.5,
#             alpha = 1,
#             show.legend = FALSE
#           ) +
#           geom_node_point(
#             aes(fill = class,
#                 size = Count),
#             shape = 21,
#             alpha = 1,
#             show.legend = TRUE
#           ) +
#           geom_node_text(
#             aes(
#               x = x,
#               y = y,
#               label = ifelse(class == "Molecule", NA, annotation)
#             ),
#             hjust = 1,
#             angle = 90,
#             size = text_size,
#             show.legend = FALSE
#           ) +
#           scale_fill_manual(values = node_color) +
#           ggraph::scale_edge_color_manual(values = edge_color) +
#           scale_size_continuous(range = c(1, 8)) +
#           ggraph::theme_graph() +
#           theme(
#             plot.background = element_rect(fill = "transparent", color = NA),
#             panel.background = element_rect(fill = "transparent", color = NA),
#             legend.position = "right",
#             legend.background = element_rect(fill = "transparent", color = NA)
#           )
#
#       }
#
#
#     }
#     plot
#   }
