# library(tidyverse)
# library(plyr)
# library(igraph)
#
# setwd(r4projects::get_project_wd())
# source("R/6-utils.R")
# source("R/8-functional_module_class.R")
# source("R/9-enrich_pathway.R")
# source("R/10-merge_pathways.R")
# source("R/11-merge_modules.R")
# source("R/12-data_visualization_plot_pathway_bar.R")
# source("R/13-data_visualization_plot_module_info.R")
# source("R/14-data_visualization_plot_similarity_network.R")
# source("R/15-data_visualization_plot_relationship_network.R")
# source("R/16-export_functional_module.R")
# source("R/30-dplyr_filter.R")
#
# setwd("demo_data/")
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("demo_data.rda")
#
# variable_info <-
#   demo_data %>%
#   massdataset::activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0) %>%
#   massdataset::extract_variable_info()
# #
# # enriched_pathways <-
# #   enrich_pathway(
# #     variable_info = variable_info,
# #     save_to_local = FALSE,
# #     path = "result",
# #     OrgDb = org.Hs.eg.db,
# #     organism = "hsa",
# #     database = c("go", "reactome", "kegg"),
# #     ont = "ALL",
# #     pvalueCutoff = 0.05,
# #     pAdjustMethod = "BH",
# #     qvalueCutoff = 0.2,
# #     minGSSize = 10,
# #     maxGSSize = 500,
# #     readable = FALSE,
# #     pool = FALSE
# #   )
# #
# # save(enriched_pathways, file = "enriched_pathways.rda")
# load("enriched_pathways.rda")
#
# # enriched_modules <-
# #   merge_pathways(
# #     object = enriched_pathways,
# #     p.adjust.cutoff.go = 0.05,
# #     p.adjust.cutoff.kegg = 0.05,
# #     p.adjust.cutoff.reactome = 0.05,
# #     count.cutoff.go = 5,
# #     count.cutoff.kegg = 5,
# #     count.cutoff.reactome = 5,
# #     sim.cutoff.go = 0.5,
# #     sim.cutoff.kegg = 0.5,
# #     sim.cutoff.reactome = 0.5,
# #     measure.method.go = "Wang",
# #     measure.method.kegg = "jaccard",
# #     measure.method.reactome = "jaccard",
# #     path = "result",
# #     save_to_local = FALSE
# #   )
# #
# # save(enriched_modules, file = "enriched_modules.rda")
# load("enriched_modules.rda")
#
# # enriched_functional_module <-
# #   merge_modules(
# #     object = enriched_modules,
# #     sim.cutoff = 0.5,
# #     measure_method = c("jaccard"),
# #     path = "result",
# #     save_to_local = FALSE
# #   )
# #
# # save(enriched_functional_module, file = "enriched_functional_module.rda")
# load("enriched_functional_module.rda")
#
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 20,
#   level = "pathway",
#   database = "go"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 20,
#   level = "pathway",
#   database = "kegg"
# )
#
# plot_pathway_bar(
#   object = enriched_functional_module,
#   top_n = 20,
#   level = "pathway",
#   database = "reactome"
# )
#
# plot_pathway_bar(object = enriched_functional_module,
#                  top_n = 20,
#                  level = "module")
#
# plot_pathway_bar(object = enriched_functional_module,
#                  top_n = 20,
#                  level = "functional_module")
#
#
# plot <-
#   plot_module_info(
#     object = enriched_functional_module,
#     level = "module",
#     database = "go",
#     module_id = "go_Module_3"
#   )
#
# plot_module_info(
#   object = enriched_functional_module,
#   level = "module",
#   database = "kegg",
#   module_id = "kegg_Module_15"
# )
#
#
# enriched_functional_module@merged_module$functional_module_result$module
#
# plot_module_info(object = enriched_functional_module,
#                  level = "functional_module",
#                  module_id = "Functional_module_17")
#
# export_module_info_plot(object = object, path = "result2")
#
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   database = "go",
#   degree_cutoff = 10
# )
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   database = "go",
#   degree_cutoff = 10,
#   module_id = "go_Module_10",
#   text_all = TRUE
# )
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   degree_cutoff = 0,
#   database = "go",
#   text_all = TRUE
# )
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   degree_cutoff = 0,
#   database = "kegg",
#   text_all = TRUE
# )
#
# plot_similarity_network(
#   object = enriched_functional_module,
#   level = "module",
#   degree_cutoff = 1,
#   database = "reactome",
#   text_all = TRUE
# )
#
#
# plot_relationship_network(
#   object = enriched_functional_module,
#   include_functional_modules = TRUE,
#   include_modules = FALSE,
#   include_pathways = TRUE,
#   include_molecules = FALSE,
#   functional_module_text = TRUE,
#   module_text = TRUE,
#   pathway_text = TRUE,
#   molecule_text = FALSE,
#   circular_plot = FALSE,
#   functional_module_arrange_position = TRUE,
#   module_arrange_position = TRUE,
#   pathway_arrange_position = TRUE,
#   molecule_arrange_position = TRUE,
#   functional_module_position_limits = c(0, 1),
#   module_position_limits = c(0, 1),
#   pathway_position_limits = c(0, 1),
#   molecule_position_limits = c(0, 1)
# )
