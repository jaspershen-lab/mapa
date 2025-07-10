# setwd(r4projects::get_project_wd())
# source("11_merge_modules.R")
# source("11_merge_pathways_bioembedsim.R")
###
# For gene overlap based
# load("demo_data/pregnancy_data/results/result_overlap/up_enriched_modules.rda")
# functional_modules <- get_functional_modules(object = merged_pathways,
#                                              sim.cutoff = 0.5,
#                                              cluster_method = c("louvain"))
####
# For biotext embedding
# setwd("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/")
# load("openai_semantic_sim_matrix.rda")
# biotext_functional_modules <- get_functional_modules(object = openai_semantic_sim_matrix,
#                                                      sim.cutoff = 0.9,
#                                                      cluster_method = "louvain")
# save(biotext_functional_modules, file = "biotext_functional_modules.rda")

#' Get Functional Modules from Pathway Similarity and Pathway Enrichment Results
#'
#' @description
#' A generic function to identify functional modules from pathway enrichment results.
#' This function supports two types of input objects:
#' \itemize{
#'   \item Traditional pathway enrichment objects (functional_module class)
#'   \item Biotext embedding similarity objects (list with enriched_pathway and sim_matrix)
#' }
#' The function clusters related pathways into functional modules using similarity-based
#' clustering methods.
#'
#' @param object An object can be either:
#'   \itemize{
#'     \item A \code{functional_module} class object from traditional pathway similarity calculation and pathway enrichment analysis
#'     \item A \code{list} object with components "enriched_pathway" and "sim_matrix" from biotext embedding pathway similarity calculation and pathway enrichment analysis
#'   }
#' @param ... Additional arguments passed to specific methods.
#'
#' @return The updated object with functional modules added to the "merged_module" slot.
#'
#'
#' The function supports multiple clustering methods including binary cut, Girvan-Newman
#' community detection, and hierarchical clustering to group related pathways into
#' functional modules.
#'
#' @examples
#' \dontrun{
#' # For traditional pathway enrichment results
#' enriched_modules <- get_functional_modules(
#'   object = my_functional_module_object,
#'   sim.cutoff = 0.5,
#'   cluster_method = "louvain"
#' )
#'
#' # For biotext embedding results
#' enriched_modules <- get_functional_modules(
#'   object = biotext_results,
#'   sim.cutoff = 0.6,
#'   cluster_method = "louvain"
#' )
#' }
#'
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export

get_functional_modules <- function(object, ...) {
  UseMethod("get_functional_modules")
}


#' Get Functional Modules from Traditional Pathway Similarity Calculation
#'
#' Identifies functional modules by calculating Jaccard similarity across different databases
#' and clustering pathways using the specified method.
#'
#' @param object A `functional_module` class object processed by `merge_pathways()` function.
#' @param sim.cutoff Numeric, similarity cutoff for clustering (default: 0.5).
#' @param measure_method Character, similarity measure. Currently supports "jaccard" (default).
#' @param cluster_method Character, clustering method options:  "louvain" (default), "hierarchical",
#'   "binary_cut", "walktrap", "infomap", "edge_betweenness", "fast_greedy", "label_prop",
#'   "leading_eigen", "optimal".
#' @param hclust.method Character, agglomeration method for hierarchical clustering including:
#'   "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#'   Only used when `cluster_method = "hierarchical"`.
#' @param path Character, directory path to save results (default: "result").
#' @param save_to_local Logical, save results to local files (default: FALSE).
#' @param ... Additional arguments (currently unused).
#'
#' @return Updated object with functional modules in the "merged_module" slot.
#'
#' @examples
#' \dontrun{
#' result <- get_functional_modules(
#'   object = my_enrichment_object,
#'   sim.cutoff = 0.5,
#'   cluster_method = "girvan newman"
#' )
#' }
#'
#' @method get_functional_modules functional_module
#' @export

get_functional_modules.functional_module <- function(object,
                                                     sim.cutoff = 0.5,
                                                     measure_method = c("jaccard"),
                                                     cluster_method = "louvain",
                                                     hclust.method = NULL,
                                                     path = "result",
                                                     save_to_local = FALSE,
                                                     ...) {

  message("Get functional modules from traditional pathway similarity calculation results ...")

  # Call your existing merge_modules function
  merge_modules(
    object = object,
    sim.cutoff = sim.cutoff,
    measure_method = measure_method,
    cluster_method = cluster_method,
    hclust.method = hclust.method,
    path = path,
    save_to_local = save_to_local
  )
}


#' Get Functional Modules from Biotext Embedding Similarity
#'
#' Identifies functional modules using pre-computed semantic similarity scores
#' from biotext embeddings.
#'
#' @param object A \code{list} containing "enriched_pathway" and "sim_matrix" components.
#' @param sim.cutoff Numeric, similarity cutoff for clustering (default: 0.5).
#' @param cluster_method Character, clustering method options:  "louvain" (default), "hierarchical",
#'   "binary_cut", "walktrap", "infomap", "edge_betweenness", "fast_greedy", "label_prop",
#'   "leading_eigen", "optimal".
#' @param hclust.method Character, agglomeration method for hierarchical clustering including:
#'   "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#'   Only used when `cluster_method = "hierarchical"`.
#' @param save_to_local Logical, save results to local files (default: FALSE).
#' @param path Character, directory path to save results (default: "result").
#' @param ... Additional arguments (currently unused).
#'
#' @return Updated object with functional modules.
#'
#' @examples
#' \dontrun{
#' biotext_results <- list(
#'   enriched_pathway = my_pathway_results,
#'   sim_matrix = embedding_similarity_matrix
#' )
#' result <- get_functional_modules(biotext_results, sim.cutoff = 0.6)
#' }
#'
#' @method get_functional_modules list
#' @export

get_functional_modules.list <- function(object,
                                        sim.cutoff = 0.5,
                                        cluster_method = "louvain",
                                        hclust.method = NULL,
                                        save_to_local = FALSE,
                                        path = "result",
                                        ...) {

  # Check if this is a bioembedding similarity object
  if (all(c("enriched_pathway", "sim_matrix") %in% names(object))) {
    message("Get functional modules from Biotext embedding results ...")

    # Call your existing merge_pathways_bioembedsim function
    merge_pathways_bioembedsim(
      object = object,
      sim.cutoff = sim.cutoff,
      cluster_method = cluster_method,
      hclust.method = hclust.method,
      save_to_local = save_to_local,
      path = path
    )
  } else {
    stop("List object must contain 'enriched_pathway' and 'sim_matrix' components.")
  }
}

