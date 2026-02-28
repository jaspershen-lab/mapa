# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/mnet_obj.rda")
# mnet <- mnet_obj
# multiplex_het <- build_MultiplexHet(mnet = mnet)
# save(multiplex_het, file = "demo_data/demo_multi-omics/multiplex_het.rda")

#' Convert a multi_omics_functional_module to a MultiplexHet object
#'
#' Transforms a [multi_omics_functional_module] S4 object into a
#' \code{MultiplexHet} object (from RandomWalkRestartMH) that is ready for
#' Random Walk with Restart on Multiplex Heterogeneous networks (RWR-MH).
#'
#' The resulting heterogeneous network is composed of two sub-networks:
#' \describe{
#'   \item{Multiplex 1 (molecules)}{A multi-layer molecular network built from
#'     each entry in \code{mol_layers_edge_weight}. Every layer is converted
#'     to an undirected, weighted igraph object.}
#'   \item{Multiplex 2 (pathways)}{A single-layer pathway similarity network
#'     built from \code{path_layer_edge_weight}. Edge weights are cosine
#'     similarities between pathway text embeddings.}
#'   \item{Bipartite relations}{Pathway–molecule associations derived from
#'     \code{pathway_mol_edge_weight}.}
#' }
#'
#' @param mnet A \code{multi_omics_functional_module} S4 object, as produced by
#'   [build_MNetwork()].
#' @param min_path_sim Numeric scalar. Minimum cosine similarity threshold for
#'   retaining pathway–pathway edges. Defaults to \code{0} (keep all positive
#'   similarities).
#' @param min_bipartite_weight Numeric scalar. Minimum IDF weight for retaining
#'   pathway–molecule bipartite edges. Defaults to \code{0} (keep all non-zero
#'   weights).
#'
#' @return A \code{MultiplexHet} object ready for use with
#'   [compute.transition.matrix()] and [Random.Walk.Restart.MultiplexHet()].
#'
#' @details
#' **Multiplex 1 – Molecular layers**
#'
#' Each element of \code{mol_layers_edge_weight} (columns \code{from},
#' \code{to}, \code{weight}) is converted to a weighted igraph via
#' [igraph::graph_from_data_frame()]. All molecular nodes present in
#' \code{mol_nodes} are added as vertices even if they have no edges in a
#' given layer, ensuring every layer shares the same node pool.
#'
#' **Multiplex 2 – Pathway layer**
#'
#' \code{path_layer_edge_weight} (columns \code{from}, \code{to},
#' \code{weight}) is converted to a single-layer weighted igraph. Edges with
#' cosine similarity below \code{min_path_sim} are discarded.
#'
#' **Bipartite relations**
#'
#' \code{pathway_mol_edge_weight} (columns \code{pathway_id}, \code{mol_id},
#' \code{weight}) is filtered to non-zero weights (and optionally by
#' \code{min_bipartite_weight}). The resulting data frame is passed as
#' \code{Nodes_relations} to [create.multiplexHet()], with column order
#' \code{[mol_id, pathway_id, weight]} (Multiplex 1 nodes first).
#'
#' @seealso [build_MNetwork()], [create.multiplex()],
#'   [create.multiplexHet()], [compute.transition.matrix()]
#'
#' @importFrom igraph graph_from_data_frame add_vertices V
#' @export
build_MultiplexHet <- function(mnet,
                               min_path_sim = 0.55,
                               min_bipartite_weight = 0) {

  ## Input validation
  if (!methods::is(mnet, "multi_omics_functional_module")) {
    stop("`mnet` must be a multi_omics_functional_module S4 object.")
  }

  ## Build Multiplex 1: molecular layers
  message("Building Multiplex 1 (molecular layers)...")

  mol_node_ids <- mnet@mol_nodes$node_id  # full universe of molecule nodes

  mol_igraph_list <- lapply(
    names(mnet@mol_layers_edge_weight),
    function(layer_name) {
      edges <- mnet@mol_layers_edge_weight[[layer_name]]

      # Create igraph from edge list
      g <- igraph::graph_from_data_frame(
        d = edges[, c("from", "to", "weight")],
        directed = FALSE,
        vertices = data.frame(name = mol_node_ids)
      )

      g
    }
  )
  names(mol_igraph_list) <- names(mnet@mol_layers_edge_weight)

  multiplex1 <- create.multiplex(mol_igraph_list)
  message(sprintf(
    "  Multiplex 1: %d nodes, %d layers [%s]",
    multiplex1$Number_of_Nodes_Multiplex,
    multiplex1$Number_of_Layers,
    paste(names(mol_igraph_list), collapse = ", ")
  ))

  ## Build Multiplex 2: pathway layer
  message("Building Multiplex 2 (pathway layer)...")

  path_edges <- mnet@path_layer_edge_weight

  # Filter by minimum similarity threshold
  path_edges <- path_edges[
    !is.na(path_edges$from) & !is.na(path_edges$to) &
      path_edges$from != path_edges$to &
      path_edges$weight > min_path_sim, ,
    drop = FALSE
  ]

  path_node_ids <- mnet@path_nodes$node_id

  path_igraph <- igraph::graph_from_data_frame(
    d = path_edges[, c("from", "to", "weight")],
    directed = FALSE,
    vertices = data.frame(name = path_node_ids)
  )

  multiplex2 <- create.multiplex(list(pathway = path_igraph))
  message(sprintf(
    "  Multiplex 2: %d pathway nodes",
    multiplex2$Number_of_Nodes_Multiplex
  ))

  ## 3. Build bipartite relations: mol <-> pathway
  message("Building bipartite relations (pathway-molecule)...")

  bip <- mnet@pathway_mol_edge_weight

  # Keep only non-zero edges above threshold
  bip <- bip[
    !is.na(bip$weight) & bip$weight > min_bipartite_weight, ,
    drop = FALSE
  ]

  # create.multiplexHet expects: col1 = Multiplex1 nodes (mol),
  #                              col2 = Multiplex2 nodes (pathway),
  #                              col3 = weight (optional)
  nodes_relations <- data.frame(
    mol = as.character(bip$mol_id),
    pathway = as.character(bip$pathway_id),
    weight  = bip$weight,
    stringsAsFactors = FALSE
  )

  # Retain only edges whose nodes actually appear in both multiplex pools
  valid_mol  <- multiplex1$Pool_of_Nodes
  valid_path <- multiplex2$Pool_of_Nodes

  nodes_relations <- nodes_relations[
    nodes_relations$mol %in% valid_mol &
      nodes_relations$pathway %in% valid_path, ,
    drop = FALSE
  ]

  if (nrow(nodes_relations) == 0) {
    stop(paste(
      "No valid bipartite edges remain after filtering.",
      "Check that mol_id / pathway_id match the node pools."
    ))
  }

  message(sprintf("  Bipartite edges retained: %d", nrow(nodes_relations)))

  ## 4. Assemble MultiplexHet
  message("Assembling MultiplexHet object...")

  multiplex_het <- create.multiplexHet(
    Multiplex_object_1 = multiplex1,
    Multiplex_object_2 = multiplex2,
    Nodes_relations    = nodes_relations
  )

  message("Done. MultiplexHet object created successfully.")
  return(multiplex_het)
}
