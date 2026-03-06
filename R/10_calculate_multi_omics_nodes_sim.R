# setwd(r4projects::get_project_wd())
# load("demo_data/demo_multi-omics/mnet_obj.rda")

#' Compute Multi-Omics Similarity Matrix
#'
#' Runs the full pipeline from a \code{multi_omics_functional_module} object to a
#' cosine-similarity matrix: builds a MultiplexHet object, computes RWR
#' diffusion profiles, and calculates pairwise cosine similarities.
#'
#' @param mnet_obj A \code{multi_omics_functional_module} S4 object produced by
#'   [build_MNetwork()].
#' @param min_path_sim Numeric. Minimum cosine similarity for pathway-pathway
#'   edges. Default \code{0.55}.
#' @param min_bipartite_weight Numeric. Minimum IDF weight for pathway-molecule
#'   edges. Default \code{0}.
#' @param TransitionMatrix Optional pre-computed transition matrix. Default
#'   \code{NULL}.
#' @param r Numeric. Restart probability for RWR. Default \code{0.5}.
#' @param eta Numeric. Inter-layer jumping probability. Default \code{0.5}.
#' @param lambda Numeric. Multiplex-to-bipartite weight. Default \code{0.2}.
#' @param delta1 Numeric. Layer weight in multiplex 1. Default \code{0.5}.
#' @param delta2 Numeric. Layer weight in multiplex 2. Default \code{0.5}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A sparse cosine-similarity matrix (rows and columns = node IDs).
#'
#' @export
get_multi_omics_sim <- function(
    mnet_obj,
    # build_MultiplexHet args
    min_path_sim = 0.55,
    min_bipartite_weight = 0,
    # compute_diffusion_profiles args
    TransitionMatrix = NULL,
    r = 0.5, eta = 0.5, lambda = 0.2, delta1 = 0.5, delta2 = 0.5,
    verbose = TRUE
) {
  # Step 1: mnet_obj → MultiplexHet
  if (verbose) message("Building MultiplexHet object ...")
  multiplex_het <- build_MultiplexHet(
    mnet = mnet_obj,
    min_path_sim = min_path_sim,
    min_bipartite_weight = min_bipartite_weight
  )

  # Step 2: MultiplexHet → diffusion profile matrix
  if (verbose) message("\n Computing diffusion profiles ...")
  diffusion_profile <- compute_diffusion_profiles(
    MultiplexHet_Object = multiplex_het,
    TransitionMatrix = TransitionMatrix,
    r = r, eta = eta, lambda = lambda, delta1 = delta1, delta2 = delta2,
    verbose = verbose
  )

  # Step 3: diffusion profile → cosine similarity matrix
  if (verbose) message("\n Calculating cosine similarity between the diffusion profiles...")
  sim_matrix <- .calculate_sim(
    m = diffusion_profile,
    zero_diag = TRUE,
    renorm_rows = TRUE
  )

  sim_matrix
}


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
#' @noRd
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

  multiplex1 <- create.multiplex.default(mol_igraph_list)
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
      path_edges$weight >= min_path_sim, ,
    drop = FALSE
  ]

  path_node_ids <- mnet@path_nodes$node_id

  path_igraph <- igraph::graph_from_data_frame(
    d = path_edges[, c("from", "to", "weight")],
    directed = FALSE,
    vertices = data.frame(name = path_node_ids)
  )

  multiplex2 <- create.multiplex.default(list(pathway = path_igraph))
  message(sprintf(
    "  Multiplex 2: %d pathway nodes",
    multiplex2$Number_of_Nodes_Multiplex
  ))

  ## 3. Build bipartite relations: mol <-> pathway
  message("Building bipartite relations (pathway-molecule)...")

  bip <- mnet@pathway_mol_edge_weight

  # Keep only non-zero edges above threshold
  bip <- bip[
    !is.na(bip$weight) & bip$weight >= min_bipartite_weight, ,
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

  multiplex_het <- create.multiplexHet.default(
    Multiplex_object_1 = multiplex1,
    Multiplex_object_2 = multiplex2,
    Nodes_relations = nodes_relations
  )

  message("Done. MultiplexHet object created successfully.")
  return(multiplex_het)
}

#' Compute Diffusion Profiles via Random Walk with Restart
#'
#' For each seed node in the MultiplexHet network, runs RWR and collects the
#' global score vector, assembling all results into a profile matrix.
#'
#' @param MultiplexHet_Object A \code{MultiplexHet} object.
#' @param TransitionMatrix Optional pre-computed transition matrix. Default
#'   \code{NULL} (computed internally).
#' @param r Numeric. Restart probability. Default \code{0.5}.
#' @param eta Numeric. Inter-layer jumping probability. Default \code{0.5}.
#' @param lambda Numeric. Multiplex-to-bipartite weight. Default \code{0.2}.
#' @param delta1 Numeric. Layer weight in multiplex 1. Default \code{0.5}.
#' @param delta2 Numeric. Layer weight in multiplex 2. Default \code{0.5}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A numeric matrix (seeds x nodes) of RWR scores.
#'
#' @noRd
compute_diffusion_profiles <- function(
    MultiplexHet_Object,
    TransitionMatrix = NULL,
    r = 0.5,
    eta = 0.5,
    lambda = 0.2,
    delta1 = 0.5,
    delta2 = 0.5,
    verbose = TRUE
) {

  # Compute transition matrix
  if (is.null(TransitionMatrix)) {
    if (verbose) message("[1/3] Computing transition matrix ...")
    TransitionMatrix <- compute.transition.matrix(
      MultiplexHet_Object,
      lambda = lambda,
      delta1 = delta1,
      delta2 = delta2
    )
  } else {
    if (verbose) message("[1/3] Using provided transition matrix.")
  }

  n_nan <- sum(is.nan(TransitionMatrix@x))
  n_inf <- sum(is.infinite(TransitionMatrix@x))
  if (n_nan > 0 || n_inf > 0) {
    if (verbose) message(sprintf("  Cleaning TransitionMatrix: %d NaN and %d Inf replaced with 0.", n_nan, n_inf))
    TransitionMatrix@x[is.nan(TransitionMatrix@x)] <- 0
    TransitionMatrix@x[is.infinite(TransitionMatrix@x)] <- 0
  }

  seeds_to_run  <- list(mol = MultiplexHet_Object$Multiplex1$Pool_of_Nodes,
                        pathway = MultiplexHet_Object$Multiplex2$Pool_of_Nodes)

  total_seeds <- sum(sapply(seeds_to_run, length))

  # Run RWR for each seed (node) -> extract global score vector

  run_single_seed <- function(seed, which_multiplex) {

    if (which_multiplex == "mol") {
      m1_seeds <- seed
      m2_seeds <- character(0)
    } else {
      m1_seeds <- character(0)
      m2_seeds <- seed
    }

    result <- tryCatch(
      Random.Walk.Restart.MultiplexHet.default(
        x = TransitionMatrix,
        MultiplexHet_Object = MultiplexHet_Object,
        Multiplex1_Seeds = m1_seeds,
        Multiplex2_Seeds = m2_seeds,
        r = r,
        eta = eta,
        DispResults = "Alphabetic"
      ),
      error = function(e) {
        warning(sprintf("RWR failed for seed '%s': %s", seed, e$message))
        return(NULL)
      }
    )
    if (is.null(result)) return(NULL)

    # global results
    global_df <- result$RWRMH_GlobalResults
    score_vec <- setNames(global_df$Score, global_df$NodeNames)
    return(score_vec)
  }

  profile_list <- list()
  counter <- 0

  for (net_type in names(seeds_to_run)) {
    seeds_vec <- seeds_to_run[[net_type]]
    for (seed in seeds_vec) {
      counter <- counter + 1
      if (verbose && (counter %% 50 == 0 || counter == total_seeds)) {
        message(sprintf("  ... %d / %d seeds done", counter, total_seeds))
      }
      profile_list[[seed]] <- run_single_seed(seed, net_type)
    }
  }

  failed <- names(which(sapply(profile_list, is.null)))
  if (length(failed) > 0) {
    warning(sprintf("%d seeds failed and were removed: %s",
                    length(failed), paste(failed, collapse = ", ")))
    profile_list <- profile_list[!names(profile_list) %in% failed]
  }

  # build profile matrix
  if (verbose) message("[3/3] Assembling profile matrix ...")

  all_col_names <- unique(unlist(lapply(profile_list, names)))

  profile_matrix <- do.call(rbind, lapply(profile_list, function(v) {
    out <- numeric(length(all_col_names))
    names(out) <- all_col_names
    out[names(v)] <- v
    out
  }))

  if (verbose) message(sprintf(
    "Done. Profile matrix: %d rows (seeds) x %d cols.",
    nrow(profile_matrix), ncol(profile_matrix)
  ))

  # return(list(
  #   profile_matrix   = profile_matrix,
  #   TransitionMatrix = TransitionMatrix
  # ))

  profile_matrix
}


#' Calculate Cosine Similarity Matrix
#'
#' @param m Numeric matrix of diffusion profiles (rows = seeds).
#' @param zero_diag Logical. Set diagonal to zero before normalisation. Default \code{TRUE}.
#' @param renorm_rows Logical. Row-normalise before computing cosine. Default \code{TRUE}.
#'
#' @return A sparse cosine-similarity matrix.
#'
#' @noRd
.calculate_sim <- function(m, zero_diag = TRUE, renorm_rows = TRUE) {

  node_ids <- rownames(m)

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Please install Matrix package.")
  }

  if (zero_diag) {
    diag(m) <- 0
  }

  if (renorm_rows) {
    rs <- Matrix::rowSums(m)
    inv_rs <- 1 / rs
    inv_rs[!is.finite(inv_rs)] <- 0
    m <- Matrix::Diagonal(x = inv_rs) %*% m
  }

  norms <- sqrt(Matrix::rowSums(m^2))
  inv_norms <- 1 / norms
  inv_norms[!is.finite(inv_norms)] <- 0

  m_norm <- Matrix::Diagonal(x = inv_norms) %*% m
  cosine_sim <- m_norm %*% Matrix::t(m_norm)

  if (!is.null(node_ids)) {
    dimnames(cosine_sim) <- list(node_ids, node_ids)
  }

  return(cosine_sim)
}

