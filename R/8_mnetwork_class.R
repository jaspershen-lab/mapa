# setwd(r4projects::get_project_wd())
# devtools::load_all()
# load("demo_data/demo_multi-omics/T_brain_up_enrich_pathway_res.rda")
# T_input <- enrich_pathway_res
# load("demo_data/demo_multi-omics/M_functional_module_res_u_and_d_met.rda")
# M_input <- functional_module_res
# load("demo_data/demo_multi-omics/P_up_enrich_pathway_res.rda")
# P_input <- enrich_pathway_res
# network_tables <- build_network_tables(transcriptome_enrich = T_input,
#                              proteome_enrich = P_input,
#                              metabolome_enrich = M_input,
#                              reactome_dir = "demo_data/reactome_db/",
#                              input_directory = "demo_data/string_db/",
#                              taxon_id = 9606,
#                              string_score_cutoff  = 0.9,
#                              tf_confidence_levels = "A")
# mnet_obj <- build_MNetwork(network_tables = network_tables)

#' @import methods
#' @importFrom Matrix sparseMatrix
NULL

#' S4 Class: multi_omics_functional_module
#'
#' Core data structure for a multi-omics functional module network, integrating
#' molecular-layer networks, a pathway layer, and their bipartite connections.
#'
#' @slot mol_nodes `data.frame`. Molecule node table; must contain columns
#'   `node_id` and `node_type`.
#' @slot path_nodes `data.frame`. Pathway node table; must contain columns
#'   `pathway_id` and `pathway_name`.
#' @slot path_to_mol `data.frame`. Edge table linking pathways to molecules.
#' @slot mol_layers `list`. Raw edge tables for each molecular layer, one
#'   data frame per layer.
#' @slot mol_layers_m `list`. Sparse adjacency matrices (`dMatrix`) for each
#'   molecular layer; each matrix has dimensions n × n.
#' @slot path_adj `dMatrix`. Pathway-layer similarity matrix; dimensions m × m.
#' @slot pathway_mol_m `dMatrix`. Bipartite matrix linking pathways to
#'   molecules; dimensions m × n (rows = pathways, columns = molecules).
#' @slot layer_weights `numeric`. Normalised weights for each molecular layer;
#'   length equals the number of layers L.
#' @slot params `list`. Diffusion algorithm parameters. Defaults:
#'   \describe{
#'     \item{restart}{Random-walk restart probability (default 0.7).}
#'     \item{interlayer}{Inter-layer jump probability (default 0.5).}
#'     \item{gamma}{Pathway-weight scaling factor (default 1.0).}
#'     \item{tol}{Convergence tolerance (default 1e-6).}
#'     \item{max_iter}{Maximum number of iterations (default 1000).}
#'   }
#' @slot index `list`. Fast-lookup index vectors: `mol_index` (molecule ID →
#'   row position) and `path_index` (pathway ID → row position).
#'
#' @seealso [build_MNetwork()] for constructing objects from standard input tables.
#'
#' @exportClass multi_omics_functional_module
setClass(
  "multi_omics_functional_module",
  representation(
    mol_nodes = "data.frame",
    path_nodes = "data.frame",
    path_to_mol = "data.frame",
    mol_layers = "list",
    mol_layers_m = "list",
    path_adj = "dMatrix",
    pathway_mol_m = "dMatrix",
    layer_weights = "numeric",
    params = "list",
    index = "list"
  )
)

# multi_omics_functional_module Validate slot dimensions and required column names.
setValidity("multi_omics_functional_module", function(object) {
  errors <- character()

  n <- nrow(object@mol_nodes)
  m <- nrow(object@path_nodes)
  L <- length(object@mol_layers)

  # Every mol_layer must be n × n
  for (nm in names(object@mol_layers_m)) {
    mat <- object@mol_layers_m[[nm]]
    if (!inherits(mat, "dMatrix")) {
      errors <- c(errors, sprintf("mol_layers_m[['%s']] must be a dMatrix.", nm))
    } else if (!identical(dim(mat), c(n, n))) {
      errors <- c(errors,
                  sprintf("mol_layers_m[['%s']] has dim %s but expected %d x %d.",
                          nm, paste(dim(mat), collapse = "x"), n, n))
    }
  }

  # pathway_mol_m must be n × m
  if (!identical(dim(object@pathway_mol_m), c(m, n))) {
    errors <- c(errors,
                sprintf("pathway_mol_m has dim %s but expected %d x %d.",
                        paste(dim(object@pathway_mol_m), collapse = "x"), m, n))
  }

  # path_adj must be m × m
  if (!identical(dim(object@path_adj), c(m, m))) {
    errors <- c(errors,
                sprintf("path_adj has dim %s but expected %d x %d.",
                        paste(dim(object@path_adj), collapse = "x"), m, m))
  }

  # layer_weights must be length L
  if (length(object@layer_weights) != L) {
    errors <- c(errors,
                sprintf("layer_weights length %d must equal length(mol_layers) = %d.",
                        length(object@layer_weights), L))
  }

  # mol_nodes must have required columns
  req_mol <- c("node_id", "node_type")
  miss_mol <- setdiff(req_mol, colnames(object@mol_nodes))
  if (length(miss_mol) > 0) {
    errors <- c(errors, paste("mol_nodes is missing columns:", paste(miss_mol, collapse = ", ")))
  }

  # path_nodes must have required columns
  req_path <- c("node_id", "node_type")
  miss_path <- setdiff(req_path, colnames(object@path_nodes))
  if (length(miss_path) > 0) {
    errors <- c(errors, paste("path_nodes is missing columns:", paste(miss_path, collapse = ", ")))
  }

  if (length(errors) == 0) TRUE else errors
})

# multi_omics_functional_module Initialise the object and run validity checks.
setMethod("initialize", "multi_omics_functional_module", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  validObject(.Object)
  .Object
})

# multi_omics_functional_module Print a concise summary of the object.
setMethod("show", "multi_omics_functional_module", function(object) {
  n <- nrow(object@mol_nodes)
  m <- nrow(object@path_nodes)
  L <- length(object@mol_layers)

  cat("── multi_omics_functional_module ─────────────────────────────────────────────────────\n")
  cat(sprintf("  Molecule nodes (n)    : %d  [%s]\n", n,
              paste(sort(table(object@mol_nodes$node_type), decreasing = TRUE),
                    names(sort(table(object@mol_nodes$node_type), decreasing = TRUE)),
                    sep = " ", collapse = ", ")))
  cat(sprintf("  Pathway nodes (m)     : %d\n", m))
  cat(sprintf("  Molecular layers (L)  : %d  [%s]\n", L,
              paste(names(object@mol_layers), collapse = ", ")))
  cat(sprintf("  Layer weights         : %s\n",
              paste(sprintf("%s=%.3f", names(object@layer_weights),
                            object@layer_weights), collapse = ", ")))
  cat(sprintf("  pathway_mol_m edges       : %d\n", Matrix::nnzero(object@pathway_mol_m)))
  cat(sprintf("  path_adj edges        : %d\n", Matrix::nnzero(object@path_adj)))
  cat("  Params                :",
      paste(names(object@params), unlist(object@params), sep = "=", collapse = ", "), "\n")
  cat("─────────────────────────────────────────────────────────────────\n")
})

#' Replace slots of a multi_omics_functional_module object
#'
#' Each replacement method validates the updated object via [validObject()]
#' before returning it.
#'
#' @param x A `multi_omics_functional_module` object.
#' @param value New value; must match the type of the corresponding slot.
#' @return The updated `multi_omics_functional_module` object.
#' @name slot-replace
NULL

#' @rdname slot-replace
#' @export
setGeneric("mol_nodes<-", function(x, value) standardGeneric("mol_nodes<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_nodes<-", function(x, value) standardGeneric("path_nodes<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_to_mol<-", function(x, value) standardGeneric("path_to_mol<-"))

#' @rdname slot-replace
#' @export
setGeneric("mol_layers<-", function(x, value) standardGeneric("mol_layers<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_adj<-", function(x, value) standardGeneric("path_adj<-"))

#' @rdname slot-replace
#' @export
setGeneric("pathway_mol_m<-", function(x, value) standardGeneric("pathway_mol_m<-"))

#' @rdname slot-replace
#' @export
setGeneric("layer_weights<-", function(x, value) standardGeneric("layer_weights<-"))

#' @rdname slot-replace
setReplaceMethod("mol_nodes", "multi_omics_functional_module", function(x, value) {
  x@mol_nodes <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_nodes", "multi_omics_functional_module", function(x, value) {
  x@path_nodes <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_to_mol", "multi_omics_functional_module", function(x, value) {
  x@path_to_mol <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("mol_layers", "multi_omics_functional_module", function(x, value) {
  x@mol_layers <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_adj", "multi_omics_functional_module", function(x, value) {
  x@path_adj <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("pathway_mol_m", "multi_omics_functional_module", function(x, value) {
  x@pathway_mol_m <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("layer_weights", "multi_omics_functional_module", function(x, value) {
  x@layer_weights <- value; validObject(x); x
})

#' Build a multi_omics_functional_module object
#'
#' Constructs a [multi_omics_functional_module] object.
#'
#' @param network_tables `list`. Standardised network tables.
#' @param mol_layers_m `list`. Edge weight matrix for each molecule layers.
#' @param path_adj `dMatrix` or `NULL`. m × m pathway adjacency matrix.
#'   A zero matrix is used when `NULL`.
#' @param layer_weights `numeric` or `NULL`. Weight vector of length L.
#'   Equal weights (1/L) are assigned when `NULL`.
#' @param params `list`. Named list overriding default diffusion parameters:
#'   `restart`, `interlayer`, `gamma`, `tol`, `max_iter`.
#'
#' @return A validated [multi_omics_functional_module] object.
#'
#' @export
build_MNetwork <- function(network_tables,
                           mol_layers_m = list(),
                           path_adj = NULL,
                           layer_weights = NULL,
                           params = list()) {

  mol_nodes <- network_tables$node_tables$mol_nodes
  path_nodes <- network_tables$node_tables$pathway_nodes
  path_to_mol <- network_tables$edge_table$pathway_mol

  n <- nrow(mol_nodes)
  m <- nrow(path_nodes)

  # Build fast-lookup index vectors
  mol_index  <- stats::setNames(seq_len(n), mol_nodes$node_id)
  path_index <- stats::setNames(seq_len(m), path_nodes$node_id)

  mol_layers <- network_tables$edge_table[names(network_tables$edge_table) != "pathway_mol"]

  pathway_mol_pair <- network_tables$edge_table$pathway_mol

  # pathway layer and molecule layers connections
  if (nrow(pathway_mol_pair) > 0) {
    bi_i <- path_index[pathway_mol_pair$from]
    bi_j <- mol_index[pathway_mol_pair$to]
    pathway_mol_m <- Matrix::sparseMatrix(
      i = bi_i,
      j = bi_j,
      x = rep(1, nrow(pathway_mol_pair)),
      dims = c(m, n),
      dimnames = list(names(path_index), names(mol_index))
    )
  } else {
    pathway_mol_m <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                          dims = c(m, n),
                                          dimnames = list(names(path_index), names(mol_index)))
  }

  pathway_mol_m <- methods::as(pathway_mol_m, "dMatrix")

  # pathway layer edge weight
  if (!is.null(path_adj)) {
    # Re-index if rownames present
    if (!is.null(rownames(path_adj))) {
      shared <- intersect(rownames(path_adj), names(path_index))
      if (length(shared) == m) {
        path_adj <- path_adj[names(path_index), names(path_index)]
      }
    }
    path_adj_mat <- methods::as(path_adj, "dMatrix")
  } else {
    # Zero m × m placeholder
    path_adj_mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                                         dims = c(m, m),
                                         dimnames = list(names(path_index), names(path_index)))
    path_adj_mat <- methods::as(path_adj_mat, "dMatrix")
  }

  # molecule layers edge weight
  L <- length(mol_layers)
  if (is.null(layer_weights)) {
    lw <- stats::setNames(rep(1 / L, L), names(mol_layers))
  } else {
    stopifnot(length(layer_weights) == L)
    lw <- stats::setNames(as.numeric(layer_weights), names(mol_layers))
  }

  # Default diffusion params
  default_params <- list(
    restart = 0.7,
    interlayer = 0.5,
    gamma = 1.0,
    tol = 1e-6,
    max_iter = 1000
  )
  params <- utils::modifyList(default_params, params)

  methods::new(
    "multi_omics_functional_module",
    mol_nodes = mol_nodes,
    path_nodes = path_nodes,
    path_to_mol = path_to_mol,
    mol_layers = mol_layers,
    mol_layers_m = mol_layers_m,
    path_adj = path_adj_mat,
    pathway_mol_m = pathway_mol_m,
    layer_weights = lw,
    params = params,
    index = list(
      mol_index = mol_index,
      path_index = path_index
    )
  )
}
