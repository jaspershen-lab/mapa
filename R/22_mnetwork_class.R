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
# load("demo_data/demo_multi-omics/network_tables.rda")
# mnet_obj <- build_MNetwork(network_tables = network_tables,
#                            api_provider = "openai",
#                            text_embedding_model = "text-embedding-3-small",
#                            api_key = api_key)
# save(mnet_obj, file = "demo_data/demo_multi-omics/mnet_obj.rda")

#' @import methods
#' @importFrom Matrix sparseMatrix
NULL

#' S4 Class: multi_omics_functional_module
#'
#' Core data structure for a multi-omics functional module network. It
#' integrates four molecular-layer networks (TF–target, PPI,
#' metabolite–reaction, enzyme–metabolite), a pathway layer, and their
#' bipartite connections.
#'
#' @slot mol_nodes `data.frame`. Molecule node table with one row per unique
#'   molecule (gene or metabolite).
#'
#' @slot mol_layers `list`. Edge tables for each molecular layer,
#' as returned by [build_network_tables()]. Named elements:
#'   \describe{
#'     \item{`tf_target`}{TF–target edges (`from`, `to`, `mor`,
#'       `confidence`, `edge_type`). Gene symbols in both columns.}
#'     \item{`ppi`}{Protein–protein interaction edges (`from`, `to`,
#'       `combined_score`, `edge_type`). Gene symbols in both columns.}
#'     \item{`metabolite_reaction`}{Metabolite co-reaction edges (`from`,
#'       `to`, `reaction_id`, `reaction_info`, `edge_type`, `source`).
#'       KEGG compound IDs in both columns.}
#'     \item{`enzyme_metabolite`}{Enzyme–metabolite edges (`from`, `to`,
#'       `source`, `reaction_id`, `reaction_info`, `ec`,
#'       `reactome_entity_id`, `edge_type`). Gene symbol in `from`; KEGG
#'       compound ID in `to`.}
#'   }
#'
#' @slot mol_layers_edge_weight `list`. Processed, undirected edge tables (with
#'   columns `from`, `to`, `weight`) for each molecular layer, as returned by
#'   [build_mol_layers_weight()]. Named elements match those of
#'   `mol_layers`: `tf_target`, `ppi`, `metabolite_reaction`,
#'   `enzyme_metabolite`. All weights are binary (`1`) except `ppi` when
#'   `weight_ppi = TRUE`, in which case `combined_score` is used.
#'
#' @slot path_nodes `data.frame`. Pathway node table with one row per enriched
#'   pathway.
#'
#' @slot path_layer_edge_weight `data.frame`. Pathway–pathway similarity edge
#'   table derived from pairwise cosine similarities
#'   of biotext embeddings.
#'
#' @slot path_to_mol `data.frame`. Bipartite edge table linking pathways to
#'   molecules.
#'
#' @slot pathway_mol_edge_weight `data.frame`. Full pathway × molecule weight
#'   table with columns `pathway_id`, `mol_id`, and `weight`. Weight is
#'   an IDF-like score (`log((N_bg + 1) / (df_bg + 1))`) when the
#'   pathway–molecule edge is observed, and `0` otherwise.
#'   Pathways with fewer than `min_bg_genes_per_path` background genes are excluded.
#'
#' @slot layer_weights `numeric`. Named vector of normalised weights for each
#'   molecular layer; length equals `length(mol_layers)`. Equal weights
#'   (`1 / L`) are assigned when `NULL` is passed to [build_MNetwork()].
#'
#' @slot params `list`. Diffusion algorithm parameters:
#'   \describe{
#'     \item{`restart`}{Random-walk restart probability (default `0.7`).}
#'     \item{`interlayer`}{Inter-layer jump probability (default `0.5`).}
#'     \item{`gamma`}{Pathway-weight scaling factor (default `1.0`).}
#'     \item{`tol`}{Convergence tolerance (default `1e-6`).}
#'     \item{`max_iter`}{Maximum number of iterations (default `1000`).}
#'   }
#'
#' @slot index `list`. Fast-lookup index vectors for O(1) node position
#'   retrieval: `mol_index` (named integer vector mapping molecule `node_id` →
#'   row position in `mol_nodes`) and `path_index` (named integer vector
#'   mapping pathway `node_id` → row position in `path_nodes`).
#'
#' @exportClass multi_omics_functional_module

setClass(
  "multi_omics_functional_module",
  representation(
    mol_nodes = "data.frame",
    mol_layers = "list",
    mol_layers_edge_weight = "list",
    path_nodes = "data.frame",
    path_layer_edge_weight = "data.frame",
    path_to_mol = "data.frame",
    pathway_mol_edge_weight = "data.frame",
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

  # # Every mol_layer must be n × n
  # for (nm in names(object@mol_layers_edge_weight)) {
  #   mat <- object@mol_layers_edge_weight[[nm]]
  #   if (!inherits(mat, "dMatrix")) {
  #     errors <- c(errors, sprintf("mol_layers_edge_weight[['%s']] must be a dMatrix.", nm))
  #   } else if (!identical(dim(mat), c(n, n))) {
  #     errors <- c(errors,
  #                 sprintf("mol_layers_edge_weight[['%s']] has dim %s but expected %d x %d.",
  #                         nm, paste(dim(mat), collapse = "x"), n, n))
  #   }
  # }

  # # pathway_mol_edge_weight must be n × m
  # if (!identical(dim(object@pathway_mol_edge_weight), c(m, n))) {
  #   errors <- c(errors,
  #               sprintf("pathway_mol_edge_weight has dim %s but expected %d x %d.",
  #                       paste(dim(object@pathway_mol_edge_weight), collapse = "x"), m, n))
  # }
  #
  # # path_layer_edge_weight must be m × m
  # if (!identical(dim(object@path_layer_edge_weight), c(m, m))) {
  #   errors <- c(errors,
  #               sprintf("path_layer_edge_weight has dim %s but expected %d x %d.",
  #                       paste(dim(object@path_layer_edge_weight), collapse = "x"), m, m))
  # }

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
  cat(sprintf("  pathway_mol_edge_weight edges       : %d\n", Matrix::nnzero(object@pathway_mol_edge_weight)))
  cat(sprintf("  path_layer_edge_weight edges        : %d\n", Matrix::nnzero(object@path_layer_edge_weight)))
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
setGeneric("mol_layers<-", function(x, value) standardGeneric("mol_layers<-"))

#' @rdname slot-replace
#' @export
setGeneric("mol_layers_edge_weight<-", function(x, value) standardGeneric("mol_layers_edge_weight<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_nodes<-", function(x, value) standardGeneric("path_nodes<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_layer_edge_weight<-", function(x, value) standardGeneric("path_layer_edge_weight<-"))

#' @rdname slot-replace
#' @export
setGeneric("path_to_mol<-", function(x, value) standardGeneric("path_to_mol<-"))

#' @rdname slot-replace
#' @export
setGeneric("pathway_mol_edge_weight<-", function(x, value) standardGeneric("pathway_mol_edge_weight<-"))

#' @rdname slot-replace
#' @export
setGeneric("layer_weights<-", function(x, value) standardGeneric("layer_weights<-"))

#' @rdname slot-replace
setReplaceMethod("mol_nodes", "multi_omics_functional_module", function(x, value) {
  x@mol_nodes <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("mol_layers", "multi_omics_functional_module", function(x, value) {
  x@mol_layers <- value; validObject(x); x
})

setReplaceMethod("mol_layers_edge_weight", "multi_omics_functional_module", function(x, value) {
  x@mol_layers_edge_weight <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_nodes", "multi_omics_functional_module", function(x, value) {
  x@path_nodes <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_layer_edge_weight", "multi_omics_functional_module", function(x, value) {
  x@path_layer_edge_weight <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("path_to_mol", "multi_omics_functional_module", function(x, value) {
  x@path_to_mol <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("pathway_mol_edge_weight", "multi_omics_functional_module", function(x, value) {
  x@pathway_mol_edge_weight <- value; validObject(x); x
})

#' @rdname slot-replace
setReplaceMethod("layer_weights", "multi_omics_functional_module", function(x, value) {
  x@layer_weights <- value; validObject(x); x
})

#' Build a multi_omics_functional_module object
#'
#' Constructs a [multi_omics_functional_module] S4 object. Molecular-layer edge
#' weights, pathway–pathway similarity weights (via biotext embedding), and
#' pathway–molecule IDF weights are computed internally.
#'
#' @param network_tables `list`. Output of [build_network_tables()], containing:
#'   \describe{
#'     \item{`node_tables`}{A list with elements `mol_nodes` and
#'       `pathway_nodes` (each a `data.frame` with columns `node_id`,
#'       `node_type`, `node_info`).}
#'     \item{`edge_table`}{A named list of edge `data.frame`s:
#'       `tf_target`, `ppi`, `metabolite_reaction`, `enzyme_metabolite`,
#'       `pathway_mol`.}
#'   }
#' @param layer_weights `numeric` or `NULL`. Named numeric vector of length L
#'   (one value per molecular layer) supplying custom layer weights. When
#'   `NULL` (default), equal weights (`1 / L`) are assigned automatically.
#' @param params `list`. Named list of diffusion algorithm parameters that
#'   override the defaults. Recognised keys: `restart` (default `0.7`),
#'   `interlayer` (default `0.5`), `gamma` (default `1.0`), `tol` (default
#'   `1e-6`), `max_iter` (default `1000`).
#' @param api_provider `character(1)`. Embedding API provider for computing
#'   pathway–pathway text similarities. One of `"c"`, `"gemini"`, or
#'   `"siliconflow"`.
#' @param text_embedding_model `character(1)`. Name of the text embedding model
#'   to use (e.g. `"text-embedding-3-small"` for OpenAI).
#' @param api_key `character(1)`. API key for the chosen `api_provider`.
#'
#' @return A validated [multi_omics_functional_module] S4 object.
#'
#' @export
build_MNetwork <- function(network_tables,
                           api_provider = c("openai", "gemini", "siliconflow"),
                           text_embedding_model,
                           api_key,
                           layer_weights = NULL,
                           params = list()) {
  mol_nodes <- network_tables$node_tables$mol_nodes |> tibble::tibble()
  mol_layers <- network_tables$edge_table[names(network_tables$edge_table) != "pathway_mol"]
  mol_layers_edge_weight <- build_mol_layers_weight(mol_layers = mol_layers,
                                                    weight_ppi = FALSE)

  path_nodes <- network_tables$node_tables$pathway_nodes |> tibble::tibble()
  path_layer_edge_weight <- build_path_layer_weight(path_nodes = path_nodes,
                                                    api_provider = api_provider,
                                                    text_embedding_model = text_embedding_model,
                                                    api_key = api_key)

  path_to_mol <- network_tables$edge_table$pathway_mol |> tibble::tibble()
  pathway_mol_edge_weight <- build_pathway_mol_weight(network_tables,
                                                      min_bg_genes_per_path = 5)

  n <- nrow(mol_nodes)
  m <- nrow(path_nodes)

  # Build fast-lookup index vectors
  mol_index  <- stats::setNames(seq_len(n), mol_nodes$node_id)
  path_index <- stats::setNames(seq_len(m), path_nodes$node_id)

  # molecule layers edge weight
  L <- length(mol_layers)
  nonempty <- vapply(mol_layers_edge_weight, function(x) nrow(x) > 0, logical(1))
  L_active <- sum(nonempty)
  if (is.null(layer_weights)) {
    lw <- stats::setNames(
      ifelse(nonempty, if (L_active > 0) 1 / L_active else 0, 0),
      names(mol_layers)
    )
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
    mol_layers = mol_layers,
    mol_layers_edge_weight = mol_layers_edge_weight,
    path_nodes = path_nodes,
    path_layer_edge_weight = path_layer_edge_weight,
    path_to_mol = path_to_mol,
    pathway_mol_edge_weight = pathway_mol_edge_weight,
    layer_weights = lw,
    params = params,
    index = list(
      mol_index = mol_index,
      path_index = path_index
    )
  )
}
