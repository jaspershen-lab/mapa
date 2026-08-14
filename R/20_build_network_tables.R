# setwd(r4projects::get_project_wd())
# rm(list = ls())
# load("demo_data/demo_multi-omics/multi-omics/results/transcriptomics_enrichment_2026-05-18.rda")
# T_input <- enrich_result
# load("demo_data/demo_multi-omics/multi-omics/results/metabolomics_enrichment_2026-05-18.rda")
# M_input <- enrich_result
# load("demo_data/demo_multi-omics/multi-omics/results/proteomics_enrichment_2026-05-18.rda")
# P_input <- enrich_result
# network_tables <- build_network_tables(transcriptome_enrich = T_input,
#                                        proteome_enrich = P_input,
#                                        metabolome_enrich = M_input,
#                                        taxon_id = 9606,
#                                        string_score_cutoff  = 0.9,
#                                        tf_confidence_levels = "A")
# attr(network_tables, "process_info")
# save(network_tables, file = "demo_data/demo_multi-omics/multi-omics/results/network_tables.rda")

#' Build Network Node and Edge Tables from Multi-Omics Enrichment Results
#'
#' @param transcriptome_enrich A \code{functional_module} object from transcriptome
#'   analysis, or \code{NULL} when transcriptomics data are unavailable.
#' @param proteome_enrich A \code{functional_module} object from proteome analysis,
#'   or \code{NULL} when proteomics data are unavailable.
#' @param metabolome_enrich A \code{functional_module} object from metabolome
#'   analysis, or \code{NULL} when metabolomics data are unavailable. At least two
#'   of the three omics inputs must be supplied.
#' @param taxon_id Integer. NCBI taxonomy ID for STRING. Default \code{9606}.
#' @param string_score_cutoff Numeric. Minimum STRING combined score. Default \code{0.9}.
#' @param tf_confidence_levels Character. DoRothEA confidence levels. Default \code{"A"}.
#' @return A list with \code{node_tables}, \code{edge_table}, and
#'   \code{pathway_version_harmonization_report}.
#' @export
build_network_tables <- function(transcriptome_enrich = NULL,
                                 proteome_enrich = NULL,
                                 metabolome_enrich = NULL,
                                 taxon_id = 9606,
                                 string_score_cutoff  = 0.9,
                                 tf_confidence_levels = "A") {
  omics_inputs <- list(
    transcriptome_enrich = transcriptome_enrich,
    proteome_enrich = proteome_enrich,
    metabolome_enrich = metabolome_enrich
  )
  provided <- !vapply(omics_inputs, is.null, logical(1))

  if (sum(provided) < 2L) {
    stop(
      "At least two of `transcriptome_enrich`, `proteome_enrich`, and ",
      "`metabolome_enrich` must be provided.",
      call. = FALSE
    )
  }

  invalid <- names(omics_inputs)[provided & !vapply(
    omics_inputs,
    function(x) is.null(x) || methods::is(x, "functional_module"),
    logical(1)
  )]
  if (length(invalid) > 0) {
    stop(
      "Non-NULL omics inputs must be `functional_module` objects: ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  go_ids <- unique(unlist(lapply(
    omics_inputs[provided],
    function(x) {
      if (is.null(x@enrichment_go_result)) character() else
        as.character(x@enrichment_go_result@result$ID)
    }
  ), use.names = FALSE))
  pathway_version_harmonization_report <-
    .empty_pathway_version_harmonization_report()
  if (length(go_ids) > 0) {
    message("Harmonizing GO term versions ...")
    pathway_version_harmonization_report <- .harmonize_go_ids(go_ids)
    if (!is.null(transcriptome_enrich) &&
        !is.null(transcriptome_enrich@enrichment_go_result)) {
      transcriptome_enrich@enrichment_go_result@result <-
        .apply_go_harmonization(
          transcriptome_enrich@enrichment_go_result@result,
          pathway_version_harmonization_report
        )
    }
    if (!is.null(proteome_enrich) &&
        !is.null(proteome_enrich@enrichment_go_result)) {
      proteome_enrich@enrichment_go_result@result <-
        .apply_go_harmonization(
          proteome_enrich@enrichment_go_result@result,
          pathway_version_harmonization_report
        )
    }
  }

  message("Ensuring edge databases are available ...")
  mapa_ensure_edge_databases()
  message("Step 1: Extracting enrichment results ...")
  enrich_out <- extract_enrichment_result(
    transcriptome_enrich = transcriptome_enrich,
    proteome_enrich = proteome_enrich
  )

  message("Step 2: Building node table ...")
  node_tables <- build_node_tables(
    transcriptome_enrich = transcriptome_enrich,
    proteome_enrich = proteome_enrich,
    metabolome_enrich = metabolome_enrich,
    enriched_pathway = enrich_out$pathway_nodes
  )

  # Collect DE symbols / KEGG IDs for edge construction
  all_gene_symbols  <- node_tables$mol_nodes |>
    dplyr::filter(node_type == "gene") |>
    dplyr::pull(node_id) |>
    stats::na.omit() |>
    unique()

  all_met_keggids <- node_tables$mol_nodes |>
    dplyr::filter(node_type == "metabolite") |>
    dplyr::pull(node_id) |>
    stats::na.omit() |>
    unique()

  empty_edges <- .empty_raw_edge_tables()

  tf_edges <- if (length(all_gene_symbols) > 0) {
    message("Step 3: Building TF-target edges ...")
    get_tf_target_edges(
      gene_symbols = all_gene_symbols,
      confidence_levels = tf_confidence_levels
    )
  } else {
    message("Step 3: Skipped TF-target edges (no gene/protein nodes).")
    empty_edges$tf_target
  }

  ppi_edges <- if (length(all_gene_symbols) > 0) {
    message("Step 4: Building PPI edges (STRING) ...")
    get_ppi_edges(
      gene_symbols = all_gene_symbols,
      taxon_id = taxon_id,
      score_cutoff = string_score_cutoff
    )
  } else {
    message("Step 4: Skipped PPI edges (no gene/protein nodes).")
    empty_edges$ppi
  }

  em_edges <- if (length(all_gene_symbols) > 0 && length(all_met_keggids) > 0) {
    message("Step 5: Building enzyme-metabolite edges ...")
    kegg_em_edges <- get_kegg_enzyme_metabolite_edges(
      protein_symbols = all_gene_symbols,
      metabolite_kegg = all_met_keggids,
      organism = "hsa"
    )
    reactome_em_edges <- get_reactome_enzyme_metabolite_edges(
      protein_symbols = all_gene_symbols,
      metabolite_kegg = all_met_keggids,
      species_prefix = "HSA"
    )
    combine_enzyme_metabolite_edges(
      kegg_edges = kegg_em_edges,
      reactome_edges = reactome_em_edges
    )
  } else {
    missing_node_types <- c(
      if (length(all_gene_symbols) == 0) "gene/protein",
      if (length(all_met_keggids) == 0) "metabolite"
    )
    message(
      "Step 5: Skipped enzyme-metabolite edges (no ",
      paste(missing_node_types, collapse = " or "),
      " nodes)."
    )
    empty_edges$enzyme_metabolite
  }

  mmrn_edges <- if (length(all_met_keggids) > 0) {
    message("Step 6: Building metabolite-metabolite edges ...")
    get_metabolite_metabolite_edges(
      metabolite_kegg = all_met_keggids,
      organism = "hsa",
      species_prefix = "HSA"
    )
  } else {
    message("Step 6: Skipped metabolite-metabolite edges (no metabolite nodes).")
    empty_edges$metabolite_reaction
  }

  message("Step 7: Building molecule-pathway edges ...")
  mp_edges <- get_molecule_pathway_edges(
    pathway_molecule_pairs = enrich_out$pathway_molecule_pairs
  )

  message("Step 8: Assembling edge table ...")
  edge_tables <- build_edge_tables(
    tf_target_edges = tf_edges,
    ppi_edges = ppi_edges,
    enzyme_metabolite_edges = em_edges,
    metabolite_metabolite_edges = mmrn_edges,
    molecule_pathway_edges  = mp_edges
  )

  result <- list(
    node_tables = node_tables,
    edge_table = edge_tables,
    pathway_version_harmonization_report =
      pathway_version_harmonization_report
  )

  attr(result, "process_info") <- list(
    transcriptome_enrich = if (is.null(transcriptome_enrich)) NULL else transcriptome_enrich@process_info,
    proteome_enrich      = if (is.null(proteome_enrich)) NULL else proteome_enrich@process_info,
    metabolome_enrich    = if (is.null(metabolome_enrich)) NULL else metabolome_enrich@process_info,
    build_network_tables = list(
      package_name  = "mapa",
      function_name = "build_network_tables()",
      parameter     = list(
        taxon_id             = taxon_id,
        string_score_cutoff  = string_score_cutoff,
        tf_confidence_levels = tf_confidence_levels
      ),
      data_sources = list(
        ppi = list(
          edge_type = "protein_protein_interaction",
          database  = "STRING v12.0",
          taxon_id  = taxon_id,
          files     = c(
            "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
            "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
            "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
          ),
          local_cache = .mapa_string_dir()
        ),
        tf_target = list(
          edge_type        = "TF_target_regulation",
          database         = "DoRothEA",
          r_package        = "dorothea",
          package_version  = as.character(utils::packageVersion("dorothea")),
          data_object      = "dorothea_hs",
          confidence_levels = tf_confidence_levels
        ),
        enzyme_metabolite = list(
          edge_type = "enzyme_metabolite",
          sources   = list(
            reactome = list(
              database    = "Reactome v96",
              files       = c(
                "https://download.reactome.org/96/ChEBI2Reactome_PE_Reactions.txt",
                "https://download.reactome.org/96/UniProt2Reactome_PE_Reactions.txt",
                "https://download.reactome.org/96/ProteinRoleReaction.txt"
              ),
              local_cache = .mapa_reactome_dir()
            ),
            kegg = list(
              database = "KEGG",
              access   = "realtime",
              r_package = "KEGGREST",
              package_version = as.character(utils::packageVersion("KEGGREST"))
            )
          )
        ),
        metabolite_reaction = list(
          edge_type = "metabolite_metabolite_reaction",
          sources   = list(
            reactome = list(
              database    = "Reactome v96",
              files       = c(
                "https://download.reactome.org/96/ChEBI2Reactome_PE_Reactions.txt"
              ),
              local_cache = .mapa_reactome_dir()
            ),
            kegg = list(
              database = "KEGG",
              access   = "realtime",
              r_package = "KEGGREST",
              package_version = as.character(utils::packageVersion("KEGGREST"))
            )
          )
        )
      ),
      time = Sys.time()
    )
  )

  result
}

#' Extract Enrichment Results from a Transcriptome/Proteome Object
#'
#' @param enrich_obj A \code{functional_module} object with enrichment result slots.
#' @return A tibble with columns \code{pathway_id}, \code{pathway_name}, \code{BgRatio},
#'   \code{p_adjust}, \code{mapped_id}.
#' @noRd
.extract_tp_enrichment <- function(enrich_obj) {
  empty_result <- function() {
    tibble::tibble(
      pathway_id = character(),
      pathway_name = character(),
      BgRatio = character(),
      p_adjust = numeric(),
      mapped_id = character()
    )
  }

  if (is.null(enrich_obj)) return(empty_result())

  extract_slot <- function(slot_result) {
    if (is.null(slot_result)) {
      return(empty_result())
    }

    enrich_res <-
      slot_result@result |>
      dplyr::select(
        pathway_id = ID,
        pathway_name = Description,
        BgRatio,
        p_adjust,
        mapped_id = geneID
      ) |>
      dplyr::filter(p_adjust < 0.05)

    if (nrow(enrich_res) == 0) {
      return(empty_result())
    } else {
      if (slot_result@keytype != "SYMBOL") {
        enrich_res |>
          dplyr::mutate(mapped_id = list(data.frame(entrezid = strsplit(mapped_id, split = "/")[[1]]) |>
                                           dplyr::left_join(enrich_obj@variable_info |>
                                                              dplyr::select(symbol, entrezid) |>
                                                              dplyr::filter(!is.na(symbol)),
                                                            by = "entrezid"))) |>
          dplyr::mutate(mapped_id = paste(mapped_id[[1]]$symbol, collapse = "/"))
      } else {
        enrich_res
      }
    }
  }

  dplyr::bind_rows(
    extract_slot(enrich_obj@enrichment_go_result),
    extract_slot(enrich_obj@enrichment_kegg_result),
    extract_slot(enrich_obj@enrichment_reactome_result)
  ) |>
    dplyr::filter(!is.na(pathway_id))
}

#' Extract and Combine Enrichment Results from Transcriptome and Proteome
#'
#' @param transcriptome_enrich A \code{functional_module} object from transcriptome analysis.
#' @param proteome_enrich A \code{functional_module} object from proteome analysis.
#' @return A list with \code{pathway_nodes} (unique pathways) and \code{pathway_molecule_pairs}
#'   (long-format pathway–molecule mapping).
#' @noRd
extract_enrichment_result <- function(transcriptome_enrich,
                                      proteome_enrich) {
  t_enrich <- .extract_tp_enrichment(transcriptome_enrich)
  p_enrich <- .extract_tp_enrichment(proteome_enrich)
  # m_enrich <- .extract_m_enrichment(metabolome_enrich)

  all_enrichment <- dplyr::bind_rows(t_enrich, p_enrich)

  pathway_nodes <- all_enrichment |>
    dplyr::select(-mapped_id) |>
    dplyr::distinct(pathway_id, .keep_all = TRUE)
  # Long-format: one row per pathway–molecule pair
  pathway_molecule_pairs <- all_enrichment |>
    dplyr::select(-c(pathway_name, BgRatio, p_adjust)) |>
    tidyr::separate_rows(mapped_id, sep = "\\s*[/;]\\s*") |>
    dplyr::filter(!is.na(mapped_id), mapped_id != "") |>
    dplyr::distinct(pathway_id, mapped_id, .keep_all = TRUE)

  list(
    pathway_nodes = pathway_nodes,
    pathway_molecule_pairs = pathway_molecule_pairs
  )
}

#' Build Gene, Metabolite, and Pathway Node Tables
#'
#' @param transcriptome_enrich A \code{functional_module} object from transcriptome analysis.
#' @param proteome_enrich A \code{functional_module} object from proteome analysis.
#' @param metabolome_enrich A \code{functional_module} object from metabolome analysis.
#' @param enriched_pathway A data frame of enriched pathways with \code{pathway_id} and metadata.
#' @return A list with \code{mol_nodes} and \code{pathway_nodes} tibbles.
#' @noRd
build_node_tables <- function(transcriptome_enrich,
                              proteome_enrich,
                              metabolome_enrich,
                              enriched_pathway) {
  # --- Gene nodes (union of transcriptome DE genes and proteome DE proteins) ---
  t_genes <- if (is.null(transcriptome_enrich)) {
    NULL
  } else {
    transcriptome_enrich@variable_info |>
      dplyr::select(symbol, ensembl, entrezid, uniprot) |>
      dplyr::mutate(node_type = "gene", dt_src = "T")
  }

  p_genes <- if (is.null(proteome_enrich)) {
    NULL
  } else {
    proteome_enrich@variable_info |>
      dplyr::select(symbol, ensembl, entrezid, uniprot) |>
      dplyr::mutate(node_type = "gene", dt_src = "P")
  }

  gene_nodes <- dplyr::bind_rows(t_genes, p_genes) |>
    dplyr::group_by(symbol) |>
    dplyr::summarise(
      dt_src = paste(sort(unique(dt_src)), collapse = ", "),
      ensembl = first(ensembl),
      entrezid = first(entrezid),
      uniprot = first(uniprot),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      node_id = symbol,
      node_type = "gene",
      node_info = purrr::pmap(
        list(symbol = symbol, ensembl = ensembl,
             entrezid = entrezid, uniprot = uniprot,
             dt_src = dt_src),
        function(...) list(...)
      )
    ) |>
    dplyr::select(-c(symbol, ensembl, entrezid, uniprot, dt_src))


  # --- Metabolite nodes ---
  met_var_info <- if (is.null(metabolome_enrich)) {
    tibble::tibble(
      keggid = character(),
      cpd_name = character(),
      diff_metric = numeric()
    )
  } else {
    metabolome_enrich@variable_info
  }
  if (!"diff_metric" %in% colnames(met_var_info)) {
    met_var_info$diff_metric <- NA_real_
  }
  met_nodes <- met_var_info |>
    dplyr::select(keggid, cpd_name, diff_metric) |>
    dplyr::filter(!is.na(keggid), keggid != "") |>
    dplyr::distinct() |>
    dplyr::mutate(node_id = keggid,
                  node_type = "metabolite",
                  node_info = purrr::pmap(
                    list(keggid = keggid,
                         cpd_name = cpd_name,
                         diff_metric = diff_metric),
                    function(...) list(...)
                  )) |>
    dplyr::select(-c(keggid, cpd_name, diff_metric))

  # --- Pathway nodes ---
  pathway_nodes <- enriched_pathway |>
    dplyr::distinct(pathway_id, .keep_all = TRUE) |>
    dplyr::rename(node_id = pathway_id) |>
    dplyr::mutate(node_type = "pathway",
                  node_info = purrr::pmap(
                    list(pathway_name = pathway_name,
                         p_adjust = p_adjust,
                         BgRatio = BgRatio),
                    function(...) list(...)
                  )) |>
    dplyr::select(-c(pathway_name, p_adjust, BgRatio))

  list(mol_nodes = rbind(gene_nodes, met_nodes) |> dplyr::select(node_id, node_type, node_info),
       pathway_nodes = pathway_nodes |> dplyr::select(node_id, node_type, node_info))
}

#' Create Typed Empty Edge Tables
#'
#' @return A named list matching the raw edge schemas consumed by
#'   \code{build_edge_tables()}.
#' @noRd
.empty_raw_edge_tables <- function() {
  list(
    tf_target = tibble::tibble(
      tf = character(), target = character(), mor = numeric(),
      confidence = character(), edge_type = character()
    ),
    ppi = tibble::tibble(
      from = character(), to = character(), combined_score = numeric(),
      edge_type = character()
    ),
    metabolite_reaction = tibble::tibble(
      from = character(), to = character(), reaction_id = character(),
      reaction_info = character(), edge_type = character(), source = character()
    ),
    enzyme_metabolite = tibble::tibble(
      edge_type = character(), source = character(), reaction_id = character(),
      reaction_info = character(), symbol = character(), ec = character(),
      keggid = character(), reactome_entity_id = character()
    )
  )
}

#' Assemble Normalized Edge Tables
#'
#' @param tf_target_edges Data frame from \code{get_tf_target_edges()}.
#' @param ppi_edges Data frame from \code{get_ppi_edges()}.
#' @param enzyme_metabolite_edges Data frame from \code{combine_enzyme_metabolite_edges()}.
#' @param metabolite_metabolite_edges Data frame from \code{get_metabolite_metabolite_edges()}.
#' @param molecule_pathway_edges Data frame from \code{get_molecule_pathway_edges()}.
#' @return A named list with normalized edge tibbles: \code{tf_target}, \code{ppi},
#'   \code{metabolite_reaction}, \code{enzyme_metabolite}, \code{pathway_mol}.
#' @noRd
build_edge_tables <- function(tf_target_edges,
                              ppi_edges,
                              enzyme_metabolite_edges,
                              metabolite_metabolite_edges,
                              molecule_pathway_edges) {
  # Normalize TF-target: from = tf, to = target (both are gene symbols)
  tf_tbl <- tf_target_edges |>
    dplyr::rename(from = tf, to = target) |>
    dplyr::select(edge_type, from, to, mor, confidence)

  # Normalize PPI
  ppi_tbl <- ppi_edges |>
    dplyr::select(edge_type, from, to, combined_score)

  # Normalize enzyme-metabolite: from = symbol (protein), to = keggid
  em_tbl <- enzyme_metabolite_edges |>
    dplyr::rename(from = symbol, to = keggid) |>
    dplyr::select(edge_type, from, to, source,
                  reaction_id, reaction_info,
                  ec, reactome_entity_id)

  # Normalize metabolite-metabolite
  mm_rnx <- metabolite_metabolite_edges |>
    dplyr::select(edge_type, everything())

  # Normalize molecule-pathway: from = pathway_id, to = mapped_id
  mp_tbl <- molecule_pathway_edges |>
    dplyr::rename(from = pathway_id, to = mapped_id) |>
    dplyr::select(edge_type, from, to)

  list(
    tf_target = tf_tbl,
    ppi = ppi_tbl,
    metabolite_reaction = mm_rnx,
    enzyme_metabolite = em_tbl,
    pathway_mol = mp_tbl
  )
}
