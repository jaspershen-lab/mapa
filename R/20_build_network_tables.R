# setwd(r4projects::get_project_wd())
# source("R/19_extract_enrichment_results.R")
# source("R/20_build_edges.R")
# library(mapa)
# load("demo_data/demo_multi-omics/T_brain_up_enrich_pathway_res.rda")
# T_input <- enrich_pathway_res
# load("demo_data/demo_multi-omics/M_functional_module_res_u_and_d_met.rda")
# M_input <- functional_module_res
# load("demo_data/demo_multi-omics/P_up_enrich_pathway_res.rda")
# P_input <- enrich_pathway_res
# network_tables <- build_network_tables(transcriptome_enrich = T_input,
#                                        proteome_enrich = P_input,
#                                        metabolome_enrich = M_input,
#                                        reactome_dir = "demo_data/reactome_db/",
#                                        input_directory = "demo_data/string_db/",
#                                        taxon_id = 9606,
#                                        string_score_cutoff  = 0.9,
#                                        tf_confidence_levels = "A")
# save(network_tables, file = "demo_data/demo_multi-omics/network_tables.rda")

build_network_tables <- function(transcriptome_enrich,
                                 proteome_enrich,
                                 metabolome_enrich,
                                 taxon_id = 9606,
                                 reactome_dir = NA,
                                 input_directory = NA,
                                 string_score_cutoff  = 0.9,
                                 tf_confidence_levels = "A") {
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

  message("Step 3: Building TF-target edges ...")
  tf_edges <- get_tf_target_edges(
    gene_symbols = all_gene_symbols,
    confidence_levels = tf_confidence_levels
  )

  message("Step 4: Building PPI edges (STRING) ...")
  ppi_edges <- get_ppi_edges(
    gene_symbols = all_gene_symbols,
    taxon_id = taxon_id,
    score_cutoff = string_score_cutoff,
    input_directory = input_directory
  )

  message("Step 5: Building enzyme-metabolite edges ...")
  kegg_em_edges <- get_kegg_enzyme_metabolite_edges(
    protein_symbols = all_gene_symbols,
    metabolite_kegg = all_met_keggids,
    organism = "hsa"
  )
  reactome_em_edges <- get_reactome_enzyme_metabolite_edges(
    protein_symbols = all_gene_symbols,
    metabolite_kegg = all_met_keggids,
    species_prefix = "HSA",
    reactome_dir = reactome_dir
  )
  em_edges <- combine_enzyme_metabolite_edges(
    kegg_edges = kegg_em_edges,
    reactome_edges = reactome_em_edges
  )

  message("Step 6: Building metabolite-metabolite edges ...")
  mmrn_edges <- get_metabolite_metabolite_edges(
    metabolite_kegg = all_met_keggids,
    reactome_dir = reactome_dir,
    organism = "hsa",
    species_prefix = "HSA"
  )

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

  list(
    node_tables = node_tables,
    edge_table = edge_tables
  )
}

# Helper: extract GO / KEGG / Reactome enrichment from a mapa object for transcriptome and proteome
.extract_tp_enrichment <- function(enrich_obj) {
  extract_slot <- function(slot_result) {
    if (is.null(slot_result)) {
      return(tibble::tibble(
        pathway_id = NA_character_,
        pathway_name = NA_character_,
        BgRatio = NA_character_,
        p_adjust = NA_real_,
        mapped_id = NA_character_
      ))
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

  dplyr::bind_rows(
    extract_slot(enrich_obj@enrichment_go_result),
    extract_slot(enrich_obj@enrichment_kegg_result),
    extract_slot(enrich_obj@enrichment_reactome_result)
  ) |>
    dplyr::filter(!is.na(pathway_id))
}

# Build combined pathway enrichment table from T / P / M objects
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

build_node_tables <- function(transcriptome_enrich,
                              proteome_enrich,
                              metabolome_enrich,
                              enriched_pathway) {
  # --- Gene nodes (union of transcriptome DE genes and proteome DE proteins) ---
  t_genes <- transcriptome_enrich@variable_info |>
    dplyr::select(symbol, ensembl, entrezid, uniprot) |>
    dplyr::mutate(node_type = "gene")

  p_genes <- proteome_enrich@variable_info |>
    dplyr::select(symbol, ensembl, entrezid, uniprot) |>
    dplyr::mutate(node_type = "gene")

  gene_nodes <- dplyr::bind_rows(t_genes, p_genes) |>
    dplyr::distinct(symbol, .keep_all = TRUE) |>
    dplyr::mutate(
      node_id = symbol,
      node_info = purrr::pmap(
        list(symbol = symbol, ensembl = ensembl,
             entrezid = entrezid, uniprot = uniprot),
        function(...) list(...)
      )
    ) |>
    dplyr::select(-c(symbol, ensembl, entrezid, uniprot))


  # --- Metabolite nodes ---
  met_nodes <- metabolome_enrich@variable_info |>
    dplyr::select(keggid) |>
    dplyr::filter(!is.na(keggid), keggid != "") |>
    dplyr::distinct() |>
    dplyr::mutate(node_type = "metabolite",
                  node_id   = keggid,
                  node_info = purrr::pmap(
                    list(keggid = keggid),
                    function(...) list(...)
                  )) |>
    dplyr::select(-keggid)


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
