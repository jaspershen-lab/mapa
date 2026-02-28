# setwd(r4projects::get_project_wd())
# metabolome_enrich
# gene_symbols <- T_input@variable_info$symbol
# test_gg_edge <- get_tf_target_edges(gene_symbols = gene_symbols, confidence_levels = "A")
#
# gene_symbols <- P_input@variable_info$symbol
# test_ppi_edge <- get_ppi_edges(gene_symbols = gene_symbols, input_directory = "demo_data/string_db/")
#
# protein_symbols <- P_input@variable_info$symbol
# metabolite_kegg <- M_input@variable_info$keggid
# reactome_dir <- "demo_data/reactome_db/"
# test_PM_reactome <- get_reactome_enzyme_metabolite_edges(protein_symbols = protein_symbols,
#                                                          metabolite_kegg = metabolite_kegg,
#                                                          reactome_dir = reactome_dir)
#
# test_PM_kegg <- get_kegg_enzyme_metabolite_edges(protein_symbols = protein_symbols, metabolite_kegg = metabolite_kegg)
#
# test_combo_res <- combine_enzyme_metabolite_edges(kegg_edges = test_PM_kegg, reactome_edges = test_PM_reactome)
# test <- get_molecule_pathway_edges(test$pathway_molecule_pairs)

get_tf_target_edges <- function(gene_symbols,
                                confidence_levels = "A") {
  data("dorothea_hs", package = "dorothea", envir = environment())

  dorothea_hs |>
    dplyr::filter(
      tf %in% gene_symbols,
      target %in% gene_symbols,
      confidence %in% confidence_levels
    ) |>
    dplyr::distinct(tf, target, mor, confidence) |>
    dplyr::arrange(confidence, tf, target) |>
    dplyr::mutate(edge_type = "TF_target")
}

get_ppi_edges <- function(gene_symbols,
                          taxon_id = 9606,
                          score_cutoff = 0.9,
                          input_directory = NA
                          ) {
  options(timeout = 600)

  string_db <- STRINGdb::STRINGdb$new(
    version = "12.0",
    species = taxon_id,
    score_threshold = score_cutoff*1000,
    input_directory = ifelse(is.na(input_directory), "", input_directory)
  )

  symbol_df <- data.frame(symbol = unique(gene_symbols),
                          stringsAsFactors = FALSE)
  mapped <- string_db$map(symbol_df, "symbol", removeUnmappedRows = TRUE)
  ppi <- string_db$get_interactions(mapped$STRING_id)

  ppi |>
    tibble::as_tibble() |>
    dplyr::distinct(from, to, combined_score) |>
    dplyr::mutate(edge_type = "PPI", combined_score = combined_score/1000) |>
    dplyr::rename(STRING_id = from) |>
    dplyr::left_join(mapped |> dplyr::rename(from = symbol), by = "STRING_id") |>
    dplyr::select(-STRING_id) |>
    dplyr::rename(STRING_id = to) |>
    dplyr::left_join(mapped |> dplyr::rename(to = symbol), by = "STRING_id") |>
    dplyr::select(-STRING_id) |>
    dplyr::filter(from %in% unique(gene_symbols) & to %in% unique(gene_symbols))
}


# Protein–Metabolite edges via Reactome
.map_symbol_to_uniprot <- function(gene_symbols) {
  ensembl  <- biomaRt::useMart("ensembl")
  datasets <- biomaRt::listDatasets(ensembl)
  dataset_name <- datasets$dataset[grepl("hsapiens_", datasets$dataset)]
  mart <- biomaRt::useMart("ensembl", dataset = dataset_name)

  sym_attr <- if ("hgnc_symbol" %in% biomaRt::listAttributes(mart)$name) {
    "hgnc_symbol"
  } else {
    "external_gene_name"
  }

  biomaRt::getBM(
    attributes = c(sym_attr, "uniprotswissprot", "uniprotsptrembl"),
    filters = sym_attr,
    values = unique(gene_symbols),
    mart = mart
  ) |>
    dplyr::mutate(
      symbol = .data[[sym_attr]],
      uniprotkb = dplyr::coalesce(uniprotswissprot, uniprotsptrembl)
    ) |>
    dplyr::filter(!is.na(uniprotkb), uniprotkb != "") |>
    dplyr::select(symbol, uniprotkb) |>
    dplyr::distinct() |>
    dplyr::filter(symbol %in% gene_symbols)
}


# Map KEGG compound IDs to ChEBI IDs
.map_kegg_to_chebi <- function(kegg_ids,
                               chunk_size = 100,
                               sleep_sec  = 0.1) {
  .chunk <- function(x, n) split(x, ceiling(seq_along(x) / n))

  kegg_unique <- unique(kegg_ids)
  chunks      <- .chunk(kegg_unique, chunk_size)
  all_maps    <- list()

  for (i in seq_along(chunks)) {
    q <- paste0("compound:", chunks[[i]])
    res <- tryCatch(
      KEGGREST::keggConv("chebi", q, querySize = chunk_size),
      error = function(e) {
        warning(sprintf("Chunk %d/%d failed: %s", i, length(chunks), conditionMessage(e)))
        character(0)
      }
    )
    if (length(res)) all_maps[[length(all_maps) + 1]] <- res
    if (sleep_sec > 0) Sys.sleep(sleep_sec)
  }

  m <- if (length(all_maps)) unlist(all_maps, use.names = TRUE) else character(0)

  out <- data.frame(
    keggid   = sub("^cpd:", "", names(m)),
    chebi_id = as.numeric(sub("^chebi:", "", unname(m))),
    stringsAsFactors = FALSE
  )

  missing <- setdiff(kegg_unique, out$keggid)
  if (length(missing)) {
    out <- rbind(out, data.frame(keggid = missing, chebi_id = NA_real_))
  }
  out[order(out$keggid), , drop = FALSE]
}


get_reactome_enzyme_metabolite_edges <- function(protein_symbols,
                                                 metabolite_kegg,
                                                 reactome_dir,
                                                 species_prefix = "HSA") {
  rxn_pat <- paste0("^R-", species_prefix, "-")

  # --- Protein: symbol -> UniProt -> catalyst reactions ---
  sym2uni <- .map_symbol_to_uniprot(protein_symbols)

  prot_role <- readr::read_delim(
    file.path(reactome_dir, "ProteinRoleReaction.txt"),
    col_names = c("protein_id", "role", "reaction_id"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(stringr::str_detect(reaction_id, rxn_pat))

  enz_rxn <- prot_role |>
    dplyr::filter(role == "catalystActivity") |>
    dplyr::distinct(protein_id, reaction_id) |>
    dplyr::rename(uniprotkb = protein_id) |>
    dplyr::inner_join(sym2uni, by = "uniprotkb")

  # --- Metabolite: KEGG -> ChEBI -> Reactome reactions ---
  kegg2chebi <- .map_kegg_to_chebi(metabolite_kegg)

  chebi2rxn <- readr::read_delim(
    file.path(reactome_dir, "ChEBI2Reactome_PE_Reactions.txt"),
    col_names = c("chebi_id", "reactome_entity_id", "reactome_entity_name",
                  "reaction_id", "URL", "reaction_name", "evidence_code", "species"),
    show_col_types = FALSE
  ) |>
    dplyr::select(chebi_id, reactome_entity_id, reactome_entity_name,
                  reaction_id, reaction_name) |>
    dplyr::filter(stringr::str_detect(reaction_id, rxn_pat)) |>
    dplyr::distinct(chebi_id, reaction_id, .keep_all = TRUE) |>
    dplyr::mutate(chebi_id = as.numeric(chebi_id))

  met_rxn <- kegg2chebi |>
    dplyr::left_join(chebi2rxn, by = "chebi_id") |>
    dplyr::filter(!is.na(reaction_id))

  # --- Join on shared reaction_id ---
  dplyr::inner_join(enz_rxn, met_rxn, by = "reaction_id",
                    relationship = "many-to-many") |>
    dplyr::distinct(uniprotkb, chebi_id, reaction_id, .keep_all = TRUE) |>
    dplyr::select(reaction_id, reaction_name,
                  symbol, uniprotkb,
                  keggid, chebi_id,
                  reactome_entity_id, reactome_entity_name) |>
    dplyr::mutate(edge_type = "enzyme_metabolite",
                  source = "Reactome")
}

# Helper: split vector into chunks and apply keggLink, return named character vector
.kegg_link_chunked <- function(target, keys, chunk_size, sleep_sec) {
  chunks <- split(keys, ceiling(seq_along(keys) / chunk_size))
  results <- purrr::map(chunks, function(chunk) {
    Sys.sleep(sleep_sec)
    tryCatch(
      KEGGREST::keggLink(target, chunk),
      error = function(e) {
        warning(sprintf("keggLink('%s') failed for chunk: %s", target, conditionMessage(e)))
        character(0)
      }
    )
  })
  unlist(unname(results))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

get_kegg_enzyme_metabolite_edges <- function(protein_symbols,
                                             metabolite_kegg,
                                             organism = "hsa",
                                             sleep_sec = 0.2,
                                             chunk_size = 50) {

  # 1) SYMBOL -> Entrez -> KEGG gene ID
  gene_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys    = unique(protein_symbols),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ) |>
    dplyr::filter(!is.na(ENTREZID)) |>
    dplyr::distinct(SYMBOL, ENTREZID) |>
    dplyr::rename(symbol = SYMBOL, entrezid = ENTREZID) |>
    dplyr::mutate(kegg_gene = paste0(organism, ":", entrezid))

  if (nrow(gene_map) == 0) return(tibble::tibble())

  # 2) Compound -> reaction IDs (chunked)
  cpd_keys <- paste0("cpd:", unique(metabolite_kegg))
  cpd2rxn  <- .kegg_link_chunked("reaction", cpd_keys, chunk_size, sleep_sec)

  cpd_rxn_tbl <- tibble::tibble(
    keggid      = sub("^cpd:", "", names(cpd2rxn)),
    reaction_id = sub("^rn:",  "", unname(cpd2rxn))
  ) |> dplyr::distinct()

  # 3a) Gene -> EC (chunked)
  gene_keys <- unique(gene_map$kegg_gene)
  gene2ec   <- .kegg_link_chunked("enzyme", gene_keys, chunk_size, sleep_sec)

  gene2ec_tbl <- tibble::tibble(
    kegg_gene = names(gene2ec),
    ec        = unname(gene2ec)
  )

  # 3b) EC -> reaction IDs (chunked)
  ec_ids  <- unique(unname(gene2ec))
  ec2rxn  <- if (length(ec_ids) > 0) {
    .kegg_link_chunked("reaction", ec_ids, chunk_size, sleep_sec)
  } else {
    character(0)
  }

  ec2rxn_tbl <- tibble::tibble(
    ec          = names(ec2rxn),
    reaction_id = sub("^rn:", "", unname(ec2rxn))
  )

  gene_rxn_tbl <- gene2ec_tbl |>
    dplyr::left_join(ec2rxn_tbl, by = "ec", relationship = "many-to-many") |>
    dplyr::filter(!is.na(reaction_id)) |>
    dplyr::left_join(gene_map, by = "kegg_gene")

  # 4) Intersect compound and gene reactions
  edges <- dplyr::inner_join(cpd_rxn_tbl, gene_rxn_tbl, by = "reaction_id", relationship = "many-to-many") |>
    dplyr::distinct(symbol, entrezid, ec, kegg_gene, keggid, reaction_id)

  if (nrow(edges) == 0) return(edges)

  # 5) Fetch reaction annotation (chunked, max 10 per keggGet)
  rxn_ids  <- unique(edges$reaction_id)
  rxn_chunks <- split(paste0("rn:", rxn_ids),
                      ceiling(seq_along(rxn_ids) / 10))

  rxn_list <- purrr::map(rxn_chunks, function(chunk) {
    Sys.sleep(sleep_sec)
    tryCatch(
      KEGGREST::keggGet(chunk),
      error = function(e) {
        warning(sprintf("keggGet failed for chunk: %s", conditionMessage(e)))
        list()
      }
    )
  }) |> purrr::flatten()

  rxn_anno <- purrr::map_dfr(rxn_list, function(x) {
    tibble::tibble(
      reaction_id   = unname(x$ENTRY)  %||% NA_character_,
      reaction_name = x$NAME[1] %||% NA_character_,
      definition    = x$DEFINITION %||% NA_character_,
      equation      = x$EQUATION %||% NA_character_
    )
  }) |> dplyr::distinct(reaction_id, .keep_all = TRUE)

  edges |>
    dplyr::left_join(rxn_anno, by = "reaction_id") |>
    dplyr::mutate(edge_type = "enzyme_metabolite",
                  source = "KEGG")
}

combine_enzyme_metabolite_edges <- function(kegg_edges, reactome_edges) {
  shared_cols <- c(
    "edge_type", "source", "reaction_id", "reaction_info",
    "symbol", "ec",
    "keggid", "reactome_entity_id"
  )

  kegg_norm <- kegg_edges |>
    tibble::as_tibble() |>
    dplyr::mutate(reactome_entity_id = NA_character_,
                  reaction_info = paste(reaction_name, definition, equation, sep = "\n")) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::select(dplyr::any_of(shared_cols))

  reactome_norm <- reactome_edges |>
    tibble::as_tibble() |>
    dplyr::mutate(ec = NA_character_,
                  reaction_info = reaction_name) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::select(dplyr::any_of(shared_cols))

  dplyr::bind_rows(kegg_norm, reactome_norm)
}

# Metabolite–Metabolite edges
get_metabolite_metabolite_edges <- function(
    metabolite_kegg,
    reactome_dir,
    organism = "hsa",
    species_prefix = "HSA",
    sleep_sec = 0.2,
    chunk_size = 50) {
  # From KEGG database
  cpd_keys <- paste0("cpd:", unique(metabolite_kegg))

  # compound → reaction
  cpd2rxn <- .kegg_link_chunked("reaction", cpd_keys, chunk_size, sleep_sec)

  cpd_rxn_tbl <- tibble::tibble(
    keggid = sub("^cpd:", "", names(cpd2rxn)),
    reaction_id = sub("^rn:",  "", unname(cpd2rxn))
  ) |>
    dplyr::filter(keggid %in% metabolite_kegg) |>
    dplyr::distinct()

  # reactions that involve ≥2 of our metabolites → self-join to get pairs
  kegg_pairs <- dplyr::inner_join(
    cpd_rxn_tbl, cpd_rxn_tbl,
    by = "reaction_id",
    relationship = "many-to-many",
    suffix = c("_from", "_to")
  ) |>
    dplyr::filter(keggid_from < keggid_to) |>
    dplyr::distinct(reaction_id, keggid_from, keggid_to)

  kegg_edges <- tibble::tibble()

  if (nrow(kegg_pairs) > 0) {
    rxn_ids <- unique(kegg_pairs$reaction_id)
    rxn_chunks <- split(paste0("rn:", rxn_ids),
                        ceiling(seq_along(rxn_ids) / 10))

    rxn_list <- purrr::map(rxn_chunks, function(chunk) {
      Sys.sleep(sleep_sec)
      tryCatch(KEGGREST::keggGet(chunk),
               error = function(e) {
                 warning(sprintf("keggGet failed: %s", conditionMessage(e)))
                 list()
               })
    }) |> purrr::flatten()

    rxn_anno <- purrr::map_dfr(rxn_list, function(x) {
      tibble::tibble(
        reaction_id   = unname(x$ENTRY)  %||% NA_character_,
        reaction_name = x$NAME[1]        %||% NA_character_,
        definition    = x$DEFINITION     %||% NA_character_,
        equation      = x$EQUATION       %||% NA_character_
      )
    }) |> dplyr::distinct(reaction_id, .keep_all = TRUE)

    kegg_edges <- kegg_pairs |>
      dplyr::left_join(rxn_anno, by = "reaction_id") |>
      dplyr::mutate(
        reaction_info = paste(reaction_name, definition, equation, sep = "\n"),
        edge_type = "metabolite_metabolite",
        source = "KEGG"
      ) |>
      dplyr::select(
        from = keggid_from,
        to = keggid_to,
        reaction_id,
        reaction_info,
        edge_type,
        source
      )
  }

  # From Reactome database
  rxn_pat <- paste0("^R-", species_prefix, "-")

  # KEGG → ChEBI (reuse the same private helper pattern inline)
  kegg2chebi <- .map_kegg_to_chebi(metabolite_kegg)

  chebi2rxn <- readr::read_delim(
    file.path(reactome_dir, "ChEBI2Reactome_PE_Reactions.txt"),
    col_names = c("chebi_id", "reactome_entity_id", "reactome_entity_name",
                  "reaction_id", "URL", "reaction_name", "evidence_code", "species"),
    show_col_types = FALSE
  ) |>
    dplyr::select(chebi_id, reactome_entity_id, reactome_entity_name,
                  reaction_id, reaction_name) |>
    dplyr::filter(stringr::str_detect(reaction_id, rxn_pat)) |>
    dplyr::distinct(chebi_id, reaction_id, .keep_all = TRUE) |>
    dplyr::mutate(chebi_id = as.numeric(chebi_id))

  met_rxn <- kegg2chebi |>
    dplyr::left_join(chebi2rxn, by = "chebi_id") |>
    dplyr::filter(!is.na(reaction_id))

  reactome_edges <- tibble::tibble()

  if (nrow(met_rxn) > 0) {
    reactome_pairs <- dplyr::inner_join(
      met_rxn, met_rxn,
      by = "reaction_id",
      relationship = "many-to-many",
      suffix = c("_from", "_to")
    ) |>
      dplyr::filter(keggid_from < keggid_to) |>
      dplyr::distinct(reaction_id, keggid_from, keggid_to,
                      reaction_name_from,
                      reactome_entity_id_from, reactome_entity_name_from,
                      reactome_entity_id_to, reactome_entity_name_to)

    if (nrow(reactome_pairs) > 0) {
      reactome_edges <- reactome_pairs |>
        dplyr::mutate(
          reaction_info = reaction_name_from,
          edge_type = "metabolite_metabolite",
          source = "Reactome"
        ) |>
        dplyr::select(
          from = keggid_from,
          to = keggid_to,
          reaction_id,
          reaction_info,
          edge_type,
          source
        )
    }
  }

  # Combine reaction information from Reactome and KEGG
  dplyr::bind_rows(kegg_edges, reactome_edges) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::distinct(from, to, reaction_id, source, .keep_all = TRUE)
}

# Molecule–Pathway edges (from enrichment)
get_molecule_pathway_edges <- function(pathway_molecule_pairs) {
  pathway_molecule_pairs |>
    dplyr::select(mapped_id, pathway_id) |>
    dplyr::mutate(edge_type = "molecule_pathway")
}
