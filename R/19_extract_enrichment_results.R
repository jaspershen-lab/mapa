# setwd(r4projects::get_project_wd())
# transcriptome_enrich <- T_input
# proteome_enrich <- P_input
# metabolome_enrich <- M_input
# all_enrichment <- build_enrichment_tables(transcriptome_enrich, proteome_enrich, metabolome_enrich)

# Helper: extract GO / KEGG / Reactome enrichment from a mapa object for transcriptome and proteome
.extract_tp_enrichment <- function(enrich_obj) {
  extract_slot <- function(slot_result) {
    if (is.null(slot_result)) {
      return(tibble::tibble(
        pathway_id = NA_character_,
        pathway_name = NA_character_,
        p_adjust = NA_real_,
        mapped_id = NA_character_
      ))
    }

    enrich_res <-
      slot_result@result |>
      dplyr::select(
        pathway_id = ID,
        pathway_name = Description,
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

# Helper: extract HMDB / MetKEGG enrichment from a mapa object for matabolome
.extract_m_enrichment <- function(enrich_obj) {
  extract_slot <- function(slot_result) {
    if (is.null(slot_result)) {
      return(tibble::tibble(
        pathway_id   = NA_character_,
        pathway_name = NA_character_,
        p_adjust     = NA_real_,
        mapped_id    = NA_character_
      ))
    }
    slot_result@result |>
      dplyr::select(pathway_id, pathway_name, p_adjust, mapped_id) |>
      dplyr::filter(p_adjust < 0.05)
  }

  dplyr::bind_rows(
    extract_slot(enrich_obj@enrichment_hmdb_result),
    extract_slot(enrich_obj@enrichment_metkegg_result)
  ) |>
    dplyr::filter(!is.na(pathway_id))
}


# Build combined pathway enrichment table from T / P / M objects
build_enrichment_tables <- function(transcriptome_enrich,
                                    proteome_enrich,
                                    metabolome_enrich) {
  t_enrich <- .extract_tp_enrichment(transcriptome_enrich)
  p_enrich <- .extract_tp_enrichment(proteome_enrich)
  m_enrich <- .extract_m_enrichment(metabolome_enrich)

  all_enrichment <- dplyr::bind_rows(t_enrich, p_enrich, m_enrich)

  # Long-format: one row per pathway–molecule pair
  pathway_molecule_pairs <- all_enrichment |>
    tidyr::separate_rows(mapped_id, sep = "\\s*[/;]\\s*") |>
    dplyr::filter(!is.na(mapped_id), mapped_id != "") |>
    dplyr::distinct(pathway_id, mapped_id, .keep_all = TRUE)

  list(
    all_enrichment = all_enrichment,
    pathway_molecule_pairs = pathway_molecule_pairs
  )
}
