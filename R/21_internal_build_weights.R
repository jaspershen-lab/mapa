# setwd(r4projects::get_project_wd())
# devtools::load_all()
# testing
# load("demo_data/demo_multi-omics/mnet_obj.rda")
# mnetwork_obj <- mnet_obj
# nodes_relations <- build_nodes_relations(mnet_obj)
# mol_layer_weight <- build_mol_layer_weight(mnetwork_obj)
# test <- build_path_layer_weight(mnetwork_obj = mnetwork_obj,
#                                 api_provider = "openai",
#                                 text_embedding_model = "text-embedding-3-small",
#                                 api_key = api_key)


build_mol_layers_weight <- function(mol_layers,
                                    weight_ppi = FALSE) {
  # Standardise any edge table to (from, to, weight = 1)
  .binary_edges <- function(df) {
    df <- df[, c("from", "to"), drop = FALSE]
    df$from <- as.character(df$from)
    df$to   <- as.character(df$to)
    # remove self-loops
    df <- df[df$from != df$to, ]
    # deduplicate (treat as undirected: sort each pair)
    df <- unique(
      data.frame(
        from = pmin(df$from, df$to),
        to = pmax(df$from, df$to),
        weight = 1,
        stringsAsFactors = FALSE
      )
    )
    df
  }

  # Same but keep a numeric weight column from an existing column
  .weighted_edges <- function(df, weight_col) {
    df <- df[, c("from", "to", weight_col), drop = FALSE]
    colnames(df) <- c("from", "to", "weight")
    df$from <- as.character(df$from)
    df$to <- as.character(df$to)
    df$weight <- as.numeric(df$weight)
    df <- df[df$from != df$to, ]
    # for undirected symmetrisation: average weight of A->B and B->A if both exist
    key <- paste(pmin(df$from, df$to), pmax(df$from, df$to), sep = "__")
    agg <- tapply(df$weight, key, mean)
    pairs <- strsplit(names(agg), "__", fixed = TRUE)
    data.frame(
      from = vapply(pairs, `[[`, character(1), 1),
      to = vapply(pairs, `[[`, character(1), 2),
      weight = as.numeric(agg),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  ## tf_target layer
  edge_tables_list <- mol_layers
  tf_raw <- edge_tables_list$tf_target
  tf_edges <- .binary_edges(tf_raw)

  ## ppi layer
  ppi_raw <- edge_tables_list$ppi
  if (weight_ppi && "combined_score" %in% colnames(ppi_raw)) {
    ppi_edges <- .weighted_edges(ppi_raw, "combined_score")
  } else {
    ppi_edges <- .binary_edges(ppi_raw)
  }

  ## metabolite_reaction layer
  rxn_raw <- edge_tables_list$metabolite_reaction
  rxn_edges <- .binary_edges(rxn_raw)

  ## enzyme_metabolite layer
  # Expected columns: from, to  (gene symbol / KEGG compound ID)
  enz_raw <- edge_tables_list$enzyme_metabolite
  enz_edges <- .binary_edges(enz_raw)

  list(
    tf_target = tf_edges,
    ppi = ppi_edges,
    metabolite_reaction = rxn_edges,
    enzyme_metabolite = enz_edges
  )
}

build_path_layer_weight <- function(path_nodes,
                                    api_provider = c("openai", "gemini", "siliconflow"),
                                    text_embedding_model = NULL,
                                    api_key = NULL) {

    if (missing(api_provider)) {
      stop("api_provider is required.")
    }
    api_provider <- match.arg(api_provider)

    if (missing(text_embedding_model)) {
      stop("text_embedding_model is required.")
    }

    if (missing(api_key)) {
      stop("api_key is required.")
    }

    ## Collect text information from databases
    all_text_info <- list()
    pathway_nodes <- path_nodes |> tibble::tibble()
    pathway_nodes <- pathway_nodes |>
      dplyr::mutate(source = dplyr::case_when(
        startsWith(node_id, "GO") ~ "GO",
        grepl("^[a-z]{3}", node_id) ~ "KEGG",
        startsWith(node_id, "R-") ~ "Reactome"
      ))

    dbs <- unique(pathway_nodes$source)
    if ("GO" %in% dbs) {
      go_ids <- pathway_nodes |>
        dplyr::filter(source == "GO") |>
        dplyr::pull(node_id)

      message("Collecting pathway text information for GO terms from Gene Ontology database...")
      go_info <- get_go_info(go_ids)
      all_text_info <- c(all_text_info, go_info)
    }

    if ("KEGG" %in% dbs) {
      kegg_ids <- pathway_nodes |>
        dplyr::filter(source == "KEGG") |>
        dplyr::pull(node_id)

      message("Collecting pathway text information for KEGG pathways from KEGG database...")
      kegg_info <- get_kegg_pathway_info(kegg_ids)
      all_text_info <- c(all_text_info, kegg_info)
    }

    if ("Reactome" %in% dbs) {
      reactome_ids <- pathway_nodes |>
        dplyr::filter(source == "Reactome") |>
        dplyr::pull(node_id)

      message("Collecting pathway text information for Reactome pathways from Reactome database...")
      reactome_info <- get_reactome_pathway_info(reactome_ids)
      all_text_info <- c(all_text_info, reactome_info)
    }

    if (any(is.na(unlist(all_text_info)))) {
      all_text_info <- all_text_info[!is.na(all_text_info)]
    }

    all_combined_info <- combine_info(info = all_text_info)

    ## Get embedding matrix
    message("Getting pathway text embeddings ...")
    embedding_matrix <- get_embedding_matrix(text = all_combined_info,
                                             api_provider = api_provider,
                                             text_embedding_model = text_embedding_model,
                                             api_key = api_key)

    ## Calculate pairwise cosine similarity
    message("Calculating cosine similairty ...")
    sim_matrix <- calculate_cosine_sim(m = embedding_matrix)
    message("Biotext embedding and similarity calculation finished.\n")

    # Convert similarity matrix to edge table (from, to, weight)
    edge_df <- as.data.frame(as.table(sim_matrix))
    colnames(edge_df) <- c("from", "to", "weight")

    # Remove self-loops and duplicate edges (keep upper triangle only)
    edge_df <- edge_df[as.character(edge_df$from) < as.character(edge_df$to), ]

    edge_df
}

build_pathway_mol_weight <- function(network_tables,
                                     min_bg_genes_per_path = 5) {
  path_node_pair <- network_tables$edge_table$pathway_mol |> tibble::tibble()
  mol_nodes <- network_tables$node_tables$mol_nodes |> tibble::tibble()
  pathway_nodes <- network_tables$node_tables$pathway_nodes |> tibble::tibble()

  # extract BgRatio for each pathway
  bgratio_vec <- vapply(
    pathway_nodes$node_info,
    FUN = function(x) {
      if (is.list(x) && !is.null(x$BgRatio)) as.character(x$BgRatio) else NA_character_
    },
    FUN.VALUE = character(1)
  )

  pathway_bg <- dplyr::tibble(
    path_id = pathway_nodes$node_id,
    BgRatio = bgratio_vec
  ) |>
    dplyr::filter(!is.na(BgRatio))

  if (nrow(pathway_bg) == 0) stop("No BgRatio found in pathway_nodes$node_info.")

  # parse "a/b" -> df_bg = a, N_bg = b
  parts <- strsplit(pathway_bg$BgRatio, "/", fixed = TRUE)
  df_bg <- suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 1)))
  N_bg  <- suppressWarnings(as.integer(vapply(parts, `[[`, character(1), 2)))

  pathway_bg <- pathway_bg |>
    dplyr::mutate(df_bg = df_bg, N_bg = N_bg) |>
    dplyr::filter(!is.na(df_bg), !is.na(N_bg)) |>
    dplyr::filter(df_bg >= as.integer(min_bg_genes_per_path))

  if (nrow(pathway_bg) == 0) stop("No pathways left after filtering by min_bg_genes_per_path.")

  # compute IDF-like weight per pathway
  pathway_bg <- pathway_bg |>
    dplyr::mutate(idf_raw = base::log((N_bg + 1) / (df_bg + 1)))

  # if (scale_to_01) {
  #   max_idf <- max(pathway_bg$idf_raw, na.rm = TRUE)
  #   pathway_bg <- pathway_bg |>
  #     dplyr::mutate(weight_path = if (max_idf == 0) 0 else idf_raw / max_idf)
  # } else {
  #   pathway_bg <- pathway_bg |> dplyr::mutate(weight_path = idf_raw)
  # }

  pathway_bg <- pathway_bg |> dplyr::mutate(weight_path = idf_raw)

  pathway_weight <- pathway_bg |>
    dplyr::select(path_id, weight_path) |>
    dplyr::distinct(path_id, .keep_all = TRUE)

  # build full pathway x molecule grid
  all_path_ids <- pathway_weight$path_id
  all_mol_ids  <- mol_nodes$node_id

  full_grid <- tidyr::expand_grid(
    path_id = all_path_ids,
    mol_id  = all_mol_ids
  )

  # mark observed (annotated) edges
  observed_edges <- path_node_pair |>
    dplyr::rename(path_id = from, mol_id = to) |>
    dplyr::distinct(path_id, mol_id) |>
    dplyr::mutate(.observed = TRUE)

  # assign weights: IDF weight if edge exists, 0 otherwise
  nodes_relations <- full_grid |>
    dplyr::left_join(pathway_weight, by = "path_id") |>
    dplyr::left_join(observed_edges, by = c("path_id", "mol_id")) |>
    dplyr::mutate(
      weight = dplyr::if_else(!is.na(.observed), weight_path, 0)
    ) |>
    dplyr::select(pathway_id = path_id, mol_id, weight)

  tibble::tibble(nodes_relations)
}

