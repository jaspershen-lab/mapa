.make_gene_omics_object <- function(symbols, source) {
  methods::new(
    "functional_module",
    variable_info = data.frame(
      symbol = symbols,
      ensembl = paste0("ENS", symbols),
      entrezid = seq_along(symbols),
      uniprot = paste0("UP", symbols),
      stringsAsFactors = FALSE
    ),
    process_info = stats::setNames(list(list(test = TRUE)), source)
  )
}

.make_metabolome_object <- function(kegg_ids) {
  methods::new(
    "functional_module",
    variable_info = data.frame(
      keggid = kegg_ids,
      cpd_name = paste0("met_", kegg_ids),
      diff_metric = seq_along(kegg_ids),
      stringsAsFactors = FALSE
    ),
    process_info = list(M = list(test = TRUE))
  )
}

.mock_network_edge_builders <- function() {
  empty_edges <- mapa:::.empty_raw_edge_tables()

  testthat::local_mocked_bindings(
    mapa_ensure_edge_databases = function(...) invisible(NULL),
    get_tf_target_edges = function(...) empty_edges$tf_target,
    get_ppi_edges = function(...) empty_edges$ppi,
    get_kegg_enzyme_metabolite_edges = function(...) tibble::tibble(),
    get_reactome_enzyme_metabolite_edges = function(...) tibble::tibble(),
    combine_enzyme_metabolite_edges = function(...) empty_edges$enzyme_metabolite,
    get_metabolite_metabolite_edges = function(...) empty_edges$metabolite_reaction,
    .package = "mapa",
    .env = parent.frame()
  )
}

test_that("build_network_tables requires at least two omics inputs", {
  transcriptome <- .make_gene_omics_object("A", "T")

  expect_error(
    build_network_tables(transcriptome_enrich = transcriptome),
    "At least two"
  )
  expect_error(
    build_network_tables(
      transcriptome_enrich = data.frame(),
      proteome_enrich = transcriptome
    ),
    "must be `functional_module` objects"
  )
})

test_that("build_network_tables accepts every pair of omics inputs", {
  .mock_network_edge_builders()

  transcriptome <- .make_gene_omics_object(c("A", "B"), "T")
  proteome <- .make_gene_omics_object(c("B", "C"), "P")
  metabolome <- .make_metabolome_object(c("C00001", "C00002"))

  combinations <- list(
    TP = list(transcriptome, proteome, NULL),
    TM = list(transcriptome, NULL, metabolome),
    PM = list(NULL, proteome, metabolome)
  )
  expected_node_types <- list(
    TP = c(gene = 3L),
    TM = c(gene = 2L, metabolite = 2L),
    PM = c(gene = 2L, metabolite = 2L)
  )

  for (combination in names(combinations)) {
    inputs <- combinations[[combination]]
    result <- suppressMessages(build_network_tables(
      transcriptome_enrich = inputs[[1]],
      proteome_enrich = inputs[[2]],
      metabolome_enrich = inputs[[3]]
    ))

    expect_equal(
      as.integer(table(result$node_tables$mol_nodes$node_type)),
      unname(expected_node_types[[combination]]),
      info = combination
    )
    expect_named(
      result$edge_table,
      c("tf_target", "ppi", "metabolite_reaction",
        "enzyme_metabolite", "pathway_mol"),
      info = combination
    )

    process_info <- attr(result, "process_info")
    expect_identical(
      unname(vapply(
        process_info[c("transcriptome_enrich", "proteome_enrich", "metabolome_enrich")],
        is.null,
        logical(1)
      )),
      vapply(inputs, is.null, logical(1)),
      info = combination
    )
  }
})

test_that("build_network_tables reports skipped layers accurately", {
  .mock_network_edge_builders()

  transcriptome <- .make_gene_omics_object(c("A", "B"), "T")
  proteome <- .make_gene_omics_object(c("B", "C"), "P")

  expect_message(
    build_network_tables(
      transcriptome_enrich = transcriptome,
      proteome_enrich = proteome
    ),
    "Step 5: Skipped enzyme-metabolite edges (no metabolite nodes).",
    fixed = TRUE
  )
  expect_message(
    build_network_tables(
      transcriptome_enrich = transcriptome,
      proteome_enrich = proteome
    ),
    "Step 6: Skipped metabolite-metabolite edges (no metabolite nodes).",
    fixed = TRUE
  )
})

test_that("missing omics sources retain correct node annotations", {
  transcriptome <- .make_gene_omics_object(c("A", "B"), "T")
  proteome <- .make_gene_omics_object(c("B", "C"), "P")
  metabolome <- .make_metabolome_object("C00001")
  empty_pathways <- tibble::tibble(
    pathway_id = character(), pathway_name = character(),
    BgRatio = character(), p_adjust = numeric()
  )

  tp <- mapa:::build_node_tables(
    transcriptome, proteome, NULL, empty_pathways
  )
  tm <- mapa:::build_node_tables(
    transcriptome, NULL, metabolome, empty_pathways
  )
  pm <- mapa:::build_node_tables(
    NULL, proteome, metabolome, empty_pathways
  )

  tp_sources <- stats::setNames(
    vapply(tp$mol_nodes$node_info, function(x) x$dt_src, character(1)),
    tp$mol_nodes$node_id
  )
  expect_identical(tp_sources[c("A", "B", "C")], c(A = "T", B = "P, T", C = "P"))
  expect_true(all(vapply(
    tm$mol_nodes$node_info[tm$mol_nodes$node_type == "gene"],
    function(x) identical(x$dt_src, "T"),
    logical(1)
  )))
  expect_true(all(vapply(
    pm$mol_nodes$node_info[pm$mol_nodes$node_type == "gene"],
    function(x) identical(x$dt_src, "P"),
    logical(1)
  )))
})
