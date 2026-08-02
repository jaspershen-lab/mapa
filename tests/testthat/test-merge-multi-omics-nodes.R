test_that("multi_omics_num counts T, P, and M layers", {
  module_nodes <- tibble::tibble(
    module = c(
      "t_only",
      "p_only",
      "tp_gene",
      "m_only",
      "t_m",
      "t_m",
      "t_p_m",
      "t_p_m",
      "pathway_only",
      "missing_source"
    ),
    node_type = c(
      "gene",
      "gene",
      "gene",
      "metabolite",
      "gene",
      "metabolite",
      "gene",
      "metabolite",
      "pathway",
      "gene"
    ),
    node_info = list(
      list(dt_src = "T"),
      list(dt_src = "P"),
      list(dt_src = "P, T"),
      list(keggid = "C00001"),
      list(dt_src = "T"),
      list(keggid = "C00002"),
      list(dt_src = "P, T"),
      list(keggid = "C00003"),
      list(pathway_name = "Example pathway"),
      list(symbol = "GENE1")
    )
  )

  result <- mapa:::.summarise_module_omics(module_nodes)
  result <- result[match(unique(module_nodes$module), result$module), ]

  expect_equal(
    result$multi_omics_num,
    c(1L, 1L, 2L, 1L, 2L, 3L, 0L, 0L)
  )
})

test_that("omics sources are parsed as exact comma-separated tokens", {
  expect_true(
    mapa:::.node_has_omics_source(list(dt_src = "P, T"), "T")
  )
  expect_true(
    mapa:::.node_has_omics_source(list(dt_src = "P, T"), "P")
  )
  expect_false(
    mapa:::.node_has_omics_source(list(dt_src = "TP"), "T")
  )
  expect_false(
    mapa:::.node_has_omics_source(list(symbol = "GENE1"), "T")
  )
})

test_that("module omics summary validates required input columns", {
  expect_error(
    mapa:::.summarise_module_omics(tibble::tibble(module = "module_1")),
    "Missing required columns: node_type, node_info",
    fixed = TRUE
  )
})
