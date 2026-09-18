.go_term_record <- function(id, obsolete = FALSE, name = id,
                            obsoletion_text = character(),
                            secondary_ids = NULL) {
  history <- lapply(obsoletion_text, function(text) {
    list(
      action = "Added",
      category = "OBSOLETION",
      text = text
    )
  })
  list(
    id = id,
    isObsolete = obsolete,
    name = name,
    secondaryIds = secondary_ids,
    history = history
  )
}

test_that("GO harmonization uses explicit status and only unique replaced_by", {
  ids <- c("GO:0000001", "GO:0000003", "GO:0044332", "GO:0999999")
  records <- list(
    "GO:0000001" = .go_term_record(
      "GO:0000001", name = "mitochondrion inheritance"
    ),
    "GO:0000003" = .go_term_record(
      "GO:0000003",
      obsolete = TRUE,
      name = "obsolete reproduction",
      obsoletion_text = "replaced_by GO:0022414 (reproductive process)"
    ),
    "GO:0044332" = .go_term_record(
      "GO:0044332",
      obsolete = TRUE,
      obsoletion_text = c(
        "consider GO:0016055 (Wnt signaling pathway)",
        "consider GO:0009950 (dorsal/ventral axis specification)"
      )
    ),
    "GO:0999999" = .go_term_record(
      "GO:0999999",
      obsolete = TRUE,
      obsoletion_text = c(
        "replaced_by GO:0008150 (biological_process)",
        "replaced_by GO:0003674 (molecular_function)"
      )
    )
  )

  report <- mapa:::.harmonize_go_ids(ids, term_records = records)

  expect_identical(
    report$action,
    c("retained", "replaced", "excluded", "excluded")
  )
  expect_identical(
    report$resolved_id,
    c("GO:0000001", "GO:0022414", NA_character_, NA_character_)
  )
  expect_identical(report$database_status,
                   c("active", "obsolete", "obsolete", "obsolete"))
  expect_false(report$gene_set_revalidated[[2]])
  expect_match(report$decision_reason[[2]], "not been revalidated")
})

test_that("GO harmonization records active secondary-ID resolution", {
  records <- list(
    "GO:1234567" = .go_term_record(
      "GO:0008150",
      name = "biological_process",
      secondary_ids = "GO:1234567"
    )
  )

  report <- mapa:::.harmonize_go_ids("GO:1234567", term_records = records)

  expect_identical(report$action, "mapped_secondary_id")
  expect_identical(report$original_id, "GO:1234567")
  expect_identical(report$resolved_id, "GO:0008150")
})

test_that("harmonization removes excluded rows and prefers active duplicates", {
  result <- data.frame(
    ID = c("GO:0000003", "GO:0022414", "GO:0044332"),
    Description = c("obsolete reproduction", "reproductive process", "obsolete"),
    p_adjust = c(0.001, 0.01, 0.02),
    geneID = c("1/2", "1/2/3", "4"),
    stringsAsFactors = FALSE
  )
  records <- list(
    "GO:0000003" = .go_term_record(
      "GO:0000003", TRUE, "obsolete reproduction",
      "replaced_by GO:0022414 (reproductive process)"
    ),
    "GO:0022414" = .go_term_record(
      "GO:0022414", FALSE, "reproductive process"
    ),
    "GO:0044332" = .go_term_record("GO:0044332", TRUE, "obsolete")
  )
  report <- mapa:::.harmonize_go_ids(result$ID, term_records = records)

  harmonized <- mapa:::.apply_go_harmonization(result, report)

  expect_identical(harmonized$ID, "GO:0022414")
  expect_identical(harmonized$geneID, "1/2/3")
  expect_identical(harmonized$Description, "reproductive process")
})

test_that("empty harmonization report keeps a stable schema", {
  report <- mapa:::.harmonize_go_ids(character())

  expect_identical(nrow(report), 0L)
  expect_named(
    report,
    c(
      "original_id", "resolved_id", "database", "database_status",
      "is_obsolete", "replaced_by", "replacement_count", "action",
      "gene_set_revalidated", "original_name", "resolved_name",
      "decision_reason", "status_source", "checked_at"
    )
  )
})
