test_that("KEGG compound names preserve input order and use batches of ten", {
  calls <- list()
  fake_kegg_get <- function(ids) {
    calls[[length(calls) + 1L]] <<- ids
    lapply(sub("^cpd:", "", ids), function(id) {
      list(
        ENTRY = stats::setNames(id, "Compound"),
        NAME = c(paste0("name_", id, ";"), paste0("alias_", id))
      )
    })
  }

  ids <- sprintf("C%05d", 1:11)
  input <- c(ids, ids[1], "cpd:C00002", NA_character_, "invalid")

  result <- mapa:::.get_kegg_compound_names(input, fake_kegg_get)

  expect_length(calls, 2L)
  expect_lte(length(calls[[1]]), 10L)
  expect_lte(length(calls[[2]]), 10L)
  expect_equal(
    result,
    c(paste0("name_", ids), "name_C00001", "name_C00002", NA, NA)
  )
})

test_that("KEGG compound-name lookup degrades to NA when the API fails", {
  failing_kegg_get <- function(ids) stop("service unavailable")

  expect_warning(
    result <- mapa:::.get_kegg_compound_names(
      c("C00001", NA_character_),
      failing_kegg_get
    ),
    "Failed to retrieve KEGG compound names"
  )
  expect_equal(result, c(NA_character_, NA_character_))
})

test_that("KEGG compound-name lookup skips invalid IDs", {
  unexpected_kegg_get <- function(ids) {
    stop("KEGG should not be called")
  }

  result <- mapa:::.get_kegg_compound_names(
    c(NA_character_, "", "HMDB0000001", "C1234"),
    unexpected_kegg_get
  )

  expect_equal(result, rep(NA_character_, 4L))
})
