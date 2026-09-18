# Pathway-version harmonization helpers ---------------------------------------

.empty_pathway_version_harmonization_report <- function() {
  data.frame(
    original_id = character(),
    resolved_id = character(),
    database = character(),
    database_status = character(),
    is_obsolete = logical(),
    replaced_by = character(),
    replacement_count = integer(),
    action = character(),
    gene_set_revalidated = logical(),
    original_name = character(),
    resolved_name = character(),
    decision_reason = character(),
    status_source = character(),
    checked_at = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )
}

.extract_go_replacements <- function(term) {
  replacement_ids <- character()
  replacement_names <- character()

  structured_replacements <- c(term$replacements, term$replacedBy)
  if (length(structured_replacements) > 0) {
    replacement_values <- unlist(structured_replacements, recursive = TRUE,
                                 use.names = FALSE)
    replacement_ids <- c(
      replacement_ids,
      grep("^GO:[0-9]{7}$", as.character(replacement_values), value = TRUE)
    )
  }

  history <- term$history
  if (!is.null(history) && length(history) > 0) {
    for (entry in history) {
      if (!identical(toupper(entry$category %||% ""), "OBSOLETION") ||
          !identical(toupper(entry$action %||% ""), "ADDED")) {
        next
      }

      matched <- regexec(
        "^replaced_by\\s+(GO:[0-9]{7})(?:\\s+\\((.*)\\))?",
        entry$text %||% "",
        perl = TRUE
      )
      parts <- regmatches(entry$text %||% "", matched)[[1]]
      if (length(parts) >= 2) {
        replacement_ids <- c(replacement_ids, parts[[2]])
        replacement_names <- c(
          replacement_names,
          if (length(parts) >= 3) parts[[3]] else NA_character_
        )
      }
    }
  }

  replacement_ids <- unique(replacement_ids)
  names_by_id <- stats::setNames(rep(NA_character_, length(replacement_ids)),
                                 replacement_ids)
  if (length(replacement_names) > 0) {
    history_ids <- grep(
      "^GO:[0-9]{7}$",
      unlist(lapply(history, function(entry) {
        text <- entry$text %||% ""
        match <- regmatches(
          text,
          regexec("^replaced_by\\s+(GO:[0-9]{7})", text, perl = TRUE)
        )[[1]]
        if (length(match) >= 2) match[[2]] else NA_character_
      })),
      value = TRUE
    )
    for (i in seq_along(history_ids)) {
      if (i <= length(replacement_names) &&
          !is.na(replacement_names[[i]]) &&
          nzchar(replacement_names[[i]])) {
        names_by_id[[history_ids[[i]]]] <- replacement_names[[i]]
      }
    }
  }

  list(ids = replacement_ids, names = names_by_id)
}

.fetch_go_term_records <- function(go_ids) {
  go_ids <- unique(as.character(go_ids))
  go_ids <- go_ids[!is.na(go_ids) & nzchar(go_ids)]
  if (length(go_ids) == 0) return(list())

  chunks <- split(go_ids, ceiling(seq_along(go_ids) / 100))
  records <- list()

  for (chunk in chunks) {
    url <- paste0(
      "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/",
      paste(chunk, collapse = ","),
      "/history"
    )

    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_headers("Accept" = "application/json") |>
        httr2::req_retry(max_tries = 3) |>
        httr2::req_perform(),
      error = function(e) {
        stop(
          "Unable to verify GO term status with QuickGO: ", e$message,
          ". MAPA stopped before downstream similarity or interpretation ",
          "to avoid retaining terms with an unknown version status.",
          call. = FALSE
        )
      }
    )
    results <- httr2::resp_body_json(response, simplifyVector = FALSE)$results

    for (original_id in chunk) {
      matches <- Filter(function(term) {
        identical(term$id, original_id) ||
          original_id %in% unlist(term$secondaryIds %||% character(),
                                  use.names = FALSE)
      }, results)
      records[[original_id]] <- if (length(matches) == 1) matches[[1]] else NULL
    }
  }

  records
}

#' Harmonize GO pathway identifiers against the current ontology
#'
#' Reads the explicit QuickGO `isObsolete` status for every GO identifier. A
#' term is mapped only when its obsoletion history contains exactly one
#' `replaced_by` identifier. Obsolete terms without a unique replacement are
#' excluded from downstream analysis.
#'
#' @param go_ids Character vector of GO identifiers.
#' @return A data frame named `pathway_version_harmonization_report` by MAPA
#'   callers, containing the original and resolved identifiers, database status,
#'   replacement decision, and gene-set revalidation status.
#' @export
harmonize_pathway_versions <- function(go_ids) {
  .harmonize_go_ids(go_ids)
}

.harmonize_go_ids <- function(go_ids, term_records = NULL) {
  go_ids <- unique(as.character(go_ids))
  go_ids <- go_ids[!is.na(go_ids) & nzchar(go_ids)]
  if (length(go_ids) == 0) {
    return(.empty_pathway_version_harmonization_report())
  }
  invalid <- !grepl("^GO:[0-9]{7}$", go_ids)
  if (any(invalid)) {
    stop("Invalid GO identifier(s): ", paste(go_ids[invalid], collapse = ", "),
         call. = FALSE)
  }

  if (is.null(term_records)) {
    term_records <- .fetch_go_term_records(go_ids)
  }
  checked_at <- Sys.time()

  rows <- lapply(go_ids, function(original_id) {
    term <- term_records[[original_id]]
    if (is.null(term)) {
      return(data.frame(
        original_id = original_id,
        resolved_id = NA_character_,
        database = "GO",
        database_status = "not_found",
        is_obsolete = NA,
        replaced_by = NA_character_,
        replacement_count = 0L,
        action = "excluded",
        gene_set_revalidated = NA,
        original_name = NA_character_,
        resolved_name = NA_character_,
        decision_reason = "GO identifier was not returned by QuickGO",
        status_source = "QuickGO",
        checked_at = checked_at,
        stringsAsFactors = FALSE
      ))
    }

    is_obsolete <- term$isObsolete
    if (!is.logical(is_obsolete) || length(is_obsolete) != 1 ||
        is.na(is_obsolete)) {
      stop("QuickGO did not return an explicit isObsolete status for ",
           original_id, ".", call. = FALSE)
    }

    replacements <- .extract_go_replacements(term)
    replacement_count <- length(replacements$ids)
    canonical_id <- term$id %||% original_id

    if (!is_obsolete) {
      resolved_id <- canonical_id
      action <- if (identical(resolved_id, original_id)) "retained" else
        "mapped_secondary_id"
      resolved_name <- term$name %||% NA_character_
      reason <- if (action == "retained") {
        "QuickGO reports the term as active"
      } else {
        "QuickGO resolved the secondary identifier to its active primary identifier"
      }
    } else if (replacement_count == 1L) {
      resolved_id <- replacements$ids[[1]]
      action <- "replaced"
      resolved_name <- unname(replacements$names[[resolved_id]])
      reason <- paste0(
        "QuickGO reports the term as obsolete with one replaced_by target; ",
        "the original enrichment gene set has not been revalidated"
      )
    } else {
      resolved_id <- NA_character_
      action <- "excluded"
      resolved_name <- NA_character_
      reason <- if (replacement_count == 0L) {
        "QuickGO reports the term as obsolete without a replaced_by target"
      } else {
        "QuickGO reports multiple replaced_by targets; no unique mapping is safe"
      }
    }

    data.frame(
      original_id = original_id,
      resolved_id = resolved_id,
      database = "GO",
      database_status = if (is_obsolete) "obsolete" else "active",
      is_obsolete = is_obsolete,
      replaced_by = if (replacement_count == 0L) NA_character_ else
        paste(replacements$ids, collapse = ";"),
      replacement_count = as.integer(replacement_count),
      action = action,
      gene_set_revalidated = if (action == "replaced") FALSE else NA,
      original_name = term$name %||% NA_character_,
      resolved_name = resolved_name,
      decision_reason = reason,
      status_source = "QuickGO",
      checked_at = checked_at,
      stringsAsFactors = FALSE
    )
  })

  report <- do.call(rbind, rows)
  rownames(report) <- NULL
  report
}

.apply_go_harmonization <- function(result, report, id_col = "ID",
                                    description_col = "Description") {
  if (is.null(result) || nrow(result) == 0 || nrow(report) == 0) return(result)
  if (!id_col %in% names(result)) {
    stop("GO enrichment table is missing identifier column `", id_col, "`.",
         call. = FALSE)
  }

  report_index <- match(as.character(result[[id_col]]), report$original_id)
  keep <- !is.na(report_index) & report$action[report_index] != "excluded" &
    !is.na(report$resolved_id[report_index])
  result <- result[keep, , drop = FALSE]
  report_index <- report_index[keep]
  if (nrow(result) == 0) return(result)

  original_ids <- as.character(result[[id_col]])
  result[[id_col]] <- report$resolved_id[report_index]
  mapped <- original_ids != result[[id_col]]

  if (description_col %in% names(result)) {
    replacement_names <- report$resolved_name[report_index]
    update_name <- mapped & !is.na(replacement_names) & nzchar(replacement_names)
    result[[description_col]][update_name] <- replacement_names[update_name]
  }

  action <- report$action[report_index]
  priority <- ifelse(action == "retained", 1L,
                     ifelse(action == "mapped_secondary_id", 2L, 3L))
  p_adjust <- if ("p_adjust" %in% names(result)) {
    suppressWarnings(as.numeric(result$p_adjust))
  } else {
    rep(Inf, nrow(result))
  }
  p_adjust[is.na(p_adjust)] <- Inf
  order_index <- order(priority, p_adjust, seq_len(nrow(result)))
  result <- result[order_index, , drop = FALSE]
  result <- result[!duplicated(result[[id_col]]), , drop = FALSE]
  rownames(result) <- NULL
  result
}

.harmonize_go_result_object <- function(object, report = NULL) {
  if (is.null(object)) {
    return(list(object = object,
                pathway_version_harmonization_report =
                  .empty_pathway_version_harmonization_report()))
  }
  go_ids <- as.character(object@result$ID)
  if (is.null(report)) report <- .harmonize_go_ids(go_ids)
  object@result <- .apply_go_harmonization(object@result, report)
  list(object = object, pathway_version_harmonization_report = report)
}

.combine_harmonization_reports <- function(...) {
  reports <- list(...)
  reports <- reports[vapply(reports, function(x) !is.null(x) && nrow(x) > 0,
                            logical(1))]
  if (length(reports) == 0) {
    return(.empty_pathway_version_harmonization_report())
  }
  result <- do.call(rbind, reports)
  result <- result[!duplicated(result$original_id), , drop = FALSE]
  rownames(result) <- NULL
  result
}
