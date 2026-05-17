# setwd(r4projects::get_project_wd())
# mapa_download_embedding_db()

#' Download the MAPA Pathway Embedding Database
#'
#' Downloads a pre-built SQLite embedding database from the MAPA GitHub release
#' and stores it in the user's cache directory. If the database file already
#' exists, the download is skipped unless \code{overwrite = TRUE}.
#'
#' @param version Character string specifying the database version to download.
#'   Currently only \code{"v1"} is supported. Default is \code{"v1"}.
#' @param overwrite Logical. If \code{TRUE}, re-download and overwrite an
#'   existing database file. Default is \code{FALSE}.
#'
#' @return The file path to the downloaded \code{.sqlite} database (invisibly).
#'
#' @details
#' The database is saved to the platform-specific user cache directory returned
#' by \code{tools::R_user_dir("mapa", which = "cache")}. The compressed
#' \code{.gz} archive is decompressed automatically via \code{R.utils::gunzip()}.
#'
#' @examples
#' \dontrun{
#' db_path <- mapa_download_embedding_db()
#' db_path <- mapa_download_embedding_db(overwrite = TRUE)
#' }
#'
#' @export
mapa_download_embedding_db <- function(version = "v1",
                                       overwrite = FALSE) {
  db_dir <- tools::R_user_dir("mapa", which = "cache")
  dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)

  db_path <- file.path(db_dir, paste0("mapa_pathway_embedding_", version, ".sqlite"))

  if (file.exists(db_path) && !overwrite) {
    message("MAPA embedding database already exists: ", db_path)
    return(invisible(db_path))
  }

  url <- switch(
    version,
    "v1" = "https://github.com/jaspershen-lab/mapa/releases/download/v2.1.0/mapa_pathway_embedding_v1.sqlite.gz",
    stop("Unknown database version: ", version)
  )

  gz_path <- tempfile(fileext = ".sqlite.gz")
  on.exit(unlink(gz_path), add = TRUE)

  # Large file (~440 MB): raise the download timeout well above the R default
  # (60 s) and restore the original value when the function exits.
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(old_timeout, 3600L))

  message("Downloading MAPA embedding database (version ", version, ") — this may take several minutes...")
  utils::download.file(url, gz_path, mode = "wb")

  R.utils::gunzip(
    gz_path,
    destname = db_path,
    overwrite = TRUE,
    remove = FALSE
  )

  message("MAPA embedding database downloaded to: ", db_path)
  invisible(db_path)
}


#' Retrieve Pathway Embeddings from the Local SQLite Database
#'
#' Queries the pre-built pathway embedding database for a set of pathway IDs
#' and returns a named numeric matrix ready for cosine similarity computation.
#' If the database has not been downloaded yet it is fetched automatically via
#' \code{\link{mapa_download_embedding_db}}.
#'
#' Organism-specific KEGG pathway IDs (e.g. \code{"hsa04010"}) are
#' automatically mapped to the reference pathway IDs stored in the database
#' (e.g. \code{"map04010"}) before the lookup.  The original IDs are restored
#' as row names in the returned matrix so downstream code remains unaffected.
#'
#' @param pathway_ids Character vector of pathway identifiers to look up
#'   (e.g. \code{"GO:0006955"}, \code{"hsa04110"}, \code{"R-HSA-168256"}).
#' @param version Character string. Database version passed to
#'   \code{\link{mapa_download_embedding_db}} (default \code{"v1"}).
#'
#' @return A numeric matrix with one row per found pathway and one column per
#'   embedding dimension.  Row names are the original (caller-supplied) pathway
#'   IDs.  A warning is issued for any IDs absent from the database.
#'
#' @importFrom DBI dbConnect dbDisconnect dbGetQuery
#' @importFrom RSQLite SQLite
#'
#' @keywords internal
get_pathway_embeddings_from_db <- function(pathway_ids, version = "v1") {
  db_path <- mapa_download_embedding_db(version = version)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Map organism-specific KEGG IDs (e.g. "hsa04010") to reference IDs ("map04010").
  # The DB only stores reference KEGG pathways with the "map" prefix.
  kegg_re <- "^[a-z]{3}(\\d+)$"
  db_ids <- ifelse(
    grepl(kegg_re, pathway_ids),
    paste0("map", sub(kegg_re, "\\1", pathway_ids)),
    pathway_ids
  )

  # Named vector: db_id -> original caller ID (used to restore row names)
  unique_db_ids <- unique(db_ids)
  orig_for_db   <- stats::setNames(pathway_ids, db_ids)[unique_db_ids]

  placeholders <- paste(rep("?", length(unique_db_ids)), collapse = ",")
  query <- sprintf(
    "SELECT id, embedding FROM pathway_embedding WHERE id IN (%s)",
    placeholders
  )
  result <- DBI::dbGetQuery(con, query, params = as.list(unique_db_ids))

  if (nrow(result) == 0) {
    stop(
      "No embeddings found in the local database for the provided pathway IDs. ",
      "Consider using embedding_source = \"api\" instead."
    )
  }

  missing_db_ids <- setdiff(unique_db_ids, result$id)
  if (length(missing_db_ids) > 0) {
    missing_orig <- unname(orig_for_db[missing_db_ids])
    warning(
      length(missing_db_ids), " pathway ID(s) not found in the local database and will be skipped: ",
      paste(head(missing_orig, 5), collapse = ", "),
      if (length(missing_orig) > 5) ", ..." else ""
    )
  }

  # Decode float32 BLOBs (stored as little-endian 4-byte floats, not R serialize)
  embeddings <- lapply(result$embedding, function(x) {
    readBin(x, what = "double", n = length(x) %/% 4L, size = 4L, endian = "little")
  })
  mat <- do.call(rbind, embeddings)

  # Restore original caller IDs as row names
  rownames(mat) <- unname(orig_for_db[result$id])
  mat
}
