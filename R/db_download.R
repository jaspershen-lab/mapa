# setwd(r4projects::get_project_wd())
# mapa_download_embedding_db()

# ── Edge database (STRING + Reactome) caching ─────────────────────────────────
# Files are stored under:
#   tools::R_user_dir("mapa", which = "cache") / edge_database /
# Downloaded once on first use; a failed download never leaves a corrupt file.

.mapa_edge_db_dir  <- function() file.path(tools::R_user_dir("mapa", which = "cache"), "edge_database")
.mapa_string_dir   <- function() file.path(.mapa_edge_db_dir(), "string_db")
.mapa_reactome_dir <- function() file.path(.mapa_edge_db_dir(), "reactome_db")

#' Check Whether the STRING and Reactome Edge Databases Are Cached
#'
#' Returns a list with logical flags \code{string_ok} and \code{reactome_ok},
#' plus the resolved \code{string_dir} and \code{reactome_dir} paths within
#' the mapa user cache directory.
#'
#' @return A named list with elements \code{string_ok}, \code{reactome_ok},
#'   \code{string_dir}, and \code{reactome_dir}.
#' @export
mapa_db_status_check <- function() {
  string_dir   <- .mapa_string_dir()
  reactome_dir <- .mapa_reactome_dir()
  list(
    string_ok    = all(file.exists(file.path(string_dir, c(
      "9606.protein.links.v12.0.txt.gz",
      "9606.protein.info.v12.0.txt.gz",
      "9606.protein.aliases.v12.0.txt.gz"
    )))),
    reactome_ok  = all(file.exists(file.path(reactome_dir, c(
      "ChEBI2Reactome_PE_Reactions.txt",
      "UniProt2Reactome_PE_Reactions.txt",
      "ProteinRoleReaction.txt"
    )))),
    string_dir   = string_dir,
    reactome_dir = reactome_dir
  )
}

#' Download and Cache the STRING and Reactome Edge Databases
#'
#' Ensures that all required STRING (PPI) and Reactome (enzyme–metabolite)
#' database files are present in the mapa user cache directory.  Files that
#' already exist are skipped; only missing files are downloaded.  Downloads are
#' written to a \code{.download} temp file and renamed only on success, so a
#' failed or interrupted download never leaves a corrupt cached file.
#'
#' @param on_file Optional callback \code{function(label)} called just before
#'   each new file download starts (useful for progress reporting in Shiny).
#' @return A named list with elements \code{string_dir} and
#'   \code{reactome_dir} giving the local paths to the two database
#'   directories.
#' @export
mapa_ensure_edge_databases <- function(on_file = NULL) {
  if (!requireNamespace("curl", quietly = TRUE))
    stop("Package 'curl' is required to download edge databases. ",
         "Install with: install.packages('curl')")

  string_dir   <- .mapa_string_dir()
  reactome_dir <- .mapa_reactome_dir()
  dir.create(string_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(reactome_dir, showWarnings = FALSE, recursive = TRUE)

  dl <- list(
    list(
      url   = "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
      dest  = file.path(string_dir, "9606.protein.links.v12.0.txt.gz"),
      label = "STRING protein links (~400 MB)"
    ),
    list(
      url   = "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
      dest  = file.path(string_dir, "9606.protein.info.v12.0.txt.gz"),
      label = "STRING protein info (~3 MB)"
    ),
    list(
      url   = "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz",
      dest  = file.path(string_dir, "9606.protein.aliases.v12.0.txt.gz"),
      label = "STRING protein aliases (~90 MB)"
    ),
    list(
      url   = "https://download.reactome.org/96/ChEBI2Reactome_PE_Reactions.txt",
      dest  = file.path(reactome_dir, "ChEBI2Reactome_PE_Reactions.txt"),
      label = "Reactome ChEBI reactions"
    ),
    list(
      url   = "https://download.reactome.org/96/UniProt2Reactome_PE_Reactions.txt",
      dest  = file.path(reactome_dir, "UniProt2Reactome_PE_Reactions.txt"),
      label = "Reactome UniProt reactions"
    ),
    list(
      url   = "https://download.reactome.org/96/ProteinRoleReaction.txt",
      dest  = file.path(reactome_dir, "ProteinRoleReaction.txt"),
      label = "Reactome protein roles"
    )
  )

  for (f in dl) {
    if (!file.exists(f$dest)) {
      if (is.function(on_file)) on_file(f$label)
      message("Downloading ", f$label, " ...")
      dest_tmp <- paste0(f$dest, ".download")
      tryCatch({
        curl::curl_download(f$url, dest_tmp, quiet = FALSE)
        file.rename(dest_tmp, f$dest)
        message("  Done.")
      }, error = function(e) {
        if (file.exists(dest_tmp)) file.remove(dest_tmp)
        stop("Failed to download ", f$label, ": ", e$message, call. = FALSE)
      })
    }
  }

  list(string_dir = string_dir, reactome_dir = reactome_dir)
}

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


#' Retrieve Source Metadata from the Local SQLite Database
#'
#' Returns the \code{source_metadata} table from the MAPA pathway embedding
#' database, which records the version and access method used to build each
#' pathway source (GO, KEGG, Reactome).
#'
#' @param version Character string. Database version (default \code{"v1"}).
#' @return A data frame with columns \code{source_db}, \code{source_version},
#'   and \code{access_method}.
#' @keywords internal
get_db_source_metadata <- function(version = "v1") {
  db_path <- mapa_download_embedding_db(version = version)
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  DBI::dbGetQuery(
    con,
    "SELECT source_db, source_version, access_method FROM source_metadata"
  )
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
