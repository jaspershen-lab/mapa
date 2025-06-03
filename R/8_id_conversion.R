# setwd("demo_data/")
# load("demo_data.rda")
# variable_info <-
#   demo_data %>%
#   massdataset::activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0) %>%
#   massdataset::extract_variable_info()
# variable_info <- variable_info |> dplyr::select(variable_id, symbol)
# converted_dt <- convert_id(
#   data = variable_info,
#   query_type = "gene",
#   from_id_type = "SYMBOL",
#   to_id_type = c("ENSEMBL"),
#   organism = org.Hs.eg.db
# )

#' Convert Between Different ID Types for Genes or Metabolites
#'
#' This function converts between different identifier types for genes or metabolites.
#' For genes, it uses the clusterProfiler package to perform ID conversion via
#' organism-specific databases. For metabolites (human only), it converts between
#' HMDB and KEGG identifiers using the metpath package.
#'
#' @param data Data frame containing the data with identifiers to be converted.
#'   Must contain a column with the identifier type specified in \code{from_id_type}.
#' @param query_type Character vector. Type of data to convert. Must be either
#'   "gene" or "metabolite". Default is c("gene", "metabolite").
#' @param from_id_type Character string. The identifier type to convert from.
#'   For genes, this should match valid identifier types supported by the
#'   organism database (e.g., "SYMBOL", "ENTREZID", "ENSEMBL"). For metabolites,
#'   supported types are "hmdbid" and "keggid".
#' @param to_id_type Character vector. The identifier types to convert to.
#'   ENTREZID is automatically included unless the input data already contains
#'   an "entrezid" column. For genes, types should match valid identifier types
#'   supported by the organism database. For metabolites, this parameter is not
#'   used as conversion is bidirectional between HMDB and KEGG IDs. Default is c("ENTREZID").
#' @param organism Character string or OrgDb object. For genes, specify the
#'   organism database to use (e.g., "org.Hs.eg.db" for human). For metabolites,
#'   must be "hsa" (human).
#'
#' @return A data frame with the original data joined with the converted identifiers.
#'   For genes, ENTREZID is included in the results unless the input data already
#'   contains an "entrezid" column. Any additional requested identifier types are
#'   also added as new columns with lowercase names. For metabolites, both HMDB
#'   and KEGG IDs are added regardless of the input type.
#'
#'
#' @examples
#' \dontrun{
#' # Convert gene symbols to ENSEMBL IDs (ENTREZID automatically included)
#' gene_data <- data.frame(symbol = c("TP53", "BRCA1", "EGFR"))
#' converted_genes <- convert_id(
#'   data = gene_data,
#'   query_type = "gene",
#'   from_id_type = "SYMBOL",
#'   to_id_type = "ENSEMBL",
#'   organism = "org.Hs.eg.db"
#' )
#'
#' # If data already has ENTREZID, only convert to requested types
#' gene_data_with_entrez <- data.frame(
#'   symbol = c("TP53", "BRCA1", "EGFR"),
#'   entrezid = c("7157", "672", "1956")
#' )
#' converted_genes2 <- convert_id(
#'   data = gene_data_with_entrez,
#'   query_type = "gene",
#'   from_id_type = "SYMBOL",
#'   to_id_type = "ENSEMBL",  # Only ENSEMBL will be added, not ENTREZID
#'   organism = "org.Hs.eg.db"
#' )
#'
#' # Convert HMDB IDs to KEGG IDs for metabolites
#' metabolite_data <- data.frame(hmdbid = c("HMDB0000001", "HMDB0000002"))
#' converted_metabolites <- convert_id(
#'   data = metabolite_data,
#'   query_type = "metabolite",
#'   from_id_type = "hmdbid",
#'   organism = "hsa"
#' )
#' }
#'
#' @importFrom dplyr distinct rename_with left_join select rename mutate filter across
#'
#' @export
convert_id <- function(data = NULL,
                       query_type = c("gene", "metabolite"),
                       from_id_type = NULL,
                       to_id_type = c("ENTREZID"),
                       organism = NULL) {
  if (missing(query_type)){
    stop("query_type is missing")
  }
  query_type <- match.arg(query_type, c("gene", "metabolite"))

  if (query_type == "gene") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE))
      BiocManager::install("clusterProfiler")

    # Check if data already has ENTREZID column
    has_entrezid <- "entrezid" %in% tolower(names(data))

    # Only include ENTREZID if not already present
    if (!has_entrezid) {
      to_id_type <- unique(c("ENTREZID", to_id_type))
    }

    converted <- clusterProfiler::bitr(
      geneID  = data[[tolower(from_id_type)]],
      fromType = from_id_type,
      toType   = to_id_type,
      OrgDb    = organism
    ) |>
      dplyr::distinct(ENTREZID, .keep_all = TRUE) |>
      dplyr::rename_with(tolower) |>
      dplyr::left_join(data, ., by = tolower(from_id_type))

    return(converted)
  }

  if (query_type == "metabolite" && organism == "hsa") {
    if (!requireNamespace("metpath", quietly = TRUE)) {BiocManager::install("metpath")}

    id_lookup <- metpath::hmdb_compound_database@spectra.info |>
      dplyr::select(HMDB.ID, KEGG.ID) |>
      dplyr::rename(hmdbid = HMDB.ID,
                    keggid = KEGG.ID) |>
      dplyr::mutate(dplyr::across(everything(), as.character)) |>
      dplyr::distinct()

    converted <- data |>
      dplyr::filter(!is.na(.data[[tolower(from_id_type)]])) |>
      dplyr::left_join(id_lookup, by = tolower(from_id_type))

    return(converted)
  }
}
