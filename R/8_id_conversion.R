# setwd("demo_data/")
# load("demo_data.rda")
# variable_info <-
#   demo_data %>%
#   massdataset::activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0) %>%
#   massdataset::extract_variable_info()
# variable_info <- variable_info |> dplyr::select(variable_id, symbol)
### Create an AnnotationHub object to access the database
# ah <- AnnotationHub::AnnotationHub()
# ### Query the object for Macaca fascicularis
# mf_query_result <- AnnotationHub::query(ah, c("Macaca fascicularis", "OrgDb"))
# ### Get the orgDb object
# mf.orgdb <- ah[["AH119900"]] #| Taxonomy ID: 9541
# converted_dt <- convert_id(
#   data = dt_u,
#   query_type = "gene",
#   from_id_type = "ensembl",
#   organism = mf.orgdb
# )
# converted_dt <- convert_id(
#   data = dt_u,
#   query_type = "gene",
#   from_id_type = "ensembl",
#   ah_id = "AH119902"
# )
# converted_dt <- convert_id(
#   data = dt_u,
#   query_type = "gene",
#   from_id_type = "ensembl",
#   ah_id = "AH119902",
#   return_orgdb = TRUE
# )

#' Convert Between Different ID Types for Genes or Metabolites
#'
#' This function converts between different identifier types for genes or metabolites.
#' For genes, it uses the clusterProfiler package and the output will always contain
#' four columns: ensembl, entrezid, uniprot, and symbol. For metabolites (human only),
#' it converts between HMDB and KEGG identifiers using the metpath package.
#'
#' @param data Data frame containing the data with identifiers to be converted.
#'   Must contain a column with the identifier type specified in \code{from_id_type}.
#' @param query_type Character vector. Type of data to convert. Must be either
#'   "gene" or "metabolite". Default is c("gene", "metabolite").
#' @param from_id_type Character string. The identifier type to convert from.
#'   For genes, must be one of "ensembl", "entrezid", "uniprot", or "symbol".
#'   For metabolites, supported types are "hmdbid" and "keggid".
#' @param organism Character string or OrgDb object. For genes, specify the
#'   organism database to use (e.g., "org.Hs.eg.db" for human). For metabolites,
#'   must be "hsa" (human).
#' @param ah_id Character string. AnnotationHub ID to use when the orgDb package
#'   is not available. If provided, the function will retrieve the organism database
#'   from AnnotationHub using this ID. Only used for gene conversion.
#' @param return_orgdb Logical. If TRUE and ah_id is provided, returns a list containing
#'   both the converted data and the OrgDb object for downstream analysis. If FALSE,
#'   returns only the converted data. Default is FALSE.
#'
#' @return A data frame with the original data joined with the converted identifiers.
#'   For genes, the output will always contain four required columns: ensembl,
#'   entrezid, uniprot, and symbol. The column corresponding to from_id_type is
#'   preserved from the original data. If no identifier is available for conversion,
#'   NA is returned. If multiple identifiers are available for one from_id_type,
#'   only the first match is retained. For metabolites, both HMDB and KEGG IDs are
#'   added regardless of the input type. If return_orgdb = TRUE and ah_id is provided,
#'   returns a list with elements 'data' (converted data) and 'orgdb' (OrgDb object).
#'
#' @examples
#' \dontrun{
#' # Convert from gene symbols using standard org.Db package
#' gene_data <- data.frame(symbol = c("TP53", "BRCA1", "EGFR"))
#' converted_genes <- convert_id(
#'   data = gene_data,
#'   query_type = "gene",
#'   from_id_type = "symbol",
#'   organism = "org.Hs.eg.db"
#' )
#'
#' # Convert using AnnotationHub ID and return both data and OrgDb object
#' gene_data_mouse <- data.frame(symbol = c("Tp53", "Brca1", "Egfr"))
#' result_with_orgdb <- convert_id(
#'   data = gene_data_mouse,
#'   query_type = "gene",
#'   from_id_type = "symbol",
#'   ah_id = "AH123456",  # Example AnnotationHub ID
#'   return_orgdb = TRUE
#' )
#' # Access converted data: result_with_orgdb$data
#' # Access OrgDb object: result_with_orgdb$orgdb
#'
#' # Convert from ENTREZ IDs
#' gene_data_entrez <- data.frame(entrezid = c("7157", "672", "1956"))
#' converted_genes2 <- convert_id(
#'   data = gene_data_entrez,
#'   query_type = "gene",
#'   from_id_type = "entrezid",
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
#' @importFrom dplyr distinct left_join select rename mutate filter
#'
#' @export
convert_id <- function(data = NULL,
                       query_type = c("gene", "metabolite"),
                       from_id_type = NULL,
                       organism = NULL,
                       ah_id = NULL,
                       return_orgdb = FALSE) {

  # Validate inputs
  if (is.null(data)) {
    stop("data is required")
  }

  if (missing(query_type)) {
    stop("query_type is missing")
  }
  query_type <- match.arg(query_type, c("gene", "metabolite"))

  if (is.null(from_id_type)) {
    stop("from_id_type is required")
  }

  if (is.null(organism) && is.null(ah_id)) {
    stop("Either organism or ah_id is required")
  }

  # Handle gene conversion
  if (query_type == "gene") {
    # Check if from_id_type is valid for genes
    valid_gene_id_types <- c("ensembl", "entrezid", "uniprot", "symbol")
    if (!from_id_type %in% valid_gene_id_types) {
      stop("For genes, from_id_type must be one of: ", paste(valid_gene_id_types, collapse = ", "))
    }

    # Check if the required column exists in data
    if (!from_id_type %in% tolower(names(data))) {
      stop("Column '", from_id_type, "' not found in data")
    }

    # Load required packages
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
      BiocManager::install("clusterProfiler")
    }

    # Handle organism database - either from package or AnnotationHub
    if (!is.null(ah_id)) {
      if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
        BiocManager::install("AnnotationHub")
      }

      message("Retrieving organism database from AnnotationHub using ID: ", ah_id)
      ah <- AnnotationHub::AnnotationHub()
      organism <- ah[[ah_id]]

      # Display organism database information using show() method
      message("Successfully loaded organism database from AnnotationHub")
      message("Database information:")
      show(organism)
    }

    # Get available ID types for the organism database
    available_id_types <- clusterProfiler::idType(organism)

    # Define mapping between our column names and clusterProfiler ID types
    id_type_mapping <- list(
      "symbol" = "SYMBOL",
      "entrezid" = "ENTREZID",
      "ensembl" = "ENSEMBL",
      "uniprot" = "UNIPROT"
    )

    # Check which ID types are unavailable and inform user
    unavailable_types <- c()
    for (id_type in names(id_type_mapping)) {
      if (!id_type_mapping[[id_type]] %in% available_id_types) {
        unavailable_types <- c(unavailable_types, id_type)
        id_type_mapping[[id_type]] <- NA  # Mark as unavailable
      }
    }

    if (length(unavailable_types) > 0) {
      message("Note: The following ID types are not available in the organism database and will be filled with NA: ",
              paste(unavailable_types, collapse = ", "))
      message("Available ID types in database: ", paste(available_id_types, collapse = ", "))
    }

    # Get the clusterProfiler ID type for the from_id_type
    from_clusterprofiler_type <- id_type_mapping[[from_id_type]]

    # Check if from_id_type is available
    if (is.na(from_clusterprofiler_type)) {
      stop("The from_id_type '", from_id_type, "' is not available in the organism database. Available types: ",
           paste(available_id_types, collapse = ", "))
    }

    # Define the target ID types (all except the from_id_type)
    target_types <- valid_gene_id_types[valid_gene_id_types != from_id_type]
    target_clusterprofiler_types <- unlist(id_type_mapping[target_types])

    # Remove NA values (unavailable types) from target types
    target_clusterprofiler_types <- target_clusterprofiler_types[!is.na(target_clusterprofiler_types)]
    # Perform the conversion only if there are target types available
    if (length(target_clusterprofiler_types) > 0) {
      tryCatch({
        converted <- clusterProfiler::bitr(
          geneID = data[[from_id_type]],
          fromType = from_clusterprofiler_type,
          toType = target_clusterprofiler_types,
          OrgDb = organism
        )

        # Convert column names to lowercase to match our naming convention
        names(converted) <- tolower(names(converted))

        # Keep only first match for each from_id_type to avoid duplicates
        converted <- converted |>
          dplyr::distinct(!!sym(from_id_type), .keep_all = TRUE)

        # Join with original data
        result <- data |>
          dplyr::left_join(converted, by = from_id_type)

      }, error = function(e) {
        warning("ID conversion failed: ", e$message, ". Returning original data with NA values for missing ID types.")
        result <- data
      })
    } else {
      message("No target ID types available for conversion from '", from_id_type, "'. All other columns will be filled with NA.")
      result <- data
    }

    # Ensure all four required columns are present
    for (col in valid_gene_id_types) {
      if (!col %in% names(result)) {
        result[[col]] <- NA_character_
      }
    }

    # Reorder columns to have the four ID columns first
    other_cols <- setdiff(names(result), valid_gene_id_types)
    result <- result |>
      dplyr::select(all_of(valid_gene_id_types), all_of(other_cols))

    # Return based on return_orgdb parameter
    if (!is.null(ah_id) && return_orgdb) {
      return(list(data = result, orgdb = organism))
    } else {
      return(result)
    }
  }

  # Handle metabolite conversion
  if (query_type == "metabolite" && organism == "hsa") {
    # Check if from_id_type is valid for metabolites
    valid_metabolite_id_types <- c("hmdbid", "keggid")
    if (!from_id_type %in% valid_metabolite_id_types) {
      stop("For metabolites, from_id_type must be one of: ", paste(valid_metabolite_id_types, collapse = ", "))
    }

    # Check if the required column exists in data
    if (!from_id_type %in% tolower(names(data))) {
      stop("Column '", from_id_type, "' not found in data")
    }

    if (!requireNamespace("metpath", quietly = TRUE)) {
      BiocManager::install("metpath")
    }

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
