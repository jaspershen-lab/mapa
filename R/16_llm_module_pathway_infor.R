# This file is for input pathway data processing and search integration information.
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# source('R/17_llm_interpretation/modules/get_pathway_name_desc_pmid.R')

#' Preprocess Module Data
#'
#' Processes a data frame of module data, extracts pathway names, gene symbols, and gene descriptions, and organizes the data into a structured list.
#'
#' @param df A data frame containing the module data. The data frame must include the following columns:
#' \itemize{
#'   \item \code{Description}: Pathway descriptions, separated by semicolons (\code{;}).
#'   \item \code{geneID} or \code{mapped_id}: Gene IDs or metabolite IDs, separated by slashes (\code{/}).
#'   \item \code{module}: Module identifiers.
#'   \item \code{pathway_id} or \code{node}: Pathway identifiers.
#' }
#' @param orgdb An annotation database object from \code{org.Hs.eg.db} or similar packages. Defaults to \code{org.Hs.eg.db}.
#'
#' @return A named list where each element corresponds to a module. Each module contains:
#' \itemize{
#'   \item \code{PathwayNames}: A character vector of pathway names.
#'   \item \code{PathwayDescription}: A character vector of pathway descriptions.
#'   \item \code{PathwayReferencePMID}: A character vector of PMIDs for pathway references.
#'   \item For gene modules:
#'     \itemize{
#'       \item \code{GeneIDs}: A character vector of gene IDs.
#'       \item \code{GeneSymbols}: A character vector of gene symbols (mapped from gene IDs).
#'       \item \code{GeneNames_vec}: A character vector of gene descriptions (mapped from gene IDs).
#'     }
#'   \item For metabolite modules:
#'     \itemize{
#'       \item \code{MetIDs}: A character vector of metabolite IDs.
#'       \item \code{MetNames_vec}: A character vector of metabolite names.
#'     }
#' }
#'
#' @importFrom dplyr rename filter pull
#' @importFrom AnnotationDbi mapIds
#'
#' @examples
#' \dontrun{
#' # Example: Preprocess a data frame of module data
#' library(org.Hs.eg.db)
#' module_data <- data.frame(
#'   Description = c("Pathway A; Pathway B", "Pathway C"),
#'   geneID = c("ENSG000001;ENSG000002", "ENSG000003/ENSG000004"),
#'   module = c("Module1", "Module2"),
#'   pathway_id = c("GO:0001234;GO:0005678", "hsa04210")
#' )
#' processed_data <- preprocess_module(module_data, orgdb = org.Hs.eg.db)
#' print(processed_data)
#' }
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Identifies whether the input data is for genes or metabolites.
#'   \item Splits pathway descriptions and gene/metabolite IDs into vectors for each row in the data frame.
#'   \item Maps gene IDs to gene symbols and gene descriptions using the specified annotation database (\code{orgdb}).
#'   \item Retrieves pathway descriptions and references based on pathway IDs.
#'   \item Organizes the results into a named list, where each module is represented as an element containing pathway and gene/metabolite information.
#' }
#'
#' The function handles both gene and metabolite data differently, based on the columns present in the input data frame.
#'
#' @seealso \code{\link[AnnotationDbi]{mapIds}} for gene ID mapping.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
preprocess_module <- function(df,
                              orgdb = org.Hs.eg.db  # Annotation database, can be replaced based on species
                              ) {
  if ("geneID" %in% colnames(df)) {
    # functional_module_result for genes
    query_type <- "gene"
    df <- df %>%
      dplyr::rename(queryid = geneID)
    if ("node" %in% colnames(df)) {
      df <- df %>%
        dplyr::rename(pathway_id = node)
    }
  } else if ("mapped_id" %in% colnames(df)) {
    # functional_module_result for metabolites
    query_type <- "metabolite"
    df <- df %>%
      dplyr::rename(queryid = mapped_id) %>%
      dplyr::rename(pathway_id = node)
  }

  required_cols <- c("Description", "queryid", "module", "pathway_id")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Data frame is missing the following columns:", paste(missing_cols, collapse = ", ")))
  }

  result_list <- list()

  for (i in seq_len(nrow(df))) {
    # Collect module information
    pthName_raw <- df$Description[i]
    qID_raw <- df$queryid[i]
    mod_raw  <- df$module[i]
    pthID_raw <- df$pathway_id[i]

    pthName_vec <- unlist(strsplit(pthName_raw, ";"))
    pthID_vec <- unlist(strsplit(pthID_raw, ";"))
    qIDs_vec <- unlist(strsplit(qID_raw, "/"))

    pthName_vec <- trimws(pthName_vec)
    qIDs_vec <- trimws(qIDs_vec)
    pthID_vec <- trimws(pthID_vec)

    # Retrieve information from databases
    if (query_type == "gene") {
      GeneIDs_vec <- qIDs_vec
      pathway_info <- get_pathway_and_gene_info(pthID_vec)
      pathwayDescription_vec <- pathway_info$pathwayDescription_vec
      pathwayReferencePMID_vec <- pathway_info$PMID_vec

      ##这是返回pathway里所有的gene，由于结果太多废弃了
      # GeneIDs_vec <- pathway_info$annotated_entrez_vec
      # GeneSymbols_vec <- pathway_info$annotated_symbol_vec  # 获取标准基因符号
      # GeneDescriptions_vec <- AnnotatedGenename_vec
      ## No need to convert id since id conversion has been completed before this step
      # gene_conversion <- convert_gene_identifiers(query_vec, orgdb)
      # GeneSymbols_vec <- gene_conversion$GeneSymbols_vec
      # GeneNames_vec <- gene_conversion$GeneNames_vec

      suppressMessages(GeneSymbols_vec <- AnnotationDbi::mapIds(orgdb, keys = GeneIDs_vec, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first") %>% unname())
      suppressMessages(GeneNames_vec <- AnnotationDbi::mapIds(orgdb, keys = GeneIDs_vec, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first") %>% unname())

      result_list[[mod_raw]] <- list(
        PathwayNames = pthName_vec,
        PathwayDescription = pathwayDescription_vec,
        PathwayReferencePMID = unique(pathwayReferencePMID_vec),
        GeneIDs = unique(GeneIDs_vec),
        GeneSymbols = unique(GeneSymbols_vec),
        GeneNames_vec = unique(GeneNames_vec)
      )

      # if (mod_raw %in% names(result_list)) {
      #   result_list[[mod_raw]]$PathwayNames <- c(result_list[[mod_raw]]$PathwayNames, pthName_vec)
      #   result_list[[mod_raw]]$PathwayDescription <- c(result_list[[mod_raw]]$PathwayDescription, pathwayDescription_vec)
      #   result_list[[mod_raw]]$PathwayReferencePMID <- unique(c(result_list[[mod_raw]]$PathwayReferencePMID, pathwayReferencePMID_vec))
      #   result_list[[mod_raw]]$GeneIDs <- unique(c(result_list[[mod_raw]]$GeneIDs, GeneIDs_vec))
      #   result_list[[mod_raw]]$GeneSymbols <- unique(c(result_list[[mod_raw]]$GeneSymbols, GeneSymbols_vec))
      #   result_list[[mod_raw]]$GeneNames_vec <- unique(c(result_list[[mod_raw]]$GeneNames_vec, GeneNames_vec))
      # }
    } else if (query_type == "metabolite") {
      MetIDs_vec <- qIDs_vec
      pathway_info <- get_pathway_and_metabolite_info(pthID_vec)
      pathwayDescription_vec = pathway_info$pathwayDescription_vec
      pathwayReferencePMID_vec = pathway_info$PMID_vec


      MetNames_vec <- get_metabolite_name(MetIDs_vec)

      result_list[[mod_raw]] <- list(
        PathwayNames = pthName_vec,
        PathwayDescription = pathwayDescription_vec,
        PathwayReferencePMID = unique(pathwayReferencePMID_vec),
        MetIDs = unique(MetIDs_vec),
        MetNames_vec = unique(MetNames_vec)
      )
    }
  }

  return(result_list)
}

#' Get Pathway and Gene Information
#'
#' Retrieves detailed information about pathways and their associated genes from different databases
#' (GO, KEGG, Reactome) based on provided pathway IDs.
#'
#' @param pathwayID_vec A character vector of pathway IDs. IDs should follow the format conventions
#'        of their respective databases: GO IDs should start with "GO:", KEGG IDs should match the
#'        pattern "hsa\\d+", and Reactome IDs should match the pattern "R-HSA-\\d+".
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{pathwayDescription_vec}: A character vector of pathway descriptions.
#'   \item \code{PMID_vec}: A character vector of PMIDs for pathway references.
#' }
#'
#' @examples
#' \dontrun{
#' # Get information for a mix of GO, KEGG, and Reactome pathways
#' pathway_ids <- c("GO:0006915", "hsa04210", "R-HSA-109581")
#' pathway_info <- get_pathway_and_gene_info(pathway_ids)
#'
#' # Access pathway descriptions
#' print(pathway_info$pathwayDescription_vec)
#'
#' # Access reference PMIDs
#' print(pathway_info$PMID_vec)
#' }
#'
#' @details
#' This function identifies the database source of each pathway ID and calls the appropriate
#' specialized function to retrieve information:
#' \itemize{
#'   \item GO pathway information is retrieved using \code{get_go_term_info}.
#'   \item KEGG pathway information is retrieved using \code{get_kegg_info}.
#'   \item Reactome pathway information is retrieved using \code{get_reactome_info}.
#' }
#' The function merges information from all sources into a single result.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_pathway_and_gene_info <- function(pathwayID_vec) {
  # Identify different types of pathway ID
  go_ids <- pathwayID_vec[grepl("^GO:", pathwayID_vec)]
  kegg_ids <- pathwayID_vec[grepl("^[a-zA-Z]+\\d+$", pathwayID_vec)]
  reactome_ids <- pathwayID_vec[grepl("^R-[A-Z]{3,4}-\\d+$", pathwayID_vec)]

  # Initialize the structure of retrieved data
  pathwayDescription_vec <- c()
  PMID_vec <- c()
  # annotated_entrez_vec <- c()
  # annotated_symbol_vec <- c()
  # annotated_queryname_vec <- c()

  # Retrieve GO info
  if (length(go_ids) > 0) {
    go_info <- get_go_term_info(go_ids)
    for (x in go_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
      # annotated_entrez_vec <- c(annotated_entrez_vec, x$annotated_entrez)
      # annotated_symbol_vec <- c(annotated_symbol_vec, x$annotated_symbol)
      # annotated_queryname_vec <- c(annotated_queryname_vec, x$annotated_genename)
    }
  }

  # Retrieve KEGG info
  if (length(kegg_ids) > 0) {
    kegg_info <- get_kegg_info(kegg_ids)
    for (x in kegg_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
      # annotated_entrez_vec <- c(annotated_entrez_vec, x$annotated_entrez)
      # annotated_symbol_vec <- c(annotated_symbol_vec, x$annotated_symbol)
      # annotated_queryname_vec <- c(annotated_queryname_vec, x$annotated_genename)
    }
  }

  # Retrieve Reactome info
  if (length(reactome_ids) > 0) {
    reactome_info <- get_reactome_info(reactome_ids)
    for (x in reactome_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
      # annotated_entrez_vec <- c(annotated_entrez_vec, x$annotated_entrez)
      # annotated_symbol_vec <- c(annotated_symbol_vec, x$annotated_symbol)
      # annotated_queryname_vec <- c(annotated_queryname_vec, x$annotated_genename)
    }
  }

  # **返回所有合并后的向量**
  return(list(
    pathwayDescription_vec = pathwayDescription_vec,
    PMID_vec = PMID_vec
    # annotated_entrez_vec = annotated_entrez_vec,
    # annotated_symbol_vec = annotated_symbol_vec,
    # annotated_queryname_vec = annotated_queryname_vec
  ))
}

#' Get Pathway and Metabolite Information
#'
#' Retrieves detailed information about pathways and their associated metabolites from different databases
#' (KEGG, SMPDB) based on provided pathway IDs.
#'
#' @param pathwayID_vec A character vector of pathway IDs. IDs should follow the format conventions
#'        of their respective databases.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{pathwayDescription_vec}: A character vector of pathway descriptions.
#'   \item \code{PMID_vec}: A character vector of PMIDs for pathway references.
#' }
#'
#' @examples
#' \dontrun{
#' # Get information for a mix of KEGG and SMPDB pathways
#' pathway_ids <- c("hsa00010", "SMP00001")
#' pathway_info <- get_pathway_and_metabolite_info(pathway_ids)
#'
#' # Access pathway descriptions
#' print(pathway_info$pathwayDescription_vec)
#'
#' # Access reference PMIDs
#' print(pathway_info$PMID_vec)
#' }
#'
#' @details
#' This function identifies the database source of each pathway ID and calls the appropriate
#' specialized function to retrieve information:
#' \itemize{
#'   \item KEGG pathway information is retrieved using \code{get_metkegg_info}.
#'   \item SMPDB pathway information is retrieved using \code{get_smpdb_info}.
#' }
#' The function merges information from all sources into a single result.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal

get_pathway_and_metabolite_info <- function(pathwayID_vec) {
  # Identify different types of pathway ID
  # kegg_ids <- pathwayID_vec[grepl("^[a-zA-Z]+\\d+$", pathwayID_vec)]
  smpdb_ids <- pathwayID_vec[grepl("^SMP\\d+$", pathwayID_vec)]
  kegg_ids <- pathwayID_vec[!(pathwayID_vec %in% smpdb_ids)]

  # Initialize the structure of retrieved data
  pathwayDescription_vec <- c()
  PMID_vec <- c()

  # Retrieve KEGG info
  if (length(kegg_ids) > 0) {
    metkegg_info <- get_metkegg_info(kegg_ids)
    for (x in metkegg_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
    }
  }

  # Retrieve SMPDB info
  if (length(smpdb_ids) > 0) {
    smpdb_info <- get_smpdb_info(smpdb_ids)
    for (x in smpdb_info) {
      pathwayDescription_vec <- c(pathwayDescription_vec, x$term_definition)
      PMID_vec <- c(PMID_vec, x$PMID)
    }
  }

  # **返回所有合并后的向量**
  return(list(
    pathwayDescription_vec = pathwayDescription_vec,
    PMID_vec = PMID_vec
  ))
}

get_metabolite_name <- function(MetIDs_vec) {
  hmdb_MetNames_vec <-
    metpath::hmdb_compound_database@spectra.info %>%
    dplyr::filter(HMDB.ID %in% MetIDs_vec) %>%
    dplyr::pull(Compound.name)

  # kegg_MetNames_vec <-
  #   MetIDs_vec[!grepl("^HMDB", MetIDs_vec)] |>
  #   sapply(function(x){
  #     compound_info <- KEGGREST::keggGet(x)
  #     compound_name <- gsub(";", "", compound_info[[1]]$NAME[1])
  #   })
  kegg_MetNames_vec <-
    MetIDs_vec[!grepl("^HMDB", MetIDs_vec)] |>
    sapply(function(x) {
      # Skip if ID is "NA" (the string)
      if (is.na(x) || x == "NA") {
        return(NA)
      }

      # Add delay between requests (1-2 seconds recommended)
      Sys.sleep(1.5)

      # Try-catch for error handling
      tryCatch({
        compound_info <- KEGGREST::keggGet(x)
        compound_name <- gsub(";", "", compound_info[[1]]$NAME[1])
        return(compound_name)
      }, error = function(e) {
        cat("Error for ID", x, ":", e$message, "\n")
        return(NA)  # Return NA for failed requests
      })
    })

  MetNames_vec <- c(hmdb_MetNames_vec, unname(kegg_MetNames_vec))

  # return(MetNames_vec)
  return(unlist(MetNames_vec, use.names = FALSE))
}

#从不同的geneid获取symbol和description
# convert_gene_identifiers <- function(gene_vec, orgdb) {
#   gene_symbols <- c()
#   gene_descriptions <- c()
#
#   for (gene in gene_vec) {
#     if (grepl("^ENSG", gene)) {  # ENSG00000163513 -> ENSEMBL
#       symbol <- mapIds(orgdb, keys = gene, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
#       description <- mapIds(orgdb, keys = gene, column = "GENENAME", keytype = "ENSEMBL", multiVals = "first")
#     } else if (grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$", gene)) {  # P08069 -> UNIPROT
#       symbol <- mapIds(orgdb, keys = gene, column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
#       description <- mapIds(orgdb, keys = gene, column = "GENENAME", keytype = "UNIPROT", multiVals = "first")
#     } else if (grepl("^[0-9]+$", gene)) {  # 3956 -> ENTREZ ID
#       symbol <- mapIds(orgdb, keys = gene, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
#       description <- mapIds(orgdb, keys = gene, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
#     } else {
#       symbol <- gene
#       description <- gene
#     }
#
#     gene_symbols <- c(gene_symbols, ifelse(is.na(symbol), gene, symbol))
#     gene_descriptions <- c(gene_descriptions, ifelse(is.na(description), gene, description))
#   }
#
#   return(list(GeneSymbols_vec = gene_symbols, GeneDescriptions_vec = gene_descriptions))
# }

