# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# load("data/enriched_pathways.rda")
# object <- enriched_pathways
# library(dplyr)
# library(purrr)
# library(httr2)
# library(org.Hs.eg.db)
# library(reactome.db)
# library(AnnotationDbi)
# library(rvest)
# library(rentrez)
# library(KEGGREST)
# library(rbioapi)
# library(httr)
# library(httr2)
# library(curl)
# library(jsonlite)
# library(rtiktoken)

# go_line <- sample(1:nrow(object@enrichment_go_result@result), 10)
# kegg_line <- sample(1:nrow(object@enrichment_kegg_result@result), 10)
# reactome_line <- sample(1:nrow(object@enrichment_reactome_result@result), 10)
#
# demo_data_go <- object@enrichment_go_result@result[go_line,]
# demo_data_kegg <- object@enrichment_kegg_result@result[kegg_line,]
# demo_data_reactome <- object@enrichment_reactome_result@result[reactome_line,]
#
# demo_go_info <- get_go_info(go_ids = demo_data_go$ID)
# go_text_info <- get_ttl_abstr(data = demo_go_info)
# demo_kegg_info <- get_kegg_info(kegg_ids = demo_data_kegg$ID)
# kegg_text_info <- get_ttl_abstr(data = demo_kegg_info)
# demo_reactome_info <- get_reactome_info(reactome_ids = demo_data_reactome$ID)
# reactome_text_info <- get_ttl_abstr(data = demo_reactome_info)
# all_text_info <- c(go_text_info, kegg_text_info, reactome_text_info)

# go_info <- get_go_info(go_ids = object@enrichment_go_result@result$ID, include_gene_name = FALSE)
# kegg_info <- get_kegg_info(kegg_ids = object@enrichment_kegg_result@result$ID, include_gene_name = FALSE)
# reactome_info <- get_reactome_info(reactome_ids = object@enrichment_reactome_result@result$ID, include_gene_name = FALSE)
# all_text_info <- c(go_info, kegg_info, reactome_info)
# all_combined_info <- combine_info(info = all_text_info, include_gene_name = FALSE)

# go_info <- get_go_info(go_ids = object@enrichment_go_result@result$ID, include_gene_name = TRUE)
# kegg_info <- get_kegg_info(kegg_ids = object@enrichment_kegg_result@result$ID, include_gene_name = TRUE)
# reactome_info <- get_reactome_info(reactome_ids = object@enrichment_reactome_result@result$ID, include_gene_name = TRUE)
# all_text_info <- c(go_info, kegg_info, reactome_info)
# all_combined_info <- combine_info(info = all_text_info, include_gene_name = TRUE)

# For genes
# openai_semantic_sim_matrix <-
#   get_bioembedsim(object = enriched_pathways,
#                   api_provider = "openai",
#                   text_embedding_model = "text-embedding-3-small",
#                   api_key = api_key,
#                   save_to_local = TRUE,
#                   path = "~/Desktop/result")
# gemini_semantic_sim_matrix <- get_bioembedsim(object = object, api_provider = "gemini",  text_embedding_model = "text-embedding-004", api_key = api_key)

# For metabolites
# openai_sim_matrix_met <-
#   get_bioembedsim(object = enriched_pathways,
#                   api_provider = "openai",
#                   database = c("hmdb", "kegg"),
#                   text_embedding_model = "text-embedding-3-small",
#                   api_key = api_key,
#                   save_to_local = FALSE)


#' Calculate Biological Pathway Similarity Based on BioText Embeddings (biotext embedding similarity)
#'
#' @description
#' This function calculates similarity between biological pathways by using text embedding
#' models to vectorize pathway information and compute cosine similarity between them.
#' It extracts text information from GO, KEGG, HMDB and/or Reactome databases based on
#' previously performed enrichment analyses, and uses API providers like OpenAI or Gemini
#' to generate embeddings.
#'
#' @param object An object of class "functional_module", typically a result from enrich_pathway function.
#' @param api_provider Character string specifying the API provider for text embeddings.
#'   Options are "openai" or "gemini".
#' @param text_embedding_model Character string specifying the embedding model to use
#'   (e.g., "text-embedding-3-small" for OpenAI or "text-embedding-004" for Gemini)
#' @param api_key Character string of the API key for the specified provider
#' @param include_gene_name Logical indicating whether to include annotated gene names
#'   in the pathway descriptions (default: FALSE)
#' @param database Character vector of databases to include. Options are "go", "kegg", "hmdb",
#'   and/or "reactome". Multiple selections allowed.
#' @param p.adjust.cutoff.go Numeric cutoff for adjusted p-value for GO terms (default: 0.05)
#' @param p.adjust.cutoff.kegg Numeric cutoff for adjusted p-value for KEGG pathways (default: 0.05)
#' @param p.adjust.cutoff.reactome Numeric cutoff for adjusted p-value for Reactome pathways (default: 0.05)
#' @param p.adjust.cutoff.hmdb Numeric cutoff for adjusted p-value for HMDB pathways from metabolite enrichment analysis result (default: 0.05)
#' @param p.adjust.cutoff.metkegg Numeric cutoff for adjusted p-value for KEGG pathways from metabolite enrichment analysis result (default: 0.05)
#' @param count.cutoff.go Minimum gene count cutoff for GO terms (default: 5)
#' @param count.cutoff.kegg Minimum gene count cutoff for KEGG pathways (default: 5)
#' @param count.cutoff.reactome Minimum gene count cutoff for Reactome pathways (default: 5)
#' @param count.cutoff.hmdb Minimum metabolite count cutoff for HMDB pathways (default: 5)
#' @param count.cutoff.metkegg Minimum metabolite count cutoff for KEGG pathways (default: 5)
#' @param save_to_local Logical. Whether to save the resulting data to local files. Default is `FALSE`.
#' @param path Character. The directory path where intermediate results will be saved, if `save_to_local = TRUE`. Default is "result".
#'
#' @return A list containing a matrix of pairwise cosine similarity values between pathways and
#' an object of class "functional_module" with updated parameter information.
#'
#' @details
#' The function works in three main steps:
#' 1. Extract text information from the specified databases (GO, KEGG, and/or Reactome)
#' 2. Generate embeddings using the specified API provider and model
#' 3. Calculate pairwise cosine similarity between pathway embeddings
#'
#' For each database, the function filters pathways based on adjusted p-value and gene count
#' cutoffs. For Gemini embeddings with long texts, it may use LLM-generated summaries to fit
#' within token limits.
#'
#' @examples
#' \dontrun{
#' # For OpenAI embeddings
#' sim_matrix <- get_bioembedsim(
#'   object = my_enrichment_object,
#'   api_provider = "openai",
#'   text_embedding_model = "text-embedding-3-small",
#'   api_key = "your_openai_key",
#'   database = c("go", "kegg"),
#'   p.adjust.cutoff.go = 0.01,
#'   count.cutoff.go = 10
#' )
#'
#' # For Gemini embeddings
#' sim_matrix <- get_bioembedsim(
#'   object = my_enrichment_object,
#'   api_provider = "gemini",
#'   text_embedding_model = "text-embedding-004",
#'   api_key = "your_gemini_key",
#'   database = "reactome",
#'   include_gene_name = TRUE
#' )
#' }
#'
#' @importFrom dplyr filter pull
#' @importFrom httr2 request req_auth_bearer_token req_headers req_body_json req_retry req_perform resp_body_json
#' @importFrom KEGGREST keggGet
#' @importFrom purrr map map_df
#' @importFrom rtiktoken get_token_count
#'
#' @export

get_bioembedsim <-
  function(object,
           api_provider = c("openai", "gemini"),
           text_embedding_model = NULL,
           api_key = NULL,
           include_gene_name = FALSE,
           database = c("go", "kegg", "reactome", "hmdb"),
           p.adjust.cutoff.go = 0.05,
           p.adjust.cutoff.kegg = 0.05,
           p.adjust.cutoff.reactome = 0.05,
           p.adjust.cutoff.hmdb = 0.05,
           p.adjust.cutoff.metkegg = 0.05,
           count.cutoff.go = 5,
           count.cutoff.kegg = 5,
           count.cutoff.reactome = 5,
           count.cutoff.hmdb = 5,
           count.cutoff.metkegg = 5,
           save_to_local = FALSE,
           path = "result") {

    if (missing(object)) {
      stop("object is required")
    }

    if (!is(object, "functional_module")) {
      stop("object must be result from enrich_pathway function")
    }

    query_type <- object@process_info$enrich_pathway@parameter$query_type

    if (missing(api_provider)) {
      stop("api_provider is required.")
    }
    api_provider <- match.arg(api_provider)

    if (missing(text_embedding_model)) {
      stop("text_embedding_model is required.")
    }

    if (missing(api_key)) {
      stop("api_key is required.")
    }

    database <- match.arg(database, several.ok = TRUE)

    ## Collect text information from databases
    all_text_info <- list()

    if (query_type == "gene") {
      if ("go" %in% database) {
        if (is.null(object@enrichment_go_result)) {
          stop("Please perform pathway enrichment based on GO database at first.")
        } else {
          go_info <-
            object@enrichment_go_result@result %>%
            dplyr::filter(p.adjust < p.adjust.cutoff.go) %>%
            dplyr::filter(Count > count.cutoff.go) %>%
            dplyr::pull(ID) %>%
            get_go_info(include_gene_name = include_gene_name)
          all_text_info <- c(all_text_info, go_info)
        }
      }

      if ("kegg" %in% database) {
        if (is.null(object@enrichment_kegg_result)) {
          stop("Please perform pathway enrichment based on KEGG database at first.")
        } else {
          kegg_info <-
            object@enrichment_kegg_result@result %>%
            dplyr::filter(p.adjust < p.adjust.cutoff.kegg) %>%
            dplyr::filter(Count > count.cutoff.kegg) %>%
            dplyr::pull(ID) %>%
            get_kegg_info(include_gene_name = include_gene_name)
          all_text_info <- c(all_text_info, kegg_info)
        }
      }

      if ("reactome" %in% database) {
        if (is.null(object@enrichment_reactome_result)) {
          stop("Please perform pathway enrichment based on Reactome database at first.")
        } else {
          reactome_info <-
            object@enrichment_reactome_result@result %>%
            dplyr::filter(p.adjust < p.adjust.cutoff.reactome) %>%
            dplyr::filter(Count > count.cutoff.reactome) %>%
            dplyr::pull(ID) %>%
            get_reactome_info(include_gene_name = include_gene_name)
          all_text_info <- c(all_text_info, reactome_info)
        }
      }
    } else if (query_type == "metabolite") {
      if ("kegg" %in% database) {
        if (is.null(object@enrichment_metkegg_result)) {
          stop("Please perform pathway enrichment based on KEGG database at first.")
        } else {
          metkegg_info <- list()
          metkegg_enrichment_result <-
            object@enrichment_metkegg_result@result %>%
            dplyr::filter(p_value_adjust < p.adjust.cutoff.metkegg) %>%
            dplyr::filter(mapped_number > count.cutoff.metkegg)
          for (i in 1:nrow(metkegg_enrichment_result)) {
            entry <- metkegg_enrichment_result[i,]
            all_info <- list(
              "id" = entry$pathway_id,
              "term_name" = entry$pathway_name,
              "term_definition" = entry$describtion
            )
            metkegg_info <- c(metkegg_info, list(all_info))
          }
          all_text_info <- c(all_text_info, metkegg_info)
        }
      }

      if ("hmdb" %in% database) {
        if (is.null(object@enrichment_hmdb_result)) {
          stop("Please perform pathway enrichment based on KEGG database at first.")
        } else {
          hmdb_info <- list()
          hmdb_enrichment_result <-
            object@enrichment_hmdb_result@result %>%
            dplyr::filter(p_value_adjust < p.adjust.cutoff.hmdb) %>%
            dplyr::filter(mapped_number > count.cutoff.hmdb)
          for (i in 1:nrow(hmdb_enrichment_result)) {
            entry <- hmdb_enrichment_result[i,]
            all_info <- list(
              "id" = entry$pathway_id,
              "term_name" = entry$pathway_name,
              "term_definition" = entry$describtion
            )
            hmdb_info <- c(hmdb_info, list(all_info))
          }
          all_text_info <- c(all_text_info, hmdb_info)
        }
      }
    }

    all_combined_info <- combine_info(info = all_text_info, include_gene_name = include_gene_name)


    ## Get embedding matrix
    embedding_matrix <- get_embedding_matrix(text = all_combined_info,
                                             include_gene_name = include_gene_name,
                                             api_provider = api_provider,
                                             text_embedding_model = text_embedding_model,
                                             api_key = api_key)


    ## Calculate pairwise cosine similarity
    sim_matrix <- calculate_cosine_sim(m = embedding_matrix)

    ## Store parameters
    if (query_type == "gene") {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "get_bioembedsim()",
        parameter = list(
          query_type = "gene",
          p.adjust.cutoff.go = p.adjust.cutoff.go,
          p.adjust.cutoff.kegg = p.adjust.cutoff.kegg,
          p.adjust.cutoff.reactome = p.adjust.cutoff.reactome,
          count.cutoff.go = count.cutoff.go,
          count.cutoff.kegg = count.cutoff.kegg,
          count.cutoff.reactome = count.cutoff.reactome
        ),
        time = Sys.time()
      )
    } else if (query_type == "metabolite") {
      parameter = new(
        Class = "tidymass_parameter",
        pacakge_name = "mapa",
        function_name = "get_bioembedsim()",
        parameter = list(
          query_type = "metabolite",
          p.adjust.cutoff.hmdb = p.adjust.cutoff.hmdb,
          p.adjust.cutoff.metkegg = p.adjust.cutoff.metkegg,
          count.cutoff.hmdb = count.cutoff.hmdb,
          count.cutoff.metkegg = count.cutoff.metkegg
        ),
        time = Sys.time()
      )
    }

    process_info <-
      slot(object, "process_info")

    process_info$merge_pathways <-
      parameter

    slot(object, "process_info") <-
      process_info

    message("Biotext embedding and similarity calculation finished")

    ## Save similarity matrix as intermediate data
    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      dir.create(
        file.path(path, "intermediate_data"),
        showWarnings = FALSE,
        recursive = TRUE
      )
      save(sim_matrix, file = file.path(path, "intermediate_data/sim_matrix.RData"))
    }

    return(list(sim_matrix = sim_matrix, enriched_pathway = object))
}

# Step1: Extract text info =====
## 1.1 Extract GO info (go_id, annotated gene IDs, name, definition, PMID)====
### test: GO:0031954 -> 0PMID; GO:1902074 -> 1 PMID; GO:0000001 -> 2PMID

get_go_info <- function(go_ids, include_gene_name = FALSE) {
  #### Generate GO terms and annotated Entrez Gene identifiers dict
  go2egs <- NULL
  eg2genename <- NULL

  # Only retrieve gene-related information if include_gene_name is TRUE
  if (include_gene_name) {
    go2egs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, columns = "ENTREZID", keytype = "GOALL"))
    #### Generate entrez ID to gene name dict
    eg2genename <-
      suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = go2egs$ENTREZID, columns = "GENENAME")) %>%
      dplyr::distinct()
  }

  chunk_size <- 100
  chunks <- split(go_ids, ceiling(seq_along(go_ids) / chunk_size))
  go_info <- list()
  for (i in 1:length(chunks)) {
    sub_go_info <-
      quickgo_api(go_ids = chunks[[i]]) %>%
      purrr::map(
        function(x) {
          # Initialize empty gene name string
          go2genename <- ""

          # Only process gene names if include_gene_name is TRUE
          if (include_gene_name) {
            # Get GENENAME: GO ID -> ENTREZID -> GENENAME
            go2genename <- go2egs %>%
              dplyr::filter(GOALL == x$id) %>%
              dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
              dplyr::left_join(eg2genename, by = "ENTREZID") %>%
              dplyr::pull(GENENAME) %>%
              paste0(collapse = ", ")
          }

          # # Get PMID
          # if ("xrefs" %in% names(x$definition)){
          #   # Extract 'dbId' from each element in the list
          #   pmid_list <- lapply(x$definition$xrefs,
          #                       function(r) {
          #                         if(r$dbCode == "PMID") {
          #                           pmid <- r$dbId
          #                         } else {
          #                           pmid <- NULL
          #                         }})
          #   # Convert the list of dbId values into a string
          #   pmid <- paste(unlist(pmid_list), collapse = ",")
          # } else {
          #   pmid <- ""
          # }

          # Collect all info
          all_info <- list(
            "id" = x$id,
            "term_name" = x$name,
            "term_definition" = x$definition$text
          )

          # Only add annotated_genename field if include_gene_name is TRUE
          if (include_gene_name) {
            all_info$annotated_genename <- go2genename
          }

          return(all_info)
        }
      )
    go_info <- c(go_info, sub_go_info)
  }
  return(go_info)
}

# Get core information about a list of terms based on their ids
quickgo_api <- function(go_ids) {
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste(go_ids, collapse = ","))
  tryCatch(
    expr = {
      # Create a request object
      req <- httr2::request(url)

      resp <- req %>%
        req_headers("Accept" = "application/json") %>%
        req_retry(max_tries = 3,
                  max_seconds = 60,
                  # Condition for when to retry
                  after = \(resp) is.null(resp) && resp$status_code != 200
        ) %>%
        req_perform()

      # Parse json to get a list
      info <- resp_body_json(resp)$result
    },
    error = function(e) {
      warning("Failed to get info from QuickGO after 3 retries:", e$message)
      NULL
    })
}

## 1.2 Extract KEGG info (pathway_id, name, definition(description), PMID) ====
get_kegg_info <- function(kegg_ids, include_gene_name = FALSE){
  chunk_size <- 10
  chunks <- split(kegg_ids, ceiling(seq_along(kegg_ids) / chunk_size))
  kegg_info <- list()

  for (i in 1:length(chunks)) {
    sub_kegg_info <-
      tryCatch({
        # Try to process the chunk as a batch first
        KEGGREST::keggGet(dbentries = chunks[[i]]) %>%
          purrr::map(
            function(x) {
              # Initialize kegg2genename as empty
              kegg2genename <- NA

              # Only get annotated Entrez ID if include_gene_name is TRUE
              if (include_gene_name && "GENE" %in% names(x)) {
                kegg2genename <-
                  x$GENE[seq(2, length(x$GENE), 2)] %>%
                  stringr::str_match(pattern = ";\\s*(.*?)\\s*\\[") %>%
                  as.data.frame() %>%
                  dplyr::pull(V2) %>%
                  paste0(collapse = ",")
              }

              # # Get PMID
              # if ("REFERENCE" %in% names(x)) {
              #   pmid_list <- lapply(x$REFERENCE,
              #                       function(r) {
              #                         if (grepl("^PMID", r$REFERENCE)) {
              #                           pmid <- gsub("[^:]+:", "", r$REFERENCE)
              #                         } else {
              #                           pmid <- ""
              #                         }
              #                       })
              #   pmid <- paste(unlist(pmid_list), collapse = ",")
              # } else {
              #   pmid <- ""
              # }

              # Collect base info (always included)
              all_info <- list(
                "id" = unname(x$ENTRY),
                "term_name" = sub(" - Homo sapiens \\(human\\)$", "", x$NAME),
                "term_definition" = paste(x$DESCRIPTION, collapse = " ")
              )

              # Only add gene name information if requested
              if (include_gene_name) {
                all_info$annotated_genename <- kegg2genename
              }

              return(all_info)
            })
      }, error = function(e) {
        # If batch processing fails, process each ID individually
        message("Batch processing failed, trying individual IDs: ", e$message)

        # Process each ID in the current chunk individually
        individual_results <- list()
        for (kegg_id in chunks[[i]]) {
          single_result <- tryCatch({
            # Get a single KEGG entry
            entry <- KEGGREST::keggGet(dbentries = kegg_id)[[1]]

            # Initialize kegg2genename as empty
            kegg2genename <- NA

            # Only get annotated Entrez ID if include_gene_name is TRUE
            if (include_gene_name && "GENE" %in% names(entry)) {
              kegg2genename <-
                entry$GENE[seq(2, length(entry$GENE), 2)] %>%
                stringr::str_match(pattern = ";\\s*(.*?)\\s*\\[") %>%
                as.data.frame() %>%
                dplyr::pull(V2) %>%
                paste0(collapse = ",")
            }

            # Collect base info (always included)
            all_info <- list(
              "id" = unname(entry$ENTRY),
              "term_name" = sub(" - Homo sapiens \\(human\\)$", "", entry$NAME),
              "term_definition" = paste(entry$DESCRIPTION, collapse = " ")
            )

            # Only add gene name information if requested
            if (include_gene_name) {
              all_info$annotated_genename <- kegg2genename
            }

            all_info
          }, error = function(e2) {
            message("  Failed to retrieve KEGG ID: ", kegg_id, " (", e2$message, ")")
            return(NULL)
          })

          # Add non-NULL results to the list
          if (!is.null(single_result)) {
            individual_results <- c(individual_results, list(single_result))
          }

          # Brief pause to avoid overwhelming the API
          Sys.sleep(0.1)
        }

        return(individual_results)
      })

    kegg_info <- c(kegg_info, sub_kegg_info)
  }
  return(kegg_info)
}

## 1.3 Extract Reactome info (pathway_id, name, definition, PMID) ====
get_reactome_info <- function(reactome_ids, include_gene_name = FALSE) {
  #### Initialize gene-related variables
  reactome2egs <- NULL
  eg2genename <- NULL

  #### Only retrieve gene-related information if include_gene_name is TRUE
  if (include_gene_name) {
    #### Generate Reactome terms and annotated Entrez Gene identifiers dict
    suppressMessages(reactome2egs <- AnnotationDbi::select(reactome.db, keys = reactome_ids, columns = "ENTREZID", keytype = "PATHID"))
    #### Generate entrez ID to gene name dict
    eg2genename <-
      suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = reactome2egs$ENTREZID, columns = "GENENAME")) %>%
      dplyr::distinct()
  }

  chunk_size <- 20
  chunks <- split(reactome_ids, ceiling(seq_along(reactome_ids) / chunk_size))
  reactome_info <- list()

  for (i in 1:length(chunks)) {
    reactome_api_result <- suppressMessages(rbioapi::rba_reactome_query(ids = chunks[[i]]))
    if (length(chunks[[i]]) == 1) {reactome_api_result <- list(reactome_api_result)}

    sub_reactome_info <-
      reactome_api_result %>%
      purrr::map(
        function(x) {
          # Initialize empty gene name string
          reactome2genename <- ""

          # Only process gene names if include_gene_name is TRUE
          if (include_gene_name) {
            # Get GENENAME: Reactome ID -> ENTREZID -> GENENAME
            reactome2genename <- reactome2egs %>%
              dplyr::filter(PATHID == x$stId) %>%
              dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
              dplyr::left_join(eg2genename, by = "ENTREZID") %>%
              dplyr::pull(GENENAME) %>%
              paste0(collapse = ", ")
          }

          # # Get PMID
          # if ("literatureReference" %in% names(x)) {
          #   pmid_list <- lapply(x$literatureReference,
          #                       function(r) {
          #                         if ("pubMedIdentifier" %in% names(r)) {
          #                           pmid <- r$pubMedIdentifier
          #                         } else {
          #                           pmid <- ""
          #                         }
          #                       })
          #   pmid <- paste(unlist(pmid_list), collapse = ",")
          # } else {
          #   pmid <- ""
          # }

          # Collect base info (always included)
          all_info <- list(
            "id" = x$stId,
            "term_name" = x$displayName,
            "term_definition" = gsub("(<BR>|<br>)", "", x$summation[[1]]$text)
          )

          # Only add gene name information if requested
          if (include_gene_name) {
            all_info$annotated_genename <- reactome2genename
          }

          return(all_info)
        }
      )
    reactome_info <- c(reactome_info, sub_reactome_info)
  }
  return(reactome_info)
}

## 1.4 Combine information (term name, term description, gene names) into a string ====
combine_info <- function(info, include_gene_name = FALSE) {
  info %>%
    purrr::map(
      function(x) {
        if (include_gene_name) {
          text_info <- sprintf("Pathway name: %s\nDefinition: %s\nAnnotated gene names: %s", x$term_name, x$term_definition, x$annotated_genename)
        } else {
          text_info <- sprintf("%s: %s", x$term_name, x$term_definition)
        }

        text <- list(
          "id" = x$id,
          "text_info" = text_info
        )

        return(text)
      }
    )
}

# # (Remove)Step2: Retrieve title and abstract from PubMed using PMID and get data from embedding ====
# get_ttl_abstr <- function(data) {
#   data %>%
#     purrr::map(
#       function(x) {
#         # Format text info
#         if (nchar(x$PMID) == 0){
#           text_info <- sprintf("ID: %s\nName: %s\nDefinition: %s\nAnnotated Gene names: %s", x$id, x$term_name, x$term_definition, x$annotated_genename)
#
#         } else {
#           ttls_abstrs <-
#             stringr::str_split(x$PMID, ",")[[1]] %>%
#             purrr::map_chr(~ retrieve_abstr_ttl(pmid = .x)) %>%
#             paste(., collapse = "\n")
#           text_info <- sprintf("ID: %s\nName: %s\nDefinition: %s\nAnnotated Gene names: %s\nReference:\n%s", x$id, x$term_name, x$term_definition, x$annotated_genename, ttls_abstrs)
#         }
#
#         text <- list(
#           "id" = x$id,
#           "text_info" = text_info
#         )
#       })
# }
#
# retrieve_abstr_ttl <- function(pmid, retries = 2, pause = 1) {
#   attempt <- 1
#   while (attempt <= retries) {
#     tryCatch({
#       # Fetch html data containing title and abstract
#       html_result <-
#         rentrez::entrez_fetch(
#           db = "pubmed",
#           id = pmid,
#           rettype = "abstract",
#           retmode = "html"
#         )
#
#       # Parse html
#       parsed_html <- rvest::read_html(html_result)
#
#       # Remove Review articles
#       publication_type <- parsed_html %>%
#         rvest::html_nodes("publicationtype") %>%
#         rvest::html_text()
#       if ("Review" %in% publication_type) {
#         return("")
#       }
#
#       # Get title
#       title <- parsed_html %>%
#         rvest::html_node("articletitle") %>%
#         rvest::html_text(trim = TRUE)
#
#       # Get abstract
#       abstract <- parsed_html %>%
#         rvest::html_node("abstract") %>%
#         rvest::html_text(trim = TRUE)
#
#       if (is.na(abstract)) {
#         return("")
#       } else {
#         return(sprintf("Title: %s\nAbstract: %s", title, abstract))
#       }
#     }, error = function(e) {
#       warning(sprintf("Attempt %d failed for PubMed IDs: %s. Error: %s",
#                       attempt, paste(pmids, collapse = ", "), e))
#     })
#
#     Sys.sleep(pause)
#     attempt <- attempt + 1
#   }
#
#   warning(sprintf("All attempts failed for PubMed IDs: %s",
#                   paste(pmids, collapse = ", ")))
#   return(lapply(pmids, function(PID) {
#     list(
#       title = sprintf("Error: Failed to retrieve title for PubMed ID %s after %d attempts",
#                       PID, retries),
#       abstract = sprintf("Error: Failed to retrieve abstract for PubMed ID %s after %d attempts",
#                          PID, retries)
#     )
#   }))
# }




# Step2: Get embeddings ====
## Get the number of tokens in the text info ====
# # token_num <- rtiktoken::get_token_count(all_text_info$text_info, model = "text-embedding-3-small")
# # barplot(token_num, ylim = c(0, 8500))
# # abline(h = 8191, col = "red", lwd = 2)

# Test openai embedding model
# openai_embedding_model <- "text-embedding-3-small"
# M <- get_embedding_matrix(text = all_combined_info[c(1:10)],
#                           api_provider = "openai",
#                           text_embedding_model = openai_embedding_model,
#                           api_key = api_key)
#
# # Test google gemini embedding model
# gemini_embedding_model <- "text-embedding-004"
# M <- get_embedding_matrix(text = all_combined_info[c(1:10)],
#                           api_provider = "gemini",
#                           text_embedding_model = gemini_embedding_model,
#                           api_key = api_key)

get_embedding_matrix <-
  function(
    text,
    include_gene_name = FALSE,
    api_provider = c("openai", "gemini"),
    text_embedding_model,
    api_key) {

  if (missing(api_provider)) {
    stop("api_provider is required.")
  }
  api_provider <- match.arg(api_provider)

  if (missing(text_embedding_model)) {
    stop("text_embedding_model is required.")
  }

  if (missing(api_key)) {
    stop("api_key is required.")
  }

  if (api_provider == "openai") {
    embedding_matrix <-
      text %>%
      purrr::map_df(function(x) {
        token_num <-
          rtiktoken::get_token_count(text = x$text_info, model = "text-embedding-3-small")
        # max input for openai embedding model is 8191 (March 05, 2025)
        if (include_gene_name == TRUE & token_num > 8000) {
          x$text_info <- sub(pattern = "\nAnnotated gene names.*$", "", x$text_info)
          message(paste0("Text information of ", x$id, " exceeded token limit and annotated gene names were removed."))
        }
        embedding <-
          get_openai_embedding_internal(input_text = x$text_info, text_embedding_model = text_embedding_model, api_key = api_key) %>%
          t() %>%
          as.data.frame()
        rownames(embedding) <- x$id

        embedding
      }) %>%
      as.matrix()
  } else if (api_provider == "gemini") {
    message("Use the tokenizer(i.e. Cl100kBase) for OpenAI's embedding model to approximate the token number for Gemini model.")
    embedding_matrix <-
      text %>%
      purrr::map_df(function(x) {
        token_num <-
          rtiktoken::get_token_count(text = x$text_info, model = "text-embedding-3-small")
        if (token_num > 1500) {
          query <- paste0("According to the following description: ", x$text_info, " Please generate a short summary about 100 words to extract the core function information of the pathway described.")
          llm_summary <- gemini_api_call(
            input_text = query,
            api_key = api_key
          )
          # Considering the token limitation (2048 tokens) of text embedding model, use llm to summarize the pathway description
          message(paste0("For pathway ID: ", x$id, ", use the summary generated by gemini-1.5-flash to do text embedding since pathway information exceeded token limitation."))
          x$text_info <- paste0("Pathway name:", x$id, "\n", llm_summary)
        }

        embedding <-
          get_gemini_embedding_internal(input_text = x$text_info, text_embedding_model = text_embedding_model, api_key = api_key) %>%
          t() %>%
          as.data.frame()
        rownames(embedding) <- x$id

        embedding
      }) %>%
      as.matrix()
  }

  return(embedding_matrix)
}


get_openai_embedding_internal <-
  function(input_text, text_embedding_model = "text-embedding-3-small", api_key) {

    url <- "https://api.openai.com/v1/embeddings"

    # Body specifying model and text
    data <- list(
      model = text_embedding_model,  # Use the specified model or default
      input = input_text  # Input text
    )

    embedding <- tryCatch(
      expr = {
        # Create a request object
        req <- httr2::request(url)

        resp <- req %>%
          req_auth_bearer_token(token = api_key) %>%
          req_body_json(data = data) %>%
          req_retry(max_tries = 3,
                    max_seconds = 60,
                    after = \(resp) is.null(resp) && resp$status_code != 200 # Condition for when to retry
          ) %>%
          req_perform()

        embedding <-
          # Parsed JSON -> an embedding, a list
          resp_body_json(resp) %>%
          {.$data[[1]]$embedding} %>%
          unlist()
      },
      error = function(e) {
        warning("Failed to get embedding after 3 retries:", e$message)
        NULL
      })

    return(embedding)
}

get_gemini_embedding_internal <-
  function(input_text, text_embedding_model = "text-embedding-004", api_key) {

  url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/", text_embedding_model, ":embedContent?key=", api_key)

  embedding <- tryCatch(
    expr = {
      # Create a request object
      req <- httr2::request(url)

      # Get response
      resp <- req %>%
        req_headers("Content-Type" = "application/json") %>%
        req_body_json(
          list(
            model = paste0("models/", text_embedding_model),
            content = list(parts = list(list(text = input_text)))
          )
        ) %>%
        req_retry(max_tries = 3,
                  max_seconds = 60,
                  after = \(resp) is.null(resp) && resp$status_code != 200 # Condition for when to retry
        ) %>%
        req_perform()

      # Parse the JSON response
      embedding <- resp %>%
        resp_body_json() %>%
        {.$embedding$values} %>%
        unlist()
    },
    error = function(e) {
      warning("Failed to get embedding after 3 retries:", e$message)
      NULL
    }
  )

  return(embedding)
}

gemini_api_call <-
  function(input_text,
           model = "gemini-1.5-flash",
           api_key,
           temperature = 1,
           topK = 40,
           topP = 0.95,
           maxOutputTokens = 8192,
           responseMimeType = "text/plain") {
    # Construct the full URL with the API key and model as query parameters
    url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/", model, ":generateContent?key=", api_key)

    # Create the request body
    request_body <- list(
      contents = list(
        list(
          role = "user",
          parts = list(list(text = input_text))
        )
      ),
      generationConfig = list(
        temperature = temperature,
        topK = topK,
        topP = topP,
        maxOutputTokens = maxOutputTokens,
        responseMimeType = responseMimeType
      )
    )

    # Create and send the request using httr2
    response <- request(url) %>%
      req_headers("Content-Type" = "application/json") %>%
      req_body_json(request_body) %>%
      req_perform()

    # Check for errors
    if (resp_is_error(response)) {
      stop(paste("API request failed with status:", resp_status(response)))
    }

    # Parse and return the JSON response
    result <- response %>% resp_body_json()

    # Extract generated text
    return(result$candidates[[1]]$content$parts[[1]]$text)
  }

# Step3: Calculate pairwise cosine similarity ====
calculate_cosine_sim <- function(m){
  dot_product <- m %*% t(m)
  norm_product <- sqrt(rowSums(m^2)) %*% t(sqrt(rowSums(m^2)))
  cosine_sim <- dot_product / norm_product
  return(cosine_sim)
}

