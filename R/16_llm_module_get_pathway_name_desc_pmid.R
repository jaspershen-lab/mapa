# library(dplyr)
# library(purrr)
# library(httr2)
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# library(reactome.db)
# library(rvest)
# library(rentrez)
# library(KEGGREST)
# library(rbioapi)
# library(httr)
# library(httr2)
# library(curl)
# library(jsonlite)
# library(rtiktoken)

#' Extract information from GO terms using GO IDs
#'
#' This internal function retrieves term name, definition, and PMIDs for given GO IDs
#' using the QuickGO API.
#'
#' @param go_ids A character vector of GO IDs (e.g., "GO:0008150")
#'
#' @return A list of GO term information, where each element contains:
#'   \item{id}{GO ID}
#'   \item{term_name}{Name of the GO term}
#'   \item{term_definition}{Definition of the GO term}
#'   \item{PMID}{Associated PubMed IDs}
#'
#' @importFrom purrr map
#' @importFrom httr2 request req_headers req_retry req_perform resp_body_json
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_go_term_info <- function(go_ids) {

  # # 选取一个已知的有效 GO ID
  # valid_go_id <- "GO:0008150"  # 该 GO ID 代表 "biological process"，通常有效
  #
  # # 记录是否需要删除该 GO ID
  # need_remove <- !(valid_go_id %in% go_ids)
  #
  # # 在 `go_ids` 里添加该有效 GO ID（如果原本不存在）
  # if (need_remove) {
  #   go_ids <- c(go_ids, valid_go_id)
  # }
  #
  # #### Get ENTIREZID: GO ID -> EntrezID
  # suppressMessages(go2egs <- AnnotationDbi::select(org.Hs.eg.db, keys = go_ids, columns = c("ENTREZID", "SYMBOL", "GENENAME"), keytype = "GO"))
  #
  # #### Generate Entrez ID to Gene Name and Symbol mapping
  # eg2genename <- go2egs %>%
  #   dplyr::distinct(ENTREZID, SYMBOL, GENENAME)

  chunk_size <- 100
  chunks <- split(go_ids, ceiling(seq_along(go_ids) / chunk_size))

  go_info <- list()

  for (i in 1:length(chunks)) {
    sub_go_info <-
      quickgo_api(go_ids = chunks[[i]]) %>%
      purrr::map(
        function(x) {

          # # Get ENTIREZID and GENENAME: GO ID -> ENTREZID -> SYMBOL & GENENAME
          # go2gene_info <- go2egs %>%
          #   dplyr::filter(GO == x$id) %>%
          #   dplyr::distinct(ENTREZID, SYMBOL, GENENAME, .keep_all = TRUE)
          #
          # # Extract ENTIREZID, SYMBOL, and GENENAME
          # go2entrez <- go2gene_info$ENTREZID  # 返回一个向量
          # go2symbol <- go2gene_info$SYMBOL    # 返回一个向量
          # go2genename <- go2gene_info$GENENAME  # 这里去掉 paste()，让它返回列表

          # Get PMID
          if ("xrefs" %in% names(x$definition)){
            # Extract 'dbId' from each element in the list
            pmid_list <- lapply(x$definition$xrefs,
                                function(r) {
                                  if(r$dbCode == "PMID") {
                                    pmid <- r$dbId
                                  } else {
                                    pmid <- NULL
                                  }})
            # Convert the list of dbId values into a string
            pmid <- unlist(pmid_list)
          } else {
            pmid <- ""
          }

          # Collect all info
          all_info <- list(
            "id" = x$id,
            "term_name" = x$name,
            "term_definition" = x$definition$text,
            # "annotated_entrez" = go2entrez,  # 新增：返回 Entrez Gene ID
            # "annotated_symbol" = go2symbol,  # 新增：返回 Gene Symbol
            # "annotated_genename" = go2genename,
            "PMID" = pmid
          )
        }
      )

    go_info <- c(go_info, sub_go_info)
  }
  # # 删除之前添加的有效 GO ID（如果原本不在输入中）
  # if (need_remove) {
  #   go_info <- go_info[!sapply(go_info, function(x) x$id == valid_go_id)]
  # }
  return(go_info)
}


#' Query the QuickGO API for GO term information
#'
#' This internal function sends a request to the QuickGO API and retrieves
#' information about GO terms.
#'
#' @param go_ids A character vector of GO IDs
#'
#' @return A list containing the query results from the QuickGO API
#'
#' @importFrom httr2 request req_headers req_retry req_perform resp_body_json
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
quickgo_api <- function(go_ids) {
  url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", paste(go_ids, collapse = ","))
  tryCatch(
    expr = {
      # Create a request object
      req <- httr2::request(url)

      resp <- req %>%
        httr2::req_headers("Accept" = "application/json") %>%
        httr2::req_retry(max_tries = 3,
                         max_seconds = 60,
                         # Condition for when to retry
                         after = \(resp) is.null(resp) && resp$status_code != 200
        ) %>%
        httr2::req_perform()

      # Parse json to get a list
      info <- httr2::resp_body_json(resp)$result
    },
    error = function(e) {
      warning("Failed to get info from QuickGO after 3 retries:", e$message)
      NULL
    })
}

#' Extract information from KEGG pathways using KEGG IDs
#'
#' This internal function retrieves pathway name, definition, and PMIDs for given KEGG IDs
#' using the KEGGREST package.
#'
#' @param kegg_ids A character vector of KEGG pathway IDs (e.g., "hsa04010")
#'
#' @return A list of KEGG pathway information, where each element contains:
#'   \item{id}{KEGG pathway ID}
#'   \item{term_name}{Name of the pathway}
#'   \item{term_definition}{Description of the pathway}
#'   \item{PMID}{Associated PubMed IDs}
#'
#' @importFrom purrr map
#' @importFrom KEGGREST keggGet
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_kegg_info <- function(kegg_ids) {
  # Process each KEGG ID individually using purrr::map
  kegg_info <-
    kegg_ids %>%
    purrr::map(function(kegg_id) {
      # Get KEGG data for this ID
      x <- KEGGREST::keggGet(dbentries = kegg_id)[[1]]

      # if ("GENE" %in% names(x)) {
      #   # Get annotated Entrez Gene IDs
      #   entrez_ids <- x$GENE[seq(1, length(x$GENE), 2)]
      #
      #   # Get SYMBOL and GENENAME
      #   gene_symbols <- stringr::str_match(x$GENE[seq(2, length(x$GENE), 2)], "^(.*?);")[,2]
      #   gene_names <- stringr::str_match(x$GENE[seq(2, length(x$GENE), 2)], ";\\s*(.*?)\\s*\\[")[,2]
      # } else {
      #   entrez_ids <- NA_character_
      #   gene_symbols <- NA_character_
      #   gene_names <- NA_character_
      # }

      # Extract PMIDs
      pmid <- if ("REFERENCE" %in% names(x)) {
        pmid_list <- purrr::map(x$REFERENCE, function(r) {
          if (grepl("^PMID", r$REFERENCE)) {
            gsub("[^:]+:", "", r$REFERENCE)
          } else {
            NULL
          }
        })
        unlist(pmid_list)
      } else {
        NA_character_
      }

      # Return compiled information
      list(
        "id" = unname(x$ENTRY),
        "term_name" = x$NAME,
        "term_definition" = paste(x$DESCRIPTION, collapse = " "),
        # "annotated_entrez" = entrez_ids,
        # "annotated_symbol" = gene_symbols,
        # "annotated_genename" = gene_names,
        "PMID" = pmid
      )
    })

  return(kegg_info)
}

#' Extract information from Reactome pathways using Reactome IDs
#'
#' This internal function retrieves pathway name, definition, and PMIDs for given Reactome IDs
#' using the rbioapi package.
#'
#' @param reactome_ids A character vector of Reactome pathway IDs
#'
#' @return A list of Reactome pathway information, where each element contains:
#'   \item{id}{Reactome pathway ID}
#'   \item{term_name}{Name of the pathway}
#'   \item{term_definition}{Description of the pathway}
#'   \item{PMID}{Associated PubMed IDs}
#'
#' @importFrom purrr map
#' @importFrom rbioapi rba_reactome_query
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_reactome_info <- function(reactome_ids) {
  # #### 1️⃣ 获取 Reactome ID 到 Entrez Gene ID 的映射
  # suppressMessages(
  #   reactome2egs <- AnnotationDbi::select(reactome.db,
  #                                         keys = reactome_ids,
  #                                         columns = "ENTREZID",
  #                                         keytype = "PATHID")
  # )
  #
  # #### 2️⃣ 获取 Entrez Gene ID 对应的 Gene Symbol 和 Gene Name
  # if (!is.null(reactome2egs$ENTREZID)) {
  #   suppressMessages(
  #     eg2info <- AnnotationDbi::select(orgdb,
  #                                      keys = reactome2egs$ENTREZID,
  #                                      columns = c("SYMBOL", "GENENAME"),
  #                                      keytype = "ENTREZID")
  #   )
  #   eg2info <- dplyr::distinct(eg2info)
  # } else {
  #   eg2info <- data.frame(ENTREZID = NA_character_, SYMBOL = NA_character_, GENENAME = NA_character_)
  # }

  #### 3️⃣ 分批查询 Reactome API
  chunk_size <- 20
  chunks <- split(reactome_ids, ceiling(seq_along(reactome_ids) / chunk_size))

  reactome_info <- list()

  for (i in 1:length(chunks)) {
    reactome_api_result <- suppressMessages(rbioapi::rba_reactome_query(ids = chunks[[i]]))

    if (length(chunks[[i]]) == 1) {
      reactome_api_result <- list(reactome_api_result)
    }

    sub_reactome_info <- purrr::map(
      reactome_api_result,
      function(x) {
        # # **获取 Entrez Gene ID, Gene Symbol, Gene Name**
        # reactome_gene_info <- reactome2egs %>%
        #   dplyr::filter(PATHID == x$stId) %>%
        #   dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
        #   dplyr::left_join(eg2info, by = "ENTREZID")
        #
        # # **修改：直接返回字符向量**
        # annotated_entrez <- reactome_gene_info$ENTREZID
        # annotated_symbol <- reactome_gene_info$SYMBOL
        # annotated_genename <- reactome_gene_info$GENENAME

        # **修改：PMID 直接返回字符向量**
        if ("literatureReference" %in% names(x)) {
          pmid_list <- lapply(x$literatureReference, function(r) {
            if ("pubMedIdentifier" %in% names(r)) {
              return(as.character(r$pubMedIdentifier))
            } else {
              return(NULL)
            }
          })
          pmid <- unlist(pmid_list)  # 直接返回字符向量
        } else {
          pmid <- NA_character_  # 空值返回字符向量
        }

        # **整理输出**
        all_info <- list(
          "id" = x$stId,
          "term_name" = x$displayName,
          "term_definition" = ifelse(!is.null(x$summation), x$summation[[1]]$text, ""),
          # "annotated_entrez" = annotated_entrez,  # ✅ 直接返回字符向量
          # "annotated_symbol" = annotated_symbol,  # ✅ 直接返回字符向量
          # "annotated_genename" = annotated_genename,  # ✅ 直接返回字符向量
          "PMID" = pmid  # ✅ 直接返回字符向量
        )
        return(all_info)
      }
    )

    reactome_info <- c(reactome_info, sub_reactome_info)
  }

  return(reactome_info)
}

#' Extract metabolic pathway information from KEGG
#'
#' This internal function retrieves metabolic pathway names, descriptions, and PMIDs
#' for given KEGG IDs.
#'
#' @param kegg_ids A character vector of KEGG metabolic pathway IDs
#'
#' @return A list of metabolic pathway information, where each element contains:
#'   \item{id}{KEGG pathway ID}
#'   \item{term_name}{Name of the pathway}
#'   \item{term_definition}{Description of the pathway}
#'   \item{PMID}{Associated PubMed IDs}
#'
#' @importFrom KEGGREST keggGet
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_metkegg_info <- function(kegg_ids) {
  metkegg_info <- list()

  for (i in 1:length(kegg_ids)) {
    kegg_id <- kegg_ids[i]
    ### Get a single KEGG entry
    entry <- KEGGREST::keggGet(dbentries = kegg_id)[[1]]

    ### Get PMID for reference
    if ("REFERENCE" %in% names(entry)) {
      pmid_list <- lapply(entry$REFERENCE, function(r) {
        if (grepl("^PMID", r$REFERENCE)) {
          pmid <- gsub("[^:]+:", "", r$REFERENCE)
        } else {
          pmid <- NULL
        }
      })
      pmid <- unlist(pmid_list)
    } else {
      pmid <- ""
    }

    # Collect base info (always included)
    all_info <- list(
      "id" = unname(entry$ENTRY),
      # "term_name" = sub(" - Homo sapiens \\(human\\)$", "", entry$NAME),
      "term_name" = sub(" - [^(]+\\([^)]+\\)$", "", entry$NAME),
      "term_definition" = paste(entry$DESCRIPTION, collapse = " "),
      "PMID" = pmid
    )

    metkegg_info <- c(list(all_info), metkegg_info)
  }

  return(metkegg_info)
}

#' Extract pathway information from SMPDB (small molecule pathway database)
#'
#' This internal function retrieves pathway names, descriptions and related information
#' for given SMPDB IDs using the metpath package.
#'
#' @param smpdb_ids A character vector of SMPDB pathway IDs
#'
#' @return A list of SMPDB pathway information, where each element contains:
#'   \item{id}{SMPDB pathway ID}
#'   \item{term_name}{Name of the pathway}
#'   \item{term_definition}{Description of the pathway}
#'   \item{PMID}{Associated PubMed IDs (NA for SMPDB)}
#'
#' @importFrom purrr map
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_smpdb_info <- function(smpdb_ids) {
  smpdb_pathway_info <- metpath::hmdb_pathway
  smpdb_info <-
    smpdb_ids %>%
    purrr::map(function(x) {
      indx <- which(smpdb_pathway_info@pathway_id == x)
      term_name <- smpdb_pathway_info@pathway_name[indx]
      term_definition <- smpdb_pathway_info@describtion[indx][[1]]
      annotated_compound_name <- smpdb_pathway_info@compound_list[indx][[1]]$Compound.name
      list(
        "id" = x,
        "term_name" = term_name,
        "term_definition" = term_definition,
        "PMID" = NA_character_
      )
    })

  return(smpdb_info)
}

#' Retrieve titles and abstracts from a list of pathway information
#'
#' This internal function extracts titles and abstracts from PubMed for each PMID
#' associated with the pathway data.
#'
#' @param data A list of pathway information containing PMIDs
#'
#' @return A list with pathway information and formatted text information including
#' titles and abstracts from PubMed
#'
#' @importFrom purrr map map_chr
#' @importFrom stringr str_split
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
get_ttl_abstr <- function(data) {
  data %>%
    purrr::map(
      function(x) {
        # Format text info
        if (nchar(x$PMID) == 0){
          text_info <- sprintf("ID: %s\nName: %s\nDefinition: %s\nAnnotated Gene names: %s", x$id, x$term_name, x$term_definition, x$annotated_genename)

        } else {
          ttls_abstrs <-
            stringr::str_split(x$PMID, ",")[[1]] %>%
            purrr::map_chr(~ retrieve_abstr_ttl(pmid = .x)) %>%
            paste(., collapse = "\n")
          text_info <- sprintf("ID: %s\nName: %s\nDefinition: %s\nAnnotated Gene names: %s\nReference:\n%s", x$id, x$term_name, x$term_definition, x$annotated_genename, ttls_abstrs)
        }

        text <- list(
          "id" = x$id,
          "text_info" = text_info
        )
      })
}

#' Retrieve title and abstract for a PubMed ID
#'
#' This internal function retrieves the title and abstract for a given PubMed ID.
#' It also filters out review articles.
#'
#' @param pmid A character string with the PubMed ID
#' @param retries Number of retries when fetching fails (default: 2)
#' @param pause Time in seconds to pause between retries (default: 1)
#'
#' @return A formatted string containing the title and abstract, or an empty string if retrieval fails
#'
#' @importFrom rentrez entrez_fetch
#' @importFrom rvest read_html html_nodes html_node html_text
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
retrieve_abstr_ttl <- function(pmid, retries = 2, pause = 1) {
  attempt <- 1
  while (attempt <= retries) {
    tryCatch({
      # Fetch html data containing title and abstract
      html_result <-
        rentrez::entrez_fetch(
          db = "pubmed",
          id = pmid,
          rettype = "abstract",
          retmode = "html"
        )

      # Parse html
      parsed_html <- rvest::read_html(html_result)

      # Remove Review articles
      publication_type <- parsed_html %>%
        rvest::html_nodes("publicationtype") %>%
        rvest::html_text()
      if ("Review" %in% publication_type) {
        return("")
      }

      # Get title
      title <- parsed_html %>%
        rvest::html_node("articletitle") %>%
        rvest::html_text(trim = TRUE)

      # Get abstract
      abstract <- parsed_html %>%
        rvest::html_node("abstract") %>%
        rvest::html_text(trim = TRUE)

      if (is.na(abstract)) {
        return("")
      } else {
        return(sprintf("Title: %s\nAbstract: %s", title, abstract))
      }
    }, error = function(e) {
      warning(sprintf("Attempt %d failed for PubMed IDs: %s. Error: %s",
                      attempt, paste(pmids, collapse = ", "), e))
    })

    Sys.sleep(pause)
    attempt <- attempt + 1
  }

  warning(sprintf("All attempts failed for PubMed IDs: %s",
                  paste(pmids, collapse = ", ")))
  return(lapply(pmids, function(PID) {
    list(
      title = sprintf("Error: Failed to retrieve title for PubMed ID %s after %d attempts",
                      PID, retries),
      abstract = sprintf("Error: Failed to retrieve abstract for PubMed ID %s after %d attempts",
                         PID, retries)
    )
  }))
}









