#' Search PubMed for All Modules and Add PubMed IDs
#'
#' This internal function performs a PubMed search for each module in the processed data,
#' retrieves PubMed IDs based on gene symbols, descriptions, and pathways, and appends
#' the IDs to the processed data.
#'
#' @param processed_data A named list where each element corresponds to a module.
#' @param chunk_size An integer specifying the size of query chunks (default is 5).
#' @param years An integer specifying how many years to look back in the search (default is 5).
#' @param retmax An integer specifying the maximum number of results to retrieve (default is 10).
#' @param thread An integer specifying the number of parallel threads to use for processing.
#'   Default is `10` for sequential processing.
#'
#' @return A named list similar to \code{processed_data}, but with an added \code{PubmedIDs}
#'   field for each module, containing the retrieved PubMed IDs.
#'
#' @importFrom parallel makeCluster clusterExport clusterEvalQ parLapply stopCluster detectCores
#' @importFrom pbmcapply pbmclapply
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
pubmed_search <- function(processed_data, phenotype, chunk_size = 5, years = 5, retmax = 10, thread = 10) {
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(thread)  # Creates clusters based on available cores
    parallel::clusterExport(cl, varlist = c("process_module", "safe_entrez_search", "perform_query","test_siliconflow_url"))
    parallel::clusterExport(cl, varlist = c("phenotype", "chunk_size", "years", "retmax", "thread"), envir = environment())
    parallel::clusterEvalQ(cl, {
      library(rentrez)
      library(curl)
    })

    results <- parallel::parLapply(cl, names(processed_data), function(module_name) {
      module <- processed_data[[module_name]]
      result <- process_module(module_name, module, phenotype, chunk_size, years, retmax)
      return(result)
    })

    parallel::stopCluster(cl)
  } else {
    results <- pbmcapply::pbmclapply(names(processed_data), function(module_name) {
      module <- processed_data[[module_name]]
      result <- process_module(module_name, module, phenotype, chunk_size, years, retmax)
      return(result)
    }, mc.cores = thread)
  }

  for (result in results) {
    processed_data[[result$module_name]]$PubmedIDs <- result$PubmedIDs
  }

  return(processed_data)
}


#' Process a Single Module to Retrieve PubMed IDs
#'
#' This internal function processes a module to query PubMed IDs based on either gene information
#' or metabolite information, combined with pathway data. It detects the module type automatically
#' based on the module's structure.
#'
#' @param module_name A character string specifying the name of the module.
#' @param module A list containing module information:
#'   \itemize{
#'     \item For gene modules (length 6): Contains PathwayNames, GeneSymbols, and GeneNames_vec
#'     \item For metabolite modules (length 5): Contains PathwayNames and MetNames_vec
#'   }
#' @param chunk_size An integer specifying the size of query chunks (default is 5).
#' @param years An integer specifying how many years to look back in the search (default is 5).
#' @param retmax An integer specifying the maximum number of results to retrieve (default is 10).
#'
#' @return A list containing:
#'   \item{module_name}{The name of the module.}
#'   \item{PubmedIDs}{A character vector of unique PubMed IDs retrieved for the module.}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
process_module <- function(module_name, module, phenotype = NULL, chunk_size = 5, years = 5, retmax = 10) {

  if (length(module) == 7) { # for multi_omics module
    query <- build_pubmed_query(pathway_names = module$PathwayNames,
                                gene_symbols = module$GeneIDs,
                                gene_names = module$GeneNames_vec,
                                met_names = module$MetNames_vec,
                                phenotype = phenotype)

    pmids <- perform_query(query = query, years = years, retmax = retmax, chunk_size = chunk_size)
  } else if (length(module) == 6) { # For gene
    pathway_names <- module$PathwayNames
    gene_symbols <- module$GeneSymbols
    gene_names <- module$GeneNames_vec

    gene_symbol_query <- build_pubmed_query(pathway_names = module$PathwayNames,
                                            gene_symbols = module$GeneSymbols,
                                            gene_names = NA,
                                            met_names = NA,
                                            phenotype = phenotype)
    gene_symbol_ids <- perform_query(query = query, years = years, retmax = retmax, chunk_size = chunk_size)

    gene_name_query <- build_pubmed_query(pathway_names = module$PathwayNames,
                                          gene_symbols = NA,
                                          gene_names = module$GeneNames_vec,
                                          met_names = NA,
                                          phenotype = phenotype)
    gene_name_ids <- perform_query(query = query, years = years, retmax = retmax, chunk_size = chunk_size)
    pmids <- unique(c(gene_symbol_ids, gene_name_ids))

  } else if (length(module) == 5) { # For metabolites
    query <- build_pubmed_query(pathway_names = module$PathwayNames,
                                gene_symbols = NA,
                                gene_names = NA,
                                met_names = module$MetNames_vec,
                                phenotype = phenotype)

    pmids <- perform_query(query = query, years = years, retmax = retmax, chunk_size = chunk_size)
  }

  return(list(module_name = module_name, PubmedIDs = pmids))
}

#' Perform PubMed Query with Terms and Pathway Information
#'
#' This internal function executes a PubMed search using query terms and pathway information,
#' with fallback strategies for handling failed queries. It implements a three-level search
#' strategy: first attempting a full query with all terms, then breaking into chunks if that
#' fails, and finally trying individual terms if chunk queries also fail.
#'
#' @param query_terms A character vector of query terms (e.g., gene symbols, gene names, metabolite names).
#' @param pathway_query A character string representing the pathway part of the query.
#' @param years An integer specifying how many years to look back in the search.
#' @param retmax An integer specifying the maximum number of results to retrieve.
#' @param chunk_size An integer specifying the size of query chunks.
#'
#' @return A character vector of unique PubMed IDs retrieved from the search.
#'
#' @importFrom rentrez entrez_search
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @keywords internal
perform_query <- function(query,
                          years,
                          retmax,
                          chunk_size) {
  search_ids <- c()

  result <- safe_entrez_search(db = "pubmed", term = query, retmax = retmax, years = years)

  if (!is.null(result)) {
    search_ids <- c(search_ids, result$ids)
  } else {
    message(sprintf("Full query failed. Attempting to split into chunks."))
    for (i in seq(1, length(query_terms), by = chunk_size)) {
      term_chunk <- query_terms[i:min(i + chunk_size - 1, length(query_terms))]
      term_query <- paste(term_chunk, collapse = " OR ")
      chunk_query <- paste("(", term_query, ")", "AND", "(", pathway_query, ")", sep = " ")

      chunk_result <- safe_entrez_search(db = "pubmed", term = chunk_query, retmax = retmax, years = years)
      if (!is.null(chunk_result)) {
        search_ids <- c(search_ids, chunk_result$ids)
      } else {
        message(sprintf("Chunk query failed. Attempting to split into individual terms."))
        for (sub_term in term_chunk) {
          sub_query <- paste("(", sub_term, ")", "AND", "(", pathway_query, ")", sep = " ")
          sub_result <- safe_entrez_search(db = "pubmed", term = sub_query, retmax = retmax, years = years)
          if (!is.null(sub_result)) {
            search_ids <- c(search_ids, sub_result$ids)
          }
        }
      }
    }
  }

  return(unique(search_ids))
}

#' Safely Perform an Entrez Search with Retries
#'
#' This internal function performs a PubMed search using Entrez with a specified number
#' of retries and pauses between attempts. It also handles date filtering using the years
#' parameter to limit search results to recent publications.
#'
#' @param db A character string specifying the Entrez database to search (e.g., "pubmed").
#' @param term A character string specifying the search term or query to execute.
#' @param retmax An integer specifying the maximum number of results to retrieve (default is 10).
#' @param retries An integer specifying the number of retry attempts if the search fails (default is 3).
#' @param pause A numeric value specifying the time in seconds to pause between retry attempts (default is 5).
#' @param years An integer specifying how many years to look back in the search (default is 5).
#'
#' @return A list containing the search results if successful, or \code{NULL} if all retries fail.
#'
#' @importFrom rentrez entrez_search
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal

safe_entrez_search <- function(db, term, retmax = 10, retries = 3, pause = 5, years = 5) {
  # Validate years if provided
  if (!is.null(years)) {
    if (!is.numeric(years) || years < 0) {
      stop("'years' must be a non-negative numeric value.")
    }
    start_year <- as.integer(format(Sys.Date(), "%Y")) - years
    date_filter <- paste0(start_year, "[PDAT] : ", format(Sys.Date(), "%Y"), "[PDAT]")
    term <- paste0("(", term, ") AND (", date_filter, ")")  # Fix: Proper AND and parentheses
  }

  attempt <- 1
  while (attempt <= retries) {
    result <- tryCatch({
      rentrez::entrez_search(db = db, term = term, retmax = retmax)
    }, error = function(e) {
      warning(sprintf("Attempt %s failed: %s", attempt, e$message))
      NULL
    })
    if (!is.null(result)) {
      return(result)
    }
    Sys.sleep(pause)
    attempt <- attempt + 1
  }
  return(NULL)
}

# Desired structure:
# ( <ANCHOR> ) AND ( (<GENE_OR_LIST>) OR (<METABOLITE_OR_LIST>) OR (<PATHWAY_OR_LIST>) )
build_pubmed_query <- function(
    pathway_names,
    gene_symbols,
    gene_names,
    met_names,
    phenotype = NULL,
    field = "tiab",
    add_mesh_for_anchor = TRUE
) {
  clean_terms <- function(x) {
    if (is.null(x)) return(character(0))
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- trimws(x)
    x <- x[nzchar(x)]
    x <- gsub('"', "", x, fixed = TRUE)   # avoid breaking quotes in query
    unique(x)
  }

  fmt_term <- function(term, field) {
    # Quote phrases with spaces; keep simple.
    if (grepl("\\s", term)) sprintf('"%s"[%s]', term, field) else sprintf('%s[%s]', term, field)
  }

  fmt_or_block <- function(terms, field) {
    terms <- clean_terms(terms)
    if (length(terms) == 0) return(NULL)
    paste0("(", paste(vapply(terms, fmt_term, character(1), field = field), collapse = " OR "), ")")
  }

  # Anchor (phenotype)
  anchor_block <- NULL
  if (!is.null(phenotype)) {
    phenotype <- clean_terms(phenotype)
    if (length(phenotype) > 0) {
      # phenotype is a single string; if user passes vector, we OR them
      anchor_tiab <- paste(vapply(phenotype, fmt_term, character(1), field = field), collapse = " OR ")

      if (add_mesh_for_anchor) {
        # Mesh term is typically quoted
        phen_mesh <- paste(sprintf('"%s"[Mesh]', phenotype), collapse = " OR ")
        anchor_block <- paste0("(", phen_mesh, " OR ", anchor_tiab, ")")
      } else {
        anchor_block <- paste0("(", anchor_tiab, ")")
      }
    }
  }

  # Entity blocks
  gene_terms <- c(clean_terms(gene_symbols), clean_terms(gene_names))
  gene_block <- fmt_or_block(gene_terms, field = field)

  met_block  <- fmt_or_block(met_names, field = field)
  path_block <- fmt_or_block(pathway_names, field = field)

  # Combine entity blocks with OR (skip NULLs)
  entity_blocks <- Filter(Negate(is.null), list(gene_block, met_block, path_block))
  entity_part <- if (length(entity_blocks) == 0) NULL else paste0("(", paste(entity_blocks, collapse = " OR "), ")")

  # Final assembly
  if (is.null(anchor_block) && is.null(entity_part)) {
    stop("No valid terms to build query: phenotype is NULL/empty AND all entity lists are empty/NA.")
  } else if (is.null(anchor_block)) {
    return(entity_part)
  } else if (is.null(entity_part)) {
    return(anchor_block)
  } else {
    return(paste(anchor_block, "AND", entity_part))
  }
}
