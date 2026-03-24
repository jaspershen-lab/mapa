#' Search PubMed for All Modules and Add PubMed IDs
#'
#' This internal function performs a PubMed search for each module in the processed data,
#' retrieves PubMed IDs based on gene symbols, descriptions, and pathways, and appends
#' the IDs to the processed data.
#'
#' @param processed_data A named list where each element corresponds to a module.
#' @param phenotype Character or NULL. Phenotype/disease to focus the search on.
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
#' @noRd
pubmed_search <- function(processed_data, phenotype = NULL, chunk_size = 5, years = 5, retmax = 10, thread = 10) {
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(thread)  # Creates clusters based on available cores
    parallel::clusterExport(cl, varlist = c("process_module", "safe_entrez_search", "perform_query", "build_anchor_block", "test_siliconflow_url"))
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
      result <- withCallingHandlers(
        process_module(module_name, module, phenotype, chunk_size, years, retmax),
        warning = function(w) invokeRestart("muffleWarning")
      )
      return(result)
    }, mc.cores = thread)

    # Unwrap any warning-wrapped results from pbmclapply
    results <- lapply(results, function(r) {
      if (is.list(r) && !is.null(r$value)) r$value else r
    })
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
#' @param phenotype Character or NULL. Phenotype/disease to focus the search on.
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
#' @noRd
process_module <- function(module_name, module, phenotype = NULL, chunk_size = 5, years = 5, retmax = 10) {

  anchor_block <- build_anchor_block(phenotype)

  clean_vec <- function(x) {
    if (length(x) == 0 || all(is.na(x))) return(character(0))
    x <- as.character(x)
    x <- x[!is.na(x)]
    x <- trimws(x)
    unique(x[nzchar(x)])
  }

  if (length(module) == 7) { # for multi_omics module
    entity_terms <- unique(c(
      clean_vec(module$GeneIDs),
      clean_vec(module$GeneNames_vec),
      clean_vec(module$MetNames_vec),
      clean_vec(module$PathwayNames)
    ))
  } else if (length(module) == 6) { # For gene
    entity_terms <- unique(c(
      clean_vec(module$GeneSymbols),
      clean_vec(module$GeneNames_vec),
      clean_vec(module$PathwayNames)
    ))
  } else if (length(module) == 5) { # For metabolites
    entity_terms <- unique(c(
      clean_vec(module$MetNames_vec),
      clean_vec(module$PathwayNames)
    ))
  }

  pmids <- perform_query(entity_terms = entity_terms, anchor_block = anchor_block,
                         years = years, retmax = retmax, chunk_size = chunk_size)

  return(list(module_name = module_name, PubmedIDs = pmids))
}


#' Build Phenotype Anchor Block for PubMed Query
#'
#' @param phenotype Character or NULL. Phenotype/disease anchor term(s).
#' @param field Character. PubMed field tag. Default \code{"tiab"}.
#' @param add_mesh Logical. Add MeSH term for phenotype. Default \code{TRUE}.
#' @return A character string for the anchor block, or \code{NULL}.
#' @noRd
build_anchor_block <- function(phenotype, field = "tiab", add_mesh = TRUE) {
  if (is.null(phenotype)) return(NULL)
  phenotype <- as.character(phenotype)
  phenotype <- phenotype[!is.na(phenotype)]
  phenotype <- trimws(phenotype)
  phenotype <- gsub('"', "", phenotype, fixed = TRUE)
  phenotype <- unique(phenotype[nzchar(phenotype)])
  if (length(phenotype) == 0) return(NULL)

  fmt_term <- function(term) {
    if (grepl("\\s", term)) sprintf('"%s"[%s]', term, field) else sprintf('%s[%s]', term, field)
  }

  anchor_tiab <- paste(vapply(phenotype, fmt_term, character(1)), collapse = " OR ")

  if (add_mesh) {
    phen_mesh <- paste(sprintf('"%s"[Mesh]', phenotype), collapse = " OR ")
    paste0("(", phen_mesh, " OR ", anchor_tiab, ")")
  } else {
    paste0("(", anchor_tiab, ")")
  }
}

#' Perform PubMed Query with Automatic Chunking on Failure
#'
#' This internal function executes a PubMed search. If the full query fails (e.g. HTTP 414),
#' it falls back to chunking the entity terms and retrying each chunk against the fixed anchor
#' block. A third level retries individual terms if a chunk also fails.
#'
#' @param entity_terms Character vector of entity terms (genes, metabolites, pathways) to OR together.
#' @param anchor_block Character string for the pre-built phenotype anchor block, or \code{NULL}.
#' @param years An integer specifying how many years to look back in the search.
#' @param retmax An integer specifying the maximum number of results to retrieve.
#' @param chunk_size An integer specifying the size of entity term chunks.
#' @param field Character. PubMed field tag used when formatting terms. Default \code{"tiab"}.
#'
#' @return A character vector of unique PubMed IDs retrieved from the search.
#'
#' @importFrom rentrez entrez_search
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @noRd
perform_query <- function(entity_terms,
                          anchor_block = NULL,
                          years,
                          retmax,
                          chunk_size,
                          field = "tiab") {

  fmt_term <- function(term) {
    term <- gsub('"', "", term, fixed = TRUE)
    if (grepl("\\s", term)) sprintf('"%s"[%s]', term, field) else sprintf('%s[%s]', term, field)
  }

  build_entity_block <- function(terms) {
    terms <- trimws(terms)
    terms <- terms[!is.na(terms) & nzchar(terms)]
    if (length(terms) == 0) return(NULL)
    paste0("(", paste(vapply(terms, fmt_term, character(1)), collapse = " OR "), ")")
  }

  assemble_query <- function(terms) {
    eb <- build_entity_block(terms)
    if (is.null(eb) && is.null(anchor_block)) return(NULL)
    if (is.null(eb)) return(anchor_block)
    if (is.null(anchor_block)) return(eb)
    paste(anchor_block, "AND", eb)
  }

  entity_terms <- entity_terms[!is.na(entity_terms) & nzchar(trimws(entity_terms))]
  search_ids <- character(0)

  # Level 1: try full query with all entity terms
  full_query <- assemble_query(entity_terms)
  if (!is.null(full_query)) {
    result <- safe_entrez_search(db = "pubmed", term = full_query, retmax = retmax, years = years)
    if (!is.null(result)) {
      return(unique(result$ids))
    }
  }

  # Level 2: chunk entity_terms and retry each chunk against the fixed anchor_block
  if (length(entity_terms) == 0) return(unique(search_ids))
  message("Full query failed (possibly HTTP 414). Splitting entity terms into chunks of ", chunk_size, ".")

  for (i in seq(1, length(entity_terms), by = chunk_size)) {
    chunk <- entity_terms[i:min(i + chunk_size - 1, length(entity_terms))]
    chunk_query <- assemble_query(chunk)
    if (is.null(chunk_query)) next

    chunk_result <- safe_entrez_search(db = "pubmed", term = chunk_query, retmax = retmax, years = years)
    if (!is.null(chunk_result)) {
      search_ids <- c(search_ids, chunk_result$ids)
    } else {
      # Level 3: individual terms within this chunk
      message("Chunk query failed. Attempting individual terms.")
      for (term in chunk) {
        term_query <- assemble_query(term)
        if (is.null(term_query)) next
        term_result <- safe_entrez_search(db = "pubmed", term = term_query, retmax = retmax, years = years)
        if (!is.null(term_result)) {
          search_ids <- c(search_ids, term_result$ids)
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
#' @noRd
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

#' Build a PubMed Query String
#'
#' @param pathway_names Character vector of pathway names.
#' @param gene_symbols Character vector of gene symbols, or \code{NA}.
#' @param gene_names Character vector of gene full names, or \code{NA}.
#' @param met_names Character vector of metabolite names, or \code{NA}.
#' @param phenotype Character or NULL. Phenotype/disease anchor term.
#' @param field Character. PubMed field tag. Default \code{"tiab"}.
#' @param add_mesh_for_anchor Logical. Add MeSH term for phenotype anchor. Default \code{TRUE}.
#' @return A character string containing the assembled PubMed query.
#' @noRd
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
