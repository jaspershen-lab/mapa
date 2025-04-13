# source("R/17_llm_interpretation/modules/utils.R")
# library(pdftools)
# library(stringr)
# library(httr)
# library(jsonlite)
# library(dplyr)
# library(data.table)
# library(rentrez)
# library(parallel)
# library(pbmcapply)


#' Process Local PDF Corpus for Embedding
#'
#' Processes all PDF files in a specified local corpus directory by generating embeddings and saving them to the database.
#'
#' @param api_key A character string containing the API key for the embedding service.
#' @param local_corpus_dir A character string specifying the name of the directory where the local files
#'        provided by users are saved. Defaults to "local_corpus".
#' @param embedding_output_dir A character string specifying the directory where the intermediate embedding results
#'        will be saved. Defaults to "embedding_output".
#' @param save_dir_local_corpus_embed A character string specifying the directory where the embedding data will be saved.
#'        Defaults to "local".
#'
#' @return Nothing. The function processes and saves the embeddings to the specified save directory.
#'
#' @importFrom tools file_ext
#'
#' @examples
#' \dontrun{
#' # Example: Process all local PDFs in the "local_corpus" folder
#' api_key <- "your_api_key"
#' embedding_local_corpus(api_key)
#'
#' # Example: Process PDFs in a custom directory and save to a custom location
#' embedding_local_corpus(api_key, local_corpus_dir = "my_pdfs",
#'                        embedding_output_dir = "embed_output",
#'                        save_dir_local_corpus_embed = "custom_save")
#' }
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Checks the specified folder in the current working directory for PDF files.
#'   \item Validates that all files in the folder are in PDF format.
#'   \item Processes each PDF file using \code{\link{embedding_single_pdf}} to generate embeddings.
#'   \item Saves the generated embeddings to the specified directory using \code{\link{save_embedding}}.
#' }
#' If non-PDF files are found or if the folder is empty, the function throws an error.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
embedding_local_corpus <-
  function(api_key = NULL,
           local_corpus_dir = "local_corpus",
           embedding_output_dir = "embedding_output",
           save_dir_local_corpus_embed = "local") {
  if (is.null(api_key)) {
    stop("Please provide api key to do embedding.")
  }

  current_dir <- getwd()
  ## get file path where user store their uploaded pdf file.
  local_corpus_dir <- file.path(current_dir, local_corpus_dir)
  files <- list.files(local_corpus_dir, full.names = TRUE)

  if (length(files) == 0) {
    stop(paste("Error: No files found in", local_corpus_dir, "- please check the path is correct."))
  }

  pdf_files <- files[tolower(tools::file_ext(files)) == "pdf"]

  if (length(pdf_files) != length(files)) {
    stop("Error: Non-PDF files found in the 'local_corpus' folder. Please ensure all files are in PDF format.")
  }

  for (pdf_path in pdf_files) {
    print(paste("Embedding in progress for file:", pdf_path))
    pdf_embeddings <- embedding_single_pdf(pdf_path, api_key)
    save_embedding(pdf_embeddings, embedding_output_dir = embedding_output_dir, save_dir = save_dir_local_corpus_embed)
  }
}

#' Process Single PDF for Embedding
#'
#' Processes a single PDF file by extracting text, splitting it into manageable chunks, and generating embeddings for each chunk.
#'
#' @param pdf_path A character string specifying the path to the PDF file.
#' @param api_key A character string containing the API key for the embedding service.
#'
#' @return A data frame containing the following columns:
#' \item{paper_title}{The title of the PDF (derived from the file name).}
#' \item{chunks}{The extracted text chunks from the PDF.}
#' \item{embedding_vector}{A list of embedding vectors for each chunk.}
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom pdftools pdf_text
#' @importFrom parallel makeCluster clusterExport clusterEvalQ parLapply stopCluster detectCores
#' @importFrom pbmcapply pbmclapply
#'
#' @examples
#' \dontrun{
#' # Example: Process a single PDF to generate embeddings
#' pdf_path <- "example.pdf"
#' api_key <- "your_api_key"
#' embedding_data <- embedding_single_pdf(pdf_path, api_key)
#' print(embedding_data)
#' }
#' @details
#' This function performs the following:
#' \enumerate{
#'   \item Extracts text from the PDF using \code{\link{pdf_text}}.
#'   \item Splits the text into smaller chunks using \code{\link{split_into_chunks}}.
#'   \item Generates embeddings for each chunk using a parallelized approach (\code{\link{parallel::parLapply}} or \code{\link{pbmclapply}}).
#' }
#' The resulting data frame contains the title of the PDF, the processed chunks, and their corresponding embedding vectors.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
embedding_single_pdf <- function(pdf_path, api_key){
  paper_title <- tools::file_path_sans_ext(basename(pdf_path))
  text_pages <- pdftools::pdf_text(pdf_path)
  full_text <- paste(text_pages, collapse = "\n")
  chunks <- split_into_chunks(full_text)

  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(min(detectCores()-1, 10)) # Creates four clusters
    parallel::clusterExport(cl, varlist = c("chunks", "api_key", "get_embedding"), envir = environment()) # Export the variables into the global environment of newly created clusters so that they can use them
    parallel::clusterEvalQ(cl, {
      library(httr2)
    })
    embeddings <- parallel::parLapply(cl, chunks,
                                      function(chunk) {
                                        Sys.sleep(1)
                                        embedding <- get_embedding(chunk, api_key)
                                        return(embedding)
                                      })
    parallel::stopCluster(cl) # Stop the clusters and close the parallel backend.
  } else if (.Platform$OS.type == "unix") {
    embeddings <- pbmclapply(chunks, function(chunk) {
      Sys.sleep(1)
      embedding <- get_embedding(chunk, api_key)
      return(embedding)
    }, mc.cores = detectCores()-1)
  }
  chunks_matrix <- do.call(rbind, chunks)

  df <- data.frame(
    paper_title = rep(paper_title, nrow(chunks_matrix)),
    chunks = I(chunks_matrix),
    stringsAsFactors = FALSE
  )

  df$embedding_vector <- embeddings
  return(df)
}

#' Split Full Text into Manageable Chunks
#'
#' Splits a given full text into smaller chunks based on specified parameters such as maximum chunk length and minimum word count.
#'
#' @param full_text A character string representing the full text to be split into chunks.
#' @param max_chunk_length An integer specifying the maximum number of words allowed in a single chunk (default is 800).
#' @param min_words An integer specifying the minimum number of words a line must contain to be included in a chunk (default is 5).
#' @param space_threshold An integer specifying the minimum number of spaces in a line to consider it as a double-column text (default is 5).
#' @param min_chunk_length An integer specifying the minimum number of words required in a chunk to be included in the result (default is 50).
#'
#' @return A character vector where each element is a text chunk.
#'
#' @importFrom stringr str_split str_count str_trim str_c str_detect
#'
#' @examples
#' \dontrun{
#' # Example: Split a full text into chunks
#' full_text <- "This is the first line.\nThis is the second line, which has more words.\n\nReferences\nThis line will be excluded."
#' chunks <- split_into_chunks(full_text)
#' print(chunks)
#' }
#' @details
#' This function handles both single-column and double-column texts. Double-column texts are identified based on the number of spaces in a line (\code{space_threshold}). The function also filters out chunks containing references or URLs and ensures that only chunks with at least \code{min_chunk_length} words are included.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
split_into_chunks <- function(full_text, max_chunk_length = 800, min_words = 5, space_threshold = 5, min_chunk_length = 50) {
  lines <- unlist(str_split(full_text, "\n"))
  lines <- lines[str_count(lines, "\\S+") >= min_words]

  chunks <- list()
  current_chunk <- ""
  in_double_column <- FALSE

  for (i in seq_along(lines)) {
    line <- lines[i]
    space_count <- str_count(line, " ")
    if (space_count >= space_threshold) {
      columns <- unlist(str_split(line, "\\s{2,}"))
      if (!in_double_column) {
        if (str_count(current_chunk, "\\S+") > 0) {
          chunks <- append(chunks, current_chunk)
        }
        current_chunk <- ""
        in_double_column <- TRUE
      }

      for (col in columns) {
        col <- str_trim(col)
        if (str_count(col, "\\S+") >= min_words) {
          if (str_count(str_trim(current_chunk), "\\S+") + str_count(str_trim(col), "\\S+") > max_chunk_length) {
            chunks <- append(chunks, current_chunk)
            current_chunk <- col
          } else {
            current_chunk <- str_c(current_chunk, col, sep = " ")
          }
        }
      }
    } else {
      if (in_double_column) {
        chunks <- append(chunks, current_chunk)
        current_chunk <- ""
        in_double_column <- FALSE
      }

      if (str_count(str_trim(current_chunk), "\\S+") + str_count(str_trim(line), "\\S+") > max_chunk_length) {
        chunks <- append(chunks, current_chunk)
        current_chunk <- line
      } else {
        current_chunk <- str_c(current_chunk, line, sep = " ")
      }
    }
  }

  if (str_trim(current_chunk) != "") {
    chunks <- append(chunks, current_chunk)
  }

  chunks <- chunks[!str_detect(chunks, regex("^References", ignore_case = TRUE))]
  chunks <- chunks[!str_detect(chunks, regex("http[s]?://|www\\.", ignore_case = TRUE))]
  chunks <- chunks[sapply(chunks, function(chunk) str_count(str_trim(chunk), "\\S+") >= min_chunk_length)]

  return(chunks)
}

#' Save Embedding Data
#'
#' Saves embedding data (including titles, chunks, and embedding vectors) to a specified directory.
#'
#' @param embedding_df A data frame with three columns: \code{paper_title}, \code{chunk}, and \code{embedding_vector}.
#'        \code{embedding_vector} should be a list of numeric vectors.
#' @param embedding_output_dir A character string specifying the path to the directory where all the intermediate embedding results to be saved.
#'        Defaults is "embedding_output".
#' @param save_dir A character string specifying the directory where the embedding data will be saved.
#'        Defaults to "local".
#'
#' @return Nothing. The function saves the data to files in the specified directory.
#'
#' @importFrom data.table as.data.table fwrite
#'
#' @examples
#' \dontrun{
#' # Example: Save embeddings to a directory
#' embedding_df <- data.frame(
#'   paper_title = c("Title 1", "Title 2"),
#'   chunks = c("Chunk 1", "Chunk 2"),
#'   embedding_vector = I(list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)))
#' )
#' save_embedding(embedding_df, embedding_output_dir = "embed_output", save_dir = "example_dir")
#' }
#' @details
#' This function performs the following:
#' \itemize{
#'   \item Cleans the \code{paper_title} and \code{chunks} by removing extra spaces, line breaks, and tabs.
#'   \item Writes the titles and chunks to \code{paper_title.txt} and \code{chunk.txt}, respectively.
#'   \item Saves the embedding vectors to \code{embedding_vector.csv.gz} in compressed format.
#'   \item Appends to existing files if they already exist in the directory.
#' }
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
save_embedding <- function(embedding_df, embedding_output_dir = "embedding_output", save_dir = "local") {
  if (ncol(embedding_df) != 3) {
    stop("Data must contain three columns: paper_title, chunk, embedding_vector")
  }

  output_dir <- file.path(embedding_output_dir, save_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  embedding_df$paper_title <- gsub("\\s+", " ", embedding_df$paper_title)
  embedding_df$chunks <- gsub("\\s+", " ", embedding_df$chunks)
  embedding_df$paper_title <- gsub("[\r\n\t]", "", embedding_df$paper_title)
  embedding_df$chunks <- gsub("[\r\n\t]", "", embedding_df$chunks)

  paper_title_file <- file.path(output_dir, "paper_title.txt")
  chunk_file <- file.path(output_dir, "chunk.txt")
  embedding_vector_file <- file.path(output_dir, "embedding_vector.csv.gz")

  write(embedding_df$paper_title, file = paper_title_file, append = TRUE, sep = "\n")
  write(embedding_df$chunks, file = chunk_file, append = TRUE, sep = "\n")

  embedding_matrix <- do.call(rbind, embedding_df$embedding_vector)
  dt_embeddings <- as.data.table(embedding_matrix)

  if (!file.exists(embedding_vector_file)) {
    fwrite(dt_embeddings,
           file = embedding_vector_file,
           compress = "gzip",
           nThread = 4)
  } else {
    fwrite(dt_embeddings,
           file = embedding_vector_file,
           append = TRUE,
           compress = "gzip",
           nThread = 4)
  }
}


#' Generate Embeddings for PubMed Search Results Across Multiple Modules
#'
#' @description
#' Processes PubMed search results for multiple modules by generating text embeddings
#' for abstracts and titles, then saving these embeddings to specified directories.
#'
#' @param pubmed_result A named list where each element corresponds to a module.
#'   Each module contains a vector of PubMed IDs under the `PubmedIDs` field.
#' @param api_key Character string containing the API key for the embedding service.
#' @param embedding_output_dir Character string specifying the base directory where
#'   embedding results will be saved.
#'
#' @return No return value; function saves embeddings to disk in the specified output directory.
#'
#' @details
#' For each module in `pubmed_result`, this function:
#' \itemize{
#'   \item Extracts PubMed IDs for the module
#'   \item Processes these IDs by calling `embedding_single_module_pubmed_search()`
#'   \item Displays progress information to the console
#' }
#'
#' @examples
#' \dontrun{
#' # Process PubMed search results for multiple modules
#' pubmed_result <- list(
#'   immune_response = list(PubmedIDs = c("12345678", "23456789")),
#'   cancer_genetics = list(PubmedIDs = c("34567890", "45678901"))
#' )
#' api_key <- "your_openai_api_key"
#' embedding_output_dir <- "path/to/embeddings"
#' embedding_pubmed_search(pubmed_result, api_key, embedding_output_dir)
#' }
#'
#' @seealso \code{\link{embedding_single_module_pubmed_search}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal

embedding_pubmed_search <- function(pubmed_result, api_key, embedding_output_dir){
  module_names <- names(pubmed_result)

  for (module_name in module_names) {
    module <- pubmed_result[[module_name]]
    PID_list <- module$PubmedIDs
    cat(sprintf("Processing module: %s\n", module_name))
    cat(sprintf("Including PID number: %s\n", length(PID_list)))
    embedding_single_module_pubmed_search(module_name, PID_list,api_key, embedding_output_dir)
  }
}

#' Generate Embeddings for PubMed Abstracts and Titles for a Single Module
#'
#' @description
#' Retrieves abstracts and titles for a set of PubMed IDs, generates text embeddings
#' for each abstract, and saves the resulting data to a specified directory.
#'
#' @param module_name Character string specifying the name of the module (used for
#'   organizing output files).
#' @param PID_list Character vector containing PubMed IDs to process.
#' @param api_key Character string containing the API key for the embedding service.
#' @param embedding_output_dir Character string specifying the base directory where
#'   embedding results will be saved.
#'
#' @return No return value; function saves embeddings to disk in module-specific directory.
#'
#' @details
#' This function performs these steps:
#' \itemize{
#'   \item Splits PubMed IDs into batches of 5 to avoid API rate limits
#'   \item Retrieves abstracts and titles for each batch using `query_abstracts_and_titles_by_pubmed_ids()`
#'   \item Generates embeddings for each abstract using parallel processing
#'   \item Creates a data frame with paper titles, abstracts, and embedding vectors
#'   \item Saves the results using `save_embedding()` in the specified directory
#' }
#'
#' The function uses different parallel processing methods based on the operating system:
#' \itemize{
#'   \item On Windows: Uses `parallel::parLapply()` with explicit cluster management
#'   \item On other platforms: Uses `pbmclapply()` for progress-tracked parallel processing
#' }
#'
#' @examples
#' \dontrun{
#' # Process PubMed IDs for a single module
#' module_name <- "immune_response"
#' PID_list <- c("12345678", "23456789", "34567890")
#' api_key <- "your_openai_api_key"
#' embedding_output_dir <- "path/to/embeddings"
#' embedding_single_module_pubmed_search(module_name, PID_list, api_key, embedding_output_dir)
#' }
#'
#' @importFrom parallel makeCluster clusterExport clusterEvalQ parLapply stopCluster detectCores
#' @importFrom httr2 req_perform req_body_json
#'
#' @seealso
#' \code{\link{query_abstracts_and_titles_by_pubmed_ids}},
#' \code{\link{get_embedding}},
#' \code{\link{save_embedding}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
embedding_single_module_pubmed_search <- function(module_name, PID_list, api_key, embedding_output_dir) {
  # 将PID_list分成每组5个的批次
  batch_size <- 5
  PID_batches <- split(PID_list, ceiling(seq_along(PID_list)/batch_size))

  # 按批次获取摘要和标题
  abstracts_and_titles_list <- list()
  for(batch in PID_batches) {
    batch_results <- query_abstracts_and_titles_by_pubmed_ids(batch)
    abstracts_and_titles_list <- c(abstracts_and_titles_list, batch_results)
    Sys.sleep(1)  # 添加短暂延迟以避免API限制
  }

  # 提取所有标题和摘要
  titles <- sapply(abstracts_and_titles_list, function(x) x$title)
  abstracts <- sapply(abstracts_and_titles_list, function(x) x$abstract)

  # 并行处理embeddings生成
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(min(parallel::detectCores()-1, 10))
    parallel::clusterExport(cl, varlist = c("get_embedding", "api_key"))
    parallel::clusterEvalQ(cl, {
      library(httr2)
    })

    embeddings <- parallel::parLapply(cl, abstracts, function(abstract) {
      get_embedding(abstract, api_key)
    })

    parallel::stopCluster(cl)
  } else {
    embeddings <- pbmclapply(abstracts, function(abstract) {
      get_embedding(abstract, api_key)
    }, mc.cores = parallel::detectCores() - 1)
  }

  embedding_df <- data.frame(
    paper_title = titles,
    chunks = abstracts,
    embedding_vector = I(embeddings),
    stringsAsFactors = FALSE
  )

  save_embedding(embedding_df, embedding_output_dir = embedding_output_dir, save_dir = module_name)
}

#' Retrieve Titles and Abstracts by PubMed IDs
#'
#' @description
#' Queries PubMed for titles and abstracts using a list of PubMed IDs (PIDs).
#'
#' @param PID_list A character vector containing PubMed IDs to query.
#' @param retries An integer specifying the number of retry attempts for failed queries (default is 2).
#' @param pause A numeric value specifying the time in seconds to pause between retries (default is 1).
#'
#' @return A named list where each element corresponds to a PubMed ID. Each element is a list containing:
#' \item{title}{The retrieved title for the PubMed ID. If retrieval fails, an error message is included.}
#' \item{abstract}{The retrieved abstract for the PubMed ID. If retrieval fails, an error message is included.}
#'
#' @details
#' This function performs the following steps for each PubMed ID:
#' \itemize{
#'   \item Queries PubMed for the title and abstract using \code{\link{retrieve_abstr_ttl}}.
#'   \item Handles errors and retries the query up to \code{retries} times in case of failure.
#'   \item Returns a default error message if all retry attempts fail.
#' }
#'
#' The function ensures robust error handling and provides meaningful messages in case of retrieval failures.
#'
#' @examples
#' \dontrun{
#' # Example: Query titles and abstracts for a list of PubMed IDs
#' PID_list <- c("12345678", "23456789")
#' results <- query_abstracts_and_titles_by_pubmed_ids(PID_list)
#' print(results)
#' }
#'
#' @seealso \code{\link{retrieve_abstr_ttl}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
query_abstracts_and_titles_by_pubmed_ids <- function(PID_list, retries = 2, pause = 1) {
  if (is.null(PID_list) || length(PID_list) == 0) {
    stop("Invalid PubMed ID list provided.")
  }

  # 将PID_list逐个处理
  results <- list()
  for (pmid in PID_list) {
    attempt <- 1
    while (attempt <= retries) {
      tryCatch({
        # 调用 retrieve_abstr_ttl 函数获取标题和摘要
        abstr_ttl <- retrieve_abstr_ttl(pmid, retries, pause)

        if (abstr_ttl == "") {
          results[[pmid]] <- list(
            title = sprintf("Error: No abstract found for PubMed ID %s", pmid),
            abstract = sprintf("Error: No abstract found for PubMed ID %s", pmid)
          )
        } else {
          # 拆分标题和摘要
          title <- sub("^Title: (.+?)\\nAbstract:.*$", "\\1", abstr_ttl)
          abstract <- sub("^Title: .+?\\nAbstract: (.+)$", "\\1", abstr_ttl)

          results[[pmid]] <- list(
            title = title,
            abstract = abstract
          )
        }
        break
      }, error = function(e) {
        warning(sprintf("Attempt %d failed for PubMed ID: %s. Error: %s",
                        attempt, pmid, e))
        Sys.sleep(pause)
      })
      attempt <- attempt + 1
    }

    if (attempt > retries) {
      results[[pmid]] <- list(
        title = sprintf("Error: Failed to retrieve title for PubMed ID %s after %d attempts", pmid, retries),
        abstract = sprintf("Error: Failed to retrieve abstract for PubMed ID %s after %d attempts", pmid, retries)
      )
    }
  }

  return(results)
}

#' Get Dataset Dimensions
#'
#' @description
#' Retrieves the dimensions of a specified dataset based on provided dataset information.
#'
#' @param info A data frame containing dataset metadata. It should include columns such as \code{name},
#'   \code{group}, and \code{dim}.
#' @param dataset_path A character string specifying the dataset's path in the format "group/dataset_name".
#'
#' @return A numeric vector representing the dimensions of the dataset.
#'
#' @details
#' The function parses the dataset path to extract group and name components, then searches for a matching
#' dataset in the provided metadata. The dimensions are extracted from the 'dim' column, which should be
#' formatted as "rowsxcolumns" (e.g., "10x20").
#'
#' @examples
#' \dontrun{
#' # Example: Retrieve dimensions for a specific dataset
#' dataset_info <- data.frame(
#'   name = c("dataset1", "dataset2"),
#'   group = c("/group1", "/group2"),
#'   dim = c("10x20", "30x40")
#' )
#' dims <- get_dataset_dim(dataset_info, "group1/dataset1")
#' print(dims) # Output: c(10, 20)
#' }
#'
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
get_dataset_dim <- function(info, dataset_path) {
  path_parts <- strsplit(dataset_path, "/")[[1]]
  group <- paste0("/", paste(path_parts[-length(path_parts)], collapse = "/"))
  name <- path_parts[length(path_parts)]
  dataset_info <- info[info$name == name & info$group == group, ]

  if (nrow(dataset_info) == 0) {
    stop(paste("Dataset", dataset_path, "not found"))
  }

  dims <- as.numeric(strsplit(dataset_info$dim, "x")[[1]])
  return(dims)
}
