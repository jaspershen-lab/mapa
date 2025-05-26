#This file is for RAG strategy that is retrieve related chunks in database
# source("R/17_llm_interpretation/modules/utils.R")
# library(data.table)
# library(Matrix)
# library(parallel)
# library(pbmcapply)
# library(roxygen2)


#' Read Embeddings from Storage
#'
#' @description
#' Reads embedding vectors from a specified directory. This function supports batch reading
#' of data, allowing partial reads from large compressed files to optimize memory usage.
#'
#' @param save_dir A character string specifying the output directory where the embedding file is stored.
#' @param start_row An integer indicating the starting row to read from (default is 1).
#' @param n_rows An integer specifying the number of rows to read. Use -1 to read all rows
#'   (default is -1).
#'
#' @return A data.table containing the embedding vectors loaded from the compressed CSV file.
#'
#' @examples
#' \dontrun{
#' # Example: Read all embeddings from a module directory
#' all_embeddings <- read_embeddings("module1")
#'
#' # Example: Read only the first 100 rows of embeddings
#' sample_embeddings <- read_embeddings("module2", start_row = 1, n_rows = 100)
#'
#' }
#'
#' @importFrom data.table fread
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
#从embedding储存文件中读取数据，到时候分批读取并计算similarity。
read_embeddings <- function(save_dir, start_row = 1, n_rows = -1) {
  # 拼接路径，确保路径前缀是 "embedding_output/"
  file_path <- paste0(save_dir, "/embedding_vector.csv.gz")

  # 读取文件
  embeddings <- data.table::fread(file_path,
                                  nrows = n_rows,
                                  skip = start_row - 1,
                                  nThread = 4)

  return(embeddings)
}


#' Read Titles
#'
#' Reads paper titles from a specified directory. This function supports batch reading of data.
#'
#' @param save_dir A character string specifying the output directory where the titles file is stored.
#' @param start_row An integer indicating the starting row to read from (default is 1).
#' @param n_rows An integer specifying the number of rows to read. Use -1 to read all rows (default is -1).
#' @return A character vector containing the titles.
#'
#' @examples
#' \dontrun{
#' # Example: Read the first 50 titles
#' titles <- read_titles("example_dir", start_row = 1, n_rows = 50)
#' }
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
# 读取title
read_titles <- function(save_dir, start_row = 1, n_rows = -1) {
  # 拼接路径，确保路径前缀是 "embedding_output/"
  file_path <- paste0(save_dir, "/paper_title.txt")

  # 读取文件
  titles <- readr::read_lines(file = file_path, skip = start_row - 1, n_max = n_rows)

  return(titles)
}

#' Read Chunks
#'
#' Reads abstract or chunk data from a specified directory. This function supports batch reading of data.
#'
#' @param save_dir A character string specifying the output directory where the chunks file is stored.
#' @param start_row An integer indicating the starting row to read from (default is 1).
#' @param n_rows An integer specifying the number of rows to read. Use -1 to read all rows (default is -1).
#'
#' @return A character vector containing the chunks.
#'
#' @examples
#' \dontrun{
#' # Example: Read the first 30 chunks
#' chunks <- read_chunks("example_dir", start_row = 1, n_rows = 30)
#' }
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal

read_chunks <- function(save_dir, start_row = 1, n_rows = -1) {
  # 拼接路径，确保路径前缀是 "embedding_output/"
  file_path <- paste0(save_dir, "/chunk.txt")

  # 读取文件
  chunks <- readr::read_lines(file = file_path, skip = start_row - 1, n_max = n_rows)

  return(chunks)
}


#' Generate Module Embedding
#'
#' Generates an embedding for a given module based on its pathway and gene descriptions.
#'
#' @param module_list A list or data frame containing module information with at least:
#'   \code{PathwayNames} (character vector) and one of \code{GeneNames_vec} or \code{MetNames_vec}
#'   (character vectors defining molecules of interest).
#' @param embedding_model A string specifying the embedding model to use (default is `"text-embedding-3-small"`).
#' @param api_key A character string containing the API key for the embedding generation service.
#'
#' @return A numeric vector representing the generated embedding for the module.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
#基于pathway和gene为module生成embedding
get_module_embedding <- function(module_list, api_key, embedding_model = "text-embedding-3-small"){
  if ("GeneNames_vec" %in% names(module_list)) {
    module_str <- paste0(
      paste(module_list$PathwayNames, collapse = " "),
      " ",
      paste(module_list$GeneNames_vec, collapse = " ")
    )
  } else if ("MetNames_vec" %in% names(module_list)) {
    module_str <- paste0(
      paste(module_list$PathwayNames, collapse = " "),
      " ",
      paste(module_list$MetNames_vec, collapse = " ")
    )
  }

  get_embedding(module_str, api_key, model_name = embedding_model)
}


#' Calculate Similarity Between Embeddings
#'
#' Calculates the cosine similarity between a target embedding and a list of embeddings.
#'
#' @param target_embeddings_list A matrix or data frame where each row represents an embedding.
#' @param module_embedding A numeric vector representing the module embedding to compare against.
#'
#' @return A numeric vector containing similarity scores for each row in `target_embeddings_list`.
#'
#' @examples
#' \dontrun{
#' # Example: Calculate similarity scores
#' target_embeddings <- matrix(runif(100), nrow = 10) # 10 embeddings with random values
#' module_embedding <- runif(10) # A random module embedding
#' similarity_scores <- calculate_similarity(target_embeddings, module_embedding)
#' }
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
#计算embedding间的相似性
calculate_similarity <- function(target_embeddings_list, module_embedding) {
  # 检查输入的维度是否匹配
  if (ncol(target_embeddings_list) != length(module_embedding)) {
    stop("Dimension mismatch between target embeddings and module embedding.")
  }

  # 计算余弦相似性
  cosine_similarity <- function(vec1, vec2) {
    return(sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2))))
  }

  # 遍历 target_embeddings_list 中的每一行
  similarity_scores <- apply(target_embeddings_list, 1, function(row) {
    cosine_similarity(row, module_embedding)
  })

  # 返回相似性分数的列表
  return(similarity_scores)
}



#' Process Chunks with GPT for Relevance Scoring and Text Cleaning
#'
#' @description
#' Processes a list of text chunks using GPT to evaluate their relevance to a specified module
#' and cleans the text by removing unrelated information. The function assigns relevance scores
#' to each chunk based on how well it aligns with the pathways and molecules defined in the module.
#'
#' @param chunks A character vector where each element is a text chunk (e.g., abstracts or articles) to process.
#' @param module_list A list or data frame containing module information with at least:
#'   \code{PathwayNames} (character vector) and one of \code{GeneNames_vec} or \code{MetNames_vec}
#'   (character vectors defining molecules of interest).
#' @param api_key A character string containing the API key for the GPT service.
#' @param model A string specifying the GPT model to use (default is `"gpt-4o-mini-2024-07-18"`).
#'
#' @return A list of results where each element is a list containing:
#' \item{relevance_score}{A numeric value between 0 and 1, indicating the relevance of the chunk.}
#' \item{cleaned_text}{A character string with unrelated information (such as author names,
#'   affiliations, and non-relevant metadata) removed.}
#' The list is sorted in descending order of \code{relevance_score}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts pathway and molecule information from the module_list
#'   \item Processes each chunk in parallel using either \code{parallel::parLapply} (Windows)
#'     or \code{pbmclapply} (other platforms)
#'   \item For each chunk, calls the GPT API with a carefully constructed prompt
#'   \item Validates and potentially retries the API call if the response format is incorrect
#'   \item Parses the JSON response to extract relevance scores and cleaned text
#'   \item Returns the results sorted by relevance score
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Process a set of scientific abstracts
#' abstracts <- c(
#'   "Abstract 1: This study investigates Pathway1 and its relationship with Gene1...",
#'   "Abstract 2: Recent findings on Pathway2 suggest that Gene2 plays a crucial role..."
#' )
#' module_info <- list(
#'   PathwayNames = c("Pathway1", "Pathway2"),
#'   GeneNames_vec = c("Gene1", "Gene2", "Gene3")
#' )
#' api_key <- "your_openai_api_key"
#' results <- GPT_process_chunk(abstracts, module_info, api_key, model)
#'
#' # Access the most relevant result
#' top_result <- results[[1]]
#' print(paste("Top score:", top_result$relevance_score))
#' print(paste("Cleaned text:", top_result$cleaned_text))
#' }
#'
#' @importFrom parallel makeCluster clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom jsonlite fromJSON
#'
#' @seealso
#' \code{\link{gpt_api_call}} used internally for API communication
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
GPT_process_chunk <- function(chunks, module_list, api_key, model = "gpt-4o-mini-2024-07-18") {
  # 初始化结果列表
  reranked_results <- list()

  # 提取模块相关信息
  pathways <- paste(module_list[["PathwayNames"]], collapse = ", ")
  if ("GeneNames_vec" %in% names(module_list)) {
    molecules <- paste(module_list[["GeneNames_vec"]], collapse = ", ")
  } else if ("MetNames_vec" %in% names(module_list)) {
    molecules <- paste(module_list[["MetNames_vec"]], collapse = ", ")
  }

  # 根据操作系统选择不同的并行处理方式
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(min(detectCores()-1, 10))

    # 导出所有必需的变量到集群
    parallel::clusterExport(cl, c("chunks", "pathways", "molecules", "api_key", "model"),
                           envir = environment())

    # 导出工具函数
    parallel::clusterExport(cl, c("process_chunk", "check_json_format",
                                 "modify_prompt_for_format"),
                           envir = environment())

    # 加载必要的包
    parallel::clusterEvalQ(cl, {
      library(jsonlite)
      library(curl)
      source("R/17_llm_module_utils.R")  # 导入包含 gpt_api_call 的文件
    })

    # 执行并行处理
    reranked_results <- parallel::parLapply(cl, chunks,
                                          function(chunk) {
                                            process_chunk(chunk, pathways, molecules, api_key, model = model)
                                          })
    parallel::stopCluster(cl)
  } else {
    reranked_results <- pbmclapply(chunks, process_chunk,
                                  pathways = pathways,
                                  molecules = molecules,
                                  api_key = api_key,
                                  model = model,
                                  mc.cores = detectCores()-1)
  }

  # 对结果按 relevance_score 进行降序排序
  reranked_results <- reranked_results[order(sapply(reranked_results, function(x) x$relevance_score), decreasing = TRUE)]

  return(reranked_results)
}

# 定义验证 JSON 格式的函数
check_json_format <- function(response) {
  tryCatch({
    result <- jsonlite::fromJSON(response)
    if (!is.null(result$relevance_score) && !is.null(result$cleaned_text)) {
      return(TRUE)
    }
  }, error = function(e) {
    return(FALSE)
  })
  return(FALSE)
}

# 定义修改格式的逻辑
modify_prompt_for_format <- function(original_prompt) {
  # 将原始 prompt 和附加信息封装为 messages 格式
  messages <- list(
    list(role = "system", content = "You are an AI tasked with ensuring responses follow a strict JSON format."),
    list(role = "user", content = paste0(
      original_prompt, "\n\n",
      "If your response does not strictly follow the JSON format, please fix the format and make sure to return a valid JSON structure like this:\n",
      "{\n",
      "  \"relevance_score\": <score>,\n",
      "  \"cleaned_text\": \"<cleaned text>\"\n",
      "}"
    ))
  )

  return(messages)
}

# 遍历 chunks，调用 GPT API
process_chunk <- function(chunk,
                          pathways,
                          molecules,
                          api_key,
                          model = "gpt-4o-mini-2024-07-18") {
  # 构建 GPT API 的 prompt
  messages <- list(
    list(role = "system", content = "You are an AI tasked with identifying the most relevant and valuable articles for the given module."),
    list(role = "user", content = paste0(
      "The module is defined by the following pathways: ", pathways, ". ",
      "The module also focuses on the following molecules: ", molecules, ". ",
      "Below is an abstract or chunk of text. Please:\n",
      "1. Rank its relevance to the module in terms of how well it aligns with the pathways and molecules.\n",
      "2. Remove unrelated information such as author names, affiliations, and non-relevant metadata.\n",
      "3. Return the cleaned text and assign it a relevance score (from 0 to 1, where 1 is highly relevant).\n\n",
      "Text:\n", chunk, "\n\n",
      "Please respond in the following JSON format:\n",
      "{\n",
      "  \"relevance_score\": <score>,\n",
      "  \"cleaned_text\": \"<cleaned text>\"\n",
      "}"
    ))
  )


  # 初次调用 GPT API
  gpt_response <- gpt_api_call(messages, api_key, model = model)

  # 检查格式是否正确
  if (!check_json_format(gpt_response)) {
    # 格式不对，重复调用两次
    gpt_response <- gpt_api_call(messages, api_key, model = model)
    if (!check_json_format(gpt_response)) {
      # 第二次格式依然不对，再次修改 prompt
      messages <- modify_prompt_for_format(messages)
      gpt_response <- gpt_api_call(messages, api_key, model = model)
      if (!check_json_format(gpt_response)) {
        # 三次调用仍然失败，返回一个固定格式的默认失败结果
        gpt_response <- '{"relevance_score": 0, "cleaned_text": "Unable to process the request."}'
      }
    }
  }

  # 解析 GPT 返回的 JSON 格式结果
  result <- jsonlite::fromJSON(gpt_response)

  # 返回解析后的结果
  return(list(relevance_score = result$relevance_score, cleaned_text = result$cleaned_text))
}

#' Retrieve and Rerank Strategy for Modules
#'
#' @description
#' Retrieves embeddings for PubMed and optionally local documents for each module,
#' calculates similarity scores, and uses GPT to rerank and clean the text.
#'
#' @param pubmed_result A named list where each element corresponds to a module,
#'   containing \code{PathwayNames} and one of \code{GeneNames_vec} or \code{MetNames_vec}.
#' @param model A string specifying the GPT model to use (default is `"gpt-4o-mini-2024-07-18"`).
#' @param embedding_model A string specifying the embedding model to use (default is `"text-embedding-3-small"`).
#' @param api_key A character string containing the API key for the GPT service.
#' @param similarity_filter_num An integer specifying the number of top documents to keep
#'   after similarity filtering (default is 15).
#' @param GPT_filter_num An integer specifying the number of top documents to keep
#'   after GPT reranking (default is 5).
#' @param local_corpus A logical value indicating whether to include local documents
#'   in the analysis (default is FALSE).
#' @param embedding_output_dir A character string specifying the directory for embeddings.
#' @param save_dir_local_corpus_embed A character string specifying the subdirectory for local
#'   corpus embeddings (required if \code{local_corpus} is TRUE).
#'
#' @return A named list where each element corresponds to a module. Each module contains
#'   a list of results with the following components:
#' \item{title}{The title of the filtered document.}
#' \item{relevance_score}{The relevance score assigned by GPT.}
#' \item{cleaned_text}{The cleaned text returned by GPT.}
#' If no documents are found or processed for a module, the module's entry will be \code{NULL}.
#'
#' @details
#' This function performs the following steps for each module:
#' \enumerate{
#'   \item Retrieves module embeddings using \code{\link{get_module_embedding}}.
#'   \item Calculates similarity scores for PubMed and optionally local documents
#'     using \code{\link{calculate_similarity}}.
#'   \item Filters the top documents based on similarity scores.
#'   \item Reads content and titles for the top documents.
#'   \item Uses \code{\link{GPT_process_chunk}} to rerank and clean the top documents.
#'   \item Merges and returns the top results, sorted by relevance scores.
#' }
#'
#' Progress updates are printed to the console throughout the processing of each module.
#'
#'
#' @seealso
#' \code{\link{get_module_embedding}},
#' \code{\link{read_embeddings}},
#' \code{\link{calculate_similarity}},
#' \code{\link{read_titles}},
#' \code{\link{read_chunks}},
#' \code{\link{GPT_process_chunk}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
#对于每个module进行检索，检索后通过GPT进行rerank然后顺便清洗一下
retrieve_strategy <- function(pubmed_result,
                              model = "gpt-4o-mini-2024-07-18",
                              embedding_model = "text-embedding-3-small",
                              api_key,
                              similarity_filter_num = 15,
                              GPT_filter_num = 5,
                              local_corpus = FALSE,
                              embedding_output_dir = NULL,
                              save_dir_local_corpus_embed = NULL) {
  result <- list()

  cat("Starting to process all modules...\n")
  total_modules <- length(names(pubmed_result))
  current_module <- 0


  for (module_name in names(pubmed_result)) {

    current_module <- current_module + 1

    cat(sprintf("\nProcessing module %d/%d: %s\n", current_module, total_modules, module_name))

    module_list <- pubmed_result[[module_name]]

    if (length(module_list$PubmedIDs) == 0) {
      cat("- No PubMed id available. \n")
      final_result <- list()
    } else {
      cat("- Generating module embeddings...\n")
      module_embedding <- get_module_embedding(module_list = module_list, embedding_model = embedding_model, api_key = api_key)

      cat("- Reading and calculating PubMed similarity...\n")
      pubmed_embeddings <- read_embeddings(save_dir = file.path(embedding_output_dir, module_name))
      pubmed_similarity <- calculate_similarity(pubmed_embeddings, module_embedding)

      # 从内存中释放 pubmed_embeddings
      rm(pubmed_embeddings)
      gc()  # 强制垃圾回收

      if (local_corpus) {
        # 检查 local_path 是否存在
        local_path <- file.path(embedding_output_dir, save_dir_local_corpus_embed)
        if (dir.exists(local_path)) {
          cat("- Reading and calculating local document similarity...\n")
          local_embeddings <- read_embeddings(save_dir = local_path)
          local_similarity <- calculate_similarity(local_embeddings, module_embedding)

          rm(local_embeddings)
          gc()
        }
      } else {
        cat("- Local documents not found, skipping local similarity calculation\n")
        local_similarity <- NULL
      }

      # 合并 PubMed 和 Local 相似性分数
      all_similarities <- c(pubmed_similarity, local_similarity)

      # 获取相似性 Top N 的索引（由 similarity_filter_num 决定）
      if (!is.null(all_similarities)) {
        top_indices <- order(all_similarities, decreasing = TRUE)[1:similarity_filter_num]
      } else {
        top_indices <- NULL
      }

      # 根据 Top N 的索引区分 PubMed 和 Local 的索引
      if (!is.null(top_indices)) {
        pubmed_count <- length(pubmed_similarity)
        pubmed_top_indices <- top_indices[top_indices <= pubmed_count]
        local_top_indices <- top_indices[top_indices > pubmed_count] - pubmed_count
      }

      # 读取 PubMed 和 Local 的标题
      pubmed_selected_titles <- read_titles(save_dir = file.path(embedding_output_dir, module_name))
      pubmed_filtered_titles <- pubmed_selected_titles[pubmed_top_indices]
      local_selected_titles <- NULL
      local_filtered_titles <- NULL

      cat("- Reading document content...\n")
      pubmed_selected_chunks <- read_chunks(save_dir = file.path(embedding_output_dir, module_name))
      pubmed_filtered_chunks <- pubmed_selected_chunks[pubmed_top_indices]
      rm(pubmed_selected_chunks)
      gc()

      # 只在local文件存在时读取local内容
      if (local_corpus) {
        if (dir.exists(local_path)) {
          cat("- Reading local document content...\n")
          local_selected_chunks <- read_chunks(save_dir = local_path)
          local_filtered_chunks <- local_selected_chunks[local_top_indices]
          rm(local_selected_chunks)
          gc()
          local_selected_titles <- read_titles(save_dir = local_path)
          local_filtered_titles <- local_selected_titles[local_top_indices]

          # 合并 PubMed 和 Local 的数据
          combined_chunks <- c(pubmed_filtered_chunks, local_filtered_chunks)
          combined_titles <- c(pubmed_filtered_titles, local_filtered_titles)
        }
      } else {
        # 如果没有local数据，只使用PubMed数据
        combined_chunks <- pubmed_filtered_chunks
        combined_titles <- pubmed_filtered_titles
      }

      cat("- Processing document content using GPT...\n")
      combined_GPT_result <- GPT_process_chunk(combined_chunks, module_list, api_key, model = model)

      cat("- Filtering the most relevant results...\n")
      # 筛选 GPT 结果 Top N（由 GPT_filter_num 决定）
      if (!is.null(combined_GPT_result)) {
        # 提取 relevance_score 列表
        relevance_scores <- sapply(combined_GPT_result, function(x) x$relevance_score)

        # 对 relevance_scores 排序，获取 Top N 索引
        GPT_top_indices <- order(relevance_scores, decreasing = TRUE)[1:GPT_filter_num]

        # 根据排序后的索引提取对应的结果
        filtered_GPT_result <- combined_GPT_result[GPT_top_indices]

        # 提取对应的标题
        filtered_titles <- combined_titles[GPT_top_indices]
      } else {
        filtered_GPT_result <- NULL
        filtered_titles <- NULL
      }

      # 如果有筛选结果，按项合并标题和 GPT 结果
      if (!is.null(filtered_GPT_result) && !is.null(filtered_titles)) {
        # 按项合并，每个元素包含 title 和对应的 GPT_result
        final_result <- Map(function(gpt, title) {
          list(
            title = title,
            relevance_score = gpt$relevance_score,
            cleaned_text = gpt$cleaned_text
          )
        }, filtered_GPT_result, filtered_titles)
      } else {
        # 如果没有结果，返回 NULL
        final_result <- list()
      }
    }
    # 将结果存储到 result 列表中
    result[[module_name]] <- final_result
    cat(sprintf("%s processing completed\n", module_name))
  }

  cat("\nAll modules processing completed!\n")
  return(result)
}

