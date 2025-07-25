# #This file includes general utils
# library(httr)
# library(httr2)
# library(curl)
# library(jsonlite)
# library(roxygen2)

#' gpt_api_call: Function for calling GPT model APIs (supports OpenAI and Gemini providers)
#'
#' This internal function is designed to interact with GPT-based APIs (such as OpenAI's GPT or the Gemini API).
#' It handles HTTP requests, retries failed attempts, and parses responses.
#' This function is not intended for direct use by package users.
#'
#' @param messages A list of messages to send to the GPT model, usually representing the conversation history.
#'                  Each message should be a list with `role` (e.g., "system", "user", or "assistant") and `content`.
#' @param api_key A string containing the API key required for authentication.
#' @param model A string specifying the GPT model to use (default is `"gpt-4o-mini-2024-07-18"`).
#' @param max_tokens An integer indicating the maximum number of tokens to generate in the response (default is `1000`).
#' @param temperature A numeric value between 0 and 1 to control the randomness of the response
#'                    (default is `0.7`, where lower values produce more deterministic results).
#' @param retry_attempts An integer specifying the maximum number of retry attempts if the API call fails (default is `3`).
#' @param api_provider A string indicating the API provider, either `"openai"`, `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param thinkingBudget An integer for the "thinking budget" parameter specific to the Gemini API (default is `0`).
#'
#' @return The generated text from the GPT model if the call is successful. If the call fails after the specified number
#'          of retries, the function returns `NULL`.
#'
#' @note Ensure you have a valid API key and the correct endpoint for your chosen provider. For Gemini, update the URL
#'       in the code if necessary.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'

#' @keywords internal
gpt_api_call <- function(
    messages,
    api_key,
    model = NULL,
    max_tokens = 1000,
    temperature = 0.7,
    retry_attempts = 3,
    api_provider = "openai",
    thinkingBudget = 0
) {
  # 1. 设置默认模型
  if (is.null(model)) {
    model <- switch(
      api_provider,
      "openai"      = "gpt-4o-mini-2024-07-18",
      "gemini"      = "gemini-2.5-flash",
      "siliconflow" = "Qwen/Qwen3-32B",
      stop("Invalid api_provider. Choose 'openai', 'gemini', or 'siliconflow'.")
    )
  }

  # 2. URL

  api_url <- switch(
    api_provider,
    "openai"      = "https://api.openai.com/v1/chat/completions",
    "gemini"      = paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                           model, ":generateContent?key=", api_key),
    "siliconflow" = "https://api.siliconflow.cn/v1/chat/completions",
    stop("Invalid api_provider. Choose 'openai', 'gemini', or 'siliconflow'.")
  )

  # 3. 构建请求体
  if (api_provider %in% c("openai", "siliconflow")) {
    request_body <- jsonlite::toJSON(list(
      model = model,
      messages = messages,
      max_tokens = max_tokens,
      temperature = temperature,
      n = 1
    ), auto_unbox = TRUE)
  } else if (api_provider == "gemini") {
    gemini_contents <- lapply(messages, function(msg) {
      if (msg$role == "user") {
        list(parts = list(list(text = msg$content)))
      } else if (msg$role == "assistant") {
        list(role = "model", parts = list(list(text = msg$content)))
      } else {
        NULL
      }
    })
    gemini_contents <- gemini_contents[!sapply(gemini_contents, is.null)]

    gen_cfg <- list(
      maxOutputTokens = max_tokens,
      temperature = temperature
    )
    if (thinkingBudget > 0) {
      gen_cfg$thinkingConfig <- list(thinkingBudget = thinkingBudget)
    }

    request_body <- jsonlite::toJSON(list(
      contents = gemini_contents,
      generationConfig = gen_cfg
    ), auto_unbox = TRUE, null = "null")
  } else {
    stop("Invalid api_provider.")
  }

  attempt <- 0
  repeat {
    attempt <- attempt + 1
    handle <- curl::new_handle()

    if (api_provider %in% c("openai", "siliconflow")) {
      curl::handle_setheaders(
        handle,
        Authorization = paste("Bearer", api_key),
        `Content-Type` = "application/json"
      )
    } else if (api_provider == "gemini") {
      curl::handle_setheaders(handle, `Content-Type` = "application/json")
    }

    curl::handle_setopt(
      handle,
      url = api_url,
      postfields = request_body,
      customrequest = "POST",
      timeout = 100
    )

    response <- tryCatch(
      curl::curl_fetch_memory(api_url, handle),
      error = function(e) {
        message(sprintf("Attempt %d failed: %s", attempt, e$message))
        NULL
      }
    )

    if (is.null(response)) {
      if (attempt >= retry_attempts) {
        message("Max retries reached, returning NULL.")
        return(NULL)
      }
      Sys.sleep(1); next
    }

    if (response$status_code == 200) {
      resp <- jsonlite::fromJSON(rawToChar(response$content))
      llm_output <- switch(
        api_provider,
        "openai" = resp$choices$message$content,
        "siliconflow" = resp[["choices"]][["message"]][["content"]],
        "gemini" = resp[["candidates"]][["content"]][["parts"]][[1]][["text"]],
        stop("Unexpected provider in parsing.")
      )
      return(llm_output)
    } else {
      message(sprintf(
        "Attempt %d failed with status %s: %s",
        attempt, response$status_code, rawToChar(response$content)
      ))
      if (attempt >= retry_attempts) {
        warning("Failed after ", retry_attempts, " attempts.")
        return(NULL)
      }
      Sys.sleep(1)
    }
  }
}

#' get_embedding: Internal function to retrieve embeddings from an API
#'
#' This internal function retrieves text embeddings from a specified API provider (e.g., OpenAI, Gemini, or SiliconFlow).
#' It supports error handling, retries, and flexible model configurations. The function returns the embedding
#' vector for a given input text chunk.
#'
#' @param chunk A character string representing the input text for which the embedding is required.
#' @param api_key A string containing the API key required for authentication with the provider.
#' @param model_name A string specifying the embedding model to use.
#'   - For OpenAI, default is `"text-embedding-3-small"`.
#'   - For Gemini, default is `"models/gemini-embedding-exp-03-07"`.
#'   - For SiliconFlow, default is `"BAAI/bge-large-zh-v1.5"`.
#' @param api_provider A string indicating the API provider, either `"openai"`, `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param task_type A string specifying the task type for Gemini API (default is `"SEMANTIC_SIMILARITY"`). Only used when `api_provider` is "gemini".
#'
#' @return A numeric vector representing the embedding for the input text. If the API call fails after retries,
#'         the function returns `NULL`.
#'
#' @note This function is for internal use only and supports retry logic for API failures. Ensure you have
#'       valid API credentials and the provider's URL is correctly set.
#'
#' @author Feifan Zhang <FEIFAN004@e.ntu.edu.sg>
#'
#' @export
get_embedding <- function(chunk, api_key, model_name = NULL, api_provider = "openai", task_type = "SEMANTIC_SIMILARITY") {

  # 设置默认模型名称
  if (is.null(model_name)) {
    model_name <- switch(
      api_provider,
      "openai" = "text-embedding-3-small",
      "gemini" = "models/gemini-embedding-exp-03-07",
      "siliconflow" = "Qwen/Qwen3-Embedding-8B", # 新增: SiliconFlow 的默认模型
      stop("Invalid API provider. Please choose 'openai', 'gemini', or 'siliconflow'.")
    )
  }

  # 根据 API 提供者选择 URL 和构建请求体
  if (api_provider %in% c("openai", "siliconflow")) { # 修改: 将 openai 和 siliconflow 分组
    url <- switch(
      api_provider,
      "openai" = "https://api.openai.com/v1/embeddings",
      "siliconflow" = "https://api.siliconflow.cn/v1/embeddings" # 新增: SiliconFlow URL
    )
    data <- list(
      model = model_name,
      input = chunk
    )
  } else if (api_provider == "gemini") {
    url <- paste0("https://generativelanguage.googleapis.com/v1beta/models/",
                  gsub("^models/", "", model_name), ":embedContent")
    data <- list(
      model = model_name,
      content = list(
        parts = list(
          list(text = chunk)
        )
      ),
      taskType = task_type
    )
  } else {
    stop("Invalid API provider. Please choose 'openai', 'gemini', or 'siliconflow'.")
  }

  # 请求嵌入向量
  embedding <- tryCatch(
    expr = {
      # 创建请求对象
      req <- httr2::request(url)

      # 根据API提供者设置认证方式
      if (api_provider %in% c("openai", "siliconflow")) { # 修改: 为 siliconflow 使用 Bearer Token
        req <- req %>% httr2::req_auth_bearer_token(token = api_key)
      } else if (api_provider == "gemini") {
        req <- req %>% httr2::req_headers(`x-goog-api-key` = api_key)
      }

      # 执行请求
      resp <- req %>%
        httr2::req_body_json(data = data) %>%
        httr2::req_retry(
          max_tries = 3,
          max_seconds = 60,
          is_transient = \(resp) httr2::resp_status(resp) %in% c(429, 500, 503), # 改进重试条件
          after = \(resp) httr2::resp_retry_after(resp)
        ) %>%
        httr2::req_perform()

      # 根据API提供者解析响应
      if (api_provider %in% c("openai", "siliconflow")) { # 修改: 为 siliconflow 使用相同的解析逻辑
        embedding <- httr2::resp_body_json(resp)$data[[1]]$embedding
      } else if (api_provider == "gemini") {
        embedding <- httr2::resp_body_json(resp)$embedding$values
      }

      unlist(embedding)
    },
    error = function(e) {
      warning("Failed to get embedding after retries for provider '", api_provider, "': ", e$message)
      NULL
    }
  )

  return(embedding)
}

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------


#' Clear the embedding_output directory
#'
#' This internal function checks for the existence of a specified directory.
#' If the directory exists, it recursively deletes all files and subdirectories within it,
#' and then removes the empty directory itself. If the directory does not exist, it provides
#' a message indicating so.
#'
#' @param output_dir A character string specifying the directory to be checked and cleared.
#'
#' @return This function does not return a value. It prints a success or failure message
#'         to the console based on whether the directory was found and cleared.
#'
#' @note This function is designed for internal use to manage temporary or output files.
#'       Ensure that the specified directory does not contain critical data before clearing.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
clear_output_dir <- function(output_dir = NULL) {
  # Check if directory exists
  if (file.exists(output_dir)) {
    # Get all files and subdirectories in the directory
    files <- list.files(output_dir, full.names = TRUE)

    if (length(files) != 0) {
      # Recursively delete files and subdirectories
      unlink(files, recursive = TRUE)
      cat("Successfully cleared the directory:", output_dir, "\n")
    }

    # # Delete the empty directory
    # unlink(output_dir, recursive = TRUE)
  } else {
    stop("Directory does not exist:", output_dir, "\n")
  }
}
