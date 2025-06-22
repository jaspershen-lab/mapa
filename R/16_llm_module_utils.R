# #This file includes general utils
# library(httr)
# library(httr2)
# library(curl)
# library(jsonlite)
# library(roxygen2)

#' gpt_api_call: Function for calling GPT model APIs (supports OpenAI and Gemini providers)
#'
#' This internal function is designed to interact with GPT-based APIs (such as OpenAI's GPT or a hypothetical Gemini API).
#' It handles HTTP requests, retries failed attempts, and parses responses.
#' This function is not intended for direct use by package users.
#'
#' @param messages A list of messages to send to the GPT model, usually representing the conversation history.
#'                 Each message should be a list with `role` (e.g., "system", "user", or "assistant") and `content`.
#' @param api_key A string containing the API key required for authentication.
#' @param model A string specifying the GPT model to use (default is `"gpt-4o-mini-2024-07-18"`).
#' @param max_tokens An integer indicating the maximum number of tokens to generate in the response (default is `1000`).
#' @param temperature A numeric value between 0 and 1 to control the randomness of the response
#'                    (default is `0.7`, where lower values produce more deterministic results).
#' @param retry_attempts An integer specifying the maximum number of retry attempts if the API call fails (default is `3`).
#' @param api_provider A string indicating the API provider, either `"openai"` or `"gemini"` (default is `"openai"`).
#'
#' @return The generated text from the GPT model if the call is successful. If the call fails after the specified number
#'         of retries, the function returns `NULL`.
#'
#' @note Ensure you have a valid API key and the correct endpoint for your chosen provider. For Gemini, update the URL
#'       in the code if necessary.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
gpt_api_call <- function(messages, api_key, model = "gpt-4o-mini-2024-07-18", max_tokens = 1000,
                         temperature = 0.7, retry_attempts = 3, api_provider = "openai") {
  # 根据 API 提供者选择 URL
  api_url <- switch(
    api_provider,
    "openai" = "https://api.openai.com/v1/chat/completions",
    "gemini" = paste0("https://generativelanguage.googleapis.com/v1beta/models/", model, ":generateContent?key=", api_key),
    stop("Invalid API provider. Please choose 'openai' or 'gemini'.")
  )
  request_body <- NULL
  if (api_provider == "openai") {
    request_body <- jsonlite::toJSON(list(
      model = model,
      messages = messages,
      max_tokens = max_tokens,
      temperature = temperature,
      n = 1
    ), auto_unbox = TRUE)
  } else if (api_provider == "gemini") {
    # Gemini API 的 messages 结构与 OpenAI 不同
    # 它需要一个 'contents' 数组，每个元素包含 'parts' 数组，其中包含 'text'
    # 假设 messages 是一个列表，每个元素是 {role: ..., content: ...}
    # 我们需要将其转换为 Gemini 的格式
    gemini_contents <- lapply(messages, function(msg) {
      if (msg$role == "user") { # Gemini API user role 是 parts 中的 text
        return(list(parts = list(list(text = msg$content))))
      } else if (msg$role == "assistant") { # Gemini API assistant role 是 model 中的 text
        return(list(role = "model", parts = list(list(text = msg$content))))
      } else {
        # 对于其他角色，可以根据实际情况进行处理或忽略
        return(NULL)
      }
    })
    # 过滤掉 NULL 元素（如果存在）
    gemini_contents <- gemini_contents[!sapply(gemini_contents, is.null)]
    
    request_body <- jsonlite::toJSON(list(
      contents = gemini_contents,
      generationConfig = list( # Gemini 的额外配置
        maxOutputTokens = max_tokens,
        temperature = temperature,
        thinkingConfig = list(
          thinkingBudget = 0 # 设置思考预算为 0
        )
      )
    ), auto_unbox = TRUE, null = "null") # auto_unbox = TRUE 用于将单元素向量转换为标量，null = "null" 处理空值
    
  } else {
    stop("Invalid API provider. Request body cannot be constructed.")
  }
  
  # 初始化重试计数
  attempt <- 0
  
  # 循环重试逻辑
  repeat {
    attempt <- attempt + 1
    
    # 使用 curl 发送请求
    handle <- curl::new_handle()
    
    # Gemini API 的 key 直接在 URL 中，OpenAI 的 key 在 Authorization header 中
    if (api_provider == "openai") {
      curl::handle_setheaders(handle,
                              Authorization = paste("Bearer", api_key),
                              `Content-Type` = "application/json")
    } else if (api_provider == "gemini") {
      curl::handle_setheaders(handle, `Content-Type` = "application/json")
      # API Key 已经在 URL 中，所以这里不需要 Authorization header
    }
    
    curl::handle_setopt(handle,
                        url = api_url,
                        postfields = request_body,
                        customrequest = "POST",
                        timeout = 100)
    
    response <- tryCatch({
      curl::curl_fetch_memory(api_url, handle)
    }, error = function(e) {
      # 捕获可能的网络请求错误
      message(paste("Attempt", attempt, "failed with error:", e$message))
      NULL
    })
    
    # 如果请求失败，判断是否需要重试
    if (is.null(response)) {
      if (attempt >= retry_attempts) {
        message("Max retries reached, returning error response.")
        return(NULL)
      }
      Sys.sleep(1) # 短暂等待后重试
      next
    }
    
    # 检查 HTTP 响应状态
    if (response$status_code == 200) {
      # 如果成功返回，解析 JSON 内容
      response_content <- jsonlite::fromJSON(rawToChar(response$content))
      
      # 根据不同的 API 提供者解析响应
      llm_output <- switch(
        api_provider,
        "openai" = response_content[["choices"]][["message"]][["content"]],
        "gemini" = response_content[["candidates"]][["content"]][["parts"]][[1]][["text"]],
        stop("Invalid API provider. Cannot parse response.")
      )
      return(llm_output)
    } else {
      # 打印错误信息
      message(paste("Attempt", attempt, "failed with status code:", response$status_code,
                    "and message:", rawToChar(response$content)))
      
      # 如果超过最大重试次数，停止执行
      if (attempt >= retry_attempts) {
        warning("Failed to call LLM API after ", retry_attempts, " attempts. Last error status code: ", response$status_code)
        return(NULL)
      }
      Sys.sleep(1) # 短暂等待后重试
    }
  }
}


#' get_embedding: Internal function to retrieve embeddings from an API
#'
#' This internal function retrieves text embeddings from a specified API provider (e.g., OpenAI or Gemini).
#' It supports error handling, retries, and flexible model configurations. The function returns the embedding
#' vector for a given input text chunk.
#'
#' @param chunk A character string representing the input text for which the embedding is required.
#' @param api_key A string containing the API key required for authentication with the provider.
#' @param model_name A string specifying the embedding model to use (default is `"text-embedding-3-small"`).
#' @param api_provider A string indicating the API provider, either `"openai"` or `"gemini"` (default is `"openai"`).
#'
#' @return A numeric vector representing the embedding for the input text. If the API call fails after retries,
#'         the function returns `NULL`.
#'
#'
#' @note This function is for internal use only and supports retry logic for API failures. Ensure you have
#'       valid API credentials and the provider's URL is correctly set.
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
get_embedding <- function(chunk, api_key, model_name = "text-embedding-3-small", api_provider = "openai") {
  # 根据 API 提供者选择 URL
  url <- switch(
    api_provider,
    "openai" = "https://api.openai.com/v1/embeddings",
    "gemini" = "https://api.gemini.com/v1/embeddings",
    stop("Invalid API provider. Please choose 'openai' or 'gemini'.")
  )

  # Body specifying model and text
  data <- list(
    model = model_name,  # 可配置的模型名称
    input = chunk        # 输入文本
  )

  # 请求嵌入向量
  embedding <- tryCatch(
    expr = {
      # 创建请求对象
      req <- httr2::request(url)

      resp <- req %>%
        httr2::req_auth_bearer_token(token = api_key) %>%
        httr2::req_body_json(data = data) %>%
        httr2::req_retry(
          max_tries = 3,
          max_seconds = 60,
          after = \(resp) is.null(resp) && resp$status_code != 200  # 重试条件
        ) %>%
        httr2::req_perform()

      # 提取嵌入向量
      embedding <- httr2::resp_body_json(resp)$data[[1]]$embedding   # 解析 JSON 响应
      unlist(embedding)
    },
    error = function(e) {
      warning("Failed to get embedding after 3 retries:", e$message)
      NULL
    }
  )

  return(embedding)
}



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
