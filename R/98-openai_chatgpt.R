#' Set ChatGPT API Key
#'
#' This function sets the ChatGPT API key as an environment variable.
#' It's used to authenticate and interact with ChatGPT API services.
#'
#' @param api_key A string representing the ChatGPT API key to be set.
#'
#' @return Invisible NULL. The function is used for its side effect
#' of setting an environment variable and does not return a value.
#'
#' @examples
#' \dontrun{
#' # Set the ChatGPT API key (example key used here)
#' set_chatgpt_api_key("your_api_key_here")
#'}
#' @export

set_chatgpt_api_key <- function(api_key) {
  api_key <- Sys.setenv(chatgpt_api_key = api_key)
}

#' Request Response from ChatGPT
#'
#' This function sends a request to OpenAI's ChatGPT service using the provided prompt.
#' It requires that a ChatGPT API key is already set using `set_chatgpt_api_key()`.
#' The function makes a POST request to the ChatGPT API and returns the generated response.
#'
#' @param prompt A string representing the prompt to be sent to the ChatGPT service.
#' Default prompt is "what is metabolomics?".
#'
#' @return A string containing the response from ChatGPT.
#' If the API key is not set, the function stops with an error message.
#'
#' @examples
#' \dontrun{
#'   # Set the ChatGPT API key first
#'   set_chatgpt_api_key("your_api_key_here")
#'
#'   # Request a response from ChatGPT service
#'   chatgpt_response <- request_chatgpt_response(prompt = "Tell me about R programming")
#' }
#'
#' @export
#' @importFrom httr POST add_headers content_type_json content

request_chatgpt_response <-
  function(prompt = "what is metabolomics?") {
    api_key <-
      Sys.getenv("chatgpt_api_key")
    if (nchar(api_key) == 0) {
      stop("Please set the API key first using the function `set_chatgpt_api_key()`")
    }
    model_query <- "chatgpt-pro:generateContent"

    response <-
      httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        httr::add_headers(Authorization = paste("Bearer", api_key)),
        httr::content_type_json(),
        encode = "json",
        body = list(model = "gpt-4-1106-preview",
                    messages = list(list(
                      role = "user", content = prompt
                    )))
      )

    choices <- httr::content(response)$choices[[1]]
    outputs <- choices$message$content
    return(outputs)
  }
