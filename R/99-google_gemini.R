#' Set Gemini API Key
#'
#' This function sets the Gemini API key as an environment variable.
#' It's used to authenticate and interact with Gemini API services.
#'
#' @param api_key A string representing the Gemini API key to be set.
#'
#' @return Invisible NULL. The function is used for its side effect
#' of setting an environment variable and does not return a value.
#'
#' @examples
#' \dontrun{
#' # Set the Gemini API key (example key used here)
#' set_gemini_api_key("your_api_key_here")
#'}
#' @export

set_gemini_api_key <- function(api_key) {
  api_key <- Sys.setenv(gemini_api_key = api_key)
}

#' Request Response from Gemini Service
#'
#' This function sends a request to the Gemini service using the provided prompt.
#' It requires that a Gemini API key is already set using `set_gemini_api_key()`.
#' The function makes a POST request and returns the response from Gemini.
#'
#' @param prompt A string representing the prompt to be sent to the Gemini service.
#'
#' @return A character vector containing the response outputs from the Gemini service.
#' If the API key is not set, the function stops with an error message.
#'
#' @examples
#' \dontrun{
#'   # Set the API key first
#'   set_gemini_api_key("your_api_key_here")
#'
#'   # Request a response from Gemini service
#'   response <- request_gemini_response("Your prompt here")
#' }
#'
#' @export
#' @importFrom httr POST content_type_json content

request_gemini_response <-
  function(prompt = "what is metabolomics?") {
    api_key <-
      Sys.getenv("gemini_api_key")
    if (nchar(api_key) == 0) {
      stop("Please set the API key first using the function `set_gemini_api_key()`")
    }
    model_query <- "gemini-pro:generateContent"

    response <- httr::POST(
      url = paste0(
        "https://generativelanguage.googleapis.com/v1beta/models/",
        model_query
      ),
      query = list(key = api_key),
      httr::content_type_json(),
      encode = "json",
      body = list(
        contents = list(parts = list(list(text = prompt))),
        generationConfig = list(temperature = 0.5,
                                maxOutputTokens = 1024)
      )
    )
    candidates <- httr::content(response)$candidates
    outputs <-
      unlist(lapply(candidates, function(candidate)
        candidate$content$parts))

    names(outputs) <- NULL

    return(outputs)
  }
