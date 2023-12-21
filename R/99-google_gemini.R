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


#' Translate Text into Different Languages
#'
#' This function translates a given text into a specified language using the Gemini translation engine.
#' It utilizes the Gemini large language models for translation.
#'
#' @param text A string representing the text to be translated.
#' Default is a description of the Gemini model.
#' @param engine A string specifying the translation engine to be used.
#' Currently, only 'gemini' is supported.
#' @param to A string indicating the target language for translation.
#' Supported languages include 'chinese', 'spanish', 'english', 'french', 'german',
#' 'italian', 'japanese', 'korean', 'portuguese', 'russian', and 'spanish'.
#'
#' @return A character vector containing the translated text.
#'
#' @examples
#' \dontrun{
#'   # Translate text to French using Gemini
#'   translated_text <- translate_language(text = "Hello World", engine = "gemini", to = "french")
#' }
#'
#' @export

translate_language <-
  function(text = "Gemini is a family of multimodal large language models
           developed by Google DeepMind,
           serving as the successor to LaMDA and PaLM 2.
           Comprising Gemini Ultra, Gemini Pro, and Gemini Nano,
           it was announced on December 6, 2023,
           positioned as a contender to OpenAI's GPT-4.",
           engine = c("gemini"),
           to = c(
             "chinese",
             "spanish",
             "english",
             "french",
             "german",
             "italian",
             "japanese",
             "korean",
             "portuguese",
             "russian",
             "spanish"
           )) {
    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    if (engine == "gemini") {
      prompt <- paste0("Translate to ", to, ":\n", text)
      request_gemini_response(prompt = prompt)
    }
  }
