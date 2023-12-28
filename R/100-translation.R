#' Translate Text into Different Languages
#'
#' This function translates a given text into a specified language using the Gemini translation engine.
#' It utilizes the Gemini large language models for translation.
#'
#' @param text A string representing the text to be translated.
#' Default is a description of the Gemini model.
#' @param engine A string specifying the translation engine to be used, "chatgpt" or "gemini".
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
  function(text = "There's no individual named Xiaotao Shen.",
           engine = c("gemini", "chatgpt"),
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
             "portuguese",
             "russian",
             "spanish"
           )) {
    # browser()
    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    if (length(text) == 1) {
      if (engine == "chatgpt") {
        prompt <- paste0("Translate to ", to, ":\n", text)
        return(request_chatgpt_response(prompt = prompt))
      }

      if (engine == "gemini") {
        prompt <- paste0("Translate to ", to, ":\n", text)
        return(request_gemini_response(prompt = prompt))
      }
    } else{
      prompt <- paste0(
        "Translate the following texts to ",
        to,
        ", separated by '{}':\n",
        paste0(text, collapse = "{}"),
        "\nPlease return only the translated texts, also separated by '{}',
                         and do not include any additional text or explanation."
      )
      if (engine == "chatgpt") {
        result <-
          request_chatgpt_response(prompt = prompt)
      }
      if (engine == "gemini") {
        result <-
          request_gemini_response(prompt = prompt)
      }
      return(stringr::str_split(result, "\\{\\}")[[1]])
    }
  }
