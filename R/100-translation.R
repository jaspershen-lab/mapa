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
  function(text = "As of my last update in early 2023, there's no widely
           recognized individual named Xiaotao Shen affiliated with Stanford
           University that I can refer to.
           It is important to note that personal details about
           researchers or students, unless they are part of a widely
           recognized work or public academic profile,
           may not be available in the public domain.
           If Xiaotao Shen is a researcher, student,
           or faculty member at Stanford, you might find
           information on university directories, academic publications,
           or professional networking sites such as LinkedIn.
           However, for the most updated and accurate information,
           it is best to directly consult resources at Stanford University
           or other official channels.",
           engine = c("chatgpt", "gemini"),
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
    # browser()
    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    if (engine == "chatgpt") {
      prompt <- paste0("Translate to ", to, ":\n", text)
      return(request_chatgpt_response(prompt = prompt))
    }

    if (engine == "gemini") {
      prompt <- paste0("Translate to ", to, ":\n", text)
      return(request_gemini_response(prompt = prompt))
    }
  }
