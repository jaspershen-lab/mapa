#' Translate Language
#'
#' This function serves as a generic method for translating text using different engines.
#' It requires a specific method to be defined for each type of input (e.g., character, functional_module).
#'
#' @param text Text to be translated.
#' @param engine Translation engine to use, defaults to 'gemini' or 'chatgpt'.
#' @param to Target language for translation, options include 'chinese', 'spanish', 'english',
#'   'french', 'german', 'italian', 'japanese', 'korean', 'portuguese', 'russian', and 'spanish'.
#' @return An object of the same class as `text` with the translated content.
#' @examples
#' \dontrun{
#'   # Translate text to French using Gemini
#'   translated_text <- translate_language(text = "Hello World", engine = "gemini", to = "french")
#' }
#' @export

translate_language <-
  function(text,
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
    if (missing(text)) {
      stop("Please provide a text to translate.")
    }
    UseMethod("translate_language")
  }



#' @method translate_language character
#' @rdname translate_language
#' @export

translate_language.character <-
  function(text,
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

    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    result <-
      translate_language_internal(text = text,
                                  engine = engine,
                                  to = to)
    result
  }



#' @method translate_language functional_module
#' @docType methods
#' @export

translate_language.functional_module <-
  function(text,
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

    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    ####GO------------------------------------------------
    ####GO pathway result
    if (length(text@enrichment_go_result) > 0) {
      message("GO...")
      result_go <-
        text@enrichment_go_result@result %>%
        dplyr::select(ID, Description)

      if (nrow(result_go) > 100) {
        Description_trans <-
          lapply(seq(
            from = 1,
            to = nrow(result_go),
            by = 100
          ),
          function(i)
            c(1:nrow(result_go))[i:min(i + 99, nrow(result_go))]) %>%
          lapply(function(index) {
            result <-
              tryCatch(
                translate_language_internal(
                  text = result_go$Description[index],
                  engine = engine,
                  to = to
                ),
                error = function(e) {
                  return(result_go$Description[index])
                }
              )
            return(result)
          }) %>%
          unlist()
      } else{
        Description_trans <-
          tryCatch(
            translate_language_internal(
              text = result_go$Description,
              engine = engine,
              to = to
            ),
            error = function(e) {
              return(result_go$Description)
            }
          )
      }

      result_go$Description_trans <-
        Description_trans

      text@enrichment_go_result@result <-
        text@enrichment_go_result@result %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_go[, c("ID", "Description_trans")],
                         by = "ID")

    } else{
      result_go <-
        NULL
    }

    ####GO module result
    if (length(text@merged_pathway_go) > 0) {
      ###module result
      module_result_go <-
        text@merged_pathway_go$module_result

      module_result_go <-
        module_result_go %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_go[, c("Description", "Description_trans")],
                         by = c("module_annotation" = "Description")) %>%
        dplyr::rename(module_annotation_trans = Description_trans)

      module_result_go$Description_trans <-
        lapply(module_result_go$node,
               function(x) {
                 result_go$Description_trans[match(stringr::str_split(x, ";")[[1]],
                                                   result_go$ID)] %>%
                   paste(collapse = ";")
               }) %>%
        unlist()

      text@merged_pathway_go$module_result <-
        module_result_go


      ###result_with_module
      result_with_module_go <-
        text@merged_pathway_go$result_with_module

      result_with_module_go <-
        result_with_module_go %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_go[, c("ID", "Description_trans")],
                         by = c("node" = "ID"))

      text@merged_pathway_go$result_with_module <-
        result_with_module_go
    } else{
      module_result_go <- NULL
    }

    ####kegg------------------------------------------------
    ####kegg pathway result
    if (length(text@enrichment_kegg_result) > 0) {
      message("KEGG...")
      result_kegg <-
        text@enrichment_kegg_result@result %>%
        dplyr::select(ID, Description)

      if (nrow(result_kegg) > 100) {
        Description_trans <-
          lapply(seq(
            from = 1,
            to = nrow(result_kegg),
            by = 100
          ),
          function(i)
            c(1:nrow(result_kegg))[i:min(i + 99, nrow(result_kegg))]) %>%
          lapply(function(index) {
            result <-
              tryCatch(
                translate_language_internal(
                  text = result_kegg$Description[index],
                  engine = engine,
                  to = to
                ),
                error = function(e) {
                  return(result_kegg$Description[index])
                }
              )
            return(result)
          }) %>%
          unlist()
      } else{
        Description_trans <-
          tryCatch(
            translate_language_internal(
              text = result_kegg$Description,
              engine = engine,
              to = to
            ),
            error = function(e) {
              return(result_kegg$Description)
            }
          )
      }

      Description_trans <-
        Description_trans[Description_trans != ""]

      result_kegg$Description_trans <-
        Description_trans

      text@enrichment_kegg_result@result <-
        text@enrichment_kegg_result@result %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_kegg[, c("ID", "Description_trans")],
                         by = "ID")

    } else{
      result_kegg <-
        NULL
    }

    ####kegg module result
    if (length(text@merged_pathway_kegg) > 0) {
      ###module result
      module_result_kegg <-
        text@merged_pathway_kegg$module_result

      module_result_kegg <-
        module_result_kegg %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_kegg[, c("Description", "Description_trans")],
                         by = c("module_annotation" = "Description")) %>%
        dplyr::rename(module_annotation_trans = Description_trans)

      module_result_kegg$Description_trans <-
        lapply(module_result_kegg$node,
               function(x) {
                 result_kegg$Description_trans[match(stringr::str_split(x, ";")[[1]],
                                                     result_kegg$ID)] %>%
                   paste(collapse = ";")
               }) %>%
        unlist()

      text@merged_pathway_kegg$module_result <-
        module_result_kegg


      ###result_with_module
      result_with_module_kegg <-
        text@merged_pathway_kegg$result_with_module

      result_with_module_kegg <-
        result_with_module_kegg %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_kegg[, c("ID", "Description_trans")],
                         by = c("node" = "ID"))

      text@merged_pathway_kegg$result_with_module <-
        result_with_module_kegg
    } else{
      module_result_kegg <- NULL
    }

    ####reactome------------------------------------------------
    ####reactome pathway result
    if (length(text@enrichment_reactome_result) > 0) {
      message("Reactome...")
      result_reactome <-
        text@enrichment_reactome_result@result %>%
        dplyr::select(ID, Description)

      if (nrow(result_reactome) > 100) {
        Description_trans <-
          lapply(seq(
            from = 1,
            to = nrow(result_reactome),
            by = 100
          ),
          function(i)
            c(1:nrow(result_reactome))[i:min(i + 99, nrow(result_reactome))]) %>%
          lapply(function(index) {
            result <-
              tryCatch(
                translate_language_internal(
                  text = result_reactome$Description[index],
                  engine = engine,
                  to = to
                ),
                error = function(e) {
                  return(result_reactome$Description[index])
                }
              )
            return(result)
          }) %>%
          unlist()
      } else{
        Description_trans <-
          tryCatch(
            translate_language_internal(
              text = result_reactome$Description,
              engine = engine,
              to = to
            ),
            error = function(e) {
              return(result_reactome$Description)
            }
          )
      }

      result_reactome$Description_trans <-
        Description_trans

      text@enrichment_reactome_result@result <-
        text@enrichment_reactome_result@result %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_reactome[, c("ID", "Description_trans")],
                         by = "ID")

    } else{
      result_reactome <-
        NULL
    }

    ####reactome module result
    if (length(text@merged_pathway_reactome) > 0) {
      ###module result
      module_result_reactome <-
        text@merged_pathway_reactome$module_result

      module_result_reactome <-
        module_result_reactome %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_reactome[, c("Description", "Description_trans")],
                         by = c("module_annotation" = "Description")) %>%
        dplyr::rename(module_annotation_trans = Description_trans)

      module_result_reactome$Description_trans <-
        lapply(module_result_reactome$node,
               function(x) {
                 result_reactome$Description_trans[match(stringr::str_split(x, ";")[[1]],
                                                         result_reactome$ID)] %>%
                   paste(collapse = ";")
               }) %>%
        unlist()

      text@merged_pathway_reactome$module_result <-
        module_result_reactome


      ###result_with_module
      result_with_module_reactome <-
        text@merged_pathway_reactome$result_with_module

      result_with_module_reactome <-
        result_with_module_reactome %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(result_reactome[, c("ID", "Description_trans")],
                         by = c("node" = "ID"))

      text@merged_pathway_reactome$result_with_module <-
        result_with_module_reactome
    } else{
      module_result_reactome <-
        NULL
    }

    ####merged module
    if (length(text@merged_module) > 0) {
      functional_module_result <-
        text@merged_module$functional_module_result

      result_with_module <-
        text@merged_module$result_with_module

      if (!is.null(module_result_go)) {
        module_result_go <-
          module_result_go[, c(
            "module",
            "module_annotation",
            "module_annotation_trans",
            "Description_trans"
          )]
      }

      if (!is.null(module_result_kegg)) {
        module_result_kegg <-
          module_result_kegg[, c(
            "module",
            "module_annotation",
            "module_annotation_trans",
            "Description_trans"
          )]
      }

      if (!is.null(module_result_reactome)) {
        module_result_reactome <-
          module_result_reactome[, c(
            "module",
            "module_annotation",
            "module_annotation_trans",
            "Description_trans"
          )]
      }

      module_result <-
        rbind(module_result_reactome,
              module_result_go,
              module_result_kegg)

      result <-
        rbind(result_reactome,
              result_go,
              result_kegg)

      functional_module_result <-
        functional_module_result %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(module_result[, c("module_annotation",
                                           "module_annotation_trans")],
                         by = c("module_annotation"))

      functional_module_result$Description_trans <-
        lapply(functional_module_result$Description,
               function(x) {
                 result$Description_trans[match(stringr::str_split(x, ";")[[1]],
                                                result$Description)] %>%
                   paste(collapse = ";")
               }) %>%
        unlist()

      text@merged_module$functional_module_result <-
        functional_module_result

      result_with_module <-
        result_with_module %>%
        dplyr::select(-dplyr::contains("_trans")) %>%
        dplyr::left_join(
          result[, c("Description", "Description_trans")] %>%
            dplyr::rename(
              module_annotation = Description,
              module_annotation_trans = Description_trans
            ),
          by = "module_annotation"
        )

      result_with_module$Description_trans <-
        lapply(result_with_module$Description,
               function(x) {
                 result$Description_trans[match(stringr::str_split(x, ";")[[1]],
                                                result$Description)] %>%
                   paste(collapse = ";")
               }) %>%
        unlist()

      text@merged_module$result_with_module <-
        result_with_module
    }

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "translate_language()",
      parameter = list(engine = engine,
                       to = to),
      time = Sys.time()
    )

    process_info <-
      slot(text, "process_info")

    process_info$translate_language <-
      parameter

    slot(text, "process_info") <-
      process_info

    message("Done")
    return(text)
  }



#' Internal Language Translation Function
#'
#' An internal helper function for `translate_language` to handle the actual translation logic.
#' It communicates with the specified translation engine and performs the translation.
#'
#' @param text Text to be translated, defaults to a specific string if not provided.
#' @param engine Translation engine to use, defaults to 'gemini' or 'chatgpt'.
#' @param to Target language for translation, with multiple language options.
#' @return Depending on the input, either a translated string or a data frame mapping original
#'   texts to their translations.
#' @noRd
translate_language_internal <-
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

    engine <-
      match.arg(engine)

    to <-
      match.arg(to)

    if (length(text) == 1) {
      prompt <-
        paste0(
          "I will give you a text.",
          " Please translate it to ",
          to,
          ". Please return only the translated text,
               do not include any additional text or explanation.\n",
          "The text is below: \n",
          text
        )

      if (engine == "chatgpt") {
        prompt <- paste0("Translate to ", to, ":\n", text)
        return(request_chatgpt_response(prompt = prompt))
      }

      if (engine == "gemini") {
        prompt <- paste0("Translate to ", to, ":\n", text)
        return(request_gemini_response(prompt = prompt))
      }
    } else{
      prompt <-
        paste0(
          "I have a string formatted with segments separated by '{}', like 'apple{}title{car}'.
          Please follow these instructions:",
          paste0("1. Translate each segment into ", to, "\n"),
          "2. Ensure that the length of each original segment is exactly the same as its translated segment\n",
          "3. Format the output as 'original=translation'.\n",
          "4. If a segment cannot be translated or you're unsure about it,
          return it as 'original=original',
          with the length matched as in point 2.",
          "Keep the '{}' separators between each segment in the output.\n",
          "5. Please return only the translated texts,
                     do not include any additional text or explanation.\n",
          "Example: For an input 'apple{}bana',
          if 'bana' can't be translated or is unknown,
          and ensuring length match,
          the output should be 'apple=苹果{}bana=bana'.\n",
          "The string is below:\n",
          paste(text, collapse = "{}")
        )


      if (engine == "chatgpt") {
        result <-
          tryCatch(
            request_chatgpt_response(prompt = prompt),
            error = function(e) {
              return(NULL)
            }
          )
      }
      if (engine == "gemini") {
        result <-
          tryCatch(
            request_gemini_response(prompt = prompt),
            error = function(e) {
              return(NULL)
            }
          )
      }

      ###if there is error of translation, return original text
      if (is.null(result)) {
        return(text)
      }

      result <-
        stringr::str_split(result, "\\{\\}")[[1]]
      result <- result[result != ""]

      result <-
        as.data.frame(do.call(rbind, stringr::str_split(result, "=")))

      colnames(result) <- c("original", "translation")


      result <-
        tryCatch(
          data.frame(original = text) %>%
            dplyr::left_join(result, by = "original") %>%
            pull(translation),
          error = function(e) {
            return(NULL)
          }
        )

      if (is.null(result)) {
        return(text)
      }

      result[is.na(result)] <-
        text[is.na(result)]

      return(result)
    }
  }
