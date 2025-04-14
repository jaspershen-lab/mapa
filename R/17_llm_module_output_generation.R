# This file is for final name and current study summary
# library(roxygen2)


#' Single Module Generation
#'
#' This internal function generates a biological module name and a research summary
#' for a specific module based on provided pathway, gene, and related article information.
#' It interacts with an AI API to perform text analysis and ensure the response is in JSON format.
#'
#' @param module_related_paper A list of related papers with titles and cleaned texts for the module.
#' @param module_info A list containing pathway names and gene symbols relevant to the module.
#' @param api_key A string containing the API key required to access the AI API.
#'
#' @return A list containing two elements: \code{module_name} (the generated biological module name)
#' and \code{summary} (the research summary).
#'
#' @importFrom jsonlite fromJSON
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal
# Generate a biological module name and research summary using GPT
single_module_generation <- function(module_related_paper,
                                     module_info,
                                     api_key) {
  pathway_info <- paste(module_info[["PathwayNames"]], "(", module_info[["PathwayDescription"]], ")", collapse = "; ")
  if ("GeneNames_vec" %in% names(module_info)) {
    gene_names <- paste(module_info[["GeneNames_vec"]], collapse = ", ")
  } else if ("MetNames_vec" %in% names(module_info)) {
    met_names <- paste(module_info[["MetNames_vec"]], collapse = ", ")
  }

  titles <- sapply(module_related_paper, function(x) x[["title"]])
  cleaned_texts <- sapply(module_related_paper, function(x) x[["cleaned_text"]])

  title_text_pairs <- mapply(function(title, cleaned_text) {
    paste("Title: ", title, "\nText: ", cleaned_text, sep = "")
  }, titles, cleaned_texts)
  combined_texts <- paste(title_text_pairs, collapse = "\n\n")

  prompt_text <- readLines("R/17_llm_prompt.md", warn = FALSE)
  prompt_text <- paste(prompt_text, collapse = "\n")
  prompt_text <- gsub("\\{pathway_info\\}", pathway_info, prompt_text)

  if ("GeneNames_vec" %in% names(module_info)) {
    prompt_text <- gsub("\\{query_molecule_names\\}", gene_names, prompt_text)
    prompt_text <- gsub("\\{query_molecules\\}", "genes", prompt_text)
    prompt_text <- gsub("\\{query_product\\}", "protein", prompt_text)
  } else if ("MetNames_vec" %in% names(module_info)) {
    prompt_text <- gsub("\\{query_molecule_names\\}", met_names, prompt_text)
    prompt_text <- gsub("\\{query_molecules\\}", "compounds", prompt_text)
    prompt_text <- gsub("\\{query_product\\}", "compound", prompt_text)
  }

  prompt_text <- gsub("\\{combined_texts\\}", combined_texts, prompt_text)

  messages <- list(
    list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist."),
    list(
      role = "user",
      content = prompt_text
    )
  )


  gpt_response <- gpt_api_call(messages, api_key)

  if (!check_json_format(gpt_response)) {
    gpt_response <- gpt_api_call(messages, api_key)
    if (!check_json_format(gpt_response)) {
      prompt <- modify_prompt_for_format(messages)
      gpt_response <- gpt_api_call(prompt, api_key)
      if (!check_json_format(gpt_response)) {
        print(gpt_response)
        gpt_response <- '{"module_name": "Default Module Name", "summary": "Unable to process the request."}'
      }
    }
  }

  result <- jsonlite::fromJSON(gpt_response)

  return(list(
    module_name = result$module_name,
    summary = result$summary,
    confidence_score = result$confidence_score
  ))
}

check_json_format <- function(response) {
  tryCatch({
    result <- jsonlite::fromJSON(response)
    if (!is.null(result$module_name) && !is.null(result$summary)) {
      return(TRUE)
    }
  }, error = function(e) {
    return(FALSE)
  })
  return(FALSE)
}

modify_prompt_for_format <- function(messages) {
  # 找到用户输入消息
  user_message_index <- which(sapply(messages, function(msg) msg$role == "user"))

  if (length(user_message_index) > 0) {
    # 修改用户输入消息的内容
    user_message_index <- user_message_index[1]  # 假定只有一个用户输入消息
    original_content <- messages[[user_message_index]]$content

    # 添加格式说明到用户输入消息
    messages[[user_message_index]]$content <- paste0(
      original_content,
      "\n\nIf your response does not strictly follow the JSON format, please fix the format and make sure to return a valid JSON structure like this:\n",
      "{\n",
      "  \"module_name\": \"<biological module name>\",\n",
      "  \"summary\": \"<summary of the current research>\"\n",
      "}"
    )
  } else {
    warning("No user message found in the provided messages.")
  }

  return(messages)
}

#' Module Name Generation
#'
#' This internal function processes a list of modules to generate final biological module names
#' and summaries for all modules by calling \code{single_module_generation}.
#'
#' @param paper_result A list containing module information, including related papers and PubMed results.
#' @param api_key A string containing the API key required to access the AI API.
#'
#' @return A list of results for each module, where each element is a list containing
#' \code{module_name} and \code{summary}.
#'
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{single_module_generation}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @keywords internal

# Generate final module names and summaries for all modules
module_name_generation <- function(paper_result, api_key) {
  for (module_index in seq_along(paper_result)) {
    module_list <- paper_result[[module_index]]

    module_related_paper <- module_list[["related_paper"]]
    module_info <- module_list[["pubmed_result"]]

    final_result <- single_module_generation(module_related_paper = module_related_paper,
                                             module_info = module_info,
                                             api_key = api_key)

    # 将结果直接存入 paper_result
    paper_result[[module_index]][["generated_name"]] <- final_result
  }

  return(paper_result)
}

