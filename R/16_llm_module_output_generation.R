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
#' @param phenotype Character string. Phenotype or disease to focus on. Default is NULL.
#' @param model A string specifying the GPT model to use. Default is `"gpt-4o-mini-2024-07-18"`.
#' @param api_key A string containing the API key required to access the AI API.
#' @param output_prompt Logical. Whether to output prompt in final annotation result. Default is TRUE.
#' @param api_provider A string indicating the API provider, either `"openai"`, `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param thinkingBudget An integer for the "thinking budget" parameter specific to the Gemini API (default is `0`).
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
                                     phenotype = NULL,
                                     model = "gpt-4o-mini-2024-07-18",
                                     api_key,
                                     output_prompt = TRUE,
                                     api_provider = "openai",
                                     thinkingBudget = 0) {
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

  if (is.null(phenotype)) {
    prompt_path <- system.file("prompts_template","16_llm_prompt.md", package = "mapa")
    prompt_text <- readLines(prompt_path, warn = FALSE)
  } else {
    prompt_path <- system.file("prompts_template","16_llm_prompt_with_phenotype.md", package = "mapa")
    prompt_text <- readLines(prompt_path, warn = FALSE)
    prompt_text <- gsub("\\{phenotype\\}", phenotype, prompt_text)
  }

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

  if (nchar(combined_texts) == 0) {
    prompt_text <- gsub("Below are related articles: \\{combined_texts\\}", combined_texts, prompt_text)
  } else {
    prompt_text <- gsub("\\{combined_texts\\}", combined_texts, prompt_text)
  }

  messages <- list(
    list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist."),
    list(
      role = "user",
      content = prompt_text
    )
  )


  gpt_response <- gpt_api_call(messages, api_key, model = model, api_provider = api_provider, thinkingBudget = thinkingBudget)

  if (!check_json_format_output_generation(gpt_response)) {
    gpt_response <- gpt_api_call(messages, api_key, model = model, api_provider = api_provider, thinkingBudget = thinkingBudget)
    if (!check_json_format_output_generation(gpt_response)) {
      prompt <- modify_prompt_for_format_output_generation(gpt_response)
      gpt_response <- gpt_api_call(prompt, api_key, model = model, api_provider = api_provider, thinkingBudget = thinkingBudget)
      if (!check_json_format_output_generation(gpt_response)) {
        print(gpt_response)
        gpt_response <- '{"module_name": "Default Module Name", "summary": "Unable to process the request."}'
      }
    }
  }

  result <- jsonlite::fromJSON(gpt_response)

  if (is.null(phenotype)) {
    if (output_prompt) {
      return(list(
        module_name = result$module_name,
        summary = result$summary,
        confidence_score = result$confidence_score,
        prompt = messages
      ))
    } else {
      return(list(
        module_name = result$module_name,
        summary = result$summary,
        confidence_score = result$confidence_score
      ))
    }
  } else {
    if (output_prompt) {
      return(list(
        module_name = result$module_name,
        summary = result$summary,
        phenotype_analysis = result$phenotype_analysis,
        confidence_score = result$confidence_score,
        prompt = messages
      ))
    } else {
      return(list(
        module_name = result$module_name,
        summary = result$summary,
        phenotype_analysis = result$phenotype_analysis,
        confidence_score = result$confidence_score
      ))
    }
  }
}

#' Check JSON Format Output Generation
#'
#' Validates if a response string contains valid JSON with required 'module_name' and 'summary' fields.
#'
#' @param response Character string containing JSON response to validate
#' @return Logical. TRUE if valid JSON with both required fields, FALSE otherwise
#' @export

check_json_format_output_generation <- function(response) {
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

#' Modify Prompt for Format Output Generation
#'
#' Creates structured messages to convert GPT responses into required JSON format.
#'
#' @param gpt_response Character string containing original GPT response to convert
#' @return List of message objects formatted for chat-based AI models
#' @export

modify_prompt_for_format_output_generation <- function(gpt_response) {
  messages <- list(
    list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist."),
    list(role = "user", content = paste0(
      "Please convert the following response to the required JSON format:\n\n",
      gpt_response,
      "\n\nReturn a valid JSON structure like this:\n",
      "{\n",
      "  \"module_name\": \"<biological module name>\",\n",
      "  \"summary\": \"<summary of the current research>\"\n",
      "}"
    ))
  )
  return(messages)
}

#' Module Name Generation
#'
#' This internal function processes a list of modules to generate final biological module names
#' and summaries for all modules by calling \code{single_module_generation}.
#'
#' @param paper_result A list containing module information, including related papers and PubMed results.
#' @param phenotype Character string. Phenotype or disease to focus on. Default is NULL.
#' @param model A string specifying the GPT model to use. Default is `"gpt-4o-mini-2024-07-18"`.
#' @param api_key A string containing the API key required to access the AI API.
#' @param output_prompt Logical. Whether to output prompt in final annotation result. Default is TRUE.
#' @param api_provider A string indicating the API provider, either `"openai"`, `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param thinkingBudget An integer for the "thinking budget" parameter specific to the Gemini API (default is `0`).
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
module_name_generation <- function(paper_result,
                                   phenotype = NULL,
                                   model = "gpt-4o-mini-2024-07-18",
                                   api_key,
                                   output_prompt = TRUE,
                                   api_provider = "openai",
                                   thinkingBudget = 0) {
  for (module_index in seq_along(paper_result)) {
    module_list <- paper_result[[module_index]]

    module_related_paper <- module_list[["related_paper"]]
    module_info <- module_list[["pubmed_result"]]

    final_result <- single_module_generation(module_related_paper = module_related_paper,
                                             module_info = module_info,
                                             phenotype = phenotype,
                                             model = model,
                                             api_key = api_key,
                                             output_prompt = output_prompt, api_provider = api_provider, thinkingBudget = thinkingBudget)

    # 将结果直接存入 paper_result
    paper_result[[module_index]][["generated_name"]] <- final_result
  }

  return(paper_result)
}

