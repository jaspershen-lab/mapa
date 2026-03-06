#' Single Module Generation
#'
#' This internal function generates a biological module name and a research summary
#' for a specific module based on provided pathway, gene, and related article information.
#' It interacts with an AI API to perform text analysis and ensure the response is in JSON format.
#'
#' @param module_related_paper A list of related papers with titles and cleaned texts for the module.
#' @param module_info A list containing pathway names and gene symbols relevant to the module.
#'   When \code{multi_omics = TRUE}, \code{module_info} must contain the fields:
#'   \code{GeneIDs}, \code{GeneNames_vec}, \code{MetIDs}, \code{MetNames_vec},
#'   \code{PathwayNames}, \code{PathwayDescription}, and \code{PathwayReferencePMID}.
#' @param phenotype Character string. Phenotype or disease to focus on. Default is NULL.
#' @param model A string specifying the GPT model to use. Default is `"gpt-4o-mini-2024-07-18"`.
#' @param api_key A string containing the API key required to access the AI API.
#' @param output_prompt Logical. Whether to output prompt in final annotation result. Default is TRUE.
#' @param api_provider A string indicating the API provider, either `"openai"`, `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param thinkingBudget An integer for the "thinking budget" parameter specific to the Gemini API (default is `0`).
#' @param multi_omics Logical. If TRUE, use multi-omics prompt templates that integrate genes,
#'   metabolites, and pathways together. Default is FALSE.
#'
#' @return A list containing: \code{module_name}, \code{summary}, \code{confidence_score},
#'   optionally \code{phenotype_analysis} (when \code{phenotype} is provided),
#'   and optionally \code{prompt} (when \code{output_prompt} is TRUE).
#'
#' @importFrom jsonlite fromJSON
#'
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @noRd
single_module_generation <- function(module_related_paper,
                                     module_info,
                                     phenotype = NULL,
                                     model = "gpt-4o-mini-2024-07-18",
                                     api_key,
                                     output_prompt = TRUE,
                                     api_provider = "openai",
                                     thinkingBudget = 0,
                                     multi_omics = FALSE) {

  # Build combined literature text
  titles <- sapply(module_related_paper, function(x) x[["title"]])
  cleaned_texts <- sapply(module_related_paper, function(x) x[["cleaned_text"]])

  title_text_pairs <- mapply(function(title, cleaned_text) {
    paste("Title: ", title, "\nText: ", cleaned_text, sep = "")
  }, titles, cleaned_texts)
  combined_texts <- paste(title_text_pairs, collapse = "\n\n")

  # multi-omics mode
  if (multi_omics) {

    # Helper: collapse a field to "NULL" if missing/empty/NA
    is_empty_field <- function(x) {
      is.null(x) || length(x) == 0 || all(is.na(x)) ||
        all(nchar(trimws(as.character(x))) == 0)
    }

    # Gene names: comma-separated, or "NULL" if missing
    gene_names_str <- if (is_empty_field(module_info[["GeneNames_vec"]])) {
      "NULL"
    } else {
      paste(module_info[["GeneNames_vec"]], collapse = ", ")
    }

    # Metabolite names: comma-separated, or "NULL" if missing
    met_names_str <- if (is_empty_field(module_info[["MetNames_vec"]])) {
      "NULL"
    } else {
      paste(module_info[["MetNames_vec"]], collapse = ", ")
    }

    # Pathway text: one entry per line pair, or a "not provided" note
    pathway_names <- module_info[["PathwayNames"]]
    pathway_desc  <- module_info[["PathwayDescription"]]

    combined_pathway_text <- if (is_empty_field(pathway_names) || is_empty_field(pathway_desc)) {
      "No pathway enrichment terms were provided for this module."
    } else {
      paste(
        mapply(function(nm, desc) paste0(nm, " (", desc, ")"),
               pathway_names, pathway_desc),
        collapse = "\n\n"
      )
    }

    # Select prompt template
    if (is.null(phenotype)) {
      prompt_path <- system.file("prompts_template", "17_llm_prompt_multiomics_module.md",
                                 package = "mapa")
    } else {
      prompt_path <- system.file("prompts_template",
                                 "17_llm_prompt_multiomics_module_with_phenotype.md",
                                 package = "mapa")
    }

    prompt_text <- readLines(prompt_path, warn = FALSE)
    prompt_text <- paste(prompt_text, collapse = "\n")

    # Fill in placeholders
    prompt_text <- gsub("\\{GeneNames_vec\\}", gene_names_str, prompt_text)
    prompt_text <- gsub("\\{MetNames_vec\\}", met_names_str, prompt_text)
    prompt_text <- gsub("\\{combined_pathway_text\\}", combined_pathway_text, prompt_text)

    if (!is.null(phenotype)) {
      prompt_text <- gsub("\\{phenotype\\}", phenotype, prompt_text)
    }

    if (nchar(combined_texts) == 0) {
      prompt_text <- gsub("\\{combined_texts\\}", "No related articles provided.", prompt_text)
    } else {
      prompt_text <- gsub("\\{combined_texts\\}", combined_texts, prompt_text)
    }

    # original single-omics mode (unchanged)
  } else {

    pathway_info <- paste(module_info[["PathwayNames"]], "(", module_info[["PathwayDescription"]], ")", collapse = "; ")
    if ("GeneNames_vec" %in% names(module_info)) {
      gene_names <- paste(module_info[["GeneNames_vec"]], collapse = ", ")
    } else if ("MetNames_vec" %in% names(module_info)) {
      met_names <- paste(module_info[["MetNames_vec"]], collapse = ", ")
    }

    if (is.null(phenotype)) {
      prompt_path <- system.file("prompts_template", "16_llm_prompt.md", package = "mapa")
      prompt_text <- readLines(prompt_path, warn = FALSE)
    } else {
      prompt_path <- system.file("prompts_template", "16_llm_prompt_with_phenotype.md", package = "mapa")
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
  }

  # Build messages and call API
  messages <- list(
    list(role = "system", content = "You are an efficient and insightful assistant to a molecular biologist."),
    list(role = "user",   content = prompt_text)
  )

  gpt_response <- gpt_api_call(messages, api_key, model = model,
                               api_provider = api_provider,
                               thinkingBudget = thinkingBudget)

  if (!check_json_format_output_generation(gpt_response)) {
    gpt_response <- gpt_api_call(messages, api_key, model = model,
                                 api_provider = api_provider,
                                 thinkingBudget = thinkingBudget)
    if (!check_json_format_output_generation(gpt_response)) {
      prompt <- modify_prompt_for_format_output_generation(gpt_response)
      gpt_response <- gpt_api_call(prompt, api_key, model = model,
                                   api_provider = api_provider,
                                   thinkingBudget = thinkingBudget)
      if (!check_json_format_output_generation(gpt_response)) {
        print(gpt_response)
        gpt_response <- '{"module_name": "Default Module Name", "summary": "Unable to process the request."}'
      }
    }
  }

  result <- jsonlite::fromJSON(gpt_response)

  # Assemble return value
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
#' @param multi_omics Logical. If TRUE, use multi-omics prompt templates. Default is FALSE.
#' @return A list of results for each module, where each element is a list containing
#' \code{module_name} and \code{summary}.
#'
#' @importFrom jsonlite fromJSON
#'
#' @seealso \code{\link{single_module_generation}}
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#'
#' @noRd
module_name_generation <- function(paper_result,
                                   phenotype = NULL,
                                   model = "gpt-4o-mini-2024-07-18",
                                   api_key,
                                   output_prompt = TRUE,
                                   api_provider = "openai",
                                   thinkingBudget = 0,
                                   multi_omics = FALSE) {
  for (module_index in seq_along(paper_result)) {
    module_list <- paper_result[[module_index]]

    module_related_paper <- module_list[["related_paper"]]
    module_info <- module_list[["pubmed_result"]]

    final_result <- single_module_generation(module_related_paper = module_related_paper,
                                             module_info = module_info,
                                             phenotype = phenotype,
                                             model = model,
                                             api_key = api_key,
                                             output_prompt = output_prompt,
                                             api_provider = api_provider,
                                             thinkingBudget = thinkingBudget,
                                             multi_omics = multi_omics)

    # 将结果直接存入 paper_result
    paper_result[[module_index]][["generated_name"]] <- final_result
  }

  return(paper_result)
}

