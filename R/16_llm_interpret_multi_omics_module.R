# load("demo_data/demo_multi-omics/multi_omics_modules.rda")
# {
#   object = multi_omics_modules
#   module_content_number_cutoff = 1
#   llm_model = "gpt-4o-mini-2024-07-18"
#   embedding_model = "text-embedding-3-small"
#   # api_key
#   embedding_output_dir = "demo_data/demo_multi-omics/embedding_output/"
#   local_corpus_dir = NULL
#   phenotype = NULL
#   chunk_size = 5
#   years = 5
#   retmax = 10
#   similarity_filter_num = 20
#   GPT_filter_num = 5
#   orgdb = org.Mmu.eg.db
#   output_prompt = TRUE
#   api_provider = "openai"
#   thinkingBudget = 0
#   thread = 10
# }
#
# object$functional_module_result <- object$functional_module_result |> dplyr::filter((multi_omics_num == 3 | multi_omics_num == 2) & module_content_number < 50)
#
# llm_interpreted_object <- llm_interpret_multi_omics_module(object = object,
#                                                            module_content_number_cutoff = 1,
#                                                            llm_model = "gpt-4o-mini-2024-07-18",
#                                                            embedding_model = "text-embedding-3-small",
#                                                            api_key = api_key,
#                                                            embedding_output_dir = "demo_data/demo_multi-omics/embedding_output/",
#                                                            local_corpus_dir = NULL,
#                                                            phenotype = NULL,
#                                                            chunk_size = 5,
#                                                            years = 5,
#                                                            retmax = 20,
#                                                            similarity_filter_num = 20,
#                                                            GPT_filter_num = 5,
#                                                            orgdb = org.Mmu.eg.db,
#                                                            output_prompt = FALSE,
#                                                            api_provider = "openai",
#                                                            thinkingBudget = 0,
#                                                            thread = 10)


#' Interpret Multi-Omics Functional Modules using LLM
#'
#' Internal implementation for the multi-omics branch of [llm_interpret_module()].
#' Retrieves relevant literature via a RAG strategy and uses an LLM to generate
#' biological names and summaries for each functional module.
#'
#' @param object A \code{list} produced by [merge_multi_omics_nodes()], containing
#'   \code{graph_data}, \code{functional_module_result}, and \code{result_with_module}.
#' @param module_content_number_cutoff Integer. Only modules with content number
#'   greater than this value are processed. Default \code{1}.
#' @param llm_model Character. LLM model identifier. Default
#'   \code{"gpt-4o-mini-2024-07-18"}.
#' @param embedding_model Character. Embedding model identifier. Default
#'   \code{"text-embedding-3-small"}.
#' @param api_key Character. API key for the chosen provider.
#' @param embedding_output_dir Character. Directory for embedding cache files.
#' @param local_corpus_dir Character or \code{NULL}. Path to user-supplied local
#'   documents. Default \code{NULL}.
#' @param phenotype Character or \code{NULL}. Phenotype/disease to focus on.
#'   Default \code{NULL}.
#' @param chunk_size Integer. PubMed query chunk size. Default \code{5}.
#' @param years Integer. Years to look back in PubMed. Default \code{5}.
#' @param retmax Integer. Max PubMed records per query. Default \code{10}.
#' @param similarity_filter_num Integer. Top-N documents after embedding
#'   similarity filtering. Default \code{20}.
#' @param GPT_filter_num Integer. Top-N documents after LLM re-ranking.
#'   Default \code{5}.
#' @param orgdb Organism database object for gene annotation. Default
#'   \code{org.Hs.eg.db}.
#' @param output_prompt Logical. Include the LLM prompt in results. Default
#'   \code{TRUE}.
#' @param api_provider Character. API provider: \code{"openai"},
#'   \code{"gemini"}, or \code{"siliconflow"}. Default \code{"openai"}.
#' @param thinkingBudget Integer. Gemini thinking budget. Default \code{0}.
#' @param thread Integer. Number of parallel threads. Default \code{10}.
#'
#' @return The updated \code{object} list with an added \code{llm_module_interpretation}
#'   element and an updated \code{functional_module_result} table containing a
#'   \code{llm_module_name} column.
#'
#' @noRd
llm_interpret_multi_omics_module <- function(object,
                                             module_content_number_cutoff = 1,
                                             llm_model = "gpt-4o-mini-2024-07-18",
                                             embedding_model = "text-embedding-3-small",
                                             api_key,
                                             embedding_output_dir,
                                             local_corpus_dir = NULL,
                                             phenotype = NULL,
                                             chunk_size = 5,
                                             years = 5,
                                             retmax = 10,
                                             similarity_filter_num = 20,
                                             GPT_filter_num = 5,
                                             orgdb = org.Hs.eg.db,
                                             output_prompt = TRUE,
                                             api_provider = "openai",
                                             thinkingBudget = 0,
                                             thread = 10) {
  # convert character orgdb into orgdb object
  if (is(orgdb, "character")) {
    # Check if package is installed
    if (!requireNamespace(orgdb, quietly = TRUE)) {
      stop("Package '", orgdb, "' is not installed. ",
           "Please install it using:\n",
           "  BiocManager::install('", orgdb, "')")
    }
    # Load the orgdb object from the package
    orgdb <- utils::getFromNamespace(orgdb, orgdb)
  }

  # 1. Collect functional module result ====
  if (module_content_number_cutoff >= max(object$functional_module_result$module_content_number)) {
    stop("module_content_number_cutoff should be smaller than the maximum of all module content numbers in your functional module result.")
  }

  functional_module_result <-
    object$functional_module_result |>
    dplyr::filter(module_content_number > module_content_number_cutoff)

  # # Check if local corpus parameters are provided when local_corpus is TRUE
  # if (local_corpus) {
  #   if (is.null(local_corpus_dir)) {
  #     stop("When local_corpus is TRUE, local_corpus_dir must be provided.")
  #   }
  #   save_dir_local_corpus_embed = "local"
  # }
  if (!is.null(local_corpus_dir)) {
    local_corpus <- TRUE
    save_dir_local_corpus_embed <- file.path(local_corpus_dir, "local_corpus_embedding_output")
  } else {
    local_corpus <- FALSE
  }

  # 2. Create vector database for local corpus uploaded by user ====
  clear_output_dir(output_dir = embedding_output_dir) # Clear output directory

  if (!is.null(local_corpus_dir)) {
    embedding_local_corpus(embedding_model = embedding_model,
                           api_provider =  api_provider,
                           api_key = api_key,
                           local_corpus_dir = local_corpus_dir,
                           embedding_output_dir = embedding_output_dir,
                           save_dir_local_corpus_embed = save_dir_local_corpus_embed,
                           thread = thread)
  }

  # 3. Extract pathway description, molecule name, and pmid of references
  processed_data <- preprocess_multi_omics_module(df = functional_module_result,
                                                  orgdb = orgdb)

  # 4. Retrieve PubMedIDs of related articles
  pubmed_result <- pubmed_search(processed_data = processed_data,
                                 phenotype = phenotype,
                                 chunk_size = chunk_size,
                                 years = years,
                                 retmax = retmax,
                                 thread = thread)

  # Reference paper PMID save to PubmedIDs
  for (module in names(pubmed_result)) {
    # Extract PathwayReferencePMID
    pmids <- pubmed_result[[module]][["PathwayReferencePMID"]]

    # Ensure pmids exists, and remove empty strings
    if (!is.null(pmids)) {
      pmids <- pmids[pmids != "" & !is.na(pmids)]

      # If PubmedIDs already exists, then merge, otherwise create
      if (!is.null(pubmed_result[[module]][["PubmedIDs"]])) {
        pubmed_result[[module]][["PubmedIDs"]] <- unique(c(pubmed_result[[module]][["PubmedIDs"]], pmids))
      } else {
        pubmed_result[[module]][["PubmedIDs"]] <- pmids
      }
    }
  }

  # 5. Save search results for each module as CSV.GZ files using embedding database
  embedding_pubmed_search(pubmed_result = pubmed_result,
                          embedding_model = embedding_model,
                          api_provider =  api_provider,
                          api_key = api_key,
                          embedding_output_dir = embedding_output_dir,
                          thread = thread)

  # 6. Retrieve and rank related papers using RAG strategy
  related_paper <- retrieve_strategy(pubmed_result = pubmed_result,
                                     model = llm_model,
                                     embedding_model = embedding_model,
                                     api_key = api_key,
                                     api_provider = api_provider,
                                     similarity_filter_num = similarity_filter_num,
                                     GPT_filter_num = GPT_filter_num,
                                     local_corpus = local_corpus,
                                     embedding_output_dir = embedding_output_dir,
                                     save_dir_local_corpus_embed = save_dir_local_corpus_embed,
                                     thread = thread,
                                     multi_omics = TRUE)

  paper_result <- Map(function(x, y) {
    # Customize operations for each pair of related data
    return(list(related_paper = x, pubmed_result = y))
  }, related_paper, pubmed_result)

  # 7. Generate module names and study summaries
  message("Start to generate name and summary for functional modules ...")
  final_result <- module_name_generation(paper_result = paper_result,
                                         phenotype = phenotype,
                                         model = llm_model,
                                         api_key = api_key,
                                         output_prompt = output_prompt,
                                         api_provider = api_provider,
                                         thinkingBudget = thinkingBudget,
                                         multi_omics = TRUE)

  # 8. Store the final result in the object's llm_module_interpretation slot
  object[["llm_module_interpretation"]] <- final_result

  # 9. Update graph_data and functional_module result according to llm interpretation
  llm_module_name_df <- data.frame()

  for (i in 1:length(object$llm_module_interpretation)){
    module <- names(object$llm_module_interpretation[i])
    llm_module_name <- object$llm_module_interpretation[[i]]$generated_name$module_name
    llm_module_name_df[i, 1] <- module
    llm_module_name_df[i, 2] <- llm_module_name
  }
  colnames(llm_module_name_df) <- c("module", "llm_module_name")
  object$functional_module_result <-
    object$functional_module_result |>
    dplyr::left_join(llm_module_name_df, by = "module")

  message("Done")

  # 11. Return the updated object
  return(object)
}
