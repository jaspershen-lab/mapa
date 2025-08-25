# source("R/16_llm_module_pathway_infor.R")
# source("R/16_llm_module_online_retrieval.R")
# source("R/16_llm_module_embedding_database.R")
# source("R/16_llm_module_RAG_strategy.R")
# source("R/16_llm_module_utils.R")
# source('R/16_llm_module_output_generation.R')

# Metabolite
# llm_interpreted_enriched_functional_module <-
#   llm_interpret_module(
#     object = object,
#     llm_model = "gpt-4o-mini-2024-07-18",
#     embedding_model = "text-embedding-3-small",
#     api_key = api_key,
#     embedding_output_dir = "demo_data/debug_yijiang/embedding_ouput/"
# )
# load("demo_data/test.rda")
# module_content_number_cutoff = 0
# api_provider = "openai"
# llm_model = "gpt-4o-mini-2024-07-18"
# embedding_model = "text-embedding-3-small"
# orgdb = org.Hs.eg.db
# embedding_output_dir = "demo_data/debug_yijiang/embedding_ouput/"
# test_llm_met <- llm_interpret_module(object = object,
#                                      module_content_number_cutoff = 0,
#                                      api_key = api_key,
#                                      embedding_output_dir = "demo_data/debug_yijiang/embedding_ouput/")

# Gene
# load("demo_data/updated_object_results_for_genes_ora/biotext_sim_result/biotext_functional_modules.rda")
# llm_interpreted_functional_module <-
#   llm_interpret_module(
#     object = biotext_functional_modules,
#     module_content_number_cutoff = 1,
#     api_provider = "openai",
#     api_key = api_key,
#     llm_model = "gpt-4o-mini-2024-07-18",
#     embedding_model = "text-embedding-3-small",
#     orgdb = org.Hs.eg.db,
#     embedding_output_dir = "demo_data/updated_object_results_for_genes_ora/biotext_sim_result/embedding_output/"
#   )

# ah <- AnnotationHub::AnnotationHub()
# mf.orgdb <- ah[["AH119900"]] #| Taxonomy ID: 9541
#
# llm_interpreted_functional_module <-
#   llm_interpret_module(
#     object = functional_module_res,
#     module_content_number_cutoff = 1,
#     api_provider = "openai",
#     api_key = api_key,
#     llm_model = "gpt-4o-mini-2024-07-18",
#     embedding_model = "text-embedding-3-small",
#     orgdb = mf.orgdb,
#     embedding_output_dir = "demo_data/embedding_output/"
#   )
#
# llm_interpreted_functional_module <-
#   llm_interpret_module(
#     object = biotext_functional_modules,
#     module_content_number_cutoff = 2,
#     api_provider = "gemini",
#     api_key = api_key,
#     llm_model = "gemini-1.5-flash",
#     embedding_model = "text-embedding-004",
#     orgdb = org.Hs.eg.db,
#     embedding_output_dir = "demo_data/updated_object_results_for_genes_ora/biotext_sim_result/embedding_output/"
#   )
# save(llm_interpreted_functional_module, file = "demo_data/updated_object_results_for_genes_ora/biotext_sim_result/llm_interpreted_functional_module_cutoff_2.rda")

#' Interpret Functional Modules using LLM with RAG Strategy
#'
#' @description
#' This function processes functional module results by retrieving relevant scientific
#' literature using a Retrieval-Augmented Generation (RAG) strategy. The process includes
#' pathway description extraction, PubMed searching, embedding of search results, and
#' retrieving related papers with ranking. Finally, functional module information and
#' abstracts/titles of retrieved papers are input into an LLM to generate meaningful
#' names and summaries for each functional module.
#'
#' @param object A `functional_module` class object that has been processed with
#'   `get_functional_modules()` function.
#' @param module_content_number_cutoff Integer. Only modules with content number
#'   greater than this value will be processed. Must be smaller than the maximum
#'   module content number in the results. Default is `1`.
#' @param llm_model Character string. The LLM model to use for text generation.
#'   Default is `"gpt-4o-mini-2024-07-18"`.
#' @param embedding_model Character string. The embedding model to use for text
#'   embeddings. Default is `"text-embedding-3-small"`.
#' @param api_key Character string. API key for OpenAI or other embedding service.
#'   Currently, only OpenAI API keys are supported.
#' @param api_provider Character string. The API provider to use, one of `"openai"`,
#'   `"gemini"`, or `"siliconflow"` (default is `"openai"`).
#' @param embedding_output_dir Character string. Directory path where embedding
#'   results will be saved. This directory will be cleared before processing.
#' @param local_corpus_dir Character string. Optional directory path containing
#'   local files provided by users for additional context. If provided, these files
#'   will be embedded and used in the RAG process. Default is `NULL`.
#' @param phenotype Character string. Optional phenotype or disease name to focus
#'   the literature search and interpretation on. Default is `NULL`.
#' @param chunk_size Integer. Chunk size for processing data in batches. Default is `5`.
#' @param years Integer. Number of recent years to include in PubMed search.
#'   Default is `5`.
#' @param retmax Integer. Maximum number of records to return from each PubMed
#'   search query. Default is `10`.
#' @param similarity_filter_num Integer. Number of papers to filter based on
#'   embedding similarity scores. Default is `20`.
#' @param GPT_filter_num Integer. Number of papers to retain after LLM-based
#'   filtering for relevance. Default is `5`.
#' @param orgdb Object. Organism database for gene annotation. Default is
#'   `org.Hs.eg.db`. Only used for gene enrichment results.
#' @param output_prompt Logical. Whether to include the LLM prompt in the final
#'   annotation results. Default is `TRUE`.
#' @param thinkingBudget Integer. The "thinking budget" parameter specific to
#'   the Gemini API, controlling the depth of reasoning. Default is `0`.
#'
#' @return
#' A `functional_module` class object with updated slots:
#' \describe{
#'   \item{`llm_module_interpretation`}{List containing the final results with
#'     LLM-generated module names and study summaries for each processed module}
#'   \item{`merged_module$functional_module_result`}{Updated data frame with a new
#'     `llm_module_name` column containing LLM-generated module names}
#'   \item{`process_info`}{Updated with parameters and timestamp for this operation}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Basic usage with OpenAI
#' result <- llm_interpret_module(
#'   object = my_functional_module,
#'   api_key = "your_openai_api_key",
#'   embedding_output_dir = "/path/to/embedding/output"
#' )
#'
#' # Advanced usage with local corpus and specific phenotype
#' result <- llm_interpret_module(
#'   object = my_functional_module,
#'   module_content_number_cutoff = 2,
#'   llm_model = "gpt-4o",
#'   api_key = "your_openai_api_key",
#'   embedding_output_dir = "/path/to/embedding/output",
#'   local_corpus_dir = "/path/to/local/papers",
#'   phenotype = "Alzheimer's disease",
#'   years = 3,
#'   retmax = 15,
#'   similarity_filter_num = 30,
#'   GPT_filter_num = 8
#' )
#'}
#'
#' @importFrom dplyr filter mutate left_join
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export

llm_interpret_module <- function(object,
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
                                 thinkingBudget = 0) {

  # 1. Collect functional module result
  if (!is(object, "functional_module")) {
    stop("object should be functional_module class")
  }
  if (all(names(object@process_info) != "merge_modules")) {
    stop("Please use the merge_modules() function to process first")
  }

  if (module_content_number_cutoff >= max(object@merged_module$functional_module_result$module_content_number)) {
    stop("module_content_number_cutoff should be smaller than the maximum of all module content numbers in your functional module result.")
  }

  functional_module_result <-
    object@merged_module$functional_module_result |>
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
    save_dir_local_corpus_embed = "local"
  } else {
    local_corpus <- FALSE
  }

  # 2. Create vector database for local corpus uploaded by user
  clear_output_dir(output_dir = embedding_output_dir) # Clear output directory

  if (!is.null(local_corpus_dir)) {
    embedding_local_corpus(embedding_model = embedding_model,
                           api_provider =  api_provider,
                           api_key = api_key,
                           local_corpus_dir = local_corpus_dir,
                           embedding_output_dir = embedding_output_dir,
                           save_dir_local_corpus_embed = save_dir_local_corpus_embed)
  }

  # 3. Extract pathway description, molecule name, and pmid of references
  processed_data <- preprocess_module(df = functional_module_result,
                                      orgdb = orgdb)

  # 4. Retrieve PubMedIDs of related articles
  pubmed_result <- pubmed_search(processed_data = processed_data,
                                 chunk_size = chunk_size,
                                 years = years,
                                 retmax = retmax)

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
                          embedding_output_dir = embedding_output_dir)

  # 6. Retrieve and rank related papers using RAG strategy
  related_paper <- retrieve_strategy(pubmed_result = pubmed_result,
                                     model = llm_model,
                                     embedding_model = embedding_model,
                                     api_key = api_key,
                                     similarity_filter_num = similarity_filter_num,
                                     GPT_filter_num = GPT_filter_num,
                                     local_corpus = local_corpus,
                                     embedding_output_dir = embedding_output_dir,
                                     save_dir_local_corpus_embed = save_dir_local_corpus_embed)

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
                                         thinkingBudget = thinkingBudget)

  # 8. Store the final result in the object's llm_module_interpretation slot
  object@llm_module_interpretation <- final_result

  # 9. Update graph_data and functional_module result according to llm interpretation
  llm_module_name_df <- data.frame()
  for (i in 1:length(object@llm_module_interpretation)){
    module <- names(object@llm_module_interpretation[i])
    llm_module_name <- object@llm_module_interpretation[[i]]$generated_name$module_name
    llm_module_name_df[i, 1] <- module
    llm_module_name_df[i, 2] <- llm_module_name
  }
  colnames(llm_module_name_df) <- c("module", "llm_module_name")
  object@merged_module$functional_module_result <-
    object@merged_module$functional_module_result |>
    dplyr::left_join(llm_module_name_df, by = "module")

  object@merged_module$functional_module_result <-
    object@merged_module$functional_module_result |>
    mutate(llm_module_name = ifelse(is.na(llm_module_name), yes = module_annotation, no = llm_module_name))

  # 9. Create a process_info entry for this operation using the tidymass_parameter class
  parameter = new(
    Class = "tidymass_parameter",
    pacakge_name = "mapa",
    function_name = "llm_interpret_module()",
    parameter = list(
      module_content_number_cutoff = module_content_number_cutoff,
      llm_model = llm_model,
      embedding_model = embedding_model,
      phenotype = phenotype,
      years = years,
      retmax = retmax,
      similarity_filter_num = similarity_filter_num,
      output_prompt = output_prompt,
      local_corpus = local_corpus
    ),
    time = Sys.time()
  )

  # 10. Add process_info entry to the object
  process_info <- slot(object, "process_info")
  process_info$llm_interpret_module <- parameter
  slot(object, "process_info") <- process_info

  message("Done")

  # 11. Return the updated object
  return(object)
}
