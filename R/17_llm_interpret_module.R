# source("R/17_llm_module_pathway_infor.R")
# source("R/17_llm_module_online_retrieval.R")
# source("R/17_llm_module_embedding_database.R")
# source("R/17_llm_module_RAG_strategy.R")
# source("R/17_llm_module_utils.R")
# source('R/17_llm_module_output_generation.R')

# Metabolite
# functional_module_annotation <-
#   llm_interpret_module(
#     object = enriched_functional_module,
#     phenotype = "pregnancy",
#     api_key = api_key,
#     embedding_output_dir = "demo_data/pregnancy_data/results/results_biotext/embedding_output/",
#     local_corpus = FALSE
#   )

# Gene
# functional_module_annotation <-
#   llm_interpret_module(
#     object = enriched_functional_module,
#     api_key = api_key,
#     embedding_output_dir = "demo_data/updated_object_results_for_genes/gene_overlap_result/embedding_output/",
#     local_corpus = FALSE
#   )

#' Interpret Functional Module using LLM Integrated with RAG Strategy
#'
#' @description This function processes functional module results by retrieving
#'   relevant papers using a Retrieval-Augmented Generation (RAG) strategy. It includes pathway
#'   description extraction, PubMed searching, embedding of search results, and retrieving
#'   related papers with ranking. Finally, input the functional module information and
#'   abstract and title of retrieved and filtered papers into LLM to generate a name and a summary
#'   for each functional module.
#'
#' @param object A functional_module class object.
#' @param api_key Character string. API key for OpenAI or other embedding service. (Currently, only API key for OpenAI can be used)
#' @param embedding_output_dir Character string. Directory where embedding results will be saved.
#' @param local_corpus Logical. Whether to use local files. Default is FALSE.
#' @param local_corpus_dir Character string. Directory containing local files provided by users.
#'   Required if local_corpus is TRUE.
#' @param save_dir_local_corpus_embed Character string. Directory to save embedded local corpus.
#'   Required if local_corpus is TRUE.
#' @param phenotype Character string. Phenotype or disease to focus on. Default is NULL.
#' @param chunk_size Integer. Chunk size for processing data. Default is 5.
#' @param years Integer. Number of recent years to search in PubMed. Default is 5.
#' @param retmax Integer. Maximum number of records to return from PubMed. Default is 10.
#' @param similarity_filter_num Integer. Number of papers to filter based on similarity. Default is 20.
#' @param GPT_filter_num Integer. Number of papers to filter using GPT. Default is 5.
#' @param orgdb Object. Organism database for gene annotation, default is org.Hs.eg.db. Only used for gene enrichment results.
#'
#' @return A list containing the final results with module names and study summaries.
#'
#' @importFrom methods is
#' @importFrom dplyr filter mutate
#' @import org.Hs.eg.db
#'
#' @author Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export

llm_interpret_module <- function(object,
                                 api_key,
                                 embedding_output_dir,
                                 local_corpus = FALSE,
                                 local_corpus_dir = NULL,
                                 save_dir_local_corpus_embed = NULL,
                                 phenotype = NULL,
                                 chunk_size = 5,
                                 years = 5,
                                 retmax = 10,
                                 similarity_filter_num = 20,
                                 GPT_filter_num = 5,
                                 orgdb = org.Hs.eg.db) {

  # 1. Collect functional module result
  if (!is(object, "functional_module")) {
    stop("object should be functional_module class")
  }
  if (all(names(object@process_info) != "merge_modules")) {
    stop("Please use the merge_modules() function to process first")
  }
  functional_module_result <- object@merged_module$functional_module_result

  # Check if local corpus parameters are provided when local_corpus is TRUE
  if (local_corpus) {
    if (is.null(local_corpus_dir) || is.null(save_dir_local_corpus_embed)) {
      stop("When local_corpus is TRUE, both local_corpus_dir and save_dir_local_corpus_embed must be provided")
    }
  }

  # 2. Create vector database for local corpus uploaded by user
  clear_output_dir(output_dir = embedding_output_dir) # Clear output directory

  if (local_corpus) {
    embedding_local_corpus(api_key = api_key,
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
                          api_key = api_key,
                          embedding_output_dir = embedding_output_dir)

  # 6. Retrieve and rank related papers using RAG strategy
  related_paper <- retrieve_strategy(pubmed_result = pubmed_result,
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
  final_result <- module_name_generation(paper_result = paper_result,
                                         phenotype = phenotype,
                                         api_key = api_key)

  return(final_result)
}
