#' Demo Dataset for Over-Representation Analysis
#'
#' Significantly downregulated proteins from the muscle of aging mice
#' (6 vs 30 months, male C57BL/6). Contains 66 proteins with |log2FC| ≥ 0.5
#' and FDR < 0.05, suitable for pathway enrichment analysis.
#'
#' @format A tibble with 66 rows and 3 columns:
#' \describe{
#'   \item{symbol}{Gene symbols}
#'   \item{log2FC (6 vs 30mo)}{Log2 fold changes (all negative)}
#'   \item{FDR (6 vs 30mo)}{False discovery rates (all < 0.05)}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
#' @examples
#' data(demo_data_ora)
#' head(demo_data_ora)
#'
#' # Extract gene symbols for ORA
#' genes <- demo_data_ora$symbol
#'
"demo_data_ora"

#' Demo Dataset for Gene Set Enrichment Analysis
#'
#' Complete proteomics dataset from liver of aging mice (6 vs 30 months,
#' male C57BL/6). Contains 5,290 proteins with fold changes and adjusted
#' p-values, suitable for gene set enrichment analysis (GSEA).
#'
#' @format A tibble with 5,290 rows and 3 columns:
#' \describe{
#'   \item{symbol}{Gene symbols}
#'   \item{fc}{Fold changes (6 vs 30 months)}
#'   \item{p_value_adjust}{Adjusted p-values}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
#' @examples
#' data(demo_data_gsea)
#' head(demo_data_gsea)
#'
#' # Create ranked gene list for GSEA
#' ranked_genes <- setNames(demo_data_gsea$fc, demo_data_gsea$symbol)
#' ranked_genes <- sort(ranked_genes, decreasing = TRUE)
#'
"demo_data_gsea"

#' #' Enriched Functional Module Dataset generated from embedding
#' #'
#' #' A dataset containing results from pathway enrichment analysis and functional module detection.
#' #' The dataset represents a formally structured S4 object of class 'functional_module' from the mapa package,
#' #' containing detailed information about functional modules, pathway enrichment, and their relationships.
#' #'
#' #' @format An S4 object of class 'functional_module'
#' #' @source Generated using the mapa package's pathway enrichment and functional module detection functions.
#' #' @examples
#' #' # Access basic information about the genes
#' #' head(enriched_functional_module@variable_info)
#' #'
#' #' # Explore GO enrichment results
#' #' head(enriched_functional_module@enrichment_go_result@result)
#' #'
#' #' # View KEGG pathway enrichment results
#' #' head(enriched_functional_module@enrichment_kegg_result@result)
#' #'
#' #' # Get functional module information
#' #' head(enriched_functional_module@merged_module$functional_module_result)
#' "embedding_enriched_functional_module"
#'
#' #' Enriched Functional Module Dataset generated mainly by gene overlap
#' #'
#' #' A dataset containing results from pathway enrichment analysis and functional module detection.
#' #' The dataset represents a formally structured S4 object of class 'functional_module' from the mapa package,
#' #' containing detailed information about functional modules, pathway enrichment, and their relationships.
#' #'
#' #' @format An S4 object of class 'functional_module'
#' #' @source Generated using the mapa package's pathway enrichment and functional module detection functions.
#' #' @examples
#' #' # Access basic information about the genes
#' #' head(gene_overlap_enriched_functional_module@variable_info)
#' #'
#' #' # Explore GO enrichment results
#' #' head(gene_overlap_enriched_functional_module@enrichment_go_result@result)
#' #'
#' #' # View KEGG pathway enrichment results
#' #' head(gene_overlap_enriched_functional_module@enrichment_kegg_result@result)
#' #'
#' #' # Get functional module information
#' #' head(gene_overlap_enriched_functional_module@merged_module$functional_module_result)
#' "gene_overlap_enriched_functional_module"
#'
#' #' Enriched Functional Module Dataset generated mainly by metabolite overlap
#' #'
#' #' A dataset containing results from pathway enrichment analysis and functional module detection.
#' #' The dataset represents a formally structured S4 object of class 'functional_module' from the mapa package,
#' #' containing detailed information about functional modules, pathway enrichment, and their relationships.
#' #'
#' #' @format An S4 object of class 'functional_module'
#' #' @source Generated using the mapa package's pathway enrichment and functional module detection functions.
#' #' @examples
#' #' # Access basic information about the metabolites
#' #' head(enriched_functional_module_met@variable_info)
#' #'
#' #' # Explore HMDB enrichment results
#' #' head(enriched_functional_module_met@enrichment_hmdb_result@result)
#' #'
#' #' # View KEGG pathway enrichment results
#' #' head(enriched_functional_module_met@enrichment_metkegg_result@result)
#' #'
#' #' # Get functional module information
#' #' head(enriched_functional_module_met@merged_module$functional_module_result)
#' "enriched_functional_module_met"


#' #' @title enriched_pathways
#' #'
#' #' @description enriched_pathways
#' #' @name enriched_pathways
#' #' @format An enriched_functional_module object.
#' #' @docType data
#' #' @keywords data
#' #' @examples
#' #' data(enriched_pathways)
#' "enriched_pathways"
