#' @title Demo data
#'
#' @description Demo data
#' @name demo_data
#' @format A mass_dataset object.
#' @docType data
#' @keywords data
#' @examples
#' data(demo_data)
"demo_data"

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
