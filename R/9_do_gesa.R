# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")

# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("covid_data/covid_data.RData")

# variable_info <-
#   covid_data %>%
#   massdataset::activate_mass_dataset(what = "variable_info") %>%
#   massdataset::extract_variable_info() #%>%
#   #dplyr::arrange(desc(fc))
#
# plot(log(variable_info$fc, 2))
#
# sum(is.na(log(variable_info$fc, 2)))
# sum(is.infinite(log(variable_info$fc, 2)))

# gsea_pathways <-
#   do_gsea(
#     variable_info = variable_info,
#     order_by = "fc",
#     save_to_local = FALSE,
#     path = "result",
#     OrgDb = org.Hs.eg.db,
#     organism = "hsa",
#     database = c("go", "kegg", "reactome"),
#     ont = "ALL",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.2,
#     minGSSize = 10,
#     maxGSSize = 500,
#     readable = FALSE,
#     pool = FALSE
#   )

# gsea_pathways
#
# dir.create("covid_data/result", showWarnings = FALSE)
#
# save(gsea_pathways, file = "covid_data/result/gsea_pathways")

#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs Gene Set Enrichment Analysis using specified databases (GO, KEGG, Reactome)
#' on provided variable information.
#'
#' @param variable_info A data frame containing variable information, including gene identifiers
#'   and associated statistics.
#' @param order_by A character string specifying the column name in `variable_info` to order genes by.
#'   Default is `"fc"`.
#' @param database A character vector specifying the databases to use for enrichment analysis.
#'   Options are `"go"`, `"kegg"`, and `"reactome"`.
#' @param save_to_local Logical value indicating whether to save the results locally.
#'   Default is `FALSE`.
#' @param path A character string specifying the directory path to save results if `save_to_local` is `TRUE`.
#'   Default is `"result"`.
#' @param OrgDb An `OrgDb` object for the organism of interest. Required if `"go"` is included in `database`.
#' @param organism A character string specifying the organism code (e.g., `"hsa"` for human).
#'   Default is `"hsa"`.
#' @param keyType A character string specifying the type of gene identifiers used.
#'   Default is `"ENTREZID"`.
#' @param exponent Numeric value specifying the exponent in the GSEA algorithm.
#'   Default is `1`.
#' @param eps Numeric value specifying the tolerance in the GSEA calculation.
#'   Default is `1e-10`.
#' @param verbose Logical value indicating whether to print messages during the analysis.
#'   Default is `TRUE`.
#' @param seed Logical value indicating whether to set a random seed for reproducibility.
#'   Default is `FALSE`.
#' @param by Character string specifying the method to combine p-values.
#'   Default is `"fgsea"`.
#' @param use_internal_data Logical value indicating whether to use internal data for KEGG analysis.
#'   Default is `FALSE`.
#' @param ont Character string specifying the ontology to use in GO analysis.
#'   Default is `"ALL"`.
#' @param pvalueCutoff Numeric value specifying the p-value cutoff for significance.
#'   Default is `0.05`.
#' @param pAdjustMethod Character string specifying the method for p-value adjustment.
#'   Default is `"BH"`.
#' @param qvalueCutoff Numeric value specifying the q-value cutoff for significance.
#'   Default is `0.2`.
#' @param minGSSize Integer specifying the minimum size of gene sets to include in the analysis.
#'   Default is `10`.
#' @param maxGSSize Integer specifying the maximum size of gene sets to include in the analysis.
#'   Default is `500`.
#' @param readable Logical value indicating whether to convert gene IDs to gene Symbols in the results.
#'   Default is `FALSE`.
#' @param ... Additional arguments passed to the enrichment functions.
#'
#' @return An object of class `functional_module` containing the enrichment analysis results.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' variable_info <- data.frame(
#'   entrezid = c("1", "2", "3"),
#'   uniprot = c("P12345", "Q67890", "A1B2C3"),
#'   fc = c(2.5, -1.3, 0.8)
#' )
#'
#' library(org.Hs.eg.db)
#' result <- do_gsea(
#'   variable_info = variable_info,
#'   order_by = "fc",
#'   database = c("go", "kegg"),
#'   OrgDb = org.Hs.eg.db,
#'   organism = "hsa"
#' )
#' }
#'
#' @export

do_gsea <-
  function(variable_info,
           order_by = "fc",
           database = c("go", "kegg", "reactome"),
           save_to_local = FALSE,
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           exponent = 1,
           eps = 1e-10,
           verbose = TRUE,
           seed = FALSE,
           by = "fgsea",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           ...) {
    ###check the id of markers
    if (readable) {
      readable <- FALSE
    }
    if (missing(database)) {
      stop("database is required")
    }

    if (all(!database %in% c("go", "kegg", "reactome"))) {
      stop("database should contains go, kegg, and reactome")
    }

    if (missing(order_by)) {
      stop("order_by is required")
    }

    if (missing(variable_info)) {
      stop("variable_info is required")
    }

    ###chang all the columns to low cases
    colnames(variable_info) <- tolower(colnames(variable_info))

    #####check variable_info
    check_variable_info(variable_info = variable_info, query_type = "gene" ,order_by = order_by)

    variable_info <-
      variable_info %>%
      dplyr::filter(!is.na(entrezid)) %>%
      dplyr::group_by(entrezid) %>%
      dplyr::mutate(p_adjust = min(p_value_adjust)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(entrezid, .keep_all = TRUE)
    message("Filter variable_info to unique 'entrezid' entries with minimum adjusted p-values.")

    gene_list <-
      log(variable_info[[order_by]], 2)

    names(gene_list) <-
      variable_info$entrezid

    ###go enrichment
    if ("go" %in% database) {
      message("GO database...")
      if (missing(OrgDb)) {
        stop("OrgDb is required")
      }

      gsea_go_result <-
        clusterProfiler::gseGO(
          geneList = sort(gene_list, decreasing = TRUE),
          ont = ont,
          OrgDb = OrgDb,
          keyType = keyType,
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          seed = seed,
          by = by
        )
      gsea_go_result@result$Count <- purrr::map_int(gsea_go_result@result$core_enrichment,
                                                   function(x) {
                                                     length(stringr::str_split(x, pattern = "/")[[1]])
                                                   })
    } else{
      gsea_go_result <- NULL
    }

    ###KEGG
    if ("kegg" %in% database) {
      message("kegg database...")

      gsea_kegg_result <-
        clusterProfiler::gseKEGG(
          geneList = sort(gene_list, decreasing = TRUE),
          organism = organism,
          keyType = "kegg",
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          use_internal_data = use_internal_data,
          seed = seed,
          by = by
        )
      gsea_kegg_result@result$Count <- purrr::map_int(gsea_kegg_result@result$core_enrichment,
                                                   function(x) {
                                                     length(stringr::str_split(x, pattern = "/")[[1]])
                                                   })
    } else{
      gsea_kegg_result <- NULL
    }

    ###Reactome
    if ("reactome" %in% database) {
      ##change the organism name
      organism2 <-
        switch(
          organism,
          hsa = "human",
          rno = "rat",
          mmu = "mouse",
          cel = "celegans",
          sce = "yeast",
          dre = "zebrafis",
          dme = "fly",
          mde = "fly"
        )
      message("Reactome database...")

      gsea_reactome_result <-
        ReactomePA::gsePathway(
          geneList = sort(gene_list, decreasing = TRUE),
          organism = organism2,
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          seed = seed,
          by = by
        )
      gsea_reactome_result@result$Count <- purrr::map_int(gsea_reactome_result@result$core_enrichment,
                                                   function(x) {
                                                     length(stringr::str_split(x, pattern = "/")[[1]])
                                                     })
    } else{
      gsea_reactome_result <- NULL
    }

    ####save results to local?
    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(path, "go"),
                 showWarnings = FALSE,
                 recursive = TRUE)
      dir.create(file.path(path, "kegg"),
                 showWarnings = FALSE,
                 recursive = TRUE)
      dir.create(file.path(path, "reactome"),
                 showWarnings = FALSE,
                 recursive = TRUE)
      save(gsea_go_result, file = file.path(path, "go/gsea_go_result"))
      save(gsea_kegg_result, file = file.path(path, "kegg/gsea_kegg_result"))
      save(gsea_reactome_result,
           file = file.path(path, "reactome/gsea_reactome_result"))
    }

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "do_gsea()",
      parameter = list(
        order_by = order_by,
        database = database,
        save_to_local = save_to_local,
        path = path,
        organism = organism,
        keyType = keyType,
        exponent = exponent,
        eps = eps,
        verbose = verbose,
        seed = seed,
        by = by,
        use_internal_data = use_internal_data,
        ont = ont,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = qvalueCutoff,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        readable = readable
      ),
      time = Sys.time()
    )

    process_info <- list()

    process_info$do_gsea <-
      parameter

    result <-
      new(
        "functional_module",
        variable_info = variable_info,
        enrichment_go_result = gsea_go_result,
        enrichment_kegg_result = gsea_kegg_result,
        enrichment_reactome_result = gsea_reactome_result,
        process_info = process_info
      )

    message("Done.")
    result
  }
