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
#     database = c("go", "kegg", "reactome"),
#     save_to_local = FALSE,
#     path = "result",
#     go.orgdb = org.Hs.eg.db,
#     go.ont = "ALL",
#     kegg.organism = "hsa",
#     reactome.organism = "human",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH"
# )

# gsea_pathways
#
# dir.create("covid_data/result", showWarnings = FALSE)
#
# save(gsea_pathways, file = "covid_data/result/gsea_pathways")

#' Perform Gene Set Enrichment Analysis (GSEA)
#'
#' @param variable_info A data frame containing variable information.
#' @param query_type Character, the category of biological entity to query. Must be one of:
#'   - "gene": for gene/protein-based enrichment
#'   - "metabolite": for metabolite-based enrichment (currently unusable)
#' @param order_by Column name to order genes by. Default is "fc".
#' @param database Databases to use: "go", "kegg", "reactome".
#' @param save_to_local Save results locally. Default is FALSE.
#' @param path Directory to save results. Default is "result".
#'
#' @param go.orgdb OrgDb object for GO analysis.
#' @param go.keytype Key type for GO. Default is "ENTREZID".
#' @param go.ont GO ontology to use. Default is "ALL".
#'
#' @param kegg.organism Organism code for KEGG. Default is "hsa".
#' @param kegg.keytype Key type for KEGG. Default is "kegg".
#'
#' @param reactome.organism Organism for Reactome. Default is "human".
#'   Supported values: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
#'
#' @param exponent GSEA ranking exponent. Default is 1.
#' @param eps Tolerance value. Default is 1e-10.
#' @param verbose Print messages. Default is TRUE.
#' @param seed Set random seed. Default is FALSE.
#' @param pvalueCutoff P-value cutoff. Default is 0.05.
#' @param pAdjustMethod P-value adjustment method. Default is "BH".
#' @param qvalueCutoff Q-value cutoff. Default is 0.2.
#' @param minGSSize Minimum gene set size. Default is 10.
#' @param maxGSSize Maximum gene set size. Default is 500.
#' @param readable Convert IDs to symbols. Default is FALSE.
#' @param ... Additional arguments
#'
#' @return A functional_module object containing enrichment results
#'
#' @export
do_gsea <-
  function(variable_info,
           query_type = "gene",
           order_by = "fc",
           database = c("go", "kegg", "reactome"),
           save_to_local = FALSE,
           path = "result",
           # GO parameters
           go.orgdb = NULL,
           go.keytype = "ENTREZID",
           go.ont = "ALL",
           # KEGG parameters
           kegg.organism = "hsa",
           kegg.keytype = "kegg",
           # Reactome parameters
           reactome.organism = "human",
           # Common parameters
           exponent = 1,
           eps = 1e-10,
           verbose = TRUE,
           seed = FALSE,
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           ...) {

    ###check parameters
    if (readable) {
      readable <- FALSE
    }

    if (missing(database)) {
      stop("database is required")
    }

    if (all(!database %in% c("go", "kegg", "reactome"))) {
      stop("database should contain go, kegg, or reactome")
    }

    if (missing(order_by)) {
      stop("order_by is required")
    }

    if (missing(variable_info)) {
      stop("variable_info is required")
    }

    ###change all columns to lower case
    colnames(variable_info) <- tolower(colnames(variable_info))

    #####check variable_info
    check_variable_info(variable_info = variable_info,
                       query_type = "gene",
                       order_by = order_by)

    variable_info <-
      variable_info %>%
      dplyr::filter(!is.na(entrezid)) %>%
      dplyr::group_by(entrezid) %>%
      dplyr::mutate(p_adjust = min(p_value_adjust)) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(entrezid, .keep_all = TRUE)

    message("Filter variable_info to unique 'entrezid' entries with minimum adjusted p-values.")
    message("'ENTREZID' are used for GSEA.")

    gene_list <-
      log(variable_info[[order_by]], 2)

    names(gene_list) <-
      variable_info$entrezid

    ###GO enrichment
    if ("go" %in% database) {
      message("GO database...")
      if (is.null(go.orgdb)) {
        stop("go.orgdb is required for GO analysis")
      }

      gsea_go_result <-
        clusterProfiler::gseGO(
          geneList = sort(gene_list, decreasing = TRUE),
          ont = go.ont,
          OrgDb = go.orgdb,
          keyType = go.keytype,
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          seed = seed,
          by = "fgsea"
        )

      if (!is.null(gsea_go_result)) {
        gsea_go_result@result$Count <-
          purrr::map_int(gsea_go_result@result$core_enrichment,
                        function(x) {
                          length(stringr::str_split(x, pattern = "/")[[1]])
                        })
        data.table::setnames(gsea_go_result@result, "p.adjust", "p_adjust")
      }
    } else {
      gsea_go_result <- NULL
    }

    ###KEGG
    if ("kegg" %in% database) {
      message("KEGG database...")

      gsea_kegg_result <-
        clusterProfiler::gseKEGG(
          geneList = sort(gene_list, decreasing = TRUE),
          organism = kegg.organism,
          keyType = kegg.keytype,
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          use_internal_data = FALSE,
          seed = seed,
          by = "fgsea"
        )

      if (!is.null(gsea_kegg_result)) {
        gsea_kegg_result@result$Count <-
          purrr::map_int(gsea_kegg_result@result$core_enrichment,
                        function(x) {
                          length(stringr::str_split(x, pattern = "/")[[1]])
                        })
        data.table::setnames(gsea_kegg_result@result, "p.adjust", "p_adjust")
      }
    } else {
      gsea_kegg_result <- NULL
    }

    ###Reactome
    if ("reactome" %in% database) {
      message("Reactome database...")

      gsea_reactome_result <-
        ReactomePA::gsePathway(
          geneList = sort(gene_list, decreasing = TRUE),
          organism = reactome.organism,
          exponent = exponent,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          eps = eps,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          verbose = verbose,
          seed = seed,
          by = "fgsea"
        )

      if (!is.null(gsea_reactome_result)) {
        gsea_reactome_result@result$Count <-
          purrr::map_int(gsea_reactome_result@result$core_enrichment,
                        function(x) {
                          length(stringr::str_split(x, pattern = "/")[[1]])
                        })
        data.table::setnames(gsea_reactome_result@result, "p.adjust", "p_adjust")
      }
    } else {
      gsea_reactome_result <- NULL
    }

    ####save results to local
    if (save_to_local) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)

      if (!is.null(gsea_go_result)) {
        dir.create(file.path(path, "go"), showWarnings = FALSE, recursive = TRUE)
        save(gsea_go_result, file = file.path(path, "go/gsea_go_result"))
      }

      if (!is.null(gsea_kegg_result)) {
        dir.create(file.path(path, "kegg"), showWarnings = FALSE, recursive = TRUE)
        save(gsea_kegg_result, file = file.path(path, "kegg/gsea_kegg_result"))
      }

      if (!is.null(gsea_reactome_result)) {
        dir.create(file.path(path, "reactome"), showWarnings = FALSE, recursive = TRUE)
        save(gsea_reactome_result, file = file.path(path, "reactome/gsea_reactome_result"))
      }
    }

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "do_gsea()",
      parameter = list(
        query_type = query_type,
        order_by = order_by,
        by = by,
        database = database,
        save_to_local = save_to_local,
        path = path,
        go.orgdb = deparse(substitute(go.orgdb)),
        go.keytype = go.keytype,
        go.ont = go.ont,
        kegg.organism = kegg.organism,
        kegg.keytype = kegg.keytype,
        kegg.use_internal_data = FALSE,
        reactome.organism = reactome.organism,
        exponent = exponent,
        eps = eps,
        verbose = verbose,
        seed = seed,
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
    process_info$do_gsea <- parameter

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
