# following code: Yuchen's testing code
#
# library(metid)
# library(metpath)
# library(tidymass)
#
# load("/do_gsea_demoData/meta.rda")
#
# plot(log(meta_info$score, 2))
#
# gsea_pathways <- do_gsea(
#   variable_info = meta_info,
#   database = c("kegg", "hmdb"),
#   order_by = "score",
#   analysis_type = "meta",
#   pvalueCutoff = 1,
#   minGSSize = 0,
#   maxGSSize = 1000
# )
#
# meta_kegg_gsea_res <-
#   as.data.frame(gsea_pathways@enrichment_metkegg_result@result)
#
# meta_hmdb_gsea_res <-
#   as.data.frame(gsea_pathways@enrichment_hmdb_result@result)
#
#
# ##### test: transcript gsea ######
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("/do_gsea_demoData/transc.rda")
#
# plot(log(transc_info$score, 2))
#
# sum(is.na(log(transc_info$score, 2)))
# sum(is.infinite(log(transc_info$score, 2)))
#
# gsea_pathways <-
#   do_gsea(
#     variable_info = transc_info,
#     analysis_type = "transc",
#     order_by = "score",
#     OrgDb = org.Hs.eg.db,
#     database = c("go", "kegg", "reactome"),
#     ont = "ALL",
#     pool = FALSE,
#     pvalueCutoff = 1
#   )
#
# transc_go_gsea_res <-
#   as.data.frame(gsea_pathways@enrichment_go_result@result)
#
# transc_kegg_gsea_res <-
#   as.data.frame(gsea_pathways@enrichment_kegg_result@result)
#
# transc_reactome_gsea_res <-
#   as.data.frame(gsea_pathways@enrichment_reactome_result@result)

# following code: Yifei's initial testing code
#
# library(r4projects)
# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("covid_data/covid_data.RData")
#
# variable_info <-
#   covid_data |>
#   massdataset::activate_mass_dataset(what = "variable_info") |>
#   massdataset::extract_variable_info() #|>
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

#' Perform Gene Set Enrichment Analysis (GSEA) for multi-omics data
#'
#' @description
#' This integrated function performs GSEA for both transcriptomics and metabolomics data:
#' \itemize{
#'   \item \strong{Transcriptomics} (analysis_type = "transc"):
#'   Supports GO/KEGG/Reactome databases. Requires Entrez IDs and fold change values.
#'
#'   \item \strong{Metabolomics} (analysis_type = "meta"):
#'   Supports KEGG/HMDB metabolic pathways. Requires KEGG IDs, HMDB IDs and fold change values.
#' }
#'
#' @param variable_info A data.frame containing feature information with required identifiers:
#'                     - For metabolomics: 'keggid' (KEGG compound IDs) and/or 'hmdbid' (HMDB IDs)
#'                     - For transcriptomics: 'entrezid' (Entrez Gene IDs)
#' @param order_by Character specifying the column in variable_info to use for ranking features.
#'                Typically contains fold change values. Must be specified. Default is 'fc'.
#' @param analysis_type Analysis type: "transc" for transcriptomics or "meta" for metabolomics.
#'                     Mandatory parameter with no default.
#' @param database Character vector specifying databases to use:
#'                - For "transc": c("go", "kegg", "reactome")
#'                - For "meta": c("kegg", "hmdb")
#' @param save_to_local Logical indicating whether to save results to local files.
#' Default is `FALSE`.
#' @param path Output directory path when save_to_local=TRUE.
#' Default is `"result"`.
#' @param OrgDb Annotation database (e.g., org.Hs.eg.db). Required for transcriptomics GO analysis. Required if `"go"` is included in `database`.
#' @param organism Organism code. Supported: "hsa" (human), "rno" (rat), "mmu" (mouse), etc.
#'               Default is `"hsa"`.
#' @param keyType ID type for transcriptomics (default is `"ENTREZID"`).
#' @param exponent Weighting exponent for GSEA (default is `1`).
#' @param eps Numeric threshold for p-value adjustment (default is `1e-10`).
#' @param verbose Logical controlling progress messages (default is `TRUE`).
#' @param seed Logical/numeric for reproducibility (default is `FALSE`).
#' @param by Algorithm selection (`"fgsea"` for fast implementation, default).
#' @param use_internal_data Logical to use KEGG internal data (default is `FALSE`).
#' @param ont GO ontology type (`"BP"`, `"MF"`, or `"CC"`). Required for GO analysis.
#' @param pvalueCutoff Adjusted p-value cutoff (default is `0.05`).
#' @param pAdjustMethod P-value adjustment method (default is `"BH"`).
#' @param qvalueCutoff q-value cutoff (default is `0.2`).
#' @param minGSSize Minimum gene set size (default is `10`).
#' @param maxGSSize Maximum gene set size (default is `500`).
#' @param readable Logical for converting IDs to gene symbols (default is `FALSE`).
#' @param ... Additional parameters passed to internal functions.
#'
#' @return An object of class `functional_module` containing the enrichment analysis results.
#'
#' @author Xiaotao Shen \email{xiaotao.shen@outlook.com}
#'
#' @examples
#' \dontrun{
#' # Transcriptomics example
#' library(org.Hs.eg.db)
#' data(gene_data) # Contains entrezid and fc columns
#' result_transc <- do_gsea(
#'   analysis_type = "transc",
#'   variable_info = gene_data,
#'   database = c("go", "kegg"),
#'   OrgDb = org.Hs.eg.db,
#'   organism = "hsa"
#' )
#'
#' # Metabolomics example
#' data(metabo_data) # Contains keggid and fc columns
#' result_meta <- do_gsea(
#'   analysis_type = "meta",
#'   variable_info = metabo_data,
#'   database = "kegg",
#'   organism = "hsa"
#' )
#'
#' # Extract KEGG pathway results
#' kegg_res <- result_meta@meta_gsea_kegg_res
#' }
#'
#' @export

do_gsea <-
  function(variable_info,
           order_by = "fc",
           analysis_type = c("transc", "meta"),
           database,
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
           ont,
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           ...) {
    ### complete this argument
    analysis_type <- match.arg(analysis_type)

    ### check if type are clarified
    if (!missing(analysis_type)) {
      if (!analysis_type %in% c("transc", "meta")) {
        stop("Invalid analysis_type detected. Choose 'transc' for transcriptomics gsea or 'meta' for metabolomics gsea.")
      }
      if (length(analysis_type) != 1) {
        stop("Incorrect analysis_type input. Please choose one and only one from 'transc' or 'meta'.")
      } else if (analysis_type == "transc") {
        do_gsea_transc(
          variable_info = variable_info,
          order_by = order_by,
          database = database,
          save_to_local = save_to_local,
          path = path,
          OrgDb,
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
          readable = readable,
          ...
        )
      } else if (analysis_type == "meta") {
        do_gsea_meta(
          variable_info = variable_info,
          order_by = order_by,
          database = database,
          save_to_local = save_to_local,
          path = path,
          organism = organism,
          exponent = exponent,
          eps = eps,
          verbose = verbose,
          seed = seed,
          by = by,
          use_internal_data = use_internal_data,
          pvalueCutoff = pvalueCutoff,
          pAdjustMethod = pAdjustMethod,
          qvalueCutoff = qvalueCutoff,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          readable = readable,
          ...
        )
      }
    } else {
      stop("Please specify the type of omic data for GSEA analysis: use 'transc' for transcriptomics or 'meta' for metabolomics.")
    }
  }

do_gsea_transc <- function(variable_info = variable_info,
                           order_by = order_by,
                           database = database,
                           save_to_local = save_to_local,
                           path = path,
                           OrgDb,
                           organism = organism,
                           keyType = "ENTREZID",
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
                           readable = readable,
                           ...) {
  ### check the id of markers
  if (readable) {
    readable <- FALSE
  }

  if (missing(database)) {
    stop("Please specify which database you intend to use.")
  }

  ### check valid db
  valid_db <- c("go", "kegg", "reactome")
  db_error_msg <- "Transcriptomics database options: go/kegg/reactome. You are welcome to select one or multiple of them."

  if (any(!database %in% valid_db)) {
    stop(paste("Invalid database for", analysis_type, ":", db_error_msg))
  }

  if (missing(order_by)) {
    stop("order_by must be specified.")
  }

  if (missing(variable_info)) {
    stop("variable_info must be specified.")
  }

  ### change all the columns to low cases
  colnames(variable_info) <- tolower(colnames(variable_info))

  ##### check variable_info
  check_variable_info(
    variable_info = variable_info,
    query_type = "gene",
    order_by = order_by
  )

  ### query cleaning
  variable_info <-
    variable_info |>
    dplyr::filter(!is.na(entrezid)) |>
    dplyr::group_by(entrezid) |>
    dplyr::mutate(p_adjust = min(p_value_adjust)) |>
    dplyr::ungroup() |>
    dplyr::distinct(entrezid, .keep_all = TRUE)

  message("Filter variable_info to unique 'entrezid' entries with minimum adjusted p-values.")

  gene_list <-
    log(variable_info[[order_by]], 2)

  ### final to-enrich object: gene_list
  names(gene_list) <-
    variable_info$entrezid

  gene_list <-
    gene_list[is.finite(gene_list)]

  ### go enrichment
  if ("go" %in% database) {
    message(
      "\n",
      "==========================\n",
      "Now running GO database...\n",
      "==========================\n"
    )

    if (missing(OrgDb)) {
      stop("OrgDb is required")
    }

    gene_gsea_go_res <-
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

    gene_gsea_go_res@result$count <-
      purrr::map_int(
        gene_gsea_go_res@result$core_enrichment,
        function(x) {
          length(stringr::str_split(x, pattern = "/")[[1]])
        }
      )
  } else {
    gene_gsea_go_res <- NULL
  }

  ### KEGG
  if ("kegg" %in% database) {
    message(
      "\n",
      "============================\n",
      "Now running KEGG database...\n",
      "============================\n"
    )

    gene_gsea_kegg_res <-
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

    gene_gsea_kegg_res@result$count <-
      purrr::map_int(
        gene_gsea_kegg_res@result$core_enrichment,
        function(x) {
          length(stringr::str_split(x, pattern = "/")[[1]])
        }
      )
  } else {
    gene_gsea_kegg_res <- NULL
  }

  ### Reactome
  if ("reactome" %in% database) {
    ## change the organism name
    organism2 <-
      switch(organism,
        hsa = "human",
        rno = "rat",
        mmu = "mouse",
        cel = "celegans",
        sce = "yeast",
        dre = "zebrafis",
        dme = "fly",
        mde = "fly"
      )

    message(
      "\n",
      "================================\n",
      "Now running Reactome database...\n",
      "================================\n"
    )


    gene_gsea_reactome_res <-
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

    gene_gsea_reactome_res@result$count <-
      purrr::map_int(
        gene_gsea_reactome_res@result$core_enrichment,
        function(x) {
          length(stringr::str_split(x, pattern = "/")[[1]])
        }
      )
  } else {
    gene_gsea_reactome_res <- NULL
  }

  #### save results to local?
  if (save_to_local) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(path, "go"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(file.path(path, "kegg"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    dir.create(file.path(path, "reactome"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    save(gene_gsea_go_res, file = file.path(path, "go/gene_gsea_go_res"))
    save(gene_gsea_kegg_res, file = file.path(path, "kegg/gene_gsea_kegg_res"))
    save(gene_gsea_reactome_res,
      file = file.path(path, "reactome/gene_gsea_reactome_res")
    )
  }

  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "mapa",
    function_name = "do_gsea()",
    parameter = list(
      variable_info = variable_info,
      order_by = order_by,
      database = database,
      save_to_local = save_to_local,
      path = path,
      OrgDb,
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
      enrichment_go_result = gene_gsea_go_res,
      enrichment_kegg_result = gene_gsea_kegg_res,
      enrichment_reactome_result = gene_gsea_reactome_res,
      process_info = process_info
    )

  message("Done.")
  result
}

do_gsea_meta <- function(variable_info = variable_info,
                         order_by = order_by,
                         database = database,
                         save_to_local = save_to_local,
                         path = path,
                         organism = organism,
                         exponent = exponent,
                         eps = eps,
                         verbose = verbose,
                         seed = seed,
                         by = by,
                         use_internal_data = use_internal_data,
                         pvalueCutoff = pvalueCutoff,
                         pAdjustMethod = pAdjustMethod,
                         minGSSize = minGSSize,
                         maxGSSize = maxGSSize,
                         readable = readable,
                         ...) {
  ### check the id of markers
  if (readable) {
    readable <- FALSE
  }

  if (missing(database)) {
    stop("Please specify which database you intend to use.")
  }

  ### check valid db
  valid_db <- c("kegg", "hmdb")

  db_error_msg <- "Metabolomics database options: kegg/hmdb. You can also select multiple of them."

  if (any(!database %in% valid_db)) {
    stop(paste("Invalid database for", analysis_type, ":", db_error_msg))
  }

  if (missing(order_by)) {
    stop("order_by must be specified.")
  }

  if (missing(variable_info)) {
    stop("variable_info cannot be empty.")
  }

  check_variable_info(
    variable_info = variable_info,
    query_type = "metabolite",
    order_by = order_by
  )

  ### initialize
  meta_gsea_kegg_res <- NULL
  meta_gsea_hmdb_res <- NULL

  get_kegg_pathway <- function() {
    data("kegg_hsa_pathway", envir = environment())
    message(
      crayon::yellow(
        "This database is downloaded in",
        kegg_hsa_pathway@database_info$version
      )
    )
    cat("\n")
    return(kegg_hsa_pathway)
  }

  get_hmdb_pathway <- function() {
    data("hmdb_pathway", envir = environment())
    message(
      crayon::yellow(
        "This database is downloaded in",
        hmdb_pathway@database_info$version
      )
    )
    cat("\n")
    return(hmdb_pathway)
  }

  if ("kegg" %in% database) {
    message(
      "\n",
      "============================\n",
      "Now running KEGG database...\n",
      "============================\n"
    )

    # process variable_info
    variable_info <- variable_info |>
      dplyr::filter(!is.na(keggid)) |>
      dplyr::group_by(keggid) |>
      dplyr::mutate(p_adjust = min(p_value_adjust)) |>
      dplyr::ungroup() |>
      dplyr::distinct(keggid, .keep_all = TRUE)

    # check KEGG ID format
    if (!any(grepl("^C\\d{5}$", variable_info$keggid))) {
      stop("All KEGG ID does not match compound format (e.g. C00001)")
    } else {
      # load KEGG pathways from metpath
      kegg_hsa_pathway <-
        get_kegg_pathway()

      # Arrange database
      kegg_ids <- sapply(kegg_hsa_pathway@compound_list, function(x) x$KEGG.ID)
      kegg_pathways <- kegg_hsa_pathway@pathway_name
      metabolite_strs <- sapply(kegg_ids, function(x) paste(x, collapse = ","))

      term2metabolite_kegg <- data.frame(
        pathway_name = kegg_pathways,
        metabolite_ids = metabolite_strs,
        stringsAsFactors = FALSE
      )

      term2metabolite_kegg <- term2metabolite_kegg |>
        tidyr::separate_rows(metabolite_ids, sep = ",") |>
        dplyr::mutate(metabolite_ids = trimws(metabolite_ids)) |>
        dplyr::filter(!is.na(metabolite_ids) & metabolite_ids != "") |>
        dplyr::distinct()

      meta_list <- log(variable_info[[order_by]], 2)
      names(meta_list) <- variable_info$keggid

      meta_gsea_kegg_res <- clusterProfiler::GSEA(
        geneList = sort(meta_list, decreasing = TRUE),
        TERM2GENE = term2metabolite_kegg,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        eps = eps,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        verbose = verbose,
        seed = seed,
        by = by
        # GSEA function in clusterProfiler
        # does not have `qvalueCutoff` parameter
      )

      meta_gsea_kegg_res@result$count <-
        purrr::map_int(
          meta_gsea_kegg_res@result$core_enrichment,
          function(x) {
            length(stringr::str_split(x, pattern = "/")[[1]])
          }
        )

      meta_gsea_kegg_res@result <-
        subset(meta_gsea_kegg_res@result, select = -ID)
      meta_gsea_kegg_res@result <-
        meta_gsea_kegg_res@result |>
        rename(pathway = Description)
    }
  }

  if ("hmdb" %in% database) {
    message(
      "\n",
      "============================\n",
      "Now running HMDB database...\n",
      "============================\n"
    )

    # process variable_info
    variable_info <- variable_info |>
      dplyr::filter(!is.na(hmdbid)) |>
      dplyr::group_by(hmdbid) |>
      dplyr::mutate(p_adjust = min(p_value_adjust)) |>
      dplyr::ungroup() |>
      dplyr::distinct(hmdbid, .keep_all = TRUE)

    # check HMDB ID format
    if (!any(grepl("^HMDB\\d{7}$", variable_info$hmdbid))) {
      stop("All HMDB ID does not match compound format (e.g. HMDB0000001)")
    } else {
      # load HMDB pathways from metpath
      hmdb_hsa_pathway <-
        get_hmdb_pathway()

      # Arrange database
      hmdb_ids <- sapply(hmdb_hsa_pathway@compound_list, function(x) x$HMDB.ID)
      hmdb_pathways <- hmdb_hsa_pathway@pathway_name
      metabolite_strs <- sapply(hmdb_ids, function(x) paste(x, collapse = ","))

      term2metabolite_hmdb <- data.frame(
        pathway_name = hmdb_pathways,
        metabolite_ids = metabolite_strs,
        stringsAsFactors = FALSE
      )

      term2metabolite_hmdb <- term2metabolite_hmdb |>
        tidyr::separate_rows(metabolite_ids, sep = ",") |>
        dplyr::mutate(metabolite_ids = trimws(metabolite_ids)) |>
        dplyr::filter(!is.na(metabolite_ids) & metabolite_ids != "") |>
        dplyr::distinct()

      meta_list <- log(variable_info[[order_by]], 2)
      names(meta_list) <- variable_info$hmdbid

      meta_gsea_hmdb_res <- clusterProfiler::GSEA(
        geneList = sort(meta_list, decreasing = TRUE),
        TERM2GENE = term2metabolite_hmdb,
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

      meta_gsea_hmdb_res@result$count <-
        purrr::map_int(
          meta_gsea_hmdb_res@result$core_enrichment,
          function(x) {
            length(stringr::str_split(x, pattern = "/")[[1]])
          }
        )

      meta_gsea_hmdb_res@result <-
        subset(meta_gsea_hmdb_res@result, select = -ID)
      meta_gsea_hmdb_res@result <-
        meta_gsea_hmdb_res@result |>
        rename(pathway = Description)
    }
  }

  ### save results to local?
  if (save_to_local) {
    dir.create(file.path(path, "kegg"), showWarnings = FALSE)
    dir.create(file.path(path, "hmdb"), showWarnings = FALSE)

    save(meta_gsea_kegg_res,
      file = file.path(path, "kegg/meta_gsea_kegg_res")
    )

    save(meta_gsea_hmdb_res,
      file = file.path(path, "hmdb/meta_gsea_hmdb_res")
    )
  }

  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "mapa",
    function_name = "do_gsea()",
    parameter = list(
      variable_info = variable_info,
      order_by = order_by,
      database = database,
      save_to_local = save_to_local,
      path = path,
      organism = organism,
      exponent = exponent,
      eps = eps,
      verbose = verbose,
      seed = seed,
      by = by,
      use_internal_data = use_internal_data,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = readable
    ),
    time = Sys.time()
  )

  process_info <- list()

  process_info$do_gsea <- parameter

  result <- new(
    "functional_module",
    variable_info = variable_info,
    enrichment_metkegg_result = meta_gsea_kegg_res,
    enrichment_hmdb_result = meta_gsea_hmdb_res,
    process_info = process_info
  )

  message("Done.")
  result
}
