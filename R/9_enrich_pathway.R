# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")
# setwd("demo_data/")
#
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
# library(metpath)
#
# load("demo_data.rda")
#
# variable_info <-
#   demo_data %>%
#   massdataset::activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0) %>%
#   massdataset::extract_variable_info()
#
# enriched_pathways <-
#   enrich_pathway(
#     variable_info = variable_info,
#     query_type = "gene",
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
# dir.create("result", showWarnings = FALSE)
# save(enriched_pathways,
#      file = "result/enriched_pathways")

# Metabolite enrichment analysis
# data("query_id_hmdb", package = "metpath")
# data("query_id_kegg", package = "metpath")
# variable_info <- data.frame("hmdbid" = query_id_hmdb, "keggid" = query_id_kegg)

# variable_info <- variable_info_down
# variable_info <- variable_info_up
# variable_info <-
#   variable_info %>%
#   dplyr::rename(hmdbid = HMDB.ID)
# id_conversion <-
#   metpath::hmdb_compound_database@spectra.info %>%
#   dplyr::select(HMDB.ID, KEGG.ID) %>%
#   dplyr::rename(hmdbid = HMDB.ID,
#                 keggid = KEGG.ID)
# variable_info <-
#   variable_info %>%
#   left_join(id_conversion)
#
# enriched_pathways <-
#   enrich_pathway(
#     variable_info = variable_info,
#     query_type = "metabolite",
#     database = c("hmdb", "kegg"),
#     use_internal_data = TRUE,
#     save_to_local = FALSE,
#     path = "result",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH"
#   )


#' Enrich Pathways
#'
#' Performs pathway enrichment analysis using multiple databases (GO, KEGG, Reactome, HMDB). This
#' function acts as a generic interface that dispatches to specific methods based on the type of query: gene-based or metabolite-based.
#' For gene queries, GO, KEGG, and Reactome enrichment are performed;
#' for metabolite queries, HMDB and KEGG enrichment are performed.
#'
#' @param variable_info A data.frame containing relevant gene information or
#' an object that has a specific method for `enrich_pathway`.
#' @param query_type Character, the category of biological entity to query ("gene", "metabolite").
#' @param database Character vector, specify which database(s) to use for enrichment ('go', 'kegg', 'reactome', 'hmdb').
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#' @param path Character, the directory where to save the results if save_to_local is TRUE.
#' @param OrgDb Object, the OrgDb object required for GO enrichment.
#' @param organism Character, the organism code.
#' @param keyType Character, the type of key to be used.
#' @param use_internal_data Logical, whether to use internal database for KEGG and HMDB enrichment.
#' @param ont Character, the ontology to be used for GO enrichment ('ALL', 'BP', 'CC', 'MF').
#' @param pvalueCutoff Numeric, the p-value cutoff for enrichment.
#' @param pAdjustMethod Character, the method for adjusting p-values.
#' @param universe_go Numeric vector, the universe for GO enrichment.
#' @param universe_kegg Numeric vector, the universe for KEGG enrichment.
#' @param universe_reactome Numeric vector, the universe for Reactome enrichment.
#' @param qvalueCutoff Numeric, the q-value cutoff for enrichment.
#' @param minGSSize Numeric, the minimum gene set size for enrichment.
#' @param maxGSSize Numeric, the maximum gene set size for enrichment.
#' @param readable Logical, whether to make the results more human-readable.
#' @param pool Logical, whether to use clusterProfiler's pooling feature.
#' @param ... Additional arguments passed to the underlying methods.
#'
#' @return An object containing the enrichment results and parameters.
#'
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom ReactomePA enrichPathway
#' @importFrom metpath enrich_hmdb enrich_kegg
#' @importFrom methods new
#' @import org.Hs.eg.db
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @export
#'
#' @examples
#' \dontrun{
#' data(demo_data)
#'
#' variable_info <-
#'   demo_data %>%
#'   massdataset::activate_mass_dataset(what = "variable_info") %>%
#'   dplyr::filter(fdr < 0.05 & score > 0) %>%
#'   massdataset::extract_variable_info()
#' require(org.Hs.eg.db)
#' enriched_pathways <-
#'   enrich_pathway(
#'     variable_info = variable_info,
#'     save_to_local = FALSE,
#'     path = "result",
#'     OrgDb = org.Hs.eg.db,
#'     organism = "hsa",
#'     database = c("kegg", "reactome"),
#'     ont = "ALL",
#'     pvalueCutoff = 0.05,
#'     pAdjustMethod = "BH",
#'     qvalueCutoff = 0.2,
#'     minGSSize = 10,
#'     maxGSSize = 500,
#'     readable = FALSE,
#'     pool = FALSE
#'   )
#'   }


enrich_pathway <-
  function(variable_info,
           query_type = c("gene", "metabolite"),
           database = c("go", "kegg", "reactome", "hmdb"),
           save_to_local = FALSE,
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe_go = NULL,
           universe_kegg = NULL,
           universe_reactome = NULL,
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           pool = FALSE,
           ...) {
    ###check the id of markers
    if (readable) {
      readable <- FALSE
    }

    if (missing(query_type)) {
      stop("query_type is required")
    }
    query_type <- match.arg(query_type)

    if (missing(database)) {
      stop("database is required")
    }
    if (all(!database %in% c("go", "kegg", "reactome", "hmdb"))) {
      stop("database should contains go, kegg, reactome and hmdb")
    }

    #####check variable_info
    check_variable_info(variable_info = variable_info, query_type = query_type)

    if (query_type == "gene") {
      ###go enrichment
      if ("go" %in% database) {
        message("GO database...")
        if (missing(OrgDb)) {
          stop("OrgDb is required")
        }
        if (is.null(universe_go)) {
          enrichment_go_result <-
            clusterProfiler::enrichGO(
              gene = variable_info$ensembl[!is.na(variable_info$ensembl)],
              OrgDb = OrgDb,
              keyType = "ENSEMBL",
              ont = "ALL",
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              pool = pool
            )
        } else{
          enrichment_go_result <-
            clusterProfiler::enrichGO(
              gene = variable_info$ensembl[!is.na(variable_info$ensembl)],
              OrgDb = OrgDb,
              keyType = "ENSEMBL",
              ont = "ALL",
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = universe_go,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              pool = pool
            )
        }
      } else{
        enrichment_go_result <- NULL
      }

      ###KEGG
      if ("kegg" %in% database) {
        message("KEGG database...")
        if (is.null(universe_kegg)) {
          enrichment_kegg_result <-
            clusterProfiler::enrichKEGG(
              gene = variable_info$uniprot[!is.na(variable_info$uniprot)],
              organism = organism,
              keyType = "uniprot",
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              use_internal_data = use_internal_data
            )
        } else{
          enrichment_kegg_result <-
            clusterProfiler::enrichKEGG(
              gene = variable_info$uniprot[!is.na(variable_info$uniprot)],
              organism = organism,
              keyType = "uniprot",
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = universe_kegg,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              use_internal_data = use_internal_data
            )
        }
      } else{
        enrichment_kegg_result <- NULL
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
        if (is.null(universe_reactome)) {
          enrichment_reactome_result <-
            ReactomePA::enrichPathway(
              gene = variable_info$entrezid[!is.na(variable_info$entrezid)],
              organism = organism2,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable
            )
        } else{
          enrichment_reactome_result <-
            ReactomePA::enrichPathway(
              gene = variable_info$entrezid[!is.na(variable_info$entrezid)],
              organism = organism2,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              universe = universe_reactome
            )
        }
      } else{
        enrichment_reactome_result <- NULL
      }
    } else if (query_type == "metabolite") {
      #### Get HMDB and KEGG databases
      dbs <- get_met_pathway_db(use_internal_data = use_internal_data)

      #### HMDB pathway enrichment
      if ("hmdb" %in% database) {
        message("HMDB enrichment ...")
        enrichment_hmdb_result <-
          metpath::enrich_hmdb(
            query_id = variable_info$hmdbid,
            query_type = "compound",
            id_type = "HMDB",
            pathway_database = dbs$hmdb_pathway,
            only_primary_pathway = TRUE,
            p_cutoff = pvalueCutoff,
            p_adjust_method = pAdjustMethod,
            threads = 3
          )
      } else {
        enrichment_hmdb_result <- NULL
      }

      #### KEGG pathway enrichment
      if ("kegg" %in% database) {
        message("KEGG enrichment ...")
        enrichment_metkegg_result <-
          metpath::enrich_kegg(
            query_id = variable_info$keggid,
            query_type = "compound",
            id_type = "KEGG",
            pathway_database = dbs$kegg_pathway,
            p_cutoff = pvalueCutoff,
            p_adjust_method = pAdjustMethod,
            threads = 3
          )
      } else {
        enrichment_metkegg_result <- NULL
      }
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
      save(enrichment_go_result,
           file = file.path(path, "go/enrichment_go_result"))
      save(enrichment_kegg_result,
           file = file.path(path, "kegg/enrichment_kegg_result"))
      save(
        enrichment_reactome_result,
        file = file.path(path, "reactome/enrichment_reactome_result")
      )
    }

    parameter = new(
      Class = "tidymass_parameter",
      pacakge_name = "mapa",
      function_name = "enrich_pathway()",
      parameter = list(
        query_type = query_type,
        database = database,
        save_to_local = save_to_local,
        path = path,
        organism = organism,
        keyType = keyType,
        use_internal_data = use_internal_data,
        ont = ont,
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = pAdjustMethod,
        qvalueCutoff = qvalueCutoff,
        universe_go = universe_go,
        universe_kegg = universe_kegg,
        universe_reactome = universe_reactome,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        readable = readable,
        pool = pool
      ),
      time = Sys.time()
    )

    process_info <- list()

    process_info$enrich_pathway <-
      parameter

    if (query_type == "gene") {
      enrichment_hmdb_result <- NULL
      enrichment_metkegg_result <- NULL
    } else {
      enrichment_go_result <- NULL
      enrichment_kegg_result <- NULL
      enrichment_reactome_result <- NULL
    }

    result <-
      new(
        "functional_module",
        variable_info = variable_info,
        enrichment_go_result = enrichment_go_result,
        enrichment_kegg_result = enrichment_kegg_result,
        enrichment_reactome_result = enrichment_reactome_result,
        enrichment_hmdb_result = enrichment_hmdb_result,
        enrichment_metkegg_result = enrichment_metkegg_result,
        process_info = process_info
      )

    message("Done.")
    result
  }
