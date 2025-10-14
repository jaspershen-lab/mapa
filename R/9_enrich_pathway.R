# setwd(r4projects::get_project_wd())
# source("R/6_utils.R")
# source("R/8_functional_module_class.R")

# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
# library(metpath)

# setwd("demo_data/")
# load("demo_data.rda")
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
#     database = c("go", "kegg", "reactome"),
#     save_to_local = FALSE,
#     go.orgdb = org.Hs.eg.db,
#     go.keytype = "ENTREZID",
#     go.ont = "ALL",
#     kegg.organism = "hsa",
#     kegg.keytype = "kegg",
#     reactome.organism = "human",
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

# enriched_pathways <-
#   enrich_pathway(
#     variable_info = variable_info_up,
#     query_type = "metabolite",
#     met_organism = "hsa",
#     database = c("metkegg", "hmdb"),
#     save_to_local = FALSE,
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH"
#   )

#' Enrich Pathways
#'
#' @description
#' Performs comprehensive pathway enrichment analysis using multiple databases (GO, KEGG, Reactome, HMDB).
#' This function serves as a unified interface that dispatches to specific methods based on the query type:
#' - For gene queries: GO, KEGG, and Reactome enrichment
#' - For metabolite queries: HMDB and KEGG enrichment
#'
#' @param variable_info A data.frame containing biological entity information for enrichment.
#'   For gene-based queries, it should contain columns like 'ensembl', 'uniprot', and 'entrezid'.
#'   For metabolite-based queries, it should contain columns like 'hmdbid' and 'keggid'.
#' @param query_type Character, the category of biological entity to query. Must be one of:
#'   - "gene": for gene/protein-based enrichment
#'   - "metabolite": for metabolite-based enrichment
#' @param database Character vector, specify which database(s) to use for enrichment:
#'   - For genes: 'go', 'kegg', 'reactome'
#'   - For metabolites: 'hmdb', 'metkegg'. 'hmdb' is only for human.
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#' @param path Character, the directory where to save the results if save_to_local is TRUE.
#'
#' @param go.orgdb  OrgDb object or character string naming the *OrgDb* annotation package required for GO enrichment (e.g., org.Hs.eg.db or "org.Hs.eg.db").
#'   See Bioconductor OrgDb packages for available organisms.
#' @param go.keytype Character, the key type to be used for GO enrichment. Default is "ENTREZID".
#'   See available keytypes in your OrgDb with `keytypes(OrgDb)`.
#' @param go.ont Character, the GO ontology to be used:
#'   - "ALL": all ontologies (default)
#'   - "BP": biological process
#'   - "CC": cellular component
#'   - "MF": molecular function
#'   See clusterProfiler::enrichGO() for details.
#' @param go.universe Numeric vector, the background universe for GO enrichment.
#'   If NULL (default), all genes in the OrgDb will be used.
#' @param go.pool Logical, if ont='ALL', whether pool 3 GO sub-ontologies.
#'
#' @param kegg.organism Character, the organism code (e.g., "hsa" for human) required for KEGG enrichment.
#'   See 'https://www.genome.jp/kegg/catalog/org_list.html' for supported organism codes.
#' @param kegg.keytype Character, the type of key to be used for KEGG. Default is "kegg".
#'   One of "kegg", "ncbi-geneid", "ncbi-proteinid", or "uniprot".
#'   Note: 'kegg' ID is typically entrezgene ID for eukaryotes and Locus ID for prokaryotes.
#' @param kegg.universe Numeric vector, the background universe for KEGG enrichment.
#'   If NULL (default), all genes in KEGG database will be used.
#'
#' @param reactome.organism Character, the organism name for Reactome (e.g., "human").
#'   Supported values: "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
#' @param reactome.universe Numeric vector, the background universe for Reactome enrichment.
#'   If NULL (default), all genes in Reactome database will be used.
#'
#' @param met_organism Character, the organism code for metabolite enrichment.
#'   Default is "hsa" (human). See 'https://www.genome.jp/kegg/catalog/org_list.html' for supported organism codes.
#'
#' @param pvalueCutoff Numeric, the p-value cutoff for enrichment.
#' @param pAdjustMethod Character, the method for adjusting p-values.
#'   See stats::p.adjust() for available methods ("BH", "bonferroni", etc.).
#' @param qvalueCutoff Numeric, the q-value cutoff for enrichment. Default is 0.2.
#'   q-values are specifically those adjusted p-values that incorporate an empirical estimate of
#'   the fraction of true nulls, aiming for maximum power at a fixed FDR.
#' @param minGSSize Numeric, the minimum gene set size for enrichment.
#'   Sets with fewer genes than this will be excluded.
#' @param maxGSSize Numeric, the maximum gene set size for enrichment.
#'   Sets with more genes than this will be excluded.
#' @param readable Logical, whether to make the gene identifiers in results more human-readable.
#'   See clusterProfiler documentation for details.
#' @param ... Additional arguments passed to the underlying enrichment functions:
#'   clusterProfiler::enrichGO(), clusterProfiler::enrichKEGG(),
#'   ReactomePA::enrichPathway(), metpath::enrich_hmdb(), or metpath::enrich_kegg().
#'
#' @return An object of class "functional_module" containing:
#'   - variable_info: The input variable information
#'   - enrichment_go_result: GO enrichment results (for gene queries)
#'   - enrichment_kegg_result: KEGG enrichment results (for gene queries)
#'   - enrichment_reactome_result: Reactome enrichment results (for gene queries)
#'   - enrichment_hmdb_result: HMDB enrichment results (for metabolite queries)
#'   - enrichment_metkegg_result: KEGG enrichment results (for metabolite queries)
#'   - process_info: Parameters and processing information
#'
#' @importFrom clusterProfiler enrichGO enrichKEGG
#' @importFrom ReactomePA enrichPathway
#' @importFrom metpath enrich_hmdb enrich_kegg
#' @importFrom methods new
#' @import org.Hs.eg.db
#'
#' @note
#' - For more details on GO and KEGG enrichment parameters, see ?clusterProfiler::enrichGO and ?clusterProfiler::enrichKEGG
#' - For Reactome enrichment parameters, see ?ReactomePA::enrichPathway
#' - For metabolite enrichment parameters, see ?metpath::enrich_hmdb and ?metpath::enrich_kegg
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Gene-based enrichment
#' enriched_pathways <-
#'   enrich_pathway(
#'     variable_info = variable_info,
#'     query_type = "gene",
#'     save_to_local = FALSE,
#'     go.orgdb = org.Hs.eg.db,
#'     go.keytype = "ENSEMBL",
#'     kegg.organism = "hsa",
#'     kegg.keytype = "uniprot",
#'     reactome.organism = "human",
#'     database = c("kegg", "reactome"),
#'     pvalueCutoff = 0.05,
#'     pAdjustMethod = "BH",
#'     qvalueCutoff = 0.2,
#'     minGSSize = 10,
#'     maxGSSize = 500,
#'     readable = FALSE
#'   )
#'
#' # Metabolite-based enrichment
#' # Prepare metabolite IDs
#' data("query_id_hmdb", package = "metpath")
#' data("query_id_kegg", package = "metpath")
#' met_variable_info <- data.frame("hmdbid" = query_id_hmdb, "keggid" = query_id_kegg)
#'
#' # Perform metabolite pathway enrichment
#' met_enriched_pathways <-
#'   enrich_pathway(
#'     variable_info = met_variable_info,
#'     query_type = "metabolite",
#'     database = c("hmdb", "metkegg"),
#'     save_to_local = FALSE,
#'     pvalueCutoff = 0.05,
#'     pAdjustMethod = "BH"
#'   )
#' }

enrich_pathway <-
  function(variable_info,
           query_type = c("gene", "metabolite"),
           database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
           save_to_local = FALSE,
           path = "result",
           # GO-specific parameters
           go.orgdb = org.Hs.eg.db,
           go.keytype = "ENTREZID",
           go.ont = "ALL",
           go.universe = NULL,
           go.pool = FALSE,
           # KEGG-specific parameters
           kegg.organism = "hsa",
           kegg.keytype = "kegg",
           kegg.universe = NULL,
           # Reactome-specific parameters
           reactome.organism = c("human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly", "bovine", "canine", "chicken"),
           ## reactome input ID should be Entrez ID
           reactome.universe = NULL,
           # Metabolite specific parameters
           met_organism = "hsa",
           # Common parameters
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

    if (missing(query_type)) {
      stop("query_type is required")
    }
    query_type <- match.arg(query_type)

    if (missing(database)) {
      stop("database is required")
    }
    database <- match.arg(database, several.ok = TRUE)

    #####check variable_info
    check_variable_info(variable_info = variable_info, query_type = query_type)

    # Initialize results to NULL
    enrichment_go_result <- NULL
    enrichment_kegg_result <- NULL
    enrichment_reactome_result <- NULL
    enrichment_hmdb_result <- NULL
    enrichment_metkegg_result <- NULL

    if (query_type == "gene") {
      ###go enrichment
      if ("go" %in% database) {
        message("GO database...")
        if (is.null(go.orgdb)) {
          stop("go.orgdb is required for GO enrichment")
        }

        # Convert keyType to lowercase for case-insensitive matching with variable_info columns
        go_keytype_lower <- tolower(go.keytype)

        # Try to find a matching column directly in variable_info
        col_match <- which(tolower(colnames(variable_info)) == go_keytype_lower)

        if (length(col_match) > 0) {
          col_name <- colnames(variable_info)[col_match]
        } else {
          stop(paste("Could not find a matching column in variable_info for keyType:", go.keytype,
                     "\nAvailable columns:", paste(colnames(variable_info), collapse=", ")))
        }

        # Extract the gene IDs from the appropriate column
        go_ids <- variable_info[[col_name]]

        enrichment_go_result <- tryCatch({
          if (is.null(go.universe)) {
            clusterProfiler::enrichGO(
              gene = go_ids[!is.na(go_ids)],
              OrgDb = go.orgdb,
              keyType = go.keytype,
              ont = go.ont,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              pool = go.pool
            )
          } else {
            clusterProfiler::enrichGO(
              gene = go_ids[!is.na(go_ids)],
              OrgDb = go.orgdb,
              keyType = go.keytype,
              ont = go.ont,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = go.universe,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              pool = go.pool
            )
          }
        }, error = function(e) {
          warning(paste("GO enrichment failed:", e$message))
          return(NULL)
        })

        if (!is.null(enrichment_go_result)) {
          data.table::setnames(enrichment_go_result@result, "p.adjust", "p_adjust")
        }
      }

      ###KEGG
      if ("kegg" %in% database) {
        message("KEGG database...")

        # Get the appropriate ID column based on kegg.keytype
        kegg_ids <- switch(kegg.keytype,
                           "uniprot" = variable_info$uniprot,
                           "kegg" = variable_info$entrezid,
                           "ncbi-geneid" = variable_info$entrezid,
                           "ncbi-proteinid" = variable_info$uniprot,
                           stop(paste("Unsupported kegg.keytype:", kegg.keytype)))

        enrichment_kegg_result <- tryCatch({
          if (is.null(kegg.universe)) {
            clusterProfiler::enrichKEGG(
              gene = kegg_ids[!is.na(kegg_ids)],
              organism = kegg.organism,
              keyType = kegg.keytype,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              use_internal_data = FALSE
            )
          } else {
            clusterProfiler::enrichKEGG(
              gene = kegg_ids[!is.na(kegg_ids)],
              organism = kegg.organism,
              keyType = kegg.keytype,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              universe = kegg.universe,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              use_internal_data = FALSE
            )
          }
        }, error = function(e) {
          warning(paste("KEGG enrichment failed:", e$message))
          return(NULL)
        })

        if (!is.null(enrichment_kegg_result)) {
          data.table::setnames(enrichment_kegg_result@result, "p.adjust", "p_adjust")
        }
      }

      ###Reactome
      if ("reactome" %in% database) {
        reactome.organism <- match.arg(reactome.organism)

        message("Reactome database...")
        enrichment_reactome_result <- tryCatch({
          if (is.null(reactome.universe)) {
            ReactomePA::enrichPathway(
              gene = variable_info$entrezid[!is.na(variable_info$entrezid)],
              organism = reactome.organism,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable
            )
          } else {
            ReactomePA::enrichPathway(
              gene = variable_info$entrezid[!is.na(variable_info$entrezid)],
              organism = reactome.organism,
              pvalueCutoff = pvalueCutoff,
              pAdjustMethod = pAdjustMethod,
              qvalueCutoff = qvalueCutoff,
              minGSSize = minGSSize,
              maxGSSize = maxGSSize,
              readable = readable,
              universe = reactome.universe
            )
          }
        }, error = function(e) {
          warning(paste("Reactome enrichment failed:", e$message))
          return(NULL)
        })

        if (!is.null(enrichment_reactome_result)) {
          data.table::setnames(enrichment_reactome_result@result, "p.adjust", "p_adjust")
        }
      }
    } else if (query_type == "metabolite") {
      if (met_organism == "hsa") {
      #### Get HMDB and KEGG databases for human
      dbs <- tryCatch({
        get_met_pathway_db(use_internal_data = TRUE)
      }, error = function(e) {
        warning(paste("Failed to get metabolite pathway databases:", e$message))
        return(list(hmdb_pathway = NULL, kegg_pathway = NULL))
      })

      #### HMDB pathway enrichment (human only)
      if ("hmdb" %in% database && !is.null(dbs$hmdb_pathway)) {
        message("HMDB enrichment ...")
        enrichment_hmdb_result <- tryCatch({
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
        }, error = function(e) {
        warning(paste("HMDB enrichment failed:", e$message))
        return(NULL)
        })

        if (!is.null(enrichment_hmdb_result)) {
        data.table::setnames(enrichment_hmdb_result@result, "p_value_adjust", "p_adjust")
        }
      }

      #### KEGG pathway enrichment for human
      if ("metkegg" %in% database && !is.null(dbs$kegg_pathway)) {
        message("KEGG enrichment for human ...")
        enrichment_metkegg_result <- tryCatch({
        metpath::enrich_kegg(
          query_id = variable_info$keggid,
          query_type = "compound",
          id_type = "KEGG",
          pathway_database = dbs$kegg_pathway,
          p_cutoff = pvalueCutoff,
          p_adjust_method = pAdjustMethod,
          threads = 3
        )
        }, error = function(e) {
        warning(paste("Metabolite KEGG enrichment failed:", e$message))
        return(NULL)
        })
      }
      } else {
      #### KEGG pathway enrichment for other organisms
      if ("metkegg" %in% database) {
        message(paste("KEGG enrichment for", met_organism, "..."))
        pathway_database <- tryCatch({
        get_kegg_pathways(organism = met_organism, local = FALSE)
        }, error = function(e) {
        warning(paste("Failed to get KEGG pathway database:", e$message))
        return(NULL)
        })

        if (!is.null(pathway_database)) {
        enrichment_metkegg_result <- tryCatch({
          metpath::enrich_kegg(
          query_id = variable_info$keggid,
          query_type = "compound",
          id_type = "KEGG",
          pathway_database = pathway_database,
          p_cutoff = pvalueCutoff,
          p_adjust_method = pAdjustMethod,
          threads = 3
          )
        }, error = function(e) {
          warning(paste("Metabolite KEGG enrichment failed:", e$message))
          return(NULL)
        })
        }
      }
      }

      if (!is.null(enrichment_metkegg_result)) {
      data.table::setnames(enrichment_metkegg_result@result, "p_value_adjust", "p_adjust")
      }
    }

    ####save results to local?
    if (save_to_local) {
      tryCatch({
        dir.create(path, recursive = TRUE, showWarnings = FALSE)

        if (!is.null(enrichment_go_result)) {
          dir.create(file.path(path, "go"), showWarnings = FALSE, recursive = TRUE)
          save(enrichment_go_result, file = file.path(path, "go/enrichment_go_result"))
        }

        if (!is.null(enrichment_kegg_result)) {
          dir.create(file.path(path, "kegg"), showWarnings = FALSE, recursive = TRUE)
          save(enrichment_kegg_result, file = file.path(path, "kegg/enrichment_kegg_result"))
        }

        if (!is.null(enrichment_reactome_result)) {
          dir.create(file.path(path, "reactome"), showWarnings = FALSE, recursive = TRUE)
          save(enrichment_reactome_result, file = file.path(path, "reactome/enrichment_reactome_result"))
        }

        if (!is.null(enrichment_hmdb_result)) {
          dir.create(file.path(path, "hmdb"), showWarnings = FALSE, recursive = TRUE)
          save(enrichment_hmdb_result, file = file.path(path, "hmdb/enrichment_hmdb_result"))
        }

        if (!is.null(enrichment_metkegg_result)) {
          dir.create(file.path(path, "metkegg"), showWarnings = FALSE, recursive = TRUE)
          save(enrichment_metkegg_result, file = file.path(path, "metkegg/enrichment_metkegg_result"))
        }
      }, error = function(e) {
        warning(paste("Failed to save results to local path:", e$message))
      })
    }

    if (!is.character(go.orgdb)) {
      go.orgdb <- deparse(substitute(go.orgdb))
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
        go.orgdb = go.orgdb,
        go.keytype = go.keytype,
        go.ont = go.ont,
        go.universe = if(is.null(go.universe)) NULL else "provided",
        go.pool = go.pool,
        kegg.organism = kegg.organism,
        kegg.keytype = kegg.keytype,
        kegg.universe = if(is.null(kegg.universe)) NULL else "provided",
        reactome.organism = reactome.organism,
        reactome.keytype = "entrezid",
        reactome.universe = if(is.null(reactome.universe)) NULL else "provided",
        met_organism = met_organism,
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
    process_info$enrich_pathway <- parameter

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
