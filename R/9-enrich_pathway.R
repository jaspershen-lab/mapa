# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("demo_data.rda")
#
# variable_info <-
#   demo_data %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0) %>%
#   extract_variable_info()
#
# enriched_pathways <-
#   enrich_pathway(
#     object = variable_info,
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

#' Enrich Pathways
#'
#' This function performs enrichment analysis on various databases such as GO, KEGG, and Reactome.
#' It's a generic function that dispatches to specific methods based on the class of the 'object' parameter.
#'
#' @param object A list containing relevant gene information or an object that has a specific method for `enrich_pathway`.
#' @param database Character vector, specify which database(s) to use for enrichment ('go', 'kegg', 'reactome').
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#' @param path Character, the directory where to save the results if save_to_local is TRUE.
#' @param OrgDb Object, the OrgDb object required for GO enrichment.
#' @param organism Character, the organism code.
#' @param keyType Character, the type of key to be used.
#' @param use_internal_data Logical, whether to use internal data for KEGG enrichment.
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
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a sample object
#' sample_object <- list(ensembl = c("ENSG000001", "ENSG000002"), uniprot = c("P12345", "P67890"))
#'
#' # Perform enrichment
#' result <- enrich_pathway(object = sample_object, database = c("go", "kegg"), OrgDb = org.Hs.eg.db)
#'}

enrich_pathway <-
  function(object,
           database = c("go", "kegg", "reactome"),
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
    UseMethod("enrich_pathway")
  }

#' @method enrich_pathway mass_dataset
#' @rdname enrich_pathway
#' @export

enrich_pathway.mass_dataset <-
  function(object,
           database = c("go", "kegg", "reactome"),
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
    variable_info <-
      massdataset::extract_variable_info(object)

    enrich_pathway_internal(
      object = variable_info,
      database = database,
      save_to_local = save_to_local,
      path = path,
      OrgDb = OrgDb,
      organism = organism,
      keyType = keyType,
      use_internal_data = use_internal_data,
      ont = ont,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = readable,
      universe_go,
      universe_kegg,
      universe_reactome,
      pool = pool,
      ...
    )

  }


#' @method enrich_pathway data.frame
#' @rdname enrich_pathway
#' @export

enrich_pathway.data.frame <-
  function(object,
           database = c("go", "kegg", "reactome"),
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
    enrich_pathway_internal(
      object = object,
      database = database,
      save_to_local = save_to_local,
      path = path,
      OrgDb = OrgDb,
      organism = organism,
      keyType = keyType,
      use_internal_data = use_internal_data,
      ont = ont,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      readable = readable,
      universe_go,
      universe_kegg,
      universe_reactome,
      pool = pool,
      ...
    )
  }


#' Internal Function for Enriching Pathways
#'
#' This is an internal function that performs enrichment analysis on various databases such as GO, KEGG, and Reactome.
#' The function takes in a list of genes and additional parameters to perform the enrichment.
#'
#' @param object A list containing relevant gene information.
#' @param database Character vector, specify which database(s) to use for enrichment ('go', 'kegg', 'reactome').
#' @param save_to_local Logical, if TRUE the results will be saved to local disk.
#' @param path Character, the directory where to save the results if save_to_local is TRUE.
#' @param OrgDb Object, the OrgDb object required for GO enrichment.
#' @param organism Character, the organism code.
#' @param keyType Character, the type of key to be used.
#' @param use_internal_data Logical, whether to use internal data for KEGG enrichment.
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
#' @param ... Additional arguments passed to the underlying functions.
#'
#' @return An object of class 'functional_module' containing the enrichment results and parameters.
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#'

enrich_pathway_internal <-
  function(object,
           database = c("go", "kegg", "reactome"),
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
    variable_info <-
      object

    if (readable) {
      readable <- FALSE
    }
    if (missing(database)) {
      stop("database is required")
    }

    if (all(!database %in% c("go", "kegg", "reactome"))) {
      stop("database should contains go, kegg, and reactome")
    }

    #####check variable_info
    check_variable_info(variable_info = variable_info)

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
      message("kegg database...")
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

    result <-
      new(
        "functional_module",
        enrichment_go_result = enrichment_go_result,
        enrichment_kegg_result = enrichment_kegg_result,
        enrichment_reactome_result = enrichment_reactome_result,
        process_info = process_info
      )

    message("Done.")
    result
  }
