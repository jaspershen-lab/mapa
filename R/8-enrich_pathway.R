# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(ReactomePA)
#
# load("demo_data")
#
# object <-
#   demo_data %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::filter(fdr < 0.05 & score > 0)
#
# enrich_pathway(object,
#                path = "result",
#                OrgDb = org.Hs.eg.db,
#                organism = "hsa")


#' Perform pathway enrichment analysis
#'
#' This function is a wrapper for enrich_pathway_internal. It dispatches the call to the
#' appropriate method based on the class of the object argument. This functions
#' is based on the clusterProfiler and ReactomePA package.
#'
#' @param object A data frame or mass_dataset object containing gene identifiers.
#' @param database A vector of databases to use for enrichment (default: c("GO", "kegg", "reactome")).
#' @param path The directory where results should be saved (default: "result").
#' @param OrgDb The OrgDb object to use for GO enrichment.
#' @param organism The organism code for KEGG and Reactome enrichment (default: "hsa").
#' @param keyType The key type for gene identifiers (default: "ENTREZID").
#' @param use_internal_data Logical, use internal data for KEGG enrichment (default: FALSE).
#' @param ont The ontology to use for GO enrichment (default: "ALL").
#' @param pvalueCutoff The p-value cutoff for enrichment (default: 0.05).
#' @param pAdjustMethod The method for p-value adjustment (default: "BH").
#' @param universe The background gene set.
#' @param qvalueCutoff The q-value cutoff for enrichment (default: 0.2).
#' @param minGSSize The minimum gene set size for enrichment (default: 10).
#' @param maxGSSize The maximum gene set size for enrichment (default: 500).
#' @param readable Logical, make the output readable (default: FALSE).
#' @param pool Logical, pool genes for GO enrichment (default: FALSE).
#' @param ... Additional arguments passed to the underlying functions.
#' @return NULL
#' @export

enrich_pathway <-
  function(object,
           database = c("go", "kegg", "reactome"),
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe,
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
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe,
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           pool = FALSE,
           ...) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    ###check the id of markers
    variable_info <-
      massdataset::extract_variable_info(object)

    if (missing(universe)) {
      enrich_pathway_internal(
        object = variable_info,
        database = database,
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
        pool = pool,
        ...
      )
    } else{
      enrich_pathway_internal(
        object = variable_info,
        database = database,
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
        universe = universe,
        pool = pool,
        ...
      )
    }
  }





#' @method enrich_pathway data.frame
#' @rdname enrich_pathway
#' @export

enrich_pathway.data.frame <-
  function(object,
           database = c("go", "kegg", "reactome"),
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe,
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           pool = FALSE,
           ...) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    if (missing(universe)) {
      enrich_pathway_internal(
        object = object,
        database = database,
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
        pool = pool,
        ...
      )
    } else{
      enrich_pathway_internal(
        object = object,
        database = database,
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
        universe = universe,
        pool = pool,
        ...
      )
    }
  }






#' Perform pathway enrichment analysis
#'
#' @param object A data frame containing gene identifiers.
#' @param database A vector of databases to use for enrichment (default: c("GO", "kegg", "reactome")).
#' @param path The directory where results should be saved (default: "result").
#' @param OrgDb The OrgDb object to use for GO enrichment. For example, human is org.Hs.eg.db.
#' @param organism The organism code for KEGG and Reactome enrichment (default: "hsa").
#' @param keyType The key type for gene identifiers (default: "ENTREZID").
#' @param use_internal_data Logical, use internal data for KEGG enrichment (default: FALSE).
#' @param ont The ontology to use for GO enrichment (default: "ALL").
#' @param pvalueCutoff The p-value cutoff for enrichment (default: 0.05).
#' @param pAdjustMethod The method for p-value adjustment (default: "BH").
#' @param universe The background gene set.
#' @param qvalueCutoff The q-value cutoff for enrichment (default: 0.2).
#' @param minGSSize The minimum gene set size for enrichment (default: 10).
#' @param maxGSSize The maximum gene set size for enrichment (default: 500).
#' @param readable Logical, make the output readable (default: FALSE).
#' @param pool Logical, pool genes for GO enrichment (default: FALSE).
#' @param ... Additional arguments passed to the underlying functions.
#' @return NULL

enrich_pathway_internal <-
  function(object,
           database = c("go", "kegg", "reactome"),
           path = "result",
           OrgDb,
           organism = "hsa",
           keyType = "ENTREZID",
           use_internal_data = FALSE,
           ont = "ALL",
           pvalueCutoff = 0.05,
           pAdjustMethod = "BH",
           universe,
           qvalueCutoff = 0.2,
           minGSSize = 10,
           maxGSSize = 500,
           readable = FALSE,
           pool = FALSE,
           ...) {
    ###check the id of markers
    variable_info <-
      object
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    if (readable) {
      readable <- FALSE
    }
    if (missing(database)) {
      stop("database is required")
    }

    if (all(!database %in% c("go", "kegg", "reactome"))) {
      stop("database should contains go, kegg, and reactome")
    }

    if (all(colnames(variable_info) != "ensembl")) {
      stop("ensembl should be in the variable_info")
    } else{
      if (all(is.na(variable_info$ensembl))) {
        stop("All ensembl column are NA")
      }
    }

    if (all(colnames(variable_info) != "symbol")) {
      stop("symbol should be in the variable_info")
    } else{
      if (all(is.na(variable_info$symbol))) {
        stop("All symbol column are NA")
      }
    }

    if (all(colnames(variable_info) != "uniprot")) {
      stop("uniprot should be in the variable_info")
    } else{
      if (all(is.na(variable_info$uniprot))) {
        stop("All uniprot column are NA")
      }
    }

    if (all(colnames(variable_info) != "entrezid")) {
      stop("entrezid should be in the variable_info")
    } else{
      if (all(is.na(variable_info$entrezid))) {
        stop("All entrezid column are NA")
      }
    }

    ###go enrichment
    if ("go" %in% database) {
      message("GO database...")
      if (missing(OrgDb)) {
        stop("OrgDb is required")
      }
      if (missing(universe)) {
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
            universe = universe,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            readable = readable,
            pool = pool
          )
      }
      save(enrichment_go_result,
           file = file.path(path, "enrichment_go_result"))
    }


    ###KEGG
    if ("kegg" %in% database) {
      message("kegg database...")
      if (missing(universe)) {
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
            universe = universe,
            minGSSize = minGSSize,
            maxGSSize = maxGSSize,
            use_internal_data = use_internal_data
          )
      }
      save(enrichment_kegg_result,
           file = file.path(path, "enrichment_kegg_result"))
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
      if (missing(universe)) {
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
            universe = universe
          )
      }
      save(enrichment_reactome_result,
           file = file.path(path, "enrichment_reactome_result"))
    }
    message("Done.")
  }
