#' Process ID conversion
id_conversion <- function(query_type = c("gene", "metabolite"),
                          data = NULL,
                          from_id_type = NULL,
                          to_id_type = NULL,
                          organism = NULL) {
  if (missing(query_type)){
    stop("query_type is missing")
  }
  query_type <- match.arg(query_type, c("gene", "metabolite"))

  if (query_type == "gene") {
    if (!requireNamespace("clusterProfiler", quietly = TRUE))
      BiocManager::install("clusterProfiler")
    library(clusterProfiler)

    converted <- clusterProfiler::bitr(
      geneID  = data[[tolower(from_id_type)]],
      fromType = from_id_type,
      toType   = to_id_type,
      OrgDb    = organism
    ) %>%
      dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
      dplyr::rename_with(tolower) %>%
      dplyr::left_join(data, ., by = tolower(from_id_type))

    conversion_code <- sprintf(
      "
      if (!requireNamespace(\"clusterProfiler\", quietly = TRUE)) {
        BiocManager::install(\"clusterProfiler\")
      }
      library(clusterProfiler)
      converted <- clusterProfiler::bitr(
        geneID  = data[[tolower(\"%s\")]],
        fromType = \"%s\",
        toType   = c(\"%s\"),
        OrgDb    = %s
       ) %%>%%
       dplyr::distinct(ENTREZID, .keep_all = TRUE) %%>%%
       dplyr::rename_with(tolower) %%>%%
       dplyr::left_join(data, ., by = tolower(\"%s\"))",
      from_id_type, from_id_type,
      paste(to_id_type, collapse = "\", \""), deparse(substitute(organism)), from_id_type)

    return(list(converted_id   = converted,
                conversion_code = conversion_code))
  }

  if (query_type == "metabolite" && organism == "hsa") {
    if (!requireNamespace("metpath", quietly = TRUE)) {BiocManager::install("metpath")}
    library(metpath)

    id_lookup <- metpath::hmdb_compound_database@spectra.info %>%
      dplyr::select(HMDB.ID, KEGG.ID) %>%
      dplyr::rename(hmdbid = HMDB.ID,
                    keggid = KEGG.ID) %>%
      dplyr::mutate(across(everything(), as.character)) %>%
      dplyr::distinct()

    converted <- data |>
      dplyr::filter(!is.na(.data[[tolower(from_id_type)]])) |>
      dplyr::left_join(id_lookup, by = tolower(from_id_type))

    conversion_code <- sprintf(
    "
    if (!requireNamespace(\"metpath\", quietly = TRUE)) {
      BiocManager::install(\"metpath\")
    }
    library(metpath)

    id_lookup <- metpath::hmdb_compound_database@spectra.info %%>%%
      dplyr::select(HMDB.ID, KEGG.ID) %%>%%
      dplyr::rename(hmdbid = HMDB.ID, keggid = KEGG.ID)

    converted <- data %%>%%
      dplyr::filter(!is.na(.data[[\"%s\"]])) %%>%%
      dplyr::left_join(id_lookup, by = \"%s\")",
    from_id_type,
    from_id_type)

    return(list(converted_id   = converted,
                conversion_code = conversion_code))
  }

  if (query_type == "metabolite" && organism != "hsa") {
    return(list(converted_id   = data,
                conversion_code = NA_character_))
  }
}

# Organism name conversion
org_kegg_2name <- c(
  "org.Hs.eg.db"  = "Human",     "org.Mm.eg.db"  = "Mouse",
  "org.Rn.eg.db"  = "Rat",       "org.Dm.eg.db"  = "Fly",
  "org.Dr.eg.db"  = "Zebrafish", "org.At.tair.db"= "Arabidopsis",
  "org.Sc.sgd.db" = "Yeast", "org.Ce.eg.db"  = "C. elegans",
  "org.Ss.eg.db" = "Pig", "org.Bt.eg.db"  = "Bovine",
  "org.Mmu.eg.db" = "Rhesus", "org.Cf.eg.db"  = "Canine",
  "org.EcK12.eg.db" = "E. coli K-12", "org.EcSakai.eg.db" = "E. coli Sakai",
  "org.Gg.eg.db"  = "Chicken",  "org.Xl.eg.db"  = "Xenopus",
  "org.Pt.eg.db"  = "Chimp",    "org.Ag.eg.db"  = "Anopheles",
  "org.Pf.plasmo.db" = "Malaria", "org.Mxanthus.db" = "Myxococcus xanthus",
  "hsa" = "Human", "mmu" = "Mouse",
  "rno" = "Rat", "dme" = "Fly",
  "dre" = "Zebrafish", "sce" = "Yeast",
  "cel" = "C. elegans", "ssc" = "Pig",
  "bta" = "Bovine", "cfa" = "Canine"
)
org2kegg <- c(
  "org.Hs.eg.db"  = "hsa",  "org.Mm.eg.db"  = "mmu",
  "org.Rn.eg.db"  = "rno",  "org.Dm.eg.db"  = "dme",
  "org.Dr.eg.db"  = "dre",  "org.At.tair.db"= "ath",
  "org.Sc.sgd.db" = "sce",  "org.Ce.eg.db"  = "cel",
  "org.Ss.eg.db"  = "ssc",  "org.Bt.eg.db"  = "bta",
  "org.Mmu.eg.db" = "mcc",  "org.Cf.eg.db"  = "cfa",
  "org.EcK12.eg.db" = "eco","org.EcSakai.eg.db" = "ecs",
  "org.Gg.eg.db"  = "gga",  "org.Xl.eg.db"  = "xla",
  "org.Pt.eg.db"  = "ptr",  "org.Ag.eg.db"  = "aga",
  "org.Pf.plasmo.db" = "pfa","org.Mxanthus.db" = "mxa"
)
org2react <- c(
  "org.Hs.eg.db"  = "human",     "org.Mm.eg.db"  = "mouse",
  "org.Rn.eg.db"  = "rat",       "org.Dm.eg.db"  = "fly",
  "org.Dr.eg.db"  = "zebrafish", "org.Sc.sgd.db" = "yeast",
  "org.Ce.eg.db"  = "celegans",  "org.Bt.eg.db"  = "bovine",
  "org.Cf.eg.db"  = "canine",    "org.Gg.eg.db"  = "chicken"
)

# Result Tab panel
# showResultTabPanel <- function(ns, tabTitle, dataTableId, buttonId, ...) {
#   tabPanel(
#     title = tabTitle,
#     shiny::dataTableOutput(ns(dataTableId)),
#     br(),
#     shinyjs::useShinyjs(),
#     downloadButton(ns(buttonId),
#                    "Download",
#                    class = "btn-primary",
#                    style = "background-color: #d83428; color: white;")
#   )
# }
