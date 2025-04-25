#' Process ID conversion

id_conversion <- function(data, from_id_type, to_id_type, orgDb) {

  # Ensure package availability
  if (!require(clusterProfiler)) {
    BiocManager::install("clusterProfiler")
    library(clusterProfiler)
  }
  # if (!require(org.Hs.eg.db)) {
  #   BiocManager::install("org.Hs.eg.db", force = TRUE)
  #   library(org.Hs.eg.db)
  # }

  # Perform ID conversion
  converted <- clusterProfiler::bitr(
    geneID = data[[tolower(from_id_type)]],
    fromType = from_id_type,
    toType = to_id_type,
    OrgDb = orgDb
  )

  # Remove duplicates according to entrezid
  converted <-
    dplyr::distinct(converted, ENTREZID, .keep_all = TRUE) %>%
    dplyr::rename_with(tolower) %>%
    dplyr::left_join(data, ., by = tolower(from_id_type))

  # Generate code for ID conversion
  conversion_code <- sprintf(
    "if (!require(clusterProfiler)) {\n  install.packages(\"clusterProfiler\")\n  library(clusterProfiler)\n}\nif (!require(%s)) {\n  install.packages(\"%s\")\n  library(%s)\n}\nconverted <- clusterProfiler::bitr(\n  geneID = data[[tolower(%s)]],\n  fromType = \"%s\",\n  toType = c(\"%s\"),\n  OrgDb = %s\n)\nconverted <- \n  dplyr::distinct(converted, ENTREZID, .keep_all = TRUE) %%>%%\n  dplyr::rename_with(tolower) %%>%%\n  dplyr::left_join(data, ., by = tolower(%s))",
    orgDb$packageName, orgDb$packageName, orgDb$packageName,
    from_id_type, from_id_type,
    paste(to_id_type, collapse = "\", \""),
    orgDb$packageName,
    from_id_type
  )

  return(list("converted_id" = converted, "conversion_code" = conversion_code))
}



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
