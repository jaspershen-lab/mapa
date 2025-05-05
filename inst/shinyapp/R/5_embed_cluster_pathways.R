embed_cluster_pathways_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "embed_cluster_pathways",
    fluidPage(
      titlePanel("Embed and Cluster Pathways"),
      fluidRow(
        column(4,
               # fileInput(
               #   ns("variable_info"),
               #   tags$h4("Choose marker information"),
               #   accept = c(
               #     "text/csv",
               #     "text/comma-separated-values,text/plain",
               #     ".csv",
               #     ".xlsx",
               #     "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
               #     "application/vnd.ms-excel"
               #   )
               # ),
               # checkboxInput(
               #   ns("use_example"),
               #   "Use example",
               #   FALSE
               # ),
               # shinyjs::hidden(
               #   div(id = ns("example_panel"),
               #       radioButtons(
               #         ns("example_choice"),
               #         tags$h4("Example dataset"),
               #         choices = c(
               #           "Pathway Enrichment Example" = "example_enrich_pathway",
               #           "GSEA Example" = "example_gsea"
               #         )
               #       )
               #   )
               # ),

               radioButtons(
                 ns("id_type"),
                 tags$h4("Input ID type"),
                 choices = list(
                   "ENSEMBL" = "ensembl",
                   "UniProt" = "uniprot",
                   "EntrezID" = "entrezid"
                 ),
                 selected = "ensembl"
               ),

               # textInput(
               #   ns("orgdb"),
               #   tags$h4("Organism Database"),
               #   value = "org.Hs.eg.db"
               # ),
               # helpText("Enter the name of an OrgDb package that is installed on your system.",
               #          "Common examples: org.Hs.eg.db (Human), org.Mm.eg.db (Mouse), org.Rn.eg.db (Rat)",
               #          "For the current list of OrgDb packages, visit: ",
               #          tags$a(
               #            href = "https://bioconductor.org/packages/release/BiocViews.html#___OrgDb",
               #            "Bioconductor OrgDb packages",
               #            target = "_blank"
               #          )),

               actionButton(
                 ns("map_id"),
                 "Submit",
                 class = "btn-primary",
                 style = "background-color: #d83428; color: white;"
               ),
               actionButton(
                 ns("go2enrich_pathways"),
                 "Next",
                 class = "btn-primary",
                 style = "background-color: #d83428; color: white;"
               ),
               actionButton(
                 ns("show_conversion_code"),
                 "Code",
                 class = "btn-primary",
                 style = "background-color: #d83428; color: white;"
               ),
               style = "border-right: 1px solid #ddd; padding-right: 20px;"
        ),
        column(8,
               tabsetPanel(
                 tabPanel(
                   title = "Marker information",
                   shiny::dataTableOutput(ns("show_variable_info")),
                   br(),
                   shinyjs::useShinyjs(),
                   downloadButton(
                     ns("download_variable_info"),
                     "Download",
                     class = "btn-primary",
                     style = "background-color: #d83428; color: white;"
                   )
                 )
               )
        )
      )
    )
  )
}

embed_cluster_pathways_server <- function(id, enriched_pathways = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ## Step1: Embedding =====

      ## Step2: Clustering ====
    }
  )
}
