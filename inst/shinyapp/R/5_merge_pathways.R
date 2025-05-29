#' Merge Pathways UI Module
#'
#' Internal UI for merging enriched pathways.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyBS bsButton bsPopover
#' @noRd

merge_pathways_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "merge_pathways",
          fluidPage(
            titlePanel("Merge Pathways"),
            fluidPage(
              fluidRow(
                column(4,
                       fluidRow(
                         column(8,
                                fileInput(inputId = ns("upload_enriched_pathways"),
                                          label = tags$h4("Upload Enrichment Result"),
                                          accept = ".rda")
                         )
                       ),

                       ### Query metabolite parameter panel ----
                       shinyjs::hidden(
                         div(
                           id = ns("parameter_panel_metabolite"),
                           h4("HMDB Network"),
                           fluidRow(
                             column(6,
                                    numericInput(
                                      ns("p.adjust.cutoff.hmdb"),
                                      "P-adjust cutoff",
                                      value = 0.05,
                                      min = 0,
                                      max = 0.5)
                             ),
                             column(6,
                                    numericInput(
                                      ns("count.cutoff.hmdb"),
                                      "Metabolite count cutoff",
                                      value = 5,
                                      min = 0,
                                      max = 1000)
                             )
                           ),
                           fluidRow(
                             column(6,
                                    selectInput(
                                      ns("measure.method.hmdb"),
                                      "Similarity method",
                                      choices = c("jaccard", "dice", "overlap", "kappa"),
                                      selected = "jaccard")
                             ),
                             column(6,
                                    numericInput(
                                      ns("sim.cutoff.hmdb"),
                                      "Similarity cutoff",
                                      value = 0.5,
                                      min = 0,
                                      max = 1)
                             )
                           ),

                           h4("KEGG Network"),
                           fluidRow(
                             column(6,
                                    numericInput(
                                      ns("p.adjust.cutoff.metkegg"),
                                      "P-adjust cutoff",
                                      value = 0.05,
                                      min = 0,
                                      max = 0.5)
                             ),
                             column(6,
                                    numericInput(
                                      ns("count.cutoff.metkegg"),
                                      "Metabolite count cutoff",
                                      value = 5,
                                      min = 0,
                                      max = 1000)
                             )
                           ),
                           fluidRow(
                             column(6,
                                    selectInput(
                                      ns("measure.method.metkegg"),
                                      "Similarity method",
                                      choices = c("jaccard", "dice", "overlap", "kappa"),
                                      selected = "jaccard")
                             ),
                             column(6,
                                    numericInput(
                                      ns("sim.cutoff.metkegg"),
                                      "Similarity cutoff",
                                      value = 0.5,
                                      min = 0,
                                      max = 1)
                             )
                           )
                         )
                       ),

                       ### Query gene parameter panel ----
                       shinyjs::hidden(
                         div(
                           id = ns("parameter_panel_gene"),

                           checkboxGroupInput(
                             ns("cluster_module_database"),
                             "Available Database",
                             choices = c(
                               "GO" = "go",
                               "KEGG" = "kegg",
                               "Reactome" = "reactome"
                             ),
                             selected = NULL
                           )
                         )
                       ),

                       shinyjs::hidden(
                         div(
                           id = ns("gene_go_network_panel"),

                           h4("GO Network"),

                           # selectInput(
                           #   ns("go_orgdb"),
                           #   label = tags$span(
                           #     "Organism Database",
                           #     shinyBS::bsButton(
                           #       ns("go_orgdb_info"),
                           #       label = "",
                           #       icon = icon("info"),
                           #       style = "info",
                           #       size = "extra_small"
                           #     )
                           #   ),
                           #   choices = c(
                           #     "Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                           #     "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
                           #     "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
                           #     "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db",
                           #     "Fly (org.Dm.eg.db)" = "org.Dm.eg.db",
                           #     "Worm (org.Ce.eg.db)" = "org.Ce.eg.db",
                           #     "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db",
                           #     "E. coli (org.EcK12.eg.db)" = "org.EcK12.eg.db"
                           #   ),
                           #   selected = "org.Hs.eg.db"
                           # ),
                           # shinyBS::bsPopover(
                           #   id = ns("go_orgdb_info"),
                           #   title = "Organism Database",
                           #   content = "OrgDb object or character string naming the OrgDb annotation package used to derive gene–GO mappings, required when database includes \"go\".",
                           #   placement = "right",
                           #   trigger = "hover",
                           #   options = list(container = "body")
                           # ),

                           fluidRow(
                             column(6,
                                    numericInput(
                                      ns("p.adjust.cutoff.go"),
                                      "P-adjust cutoff",
                                      value = 0.05,
                                      min = 0,
                                      max = 0.5)
                             ),
                             column(6,
                                    numericInput(
                                      ns("count.cutoff.go"),
                                      "Gene count cutoff",
                                      value = 5,
                                      min = 0,
                                      max = 1000
                                    )
                             )
                           ),
                           fluidRow(
                             column(6,
                                    selectInput(
                                      ns("measure.method.go"),
                                      "Similarity method",
                                      choices = c("Sim_XGraSM_2013", "Sim_Wang_2007", "Sim_Lin_1998",
                                                  "Sim_Resnik_1999", "Sim_FaITH_2010", "Sim_Relevance_2006",
                                                  "Sim_SimIC_2010", "Sim_EISI_2015", "Sim_AIC_2014",
                                                  "Sim_Zhang_2006", "Sim_universal", "Sim_GOGO_2018",
                                                  "Sim_Rada_1989", "Sim_Resnik_edge_2005", "Sim_Leocock_1998",
                                                  "Sim_WP_1994", "Sim_Slimani_2006", "Sim_Shenoy_2012",
                                                  "Sim_Pekar_2002", "Sim_Stojanovic_2001", "Sim_Wang_edge_2012",
                                                  "Sim_Zhong_2002", "Sim_AlMubaid_2006", "Sim_Li_2003",
                                                  "Sim_RSS_2013", "Sim_HRSS_2013", "Sim_Shen_2010",
                                                  "Sim_SSDD_2013", "Sim_Jiang_1997", "Sim_Kappa", "Sim_Jaccard",
                                                  "Sim_Dice",  "Sim_Overlap", "Sim_Ancestor"),
                                      selected = "Sim_XGraSM_2013"
                                    )
                             ),
                             column(6,
                                    numericInput(
                                      ns("sim.cutoff.go"),
                                      "Similarity cutoff",
                                      value = 0.5,
                                      min = 0,
                                      max = 1)
                             )
                           )
                         )
                       ),

                       shinyjs::hidden(
                         div(
                           id = ns("gene_kegg_network_panel"),

                           h4("KEGG Network"),
                           fluidRow(
                             column(6,
                                    numericInput(
                                      ns("p.adjust.cutoff.kegg"),
                                      "P-adjust cutoff",
                                      value = 0.05,
                                      min = 0,
                                      max = 0.5)
                             ),
                             column(6,
                                    numericInput(
                                      ns("count.cutoff.kegg"),
                                      "Gene count cutoff",
                                      value = 5,
                                      min = 0,
                                      max = 1000)
                             )
                           ),
                           fluidRow(
                             column(6,
                                    selectInput(
                                      ns("measure.method.kegg"),
                                      "Similarity method",
                                      choices = c("jaccard", "dice", "overlap", "kappa"),
                                      selected = "jaccard")
                             ),
                             column(6,
                                    numericInput(
                                      ns("sim.cutoff.kegg"),
                                      "Similarity cutoff",
                                      value = 0.5,
                                      min = 0,
                                      max = 1)
                             )
                           )
                         )
                       ),

                       shinyjs::hidden(
                         div(
                           id = ns("gene_reactome_network_panel"),

                           h4("Reactome Network"),
                           fluidRow(
                             column(6,
                                    numericInput(
                                      ns("p.adjust.cutoff.reactome"),
                                      "P-adjust cutoff",
                                      value = 0.05,
                                      min = 0,
                                      max = 0.5)
                             ),
                             column(6,
                                    numericInput(
                                      ns("count.cutoff.reactome"),
                                      "Gene count cutoff",
                                      value = 5,
                                      min = 0,
                                      max = 1000)
                             )
                           ),
                           fluidRow(
                             column(6,
                                    selectInput(
                                      ns("measure.method.reactome"),
                                      "Similarity method",
                                      choices = c("jaccard", "dice", "overlap", "kappa"),
                                      selected = "jaccard")
                             ),
                             column(6,
                                    numericInput(
                                      ns("sim.cutoff.reactome"),
                                      "Similarity cutoff",
                                      value = 0.5,
                                      min = 0,
                                      max = 1)
                             )
                           )
                         )
                       ),

                       actionButton(
                         ns("submit_merge_pathways"),
                         "Submit",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),

                       actionButton(
                         ns("go2merge_modules"),
                         "Next",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),

                       actionButton(
                         ns("show_merge_pathways_code"),
                         "Code",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),

                       style = "border-right: 1px solid #ddd; padding-right: 20px;"
                ),

                ### Result display ----
                column(8,
                       tabsetPanel(
                         tabPanel("Table",
                                  shinyjs::hidden(
                                    div(
                                      id = ns("table_panel_gene"),
                                      tabsetPanel(
                                        tabPanel(
                                          title = "GO",
                                          shiny::dataTableOutput(ns("merged_pathway_go")),
                                          br(),
                                          shinyjs::useShinyjs(),
                                          downloadButton(ns("download_merged_pathway_go"),
                                                         "Download",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                        ),
                                        tabPanel(
                                          title = "KEGG",
                                          shiny::dataTableOutput(ns("merged_pathway_kegg")),
                                          br(),
                                          shinyjs::useShinyjs(),
                                          downloadButton(ns("download_merged_pathway_kegg"),
                                                         "Download",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                        ),
                                        tabPanel(
                                          title = "Reactome",
                                          shiny::dataTableOutput(ns("merged_pathway_reactome")),
                                          br(),
                                          shinyjs::useShinyjs(),
                                          downloadButton(ns("download_merged_pathway_reactome"),
                                                         "Download",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                        )
                                      )
                                    )
                                  ),
                                  shinyjs::hidden(
                                    div(
                                      id = ns("table_panel_metabolite"),
                                      tabsetPanel(
                                        tabPanel(
                                          title = "HMDB",
                                          shiny::dataTableOutput(ns("merged_pathway_hmdb")),
                                          br(),
                                          shinyjs::useShinyjs(),
                                          downloadButton(ns("download_merged_pathway_hmdb"),
                                                         "Download",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                        ),
                                        tabPanel(
                                          title = "KEGG",
                                          shiny::dataTableOutput(ns("merged_pathway_metkegg")),
                                          br(),
                                          shinyjs::useShinyjs(),
                                          downloadButton(ns("download_merged_pathway_metkegg"),
                                                         "Download",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                        )
                                      )
                                    )
                                  )
                                  ),
                         tabPanel(
                           title = "Data visualization",
                           shinyjs::hidden(
                             div(
                               id = ns("plot_panel_gene"),
                               tabsetPanel(
                                 tabPanel(
                                   title = "GO",
                                   shiny::plotOutput(ns("enirched_module_go_plot")),
                                   br(),
                                   fluidRow(
                                     column(3,
                                            actionButton(ns("generate_enirched_module_plot_go"),
                                                         "Generate plot",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_go"), "Text", FALSE)
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_all_go"), "Text all", FALSE)
                                     ),
                                     column(3,
                                            numericInput(
                                              ns("enirched_module_plot_degree_cutoff_go"),
                                              "Degree cutoff",
                                              value = 0,
                                              min = 0,
                                              max = 1000)
                                     )
                                   )
                                 ),
                                 tabPanel(
                                   title = "KEGG",
                                   shiny::plotOutput(ns("enirched_module_kegg_plot")),
                                   br(),
                                   fluidRow(
                                     column(3,
                                            actionButton(ns("generate_enirched_module_plot_kegg"),
                                                         "Generate plot",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_kegg"), "Text", FALSE)
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_all_kegg"), "Text all", FALSE)
                                     ),
                                     column(3,
                                            numericInput(
                                              ns("enirched_module_plot_degree_cutoff_kegg"),
                                              "Degree cutoff",
                                              value = 0,
                                              min = 0,
                                              max = 1000
                                            )
                                     )
                                   )
                                 ),
                                 tabPanel(
                                   title = "Reactome",
                                   shiny::plotOutput(ns("enirched_module_reactome_plot")),
                                   br(),
                                   fluidRow(
                                     column(3,
                                            actionButton(ns("generate_enirched_module_plot_reactome"),
                                                         "Generate plot",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_reactome"), "Text", FALSE)
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_all_reactome"), "Text all", FALSE)
                                     ),
                                     column(3,
                                            numericInput(
                                              ns("enirched_module_plot_degree_cutoff_reactome"),
                                              "Degree cutoff",
                                              value = 0,
                                              min = 0,
                                              max = 1000)
                                     )
                                   )
                                 )
                               )
                             )
                           ),
                           shinyjs::hidden(
                             div(
                               id = ns("plot_panel_metabolite"),
                               tabsetPanel(
                                 tabPanel(
                                   title = "HMDB",
                                   shiny::plotOutput(ns("enirched_module_hmdb_plot")),
                                   br(),
                                   fluidRow(
                                     column(3,
                                            actionButton(ns("generate_enirched_module_plot_hmdb"),
                                                         "Generate plot",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_hmdb"), "Text", FALSE)
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_all_hmdb"), "Text all", FALSE)
                                     ),
                                     column(3,
                                            numericInput(
                                              ns("enirched_module_plot_degree_cutoff_hmdb"),
                                              "Degree cutoff",
                                              value = 0,
                                              min = 0,
                                              max = 1000)
                                     )
                                   )
                                 ),
                                 tabPanel(
                                   title = "KEGG",
                                   shiny::plotOutput(ns("enirched_module_metkegg_plot")),
                                   br(),
                                   fluidRow(
                                     column(3,
                                            actionButton(ns("generate_enirched_module_plot_metkegg"),
                                                         "Generate plot",
                                                         class = "btn-primary",
                                                         style = "background-color: #d83428; color: white;")
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_metkegg"), "Text", FALSE)
                                     ),
                                     column(3,
                                            checkboxInput(ns("enirched_module_plot_text_all_metkegg"), "Text all", FALSE)
                                     ),
                                     column(3,
                                            numericInput(
                                              ns("enirched_module_plot_degree_cutoff_metkegg"),
                                              "Degree cutoff",
                                              value = 0,
                                              min = 0,
                                              max = 1000
                                            )
                                     )
                                   )
                                 )
                               )
                             )
                           )
                         ),
                         tabPanel(
                           title = "R object",
                           verbatimTextOutput(ns("enriched_modules_object")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_modules_object"),
                                          label = tags$span("Download",
                                                            shinyBS::bsButton(ns("download_enriched_modules_object_info"),
                                                                              label = "",
                                                                              icon = icon("info"),
                                                                              style = "info",
                                                                              size = "extra_small")),
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;"),
                           shinyBS::bsPopover(
                             id = ns("download_enriched_modules_object_info"),
                             title = "",
                             content = "You can download the functional module file for data visualization.",
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")
                           )
                         )
                       )
                )
              )
            )
          )
          )
}

#' Merge Pathways Server Module
#'
#' Internal server logic for merging enriched pathways into database-specific modules.
#' This module handles the server-side functionality for pathway clustering and
#' merging based on similarity measures, supporting both gene and metabolite data.
#'
#' @param id Character string. Module namespace identifier.
#' @param enriched_pathways ReactiveValues object containing:
#'   \describe{
#'     \item{enriched_pathways_res}{The enriched pathways result object}
#'     \item{query_type}{Character. Type of query data ("gene" or "metabolite")}
#'     \item{organism}{Character. Organism information for analysis}
#'     \item{available_db}{Character vector. Available databases for analysis}
#'   }
#' @param enriched_modules ReactiveVal object to store the merged pathway results.
#' @param tab_switch Function to handle navigation between application tabs.
#'
#' @section Supported Databases:
#' \describe{
#'   \item{Gene queries}{GO (Gene Ontology), KEGG, Reactome}
#'   \item{Metabolite queries}{HMDB (Human Metabolome Database), KEGG}
#' }
#'
#'
#' @import shiny
#' @importFrom shinyjs toggleState useShinyjs toggleElement enable disable
#'
#' @note This function requires the clusterProfiler and ReactomePA packages
#'   to be installed and loaded. For GO analysis, appropriate organism databases
#'   (e.g., org.Hs.eg.db) must also be available.
#'
#' @examples
#' \dontrun{
#' # Used within the main server function
#' merge_pathways_server("merge_pathways_tab",
#'                       enriched_pathways = enriched_pathways,
#'                       enriched_modules = enriched_modules,
#'                       tab_switch = tab_switch)
#' }
#'
#' @noRd

merge_pathways_server <- function(id, enriched_pathways, enriched_modules, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {

      ns <- session$ns

      observeEvent(input$upload_enriched_pathways, {
        if (!is.null(input$upload_enriched_pathways$datapath)) {
          message("Loading data")
          tempEnv <- new.env()
          load(input$upload_enriched_pathways$datapath,
               envir = tempEnv)

          names <- ls(tempEnv)

          if (length(names) == 1) {
            enriched_pathways$enriched_pathways_res <- get(names[1], envir = tempEnv)
            if ("enrich_pathway" %in% names(enriched_pathways$enriched_pathways_res@process_info)) {
              enriched_pathways$available_db <- enriched_pathways$enriched_pathways_res@process_info$enrich_pathway@parameter$database
              enriched_pathways$query_type <- enriched_pathways$enriched_pathways_res@process_info$enrich_pathway@parameter$query_type
              if (enriched_pathways$query_type == "gene") {
                enriched_pathways$organism <- enriched_pathways$enriched_pathways_res@process_info$enrich_pathway@parameter$go.orgdb
              } else {
                enriched_pathways$organism <- enriched_pathways$enriched_pathways_res@process_info$enrich_pathway@parameter$met_organism
              }
            } else {
              enriched_pathways$query_type <- enriched_pathways$enriched_pathways_res@process_info$do_gsea@parameter$query_type
              enriched_pathways$available_db <- enriched_pathways$enriched_pathways_res@process_info$do_gsea@parameter$database
              enriched_pathways$organism <- enriched_pathways$enriched_pathways_res@process_info$do_gsea@parameter$go.orgdb
            }
          } else {
            message("The .rda file does not contain exactly one object.")
            showModal(
              modalDialog(
                title = "Error",
                "The uploaded file should contain exactly one object.",
                easyClose = TRUE,
                footer = modalButton("Close")
              )
            )
          }
        }
      })

      # query_type <- reactive({enriched_pathways$query_type})

      observe({
        req(enriched_pathways$query_type)
        ## For gene
        shinyjs::toggleElement(
          id = "parameter_panel_gene",
          condition = enriched_pathways$query_type == "gene"
        )
        shinyjs::toggleElement(
          id = "table_panel_gene",
          condition = enriched_pathways$query_type == "gene"
        )
        shinyjs::toggleElement(
          id = "plot_panel_gene",
          condition = enriched_pathways$query_type == "gene"
        )

        ## For metabolite
        shinyjs::toggleElement(
          id = "parameter_panel_metabolite",
          condition = enriched_pathways$query_type == "metabolite"
        )
        shinyjs::toggleElement(
          id = "table_panel_metabolite",
          condition = enriched_pathways$query_type == "metabolite"
        )
        shinyjs::toggleElement(
          id = "plot_panel_metabolite",
          condition = enriched_pathways$query_type == "metabolite"
        )
      })

      go_orgdb <- reactiveVal(NULL)
      db_choices <- c("GO" = "go", "KEGG" = "kegg", "Reactome" = "reactome")
      observeEvent(enriched_pathways$available_db,
                   {
                     # print(enriched_pathways$available_db)
                     if (enriched_pathways$query_type == "gene") {
                       updateCheckboxGroupInput(
                         session, "cluster_module_database",
                         choices  = db_choices[db_choices %in% enriched_pathways$available_db],
                         selected = enriched_pathways$available_db
                       )
                     }
                   })

      observe(
        {
          req(input$cluster_module_database)
          # Define database-panel pairs
          db_panels <- list(
            "go" = "gene_go_network_panel",
            "kegg" = "gene_kegg_network_panel",
            "reactome" = "gene_reactome_network_panel"
          )

          # Loop through each database-panel pair and toggle accordingly
          for (db in names(db_panels)) {
            shinyjs::toggleElement(
              id = db_panels[[db]],
              condition = db %in% input$cluster_module_database
            )
          }

          if ("go" %in% input$cluster_module_database) {
            go_orgdb(enriched_pathways$organism)
          }
        }
      )

      merge_pathways_code <-
        reactiveVal()

      observeEvent(
        input$submit_merge_pathways,
        {
          if (is.null(enriched_pathways$enriched_pathways_res) || length(enriched_pathways$enriched_pathways_res) == 0) {
            showModal(
              modalDialog(
                title = "Warning",
                "No enriched pathways data available. Please 'Enrich pathways' first.",
                easyClose = TRUE,
                footer = modalButton("Close")
              )
            )
          } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch({
              library(clusterProfiler)
              library(ReactomePA)

              if (!is.null(go_orgdb())) {
                # Validate input format
                if (!grepl("^org\\.[A-Za-z]+\\..+\\.db$", go_orgdb())) {
                  stop("Invalid OrgDb package name. Expected format: org.XX.eg.db")
                }
                # Check if the package is installed
                if (!requireNamespace(go_orgdb(), quietly = TRUE)) {
                  stop(paste("Package", go_orgdb(), "is not installed. Please install it using BiocManager::install('", go_orgdb(), "')"))
                }
                # Load the package
                requireNamespace(go_orgdb())
                # Get the OrgDb object
                org_db_obj <- get(go_orgdb())
              } else {
                org_db_obj <- NULL
              }

              result <-
                merge_pathways(
                  object = enriched_pathways$enriched_pathways_res,
                  database = input$cluster_module_database,
                  go.orgdb = org_db_obj,
                  p.adjust.cutoff.go = input$p.adjust.cutoff.go,
                  p.adjust.cutoff.kegg = input$p.adjust.cutoff.kegg,
                  p.adjust.cutoff.reactome = input$p.adjust.cutoff.reactome,
                  p.adjust.cutoff.hmdb = input$p.adjust.cutoff.hmdb,
                  p.adjust.cutoff.metkegg = input$p.adjust.cutoff.metkegg,
                  count.cutoff.go = input$count.cutoff.go,
                  count.cutoff.kegg = input$count.cutoff.kegg,
                  count.cutoff.reactome = input$count.cutoff.reactome,
                  count.cutoff.hmdb = input$count.cutoff.hmdb,
                  count.cutoff.metkegg = input$count.cutoff.metkegg,
                  sim.cutoff.go = input$sim.cutoff.go,
                  sim.cutoff.kegg = input$sim.cutoff.kegg,
                  sim.cutoff.reactome = input$sim.cutoff.reactome,
                  sim.cutoff.hmdb = input$sim.cutoff.hmdb,
                  sim.cutoff.metkegg = input$sim.cutoff.metkegg,
                  measure.method.go = input$measure.method.go,
                  measure.method.kegg = input$measure.method.kegg,
                  measure.method.reactome = input$measure.method.reactome,
                  measure.method.hmdb = input$measure.method.hmdb,
                  measure.method.metkegg = input$measure.method.metkegg,
                  path = "result",
                  save_to_local = FALSE
                )
            },
            error = function(e) {
              showModal(modalDialog(
                title = "Error",
                paste("Details:", e$message),
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
            })
          })

          enriched_modules(result)

          # shinyjs::hide("loading")

          ## Save code ====
          if (enriched_pathways$query_type == "gene") {

            go_params <- ""
            if ("go" %in% input$cluster_module_database) {
              go_params <- sprintf(
                '
                  p.adjust.cutoff.go = %s,
                  count.cutoff.go = %s,
                  sim.cutoff.go = %s,
                  measure.method.go = %s,
                  go.orgdb = %s,
                  ',
                input$p.adjust.cutoff.go,
                input$count.cutoff.go,
                input$sim.cutoff.go,
                paste0('"', input$measure.method.go, '"'),
                go_orgdb()
              )}

            kegg_params <- ""
            if ("kegg" %in% input$cluster_module_database) {
              kegg_params <- sprintf(
                ' p.adjust.cutoff.kegg = %s,
                  count.cutoff.kegg = %s,
                  sim.cutoff.kegg = %s,
                  measure.method.kegg = %s,
                 ',
                input$p.adjust.cutoff.kegg,
                input$count.cutoff.kegg,
                input$sim.cutoff.kegg,
                paste0('"', input$measure.method.kegg, '"')
              )}

            reactome_params <- ""
            if ("reactome" %in% input$cluster_module_database) {
              reactome_params <- sprintf(
                ' p.adjust.cutoff.reactome = %s,
                  count.cutoff.reactome = %s,
                  sim.cutoff.reactome = %s,
                  measure.method.reactome = %s,',
                input$p.adjust.cutoff.reactome,
                input$count.cutoff.reactome,
                input$sim.cutoff.reactome,
                paste0('"', input$measure.method.reactome, '"')
              )}

            db_vector <- paste0('c("', paste(input$cluster_module_database, collapse = '", "'), '")')

            merge_pathways_code_str <- sprintf(
              '
              enriched_modules <-
                merge_pathways(
                  object = enriched_pathways,
                  database = %s,%s%s%s
                  save_to_local = FALSE
                  )
              ',
              db_vector,
              go_params,
              kegg_params,
              reactome_params
            )

          } else if (enriched_pathways$query_type == "metabolite") {

            merge_pathways_code_str <- sprintf(
              '
              enriched_modules <-
                merge_pathways(
                object = enriched_pathways,
                database = c("hmdb", "kegg"),
                p.adjust.cutoff.hmdb = %s,
                p.adjust.cutoff.metkegg = %s,
                count.cutoff.hmdb = %s,
                count.cutoff.metkegg = %s,
                sim.cutoff.hmdb = %s,
                sim.cutoff.metkegg = %s,
                measure.method.hmdb = %s,
                measure.method.metkegg = %s
                )
              ',
              input$p.adjust.cutoff.hmdb,
              input$p.adjust.cutoff.metkegg,
              input$count.cutoff.hmdb,
              input$count.cutoff.metkegg,
              input$sim.cutoff.hmdb,
              input$sim.cutoff.metkegg,
              paste0('"', input$measure.method.hmdb, '"'),
              paste0('"', input$measure.method.metkegg, '"')
            )

          }

          merge_pathways_code(merge_pathways_code_str)
        }
      })

      ## Show object =====
      output$enriched_modules_object <-
        renderText({
          enriched_modules <- enriched_modules()
          captured_output1 <- capture.output(enriched_modules,
                                             type = "message")
          captured_output2 <- capture.output(enriched_modules,
                                             type = "output")
          captured_output <-
            c(captured_output1,
              captured_output2)
          paste(captured_output, collapse = "\n")
        })

      ## Show table ====
      output$merged_pathway_go <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_modules()@merged_pathway_go$module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$merged_pathway_kegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_modules()@merged_pathway_kegg$module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$merged_pathway_reactome <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_modules()@merged_pathway_reactome$module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      ### For metabolite
      output$merged_pathway_hmdb <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_modules()@merged_pathway_hmdb$module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))
      output$merged_pathway_metkegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_modules()@merged_pathway_metkegg$module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      ## Download results ====
      output$download_merged_pathway_go <-
        shiny::downloadHandler(
          filename = function() {
            "merged_pathway_go.csv"
          },
          content = function(file) {
            write.csv(enriched_modules()@merged_pathway_go$module_result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_modules()) ||
                length(enriched_modules()) == 0) {
              shinyjs::disable("download_merged_pathway_go")
            } else {
              if (length(enriched_modules()@merged_pathway_go) == 0) {
                shinyjs::disable("download_merged_pathway_go")
              } else{
                shinyjs::enable("download_merged_pathway_go")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_merged_pathway_go")
          }
        )
      })

      output$download_merged_pathway_kegg <-
        shiny::downloadHandler(
          filename = function() {
            "merged_pathway_kegg.csv"
          },
          content = function(file) {
            write.csv(enriched_modules()@merged_pathway_kegg$module_result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_modules()) ||
                length(enriched_modules()) == 0) {
              shinyjs::disable("download_merged_pathway_kegg")
            } else {
              if (length(enriched_modules()@merged_pathway_kegg) == 0) {
                shinyjs::disable("download_merged_pathway_kegg")
              } else{
                shinyjs::enable("download_merged_pathway_kegg")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_merged_pathway_kegg")
          }
        )
      })

      output$download_merged_pathway_reactome <-
        shiny::downloadHandler(
          filename = function() {
            "merged_pathway_reactome.csv"
          },
          content = function(file) {
            write.csv(enriched_modules()@merged_pathway_reactome$module_result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_modules()) ||
                length(enriched_modules()) == 0) {
              shinyjs::disable("download_merged_pathway_reactome")
            } else {
              if (length(enriched_modules()@merged_pathway_reactome) == 0) {
                shinyjs::disable("download_merged_pathway_reactome")
              } else{
                shinyjs::enable("download_merged_pathway_reactome")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_merged_pathway_reactome")
          }
        )
      })

      output$download_merged_pathway_hmdb <-
        shiny::downloadHandler(
          filename = function() {
            "merged_pathway_hmdb.csv"
          },
          content = function(file) {
            write.csv(enriched_modules()@merged_pathway_hmdb$module_result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_modules()) ||
                length(enriched_modules()) == 0) {
              shinyjs::disable("download_merged_pathway_hmdb")
            } else {
              if (length(enriched_modules()@merged_pathway_hmdb) == 0) {
                shinyjs::disable("download_merged_pathway_hmdb")
              } else{
                shinyjs::enable("download_merged_pathway_hmdb")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_merged_pathway_hmdb")
          }
        )
      })

      output$download_merged_pathway_metkegg <-
        shiny::downloadHandler(
          filename = function() {
            "merged_pathway_metkegg.csv"
          },
          content = function(file) {
            write.csv(enriched_modules()@merged_pathway_metkegg$module_result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_modules()) ||
                length(enriched_modules()) == 0) {
              shinyjs::disable("download_merged_pathway_metkegg")
            } else {
              if (length(enriched_modules()@merged_pathway_metkegg) == 0) {
                shinyjs::disable("download_merged_pathway_metkegg")
              } else{
                shinyjs::enable("download_merged_pathway_metkegg")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_merged_pathway_metkegg")
          }
        )
      })


      output$download_enriched_modules_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_modules.rda"
          },
          content = function(file) {
            em <- enriched_modules()
            save(em, file = file)
          }
        )

      observe({
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          shinyjs::disable("download_enriched_modules_object")
        } else {
          shinyjs::enable("download_enriched_modules_object")
        }
      })

      ## Data visualization ====
      # GO Plot generation logic
      enirched_module_go_plot <-
        reactiveVal()
      observeEvent(input$generate_enirched_module_plot_go, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Merge pathways' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_modules(),
                  level = "module",
                  database = "go",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_go,
                  text = input$enirched_module_plot_text_go,
                  text_all = input$enirched_module_plot_text_all_go
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              }
            )

          })

          enirched_module_go_plot(plot)
          # shinyjs::hide("loading")
        }
      })

      output$enirched_module_go_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_module_go_plot(),
            error = function(e)
              NULL
          ))
        })

      # kegg Plot generation logic
      enirched_module_kegg_plot <-
        reactiveVal()
      observeEvent(input$generate_enirched_module_plot_kegg, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Merge pathways' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_modules(),
                  level = "module",
                  database = "kegg",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_kegg,
                  text = input$enirched_module_plot_text_kegg,
                  text_all = input$enirched_module_plot_text_all_kegg
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              }
            )
          })

          enirched_module_kegg_plot(plot)

          # shinyjs::hide("loading")
        }
      })

      output$enirched_module_kegg_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_module_kegg_plot(),
            error = function(e)
              NULL
          ))
        })

      # reactome Plot generation logic
      enirched_module_reactome_plot <-
        reactiveVal()
      observeEvent(input$generate_enirched_module_plot_reactome, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Merge pathways' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_modules(),
                  level = "module",
                  database = "reactome",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_reactome,
                  text = input$enirched_module_plot_text_reactome,
                  text_all = input$enirched_module_plot_text_all_reactome
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              }
            )

          })

          enirched_module_reactome_plot(plot)
          # shinyjs::hide("loading")
        }
      })

      output$enirched_module_reactome_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_module_reactome_plot(),
            error = function(e)
              NULL
          ))
        })

      # hmdb Plot generation logic
      enirched_module_hmdb_plot <-
        reactiveVal()
      observeEvent(input$generate_enirched_module_plot_hmdb, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Merge pathways' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_modules(),
                  level = "module",
                  database = "hmdb",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_hmdb,
                  text = input$enirched_module_plot_text_hmdb,
                  text_all = input$enirched_module_plot_text_all_hmdb
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              }
            )
          })

          enirched_module_hmdb_plot(plot)

          # shinyjs::hide("loading")
        }
      })

      output$enirched_module_hmdb_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_module_hmdb_plot(),
            error = function(e)
              NULL
          ))
        })

      # metabolite KEGG Plot generation logic
      enirched_module_metkegg_plot <-
        reactiveVal()
      observeEvent(input$generate_enirched_module_plot_metkegg, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Merge pathways' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_modules(),
                  level = "module",
                  database = "metkegg",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_metkegg,
                  text = input$enirched_module_plot_text_metkegg,
                  text_all = input$enirched_module_plot_text_all_metkegg
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              }
            )
          })

          enirched_module_metkegg_plot(plot)

          # shinyjs::hide("loading")
        }
      })

      output$enirched_module_metkegg_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_module_metkegg_plot(),
            error = function(e)
              NULL
          ))
        })


      ####show code
      observeEvent(input$show_merge_pathways_code, {
        if (is.null(merge_pathways_code()) ||
            length(merge_pathways_code()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No available code",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else{
          code_content <-
            merge_pathways_code()
          code_content <-
            paste(code_content, collapse = "\n")
          showModal(modalDialog(
            title = "Code",
            tags$pre(code_content),
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        }
      })

      ##Go to merge modules tab
      ###if there is not enriched_modules, show a warning message
      observeEvent(input$go2merge_modules, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "Please merge pathways first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("merge_modules")
        }
      })
    }
  )
}
