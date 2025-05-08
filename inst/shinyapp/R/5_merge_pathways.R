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
                                          label = tags$span("Upload Enrichment Analysis Result"),
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
                           h4("GO Network"),
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
                           ),

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
                           ),

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
#' Internal server logic for merging enriched pathways.
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @param id Module id.
#' @param enriched_pathways Reactive value containing enriched pathways.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleState useShinyjs
#' @importFrom clusterProfiler merge_pathways
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA enrichPathway
#' @noRd

merge_pathways_server <- function(id, enriched_pathways = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {

      observeEvent(input$upload_enriched_pathways, {
        if (!is.null(input$upload_enriched_pathways$datapath)) {
          message("Loading data")
          tempEnv <- new.env()
          load(input$upload_enriched_pathways$datapath,
               envir = tempEnv)

          names <- ls(tempEnv)

          if (length(names) == 1) {
            enriched_pathways(get(names[1], envir = tempEnv))
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

      query_type <- reactive({
        req(enriched_pathways())
        enriched_pathways()@process_info$enrich_pathway@parameter$query_type
      })

      observe({
        req(query_type())
        ## For gene
        shinyjs::toggleElement(
          id = "parameter_panel_gene",
          condition = query_type() == "gene"
        )
        shinyjs::toggleElement(
          id = "table_panel_gene",
          condition = query_type() == "gene"
        )
        shinyjs::toggleElement(
          id = "plot_panel_gene",
          condition = query_type() == "gene"
        )

        ## For metabolite
        shinyjs::toggleElement(
          id = "parameter_panel_metabolite",
          condition = query_type() == "metabolite"
        )
        shinyjs::toggleElement(
          id = "table_panel_metabolite",
          condition = query_type() == "metabolite"
        )
        shinyjs::toggleElement(
          id = "plot_panel_metabolite",
          condition = query_type() == "metabolite"
        )
      })

      ## Define enriched_modules as a reactive value
      enriched_modules <-
        reactiveVal()

      merge_pathways_code <-
        reactiveVal()

      observeEvent(input$submit_merge_pathways, {
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
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
              library(org.Hs.eg.db)
              library(ReactomePA)
              result <-
                merge_pathways(
                  object = enriched_pathways(),
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
          if (query_type() == "gene") {

            merge_pathways_code_str <- sprintf(
              '
              enriched_modules <-
                merge_pathways(
                  object               = enriched_pathways,
                  p.adjust.cutoff.go   = %s,
                  p.adjust.cutoff.kegg = %s,
                  p.adjust.cutoff.reactome = %s,
                  count.cutoff.go      = %s,
                  count.cutoff.kegg    = %s,
                  count.cutoff.reactome = %s,
                  sim.cutoff.go        = %s,
                  sim.cutoff.kegg      = %s,
                  sim.cutoff.reactome  = %s,
                  measure.method.go    = %s,
                  measure.method.kegg  = %s,
                  measure.method.reactome = %s
                  )
              ',
              input$p.adjust.cutoff.go,
              input$p.adjust.cutoff.kegg,
              input$p.adjust.cutoff.reactome,
              input$count.cutoff.go,
              input$count.cutoff.kegg,
              input$count.cutoff.reactome,
              input$sim.cutoff.go,
              input$sim.cutoff.kegg,
              input$sim.cutoff.reactome,
              paste0('"', input$measure.method.go, '"'),
              paste0('"', input$measure.method.kegg, '"'),
              paste0('"', input$measure.method.reactome, '"')
            )

          } else if (query_type() == "metabolite") {

            merge_pathways_code_str <- sprintf(
              '
              enriched_modules <-
                merge_pathways(
                object                  = enriched_pathways,
                p.adjust.cutoff.hmdb    = %s,
                p.adjust.cutoff.metkegg = %s,
                count.cutoff.hmdb       = %s,
                count.cutoff.metkegg    = %s,
                sim.cutoff.hmdb         = %s,
                sim.cutoff.metkegg      = %s,
                measure.method.hmdb     = %s,
                measure.method.metkegg  = %s
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
                    "Please check your input parameters.",
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
                    "Please check your input parameters.",
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
                    "Please check your input parameters.",
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
                    "Please check your input parameters.",
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
                  database = "kegg",
                  degree_cutoff = input$enirched_module_plot_degree_cutoff_metkegg,
                  text = input$enirched_module_plot_text_metkegg,
                  text_all = input$enirched_module_plot_text_all_metkegg
                ),
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    "Please check your input parameters.",
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

      return(enriched_modules)
    }
  )
}
