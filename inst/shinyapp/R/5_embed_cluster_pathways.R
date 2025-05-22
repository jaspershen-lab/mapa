#' Embed and Cluster Pathways UI Module
#'
#' Internal UI for embedding pathway descriptions, clustering them into
#' functional modules, and displaying the resulting tables, plots, and
#' downloadable objects.
#'
#' @param id Module id.
#'
#' @import shiny
#' @importFrom shinyjs hidden toggleElement toggleState useShinyjs disable enable
#' @noRd

embed_cluster_pathways_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "embed_cluster_pathways",
    fluidPage(
      titlePanel("Embed and Cluster Pathways"),
      fluidRow(
        column(4,
               ### Upload enrichment analysis result ====
               fluidRow(
                 column(8,
                        fileInput(inputId = ns("upload_enriched_pathways"),
                                  label = tags$h4("Upload Enrichment Result"),
                                  accept = ".rda")
                 )
               ),

               ### Text embedding model ====
               tags$h4("Pathway Biotext Embedding"),
               fluidRow(
                 column(4,
                        selectInput(
                          ns("api_provider"),
                          tags$span(
                            class = "normal-label",
                            "API provider"),
                          choices = c("openai", "gemini"),
                          selected = "openai")
                  ),
                 column(8,
                        textInput(
                          ns("embedding_model"),
                          tags$span(
                            class = "normal-label",
                            "Embedding model"),
                          value = "text-embedding-3-small")
                 )),
               fluidRow(
                 column(
                   12,
                   textInput(ns("api_key"),
                             tags$span(
                               class = "normal-label",
                               "API key"),
                             value = "")
                 )
               ),

               ### Query metabolite parameter panel ====
               shinyjs::hidden(
                 div(
                   id = ns("parameter_panel_metabolite"),

                   checkboxGroupInput(
                     ns("cluster_module_database"),
                     "Available Database",
                     choices = c(
                       "HMDB" = "hmdb",
                       "KEGG" = "metkegg"
                     ),
                     selected = NULL
                   )
                 )
               ),

               shinyjs::hidden(
                 div(
                   id = ns("metabolite_hmdb_panel"),

                   span(tags$b("HMDB")),
                   fluidRow(
                     column(6,
                            numericInput(
                              ns("p.adjust.cutoff.hmdb"),
                              tags$span(
                                class = "normal-label",
                                "P-adjust cutoff"),
                              value = 0.05,
                              min = 0,
                              max = 0.5)
                     ),
                     column(6,
                            numericInput(
                              ns("count.cutoff.hmdb"),
                              tags$span(
                                class = "normal-label",
                                "Metabolite count cutoff"),
                              value = 5,
                              min = 0,
                              max = 1000)
                     )
                   )
                  )
                 ),

               shinyjs::hidden(
                 div(
                   id = ns("metabolite_metkegg_panel"),

                   span(tags$b("KEGG")),
                   fluidRow(
                     column(6,
                            numericInput(
                              ns("p.adjust.cutoff.metkegg"),
                              tags$span(
                                class = "normal-label",
                                "P-adjust cutoff"),
                              value = 0.05,
                              min = 0,
                              max = 0.5)
                     ),
                     column(6,
                            numericInput(
                              ns("count.cutoff.metkegg"),
                              tags$span(
                                class = "normal-label",
                                "Metabolite count cutoff"),
                              value = 5,
                              min = 0,
                              max = 1000)
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
                  id = ns("gene_go_panel"),

                  span(tags$b("GO")),
                  fluidRow(
                    column(6,
                          numericInput(
                            ns("p.adjust.cutoff.go"),
                            tags$span(
                              class = "normal-label",
                              "P-adjust cutoff"),
                            value = 0.05,
                            min = 0,
                            max = 0.5)
                    ),
                    column(6,
                           numericInput(
                             ns("count.cutoff.go"),
                             tags$span(
                               class = "normal-label",
                               "Gene count cutoff"),
                             value = 5,
                             min = 0,
                             max = 1000
                           )
                    )
                  )
                )
              ),

              shinyjs::hidden(
                div(
                  id = ns("gene_kegg_panel"),

                  span(tags$b("KEGG")),
                  fluidRow(
                    column(6,
                           numericInput(
                             ns("p.adjust.cutoff.kegg"),
                             tags$span(
                               class = "normal-label",
                               "P-adjust cutoff"),
                             value = 0.05,
                             min = 0,
                             max = 0.5)
                    ),
                    column(6,
                           numericInput(
                             ns("count.cutoff.kegg"),
                             tags$span(
                               class = "normal-label",
                               "Gene count cutoff"),
                             value = 5,
                             min = 0,
                             max = 1000
                           )
                    )
                  )
                )
              ),

              shinyjs::hidden(
                div(
                  id = ns("gene_reactome_panel"),

                  span(tags$b("Reactome")),
                  fluidRow(
                    column(6,
                           numericInput(
                             ns("p.adjust.cutoff.reactome"),
                             tags$span(
                               class = "normal-label",
                               "P-adjust cutoff"),
                             value = 0.05,
                             min = 0,
                             max = 0.5)
                    ),
                    column(6,
                           numericInput(
                             ns("count.cutoff.reactome"),
                             tags$span(
                               class = "normal-label",
                               "Gene count cutoff"),
                             value = 5,
                             min = 0,
                             max = 1000
                           )
                    )
                  )
                )
              ),

              ### Embeddings Clustering ====
              tags$h4("Embeddings Clustering"),
              fluidRow(
                column(4,
                       numericInput(
                         ns("sim_cutoff"),
                         label = tags$span(
                           class = "normal-label",
                           "Similarity cutoff"
                         ),
                         value = 0.5,
                         min = 0,
                         max = 1,
                         step = 0.05
                       )),
                column(8,
                       selectInput(
                         ns("cluster_method"),
                         tags$span(
                           class = "normal-label",
                           "Clustering method"),
                         choices = c("Binary cut" = "binary cut",
                                     "Girvan Newman" = "girvan newman",
                                     "Hierarchical" = "hierarchical"),
                         selected = "Binary cut"))
              ),

              shinyjs::hidden(
                div(
                  id = ns("hclust.method_panel"),
                  fluidRow(
                    column(12,
                           selectInput(
                             ns("hclust.method"),
                             tags$span(
                               class = "normal-label",
                               "Linkage methods"),
                             choices = c("ward.D", "ward.D2", "single",
                                         "complete", "average (UPGMA)",
                                         "mcquitty (WPGMA)", "median (WPGMC)",
                                         "centroid (UPGMC)"),
                             selected = "complete"))
                  )
                )
              ),

              actionButton(
                ns("submit_merge_pathways"),
                "Submit",
                class = "btn-primary",
                style = "background-color: #d83428; color: white;"
              ),
              actionButton(
                ns("go2llm_interpretation"),
                "Next",
                class = "btn-primary",
                style = "background-color: #d83428; color: white;"
              ),
              actionButton(
                ns("show_merge_modules_code"),
                "Code",
                class = "btn-primary",
                style = "background-color: #d83428; color: white;"
              ),

              ### Add CSS class
              tags$head(
                tags$style(HTML("
                 .normal-label {
                   font-weight: normal !important;
                 }
                "))
              ),
              style = "border-right: 1px solid #ddd; padding-right: 20px;"
        ),

        ### Result display ====
        column(8,
               tabsetPanel(
                 tabPanel(
                   title = "Table",
                   shiny::dataTableOutput(ns("enriched_functional_modules")),
                   br(),
                   shinyjs::useShinyjs(),
                   downloadButton(ns("download_enriched_functional_modules"),
                                  "Download",
                                  class = "btn-primary",
                                  style = "background-color: #d83428; color: white;")
                 ),
                 tabPanel(
                   title = "Data visualization",
                   shiny::plotOutput(ns("enirched_functional_module_plot")),
                   br(),
                   fluidRow(
                     column(3,
                            actionButton(ns("generate_enirched_functional_module"),
                                         "Generate plot",
                                         class = "btn-primary",
                                         style = "background-color: #d83428; color: white;")
                     ),
                     column(3,
                            checkboxInput(ns("enirched_functional_module_plot_text"), "Text", FALSE)
                     ),
                     column(3,
                            checkboxInput(ns("enirched_functional_module_plot_text_all"), "Text all", FALSE)
                     ),
                     column(3,
                            numericInput(
                              ns("enirched_functional_module_plot_degree_cutoff"),
                              "Degree cutoff",
                              value = 0,
                              min = 0,
                              max = 1000)
                     )
                   )
                 ),
                 tabPanel(
                   title = "R object",
                   verbatimTextOutput(ns("enriched_functional_module_object")),
                   br(),
                   shinyjs::useShinyjs(),
                   downloadButton(ns("download_enriched_functional_module_object"),
                                  "Download",
                                  class = "btn-primary",
                                  style = "background-color: #d83428; color: white;")
                 )
               )
        )
      )
    )
  )
}

#' Embed and Cluster Pathways Server Module
#'
#' Internal server logic for the embedding-based clustering step that follows
#' pathway enrichment analysis.
#'
#' @param input,output,session Internal parameters for \strong{shiny}.
#'   **Do not remove.**
#' @param id Module id.
#' @param enriched_pathways A reactive value that returns the
#'   \code{enrich_pathway} object produced by the enrichment step.
#' @param tab_switch    A callback used to programmatically switch tabs in the
#'   parent \code{shinydashboard}.
#'
#' @import shiny
#' @importFrom shinyjs toggleElement toggleState disable enable useShinyjs
#' @importFrom utils capture.output write.csv
#' @noRd

embed_cluster_pathways_server <- function(id, enriched_pathways, enriched_functional_module, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ### Step0: Load enrichment analysis result ====
      observeEvent(input$upload_enriched_pathways, {
        if (!is.null(input$upload_enriched_pathways$datapath)) {
          message("Loading data")
          tempEnv <- new.env()
          load(input$upload_enriched_pathways$datapath,
               envir = tempEnv)

          names <- ls(tempEnv)

          if (length(names) == 1) {
            enriched_pathways$enriched_pathways_res <- get(names[1], envir = tempEnv)
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

      observe({
        req(input$cluster_method)

        shinyjs::toggleElement(
          id = "hclust.method_panel",
          condition = input$cluster_method == "hierarchical"
        )
      })

      query_type <- reactive({ enriched_pathways$query_type })

      observe({
        req(query_type())
        ## For gene
        shinyjs::toggleElement(
          id = "parameter_panel_gene",
          condition = query_type() == "gene"
        )

        ## For metabolite
        shinyjs::toggleElement(
          id = "parameter_panel_metabolite",
          condition = query_type() == "metabolite"
        )

        db_choices <- c("GO" = "go",
                        "KEGG" = "kegg",
                        "Reactome" = "reactome",
                        "HMDB" = "hmdb",
                        "KEGG" = "metkegg")
        updateCheckboxGroupInput(
          session, "cluster_module_database",
          choices  = db_choices[db_choices %in% enriched_pathways$available_db],
          selected = enriched_pathways$available_db
        )
      })

      observe({
        req(input$cluster_module_database)

        # Define database-panel pairs
        db_panels <- list(
          "go" = "gene_go_panel",
          "kegg" = "gene_kegg_panel",
          "reactome" = "gene_reactome_panel",
          "hmdb" = "metabolite_hmdb_panel",
          "metkegg" = "metabolite_metkegg_panel"
        )

        # Loop through each database-panel pair and toggle accordingly
        for (db in names(db_panels)) {
          shinyjs::toggleElement(
            id = db_panels[[db]],
            condition = db %in% input$cluster_module_database
          )
        }
      })

      merge_modules_code <- reactiveVal()

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
                ### Step1: Embedding ====
                bioembed_sim_matrix <-
                  get_bioembedsim(
                    object = enriched_pathways$enriched_pathways_res,
                    api_provider = input$api_provider,
                    text_embedding_model = input$embedding_model,
                    api_key = input$api_key,
                    database = input$cluster_module_database,
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
                    save_to_local = FALSE
                  )

                ### Step2: Clustering ====
                result <-
                  merge_pathways_bioembedsim(
                    object = bioembed_sim_matrix,
                    sim.cutoff = input$sim_cutoff,
                    cluster_method = input$cluster_method,
                    hclust.method = input$hclust.method,
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

            enriched_functional_module(result)

            ### Save code ====
            if (query_type() == "gene") {

              go_params <- ""
              if ("go" %in% input$cluster_module_database) {
                go_params <- sprintf(
                  '
                  p.adjust.cutoff.go = %s,
                  count.cutoff.go = %s,
                  ',
                  input$p.adjust.cutoff.go,
                  input$count.cutoff.go
                )}

              kegg_params <- ""
              if ("kegg" %in% input$cluster_module_database) {
                kegg_params <- sprintf(
                 'p.adjust.cutoff.kegg = %s,
                  count.cutoff.kegg = %s,
                  ',
                  input$p.adjust.cutoff.kegg,
                  input$count.cutoff.kegg
                )}

              reactome_params <- ""
              if ("reactome" %in% input$cluster_module_database) {
                reactome_params <- sprintf(
                 'p.adjust.cutoff.reactome = %s,
                  count.cutoff.reactome = %s',
                  input$p.adjust.cutoff.reactome,
                  input$count.cutoff.reactome
                )}

              db_vector <- paste0('c("', paste(input$cluster_module_database, collapse = '", "'), '")')

              merge_modules_code_str <- sprintf(
                '
              bioembed_sim_matrix <-
                get_bioembedsim(
                  object = enriched_pathways,
                  api_provider = %s,
                  text_embedding_model = %s,
                  api_key = %s,
                  include_gene_name = FALSE,
                  database = %s,%s%s%s
                  save_to_local = FALSE
                  )
               enriched_functional_module <-
                 merge_pathways_bioembedsim(
                   object = bioembed_sim_matrix,
                   sim.cutoff = %s,
                   cluster_method = %s,
                   hclust.method = %s,
                   save_to_local = FALSE
                 )
              ',
                input$api_provider,
                input$embedding_model,
                input$api_key,
                db_vector,
                go_params,
                kegg_params,
                reactome_params,
                input$sim_cutoff,
                input$cluster_method,
                input$hclust.method
              )

            } else if (query_type() == "metabolite") {
              hmdb_params <- ""
              if ("hmdb" %in% input$cluster_module_database) {
                hmdb_params <- sprintf(
                  '
                  p.adjust.cutoff.hmdb = %s,
                  count.cutoff.hmdb = %s,
                  ',
                  input$p.adjust.cutoff.hmdb,
                  input$count.cutoff.hmdb
                )}

              metkegg_params <- ""
              if ("metkegg" %in% input$cluster_module_database) {
                metkegg_params <- sprintf(
                 'p.adjust.cutoff.metkegg = %s,
                  count.cutoff.metkegg = %s,',
                  input$p.adjust.cutoff.metkegg,
                  input$count.cutoff.metkegg
                )}

              db_vector <- paste0('c("', paste(input$cluster_module_database, collapse = '", "'), '")')

              merge_modules_code_str <- sprintf(
                '
              bioembed_sim_matrix <-
                get_bioembedsim(
                  object = enriched_pathways,
                  api_provider = %s,
                  text_embedding_model = %s,
                  api_key = %s,
                  include_gene_name = FALSE,
                  database = %s,%s%s
                  save_to_local = FALSE
                  )
               enriched_functional_module <-
                 merge_pathways_bioembedsim(
                   object = bioembed_sim_matrix,
                   sim.cutoff = %s,
                   cluster_method = %s,
                   hclust.method = %s,
                   save_to_local = FALSE
                 )
              ',
                input$api_provider,
                input$embedding_model,
                input$api_key,
                db_vector,
                hmdb_params,
                metkegg_params,
                input$sim_cutoff,
                input$cluster_method,
                input$hclust.method
              )

            }

            merge_modules_code(merge_modules_code_str)
          }
        })

      ####show code
      observeEvent(input$show_merge_modules_code, {
        if (is.null(merge_modules_code()) ||
            length(merge_modules_code()) == 0) {
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
            merge_modules_code()
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

      output$enriched_functional_module_object <-
        renderText({
          req(enriched_functional_module())
          enriched_functional_module <- enriched_functional_module()
          captured_output1 <-
            capture.output(enriched_functional_module,
                           type = "message")
          captured_output2 <-
            capture.output(enriched_functional_module,
                           type = "output")
          captured_output <-
            c(captured_output1,
              captured_output2)
          paste(captured_output, collapse = "\n")
        })

      output$enriched_functional_modules <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_functional_module()@merged_module$functional_module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$download_enriched_functional_modules <-
        shiny::downloadHandler(
          filename = function() {
            "merged_modules.csv"
          },
          content = function(file) {
            write.csv(
              enriched_functional_module()@merged_module$functional_module_result,
              file,
              row.names = FALSE
            )
          }
        )

      output$download_enriched_functional_module_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_functional_module.rda"
          },
          content = function(file) {
            enriched_functional_module_res <- enriched_functional_module()
            save(enriched_functional_module_res, file = file)
          }
        )

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_functional_module()) ||
                length(enriched_functional_module()) == 0) {
              shinyjs::disable("download_enriched_functional_modules")
              shinyjs::disable("download_enriched_functional_module_object")
            } else {
              if (length(enriched_functional_module()@merged_module) == 0) {
                shinyjs::disable("download_enriched_functional_modules")
                shinyjs::disable("download_enriched_functional_module_object")
              } else {
                shinyjs::enable("download_enriched_functional_modules")
                shinyjs::enable("download_enriched_functional_module_object")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_enriched_functional_modules")
            shinyjs::disable("download_enriched_functional_module_object")
          }
        )
      })

      ### Step3:  Data visualization ====
      ###define enirched_functional_module_plot
      enirched_functional_module_plot <- reactiveVal()
      observeEvent(input$generate_enirched_functional_module, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional modules data available. Please 'Merge modules' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_functional_module(),
                  level = "functional_module",
                  degree_cutoff = input$enirched_functional_module_plot_degree_cutoff,
                  text = input$enirched_functional_module_plot_text,
                  text_all = input$enirched_functional_module_plot_text_all
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

          enirched_functional_module_plot(plot)
        }
      })

      output$enirched_functional_module_plot <-
        shiny::renderPlot({
          req(tryCatch(
            enirched_functional_module_plot(),
            error = function(e)
              NULL
          ))
        })

      ### Step4: Go to llm interpretation tab ====
      ####if there is not enriched_functional_module, show a warning message
      observeEvent(input$go2llm_interpretation, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "Please merge modules first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("llm_interpretation")
        }
      })
    }
  )
}
