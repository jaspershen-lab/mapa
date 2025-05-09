#' Merge Modules UI Module
#'
#' Internal UI for merging enriched functional modules.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @importFrom shinyBS bsButton bsPopover
#' @noRd

merge_modules_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "merge_modules",
    fluidPage(titlePanel("Merge Modules"),
              fluidPage(
                fluidRow(
                  column(4,
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("measure.method.module"),
                                    "Similarity method",
                                    choices = c("jaccard"),
                                    selected = "jaccard")
                           ),
                           column(6,
                                  numericInput(
                                    ns("sim.cutoff.module"),
                                    "Similarity cutoff",
                                    value = 0.5,
                                    min = 0,
                                    max = 1)
                           )),

                         actionButton(
                           ns("submit_merge_modules"),
                           "Submit",
                           class = "btn-primary",
                           style = "background-color: #d83428; color: white;"),

                         actionButton(
                           ns("go2llm_interpretation"),
                           "Next",
                           class = "btn-primary",
                           style = "background-color: #d83428; color: white;"),

                         actionButton(
                           ns("show_merge_modules_code"),
                           "Code",
                           class = "btn-primary",
                           style = "background-color: #d83428; color: white;"),
                         style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
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
                                        ns("enirched_functional_moduleplot__degree_cutoff"),
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
                         ))
                )))
  )
}


#' Merge Modules Server Module
#'
#' Internal server logic for merging enriched modules into functional modules.
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @param id Module id.
#' @param enriched_modules Reactive value containing enriched modules.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleState useShinyjs
#' @importFrom clusterProfiler merge_modules
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @noRd

merge_modules_server <- function(id, enriched_modules = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      # Define enriched_functional_module as a reactive value
      enriched_functional_module <- reactiveVal()
      merge_modules_code <- reactiveVal()

      ### merge modules ====
      observeEvent(input$submit_merge_modules, {
        # Check if enriched_modules is available
        if (is.null(enriched_modules()) ||
            length(enriched_modules()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched modules data available. Please 'Enrich modules' first.",
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
                merge_modules(
                  object = enriched_modules(),
                  sim.cutoff = input$sim.cutoff.module,
                  measure_method = input$measure.method.module,
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

          enriched_functional_module(result)
          # shinyjs::hide("loading")

          ##save code
          merge_modules_code <-
            sprintf(
              '
            enriched_functional_module <-
            merge_modules(
            object = enriched_modules,
            sim.cutoff = %s,
            measure_method = %s)
            ',
              input$sim.cutoff.module,
              paste0('"', input$measure.method.module, '"')
            )

          merge_modules_code(merge_modules_code)
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

      observe({
        tryCatch(
          expr = {
            if (is.null(enriched_functional_module()) ||
                length(enriched_functional_module()) == 0) {
              shinyjs::disable("download_enriched_functional_modules")
            } else {
              if (length(enriched_functional_module()@merged_module) == 0) {
                shinyjs::disable("download_enriched_functional_modules")
              } else {
                shinyjs::enable("download_enriched_functional_modules")
              }
            }
          },
          error = function(e) {
            shinyjs::disable("download_enriched_functional_modules")
          }
        )
      })

      ### data visualization ====
      ###define enirched_functional_module_plot
      enirched_functional_module_plot <-
        reactiveVal()
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
          # shinyjs::show("loading")
          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_similarity_network(
                  object = enriched_functional_module(),
                  level = "functional_module",
                  degree_cutoff = input$enirched_functional_moduleplot__degree_cutoff,
                  text = input$enirched_functional_module_plot_text,
                  text_all = input$enirched_functional_module_plot_text_all
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

          enirched_functional_module_plot(plot)
          # shinyjs::hide("loading")
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


      output$download_enriched_functional_module_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_functional_module.rda"
          },
          content = function(file) {
            save(enriched_functional_module(), file = file)
          }
        )

      observe({
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          shinyjs::disable("download_enriched_functional_module_object")
        } else {
          shinyjs::enable("download_enriched_functional_module_object")
        }
      })


      output$download_enriched_functional_module_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_functional_module.rda"
          },
          content = function(file) {
            enriched_functional_module <-
              enriched_functional_module()
            save(enriched_functional_module, file = file)
          }
        )

      observe({
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          shinyjs::disable("download_enriched_functional_module_object")
        } else {
          shinyjs::enable("download_enriched_functional_module_object")
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


      ### Go to llm interpretation tab
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

      return(enriched_functional_module)
    }
  )
}
