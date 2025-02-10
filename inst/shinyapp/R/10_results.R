#' Results UI Module
#'
#' Internal UI for generating and displaying the results report.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @noRd

results_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "results",
    fluidPage(titlePanel("Results and Report"),
              fluidPage(
                column(4,
                       br(),
                       fluidRow(
                         actionButton(
                           inputId = ns("generate_report"),
                           label = "Generate report",
                           class = "btn-primary",
                           style = "background-color: #d83428; color: white;"
                         ),
                         shinyjs::useShinyjs(),
                         downloadButton(ns("download_report"),
                                        "Download",
                                        class = "btn-primary",
                                        style = "background-color: #d83428; color: white;"),
                         actionButton(
                           inputId = ns("show_report_code"),
                           label = "Code",
                           class = "btn-primary",
                           style = "background-color: #d83428; color: white;"
                         )
                       ),
                       style = "border-right: 1px solid #ddd; padding-right: 20px;"
                ),
                column(8,
                       tabsetPanel(
                         tabPanel(
                           title = "Report",
                           uiOutput(ns("mapa_report"))
                         ))
                )
              )
    )
  )
}

#' Results Server Module
#'
#' Internal server logic for generating the report, displaying its content,
#' and handling report download.
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @param id Module id.
#' @param enriched_functional_module Reactive value containing enriched functional module data.
#' @param llm_interpretation_result Reactive value containing LLM interpretation result.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs useShinyjs
#' @noRd

results_server <- function(id, enriched_functional_module = NULL, llm_interpretation_result = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {

      report_code <- reactiveVal()
      report_path <- reactiveVal()

      observeEvent(input$generate_report, {
        # Check if enriched_functional_module and llm_interpretation_result are
        #  available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional modules data available.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          if (is.null(llm_interpretation_result()) ||
              length(llm_interpretation_result()) == 0) {
            showModal(
              modalDialog(
                title = "Warning",
                "No LLM interpreatation result available.",
                easyClose = TRUE,
                footer = modalButton("Close")
              )
            )
          } else{
            # shinyjs::show("loading")

            withProgress(message = 'Analysis in progress...', {
              tryCatch({
                report_path <-
                  file.path("files",
                            paste(sample(
                              c(letters, LETTERS, 0:9),
                              30, replace = TRUE
                            ), collapse = ""))

                report_functional_module(
                  object = enriched_functional_module(),
                  interpretation_result = llm_interpretation_result(),
                  path = report_path,
                  type = "html"
                )
              },
              error = function(e) {
                showModal(
                  modalDialog(
                    title = "Error",
                    paste("Details:", e$message),
                    easyClose = TRUE,
                    footer = modalButton("Close")
                  )
                )
              })
            })

            report_path(report_path)

            # shinyjs::hide("loading")

            ##save code
            report_code <-
              sprintf(
                '
              report_functional_module(
              object = enriched_functional_module,
              interpretation_result = interpretation_result,
              path = %s,
              type = "html")
            ',
                paste0('"', report_path(), '"')
              )
            report_code(report_code)
          }
        }
      })

      # Update UI to display HTML
      output$mapa_report <-
        renderUI({
          tryCatch(
            includeHTML(file.path(
              report_path(), "Report/mapa_report.html"
            )),
            error = function(e) {
              NULL
            }
          )
        })

      ###donload the zip report file
      output$download_report <-
        shiny::downloadHandler(
          filename = function() {
            "Report.zip"
          },
          content = function(file) {
            zip_path <-
              paste0(report_path(), "/Report.zip")
            zip(zipfile = zip_path,
                files = paste0(report_path(),
                               "/Report"))
            file.copy(zip_path, file)
          }
        )

      observe({
        if (is.null(report_path()) ||
            length(report_path()) == 0) {
          shinyjs::disable("download_report")
        } else {
          shinyjs::enable("download_report")
        }
      })

      ###To delete the zip file and folder when the user closes the app
      session$onSessionEnded(function() {
        all_files_folders <-
          list.files("files", full.names = TRUE)

        folders <-
          Filter(function(x) {
            file.info(x)$isdir
          }, all_files_folders)

        regex_pattern <- "^[A-Za-z0-9]{30}$"

        report_dirs <-
          Filter(function(folder) {
            folder_name <- basename(folder)
            grepl(regex_pattern, folder_name)
          }, folders)

        unlink(report_dirs, recursive = TRUE)
      })

      ####show code
      observeEvent(input$show_report_code, {
        if (is.null(report_code()) ||
            length(report_code()) == 0) {
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
            report_code()
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
    }
  )
}
