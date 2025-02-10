#' Upload Data Module UI
#'
#' Internal UI for data upload.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs hidden
#' @noRd

upload_data_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "upload_data",
    fluidPage(
      titlePanel("Upload Data"),
      fluidRow(
        column(4,
               fileInput(
                 ns("variable_info"),
                 "Choose marker information",
                 accept = c(
                   "text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv",
                   ".xlsx",
                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                   "application/vnd.ms-excel"
                 )
               ),
               checkboxInput(
                 ns("use_example"),
                 "Use example",
                 FALSE
               ),
               shinyjs::hidden(
                 div(id = ns("example_panel"),
                     radioButtons(
                       ns("example_choice"),
                       "Example dataset",
                       choices = c(
                         "Pathway Enrichment Example" = "example_enrich_pathway",
                         "GSEA Example" = "example_gsea"
                       )
                     )
                 )
               ),

               radioButtons(
                 ns("id_type"),
                 "ID type",
                 choices = list(
                   "ENSEMBL" = "ensembl",
                   "UniProt" = "uniprot",
                   "EntrezID" = "entrezid"
                 ),
                 selected = "ensembl"
               ),

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
    ))
}

#' Upload Data Module Server
#'
#' Internal server for data upload and conversion.
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @param id Module id.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleElement toggleState disable enable hidden useShinyjs
#' @importFrom readxl read_excel
#' @noRd

upload_data_server <- function(id, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ## Upload variable info ====

      ### Toggle example selection visibility
      observe({
        shinyjs::toggleElement(
          id = "example_panel",
          condition = input$use_example
        )
      })

      ### Disable file upload when using example
      observe({
        shinyjs::toggleState(
          id = "variable_info",
          condition = !input$use_example
        )
      })

      variable_info <- reactive({
        if (input$use_example) {
          req(input$example_choice)
          tryCatch({
            switch(input$example_choice,
                   "example_enrich_pathway" = read.csv("inst/shinyapp/files/example_enrich_pathway.csv"),
                   "example_gsea" = read.csv("inst/shinyapp/files/example_gsea.csv")
            )
          }, error = function(e) {
            showNotification("Failed to load example data", type = "error")
            NULL
          })
        } else {
          req(input$variable_info)
          in_file <- input$variable_info
          tryCatch({
            if (grepl("\\.csv$", in_file$name)) {
              read.csv(in_file$datapath)
            } else {
              readxl::read_excel(in_file$datapath)
            }
          }, error = function(e) {
            showNotification("Error reading uploaded file", type = "error")
            NULL
          })
        }
      })

      ## Define new reactive values
      variable_info_new <- reactiveVal()
      conversion_code <- reactiveVal()

      ## If the user don't upload data and don't use example data,
      ### Show a warning message
      observeEvent(input$map_id, {
        if (is.null(variable_info())) {
          showModal(
            modalDialog(
              title = "Warning",
              "No data is uploaded",
              easyClose = FALSE,
              footer = modalButton("Close")
            )
          )
        }
      })

      ## ID conversion ====
      observeEvent(input$map_id, {
        req(variable_info())
        variable_info_old <- variable_info()

        conversion_params <- switch(input$id_type,
                                    "ensembl" = list(
                                      from_id_type = "ENSEMBL",
                                      to_id_type = c("UNIPROT", "ENTREZID", "SYMBOL")
                                    ),
                                    "uniprot" = list(
                                      from_id_type = "UNIPROT",
                                      to_id_type = c("ENSEMBL", "ENTREZID", "SYMBOL")
                                    ),
                                    "entrezid" = list(
                                      from_id_type = "ENTREZID",
                                      to_id_type = c("ENSEMBL", "UNIPROT", "SYMBOL")
                                    )
        )

        tryCatch({
          # Process data conversion
          converted_data <- id_conversion(
            data = variable_info_old,
            from_id_type = conversion_params$from_id_type,
            to_id_type = conversion_params$to_id_type,
            orgDb = org.Hs.eg.db
          )
        }, error = function(e) {
          showModal(modalDialog(
            title = "Error",
            paste("Conversion failed:", e$message),
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        })

        # Update reactive values
        variable_info_new(converted_data$converted_id)
        conversion_code(converted_data$conversion_code)
      })

      output$show_variable_info <-
        shiny::renderDataTable({
          req(variable_info_new())
          variable_info_new()
        },
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )

      ## Download the variable_info ====
      output$download_variable_info <-
        shiny::downloadHandler(
          filename = function() {
            "variable_info.csv"
          },
          content = function(file) {
            write.csv(variable_info_new(),
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(variable_info_new()) ||
            length(variable_info_new()) == 0) {
          shinyjs::disable("download_variable_info")
        } else {
          shinyjs::enable("download_variable_info")
        }
      })

      ## Show code ====
      observeEvent(input$show_conversion_code, {
        if (is.null(conversion_code()) ||
            length(conversion_code()) == 0) {
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
            conversion_code()
          showModal(modalDialog(
            title = "Code",
            tags$pre(code_content),
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        }
      })

      ####Click next and then go to enrich pathways tab
      observeEvent(input$go2enrich_pathways, {
        # Check if variable_info_new is available
        if (is.null(variable_info_new()) ||
            length(variable_info_new()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "Please upload data and click 'Submit' first.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # Navigate to the enrich_pathways tab
          tab_switch("enrich_pathways")
        }
      })

      return(variable_info_new)
    }
  )
}
