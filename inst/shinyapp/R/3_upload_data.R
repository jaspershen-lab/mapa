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

               radioButtons(
                 ns("query_type"),
                 "Query type",
                 choices = c(
                   "Gene" = "gene",
                   "Metabolite" = "metabolite"
                 ),
                 inline = TRUE
               ),

               shinyjs::hidden(
                 div(id = ns("gene_panel"),
                     checkboxGroupInput(
                       ns("example_choice"),
                       "Example dataset",
                       choices = c(
                         "Pathway Enrichment Example" = "example_enrich_pathway",
                         "GSEA Example" = "example_gsea"
                       ),
                       selected = character(0)
                     ),

                     selectInput(
                       ns("organism"),
                       "Organism",
                       choices = c(
                         "Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                         "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
                         "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
                         "Fly (org.Dm.eg.db)" = "org.Dm.eg.db",
                         "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db",
                         "Arabidopsis (org.At.tair.db)" = "org.At.tair.db",
                         "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db",
                         "Worm (org.Ce.eg.db)" = "org.Ce.eg.db",
                         "Pig (org.Ss.eg.db)" = "org.Ss.eg.db",
                         "Bovine (org.Bt.eg.db)" = "org.Bt.eg.db",
                         "Rhesus (org.Mmu.eg.db)" = "org.Mmu.eg.db",
                         "Canine (org.Cf.eg.db)" = "org.Cf.eg.db",
                         "E. coli strain K12(org.EcK12.eg.db)" = "org.EcK12.eg.db",
                         "E coli strain Sakai" = "org.EcSakai.eg.db",
                         "Chicken (org.Gg.eg.db)" = "org.Gg.eg.db",
                         "Xenopus (org.Xl.eg.db)" = "org.Xl.eg.db",
                         "Chimp (org.Pt.eg.db)" = "org.Pt.eg.db",
                         "Anopheles (org.Ag.eg.db)" = "org.Ag.eg.db",
                         "Malaria (org.Pf.plasmo.db)" = "org.Pf.plasmo.db",
                         "Myxococcus xanthus DK 1622" = "org.Mxanthus.db"
                       ),
                       selected = "org.Hs.eg.db"
                     ),
                     helpText("Enter the name of an OrgDb package that is installed on your system.",
                              "Common examples: org.Hs.eg.db (Human), org.Mm.eg.db (Mouse), org.Rn.eg.db (Rat)",
                              "For the current list of OrgDb packages, visit: ",
                              tags$a(
                                href = "https://bioconductor.org/packages/release/BiocViews.html#___OrgDb",
                                "Bioconductor OrgDb packages",
                                target = "_blank"
                              )),

                     selectInput(
                       ns("id_type"),
                       "Input ID type",
                       choices = list(
                         "ENSEMBL" = "ensembl",
                         "UniProt" = "uniprot",
                         "EntrezID" = "entrezid"
                       ),
                       selected = "ensembl"
                     )
                 )),

               shinyjs::hidden(
                 div(id = ns("metabolite_panel"),
                     checkboxGroupInput(
                       ns("met_example_choice"),
                       "Example dataset",
                       choices = c(
                         "Pathway Enrichment Example" = "example_enrich_pathway"
                       ),
                       selected = character(0)
                     ),

                     selectInput(
                       ns("met_organism"),
                       "Organism",
                       choices = c(
                         "Human (hsa)"       = "hsa",
                         "Mouse (mmu)"       = "mmu",
                         "Rat (rno)"         = "rno",
                         "Fly (dme)"         = "dme",
                         "Zebrafish (dre)"   = "dre",
                         "Yeast (sce)"       = "sce",
                         "Worm (cel)"        = "cel",
                         "Pig (ssc)"         = "ssc",
                         "Bovine (bta)"      = "bta",
                         "Canine (cfa)"      = "cfa"
                       ),
                       selected = "hsa"
                     ),
                     helpText(
                       "The KEGG organism code is a three or four letter abbreviation.",
                       "Examples: 'hsa' (Human), 'mmu' (Mouse), 'rno' (Rat).",
                       "For a complete list of organism codes, visit: ",
                       tags$a(
                         href = "https://www.genome.jp/kegg/catalog/org_list.html",
                         "KEGG Organism Codes",
                         target = "_blank"
                       )
                     ),
                     selectInput(
                       ns("met_id_type"),
                       "Input ID type",
                       choices = list(
                         "KEGG ID" = "keggid",
                         "HMDB ID" = "hmdbid"
                       ),
                       selected = "hmdbid"
                     )
                 )),

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
#' Internal server function for handling data upload, conversion, and processing in the Shiny application.
#' This module manages file uploads, example data selection, ID conversion for both genes and metabolites,
#' and provides download capabilities for processed data.
#'
#' @param id Character string. The module ID for namespace isolation in Shiny.
#' @param processed_info Reactive value object. Stores the processed data information to be shared across modules.
#' @param tab_switch Function. Callback function to switch between different tabs in the application.
#'
#' @details
#' The module handles the following functionalities:
#' * File upload support for CSV and Excel files
#' * Example dataset selection
#' * ID conversion for genes (ENSEMBL, UniProt, EntrezID) and metabolites
#' * Dynamic UI updates based on query type
#' * Data validation and error handling
#' * Download capabilities for processed data
#' * Conversion code display
#'
#' @return None (creates a Shiny module server)
#'
#' @import shiny
#' @importFrom shinyjs toggleElement toggleState disable enable hidden useShinyjs
#' @importFrom readxl read_excel
#'
#' @noRd

upload_data_server <- function(id, processed_info, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      # Initialize reactive values for storing data throughout the module
      data_values <- reactiveValues(
        raw_data = NULL,        # Original uploaded or example data
        converted_data = NULL,  # Data after ID conversion
        conversion_code = NULL  # The code used for conversion
      )

      ## Toggle UI panels based on query type
      observe({
        shinyjs::toggleElement(
          id = "gene_panel",
          condition = input$query_type == "gene"
        )
      })

      observe({
        shinyjs::toggleElement(
          id = "metabolite_panel",
          condition = input$query_type == "metabolite"
        )
      })

      # Update metabolite ID type when organism changes
      observeEvent(input$met_organism, {
        if (input$met_organism != "hsa") {
          updateSelectInput(
            session,
            "met_id_type",
            choices = list(
              "KEGG ID" = "keggid"
            ),
            selected = "keggid"
          )
        }
      })

      # Load data from file upload or example selection
      # Observer for gene example data selection
      observeEvent(input$example_choice, {
        if (input$query_type == "gene" && length(input$example_choice) > 0) {
          # Example data selected
          example_path <- switch(input$example_choice,
                                 "example_enrich_pathway" = "inst/shinyapp/files/example_enrich_pathway.csv",
                                 "example_gsea" = "inst/shinyapp/files/example_gsea.csv",
                                 NULL)

          if (!is.null(example_path)) {
            tryCatch({
              data_values$raw_data <- read.csv(example_path)
            }, error = function(e) {
              showNotification(paste("Failed to load example data:", e$message), type = "error")
            })
          }
        }
      }, ignoreNULL = FALSE)

      # Observer for metabolite example data selection
      observeEvent(input$met_example_choice, {
        if (input$query_type == "metabolite" && length(input$met_example_choice) > 0) {
          # Example data selected
          example_path <- switch(input$met_example_choice,
                                 "example_enrich_pathway" = "inst/shinyapp/files/example_enrich_pathway_metabolite.csv",
                                 NULL)

          if (!is.null(example_path)) {
            tryCatch({
              data_values$raw_data <- read.csv(example_path)
            }, error = function(e) {
              showNotification(paste("Failed to load example data:", e$message), type = "error")
            })
          }
        }
      }, ignoreNULL = FALSE)

      # Observer for file uploads (works for both gene and metabolite)
      observeEvent(input$variable_info, {
        if (!is.null(input$variable_info)) {
          # File uploaded
          in_file <- input$variable_info
          tryCatch({
            if (grepl("\\.csv$", in_file$name)) {
              data_values$raw_data <- read.csv(in_file$datapath)
            } else if (grepl("\\.(xlsx|xls)$", in_file$name)) {
              data_values$raw_data <- readxl::read_excel(in_file$datapath)
            } else {
              showNotification("Unsupported file format. Please upload CSV or Excel file.", type = "error")
            }
          }, error = function(e) {
            showNotification(paste("Error reading uploaded file:", e$message), type = "error")
          })
        }
      }, ignoreNULL = TRUE)

      # Process data when Submit button is clicked
      observeEvent(input$map_id, {
        if (is.null(data_values$raw_data)) {
          showModal(modalDialog(
            title = "Warning",
            "No data is available. Please upload a file or select an example dataset.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          return()
        }

        # Process based on query type
        if (input$query_type == "gene") {
          # Set up conversion parameters
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

          # Validate OrgDb format
          if (!grepl("^org\\.[A-Za-z]+\\..+\\.db$", input$organism)) {
            showNotification("Invalid OrgDb package name. Expected format: org.XX.eg.db", type = "error")
            return()
          }

          # Check if package is installed
          if (!requireNamespace(input$organism, quietly = TRUE)) {
            showModal(modalDialog(
              title = "Missing Package",
              paste0("Package ", input$organism, " is not installed. Please install it using:\n",
                     "BiocManager::install('", input$organism, "')"),
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
            return()
          }

          # Load the package and get OrgDb object
          requireNamespace(input$organism)
          org_db_obj <- get(input$organism)

          # Perform conversion
          tryCatch({
            result <- id_conversion(
              query_type = input$query_type,
              data = data_values$raw_data,
              from_id_type = conversion_params$from_id_type,
              to_id_type = conversion_params$to_id_type,
              organism = org_db_obj
            )

            data_values$converted_data <- result$converted_id
            data_values$conversion_code <- result$conversion_code

            processed_info$variable_info <- result$converted_id
            processed_info$query_type <- input$query_type
            processed_info$organism <- input$organism
            # Show success message
            showNotification("Data successfully processed", type = "message")
          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              paste("Conversion failed:", e$message),
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          })
        }
        else if (input$query_type == "metabolite") {
          # Perform metabolite conversion
          tryCatch({
            result <- id_conversion(
              query_type = input$query_type,
              data = data_values$raw_data,
              from_id_type = input$met_id_type,
              organism = input$met_organism
            )

            data_values$converted_data <- result$converted_id
            data_values$conversion_code <- result$conversion_code

            processed_info$variable_info <- result$converted_id
            processed_info$query_type <- input$query_type
            processed_info$organism <- input$met_organism

            # Show success message
            showNotification("Data successfully processed", type = "message")
          }, error = function(e) {
            showModal(modalDialog(
              title = "Error",
              paste("Conversion failed:", e$message),
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          })
        }
      })

      # Render data table with converted data
      output$show_variable_info <- shiny::renderDataTable({
        if (!is.null(data_values$converted_data)) {
          data_values$converted_data
        } else if (!is.null(data_values$raw_data)) {
          data_values$raw_data
        }
      }, options = list(pageLength = 10, scrollX = TRUE))

      # Handle download button state
      observe({
        if (is.null(data_values$converted_data)) {
          shinyjs::disable("download_variable_info")
        } else {
          shinyjs::enable("download_variable_info")
        }
      })

      # Download handler for processed data
      output$download_variable_info <- shiny::downloadHandler(
        filename = function() {
          "processed_data.csv"
        },
        content = function(file) {
          write.csv(data_values$converted_data, file, row.names = FALSE)
        }
      )

      # Show conversion code when requested
      observeEvent(input$show_conversion_code, {
        if (is.null(data_values$conversion_code)) {
          showModal(modalDialog(
            title = "Warning",
            "No conversion code available. Please process data first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        } else {
          showModal(modalDialog(
            title = "Conversion Code",
            tags$pre(data_values$conversion_code),
            easyClose = TRUE,
            size = "l",
            footer = modalButton("Close")
          ))
        }
      })

      # Handle navigation to next tab
      observeEvent(input$go2enrich_pathways, {
        if (is.null(data_values$converted_data)) {
          showModal(modalDialog(
            title = "Warning",
            "Please process data before proceeding to the next step.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        } else {
          tab_switch("enrich_pathways")
        }
      })
    }
  )
}
