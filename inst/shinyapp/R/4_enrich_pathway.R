#' Enrich Pathways UI Module
#'
#' Internal server for enrichment analysis. For genes, both pathway enrichment analysis and Gene set enrichment analysis can be performed. For metabolites, only pathway enrichment analysis is available.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs hidden toggleElement useShinyjs
#' @noRd

enrich_pathway_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "enrich_pathways",
          fluidPage(
            titlePanel("Enrich Pathways"),
            fluidPage(
              fluidRow(
                column(4,
                       radioButtons(
                         ns("query_type"),
                         tags$h4("Query type"),
                         choices = c(
                           "Gene" = "gene",
                           "Metabolite" = "metabolite"
                         ),
                         inline = TRUE
                       ),

                       ## Gene enrichment analysis panel ----
                       shinyjs::hidden(
                         div(id = ns("gene_panel"),
                             radioButtons(
                               ns("analysis_type"),
                               tags$h4("Analysis type"),
                               choices = list(
                                 "Pathway enrichment analysis" = "enrich_pathway",
                                 "Gene set enrichment analysis" = "do_gsea"
                               )
                             ),
                             shinyjs::hidden(
                               div(id = ns("gsea_order_by_panel"),
                                   radioButtons(
                                     ns("order_by"),
                                     tags$h4("Order by"),
                                     choices = c(
                                       "Fold change" = "fc",
                                       "Adjusted p value" = "p_value_adjust"
                                     )
                                   )
                               )
                             ),
                             checkboxGroupInput(
                               ns("pathway_database"),
                               tags$h4("Database"),
                               choices = c(
                                 "GO" = "go",
                                 "KEGG" = "kegg",
                                 "Reactome" = "reactome"
                               ),
                               selected = NULL
                             ),

                             # GO specific parameters
                             shinyjs::hidden(
                               div(id = ns("go_panel"),
                                   h4("GO Parameters"),
                                   textInput(
                                     ns("go_orgdb"),
                                     "Organism Database",
                                     value = "org.Hs.eg.db"
                                   ),
                                   helpText("Enter the name of an OrgDb package that is installed on your system.",
                                            "Common examples: org.Hs.eg.db (Human), org.Mm.eg.db (Mouse), org.Rn.eg.db (Rat)",
                                            "For the current list of OrgDb packages, visit: ",
                                            tags$a(
                                              href = "https://bioconductor.org/packages/release/BiocViews.html#___OrgDb",
                                              "Bioconductor OrgDb packages",
                                              target = "_blank"
                                            )),
                                   textInput(
                                     ns("go_keytype"),
                                     "GO Keytype",
                                     value = "ENSEMBL"
                                   ),
                                   selectInput(
                                     ns("go_ont"),
                                     "GO Ontology",
                                     choices = c(
                                       "All" = "ALL",
                                       "Biological Process" = "BP",
                                       "Cellular Component" = "CC",
                                       "Molecular Function" = "MF"
                                     ),
                                     selected = "ALL"
                                   )
                               )
                             ),

                             # KEGG specific parameters
                             shinyjs::hidden(
                               div(id = ns("kegg_panel"),
                                   h4("KEGG Parameters"),
                                   textInput(
                                     ns("kegg_organism"),
                                     "KEGG Organism",
                                     value = "NULL"
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
                                     ns("kegg_keytype"),
                                     "KEGG Keytype",
                                     choices = c(
                                       "KEGG/Entrez" = "kegg",
                                       "NCBI Gene ID" = "ncbi-geneid",
                                       "NCBI Protein ID" = "ncbi-proteinid",
                                       "UniProt" = "uniprot"
                                     ),
                                     selected = "uniprot"
                                   ),
                                   checkboxInput(
                                     ns("use_internal_data"),
                                     "Use internal database",
                                     value = TRUE
                                   )
                               )
                             ),

                             # Reactome specific parameters
                             shinyjs::hidden(
                               div(id = ns("reactome_panel"),
                                   h4("Reactome Parameters"),
                                   selectInput(
                                     ns("reactome_organism"),
                                     "Reactome Organism",
                                     choices = c(
                                       "Human" = "human",
                                       "Rat" = "rat",
                                       "Mouse" = "mouse",
                                       "C. elegans" = "celegans",
                                       "Yeast" = "yeast",
                                       "Zebrafish" = "zebrafish",
                                       "Fruit fly" = "fly"
                                     ),
                                     selected = "human"
                                   )
                               )
                             ),

                             # Common parameters
                             numericInput(
                               ns("p_value_cutoff"),
                               "P-value cutoff",
                               value = 0.05,
                               min = 0,
                               max = 0.5),
                             selectInput(
                               ns("p_adjust_method"),
                               "P-adjust method",
                               choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                               selected = "BH"),
                             numericInput(
                               ns("q_value_cutoff"),
                               "Q-value cutoff",
                               value = 0.2,
                               min = 0,
                               max = 1
                             ),
                             sliderInput(
                               ns("gene_set_size"),
                               "Gene set size",
                               min = 5,
                               max = 2000,
                               value = c(10, 500))
                         )
                       ),

                       ## Metabolite enrichment analysis panel ----
                       shinyjs::hidden(
                         div(id = ns("metabolite_panel"),
                             checkboxGroupInput(
                               ns("pathway_database"),
                               "Database",
                               choices = c(
                                 "HMDB" = "hmdb",
                                 "KEGG" = "kegg"
                               ),
                               selected = NULL
                             ),
                             checkboxInput(
                               ns("use_internal_data"),
                               "Use internal database",
                               value = TRUE
                             ),
                             numericInput(
                               ns("p_value_cutoff"),
                               "P-value cutoff",
                               value = 0.05,
                               min = 0,
                               max = 0.5),
                             selectInput(
                               ns("p_adjust_method"),
                               "P-adjust method",
                               choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                               selected = "BH")
                         )
                       ),

                       actionButton(
                         ns("submit_enrich_pathways"),
                         "Submit",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),
                       actionButton(
                         ns("go2merge_pathways"),
                         "Next",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),
                       actionButton(
                         ns("show_enrich_pathways_code"),
                         "Code",
                         class = "btn-primary",
                         style = "background-color: #d83428; color: white;"),

                       style = "border-right: 1px solid #ddd; padding-right: 20px;"
                ),

                ## Results panel ----
                column(8,
                       tabsetPanel(
                         id = ns("enriched_pathways_result"),

                         tabPanel(
                           title = "GO",
                           shiny::dataTableOutput(ns("enriched_pathways_go")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_go"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "KEGG",
                           shiny::dataTableOutput(ns("enriched_pathways_kegg")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_kegg"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "Reactome",
                           shiny::dataTableOutput(ns("enriched_pathways_reactome")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_reactome"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "HMDB",
                           shiny::dataTableOutput(ns("enriched_pathways_hmdb")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_hmdb"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "KEGG Metabolite",
                           shiny::dataTableOutput(ns("enriched_pathways_metkegg")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_metkegg"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "R object",
                           verbatimTextOutput(ns("enriched_pathways_object")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_enriched_pathways_object"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         )
                       )
                )
              )
            )
          ))
}


#' Enrich Pathways Server Module
#'
#' Internal server for enrichment analysis. For genes, both pathway enrichment analysis and Gene set enrichment analysis can be performed. For metabolites, only pathway enrichment analysis is available.
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @param id Module id.
#' @param variable_info Reactive variable containing uploaded variable info.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleElement disable enable useShinyjs
#' @importFrom clusterProfiler enrich_pathway do_gsea
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA enrichPathway
#' @noRd

enrich_pathway_server <- function(id, variable_info = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ## Perform gene enrichment analysis ----
      observe({
        shinyjs::toggleElement(
          id = "gene_panel",
          condition = input$query_type == "gene"
        )
      })

      observe({
        shinyjs::toggleElement(
          id = "gsea_order_by_panel",
          condition = input$analysis_type == "do_gsea"
        )
      })

      # Toggle database-specific parameter panels
      observe({
        shinyjs::toggleElement(
          id = "go_panel",
          condition = "go" %in% input$pathway_database
        )
      })

      observe({
        shinyjs::toggleElement(
          id = "kegg_panel",
          condition = "kegg" %in% input$pathway_database
        )
      })

      observe({
        shinyjs::toggleElement(
          id = "reactome_panel",
          condition = "reactome" %in% input$pathway_database
        )
      })

      ## Perform metabolite enrichment analysis ----
      observe({
        shinyjs::toggleElement(
          id = "metabolite_panel",
          condition = input$query_type == "metabolite"
        )
      })

      ### Define enriched_pathways as a reactive value
      enriched_pathways <- reactiveVal()
      enrich_pathways_code <- reactiveVal()

      ### Submit variable info for enrichment analysis
      observeEvent(input$submit_enrich_pathways, {
        ### Check if variable_info is available
        if (is.null(variable_info()) || length(variable_info()) == 0) {
          showModal(modalDialog(
            title = "Warning",
            "No data available. Please 'Upload data' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
        } else {
          if (length(input$pathway_database) == 0) {
            showModal(modalDialog(
              title = "Warning",
              "Please select at least one pathway database.",
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          } else {
            withProgress(message = 'Analysis in progress...', {
              result <- tryCatch({
                library(clusterProfiler)
                library(ReactomePA)

                if (input$analysis_type == "enrich_pathway") {
                  # Extract common parameters
                  common_params <- list(
                    variable_info = variable_info(),
                    query_type = input$query_type,
                    database = input$pathway_database,
                    save_to_local = FALSE,
                    pvalueCutoff = input$p_value_cutoff,
                    pAdjustMethod = input$p_adjust_method,
                    use_internal_data = input$use_internal_data
                  )

                  # Add gene-specific parameters if query type is gene
                  if (input$query_type == "gene") {
                    common_params$qvalueCutoff <- input$q_value_cutoff
                    # GO parameters
                    if ("go" %in% input$pathway_database) {
                      common_params$go.keytype <- input$go_keytype
                      common_params$go.ont <- input$go_ont

                      # Validate input format
                      if (!grepl("^org\\.[A-Za-z]+\\..+\\.db$", input$go_orgdb)) {
                        stop("Invalid OrgDb package name. Expected format: org.XX.eg.db")
                      }
                      # Check if the package is installed
                      if (!requireNamespace(input$go_orgdb, quietly = TRUE)) {
                        stop(paste("Package", input$go_orgdb, "is not installed. Please install it using BiocManager::install('", input$go_orgdb, "')"))
                      }
                      # Load the package
                      requireNamespace(input$go_orgdb)
                      # Get the OrgDb object
                      org_db_obj <- get(input$go_orgdb)
                      common_params$go.orgdb <- org_db_obj
                    }

                    # KEGG parameters
                    if ("kegg" %in% input$pathway_database) {
                      common_params$kegg.organism <- input$kegg_organism
                      common_params$kegg.keytype <- input$kegg_keytype
                    }

                    # Reactome parameters
                    if ("reactome" %in% input$pathway_database) {
                      common_params$reactome.organism <- input$reactome_organism
                    }

                    # Gene set size parameters
                    common_params$minGSSize <- input$gene_set_size[1]
                    common_params$maxGSSize <- input$gene_set_size[2]
                  }

                  do.call(enrich_pathway, common_params)

                } else if (input$analysis_type == "do_gsea") {
                  do_gsea(
                    variable_info(),
                    order_by = input$order_by,
                    database = input$pathway_database,
                    save_to_local = FALSE,
                    path = "result",
                    OrgDb = org.Hs.eg.db,
                    organism = input$organism,
                    ont = "ALL",
                    pvalueCutoff = input$p_value_cutoff,
                    pAdjustMethod = input$p_adjust_method,
                    qvalueCutoff = input$q_value_cutoff,
                    minGSSize = input$gene_set_size[1],
                    maxGSSize = input$gene_set_size[2],
                    readable = FALSE,
                    pool = FALSE
                  )
                }
              }, error = function(e) {
                showModal(modalDialog(
                  title = "Error",
                  paste("Details:", e$message),
                  easyClose = TRUE,
                  footer = modalButton("Close")
                ))
                return(NULL)
              })

              enriched_pathways(result)

              # shinyjs::hide("loading")

              # Save code for reproducibility
              pathway_database <-
                paste0("c(", paste(unlist(
                  lapply(paste(input$pathway_database), function(x)
                    paste0('"', x, '"'))
                ),
                collapse = ", "), ")")

              if (input$analysis_type == "enrich_pathway") {
                if (input$query_type == "gene") {
                  # Build parameter parts based on selected databases
                  go_params <- ""
                  if ("go" %in% input$pathway_database) {
                    go_params <- sprintf(
                    '
                    go.orgdb = "%s",
                    go.keytype = "%s",
                    go.ont = "%s",
                    go.universe = NULL,
                    go.pool = FALSE,
                    ',
                    input$go_orgdb,
                    input$go_keytype,
                    input$go_ont
                    )}

                  kegg_params <- ""
                  if ("kegg" %in% input$pathway_database) {
                    kegg_params <- sprintf(
                    '
                    kegg.organism = "%s",
                    kegg.keytype = "%s",
                    use_internal_data = "%s",
                    kegg.universe = NULL,
                    ',
                    input$kegg_organism,
                    input$kegg_keytype,
                    as.character(input$use_internal_data)
                    )}

                  reactome_params <- ""
                  if ("reactome" %in% input$pathway_database) {
                    reactome_params <- sprintf(
                    '
                    reactome.organism = "%s",
                    reactome.universe = NULL,
                    ',
                    input$reactome_organism
                    )}

                  code <- sprintf(
                  '
                  enriched_pathways <-
                    enrich_pathway(
                      variable_info,
                      query_type = "%s",
                      database = %s,%s%s%s
                      pvalueCutoff = %s,
                      pAdjustMethod = "%s",
                      qvalueCutoff = %s,
                      minGSSize = %s,
                      maxGSSize = %s,
                      readable = FALSE,
                      save_to_local = FALSE
                    )',
                  input$query_type,
                  pathway_database,
                  go_params,
                  kegg_params,
                  reactome_params,
                  input$p_value_cutoff,
                  input$p_adjust_method,
                  input$q_value_cutoff,
                  input$gene_set_size[1],
                  input$gene_set_size[2]
                  )

                  enrich_pathways_code(code)

                  } else { # Metabolite code
                    code <- sprintf(
                    '
                    enriched_pathways <-
                    enrich_pathway(
                      variable_info,
                      query_type = "%s",
                      database = %s,
                      use_internal_data = %s,
                      pvalueCutoff = %s,
                      pAdjustMethod = "%s"
                    )
                    ',
                    input$query_type,
                    pathway_database,
                    as.character(input$use_internal_data),
                    input$p_value_cutoff,
                    input$p_adjust_method
                    )
                    enrich_pathways_code(code)
                    }
                } else { # GSEA code
                  code <- sprintf(
                  '
                  enriched_pathways <-
                  do_gsea(
                    variable_info,
                    order_by = "%s",
                    database = %s,
                    OrgDb = org.Hs.eg.db,
                    organism = "%s",
                    ont = "%s",
                    pvalueCutoff = %s,
                    pAdjustMethod = "%s",
                    qvalueCutoff = %s,
                    minGSSize = %s,
                    maxGSSize = %s,
                    readable = FALSE
                  )
                  ',
                  input$order_by,
                  pathway_database,
                  input$kegg_organism,
                  input$go_ont,
                  input$p_value_cutoff,
                  input$p_adjust_method,
                  input$q_value_cutoff,
                  input$gene_set_size[1],
                  input$gene_set_size[2]
                  )
                  enrich_pathways_code(code)
                  }
            })
          }
        }
      })

      output$enriched_pathways_go <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways()@enrichment_go_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_kegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways()@enrichment_kegg_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_reactome <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways()@enrichment_reactome_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_hmdb <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways()@enrichment_hmdb_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enrichment_pathway_metkegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways()@enrichment_metkegg_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_object <- renderText({
        req(enriched_pathways())
        enriched_pathways <- enriched_pathways()
        captured_output1 <- capture.output(enriched_pathways,
                                           type = "message")
        captured_output2 <- capture.output(enriched_pathways,
                                           type = "output")
        captured_output <-
          c(captured_output1,
            captured_output2)
        paste(captured_output, collapse = "\n")
      })

      output$download_enriched_pathways_go <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways_go.csv"
          },
          content = function(file) {
            write.csv(enriched_pathways()@enrichment_go_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_go")
        } else {
          if (length(enriched_pathways()@enrichment_go_result) == 0) {
            shinyjs::disable("download_enriched_pathways_go")
          } else{
            shinyjs::enable("download_enriched_pathways_go")
          }
        }
      })

      output$download_enriched_pathways_kegg <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways_kegg.csv"
          },
          content = function(file) {
            write.csv(enriched_pathways()@enrichment_kegg_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_kegg")
        } else {
          if (length(enriched_pathways()@enrichment_kegg_result) == 0) {
            shinyjs::disable("download_enriched_pathways_kegg")
          } else{
            shinyjs::enable("download_enriched_pathways_kegg")
          }
        }
      })

      output$download_enriched_pathways_reactome <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways_reactome.csv"
          },
          content = function(file) {
            write.csv(enriched_pathways()@enrichment_reactome_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_reactome")
        } else {
          if (length(enriched_pathways()@enrichment_reactome_result) == 0) {
            shinyjs::disable("download_enriched_pathways_reactome")
          } else{
            shinyjs::enable("download_enriched_pathways_reactome")
          }
        }
      })

      output$download_enriched_pathways_hmdb <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways_hmdb.csv"
          },
          content = function(file) {
            write.csv(enriched_pathways()@enrichment_hmdb_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_hmdb")
        } else {
          if (length(enriched_pathways()@enrichment_hmdb_result) == 0) {
            shinyjs::disable("download_enriched_pathways_hmdb")
          } else{
            shinyjs::enable("download_enriched_pathways_hmdb")
          }
        }
      })

      output$download_enriched_pathways_metkegg <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways_metkegg.csv"
          },
          content = function(file) {
            write.csv(enriched_pathways()@enrichment_metkegg_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_metkegg")
        } else {
          if (length(enriched_pathways()@enrichment_metkegg_result) == 0) {
            shinyjs::disable("download_enriched_pathways_metkegg")
          } else{
            shinyjs::enable("download_enriched_pathways_metkegg")
          }
        }
      })

      output$download_enriched_pathways_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways.rda"
          },
          content = function(file) {
            enriched_pathways <-
              enriched_pathways()
            save(enriched_pathways, file = file)
          }
        )
      observe({
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          shinyjs::disable("download_enriched_pathways_object")
        } else {
          shinyjs::enable("download_enriched_pathways_object")
        }
      })

      ### show code
      observeEvent(input$show_enrich_pathways_code, {
        if (is.null(enrich_pathways_code()) ||
            length(enrich_pathways_code()) == 0) {
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
            enrich_pathways_code()
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

      # Go to merge pathways tab ====
      ###if there is not enriched_pathways,
      ###show a warning message
      observeEvent(input$go2merge_pathways, {
        if (is.null(enriched_pathways()) ||
            length(enriched_pathways()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "Please enrich pathways first",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # Navigate to the merge pathways tab
          tab_switch("merge_pathways")
        }
      })

      return(enriched_pathways)
    }
  )
}
