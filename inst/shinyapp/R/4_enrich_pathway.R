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
                       # radioButtons(
                       #   ns("query_type"),
                       #   tags$h4("Query type"),
                       #   choices = c(
                       #     "Gene" = "gene",
                       #     "Metabolite" = "metabolite"
                       #   ),
                       #   inline = TRUE
                       # ),

                       "Organism",
                       verbatimTextOutput(ns("organism")),
                       ## Gene enrichment analysis panel ----
                       shinyjs::hidden(
                         div(id = ns("gene_panel"),
                             radioButtons(
                               ns("analysis_type"),
                               "Analysis type",
                               choices = list(
                                 "Pathway enrichment analysis" = "enrich_pathway",
                                 "Gene set enrichment analysis" = "do_gsea"
                               )
                             ),
                             shinyjs::hidden(
                               div(id = ns("gsea_order_by_panel"),
                                     selectInput(
                                     ns("order_by"),
                                     "Order by",
                                     choices = c(
                                       "Fold change" = "fc",
                                       "Adjusted p value" = "p_value_adjust"
                                     )
                                     )
                               )
                             ),
                             checkboxGroupInput(
                               ns("pathway_database"),
                               "Available Database",
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
                                   #   selectInput(
                                   #     ns("go_orgdb"),
                                   #     "GO Organism",
                                   #     choices = c(
                                   #       "Human (org.Hs.eg.db)" = "org.Hs.eg.db",
                                   #       "Mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
                                   #       "Rat (org.Rn.eg.db)" = "org.Rn.eg.db",
                                   #       "Fly (org.Dm.eg.db)" = "org.Dm.eg.db",
                                   #       "Zebrafish (org.Dr.eg.db)" = "org.Dr.eg.db",
                                   #       "Arabidopsis (org.At.tair.db)" = "org.At.tair.db",
                                   #       "Yeast (org.Sc.sgd.db)" = "org.Sc.sgd.db",
                                   #       "Worm (org.Ce.eg.db)" = "org.Ce.eg.db",
                                   #       "Pig (org.Ss.eg.db)" = "org.Ss.eg.db",
                                   #       "Bovine (org.Bt.eg.db)" = "org.Bt.eg.db",
                                   #       "Rhesus (org.Mmu.eg.db)" = "org.Mmu.eg.db",
                                   #       "Canine (org.Cf.eg.db)" = "org.Cf.eg.db",
                                   #       "E. coli strain K12(org.EcK12.eg.db)" = "org.EcK12.eg.db",
                                   #       "E coli strain Sakai" = "org.EcSakai.eg.db",
                                   #       "Chicken (org.Gg.eg.db)" = "org.Gg.eg.db",
                                   #       "Xenopus (org.Xl.eg.db)" = "org.Xl.eg.db",
                                   #       "Chimp (org.Pt.eg.db)" = "org.Pt.eg.db",
                                   #       "Anopheles (org.Ag.eg.db)" = "org.Ag.eg.db",
                                   #       "Malaria (org.Pf.plasmo.db)" = "org.Pf.plasmo.db",
                                   #       "Myxococcus xanthus DK 1622" = "org.Mxanthus.db"
                                   #     ),
                                   #     selected = "org.Hs.eg.db"
                                   #     ),
                                   # helpText("Enter the name of an OrgDb package that is installed on your system.",
                                   #          "Common examples: org.Hs.eg.db (Human), org.Mm.eg.db (Mouse), org.Rn.eg.db (Rat)",
                                   #          "For the current list of OrgDb packages, visit: ",
                                   #          tags$a(
                                   #            href = "https://bioconductor.org/packages/release/BiocViews.html#___OrgDb",
                                   #            "Bioconductor OrgDb packages",
                                   #            target = "_blank"
                                   #          )),
                                   # textInput(
                                   #   ns("go_keytype"),
                                   #   "GO Keytype",
                                   #   value = "ENTREZID"
                                   # ),
                                   selectInput(
                                     ns("go_keytype"),
                                     "GO Keytype",
                                     choices = c(
                                       "Entrez Gene ID"  = "ENTREZID",
                                       "Gene Symbol"     = "SYMBOL",
                                       "Ensembl Gene ID" = "ENSEMBL",
                                       "UniProt ID"      = "UNIPROT"
                                     ),
                                     selected = "ENTREZID"
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
                                   #   selectInput(
                                   #   ns("kegg_organism"),
                                   #   "KEGG Organism",
                                   #   choices = c(
                                   #     "Human (hsa)" = "hsa",
                                   #     "Mouse (mmu)" = "mmu",
                                   #     "Rat (rno)" = "rno",
                                   #     "Fly (dme)" = "dme",
                                   #     "Zebrafish (dre)" = "dre",
                                   #     "Arabidopsis (ath)" = "ath",
                                   #     "Yeast (sce)" = "sce",
                                   #     "Worm (cel)" = "cel",
                                   #     "Pig (ssc)" = "ssc",
                                   #     "Bovine (bta)" = "bta",
                                   #     "Rhesus (mcc)" = "mcc",
                                   #     "Canine (cfa)" = "cfa",
                                   #     "E. coli K-12 (eco)" = "eco",
                                   #     "E. coli Sakai (ecs)" = "ecs",
                                   #     "Chicken (gga)" = "gga",
                                   #     "Xenopus (xla)" = "xla",
                                   #     "Chimp (ptr)" = "ptr",
                                   #     "Anopheles (aga)" = "aga",
                                   #     "Malaria (pfa)" = "pfa",
                                   #     "Myxococcus xanthus (mxa)" = "mxa"
                                   #   ),
                                   #   selected = "hsa"
                                   #   ),
                                   # helpText(
                                   #   "The KEGG organism code is a three or four letter abbreviation.",
                                   #   "Examples: 'hsa' (Human), 'mmu' (Mouse), 'rno' (Rat).",
                                   #   "For a complete list of organism codes, visit: ",
                                   #   tags$a(
                                   #     href = "https://www.genome.jp/kegg/catalog/org_list.html",
                                   #     "KEGG Organism Codes",
                                   #     target = "_blank"
                                   #   )
                                   # ),
                                   selectInput(
                                     ns("kegg_keytype"),
                                     "KEGG Keytype",
                                     choices = c(
                                       "Entrez" = "kegg",
                                       "NCBI Gene ID" = "ncbi-geneid",
                                       "NCBI Protein ID" = "ncbi-proteinid",
                                       "UniProt" = "uniprot"
                                     ),
                                     selected = "kegg"
                                   )
                               )
                             ),

                             # Reactome specific parameters
                             # shinyjs::hidden(
                             #   div(id = ns("reactome_panel"),
                             #       h4("Reactome Parameters"),
                             #       selectInput(
                             #         ns("reactome_organism"),
                             #         "Reactome Organism",
                             #         choices = c(
                             #           "Human" = "human",
                             #           "Rat" = "rat",
                             #           "Mouse" = "mouse",
                             #           "C. elegans" = "celegans",
                             #           "Yeast" = "yeast",
                             #           "Zebrafish" = "zebrafish",
                             #           "Fruit fly" = "fly",
                             #           "Bovine" = "bovine",
                             #           "Canine" = "canine",
                             #           "Chicken" = "chicken"
                             #         ),
                             #         selected = "human"
                             #       )
                             #   )
                             # ),

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
                               choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
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
                             # selectInput(
                             #   ns("met_organism"),
                             #   "Organism",
                             #   choices = c(
                             #     "Human (hsa)"       = "hsa",
                             #     "Mouse (mmu)"       = "mmu",
                             #     "Rat (rno)"         = "rno",
                             #     "Fly (dme)"         = "dme",
                             #     "Zebrafish (dre)"   = "dre",
                             #     "Yeast (sce)"       = "sce",
                             #     "Worm (cel)"        = "cel",
                             #     "Pig (ssc)"         = "ssc",
                             #     "Bovine (bta)"      = "bta",
                             #     "Canine (cfa)"      = "cfa"
                             #   ),
                             #   selected = "hsa"
                             # ),
                             # helpText(
                             #   "The KEGG organism code is a three or four letter abbreviation.",
                             #   "Examples: 'hsa' (Human), 'mmu' (Mouse), 'rno' (Rat).",
                             #   "For a complete list of organism codes, visit: ",
                             #   tags$a(
                             #     href = "https://www.genome.jp/kegg/catalog/org_list.html",
                             #     "KEGG Organism Codes",
                             #     target = "_blank"
                             #   )
                             # ),
                             checkboxGroupInput(
                               ns("met_pathway_database"),
                               "Available Database",
                               choices = c(
                                 "HMDB" = "hmdb",
                                 "KEGG" = "metkegg"
                               ),
                               selected = NULL
                             ),
                             # checkboxInput(
                             #   ns("use_internal_data"),
                             #   "Use internal database",
                             #   value = TRUE
                             # ),
                             numericInput(
                               ns("p_value_cutoff"),
                               "P-value cutoff",
                               value = 0.05,
                               min = 0,
                               max = 0.5),
                             selectInput(
                               ns("p_adjust_method"),
                               "P-adjust method",
                               choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
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
                       shinyjs::hidden(
                         div(
                           id = ns("gene_enrichment_res_table"),
                           tabsetPanel(
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
                             title = "R object",
                             verbatimTextOutput(ns("gene_enriched_pathways_object")),
                             br(),
                             shinyjs::useShinyjs(),
                             downloadButton(ns("gene_download_enriched_pathways_object"),
                                            "Download",
                                            class = "btn-primary",
                                            style = "background-color: #d83428; color: white;")
                           ))
                         )
                       ),
                       shinyjs::hidden(
                         div(
                           id = ns("met_enrichment_res_table"),
                           tabsetPanel(
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
                               title = "R object",
                               verbatimTextOutput(ns("met_enriched_pathways_object")),
                               br(),
                               shinyjs::useShinyjs(),
                               downloadButton(ns("met_download_enriched_pathways_object"),
                                              "Download",
                                              class = "btn-primary",
                                              style = "background-color: #d83428; color: white;")
                             )
                           )
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
#' @param processed_info Reactive variable containing uploaded variable info.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleElement disable enable useShinyjs
#' @importFrom clusterProfiler enrich_pathway do_gsea
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA enrichPathway
#' @noRd

enrich_pathway_server <- function(id, processed_info, enriched_pathways, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      ## Perform gene enrichment analysis ----
      # query_type <- reactive({ processed_info$query_type })

      variable_info <- reactive({ processed_info$variable_info })

      gene_params <- reactiveValues(
        db_choices = NULL,
        go.orgdb = NULL,
        kegg.organism = NULL,
        reactome.organism = NULL
      )

      met_params <- reactiveValues(
        db_choices = NULL,
        met_organism = NULL
      )

      skip_merge <- reactive({
        processed_info$query_type == "metabolite" &&
          !is.null(met_params$met_organism) &&
          met_params$met_organism != "hsa"
      })

      observeEvent(processed_info$organism, {
        enriched_pathways$query_type <- processed_info$query_type
        if (enriched_pathways$query_type == "gene") {
          enriched_pathways$organism <- processed_info$organism
          orgdb <- processed_info$organism
          gene_params$go.orgdb <- orgdb
          gene_params$kegg.organism <- unname(org2kegg[orgdb])
          gene_params$reactome.organism <- unname(org2react[orgdb])

          gene_params$db_choices <- c("GO" = "go")
          if (!is.na(gene_params$kegg.organism)) {gene_params$db_choices <- c(gene_params$db_choices, "KEGG" = "kegg")}
          if (!is.na(gene_params$reactome.organism)) {gene_params$db_choices <- c(gene_params$db_choices, "Reactome" = "reactome")}

          updateCheckboxGroupInput(
            session, "pathway_database",
            choices  = gene_params$db_choices,
            selected = gene_params$db_choices
          )
        } else if (enriched_pathways$query_type == "metabolite") {
          met_params$met_organism <- processed_info$organism

          if (met_params$met_organism == "hsa") {
            updateCheckboxGroupInput(
              session, "met_pathway_database",
              choices  = c("HMDB" = "hmdb", "KEGG" = "metkegg"),
              selected = c("hmdb", "metkegg")
            )
          } else {
            updateCheckboxGroupInput(
              session, "met_pathway_database",
              choices  = c("KEGG" = "metkegg"),
              selected = "metkegg"
            )
          }
        }
      })

      output$organism <- renderText({
        req(processed_info)
        unname(org_kegg_2name[processed_info$organism])
      })
      outputOptions(output, "organism",
                    suspendWhenHidden = FALSE,  # keep it running
                    priority           = 100)


      observe({
        req(processed_info$query_type)
        shinyjs::toggleElement(
          id = "gene_panel",
          condition = processed_info$query_type == "gene"
        )
        shinyjs::toggleElement(
          id = "metabolite_panel",
          condition = processed_info$query_type == "metabolite"
        )
        shinyjs::toggleElement(
          id = "gene_enrichment_res_table",
          condition = processed_info$query_type == "gene"
        )
        shinyjs::toggleElement(
          id = "met_enrichment_res_table",
          condition = processed_info$query_type == "metabolite"
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
        req(input$pathway_database)
        shinyjs::toggleElement(
          id = "go_panel",
          condition = "go" %in% input$pathway_database
        )
        shinyjs::toggleElement(
          id = "kegg_panel",
          condition = "kegg" %in% input$pathway_database
        )
        shinyjs::toggleElement(
          id = "reactome_panel",
          condition = "reactome" %in% input$pathway_database
        )
      })

      ## Perform metabolite enrichment analysis ----

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
          withProgress(message = 'Analysis in progress...', {
            result <- tryCatch({
              library(clusterProfiler)
              library(ReactomePA)

              # Extract common parameters
              common_params <- list(
                variable_info = variable_info(),
                query_type = processed_info$query_type,
                save_to_local = FALSE,
                pvalueCutoff = input$p_value_cutoff,
                pAdjustMethod = input$p_adjust_method
              )
              if (input$analysis_type == "enrich_pathway") {
                # Add gene-specific parameters if query type is gene
                if (processed_info$query_type == "gene") {
                  enriched_pathways$available_db <- input$pathway_database
                  common_params$database <- input$pathway_database
                  common_params$qvalueCutoff <- input$q_value_cutoff
                  # GO parameters
                  if ("go" %in% common_params$database) {
                    common_params$go.keytype <- input$go_keytype
                    common_params$go.ont <- input$go_ont

                    # Validate input format
                    if (!grepl("^org\\.[A-Za-z]+\\..+\\.db$", gene_params$go.orgdb)) {
                      stop("Invalid OrgDb package name. Expected format: org.XX.eg.db")
                    }
                    # Check if the package is installed
                    if (!requireNamespace(gene_params$go.orgdb, quietly = TRUE)) {
                      stop(paste("Package", gene_params$go.orgdb, "is not installed. Please install it using BiocManager::install('", gene_params$go.orgdb, "')"))
                    }
                    # Load the package
                    requireNamespace(gene_params$go.orgdb)
                    # Get the OrgDb object
                    org_db_obj <- get(gene_params$go.orgdb)
                    common_params$go.orgdb <- org_db_obj
                  }

                  # KEGG parameters
                  if ("kegg" %in% common_params$database) {
                    common_params$kegg.organism <- gene_params$kegg.organism
                    common_params$kegg.keytype <- input$kegg_keytype
                  }

                  # Reactome parameters
                  if ("reactome" %in% common_params$database) {
                    common_params$reactome.organism <- gene_params$reactome.organism
                  }

                  # Gene set size parameters
                  common_params$minGSSize <- input$gene_set_size[1]
                  common_params$maxGSSize <- input$gene_set_size[2]
                }

                if (processed_info$query_type == "metabolite") {
                  enriched_pathways$available_db <- input$met_pathway_database

                  common_params$database <- input$met_pathway_database
                  common_params$met_organism <- met_params$met_organism
                }

                do.call(enrich_pathway, common_params)

              } else if (input$analysis_type == "do_gsea") {
                enriched_pathways$available_db <- input$pathway_database

                common_params$order_by <- input$order_by
                common_params$database <- input$pathway_database
                common_params$qvalueCutoff <- input$q_value_cutoff
                # Gene set size parameters
                common_params$minGSSize <- input$gene_set_size[1]
                common_params$maxGSSize <- input$gene_set_size[2]

                # GO parameters
                if ("go" %in% common_params$database) {
                  common_params$go.keytype <- input$go_keytype
                  common_params$go.ont <- input$go_ont

                  # Validate input format
                  if (!grepl("^org\\.[A-Za-z]+\\..+\\.db$", gene_params$go.orgdb)) {
                    stop("Invalid OrgDb package name. Expected format: org.XX.eg.db")
                  }
                  # Check if the package is installed
                  if (!requireNamespace(gene_params$go.orgdb, quietly = TRUE)) {
                    stop(paste("Package", gene_params$go.orgdb, "is not installed. Please install it using BiocManager::install('", gene_params$go.orgdb, "')"))
                  }
                  # Load the package
                  requireNamespace(gene_params$go.orgdb)
                  # Get the OrgDb object
                  org_db_obj <- get(gene_params$go.orgdb)
                  common_params$go.orgdb <- org_db_obj
                }

                # KEGG parameters
                if ("kegg" %in% common_params$database) {
                  common_params$kegg.organism <- gene_params$kegg.organism
                  common_params$kegg.keytype <- input$kegg_keytype
                }

                # Reactome parameters
                if ("reactome" %in% common_params$database) {
                  common_params$reactome.organism <- gene_params$reactome.organism
                }

                do.call(do_gsea, common_params)
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

            enriched_pathways$enriched_pathways_res <- result


            observe({
              if(input$analysis_type == "do_gsea") {
                updateTextInput(session,
                                "go_keytype",
                                value = "ENTREZID")
                shinyjs::disable("go_keytype")
              } else {
                shinyjs::enable("go_keytype")
              }
            })

            observe({
              if(input$analysis_type == "do_gsea") {
                updateSelectInput(session,
                                  "kegg_keytype",
                                  selected = "kegg")
                shinyjs::disable("kegg_keytype")
              } else {
                shinyjs::enable("kegg_keytype")
              }
            })

            # shinyjs::hide("loading")

            # Save code for reproducibility
            # pathway_database <-
            #   paste0("c(", paste(unlist(
            #     lapply(paste(input$pathway_database), function(x)
            #       paste0('"', x, '"'))
            #   ),
            #   collapse = ", "), ")")

            if (input$analysis_type == "enrich_pathway") {
              if (processed_info$query_type == "gene") {
                pathway_database <-
                  paste0("c(", paste(unlist(
                    lapply(paste(input$pathway_database), function(x)
                      paste0('"', x, '"'))
                  ),
                  collapse = ", "), ")")

                # Build parameter parts based on selected databases
                go_params <- ""
                if ("go" %in% input$pathway_database) {
                  go_params <- sprintf(
                    ' go.orgdb = "%s",
                      go.keytype = "%s",
                      go.ont = "%s",
                      go.universe = NULL,
                      go.pool = FALSE,
                      ',
                    gene_params$go.orgdb,
                    input$go_keytype,
                    input$go_ont
                  )}

                kegg_params <- ""
                if ("kegg" %in% input$pathway_database) {
                  kegg_params <- sprintf(
                    ' kegg.organism = "%s",
                      kegg.keytype = "%s",
                      kegg.universe = NULL,
                      ',
                    gene_params$kegg.organism,
                    input$kegg_keytype
                  )}

                reactome_params <- ""
                if ("reactome" %in% input$pathway_database) {
                  reactome_params <- sprintf(
                    ' reactome.organism = "%s",
                      reactome.universe = NULL,',
                    gene_params$reactome.organism
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
                  processed_info$query_type,
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
                pathway_database <-
                  paste0("c(", paste(unlist(
                    lapply(paste(input$met_pathway_database), function(x)
                      paste0('"', x, '"'))
                  ),
                  collapse = ", "), ")")

                code <- sprintf(
                  '
                  enriched_pathways <-
                    enrich_pathway(
                      variable_info,
                      query_type = "%s",
                      met_organism = "%s",
                      database = %s,
                      pvalueCutoff = %s,
                      pAdjustMethod = "%s"
                    )
                    ',
                  processed_info$query_type,
                  met_params$met_organism,
                  pathway_database,
                  input$p_value_cutoff,
                  input$p_adjust_method
                )
                enrich_pathways_code(code)
              }
            } else { # GSEA code
              pathway_database <-
                paste0("c(", paste(unlist(
                  lapply(paste(input$pathway_database), function(x)
                    paste0('"', x, '"'))
                ),
                collapse = ", "), ")")

              # Build parameter parts based on selected databases
              go_params <- ""
              if ("go" %in% input$pathway_database) {
                go_params <- sprintf(
                  '   go.orgdb = "%s",
                      go.keytype = "ENTREZID",
                      go.ont = "%s",
                  ',
                  gene_params$go.orgdb,
                  input$go_ont
                )}

              kegg_params <- ""
              if ("kegg" %in% input$pathway_database) {
                kegg_params <- sprintf(
                  '   kegg.organism = "%s",
                      kegg.keytype = "kegg",
                  ',
                  gene_params$kegg.organism
                )}

              reactome_params <- ""
              if ("reactome" %in% input$pathway_database) {
                reactome_params <- sprintf(
                  '   reactome.organism = "%s",',
                  gene_params$reactome.organism
                )}

              code <- sprintf(
                '
                enriched_pathways <-
                  do_gsea(
                    variable_info,
                    order_by = "%s",
                    database = %s,%s%s%s
                    pvalueCutoff = %s,
                    pAdjustMethod = "%s",
                    qvalueCutoff = %s,
                    minGSSize = %s,
                    maxGSSize = %s,
                    save_to_local = FALSE
                    )
                    ',
                input$order_by,
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
            }
          })

        }
      })

      output$enriched_pathways_go <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways$enriched_pathways_res@enrichment_go_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_kegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways$enriched_pathways_res@enrichment_kegg_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_reactome <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways$enriched_pathways_res@enrichment_reactome_result@result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_hmdb <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways$enriched_pathways_res@enrichment_hmdb_result@result |>
              dplyr::rename(description = describtion),
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$enriched_pathways_metkegg <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_pathways$enriched_pathways_res@enrichment_metkegg_result@result |>
              dplyr::rename(description = describtion),
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$gene_enriched_pathways_object <- renderText({
        req(enriched_pathways$enriched_pathways_res)
        enriched_pathways_res <- enriched_pathways$enriched_pathways_res
        captured_output1 <- capture.output(enriched_pathways_res,
                                           type = "message")
        captured_output2 <- capture.output(enriched_pathways_res,
                                           type = "output")
        captured_output <-
          c(captured_output1,
            captured_output2)
        paste(captured_output, collapse = "\n")
      })
      output$met_enriched_pathways_object <- renderText({
        req(enriched_pathways$enriched_pathways_res)
        enriched_pathways_res <- enriched_pathways$enriched_pathways_res
        captured_output1 <- capture.output(enriched_pathways_res,
                                           type = "message")
        captured_output2 <- capture.output(enriched_pathways_res,
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
            write.csv(enriched_pathways$enriched_pathways_res@enrichment_go_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("download_enriched_pathways_go")
        } else {
          if (length(enriched_pathways$enriched_pathways_res@enrichment_go_result) == 0) {
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
            write.csv(enriched_pathways$enriched_pathways_res@enrichment_kegg_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("download_enriched_pathways_kegg")
        } else {
          if (length(enriched_pathways$enriched_pathways_res@enrichment_kegg_result) == 0) {
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
            write.csv(enriched_pathways$enriched_pathways_res@enrichment_reactome_result@result,
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("download_enriched_pathways_reactome")
        } else {
          if (length(enriched_pathways$enriched_pathways_res@enrichment_reactome_result) == 0) {
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
            write.csv(enriched_pathways$enriched_pathways_res@enrichment_hmdb_result@result |> dplyr::rename(description = describtion),
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("download_enriched_pathways_hmdb")
        } else {
          if (length(enriched_pathways$enriched_pathways_res@enrichment_hmdb_result) == 0) {
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
            write.csv(enriched_pathways$enriched_pathways_res@enrichment_metkegg_result@result |> dplyr::rename(description = describtion),
                      file,
                      row.names = FALSE)
          }
        )

      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("download_enriched_pathways_metkegg")
        } else {
          if (length(enriched_pathways$enriched_pathways_res@enrichment_metkegg_result) == 0) {
            shinyjs::disable("download_enriched_pathways_metkegg")
          } else{
            shinyjs::enable("download_enriched_pathways_metkegg")
          }
        }
      })

      output$gene_download_enriched_pathways_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways.rda"
          },
          content = function(file) {
            enriched_pathways <-
              enriched_pathways$enriched_pathways_res
            save(enriched_pathways, file = file)
          }
        )
      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("gene_download_enriched_pathways_object")
        } else {
          shinyjs::enable("gene_download_enriched_pathways_object")
        }
      })

      output$met_download_enriched_pathways_object <-
        shiny::downloadHandler(
          filename = function() {
            "enriched_pathways.rda"
          },
          content = function(file) {
            enriched_pathways <-
              enriched_pathways$enriched_pathways_res
            save(enriched_pathways, file = file)
          }
        )
      observe({
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
          shinyjs::disable("met_download_enriched_pathways_object")
        } else {
          shinyjs::enable("met_download_enriched_pathways_object")
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
        if (is.null(enriched_pathways$enriched_pathways_res) ||
            length(enriched_pathways$enriched_pathways_res) == 0) {
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
          # tab_switch("merge_pathways")
          if (skip_merge()) {
            tab_switch("embed_cluster_pathways")   # jump over “merge pathways”
          } else {
            tab_switch("merge_pathways")           # normal route
          }
        }
      })
    }
  )
}
