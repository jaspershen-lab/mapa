#' Data Visualization UI Module
#'
#' Internal UI for data visualization including barplot, module similarity network,
#' module information, and relationship network.
#'
#' @param id Module id.
#' @import shiny
#' @importFrom shinyjs hidden toggleElement useShinyjs
#' @importFrom shinyBS bsButton bsPopover
#' @importFrom shinyWidgets colorPickr
#' @noRd

data_visualization_ui <- function(id) {
  ns <- NS(id)
  tabItem(tabName = "data_visualization",
          fluidPage(
            titlePanel("Data Visualization"),
            tabsetPanel(
              ## tab1: Barplot panel ----
              tabPanel(
                title = "Barplot",
                fluidRow(
                  column(4,
                         br(),
                         fluidRow(
                           column(8,
                                  fileInput(inputId = ns("upload_enriched_functional_module"),
                                            label = tags$span("Upload functional module",
                                                              shinyBS::bsButton(ns("upload_functional_module_info"),
                                                                                label = "",
                                                                                icon = icon("info"),
                                                                                style = "info",
                                                                                size = "extra-small")),
                                            accept = ".rda"),
                                  bsPopover(
                                    id = ns("upload_functional_module_info"),
                                    title = "",
                                    content = "You can upload the functional module file here for data visualization only.",
                                    placement = "right",
                                    trigger = "hover",
                                    options = list(container = "body")
                                  )
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("barplot_level"),
                                    "Level",
                                    choices = c(
                                      "FM" = "functional_module",
                                      "Module" = "module",
                                      "Pathway" = "pathway"),
                                    selected = "pathway")
                           )
                         ),
                         fluidRow(
                           # ### level selection panel for gene ----
                           # shinyjs::hidden(
                           #   div(
                           #     id = ns("level_panel_for_gene"),
                           #     column(4,
                           #            selectInput(
                           #              inputId = ns("barplot_level_gene"),
                           #              label = "Level",
                           #              choices = c(
                           #                "FM" = "functional_module",
                           #                "Module" = "module",
                           #                "Pathway" = "pathway"),
                           #              selected = "functional_module")
                           #     )
                           #   )
                           # ),
                           # ### level selection panel for metabolite ----
                           # shinyjs::hidden(
                           #   div(
                           #     id = ns("level_panel_for_metabolite"),
                           #     column(4,
                           #            selectInput(
                           #              inputId = ns("barplot_level_metabolite"),
                           #              label = "Level",
                           #              choices = c("Pathway" = "pathway")
                           #            )
                           #            )
                           #   )
                           # ),

                           column(4,
                                  numericInput(
                                    inputId = ns("barplot_top_n"),
                                    label = "Top N",
                                    value = 5,
                                    min = 1,
                                    max = 1000)
                           ),
                           column(4,
                                  selectInput(
                                    ns("line_type"),
                                    "Line type",
                                    choices = c(
                                      "Straight" = "straight",
                                      "Meteor" = "meteor"))
                           ),
                           # column(4,
                           #        checkboxInput(ns("barplot_llm_text"), "LLM text", FALSE)
                           # ),
                           column(
                             4,
                             tags$div(
                               style = "display: flex; flex-direction: column;",
                               tags$label("LLM text", `for` = ns("barplot_llm_text")),
                               checkboxInput(ns("barplot_llm_text"), "", FALSE)
                             )
                           )
                         ),
                         fluidRow(
                           column(6,
                                  numericInput(
                                    ns("barplot_y_lable_width"),
                                    "Y label width",
                                    value = 50,
                                    min = 20,
                                    max = 100)
                           ),
                           column(6,
                                  numericInput(ns("barplot_count_cutoff"),
                                               "Count cutoff",
                                               value = 5,
                                               min = 1,
                                               max = 1000)
                           )
                         ),
                         fluidRow(
                           column(6,
                                  numericInput(ns("barplot_p_adjust_cutoff"),
                                               "P-adjust cutoff",
                                               value = 0.05,
                                               min = 0,
                                               max = 0.5)),
                           column(6,
                                  selectInput(
                                    ns("x_axis_name"),
                                    "X axis name",
                                    choices = NULL
                                  ))
                         ),

                         ### db_color_panel_for_gene ----
                         shinyjs::hidden(
                           div(
                             id = ns("db_color_panel_gene"),
                             fluidRow(
                               column(8,
                                      selectInput(ns("gene_barplot_database"),
                                                         "Database",
                                                         choices = c(
                                                           "GO" = "go",
                                                           "KEGG" = "kegg",
                                                           "Reactome" = "reactome"
                                                         ),
                                                  selected = NULL,
                                                  multiple = TRUE)
                               ) #,
                               # column(4,
                               #        checkboxInput(ns("barplot_translation"), "Translation", FALSE)
                               # )
                             ),
                             h4("Database color"),
                             fluidRow(
                               column(4,
                                      shinyWidgets::colorPickr(
                                        inputId = ns("barplot_go_color"),
                                        label = "GO",
                                        selected = "#1F77B4FF",
                                        theme = "monolith",
                                        width = "100%")
                               ),
                               column(4,
                                      shinyWidgets::colorPickr(
                                        inputId = ns("barplot_kegg_color"),
                                        label = "KEGG",
                                        selected = "#FF7F0EFF",
                                        theme = "monolith",
                                        width = "100%")
                               ),
                               column(4,
                                      shinyWidgets::colorPickr(
                                        inputId = ns("barplot_reactome_color"),
                                        label = "Reactome",
                                        selected = "#2CA02CFF",
                                        theme = "monolith",
                                        width = "100%")
                               )
                             )
                           )
                         ),

                         ### db_color_panel_for_metabolite ----
                         shinyjs::hidden(
                           div(
                             id = ns("db_color_panel_metabolite"),
                             fluidRow(
                               column(8,
                                      selectInput(ns("met_barplot_database"),
                                                  "Database",
                                                  choices = c("HMDB" = "hmdb", "KEGG" = "metkegg"),
                                                  selected = NULL,
                                                  multiple = TRUE)
                               ) #,
                               # column(4,
                               #        checkboxInput(ns("barplot_translation"), "Translation", FALSE)
                               # )
                             ),
                             h4("Database color"),
                             fluidRow(
                               column(4,
                                      shinyWidgets::colorPickr(
                                        inputId = ns("barplot_hmdb_color"),
                                        label = "HMDB",
                                        selected = "#9467BDFF",
                                        theme = "monolith",
                                        width = "100%")
                               ),
                               column(4,
                                      shinyWidgets::colorPickr(
                                        inputId = ns("barplot_metkegg_color"),
                                        label = "KEGG",
                                        selected = "#FF7F0EFF",
                                        theme = "monolith",
                                        width = "100%")
                               )
                             )
                           )
                         ),

                         fluidRow(
                           column(4,
                                  selectInput(ns("barplot_type"), "Type",
                                              choices = c("pdf", "png", "jpeg"))
                           ),
                           column(4,
                                  numericInput(ns("barplot_width"), "Width",
                                               value = 7, min = 4, max = 20)

                           ),
                           column(4,
                                  numericInput(ns("barplot_height"), "Height",
                                               value = 7, min = 4, max = 20)

                           )
                         ),
                         fluidRow(
                           column(12,
                                  actionButton(ns("generate_barplot"),
                                               "Generate plot",
                                               class = "btn-primary",
                                               style = "background-color: #d83428; color: white;"),
                                  downloadButton(ns("download_barplot"),
                                                 "Download",
                                                 class = "btn-primary",
                                                 style = "background-color: #d83428; color: white;")
                           )
                         ),
                         br(),
                         fluidRow(
                           column(12,
                                  actionButton(
                                    ns("go2results_1"),
                                    "Next",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;"),
                                  actionButton(
                                    ns("show_barplot_code"),
                                    "Code",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;")
                           )
                         ),
                         style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
                  column(8,
                         br(),
                         shiny::plotOutput(ns("barplot"))
                  )
                )
              ),

              ## tab2: Module Similarity Network panel ----
              tabPanel(
                title = "Module Similarity Network",
                fluidRow(
                  column(4,
                         br(),
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("module_similarity_network_database"),
                                    "Database",
                                    choices = c("GO" = "go",
                                                "KEGG" = "kegg",
                                                "Reactome" = "reactome"),
                                    selected = "go")
                           ),
                           column(6,
                                  numericInput(
                                    ns("module_similarity_network_degree_cutoff"),
                                    "Degree cutoff",
                                    value = 0,
                                    min = 0,
                                    max = 1000)
                           )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("module_similarity_network_level"),
                                    "Level",
                                    choices = c(
                                      "FM" = "functional_module",
                                      "Module" = "module"),
                                    selected = "functional_module")
                           )
                         ),
                         fluidRow(
                           # column(4,
                           #        checkboxInput("module_similarity_network_translation", "Translation", FALSE)
                           # ),
                           column(4,
                                  checkboxInput(ns("module_similarity_network_text"), "Text", FALSE)
                           ),
                           column(4,
                                  checkboxInput(ns("module_similarity_network_llm_text"), "LLM text", FALSE)
                           ),
                           column(4,
                                  checkboxInput(ns("module_similarity_network_text_all"), "Text all", FALSE)
                           )
                         ),
                         fluidRow(
                           column(4,
                                  selectInput(ns("module_similarity_network_type"), "Type",
                                              choices = c("pdf", "png", "jpeg"))
                           ),
                           column(4,
                                  numericInput(ns("module_similarity_network_width"), "Width",
                                               value = 7, min = 4, max = 20)),
                           column(4,
                                  numericInput(ns("module_similarity_network_height"), "Height",
                                               value = 7, min = 4, max = 20))
                         ),
                         fluidRow(
                           column(12,
                                  shinyjs::useShinyjs(),
                                  actionButton(ns("generate_module_similarity_network"),
                                               "Generate plot",
                                               class = "btn-primary",
                                               style = "background-color: #d83428; color: white;"),
                                  shinyjs::useShinyjs(),
                                  downloadButton(ns("download_module_similarity_network"),
                                                 "Download",
                                                 class = "btn-primary",
                                                 style = "background-color: #d83428; color: white;")
                           )
                         ),
                         br(),
                         fluidRow(
                           column(12,
                                  actionButton(
                                    ns("go2results_2"),
                                    "Next",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;"
                                  ),
                                  actionButton(
                                    ns("show_module_similarity_network_code"),
                                    "Code",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;")
                           )
                         ),
                         style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
                  column(8,
                         br(),
                         shiny::plotOutput(ns("module_similarity_network"))
                  )
                )
              ),

              ## tab3: Module information panel ----
              tabPanel(
                title = "Module information",
                fluidRow(
                  column(4,
                         br(),
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("module_information_level"),
                                    "Level",
                                    choices = c(
                                      "FM" = "functional_module",
                                      "Module" = "module"),
                                    selected = "functional_module"
                                    )
                           ),
                           column(6,
                                  selectInput(
                                    ns("module_information_database"),
                                    "Database",
                                    choices = c("GO" = "go",
                                                "KEGG" = "kegg",
                                                "Reactome" = "reactome"),
                                    selected = "go")
                           )
                         ),
                         fluidRow(
                           column(7,
                                  selectInput(
                                    ns("module_information_module_id"),
                                    "Module ID",
                                    choices = NULL)
                           ),
                           column(
                             4,
                             tags$div(
                               style = "display: flex; flex-direction: column;",
                               tags$label("LLM text", `for` = ns("module_information_llm_text")),
                               checkboxInput(ns("module_information_llm_text"), "", FALSE)
                             )
                           )
                           # column(4,
                           #        checkboxInput(ns("module_information_llm_text"), "LLM text", FALSE)
                           # )
                           # column(5,
                           #        checkboxInput("module_information_translation",
                           #                      "Translation", FALSE)
                           # )
                         ),
                         fluidRow(
                           column(4,
                                  selectInput(ns("module_information_type"), "Type",
                                              choices = c("pdf", "png", "jpeg"))
                           ),
                           column(4,
                                  numericInput(ns("module_information_width"), "Width",
                                               value = 7, min = 4, max = 30)),
                           column(4,
                                  numericInput(ns("module_information_height"), "Height",
                                               value = 21, min = 4, max = 20))
                         ),
                         fluidRow(
                           column(12,
                                  shinyjs::useShinyjs(),
                                  actionButton(ns("generate_module_information"),
                                               "Generate plot",
                                               class = "btn-primary",
                                               style = "background-color: #d83428; color: white;"),
                                  downloadButton(ns("download_module_information"),
                                                 "Download",
                                                 class = "btn-primary",
                                                 style = "background-color: #d83428; color: white;")
                           )
                         ),
                         br(),
                         fluidRow(
                           column(12,
                                  shinyjs::useShinyjs(),
                                  actionButton(
                                    ns("go2results_3"),
                                    "Next",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;"),
                                  actionButton(
                                    ns("show_module_information_code"),
                                    "Code",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;")
                           )
                         ),
                         style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
                  column(8,
                         br(),
                         shiny::plotOutput(ns("module_information1")),
                         shiny::plotOutput(ns("module_information2")),
                         shiny::plotOutput(ns("module_information3"))
                  )
                )
              ),

              ## tab4: Relationship network panel ----
              tabPanel(
                title = "Relationship network",
                fluidRow(
                  column(4,
                         br(),
                         fluidRow(
                           column(6,
                                  checkboxInput(ns("relationship_network_circular_plot"),
                                                "Circular layout", FALSE)
                           ),
                           # column(4,
                           #        checkboxInput(ns("relationship_network_filter"),
                           #                      "Filter", FALSE)
                           # ),
                           column(6,
                                  checkboxInput(ns("relationship_network_llm_text"), "LLM text", value = TRUE)
                           )
                           # column(4,
                           #        checkboxInput(ns("relationship_network_translation"),
                           #                      "Translation", FALSE)
                           # )
                         ),
                         fluidRow(
                           column(6,
                                  selectInput(
                                    ns("relationship_network_level"),
                                    "Filter Level",
                                    choices = c(
                                      "FM" = "functional_module",
                                      "Module" = "module"),
                                    selected = "functional_module")
                           ),
                           column(6,
                                  selectInput(
                                    ns("relationship_network_module_id"),
                                    "Module ID",
                                    choices = NULL,
                                    multiple = TRUE,
                                    selected = NULL)
                           )
                         ),
                         h4("Levels Included"),
                         fluidRow(
                           column(3,
                                  checkboxInput(ns("relationship_network_include_functional_modules"),
                                                label = tags$span("FM",
                                                                  shinyBS::bsButton(ns("functional_module_info"),
                                                                                    label = "",
                                                                                    icon = icon("info"),
                                                                                    style = "info",
                                                                                    size = "extra-small")),
                                                value = TRUE)),
                           bsPopover(
                             id = ns("functional_module_info"),
                             title = "",
                             content = "FM is functional module",
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_include_modules"),
                                                "Modules",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_include_pathways"),
                                                "Pathways",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_include_molecules"),
                                                "Molecules",
                                                value = TRUE)
                           )
                         ),
                         h4("Colors"),
                         fluidRow(
                           column(3,
                                  shinyWidgets::colorPickr(
                                    inputId = ns("relationship_network_functional_module_color"),
                                    label = "FM",
                                    selected = "#F05C3BFF",
                                    theme = "monolith",
                                    width = "100%"
                                  )
                           ),
                           column(3,
                                  shinyWidgets::colorPickr(
                                    inputId = ns("relationship_network_module_color"),
                                    label = "Module",
                                    selected = "#46732EFF",
                                    theme = "monolith",
                                    width = "100%"
                                  )
                           ),
                           column(3,
                                  shinyWidgets::colorPickr(
                                    inputId = ns("relationship_network_pathway_color"),
                                    label = "Pathway",
                                    selected = "#197EC0FF",
                                    theme = "monolith",
                                    width = "100%"
                                  )
                           ),
                           column(3,
                                  shinyWidgets::colorPickr(
                                    inputId = ns("relationship_network_molecule_color"),
                                    label = "Molecule",
                                    selected = "#3B4992FF",
                                    theme = "monolith",
                                    width = "100%"
                                  )
                           )
                         ),
                         h4("Text"),
                         fluidRow(
                           column(3,
                                  checkboxInput(ns("relationship_network_functional_module_text"),
                                                "FM",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_module_text"),
                                                "Module",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_pathway_text"),
                                                "Pathway",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_molecule_text"),
                                                "Molecules",
                                                value = TRUE)
                           )
                         ),
                         h4("Text size"),
                         fluidRow(
                           column(3,
                                  numericInput(ns("relationship_network_functional_module_text_size"),
                                               "FM",
                                               value = 3, min = 0.3, max = 10)
                           ),
                           column(3,
                                  numericInput(ns("relationship_network_module_text_size"),
                                               "Module",
                                               value = 3, min = 0.3, max = 10)
                           ),
                           column(3,
                                  numericInput(ns("relationship_network_pathway_text_size"),
                                               "Pathway",
                                               value = 3, min = 0.3, max = 10)
                           ),
                           column(3,
                                  numericInput(ns("relationship_network_molecule_text_size"),
                                               "Molecule",
                                               value = 3, min = 0.3, max = 10)
                           )
                         ),
                         h4("Arrange position"),
                         fluidRow(
                           column(3,
                                  checkboxInput(ns("relationship_network_functional_module_arrange_position"),
                                                "FM",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_module_arrange_position"),
                                                "Module",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_pathway_arrange_position"),
                                                "Pathway",
                                                value = TRUE)
                           ),
                           column(3,
                                  checkboxInput(ns("relationship_network_molecule_arrange_position"),
                                                "Molecules",
                                                value = TRUE)
                           )
                         ),
                         h4("Position limits"),
                         fluidRow(
                           column(6,
                                  sliderInput(
                                    ns("relationship_network_functional_module_position_limits"),
                                    "Functional module",
                                    min = 0, max = 1,
                                    value = c(0, 1))
                           ),
                           column(6,
                                  sliderInput(
                                    ns("relationship_network_module_position_limits"),
                                    "Module",
                                    min = 0, max = 1,
                                    value = c(0, 1))
                           )
                         ),
                         fluidRow(
                           column(6,
                                  sliderInput(
                                    ns("relationship_network_pathway_position_limits"),
                                    "Pathway",
                                    min = 0, max = 1,
                                    value = c(0, 1))
                           ),
                           column(6,
                                  sliderInput(
                                    ns("relationship_network_molecule_position_limits"),
                                    "Molecule",
                                    min = 0, max = 1,
                                    value = c(0, 1))
                           )
                         ),
                         fluidRow(
                           column(4,
                                  selectInput(ns("relationship_network_type"),
                                              "Type",
                                              choices = c("pdf", "png", "jpeg")
                                  )
                           ),
                           column(4,
                                  numericInput(ns("relationship_network_width"),
                                               "Width",
                                               value = 21, min = 4, max = 30)
                           ),
                           column(4,
                                  numericInput(ns("relationship_network_height"),
                                               "Height",
                                               value = 7, min = 4, max = 20)
                           )
                         ),
                         fluidRow(
                           column(12,
                                  actionButton(ns("generate_relationship_network"),
                                               "Generate plot",
                                               class = "btn-primary",
                                               style = "background-color: #d83428; color: white;"),
                                  shinyjs::useShinyjs(),
                                  downloadButton(ns("download_relationship_network"),
                                                 "Download",
                                                 class = "btn-primary",
                                                 style = "background-color: #d83428; color: white;")
                           )
                         ),
                         br(),
                         fluidRow(
                           column(12,
                                  actionButton(
                                    ns("go2results_4"),
                                    "Next",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;"),
                                  actionButton(
                                    ns("show_relationship_network_code"),
                                    "Code",
                                    class = "btn-primary",
                                    style = "background-color: #d83428; color: white;"
                                  )
                           )
                         ),
                         style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
                  column(8,
                         br(),
                         shiny::plotOutput(ns("relationship_network"))
                  )
                )
              )
            )
          ))
}

#' Data Visualization Server Module
#'
#' Internal server logic for data visualization including barplot, module similarity network,
#' module information, and relationship network.
#'
#' @param input,output,session Internal parameters for {shiny}. DO NOT REMOVE.
#' @param id Module id.
#' @param enriched_functional_module Reactive value containing enriched functional module data.
#' @param tab_switch Function to switch tabs.
#' @import shiny
#' @importFrom shinyjs toggleElement useShinyjs
#' @importFrom stringr str_sort
#' @noRd

data_visualization_server <- function(id, enriched_functional_module, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      query_type <- reactive({
        req(enriched_functional_module())
        if ("enrich_pathway" %in% names(enriched_functional_module()@process_info)) {
          enriched_functional_module()@process_info$enrich_pathway@parameter$query_type
        } else {
          enriched_functional_module()@process_info$do_gsea@parameter$query_type
        }
        })

      observe({
        req(query_type())
        shinyjs::toggleElement(
          id = "level_panel_for_gene",
          condition = query_type() == "gene"
        )
        shinyjs::toggleElement(
          id = "db_color_panel_gene",
          condition = query_type() == "gene"
        )
      })

      observe({
        req(query_type())
        shinyjs::toggleElement(
          id = "level_panel_for_metabolite",
          condition = query_type() == "metabolite"
        )
        shinyjs::toggleElement(
          id = "db_color_panel_metabolite",
          condition = query_type() == "metabolite"
        )
      })

      observe({
        req(enriched_functional_module())
        db_choices <- c("GO" = "go", "KEGG" = "kegg", "Reactome" = "reactome", "HMDB" = "hmdb", "KEGG" = "metkegg")
        if ((query_type() == "gene") & ("enrich_pathway" %in% names(enriched_functional_module()@process_info))) {
          all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
          available_db <- enriched_functional_module()@process_info$enrich_pathway@parameter$database
          updateSelectInput(
            session,
            "gene_barplot_database",
            choices = db_choices[db_choices %in% available_db],
            selected = available_db
          )
        } else if ((query_type() == "gene") & ("do_gsea" %in% names(enriched_functional_module()@process_info))) {
          all_choices <- c("NES")
          available_db <- enriched_functional_module()@process_info$do_gsea@parameter$database
          updateSelectInput(
            session,
            "gene_barplot_database",
            choices = db_choices[db_choices %in% available_db],
            selected = available_db
          )
        } else if (query_type() == "metabolite") {
          all_choices <- c("qscore")
          available_db <- enriched_functional_module()@process_info$enrich_pathway@parameter$database
          updateSelectInput(
            session,
            "met_barplot_database",
            choices = db_choices[db_choices %in% available_db],
            selected = available_db
          )
        }

        updateSelectInput(
          session,
          "x_axis_name",
          choices = all_choices
        )
        # Module similarity network
        updateSelectInput(
          session,
          "module_similarity_network_database",
          choices = levels(factor(enriched_functional_module()@merged_module$result_with_module$database))
        )
        # Module information
        updateSelectInput(
          session,
          "module_information_database",
          choices = levels(factor(enriched_functional_module()@merged_module$result_with_module$database))
        )
      })

      observeEvent(input$upload_enriched_functional_module, {
        if (!is.null(input$upload_enriched_functional_module$datapath)) {
          message("Loading data")
          tempEnv <- new.env()
          load(input$upload_enriched_functional_module$datapath,
               envir = tempEnv)

          names <- ls(tempEnv)

          if (length(names) == 1) {
            # If enriched_functional_module is another reactiveVal, uncomment the next line
            enriched_functional_module(get(names[1], envir = tempEnv))
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

        if ((query_type() == "gene") & ("enrich_pathway" %in% names(enriched_functional_module()@process_info))) {
          all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
        } else if ((query_type() == "gene") & ("do_gsea" %in% names(enriched_functional_module()@process_info))) {
          all_choices <- c("NES")
        } else if (query_type() == "metabolite") {
          all_choices <- c("qscore")
        }
        # barplot x_axis
        updateSelectInput(
          session,
          "x_axis_name",
          choices = all_choices
        )
        # Module similarity network
        updateSelectInput(
          session,
          "module_similarity_network_database",
          choices = levels(factor(enriched_functional_module()@merged_module$result_with_module$database))
        )
        # Module information
        updateSelectInput(
          session,
          "module_information_database",
          choices = levels(factor(enriched_functional_module()@merged_module$result_with_module$database))
        )
      })

      observeEvent(enriched_functional_module(), {

        if (query_type() == "gene") {
          if (length(c(enriched_functional_module()@merged_pathway_go,
                       enriched_functional_module()@merged_pathway_kegg,
                       enriched_functional_module()@merged_pathway_reactome)) == 0) {
            updateSelectInput(session, "barplot_level",
                              choices = c("FM" = "functional_module",
                                          "Pathway" = "pathway"),
                              selected = "pathway")
            updateSelectInput(session, "module_similarity_network_level",
                              choices = c("FM" = "functional_module"),
                              selected = "functional_module")
            updateSelectInput(session, "module_information_level",
                              choices = c("FM" = "functional_module"),
                              selected = "functional_module")
            updateSelectInput(session, "relationship_network_level", choices = c("FM" = "functional_module"), selected = "functional_module")
            updateCheckboxInput(session, "relationship_network_include_modules", value = FALSE)
            disable("relationship_network_include_modules")
            updateCheckboxInput(session, "relationship_network_module_text", value = FALSE)
            disable("relationship_network_module_text")
            updateCheckboxInput(session, "relationship_network_module_arrange_position", value = FALSE)
            disable("relationship_network_module_arrange_position")
          }
        } else if (query_type() == "metabolite") {
          if (length(c(enriched_functional_module()@merged_pathway_hmdb,
                       enriched_functional_module()@merged_pathway_metkegg)) == 0) {
            updateSelectInput(session, "barplot_level",
                              choices = c("FM" = "functional_module",
                                          "Pathway" = "pathway"),
                              selected = "pathway")
            updateSelectInput(session, "module_similarity_network_level",
                              choices = c("FM" = "functional_module"),
                              selected = "functional_module")
            updateSelectInput(session, "module_information_level",
                              choices = c("FM" = "functional_module"),
                              selected = "functional_module")
            updateSelectInput(session, "relationship_network_level", choices = c("FM" = "functional_module"), selected = "functional_module")
            updateCheckboxInput(session, "relationship_network_include_modules", value = FALSE)
            disable("relationship_network_include_modules")
            updateCheckboxInput(session, "relationship_network_module_text", value = FALSE)
            disable("relationship_network_module_text")
            updateCheckboxInput(session, "relationship_network_module_arrange_position", value = FALSE)
            disable("relationship_network_module_arrange_position")
          }
        }

        if (length(enriched_functional_module()@llm_module_interpretation) == 0) {
          updateCheckboxInput(session, "barplot_llm_text", value = FALSE)
          disable("barplot_llm_text")
          updateCheckboxInput(session, "module_similarity_network_llm_text", value = FALSE)
          disable("module_similarity_network_llm_text")
          updateCheckboxInput(session, "module_information_llm_text", value = FALSE)
          disable("module_information_llm_text")
          updateCheckboxInput(session, "relationship_network_llm_text", value = FALSE)
          disable("relationship_network_llm_text")
        }

        if (input$module_information_level == "functional_module") {
          updateSelectInput(session, "module_information_database", selected = NULL)
          disable("module_information_database")
        }
      })

      ## Barplot ----
      # Observe generate barplot button click
      barplot <-
        reactiveVal()
      barplot_code <-
        reactiveVal()

      observeEvent(input$generate_barplot, {
        message("generating barplot")
        if (is.null(enriched_functional_module())) {
          # No enriched functional module available
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional module data available. Please complete the previous steps or upload the data",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # shinyjs::show("loading")
          withProgress(message = 'Analysis in progress...', {
            tryCatch({
              if (query_type() == "gene") {
                plot <-
                  plot_pathway_bar(
                    object = enriched_functional_module(),
                    top_n = input$barplot_top_n,
                    x_axis_name = input$x_axis_name,
                    y_label_width = input$barplot_y_lable_width,
                    p.adjust.cutoff = input$barplot_p_adjust_cutoff,
                    count.cutoff = input$barplot_count_cutoff,
                    level = input$barplot_level,
                    llm_text = input$barplot_llm_text,
                    database = input$gene_barplot_database,
                    line_type = input$line_type,
                    database_color = c(
                      GO = input$barplot_go_color,
                      KEGG = input$barplot_kegg_color,
                      Reactome = input$barplot_reactome_color
                    )
                    # translation = input$barplot_translation
                  )
              } else {
                plot <-
                  plot_pathway_bar(
                    object = enriched_functional_module(),
                    top_n = input$barplot_top_n,
                    x_axis_name = input$x_axis_name,
                    y_label_width = input$barplot_y_lable_width,
                    p.adjust.cutoff = input$barplot_p_adjust_cutoff,
                    count.cutoff = input$barplot_count_cutoff,
                    level = input$barplot_level,
                    llm_text = input$barplot_llm_text,
                    database = input$met_barplot_database,
                    line_type = input$line_type,
                    database_color = c(
                      HMDB = input$barplot_hmdb_color,
                      KEGG = input$barplot_metkegg_color
                    )
                    # translation = input$barplot_translation
                  )
              }
              },
              error = function(e) {
                showModal(modalDialog(
                  title = "Error",
                  paste("Details:", e$message),
                  easyClose = TRUE,
                  footer = modalButton("Close")
                ))
              }
            )
          })

          # shinyjs::hide("loading")

          barplot(plot)

          ###save code
          if (query_type() == "gene") {
            data_color <-
              paste0("c(", paste(paste(
                c("GO", "KEGG", "Reactome"),
                c(
                  paste0('"', input$barplot_go_color, '"'),
                  paste0('"', input$barplot_kegg_color, '"'),
                  paste0('"', input$barplot_reactome_color, '"')
                ),
                sep = " = "
              ),
              collapse = ", "), ")")
            barplot_database <-
              paste0("c(", paste(unlist(lapply(paste(input$gene_barplot_database), function(x)
                paste0('"', x, '"'))),
                collapse = ", "), ")")
          } else {
            data_color <-
              paste0("c(", paste(paste(
                c("HMDB", "KEGG"),
                c(
                  paste0('"', input$barplot_hmdb_color, '"'),
                  paste0('"', input$barplot_metkegg_color, '"')
                ),
                sep = " = "
              ),
              collapse = ", "), ")")
            barplot_database <-
              paste0("c(", paste(unlist(lapply(paste(input$met_barplot_database), function(x)
                paste0('"', x, '"'))),
                collapse = ", "), ")")
          }

          barplot_code <-
            sprintf(
              "
            plot_pathway_bar(
            object = enriched_functional_module,
            top_n = %s,
            x_axis_name = %s,
            y_lable_width = %s,
            p.adjust.cutoff = %s,
            count.cutoff = %s,
            level = %s,
            llm_text = %s,
            database = %s,
            line_type = %s,
            database_color = %s)
            ",
              input$barplot_top_n,
              input$x_axis_name,
              input$barplot_y_lable_width,
              input$barplot_p_adjust_cutoff,
              input$barplot_count_cutoff,
              paste0('"', input$barplot_level, '"'),
              input$barplot_llm_text,
              barplot_database,
              paste0('"', input$line_type, '"'),
              data_color
            )

          barplot_code(barplot_code)
        }
      })


      output$barplot <-
        renderPlot({
          req(barplot())
          barplot()
        },
        res = 96)

      # output$barplot <-
      #   renderPlot({
      #     req(barplot())
      #     barplot()
      #   },
      #   width = function() {
      #     input$barplot_width_show
      #   },
      #   height = function() {
      #     input$barplot_height_show
      #   })

      output$download_barplot <-
        downloadHandler(
          filename = function() {
            paste0("pathway_barplot.", input$barplot_type)
          },
          content = function(file) {
            ggsave(
              file,
              plot = barplot(),
              width = input$barplot_width,
              height = input$barplot_height
            )
          }
        )

      observe({
        if (is.null(barplot()) ||
            length(barplot()) == 0) {
          shinyjs::disable("download_barplot")
        } else {
          shinyjs::enable("download_barplot")
        }
      })

      ####show code
      observeEvent(input$show_barplot_code, {
        if (is.null(barplot_code()) ||
            length(barplot_code()) == 0) {
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
            barplot_code()
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



      ## Module similarity network ----
      # Observe generate module_similarity_network button click
      module_similarity_network <-
        reactiveVal()

      module_similarity_network_code <-
        reactiveVal()

      observeEvent(input$generate_module_similarity_network, {
        if (is.null(enriched_functional_module())) {
          # No enriched functional module available
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional module data available. Please complete the previous steps or upload the data.",
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
                  level = input$module_similarity_network_level,
                  database = input$module_similarity_network_database,
                  degree_cutoff = input$module_similarity_network_degree_cutoff,
                  text = input$module_similarity_network_text,
                  llm_text = input$module_similarity_network_llm_text,
                  text_all = input$module_similarity_network_text_all
                  # translation = input$module_similarity_network_translation
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

          # shinyjs::hide("loading")

          module_similarity_network(plot)

          ###save code
          module_similarity_network_code <-
            sprintf(
              '
            plot_similarity_network(
            object = enriched_functional_module,
            level = %s,
            database = %s,
            degree_cutoff = %s,
            text = %s,
            llm_text = %s,
            text_all = %s)
            ',
              paste0('"', input$module_similarity_network_level, '"'),
              paste0('"', input$module_similarity_network_database, '"'),
              input$module_similarity_network_degree_cutoff,
              input$module_similarity_network_text,
              input$module_similarity_network_llm_text,
              input$module_similarity_network_text_all
            )

          module_similarity_network_code(module_similarity_network_code)

        }
      })

      output$module_similarity_network <-
        renderPlot({
          req(module_similarity_network())
          module_similarity_network()
        },
        res = 96)


      ######code for module_similarity_network
      ####show code
      observeEvent(input$show_module_similarity_network_code, {
        if (is.null(module_similarity_network_code()) ||
            length(module_similarity_network_code()) == 0) {
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
            module_similarity_network_code()
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


      output$download_module_similarity_network <-
        downloadHandler(
          filename = function() {
            paste0(
              "module_similarity_network_",
              ifelse(
                input$module_similarity_network_level == "module",
                input$module_similarity_network_database,
                "functional_module"
              ),
              ".",
              input$module_similarity_network_type
            )
          },
          content = function(file) {
            ggsave(
              file,
              plot = module_similarity_network(),
              width = input$module_similarity_network_width,
              height = input$module_similarity_network_height
            )
          }
        )

      observe({
        if (is.null(module_similarity_network()) ||
            length(module_similarity_network()) == 0) {
          shinyjs::disable("download_module_similarity_network")
        } else {
          shinyjs::enable("download_module_similarity_network")
        }
      })




      ## Module information plot ----
      # Update the module ID
      module_information_module_id <-
        reactiveVal()

      observe({
        if (!is.null(enriched_functional_module()) &
            length(enriched_functional_module()) != 0) {
          ####level is functional module
          if (input$module_information_level == "functional_module") {
            if (length(enriched_functional_module()@merged_module) > 0) {
              module_information_module_id <-
                unique(
                  enriched_functional_module()@merged_module$functional_module_result$module
                )
              module_information_module_id(module_information_module_id)
            }
          }

          ####level is module
          if (input$module_information_level == "module") {
            ####database is go
            if (input$module_information_database == "go") {
              if (length(enriched_functional_module()@merged_pathway_go) > 0) {
                module_information_module_id <-
                  unique(
                    enriched_functional_module()@merged_pathway_go$module_result$module
                  )
                module_information_module_id(module_information_module_id)
              }
            }

            ####database is kegg
            if (input$module_information_database == "kegg") {
              if (length(enriched_functional_module()@merged_pathway_kegg) > 0) {
                module_information_module_id <-
                  unique(
                    enriched_functional_module()@merged_pathway_kegg$module_result$module
                  )
                module_information_module_id(module_information_module_id)
              }
            }

            ####database is reactome
            if (input$module_information_database == "reactome") {
              if (length(enriched_functional_module()@merged_pathway_reactome) > 0) {
                module_information_module_id <-
                  unique(
                    enriched_functional_module()@merged_pathway_reactome$module_result$module
                  )
                module_information_module_id(module_information_module_id)
              }
            }

            # metabolite module information ========
            ####database is go
            if (input$module_information_database == "go") {
              if (length(enriched_functional_module()@merged_pathway_go) > 0) {
                module_information_module_id <-
                  unique(
                    enriched_functional_module()@merged_pathway_go$module_result$module
                  )
                module_information_module_id(module_information_module_id)
              }
            }

            ####database is kegg
            if (input$module_information_database == "kegg") {
              if (length(enriched_functional_module()@merged_pathway_kegg) > 0) {
                module_information_module_id <-
                  unique(
                    enriched_functional_module()@merged_pathway_kegg$module_result$module
                  )
                module_information_module_id(module_information_module_id)
              }
            }

          }

          updateSelectInput(
            session,
            "module_information_module_id",
            choices = stringr::str_sort(module_information_module_id(), numeric = TRUE),
            selected = stringr::str_sort(module_information_module_id(), numeric = TRUE)[1]
          )
        }
      })

      # Observe generate module information button click
      module_information <-
        reactiveVal()
      module_information1 <-
        reactiveVal()
      module_information2 <-
        reactiveVal()
      module_information3 <-
        reactiveVal()
      module_information_code <-
        reactiveVal()

      ####if the module_information_module_id() is null, then show warning
      observeEvent(input$generate_module_information, {
        if (is.null(enriched_functional_module())) {
          # No enriched functional module available
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional module data available. Please complete the previous steps or upload the data.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          if (is.null(module_information_module_id())) {
            # No enriched functional module available
            showModal(
              modalDialog(
                title = "Warning",
                "Select a module ID first",
                easyClose = TRUE,
                footer = modalButton("Close")
              )
            )
          } else{
            # shinyjs::show("loading")
            # browser()
            withProgress(message = 'Analysis in progress...', {
              tryCatch(
                plot <-
                  plot_module_info(
                    object = enriched_functional_module(),
                    level = input$module_information_level,
                    llm_text = input$module_information_llm_text,
                    database = input$module_information_database,
                    module_id = input$module_information_module_id
                    # translation = input$module_information_translation
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

            # shinyjs::hide("loading")
            if (is(plot, "ggplot")) {
              plot_all <-
                plot + plot + plot +
                patchwork::plot_layout(ncol = 1)

              module_information(plot_all)
              module_information1(plot)
              module_information2(plot)
              module_information3(plot)
            } else{
              plot_all <-
                plot[[1]] + plot[[2]] + plot[[3]] +
                patchwork::plot_layout(ncol = 1)

              module_information(plot_all)
              module_information1(plot[[1]])
              module_information2(plot[[2]])
              module_information3(plot[[3]])
            }


            ###save code
            module_information_code <-
              sprintf(
                '
            plot_module_info(
            object = enriched_functional_module,
            level = %s,
            llm_text = %s,
            database = %s,
            module_id = %s)
            ',
                paste0('"', input$module_information_level, '"'),
                input$module_information_llm_text,
                paste0('"', input$module_information_database, '"'),
                paste0('"', input$module_information_module_id, '"')
              )

            module_information_code(module_information_code)

          }
        }
      })

      output$module_information <-
        renderPlot({
          req(module_information())
          module_information()
        }, res = 96)

      output$module_information1 <-
        renderPlot({
          req(module_information1())
          module_information1()
        }, res = 96)

      output$module_information2 <-
        renderPlot({
          req(module_information2())
          module_information2()
        }, res = 96)

      output$module_information3 <-
        renderPlot({
          req(module_information3())
          module_information3()
        }, res = 96)

      ######code for module_information
      ####show code
      observeEvent(input$show_module_information_code, {
        if (is.null(module_information_code()) ||
            length(module_information_code()) == 0) {
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
            module_information_code()
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


      output$download_module_information <-
        downloadHandler(
          filename = function() {
            paste0(
              "module_information_",
              input$module_information_module_id,
              ".",
              input$module_information_type
            )
          },
          content = function(file) {
            ggsave(
              file,
              plot = module_information(),
              width = input$module_information_width,
              height = input$module_information_height
            )
          }
        )

      observe({
        if (is.null(module_information1()) ||
            length(module_information1()) == 0) {
          shinyjs::disable("download_module_information")
        } else {
          shinyjs::enable("download_module_information")
        }
      })


      ## Relationship network plot ----
      # Update the module ID
      observeEvent(list(enriched_functional_module(), input$relationship_network_level),{
        if (!is.null(enriched_functional_module()) &
            length(enriched_functional_module()) != 0) {
          ####level is functional module
          if (input$relationship_network_level == "functional_module") {
            if (length(enriched_functional_module()@merged_module) > 0) {
              relationship_network_module_id <-
                unique(
                  enriched_functional_module()@merged_module$functional_module_result$module
                )
            }
          }

          ####level is module
          if (input$relationship_network_level == "module") {
            if (length(enriched_functional_module()@merged_pathway_go) > 0) {
              relationship_network_module_id_go <-
                unique(enriched_functional_module()@merged_pathway_go$module_result$module)
            } else{
              relationship_network_module_id_go <- NULL
            }

            if (length(enriched_functional_module()@merged_pathway_kegg) > 0) {
              relationship_network_module_id_kegg <-
                unique(
                  enriched_functional_module()@merged_pathway_kegg$module_result$module
                )
            } else{
              relationship_network_module_id_kegg <- NULL
            }

            ####database is reactome
            if (length(enriched_functional_module()@merged_pathway_reactome) > 0) {
              relationship_network_module_id_reactome <-
                unique(
                  enriched_functional_module()@merged_pathway_reactome$module_result$module
                )
            } else{
              relationship_network_module_id_reactome <- NULL
            }

            if (length(enriched_functional_module()@merged_pathway_hmdb) > 0) {
              relationship_network_module_id_hmdb <-
                unique(enriched_functional_module()@merged_pathway_hmdb$module_result$module)
            } else{
              relationship_network_module_id_hmdb <- NULL
            }

            if (length(enriched_functional_module()@merged_pathway_metkegg) > 0) {
              relationship_network_module_id_metkegg <-
                unique(
                  enriched_functional_module()@merged_pathway_metkegg$module_result$module
                )
            } else{
              relationship_network_module_id_metkegg <- NULL
            }

            relationship_network_module_id <-
              c(
                relationship_network_module_id_go,
                relationship_network_module_id_kegg,
                relationship_network_module_id_reactome,
                relationship_network_module_id_hmdb,
                relationship_network_module_id_metkegg
              )
          }

          updateSelectInput(
            session,
            "relationship_network_module_id",
            choices = stringr::str_sort(relationship_network_module_id, numeric = TRUE),
            selected = NULL
          )
        }
      })

      # Observe generate relationship network button click
      relationship_network <-
        reactiveVal()
      relationship_network_code <-
        reactiveVal()

      ####get the filtered enriched_functional_module
      object <-
        reactiveVal()

      observeEvent(input$generate_relationship_network, {
        if (is.null(enriched_functional_module())) {
          # No enriched functional module available
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional module data available. Please complete the previous steps or upload the data.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          ####if filtered by functiobal module and modules
          object <-
            enriched_functional_module()
          if (!is.null(input$relationship_network_module_id)) {
            tryCatch({
              object <-
                filter_functional_module(
                  object,
                  level = input$relationship_network_level,
                  remain_id = input$relationship_network_module_id
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
          }

          object(object)

          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch(
              plot <-
                plot_relationship_network(
                  object = object(),
                  include_functional_modules = input$relationship_network_include_functional_modules,
                  include_modules = input$relationship_network_include_modules,
                  include_pathways = input$relationship_network_include_pathways,
                  include_molecules = input$relationship_network_include_molecules,
                  functional_module_text = input$relationship_network_functional_module_text,
                  llm_text = input$relationship_network_llm_text,
                  module_text = input$relationship_network_module_text,
                  pathway_text = input$relationship_network_pathway_text,
                  molecule_text = input$relationship_network_molecule_text,
                  circular_plot = input$relationship_network_circular_plot,
                  functional_module_color = input$relationship_network_functional_module_color,
                  module_color = input$relationship_network_module_color,
                  pathway_color = input$relationship_network_pathway_color,
                  molecule_color = input$relationship_network_molecule_color,
                  functional_module_arrange_position = input$relationship_network_functional_module_arrange_position,
                  module_arrange_position = input$relationship_network_module_arrange_position,
                  pathway_arrange_position = input$relationship_network_pathway_arrange_position,
                  molecule_arrange_position = input$relationship_network_molecule_arrange_position,
                  functional_module_position_limits = c(
                    input$relationship_network_functional_module_position_limits[1],
                    input$relationship_network_functional_module_position_limits[2]
                  ),
                  module_position_limits = c(
                    input$relationship_network_module_position_limits[1],
                    input$relationship_network_module_position_limits[2]
                  ),
                  pathway_position_limits = c(
                    input$relationship_network_pathway_position_limits[1],
                    input$relationship_network_pathway_position_limits[2]
                  ),
                  molecule_position_limits = c(
                    input$relationship_network_molecule_position_limits[1],
                    input$relationship_network_molecule_position_limits[2]
                  )
                  # translation = input$relationship_network_translation
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

          # shinyjs::hide("loading")

          relationship_network(plot)

          ###save code
          relationship_network_module_id <-
            paste0("c(",
                   paste0(
                     paste0('"',
                            input$relationship_network_module_id,
                            '"'),
                     collapse = ", "
                   ),
                   ")")

          relationship_network_code1 <-
            sprintf(
              '
            object <-
            filter_functional_module(
              object,
              level = %s,
              remain_id = %s
            )
            ',
              paste0('"', input$relationship_network_level, '"'),
              relationship_network_module_id
            )

          functional_module_position_limits <-
            paste0(
              "c(",
              paste0(
                input$relationship_network_functional_module_position_limits,
                collapse = ", "
              ),
              ")"
            )

          module_position_limits <-
            paste0(
              "c(",
              paste0(
                input$relationship_network_module_position_limits,
                collapse = ", "
              ),
              ")"
            )

          pathway_position_limits <-
            paste0(
              "c(",
              paste0(
                input$relationship_network_pathway_position_limits,
                collapse = ", "
              ),
              ")"
            )

          molecule_position_limits <-
            paste0(
              "c(",
              paste0(
                input$relationship_network_molecule_position_limits,
                collapse = ", "
              ),
              ")"
            )

          relationship_network_code2 <-
            sprintf(
              '
            plot_relationship_network(
            object = object,
            include_functional_modules = %s,
            include_modules = %s,
            include_pathways = %s,
            include_molecules = %s,
            functional_module_text = %s,
            llm_text = %s,
            module_text = %s,
            pathway_text = %s,
            molecule_text = %s,
            circular_plot = %s,
            functional_module_color = %s,
            module_color = %s,
            pathway_color = %s,
            molecule_color = %s,
            functional_module_arrange_position = %s,
            module_arrange_position = %s,
            pathway_arrange_position = %s,
            molecule_arrange_position = %s,
            functional_module_position_limits = %s,
            module_position_limits = %s,
            pathway_position_limits = %s,
            molecule_position_limits = %s)
            ',
              input$relationship_network_include_functional_modules,
              input$relationship_network_include_modules,
              input$relationship_network_include_pathways,
              input$relationship_network_include_molecules,
              input$relationship_network_functional_module_text,
              input$relationship_network_llm_text,
              input$relationship_network_module_text,
              input$relationship_network_pathway_text,
              input$relationship_network_molecule_text,
              input$relationship_network_circular_plot,
              paste0(
                '"',
                input$relationship_network_functional_module_color,
                '"'
              ),
              paste0('"', input$relationship_network_module_color, '"'),
              paste0('"', input$relationship_network_pathway_color, '"'),
              paste0('"', input$relationship_network_molecule_color, '"'),
              input$relationship_network_functional_module_arrange_position,
              input$relationship_network_module_arrange_position,
              input$relationship_network_pathway_arrange_position,
              input$relationship_network_molecule_arrange_position,
              functional_module_position_limits,
              module_position_limits,
              pathway_position_limits,
              molecule_position_limits
            )

          relationship_network_code <-
            paste0(relationship_network_code1,
                   relationship_network_code2,
                   sep = "\n")

          relationship_network_code(relationship_network_code)

        }
      })

      output$relationship_network <-
        renderPlot({
          req(relationship_network())
          relationship_network()
        },
        res = 96)

      ######code for relationship_network
      ####show code
      observeEvent(input$show_relationship_network_code, {
        if (is.null(relationship_network_code()) ||
            length(relationship_network_code()) == 0) {
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
            relationship_network_code()
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


      output$download_relationship_network <-
        downloadHandler(
          filename = function() {
            paste0(
              "relationship_network_",
              input$relationship_network_level,
              ".",
              input$relationship_network_type
            )
          },
          content = function(file) {
            ggsave(
              file,
              plot = relationship_network(),
              width = input$relationship_network_width,
              height = input$relationship_network_height
            )
          }
        )

      observe({
        if (is.null(relationship_network()) ||
            length(relationship_network()) == 0) {
          shinyjs::disable("download_relationship_network")
        } else {
          shinyjs::enable("download_relationship_network")
        }
      })


      ## Go to results tab
      ####if there is not enriched_functional_module, show a warning message
      observeEvent(input$go2results_1, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched_functional_module available",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("results")
        }
      })

      observeEvent(input$go2results_2, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched_functional_module available",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("results")
        }
      })

      observeEvent(input$go2results_3, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched_functional_module available",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("results")
        }
      })

      observeEvent(input$go2results_4, {
        # Check if enriched_functional_module is available
        if (is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched_functional_module available",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          tab_switch("results")
        }
      })
    }
  )
}
