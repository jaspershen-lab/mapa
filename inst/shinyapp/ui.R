library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(shinyWidgets)

ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "MAPA"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem(
        "Introduction",
        tabName = "introduction",
        icon = icon("info-circle")
      ),
      menuItem("Totorial", tabName = "tutorial", icon = icon("book")),
      menuItem("Upload data", tabName = "upload_data", icon = icon("upload")),
      menuItem(
        "Enrich pathways",
        tabName = "enrich_pathways",
        icon = icon("cogs")
      ),
      menuItem("Merge pathways", tabName = "merge_pathways", icon = icon("cogs")),
      menuItem("Merge modules", tabName = "merge_modules", icon = icon("cogs")),
      menuItem(
        "Data visualization",
        tabName = "data_visualization",
        icon = icon("chart-line")
      ),
      menuItem(
        "Results and report",
        tabName = "results",
        icon = icon("clipboard-list")
      )
    )
  ),
  dashboardBody(
    useShinyjs(),
    div(
      id = "loading",
      hidden = TRUE,
      class = "loading-style",
      "",
      tags$img(src = "loading.gif",
               height = "200px")
    ),

    tags$style(
      HTML(
        "
      .content-wrapper {
        padding-bottom: 120px;
      }
      .loading-style {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        text-align: center;
        z-index: 100;
      }
    "
      )
    ),

    tabItems(
      tabItem(tabName = "introduction",
              fluidPage(
                titlePanel("Introduction of MAPA"),
                includeHTML("files/introduction.html")
              )),
      tabItem(tabName = "tutorial",
              fluidPage(
                titlePanel("Tutorials of MAPA"),
                includeHTML("files/tutorials.html")
              )),
      tabItem(tabName = "upload_data",
              fluidPage(
                titlePanel("Upload data"),
                fluidRow(
                  column(
                    4,
                    fileInput(
                      "variable_info",
                      "Choose File",
                      accept = c(
                        "text/csv",
                        "text/comma-separated-values,text/plain",
                        ".csv",
                        ".xlsx",
                        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        "application/vnd.ms-excel"
                      )
                    ),
                    checkboxInput("use_example", "Use example", FALSE),
                    radioButtons(
                      "id_type",
                      "ID type:",
                      choices = list(
                        "ENSEMBL" = "ensembl",
                        "UniProt" = "uniprot",
                        "EntrezID" = "entrezid"
                      )
                    ),

                    fluidRow(column(
                      3,
                      actionButton("map_id",
                                   "Map ID",
                                   class = "btn-primary",
                                   style = "background-color: #d83428; color: white;")
                    ),
                    column(
                      3,
                      actionButton(
                        "go2enrich_pathways",
                        "Next",
                        class = "btn-primary",
                        style = "background-color: #d83428; color: white;"
                      ),
                    )),
                    br(),
                    br(),
                    actionButton(
                      "show_code_upload_data",
                      "Show/Hide Code",
                      class = "btn-primary",
                      style = "background-color: #d83428; color: white;"
                    ),

                    style = "border-right: 1px solid #ddd; padding-right: 20px;"
                  ),
                  column(
                    8,
                    shiny::dataTableOutput("contents"),
                    downloadButton("download_variable_info",
                                   "Download")
                  )
                ),
                ###show code
                br(),
                uiOutput("data_upload_code")
              )),
      tabItem(tabName = "enrich_pathways",
              fluidPage(
                titlePanel("Enrich Pathways"),
                fluidPage(
                  fluidRow(
                    column(
                      4,
                      checkboxGroupInput(
                        "pathway_database",
                        "Database:",
                        choices = c(
                          "GO" = "go",
                          "KEGG" = "kegg",
                          "Reactome" = "reactome"
                        ),
                        selected = c("go", "kegg", "reactome")
                      ),
                      selectInput(
                        "organism",
                        "Organism:",
                        choices = list(
                          "Human" = "hsa",
                          "Rat" = "rno",
                          "Mouse" = "mmu",
                          "C. elegans" = "cel",
                          "Yeast" = "sce",
                          "Zebrafish" = "dre",
                          "Fruit fly" = "dme"
                        ),
                        # Note: "Fly" is duplicated here
                        selected = "hsa"
                      ),
                      numericInput(
                        "p_value_cutoff",
                        "P-value Cutoff:",
                        value = 0.05,
                        min = 0,
                        max = 0.5
                      ),
                      selectInput(
                        "p_adjust_method",
                        "P-Adjust Method:",
                        choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                        selected = "BH"
                      ),
                      sliderInput(
                        "gene_set_size",
                        "Gene Set Size:",
                        min = 5,
                        max = 2000,
                        value = c(10, 500)
                      ),

                      fluidRow(column(
                        3,
                        actionButton(
                          "submit_enrich_pathways",
                          "Submit",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      ),
                      column(
                        3,
                        actionButton(
                          "go2merge_pathways",
                          "Next",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      )),

                      br(),
                      br(),
                      actionButton(
                        "show_code_enrich_pathways",
                        "Show/Hide Code",
                        class = "btn-primary",
                        style = "background-color: #d83428; color: white;"
                      ),
                      style = "border-right: 1px solid #ddd; padding-right: 20px;"
                    ),
                    column(
                      8,
                      tabsetPanel(
                        id = "enriched_pathways_result",
                        tabPanel(
                          title = "GO",
                          shiny::dataTableOutput("enriched_pathways_go"),
                          downloadButton("download_enriched_pathways_go",
                                         "Download")
                        ),
                        tabPanel(
                          title = "KEGG",
                          shiny::dataTableOutput("enriched_pathways_kegg"),
                          downloadButton("download_enriched_pathways_kegg",
                                         "Download")
                        ),
                        tabPanel(
                          title = "Reactome",
                          shiny::dataTableOutput("enriched_pathways_reactome"),
                          downloadButton("download_enriched_pathways_reactome",
                                         "Download")
                        )
                      )
                    )
                  ),
                  ###show code
                  br(),
                  uiOutput("enrich_pathways_code")
                )
              )),
      tabItem(tabName = "merge_pathways",
              fluidPage(
                titlePanel("Merge pathways"),
                fluidPage(
                  fluidRow(
                    column(
                      4,
                      h3("GO"),
                      fluidRow(column(
                        6,
                        numericInput(
                          "p.adjust.cutoff.go",
                          "Adjust P Value Cutoff:",
                          value = 0.05,
                          min = 0,
                          max = 0.5
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "count.cutoff.go",
                          "Gene Count Cutoff:",
                          value = 5,
                          min = 0,
                          max = 1000
                        )
                      )),
                      h4("Network"),
                      fluidRow(column(
                        6,
                        selectInput(
                          "measure.method.go",
                          "Similarity Method:",
                          choices = c("Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS", "jaccard"),
                          selected = "Wang"
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "sim.cutoff.go",
                          "Similarity Cutoff:",
                          value = 0.5,
                          min = 0,
                          max = 1
                        )
                      )),

                      h3("KEGG"),
                      fluidRow(column(
                        6,
                        numericInput(
                          "p.adjust.cutoff.kegg",
                          "Adjust P Value Cutoff:",
                          value = 0.05,
                          min = 0,
                          max = 0.5
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "count.cutoff.kegg",
                          "Gene Count Cutoff:",
                          value = 5,
                          min = 0,
                          max = 1000
                        )
                      )),
                      h4("Network"),
                      fluidRow(column(
                        6,
                        selectInput(
                          "measure.method.kegg",
                          "Similarity Method:",
                          choices = c("jaccard"),
                          selected = "jaccard"
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "sim.cutoff.kegg",
                          "Similarity Cutoff:",
                          value = 0.5,
                          min = 0,
                          max = 1
                        )
                      )),

                      h3("Reactome"),
                      fluidRow(column(
                        6,
                        numericInput(
                          "p.adjust.cutoff.reactome",
                          "Adjust P Value Cutoff:",
                          value = 0.05,
                          min = 0,
                          max = 0.5
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "count.cutoff.reactome",
                          "Gene Count Cutoff:",
                          value = 5,
                          min = 0,
                          max = 1000
                        )
                      )),
                      h4("Network"),
                      fluidRow(column(
                        6,
                        selectInput(
                          "measure.method.reactome",
                          "Similarity Method:",
                          choices = c("jaccard"),
                          selected = "jaccard"
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "sim.cutoff.reactome",
                          "Similarity Cutoff:",
                          value = 0.5,
                          min = 0,
                          max = 1
                        )
                      )),

                      fluidRow(column(
                        3,
                        actionButton(
                          "submit_merge_pathways",
                          "Submit",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      ),
                      column(
                        3,
                        actionButton(
                          "go2merge_modules",
                          "Next",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      )),

                      br(),
                      br(),
                      actionButton(
                        "show_code_merge_pathways",
                        "Show/Hide Code",
                        class = "btn-primary",
                        style = "background-color: #d83428; color: white;"
                      ),
                      style = "border-right: 1px solid #ddd; padding-right: 20px;"
                    ),
                    column(8,
                           tabsetPanel(
                             tabPanel("Table",
                                      tabsetPanel(
                                        tabPanel(
                                          title = "GO",
                                          shiny::dataTableOutput("merged_pathway_go"),
                                          downloadButton("download_merged_pathway_go",
                                                         "Download")
                                        ),
                                        tabPanel(
                                          title = "KEGG",
                                          shiny::dataTableOutput("merged_pathway_kegg"),
                                          downloadButton("download_merged_pathway_kegg",
                                                         "Download")
                                        ),
                                        tabPanel(
                                          title = "Reactome",
                                          shiny::dataTableOutput("merged_pathway_reactome"),
                                          downloadButton("download_merged_pathway_reactome",
                                                         "Download")
                                        )
                                      )),
                             tabPanel(
                               "Data Visualization",
                               tabsetPanel(
                                 tabPanel(
                                   title = "GO",
                                   shiny::plotOutput("enirched_module_go_plot"),
                                   fluidRow(
                                     column(
                                       3,
                                       actionButton("generate_enirched_module_plot_go", "Generate Plot")
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_go", "Text", FALSE)
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_all_go", "Text all", FALSE)
                                     ),
                                     column(
                                       3,
                                       numericInput(
                                         "enirched_module_plot_degree_cutoff_go",
                                         "Degree cutoff:",
                                         value = 0,
                                         min = 0,
                                         max = 1000
                                       )
                                     )
                                   )
                                 ),
                                 tabPanel(
                                   title = "KEGG",
                                   shiny::plotOutput("enirched_module_kegg_plot"),
                                   fluidRow(
                                     column(
                                       3,
                                       actionButton("generate_enirched_module_plot_kegg", "Generate Plot")
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_kegg", "Text", FALSE)
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_all_kegg", "Text all", FALSE)
                                     ),
                                     column(
                                       3,
                                       numericInput(
                                         "enirched_module_plot_degree_cutoff_kegg",
                                         "Degree cutoff:",
                                         value = 0,
                                         min = 0,
                                         max = 1000
                                       )
                                     )
                                   )
                                 ),
                                 tabPanel(
                                   title = "Reactome",
                                   shiny::plotOutput("enirched_module_reactome_plot"),
                                   fluidRow(
                                     column(
                                       3,
                                       actionButton("generate_enirched_module_plot_reactome", "Generate Plot")
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_reactome", "Text", FALSE)
                                     ),
                                     column(
                                       3,
                                       checkboxInput("enirched_module_plot_text_all_reactome", "Text all", FALSE)
                                     ),
                                     column(
                                       3,
                                       numericInput(
                                         "enirched_module_plot_degree_cutoff_reactome",
                                         "Degree cutoff:",
                                         value = 0,
                                         min = 0,
                                         max = 1000
                                       )
                                     )
                                   )
                                 )
                               )
                             )
                           ))
                  ),
                  ###show code
                  br(),
                  uiOutput("merge_pathways_code")
                )
              )),

      tabItem(
        tabName = "merge_modules",
        fluidPage(titlePanel("Merge modules"),
                  fluidPage(fluidRow(
                    column(
                      4,
                      fluidRow(column(
                        6,
                        selectInput(
                          "measure.method.module",
                          "Similarity Method:",
                          choices = c("jaccard"),
                          selected = "jaccard"
                        )
                      ),
                      column(
                        6,
                        numericInput(
                          "sim.cutoff.module",
                          "Similarity Cutoff:",
                          value = 0.5,
                          min = 0,
                          max = 1
                        )
                      )),

                      fluidRow(column(
                        3,
                        actionButton(
                          "submit_merge_modules",
                          "Submit",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      ),
                      column(
                        3,
                        actionButton(
                          "go2data_visualization",
                          "Next",
                          class = "btn-primary",
                          style = "background-color: #d83428; color: white;"
                        )
                      )),

                      br(),
                      br(),
                      actionButton(
                        "show_code_merge_modules",
                        "Show/Hide Code",
                        class = "btn-primary",
                        style = "background-color: #d83428; color: white;"
                      ),
                      style = "border-right: 1px solid #ddd; padding-right: 20px;"
                    ),
                    column(8,
                           tabsetPanel(
                             tabPanel(
                               title = "Table",
                               shiny::dataTableOutput("enriched_functional_modules"),
                               downloadButton("download_enriched_functional_modules",
                                              "Download")
                             ),
                             tabPanel(
                               title = "Data Visualization",
                               shiny::plotOutput("enirched_functional_module_plot"),
                               fluidRow(
                                 column(
                                   3,
                                   actionButton("generate_enirched_functional_module", "Generate Plot")
                                 ),
                                 column(
                                   3,
                                   checkboxInput("enirched_functional_module_plot_text", "Text", FALSE)
                                 ),
                                 column(
                                   3,
                                   checkboxInput("enirched_functional_module_plot_text_all", "Text all", FALSE)
                                 ),
                                 column(
                                   3,
                                   numericInput(
                                     "enirched_functional_moduleplot__degree_cutoff",
                                     "Degree cutoff:",
                                     value = 0,
                                     min = 0,
                                     max = 1000
                                   )
                                 )
                               )
                             )
                           ))
                  ))),
        ###show code
        br(),
        uiOutput("merge_modules_code")
      ),
      tabItem(tabName = "data_visualization",
              fluidPage(
                titlePanel("Data Visualization"),
                tabsetPanel(
                  tabPanel(
                    title = "Barplot",
                    fluidRow(
                      column(4,
                             br(),
                             fluidRow(
                               column(12,
                                      fileInput(inputId = "upload_enriched_functional_module",
                                                label = tags$span("Upload functional module",
                                                                  shinyBS::bsButton("upload_functional_module_info",
                                                                           label = "",
                                                                           icon = icon("info"),
                                                                           style = "info",
                                                                           size = "extra-small")),
                                                accept = ".rda"),
                                      bsPopover(
                                        id = "upload_functional_module_info",
                                        title = "",
                                        content = "You can upload the functional module file here for data visualization only.",
                                        placement = "right",
                                        trigger = "hover",
                                        options = list(container = "body")
                                      )

                                      )
                               ),
                             fluidRow(
                               column(4,
                                      selectInput(
                                        "barplot_level",
                                        "Level",
                                        choices = c(
                                          "Pathway" = "pathway",
                                          "Module" = "module",
                                          "Functional module" = "functional_module")
                                      )
                                      ),
                               column(4,numericInput("barplot_top_n", "Top N:",
                                                          value = 5,
                                                          min = 1,
                                                          max = 1000)
                                      ),
                               column(4,
                                      selectInput(
                                        "line_type",
                                        "Line type",
                                        choices = c(
                                          "Straight" = "straight",
                                          "Meteor" = "meteor")
                                      )
                               )
                             ),
                             fluidRow(
                               column(4,
                                      numericInput(
                                        "barplot_y_lable_width",
                                        "Y Label Width",
                                        value = 50,
                                        min = 20,
                                        max = 100
                                      )
                               ),
                               column(4,
                                      numericInput("barplot_p_adjust_cutoff",
                                                   "P-Adjust Cutoff",
                                                   value = 0.05,
                                                   min = 0,
                                                   max = 0.5),
                               ),
                               column(4,
                                      numericInput("barplot_count_cutoff",
                                                   "Count cutoff:",
                                                   value = 5,
                                                   min = 1,
                                                   max = 1000)
                               )
                             ),
                             fluidRow(
                               column(12,
                                      checkboxGroupInput("barplot_database",
                                                         "Database",
                                                         choices = c(
                                                           "GO" = "go",
                                                           "KEGG" = "kegg",
                                                           "Reactome" = "reactome"
                                                         ),
                                                         selected = c("go", "kegg", "reactome"),
                                                         inline = TRUE
                                      )
                                      )
                             ),
                             h4("Database Color"),
                             fluidRow(
                             column(4,
                                    shinyWidgets::colorPickr(
                                      inputId = "barplot_go_color",
                                      label = "GO",
                                      selected = "#1F77B4FF",
                                      theme = "monolith",
                                      width = "100%"
                                    )
                             ),
                             column(4,
                                    shinyWidgets::colorPickr(
                                      inputId = "barplot_kegg_color",
                                      label = "KEGG",
                                      selected = "#FF7F0EFF",
                                      theme = "monolith",
                                      width = "100%"
                                    )
                                    ),
                             column(4,
                                    shinyWidgets::colorPickr(
                                      inputId = "barplot_reactome_color",
                                      label = "Reactome",
                                      selected = "#2CA02CFF",
                                      theme = "monolith",
                                      width = "100%"
                                    )
                             )
                             ),
                             fluidRow(
                               column(12,
                                      downloadButton("download_barplot", "Download"))
                             ),
                             fluidRow(
                               column(4,
                                      selectInput("barplot_type", "Type",
                                                  choices = c("pdf", "png", "jpeg"))
                                      ),
                               column(4,
                                      numericInput("barplot_width", "Width",
                                                   value = 7, min = 4, max = 20)),
                               column(4,
                                      numericInput("barplot_height", "Height",
                                                   value = 7, min = 4, max = 20))
                             ),
                             style = "border-right: 1px solid #ddd; padding-right: 20px;"
                             ),
                      column(8,
                             shiny::plotOutput("barplot"),
                             br(),
                             fluidRow(
                               column(3,
                                      actionButton("generate_barplot",
                                                   "Generate Plot",
                                                   class = "btn-primary",
                                                   style = "background-color: #d83428; color: white;")
                               ),
                               column(3,
                                      actionButton(
                                        "show_barplot_code",
                                        "Show/Hide Code",
                                        class = "btn-primary",
                                        style = "background-color: #d83428; color: white;"
                                      ))
                             ),
                             br(),
                             uiOutput("barplot_code")
                             )
                    )
                  )
                )
              ))
    ),
    tags$footer(
      div(
        style = "background-color: #ecf0f4; display: flex; align-items: center; justify-content: left; padding: 10px; height: 80px; position: fixed; bottom: 0; width: 100%; z-index: 100; border-top: 1px solid #ccc;",
        tags$img(
          src = "shen_lab_logo.png",
          height = "50px",
          style = "margin-right: 15px;"
        ),
        div(
          HTML("The Shen Lab at Nanyang Technological University Singapore"),
          HTML("<br>"),
          tags$a(
            href = "http://www.shen-lab.org",
            target = "_blank",
            tags$i(class = "fa fa-house", style = "color: #e04c3c;"),
            " Shen Lab",
            style = "text-align: left; margin-left: 10px;"
          ),
          tags$a(
            href = "https://www.shen-lab.org/#contact",
            target = "_blank",
            tags$i(class = "fa fa-envelope", style = "color: #e04c3c;"),
            " Email",
            style = "text-align: left; margin-left: 10px;"
          ),
          tags$a(
            href = "https://github.com/jaspershen/mapa",
            target = "_blank",
            tags$i(class = "fa fa-github", style = "color: #e04c3c;"),
            " GitHub",
            style = "text-align: left; margin-left: 10px;"
          ),
          style = "text-align: left;"
        ),
        tags$img(
          src = "mapa_logo.png",
          height = "50px",
          style = "margin-left: 15px;"
        )
      )
    )
  )
)
