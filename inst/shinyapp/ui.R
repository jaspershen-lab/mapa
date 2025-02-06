options(shiny.maxRequestSize = 300 * 1024 ^ 2)
options(shiny.legacy.datatable = TRUE)

# setwd(r4projects::get_project_wd())
# source("inst/shinyapp/R/0_utils.R")
# source("inst/shinyapp/R/3_upload_data.R")
# source("inst/shinyapp/R/4_enrich_pathway.R")
devtools::load_all()

if (!require(tidyverse)) {
  install.packages("tidyverse")
  library(tidyverse)
}

if (!require(shiny)) {
  install.packages("shiny")
  library(shiny)
}

if (!require(shinydashboard)) {
  install.packages("shinydashboard")
  library(shinydashboard)
}
if (!require(shinyjs)) {
  install.packages("shinyjs")
  library(shinyjs)
}
if (!require(shinyBS)) {
  install.packages("shinyBS")
  library(shinyBS)
}
if (!require(patchwork)) {
  install.packages("patchwork")
  library(patchwork)
}
if (!require(shinyWidgets)) {
  install.packages("shinyWidgets")
  library(shinyWidgets)
}
if (!require(markdown)) {
  install.packages("markdown")
  library(markdown)
}
if (!require(massdataset)) {
  remotes::install_github("tidymass/massdataset")
}
if (!require(mapa)) {
  remotes::install_github("jaspershen/mapa")
  library(mapa)
}

if (!require(readxl)) {
  install.packages("readxl")
  library(readxl)
}

if (!require(extrafont)) {
  install.packages("extrafont")
  library(extrafont)
  extrafont::loadfonts()
}

# Define items ====
menu_var <- tibble::tribble(
  ~ text, ~ tabName, ~ icon,
  "Introduction", "introduction", "info-circle",
  "Totorial", "tutorial", "book",
  "Upload Data", "upload_data", "upload",
  "Enrich Pathways", "enrich_pathways", "cogs",
  "Merge Pathways", "merge_pathways", "cogs"#,
  # "Translation", "translation", "globe",
  # "Data Visualization", "data_visualization", "chart-line",
  # "LLM Interpretation", "llm_interpretation", "brain",
  # "Results and Report", "results", "clipboard-list"
)

menu_items <- purrr::pmap(
  menu_var,
  function(text, tabName, icon) {
    menuItem(
      text = text,
      tabName = tabName,
      icon = icon(icon)
    )
  }
)

intro_html_content <- readLines("inst/shinyapp/files/introduction.html")
intro_cleaned_content <- grep("<(/?(html|head|body))>", intro_html_content, invert = TRUE, value = TRUE)

tutorial_html_content <- readLines("inst/shinyapp/files/tutorials.html")
tutorial_cleaned_content <- grep("<(/?(html|head|body))>", tutorial_html_content, invert = TRUE, value = TRUE)


# Define UI ====
ui <- dashboardPage(
  skin = "red",

  dashboardHeader(title = "MAPA"),

  ## sidebar of the app ====
  dashboardSidebar(
    do.call(sidebarMenu, c(
      list(id = "tabs"),
      menu_items)
    )
  ),

  ## dashboard body code ====
  dashboardBody(
    shinyjs::useShinyjs(),
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
  ### tabitems ====
    tabItems(
      #### 1. introduction tab ====
      tabItem(tabName = "introduction",
              fluidPage(
                titlePanel("Introduction of MAPA"),
                fluidRow(
                  column(12,
                         htmltools::HTML(intro_cleaned_content)
                        )
                )
              )),
     #### 2. tutorial tab ====
      tabItem(tabName = "tutorial",
              fluidPage(
                titlePanel("Tutorials of MAPA"),
                fluidRow(
                  column(12,
                         htmltools::HTML(tutorial_cleaned_content)
                         )
                )
              )),
     #### 3. upload data tab ====
      upload_data_ui("upload_data_tab"),

     #### 4. enrich pathways tab ====
      enrich_pathway_ui("enrich_pathway_tab"),

     #### 5. merge pathways tab ====
      merge_pathways_ui("merge_pathways_tab")

  #     #### Merge modules tab ====
  #     tabItem(
  #       tabName = "merge_modules",
  #       fluidPage(titlePanel("Merge Modules"),
  #                 fluidPage(
  #                   fluidRow(
  #                     column(4,
  #                            fluidRow(
  #                              column(6,
  #                                     selectInput(
  #                                       "measure.method.module",
  #                                       "Similarity method",
  #                                       choices = c("jaccard"),
  #                                       selected = "jaccard")
  #                                     ),
  #                              column(6,
  #                                     numericInput(
  #                                       "sim.cutoff.module",
  #                                       "Similarity cutoff",
  #                                       value = 0.5,
  #                                       min = 0,
  #                                       max = 1)
  #                                     )),
  #
  #                            actionButton(
  #                              "submit_merge_modules",
  #                              "Submit",
  #                              class = "btn-primary",
  #                              style = "background-color: #d83428; color: white;"),
  #
  #                            actionButton(
  #                              "go2translation",
  #                              "Next",
  #                              class = "btn-primary",
  #                              style = "background-color: #d83428; color: white;"),
  #
  #                            actionButton(
  #                              "show_merge_modules_code",
  #                              "Code",
  #                              class = "btn-primary",
  #                              style = "background-color: #d83428; color: white;"),
  #                            style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                            ),
  #                     column(8,
  #                            tabsetPanel(
  #                              tabPanel(
  #                                title = "Table",
  #                                DT::DTOutput("enriched_functional_modules"),
  #                                br(),
  #                                shinyjs::useShinyjs(),
  #                                downloadButton("download_enriched_functional_modules",
  #                                               "Download",
  #                                               class = "btn-primary",
  #                                               style = "background-color: #d83428; color: white;")
  #                                ),
  #                              tabPanel(
  #                                title = "Data visualization",
  #                                shiny::plotOutput("enirched_functional_module_plot"),
  #                                br(),
  #                                fluidRow(
  #                                  column(3,
  #                                         actionButton("generate_enirched_functional_module",
  #                                                      "Generate plot",
  #                                                      class = "btn-primary",
  #                                                      style = "background-color: #d83428; color: white;")
  #                                         ),
  #                                  column(3,
  #                                         checkboxInput("enirched_functional_module_plot_text", "Text", FALSE)
  #                                         ),
  #                                  column(3,
  #                                         checkboxInput("enirched_functional_module_plot_text_all", "Text all", FALSE)
  #                                         ),
  #                                column(3,
  #                                       numericInput(
  #                                         "enirched_functional_moduleplot__degree_cutoff",
  #                                         "Degree cutoff",
  #                                         value = 0,
  #                                         min = 0,
  #                                         max = 1000)
  #                                       )
  #                                )
  #                                ),
  #                              tabPanel(
  #                                title = "R object",
  #                                verbatimTextOutput("enriched_functional_module_object"),
  #                                br(),
  #                                shinyjs::useShinyjs(),
  #                                downloadButton("download_enriched_functional_module_object",
  #                                               "Download",
  #                                               class = "btn-primary",
  #                                               style = "background-color: #d83428; color: white;")
  #                                )
  #                          ))
  #                 )))
  #     ),
  #
  #
  #
  #     #### Translation tab====
  #     tabItem(
  #       tabName = "translation",
  #       fluidPage(titlePanel("Translation"),
  #                 fluidPage(fluidRow(
  #                   column(
  #                     4,
  #                     fluidRow(column(
  #                       5,
  #                       selectInput(
  #                         "translation_model",
  #                         "Model",
  #                         choices = c(
  #                           "Gemini" = "gemini",
  #                           "ChatGPT" = "chatgpt"),
  #                         selected = "gemini"
  #                       )
  #                     ),
  #                     column(
  #                       7,
  #                       textInput("translation_model_ai_key",
  #                                 "AI key",
  #                                 value = "")
  #                     )),
  #                     fluidRow(column(
  #                       12,
  #                       selectInput(
  #                         "translation_to",
  #                         "Translation to",
  #                         choices = c("Chinese" = "chinese",
  #                                     "Spanish" = "spanish",
  #                                     "English" = "english",
  #                                     "French" = "french",
  #                                     "German" = "german",
  #                                     "Italian" = "italian",
  #                                     "Japanese" = "japanese",
  #                                     "Korean" = "korean",
  #                                     "Portuguese" = "portuguese",
  #                                     "Russian" = "russian",
  #                                     "Spanish" = "spanish"),
  #                         selected = "chinese"
  #                       )
  #                     )),
  #                     actionButton(
  #                       "submit_translation",
  #                       "Submit",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #
  #                     actionButton(
  #                       "skip_translation",
  #                       "Skip",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #
  #                     actionButton(
  #                       "go2data_visualization",
  #                       "Next",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #                     actionButton(
  #                       "show_translation_code",
  #                       "Code",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #                     style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                   ),
  #                   column(8,
  #                          tabsetPanel(
  #                            tabPanel(
  #                              title = "R object",
  #                              verbatimTextOutput("enriched_functional_module_object2"),
  #                              br(),
  #                              shinyjs::useShinyjs(),
  #                              downloadButton("download_enriched_functional_module_object2",
  #                                             "Download",
  #                                             class = "btn-primary",
  #                                             style = "background-color: #d83428; color: white;")
  #                            )
  #                          )
  #                   )
  #                 )))
  #     ),
  #
  #     #### Data visualization tab ====
  #     tabItem(tabName = "data_visualization",
  #             fluidPage(
  #               titlePanel("Data Visualization"),
  #               tabsetPanel(
  #                 tabPanel(
  #                   title = "Barplot",
  #                   fluidRow(
  #                     column(4,
  #                            br(),
  #                            fluidRow(
  #                              column(8,
  #                                     fileInput(inputId = "upload_enriched_functional_module",
  #                                               label = tags$span("Upload functional module",
  #                                                                 shinyBS::bsButton("upload_functional_module_info",
  #                                                                          label = "",
  #                                                                          icon = icon("info"),
  #                                                                          style = "info",
  #                                                                          size = "extra-small")),
  #                                               accept = ".rda"),
  #                                     bsPopover(
  #                                       id = "upload_functional_module_info",
  #                                       title = "",
  #                                       content = "You can upload the functional module file here for data visualization only.",
  #                                       placement = "right",
  #                                       trigger = "hover",
  #                                       options = list(container = "body")
  #                                     )
  #                                     )
  #                              ),
  #                            fluidRow(
  #                              column(4,
  #                                     selectInput(
  #                                       inputId = "barplot_level",
  #                                       label = "Level",
  #                                       choices = c(
  #                                         "FM" = "functional_module",
  #                                         "Module" = "module",
  #                                         "Pathway" = "pathway"),
  #                                       selected = "functional_module")
  #                                     ),
  #                              column(4,
  #                                     numericInput(
  #                                       inputId = "barplot_top_n",
  #                                       label = "Top N",
  #                                       value = 5,
  #                                       min = 1,
  #                                       max = 1000)
  #                                     ),
  #                              column(4,
  #                                     selectInput(
  #                                       "line_type",
  #                                       "Line type",
  #                                       choices = c(
  #                                         "Straight" = "straight",
  #                                         "Meteor" = "meteor"))
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(6,
  #                                     numericInput(
  #                                       "barplot_y_lable_width",
  #                                       "Y label width",
  #                                       value = 50,
  #                                       min = 20,
  #                                       max = 100)
  #                              ),
  #                              column(6,
  #                                     numericInput("barplot_count_cutoff",
  #                                                  "Count cutoff",
  #                                                  value = 5,
  #                                                  min = 1,
  #                                                  max = 1000)
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(6,
  #                                     numericInput("barplot_p_adjust_cutoff",
  #                                                  "P-adjust cutoff",
  #                                                  value = 0.05,
  #                                                  min = 0,
  #                                                  max = 0.5)),
  #                              column(6,
  #                                     selectInput(
  #                                       "x_axis_name",
  #                                       "X axis name",
  #                                       choices = NULL
  #                                     ))
  #                            ),
  #                            fluidRow(
  #                              column(8,
  #                                     checkboxGroupInput("barplot_database",
  #                                                        "Database",
  #                                                        choices = c(
  #                                                          "GO" = "go",
  #                                                          "KEGG" = "kegg",
  #                                                          "Reactome" = "reactome"
  #                                                        ),
  #                                                        selected = c("go", "kegg", "reactome"),
  #                                                        inline = TRUE)
  #                                     ),
  #                              column(4,
  #                                     checkboxInput("barplot_translation", "Translation", FALSE)
  #                              )
  #                            ),
  #                            h4("Database color"),
  #                            fluidRow(
  #                            column(4,
  #                                   shinyWidgets::colorPickr(
  #                                     inputId = "barplot_go_color",
  #                                     label = "GO",
  #                                     selected = "#1F77B4FF",
  #                                     theme = "monolith",
  #                                     width = "100%")
  #                            ),
  #                            column(4,
  #                                   shinyWidgets::colorPickr(
  #                                     inputId = "barplot_kegg_color",
  #                                     label = "KEGG",
  #                                     selected = "#FF7F0EFF",
  #                                     theme = "monolith",
  #                                     width = "100%")
  #                                   ),
  #                            column(4,
  #                                   shinyWidgets::colorPickr(
  #                                     inputId = "barplot_reactome_color",
  #                                     label = "Reactome",
  #                                     selected = "#2CA02CFF",
  #                                     theme = "monolith",
  #                                     width = "100%")
  #                            )
  #                            ),
  #                            fluidRow(
  #                              column(4,
  #                                     selectInput("barplot_type", "Type",
  #                                                 choices = c("pdf", "png", "jpeg"))
  #                              ),
  #                              column(4,
  #                                     numericInput("barplot_width", "Width",
  #                                                  value = 7, min = 4, max = 20)
  #
  #                              ),
  #                              column(4,
  #                                     numericInput("barplot_height", "Height",
  #                                                  value = 7, min = 4, max = 20)
  #
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(12,
  #                                     actionButton("generate_barplot",
  #                                                  "Generate plot",
  #                                                  class = "btn-primary",
  #                                                  style = "background-color: #d83428; color: white;"),
  #                                     downloadButton("download_barplot",
  #                                                    "Download",
  #                                                    class = "btn-primary",
  #                                                    style = "background-color: #d83428; color: white;")
  #                                     )
  #                            ),
  #                            br(),
  #                            fluidRow(
  #                              column(12,
  #                                     actionButton(
  #                                       "go2llm_interpretation_1",
  #                                       "Next",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;"),
  #                                     actionButton(
  #                                       "show_barplot_code",
  #                                       "Code",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;")
  #                              )
  #                            ),
  #                            style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                            ),
  #                     column(8,
  #                            br(),
  #                            shiny::plotOutput("barplot")
  #                            )
  #                   )
  #                 ),
  #
  #                 ##### Module Similarity Network =====
  #                 tabPanel(
  #                   title = "Module Similarity Network",
  #                   fluidRow(
  #                     column(4,
  #                            br(),
  #                            fluidRow(
  #                              column(6,
  #                                     selectInput(
  #                                       "module_similarity_network_database",
  #                                       "Database",
  #                                       choices = c("GO" = "go",
  #                                                   "KEGG" = "kegg",
  #                                                   "Reactome" = "reactome"),
  #                                       selected = "go")
  #                                     ),
  #                              column(6,
  #                                     numericInput(
  #                                       "module_similarity_network_degree_cutoff",
  #                                       "Degree cutoff",
  #                                       value = 0,
  #                                       min = 0,
  #                                       max = 1000)
  #                                     )
  #                            ),
  #                            fluidRow(
  #                              column(6,
  #                                     selectInput(
  #                                       "module_similarity_network_level",
  #                                       "Level",
  #                                       choices = c(
  #                                         "FM" = "functional_module",
  #                                         "Module" = "module"),
  #                                       selected = "functional_module")
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(4,
  #                                     checkboxInput("module_similarity_network_translation", "Translation", FALSE)
  #                              ),
  #                              column(4,
  #                                     checkboxInput("module_similarity_network_text", "Text", FALSE)
  #                              ),
  #                              column(4,
  #                                     checkboxInput("module_similarity_network_text_all", "Text all", FALSE)
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(4,
  #                                     selectInput("module_similarity_network_type", "Type",
  #                                                 choices = c("pdf", "png", "jpeg"))
  #                              ),
  #                              column(4,
  #                                     numericInput("module_similarity_network_width", "Width",
  #                                                  value = 7, min = 4, max = 20)),
  #                              column(4,
  #                                     numericInput("module_similarity_network_height", "Height",
  #                                                  value = 7, min = 4, max = 20))
  #                            ),
  #                            fluidRow(
  #                              column(12,
  #                                     shinyjs::useShinyjs(),
  #                                     actionButton("generate_module_similarity_network",
  #                                                  "Generate plot",
  #                                                  class = "btn-primary",
  #                                                  style = "background-color: #d83428; color: white;"),
  #                                     shinyjs::useShinyjs(),
  #                                     downloadButton("download_module_similarity_network",
  #                                                    "Download",
  #                                                    class = "btn-primary",
  #                                                    style = "background-color: #d83428; color: white;")
  #                              )
  #                            ),
  #                            br(),
  #                            fluidRow(
  #                              column(12,
  #                                     actionButton(
  #                                       "go2llm_interpretation_2",
  #                                       "Next",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;"
  #                                     ),
  #                                     actionButton(
  #                                       "show_module_similarity_network_code",
  #                                       "Code",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;")
  #                                     )
  #                              ),
  #                            style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                     ),
  #                     column(8,
  #                            br(),
  #                            shiny::plotOutput("module_similarity_network")
  #                     )
  #                   )
  #                 ),
  #                 tabPanel(
  #                   title = "Module information",
  #                   fluidRow(
  #                     column(4,
  #                            br(),
  #                            fluidRow(
  #                              column(6,
  #                                     selectInput(
  #                                       "module_information_level",
  #                                       "Level",
  #                                       choices = c(
  #                                         "FM" = "functional_module",
  #                                         "Module" = "module"),
  #                                       selected = "functional_module")
  #                                     ),
  #                              column(6,
  #                                     selectInput(
  #                                       "module_information_database",
  #                                       "Database",
  #                                       choices = c("GO" = "go",
  #                                                   "KEGG" = "kegg",
  #                                                   "Reactome" = "reactome"),
  #                                       selected = "go")
  #                                     )
  #                            ),
  #                            fluidRow(
  #                              column(7,
  #                                     selectInput(
  #                                       "module_information_module_id",
  #                                       "Module ID",
  #                                       choices = NULL)
  #                              ),
  #                              column(5,
  #                                     checkboxInput("module_information_translation",
  #                                                   "Translation", FALSE)
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(4,
  #                                     selectInput("module_information_type", "Type",
  #                                                 choices = c("pdf", "png", "jpeg"))
  #                              ),
  #                              column(4,
  #                                     numericInput("module_information_width", "Width",
  #                                                  value = 7, min = 4, max = 30)),
  #                              column(4,
  #                                     numericInput("module_information_height", "Height",
  #                                                  value = 21, min = 4, max = 20))
  #                            ),
  #                            fluidRow(
  #                              column(12,
  #                                     shinyjs::useShinyjs(),
  #                                     actionButton("generate_module_information",
  #                                                  "Generate plot",
  #                                                  class = "btn-primary",
  #                                                  style = "background-color: #d83428; color: white;"),
  #                                     downloadButton("download_module_information",
  #                                                    "Download",
  #                                                    class = "btn-primary",
  #                                                    style = "background-color: #d83428; color: white;")
  #                              )
  #                            ),
  #                            br(),
  #                            fluidRow(
  #                              column(12,
  #                                     shinyjs::useShinyjs(),
  #                                     actionButton(
  #                                       "go2llm_interpretation_3",
  #                                       "Next",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;"),
  #                                     actionButton(
  #                                       "show_module_information_code",
  #                                       "Code",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;")
  #                                     )
  #                            ),
  #                            style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                     ),
  #                     column(8,
  #                            br(),
  #                            shiny::plotOutput("module_information1"),
  #                            shiny::plotOutput("module_information2"),
  #                            shiny::plotOutput("module_information3")
  #                     )
  #                   )
  #                 ),
  #                 tabPanel(
  #                   title = "Relationship network",
  #                   fluidRow(
  #                     column(4,
  #                            br(),
  #                            fluidRow(
  #                              column(5,
  #                                     checkboxInput("relationship_network_circular_plot",
  #                                                   "Circular layout", FALSE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_filter",
  #                                                   "Filter", FALSE)
  #                                     ),
  #                              column(4,
  #                                     checkboxInput("relationship_network_translation",
  #                                                   "Translation", FALSE)
  #                              )
  #                            ),
  #                            fluidRow(
  #                              column(6,
  #                                     selectInput(
  #                                       "relationship_network_level",
  #                                       "Filter Level",
  #                                       choices = c(
  #                                         "FM" = "functional_module",
  #                                         "Module" = "module"),
  #                                       selected = "functional_module")
  #                                     ),
  #                              column(6,
  #                                     selectInput(
  #                                       "relationship_network_module_id",
  #                                       "Module ID",
  #                                       choices = NULL,
  #                                       multiple = TRUE)
  #                                     )
  #                              ),
  #                            h4("Includes"),
  #                            fluidRow(
  #                              column(3,
  #                                     checkboxInput("relationship_network_include_functional_modules",
  #                                                   label = tags$span("FM",
  #                                                                     shinyBS::bsButton("functional_module_info",
  #                                                                                       label = "",
  #                                                                                       icon = icon("info"),
  #                                                                                       style = "info",
  #                                                                                       size = "extra-small")),
  #                                                   TRUE)),
  #                              bsPopover(
  #                                id = "functional_module_info",
  #                                title = "",
  #                                content = "FM is functional module",
  #                                placement = "right",
  #                                trigger = "hover",
  #                                options = list(container = "body")
  #                              ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_include_modules",
  #                                                   "Modules", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_include_pathways",
  #                                                   "Pathways", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_include_molecules",
  #                                                   "Molecules", TRUE)
  #                                     )
  #                            ),
  #                            h4("Colors"),
  #                            fluidRow(
  #                              column(3,
  #                                     shinyWidgets::colorPickr(
  #                                       inputId = "relationship_network_functional_module_color",
  #                                       label = "FM",
  #                                       selected = "#F05C3BFF",
  #                                       theme = "monolith",
  #                                       width = "100%"
  #                                     )
  #                                     ),
  #                              column(3,
  #                                     shinyWidgets::colorPickr(
  #                                       inputId = "relationship_network_module_color",
  #                                       label = "Module",
  #                                       selected = "#46732EFF",
  #                                       theme = "monolith",
  #                                       width = "100%"
  #                                     )
  #                                     ),
  #                              column(3,
  #                                     shinyWidgets::colorPickr(
  #                                       inputId = "relationship_network_pathway_color",
  #                                       label = "Pathway",
  #                                       selected = "#197EC0FF",
  #                                       theme = "monolith",
  #                                       width = "100%"
  #                                     )
  #                                     ),
  #                              column(3,
  #                                     shinyWidgets::colorPickr(
  #                                       inputId = "relationship_network_molecule_color",
  #                                       label = "Molecule",
  #                                       selected = "#3B4992FF",
  #                                       theme = "monolith",
  #                                       width = "100%"
  #                                     )
  #                                     )
  #                              ),
  #                            h4("Text"),
  #                            fluidRow(
  #                              column(3,
  #                                     checkboxInput("relationship_network_functional_module_text",
  #                                                   "FM", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_module_text",
  #                                                   "Module", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_pathway_text",
  #                                                   "Pathway", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_molecule_text",
  #                                                   "Molecules", FALSE)
  #                                     )
  #                            ),
  #                            h4("Text size"),
  #                            fluidRow(
  #                              column(3,
  #                                     numericInput("relationship_network_functional_module_text_size",
  #                                                  "FM",
  #                                                  value = 3, min = 0.3, max = 10)
  #                                     ),
  #                              column(3,
  #                                     numericInput("relationship_network_module_text_size",
  #                                                  "Module",
  #                                                  value = 3, min = 0.3, max = 10)
  #                                     ),
  #                              column(3,
  #                                     numericInput("relationship_network_pathway_text_size",
  #                                                  "Pathway",
  #                                                  value = 3, min = 0.3, max = 10)
  #                                     ),
  #                              column(3,
  #                                     numericInput("relationship_network_molecule_text_size",
  #                                                  "Molecule",
  #                                                  value = 3, min = 0.3, max = 10)
  #                                     )
  #                              ),
  #                            h4("Arrange posision"),
  #                            fluidRow(
  #                              column(3,
  #                                     checkboxInput("relationship_network_functional_module_arrange_position",
  #                                                   "FM", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_module_arrange_position",
  #                                                   "Module", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_pathway_arrange_position",
  #                                                   "Pathway", TRUE)
  #                                     ),
  #                              column(3,
  #                                     checkboxInput("relationship_network_molecule_arrange_position",
  #                                                   "Molecules", FALSE)
  #                                     )
  #                            ),
  #                            h4("Posision limits"),
  #                            fluidRow(
  #                              column(6,
  #                                     sliderInput(
  #                                       "relationship_network_functional_module_position_limits",
  #                                       "Functional module",
  #                                       min = 0, max = 1,
  #                                       value = c(0, 1))
  #                                     ),
  #                              column(6,
  #                                     sliderInput(
  #                                       "relationship_network_module_position_limits",
  #                                       "Module",
  #                                       min = 0, max = 1,
  #                                       value = c(0, 1))
  #                                     )
  #                            ),
  #                            fluidRow(
  #                              column(6,
  #                                     sliderInput(
  #                                       "relationship_network_pathway_position_limits",
  #                                       "Pathway",
  #                                       min = 0, max = 1,
  #                                       value = c(0, 1))
  #                                     ),
  #                              column(6,
  #                                     sliderInput(
  #                                       "relationship_network_molecule_position_limits",
  #                                       "Molecule",
  #                                       min = 0, max = 1,
  #                                       value = c(0, 1))
  #                                     )
  #                            ),
  #                            fluidRow(
  #                              column(4,
  #                                     selectInput("relationship_network_type",
  #                                                 "Type",
  #                                                 choices = c("pdf", "png", "jpeg")
  #                                                 )
  #                              ),
  #                              column(4,
  #                                     numericInput("relationship_network_width",
  #                                                  "Width",
  #                                                  value = 21, min = 4, max = 30)
  #                                     ),
  #                              column(4,
  #                                     numericInput("relationship_network_height",
  #                                                  "Height",
  #                                                  value = 7, min = 4, max = 20)
  #                                     )
  #                            ),
  #                            fluidRow(
  #                              column(12,
  #                                     actionButton("generate_relationship_network",
  #                                                  "Generate plot",
  #                                                  class = "btn-primary",
  #                                                  style = "background-color: #d83428; color: white;"),
  #                                     shinyjs::useShinyjs(),
  #                                     downloadButton("download_relationship_network",
  #                                                    "Download",
  #                                                    class = "btn-primary",
  #                                                    style = "background-color: #d83428; color: white;")
  #                              )
  #                            ),
  #                            br(),
  #                            fluidRow(
  #                              column(12,
  #                                     actionButton(
  #                                       "go2llm_interpretation_4",
  #                                       "Next",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;"),
  #                                     actionButton(
  #                                       "show_relationship_network_code",
  #                                       "Code",
  #                                       class = "btn-primary",
  #                                       style = "background-color: #d83428; color: white;"
  #                                     )
  #                                     )
  #                            ),
  #                            style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                     ),
  #                     column(8,
  #                            br(),
  #                            shiny::plotOutput("relationship_network")
  #                     )
  #                   )
  #                 )
  #               )
  #             )),
  #
  #     #### LLM Interpretation tab
  #     tabItem(
  #       tabName = "llm_interpretation",
  #       fluidPage(titlePanel("LLM Interpretation"),
  #                 fluidPage(fluidRow(
  #                   column(
  #                     4,
  #                     fluidRow(column(
  #                       5,
  #                       selectInput(
  #                         "llm_model",
  #                         "LLM model",
  #                         choices = c("ChatGPT" = "chatgpt"),
  #                         selected = "chatgpt"
  #                       )
  #                     ),
  #                     column(
  #                       7,
  #                       textInput("openai_key",
  #                                 "OpenAI key",
  #                                 value = "")
  #                     )),
  #                     fluidRow(column(
  #                       12,
  #                       textInput(
  #                         "llm_interpretation_disease",
  #                         "Disease or phenotype",
  #                         value = "pregnancy",
  #                         width = "100%"
  #                       )
  #                     )),
  #                     fluidRow(
  #                       column(4,
  #                              numericInput("llm_interpretation_top_n",
  #                                           "Top N",
  #                                           value = 5,
  #                                           min = 1,
  #                                           max = 1000)
  #                       ),
  #                       column(4,
  #                              numericInput("llm_interpretation_p_adjust_cutoff",
  #                                           "P-adjust cutoff",
  #                                           value = 0.05,
  #                                           min = 0,
  #                                           max = 0.5),
  #                       ),
  #                       column(4,
  #                              numericInput("llm_interpretation_count_cutoff",
  #                                           "Count cutoff",
  #                                           value = 5,
  #                                           min = 1,
  #                                           max = 1000)
  #                       )
  #                     ),
  #                     actionButton(
  #                       "submit_llm_interpretation",
  #                       "Submit",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #
  #                     actionButton(
  #                       "go2results",
  #                       "Next",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #                     actionButton(
  #                       "show_llm_interpretation_code",
  #                       "Code",
  #                       class = "btn-primary",
  #                       style = "background-color: #d83428; color: white;"
  #                     ),
  #                     style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                   ),
  #                   column(8,
  #                          tabsetPanel(
  #                            tabPanel(
  #                              title = "LLM interpretation results",
  #                              uiOutput("llm_interpretation_result"),
  #                              br(),
  #                              shinyjs::useShinyjs(),
  #                              downloadButton("download_llm_interpretation_result",
  #                                             "Download",
  #                                             class = "btn-primary",
  #                                             style = "background-color: #d83428; color: white;")
  #                            ),
  #                            tabPanel(
  #                              title = "Functional module table 1",
  #                              DT::DTOutput("llm_enriched_functional_modules1")
  #                            ),
  #                            tabPanel(
  #                              title = "Functional module table 2",
  #                              DT::DTOutput("llm_enriched_functional_modules2")
  #                            )
  #                          )
  #                          )
  #                 )))
  #       ),
  #
  #     #### result tab
  #     tabItem(
  #       tabName = "results",
  #       fluidPage(titlePanel("Results and Report"),
  #                 fluidPage(
  #                   column(4,
  #                          br(),
  #                          fluidRow(
  #                            actionButton(
  #                              inputId = "generate_report",
  #                              label = "Generate report",
  #                              class = "btn-primary",
  #                              style = "background-color: #d83428; color: white;"
  #                            ),
  #                            shinyjs::useShinyjs(),
  #                            downloadButton("download_report",
  #                                           "Download",
  #                                           class = "btn-primary",
  #                                           style = "background-color: #d83428; color: white;"),
  #                            actionButton(
  #                              inputId = "show_report_code",
  #                              label = "Code",
  #                              class = "btn-primary",
  #                              style = "background-color: #d83428; color: white;"
  #                            )
  #                          ),
  #                          style = "border-right: 1px solid #ddd; padding-right: 20px;"
  #                          ),
  #                   column(8,
  #                          tabsetPanel(
  #                            tabPanel(
  #                              title = "Report",
  #                              uiOutput("mapa_report")
  #                              ))
  #                          )
  #                 )
  #                 )
  #     )
    ),

  ### footer ====
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

