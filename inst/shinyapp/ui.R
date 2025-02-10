options(shiny.maxRequestSize = 300 * 1024 ^ 2)
options(shiny.legacy.datatable = TRUE)

setwd(r4projects::get_project_wd())
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
  "Merge Pathways", "merge_pathways", "cogs",
  "Merge Modules", "merge_modules", "cogs",
  # "Translation", "translation", "globe",
  "Data Visualization", "data_visualization", "chart-line",
  "LLM Interpretation", "llm_interpretation", "brain",
  "Results and Report", "results", "clipboard-list"
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
      #### 1. Introduction tab ====
      tabItem(tabName = "introduction",
              fluidPage(
                titlePanel("Introduction of MAPA"),
                fluidRow(
                  column(12,
                         htmltools::HTML(intro_cleaned_content)
                        )
                )
              )),
     #### 2. Tutorial tab ====
      tabItem(tabName = "tutorial",
              fluidPage(
                titlePanel("Tutorials of MAPA"),
                fluidRow(
                  column(12,
                         htmltools::HTML(tutorial_cleaned_content)
                         )
                )
              )),
     #### 3. Upload data tab ====
      upload_data_ui("upload_data_tab"),

     #### 4. Enrich pathways tab ====
      enrich_pathway_ui("enrich_pathway_tab"),

     #### 5. Merge pathways tab ====
      merge_pathways_ui("merge_pathways_tab"),

     #### 6. Merge modules tab ====
      merge_modules_ui("merge_modules_tab"),

     #### 7. Translation tab ====

     #### 8. Data visualization tab ====
      data_visualization_ui("data_visualization_tab"),

     #### 9. LLM Interpretation tab ====
      llm_interpretation_ui("llm_interpretation_tab"),

     #### 10. Result and report tab =====
      results_ui("results_tab")
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

