options(shiny.maxRequestSize = 300 * 1024 ^ 2)
options(shiny.legacy.datatable = TRUE)

setwd(r4projects::get_project_wd())
# source("inst/shinyapp/R/0_utils.R")
# source("inst/shinyapp/R/3_upload_data.R")
# source("inst/shinyapp/R/4_enrich_pathway.R")
devtools::load_all()

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

if (!require(future)) {
  install.packages("future")
  library(future)
}

if (!require(promises)) {
  install.packages("promises")
  library(promises)
}

if (!require(htmltools)) {
  install.packages("htmltools")
  library(htmltools)
}

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
    sidebarMenu(
      id = "tabs",
      menuItem(text = "Introduction", tabName = "introduction", icon = icon("info-circle")),
      menuItem(text = "Tutorial", tabName = "tutorial", icon = icon("book")),
      menuItem(text = "Upload Data", tabName = "upload_data", icon = icon("upload")),
      menuItem(text = "Enrich Pathways", tabName = "enrich_pathways", icon = icon("cogs")),
      menuItem(text = "Pathway Clustering", tabName = NULL, icon = icon("sitemap"),
               menuItem(text = HTML("Method1:<br>Overlap / semantic → Modules"), tabName = NULL,
                        menuSubItem(text = "Step1: Merge Pathways", tabName = "merge_pathways", icon = NULL),
                        menuSubItem(text = "Step2: Merge Modules", tabName = "merge_modules", icon = NULL)
                        ),
               menuItem(text = HTML("Method2:<br>Embed → Modules"), tabName = "embed_cluster_pathways")
      ),
      menuItem(text = "LLM Interpretation", tabName = "llm_interpretation", icon = icon("brain")),
      menuItem(text = "Data Visualization", tabName = "data_visualization", icon = icon("chart-line")),
      menuItem(text = "Results & Report", tabName = "results", icon = icon("clipboard-list"))
    )
  ),

  ## dashboard body code ====
  dashboardBody(
    shinyjs::useShinyjs(),

    tags$script(HTML("
      // Custom function to toggle sidebar
      function toggleSidebar() {
        var sidebarElement = $('.main-sidebar');
        var bodyElement = $('.content-wrapper');

        if (sidebarElement.css('display') === 'none') {
          sidebarElement.show();
          bodyElement.css('margin-left', '230px');
        } else {
          sidebarElement.hide();
          bodyElement.css('margin-left', '0px');
        }
      }

      // Override the default sidebar toggle behavior
      $(document).ready(function() {
        $('.sidebar-toggle').on('click', function(e) {
          e.preventDefault();
          e.stopPropagation();
          toggleSidebar();
        });

        // Fix for submenu toggle functionality
        $('.treeview > a').on('click', function(e) {
          e.preventDefault();
          $(this).parent().toggleClass('active');
          $(this).parent().children('.treeview-menu').slideToggle('fast');
        });
      });
    ")),

    # Add CSS to make sidebar wider
    tags$head(
      tags$style(HTML("
      .main-sidebar {
        width: 230px !important;
      }
    "))
    ),

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

     #### 5-6. Pathway clustering tab ===
     #### 5a. Merge pathways tab ====
     merge_pathways_ui("merge_pathways_tab"),

     #### 6a. Merge modules tab ====
     merge_modules_ui("merge_modules_tab"),

     #### 5-6b. Embed and cluster pathways tab =====
     embed_cluster_pathways_ui("embed_cluster_pathways_tab"),

     #### 7. Translation tab ====

     #### 8. LLM Interpretation tab ====
     llm_interpretation_ui("llm_interpretation_tab"),

     #### 9. Data visualization tab ====
     data_visualization_ui("data_visualization_tab"),

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

