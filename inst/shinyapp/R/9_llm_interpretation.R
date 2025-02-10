llm_interpretation_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "llm_interpretation",
    fluidPage(titlePanel("LLM Interpretation"),
              fluidPage(fluidRow(
                column(
                  4,
                  fluidRow(column(
                    5,
                    selectInput(
                      ns("llm_model"),
                      "LLM model",
                      choices = c("ChatGPT" = "chatgpt"),
                      selected = "chatgpt"
                    )
                  ),
                  column(
                    7,
                    textInput(ns("openai_key"),
                              "OpenAI key",
                              value = "")
                  )),
                  fluidRow(column(
                    12,
                    textInput(
                      ns("llm_interpretation_disease"),
                      "Disease or phenotype",
                      value = "pregnancy",
                      width = "100%"
                    )
                  )),
                  fluidRow(
                    column(4,
                           numericInput(ns("llm_interpretation_top_n"),
                                        "Top N",
                                        value = 5,
                                        min = 1,
                                        max = 1000)
                    ),
                    column(4,
                           numericInput(ns("llm_interpretation_p_adjust_cutoff"),
                                        "P-adjust cutoff",
                                        value = 0.05,
                                        min = 0,
                                        max = 0.5),
                    ),
                    column(4,
                           numericInput(ns("llm_interpretation_count_cutoff"),
                                        "Count cutoff",
                                        value = 5,
                                        min = 1,
                                        max = 1000)
                    )
                  ),
                  actionButton(
                    ns("submit_llm_interpretation"),
                    "Submit",
                    class = "btn-primary",
                    style = "background-color: #d83428; color: white;"
                  ),

                  actionButton(
                    ns("go2results"),
                    "Next",
                    class = "btn-primary",
                    style = "background-color: #d83428; color: white;"
                  ),
                  actionButton(
                    ns("show_llm_interpretation_code"),
                    "Code",
                    class = "btn-primary",
                    style = "background-color: #d83428; color: white;"
                  ),
                  style = "border-right: 1px solid #ddd; padding-right: 20px;"
                ),
                column(8,
                       tabsetPanel(
                         tabPanel(
                           title = "LLM interpretation results",
                           uiOutput(ns("llm_interpretation_result")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_llm_interpretation_result"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "Functional module table 1",
                           shiny::dataTableOutput(ns("llm_enriched_functional_modules1"))
                         ),
                         tabPanel(
                           title = "Functional module table 2",
                           shiny::dataTableOutput(ns("llm_enriched_functional_modules2"))
                         )
                       )
                )
              )))
  )
}

llm_interpretation_server <- function(id, enriched_functional_module = NULL, tab_switch) {
  moduleServer(
    id,
    function(input, output, session) {
      # Define enriched_functional_module as a reactive value
      llm_interpretation_result <- reactiveVal("")

      llm_interpretation_code <- reactiveVal()

      openai_key <- reactiveVal()

      observeEvent(input$submit_llm_interpretation, {
        openai_key1 <-
          Sys.getenv("chatgpt_api_key")

        openai_key2 <-
          input$openai_key

        # Check if enriched_modules is available
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
          if (openai_key1 != "") {
            openai_key(openai_key1)
          } else{
            if (openai_key2 != "") {
              openai_key(openai_key2)
            } else{
              openai_key("")
            }
          }

          if (openai_key() == "") {
            showModal(
              modalDialog(
                title = "Warning",
                "No OpenAI Key provided. No interpretation will be generated.",
                easyClose = TRUE,
                footer = modalButton("Close")
              )
            )
          } else{
            set_chatgpt_api_key(api_key = openai_key())
          }

          # shinyjs::show("loading")

          withProgress(message = 'Analysis in progress...', {
            tryCatch({
              llm_interpretation_result <-
                interpret_pathways(
                  object = enriched_functional_module(),
                  p.adjust.cutoff = input$llm_interpretation_p_adjust_cutoff,
                  disease = input$llm_interpretation_disease,
                  count.cutoff = input$llm_interpretation_count_cutoff,
                  top_n = input$llm_interpretation_top_n
                )
              llm_interpretation_result(llm_interpretation_result)
            },
            error = function(e) {
              showModal(modalDialog(
                title = "Error",
                paste("Details:", e$message),
                easyClose = TRUE,
                footer = modalButton("Close")
              ))
              llm_interpretation_result("No result")
            })
          })
          # shinyjs::hide("loading")

          ##save code
          llm_interpretation_code <-
            sprintf(
              '
            llm_interpretation_result <-
            interpret_pathways(
            object = enriched_functional_module,
            p.adjust.cutoff = %s,
            disease = %s,
            count.cutoff = %s,
            top_n = %s)
            ',
              input$llm_interpretation_p_adjust_cutoff,
              paste0('"', input$llm_interpretation_disease, '"'),
              input$llm_interpretation_count_cutoff,
              input$llm_interpretation_top_n
            )

          llm_interpretation_code(llm_interpretation_code)
        }
      })

      output$llm_interpretation_result <-
        renderUI({
          shiny::HTML(markdown::markdownToHTML(llm_interpretation_result(),
                                               fragment.only = TRUE))
        })

      output$download_llm_interpretation_result <-
        shiny::downloadHandler(
          filename = function() {
            "llm_interpretation_result.md"
          },
          content = function(file) {
            writeLines(llm_interpretation_result(), file)
          }
        )

      observe({
        if (is.null(llm_interpretation_result()) ||
            length(llm_interpretation_result()) == 0) {
          shinyjs::disable("download_llm_interpretation_result")
        } else {
          shinyjs::enable("download_llm_interpretation_result")
        }
      })


      output$llm_enriched_functional_modules1 <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_functional_module()@merged_module$functional_module_result,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      output$llm_enriched_functional_modules2 <-
        shiny::renderDataTable({
          req(tryCatch(
            enriched_functional_module()@merged_module$result_with_module,
            error = function(e)
              NULL
          ))
        },
        options = list(pageLength = 10,
                       scrollX = TRUE))

      observe({
        if (is.null(llm_interpretation_result()) ||
            length(llm_interpretation_result()) == 0 ||
            llm_interpretation_result() == "") {
          shinyjs::disable("download_llm_interpretation_result")
        } else {
          shinyjs::enable("download_llm_interpretation_result")

        }
      })


      ####show code
      observeEvent(input$show_llm_interpretation_code, {
        if (is.null(llm_interpretation_code()) ||
            length(llm_interpretation_code()) == 0) {
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
            llm_interpretation_code()
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

      ###Go to result tab
      ####if there is not enriched_functional_module, show a warning message
      observeEvent(input$go2results, {
        # Check if enriched_functional_module is available
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
          tab_switch("results")
        }
      })

      return(llm_interpretation_result)
    }
  )
}
