llm_interpretation_ui <- function(id) {
  ns <- NS(id)
  tabItem(
    tabName = "llm_interpretation",
    fluidPage(titlePanel("LLM Interpretation"),
              fluidPage(
                fluidRow(
                ## Input layout ====
                column(
                  4,
                  fluidRow(
                    column(12,
                           fileInput(
                             inputId = ns("upload_interpretation_result"),
                             label   = tags$span(
                               "Upload interpretation result (.rda)",
                               shinyBS::bsButton(ns("upload_interpretation_result_info"),
                                                 label = "", icon = icon("info"),
                                                 style = "info", size = "extra-small")
                             ),
                             accept = ".rda"
                           ),
                           bsPopover(
                             id        = ns("upload_interpretation_result_info"),
                             title     = "",
                             content   = "Drop a previously saved result with LLM interpretation (.rda) here to inspect it without re-running the LLM.",
                             placement = "right", trigger = "hover",
                             options   = list(container = "body")
                           )
                    ),
                  ),
                  tags$hr(style = "border-top: 1px solid #ddd; margin-top: 10px; margin-bottom: 10px;"),
                  fluidRow(
                    column(12,
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
                             content = "You can upload the functional module result here.",
                             placement = "right",
                             trigger = "hover",
                             options = list(container = "body")
                           )
                    )
                  ),
                  # br(),
                  fluidRow(
                    column(
                      6,
                      textInput(ns("llm_model"),
                                "LLM model",
                                value = "gpt-4o-mini-2024-07-18")
                    ),
                    column(
                      6,
                      textInput(ns("embedding_model"),
                                "Embedding model",
                                value = "text-embedding-3-small")
                    )),
                  fluidRow(
                    column(
                      5,
                      numericInput(ns("years"),
                                   "Years to search",
                                   value = 5,
                                   min = 1,
                                   max = 1000)
                    ),
                    column(
                      7,
                      textInput(
                        ns("phenotype"),
                        "Disease or phenotype",
                        value = "NULL",
                        width = "100%"
                      )
                    )),
                  fluidRow(
                    column(
                      12,
                      textInput(ns("api_key"),
                                "API key",
                                value = "")
                    )
                  ),
                  tags$p("(Optional) Select a local folder with PDFs as local corpus:"),
                  fluidRow(
                    column(2,
                           shinyFiles::shinyDirButton(
                             ns("local_corpus_dir"),
                             "Browse",
                             title = "Please select the local corpus directory",
                             class = "btn-primary",
                             style = "background-color: #d83428; color: white; width: 75px;"
                           )
                    ),
                    column(10,
                           textInput(ns("corpus_path_display"), NULL, width = "100%", placeholder = "No local corpus directory selected", value = "")
                           )
                  ),
                  br(),
                  tags$p("Select a local folder to store generated embeddings:"),
                  fluidRow(
                    column(2,
                           shinyFiles::shinyDirButton(
                             ns("embedding_output_dir"),
                             "Browse",
                             title = "Please select a directory to save embeddings",
                             class = "btn-primary",
                             style = "background-color: #d83428; color: white; width: 75px;"
                           )
                           ),
                    column(10,
                           textInput(ns("embedding_path_display"), NULL, width = "100%", placeholder = "No embeddings output directory selected", value = "")
                           )
                  ),
                  br(),
                  # fluidRow(
                  #   # column(4,
                  #   #        numericInput(ns("llm_interpretation_p_adjust_cutoff"),
                  #   #                     "P-adjust cutoff",
                  #   #                     value = 0.05,
                  #   #                     min = 0,
                  #   #                     max = 0.5),
                  #   # ),
                  #   # column(4,
                  #   #        numericInput(ns("llm_interpretation_count_cutoff"),
                  #   #                     "Count cutoff",
                  #   #                     value = 5,
                  #   #                     min = 1,
                  #   #                     max = 1000)
                  #   # )
                  # ),
                  actionButton(
                    ns("submit_llm_interpretation"),
                    "Submit",
                    class = "btn-primary",
                    style = "background-color: #d83428; color: white;"
                  ),

                  actionButton(
                    ns("go2data_visualization"),
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
                ## Output layout ====
                column(8,
                       tabsetPanel(
                         tabPanel(
                           title = "Interpretation results",
                           uiOutput(ns("module_details")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_llm_interpretation_result"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "Full Prompt",
                           # h3("LLM Prompt Used"),
                           uiOutput(ns("prompt")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_llm_interpretation_prompt"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
                         ),
                         tabPanel(
                           title = "R object",
                           verbatimTextOutput(ns("llm_interpreted_functional_module_object")),
                           br(),
                           shinyjs::useShinyjs(),
                           downloadButton(ns("download_llm_interpreted_functional_module_object"),
                                          "Download",
                                          class = "btn-primary",
                                          style = "background-color: #d83428; color: white;")
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
      ns <- session$ns

      if (!is.reactive(enriched_functional_module)) {
        enriched_functional_module <- reactiveVal(NULL)
      }

      ## Section1: Load enriched_functional_module.rda and navigate to specified directory ====
      observeEvent(
        input$upload_enriched_functional_module, {
          if (!is.null(input$upload_enriched_functional_module$datapath)) {
            message("Loading data")
            tempEnv <- new.env()
            load(input$upload_enriched_functional_module$datapath,
                 envir = tempEnv)

            names <- ls(tempEnv)

            if (length(names) == 1) {
              enriched_functional_module(get(names[1], envir = tempEnv))
              # uploaded_enriched_functional_module(get(names[1], envir = tempEnv))
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
      })

      ## Display interpreted result ---------------------------------------------------------------------------
      observeEvent(input$upload_interpretation_result, {
        tmp <- new.env()
        load(input$upload_interpretation_result$datapath, envir = tmp)

        # The .rda created by mapa::llm_interpret_module() should contain exactly
        # ONE object (a named list of modules).  Guard against user mistakes:
        if (length(ls(tmp)) == 1) {
          interpreted_functional_module <- get(ls(tmp), envir = tmp)
          annotation_result(interpreted_functional_module@llm_module_interpretation)
          llm_interpreted_functional_module(interpreted_functional_module)

          showNotification("Interpretation result loaded successfully!",
                           type = "message")
        } else {
          showModal(modalDialog(
            title   = "Error",
            "The uploaded .rda must contain exactly one object
       (what `mapa::llm_interpret_module()` saves).",
            easyClose = TRUE,
            footer    = modalButton("Close")
          ))
        }
      })
      ##---------------------------------------------------------------------------

      ### Set up shinyFiles directory selection
      volumes <- shinyFiles::getVolumes()
      ### Define reactive values to store the selected directories
      selected_dirs <- reactiveValues(
        local_corpus_dir = NULL,
        embedding_output_dir = NULL
      )
      ### Set up the directory chooser for the local corpus directory
      shinyFiles::shinyDirChoose(
        input,
        "local_corpus_dir",
        roots = volumes(),
        session = session,
        restrictions = system.file(package = "base")
      )
      ### Observer for the corpus directory selection
      observeEvent(input$local_corpus_dir, {
        if (!is.null(input$local_corpus_dir)) {
          # Get the directory path
          local_corpus_path <- shinyFiles::parseDirPath(volumes(), input$local_corpus_dir)

          # Update the reactive value
          if (length(local_corpus_path) > 0) {
            selected_dirs$local_corpus_dir <- local_corpus_path

            # Optional: Display a confirmation message
            showNotification(
              paste("Selected local corpus directory:", local_corpus_path),
              type = "message"
            )
          }
        }
      })
      ### Set up the directory chooser for the embedding output directory
      shinyFiles::shinyDirChoose(
        input,
        "embedding_output_dir",
        roots = volumes(),
        session = session,
        restrictions = system.file(package = "base")
      )
      ### Observer for the embedding directory selection
      observeEvent(input$embedding_output_dir, {
        if (!is.null(input$embedding_output_dir)) {
          # Get the directory path
          embedding_output_path <- shinyFiles::parseDirPath(volumes(), input$embedding_output_dir)

          # Update the reactive value
          if (length(embedding_output_path) > 0) {
            selected_dirs$embedding_output_dir <- embedding_output_path

            # Optional: Display a confirmation message
            showNotification(
              paste("Selected embeddings out directory:", embedding_output_path),
              type = "message"
            )
          }
        }
      })

      ### show pathway
      observeEvent(selected_dirs$local_corpus_dir, {
        if (!is.null(selected_dirs$local_corpus_dir) && length(selected_dirs$local_corpus_dir) > 0) {
          updateTextInput(session, "corpus_path_display", value = selected_dirs$local_corpus_dir)
        }
      })

      observeEvent(selected_dirs$embedding_output_dir, {
        if (!is.null(selected_dirs$embedding_output_dir) && length(selected_dirs$embedding_output_dir) > 0) {
          updateTextInput(session, "embedding_path_display", value = selected_dirs$embedding_output_dir)
        }
      })

      ## Define annotation result as reactive values
      annotation_result <- reactiveVal()
      llm_interpreted_functional_module <- reactiveVal()
      llm_interpretation_code <- reactiveVal()
      module_prompt <- reactiveVal()

      # Set up the future plan - this determines how parallel tasks will run
      future::plan(future::multisession)

      observe({
        req(input$submit_llm_interpretation, enriched_functional_module())

        message("Interpreting functional modules in progress. This comprehensive analysis requires some time...")

        library(future)
        library(promises)
        library(mapa)

        object <- enriched_functional_module()
        llm_model <- input$llm_model
        embedding_model <- input$embedding_model
        api_key <- input$api_key
        embedding_output_dir <- selected_dirs$embedding_output_dir
        local_corpus_dir <- selected_dirs$local_corpus_dir
        phenotype <- if(input$phenotype == "NULL") NULL else input$phenotype
        years <- input$years

        # Show a modal with a spinner to indicate work is happening
        showModal(modalDialog(
          title = "Analysis in Progress",
          "The LLM interpretation is running in the background. Results will appear when ready.",
          footer = modalButton("Close"),
          easyClose = FALSE,
          size = "m"
        ))

        if (is.null(enriched_functional_module())) {
          removeModal()
          showModal(modalDialog(
            title = "Warning",
            "No enriched functional module data available. Please complete the previous steps or upload the data",
            easyClose = TRUE,
            footer = modalButton("Close")
          ))
          return()
        }

        # Run the interpretation asynchronously
        promises::future_promise({
          library(mapa)
          # This code runs in a separate R process
          result <-
            mapa::llm_interpret_module(
            object = object,
            llm_model = llm_model,
            embedding_model = embedding_model,
            api_key = api_key,
            embedding_output_dir = embedding_output_dir,
            local_corpus_dir = local_corpus_dir,
            phenotype = phenotype,
            years = years
          )
        },
        globals = TRUE) %...>%
          # This runs when the future completes successfully
          (function(result) {
            # Store the results
            llm_interpreted_functional_module(result)
            annotation_result(result@llm_module_interpretation)

            # # Show success notification
            showNotification("LLM interpretation completed successfully!", type = "message")
          }) %...!%
          # This runs if the future encounters an error
          (function(error) {
            # Close the modal
            removeModal()

            # Show error message
            showModal(modalDialog(
              title = "Error",
              HTML(paste("An error occurred during LLM interpretation:<br><pre>",
                         error$message, "</pre>")),
              easyClose = TRUE,
              footer = modalButton("Close")
            ))
          })
      })

      output$module_details <- renderUI({
        req(annotation_result())

        tagList(
          selectInput(
            inputId = ns("module_selector"),
            label = "Select Functional Module:",
            choices = names(annotation_result()),
            selected = names(annotation_result())[1]
          ),
          hr(),
          h3("Module Information"),
          strong("Module Name:"),
          textOutput(ns("module_name"), container = span),
          br(),
          strong("Module Summary:"),
          textOutput(ns("module_summary"), container = span),
          br(),
          strong("Association With Phenotype:"),
          textOutput(ns("association_summary"), container = span),
          br(),
          strong("Confidence Score:"),
          textOutput(ns("confidence_score"), container = span)
        )
      })

      observeEvent(                          # fire when either changes
        list(annotation_result(), input$module_selector),
        {
          req(annotation_result(), input$module_selector)

          if (input$module_selector %in% names(annotation_result())) {
            module_prompt(
              annotation_result()[[input$module_selector]]$generated_name$prompt
            )
          }
        },
        ignoreInit = TRUE                    # skip the very first (empty) run
      )

      observeEvent(annotation_result(), {
        req(annotation_result())
        ### Save code
        interpretation_code <-
          functional_module_annotation_code <-
          sprintf(
            "
            functional_module_annotation <-
              llm_interpret_module(
                object              = enriched_functional_module,
                llm_model           = %s,
                embedding_model     = %s,
                api_key             = %s,
                embedding_output_dir = %s,
                local_corpus_dir    = %s,
                phenotype           = %s,
                years               = %s
              )
            ",
            ## wrap character inputs in quotes:
            paste0('"', input$llm_model, '"'),
            paste0('"', input$embedding_model, '"'),
            paste0('"', input$api_key, '"'),
            paste0('"', selected_dirs$embedding_output_dir, '"'),
            if (is.null(selected_dirs$local_corpus_dir)) NULL else (paste0('"', selected_dirs$local_corpus_dir, '"')),
            paste0('"', input$phenotype, '"'),
            input$years
          )
        llm_interpretation_code(interpretation_code)
      })

      output$module_name <- renderText({
        req(annotation_result(), input$module_selector)
        tryCatch(
          annotation_result()[[input$module_selector]]$generated_name$module_name,
          error = function(e)
            paste("Error details:", e$message)
        )
      })

      output$module_summary <- renderText({
        req(annotation_result(), input$module_selector)
        tryCatch(
          annotation_result()[[input$module_selector]]$generated_name$summary,
          error = function(e)
            paste("Error details:", e$message)
        )
      })

      output$association_summary <- renderText({
        req(annotation_result(), input$module_selector)
        tryCatch(
          annotation_result()[[input$module_selector]]$generated_name$phenotype_analysis,
          error = function(e)
            paste("Error details:", e$message)
        )
      })

      output$confidence_score <- renderText({
        req(annotation_result(), input$module_selector)
        tryCatch(
          annotation_result()[[input$module_selector]]$generated_name$confidence_score,
          error = function(e)
            paste("Error details:", e$message)
        )
      })

      ### Show object
      output$llm_interpreted_functional_module_object <-
        renderText({
          req(llm_interpreted_functional_module())
          llm_interpreted_functional_module_obj <- llm_interpreted_functional_module()
          captured_output1 <-
            capture.output(llm_interpreted_functional_module_obj,
                           type = "message")
          captured_output2 <-
            capture.output(llm_interpreted_functional_module_obj,
                           type = "output")
          captured_output <-
            c(captured_output1,
              captured_output2)
          paste(captured_output, collapse = "\n")
        })

      ### Show code
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

      ### Format the prompt for display
      output$prompt <- renderUI({
        req(module_prompt())

        moduleID_display <- tags$div(
          tags$strong("Module ID:"),
          textOutput(ns("module_id"), inline = TRUE)
        )

        all_content <- paste(
          sapply(module_prompt(), function(item) {
            item$content
          }),
          collapse = "\n\n---\n\n"
        )

        content_html <- HTML(markdown::markdownToHTML(
          text = all_content,
          fragment.only = TRUE
        ))

        tagList(
          h3("LLM Prompt Used"),
          moduleID_display,
          br(),
          content_html
        )
      })

      output$module_id <- renderText({
        req(annotation_result(), input$module_selector)
        tryCatch(
          input$module_selector,
          error = function(e)
            paste("Error details:", e$message)
        )
      })

      observe({
        if (is.null(annotation_result()) ||
            length(annotation_result()) == 0) {
          shinyjs::disable("download_llm_interpretation_result")
          shinyjs::disable("download_llm_interpreted_functional_module_object")
        } else {
          shinyjs::enable("download_llm_interpretation_result")
          shinyjs::enable("download_llm_interpreted_functional_module_object")
        }
        if (is.null(module_prompt()) ||
            length(module_prompt()) == 0) {
          shinyjs::disable("download_llm_interpretation_prompt")
        } else {
          shinyjs::enable("download_llm_interpretation_prompt")
        }
      })

      # Download handler for llm interpreted functional modules object
      output$download_llm_interpreted_functional_module_object <-
        shiny::downloadHandler(
          filename = function() {
            "llm_interpreted_functional_module.rda"
          },
          content = function(file) {
            llm_interpreted_functional_module <-
              llm_interpreted_functional_module()
            save(llm_interpreted_functional_module, file = file)
          }
        )

      # Download handler for llm interpretation for a specific functional module ====
      output$download_llm_interpretation_result <- downloadHandler(
        filename = function() {
          # e.g. “FM3_interpretation.txt”
          paste0(input$module_selector, "_interpretation.txt")
        },
        content = function(file) {
          req(annotation_result(),
              input$module_selector)

          mod_id   <- input$module_selector
          mod_info <- annotation_result()[[mod_id]]$generated_name

          # Compose the text lines
          txt <- c(
            paste0("Module ID: ",        mod_id),
            paste0("Module Name: ",      mod_info$module_name),
            paste0("Module Summary: ",   mod_info$summary),
            paste0("Association Summary: ", mod_info$phenotype_analysis),
            paste0("Confidence Score: ", mod_info$confidence_score)
          )

          writeLines(txt, file, useBytes = TRUE)
        }
      )


      # Download handler for the prompt
      output$download_llm_interpretation_prompt <- downloadHandler(
        filename = function() {
          "llm_interpretation_prompt.txt"
        },
        content = function(file) {
          req(module_prompt())

          if (is.list(module_prompt())) {
            text_content <- paste(
              vapply(module_prompt(), function(item) item$content, character(1)),
              collapse = "\n\n---\n\n"
            )
          }

          writeLines(text_content, file, useBytes = TRUE)
        }
      )

      # # Define enriched_functional_module as a reactive value
      # llm_interpretation_result <- reactiveVal("")
      #
      # llm_interpretation_code <- reactiveVal()
      #
      # openai_key <- reactiveVal()
      #
      # observeEvent(input$submit_llm_interpretation, {
      #   openai_key1 <-
      #     Sys.getenv("chatgpt_api_key")
      #
      #   openai_key2 <-
      #     input$openai_key
      #
      #   # Check if enriched_modules is available
      #   if (is.null(enriched_functional_module()) ||
      #       length(enriched_functional_module()) == 0) {
      #     showModal(
      #       modalDialog(
      #         title = "Warning",
      #         "No enriched functional modules data available.",
      #         easyClose = TRUE,
      #         footer = modalButton("Close")
      #       )
      #     )
      #   } else {
      #     if (openai_key1 != "") {
      #       openai_key(openai_key1)
      #     } else{
      #       if (openai_key2 != "") {
      #         openai_key(openai_key2)
      #       } else{
      #         openai_key("")
      #       }
      #     }
      #
      #     if (openai_key() == "") {
      #       showModal(
      #         modalDialog(
      #           title = "Warning",
      #           "No OpenAI Key provided. No interpretation will be generated.",
      #           easyClose = TRUE,
      #           footer = modalButton("Close")
      #         )
      #       )
      #     } else{
      #       set_chatgpt_api_key(api_key = openai_key())
      #     }
      #
      #     # shinyjs::show("loading")
      #
      #     withProgress(message = 'Analysis in progress...', {
      #       tryCatch({
      #         llm_interpretation_result <-
      #           interpret_pathways(
      #             object = enriched_functional_module(),
      #             p.adjust.cutoff = input$llm_interpretation_p_adjust_cutoff,
      #             disease = input$llm_interpretation_disease,
      #             count.cutoff = input$llm_interpretation_count_cutoff,
      #             top_n = input$llm_interpretation_top_n
      #           )
      #         llm_interpretation_result(llm_interpretation_result)
      #       },
      #       error = function(e) {
      #         showModal(modalDialog(
      #           title = "Error",
      #           paste("Details:", e$message),
      #           easyClose = TRUE,
      #           footer = modalButton("Close")
      #         ))
      #         llm_interpretation_result("No result")
      #       })
      #     })
      #     # shinyjs::hide("loading")
      #
      #     ##save code
      #     llm_interpretation_code <-
      #       sprintf(
      #         '
      #       llm_interpretation_result <-
      #       interpret_pathways(
      #       object = enriched_functional_module,
      #       p.adjust.cutoff = %s,
      #       disease = %s,
      #       count.cutoff = %s,
      #       top_n = %s)
      #       ',
      #         input$llm_interpretation_p_adjust_cutoff,
      #         paste0('"', input$llm_interpretation_disease, '"'),
      #         input$llm_interpretation_count_cutoff,
      #         input$llm_interpretation_top_n
      #       )
      #
      #     llm_interpretation_code(llm_interpretation_code)
      #   }
      # })
      #
      # output$llm_interpretation_result <-
      #   renderUI({
      #     shiny::HTML(markdown::markdownToHTML(llm_interpretation_result(),
      #                                          fragment.only = TRUE))
      #   })
      #
      #
      # output$llm_enriched_functional_modules1 <-
      #   shiny::renderDataTable({
      #     req(tryCatch(
      #       enriched_functional_module()@merged_module$functional_module_result,
      #       error = function(e)
      #         NULL
      #     ))
      #   },
      #   options = list(pageLength = 10,
      #                  scrollX = TRUE))
      #
      # output$llm_enriched_functional_modules2 <-
      #   shiny::renderDataTable({
      #     req(tryCatch(
      #       enriched_functional_module()@merged_module$result_with_module,
      #       error = function(e)
      #         NULL
      #     ))
      #   },
      #   options = list(pageLength = 10,
      #                  scrollX = TRUE))
      #
      # observe({
      #   if (is.null(llm_interpretation_result()) ||
      #       length(llm_interpretation_result()) == 0 ||
      #       llm_interpretation_result() == "") {
      #     shinyjs::disable("download_llm_interpretation_result")
      #   } else {
      #     shinyjs::enable("download_llm_interpretation_result")
      #
      #   }
      # })
      #
      #

      ###Go to go2data_visualization tab
      ####if there is not enriched_functional_module, show a warning message
      observeEvent(input$go2data_visualization, {
        # Check if enriched_functional_module is available
        if ((is.null(enriched_functional_module()) ||
            length(enriched_functional_module()) == 0) &&
            is.null(llm_interpreted_functional_module())) {
          showModal(
            modalDialog(
              title = "Warning",
              "No enriched functional modules data available.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else {
          # User never pressed “Submit” in this tab
          if (is.null(llm_interpreted_functional_module()) ||
              length(llm_interpreted_functional_module()) == 0) {

            llm_interpreted_functional_module(
              enriched_functional_module()
            )
          }

          tab_switch("data_visualization")
        }
      })

      return(llm_interpreted_functional_module)
    }
  )
}
