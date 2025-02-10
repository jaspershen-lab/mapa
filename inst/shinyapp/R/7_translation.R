# Translation tab ui====
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


### Step 5 Translation ----
#
# translation_code <-
#   reactiveVal()
#
# translation_model_ai_key <-
#   reactiveVal()
#
# enriched_functional_module2 <-
#   reactiveVal()
#
# #####if the user translate the object
# observeEvent(input$submit_translation, {
#   openai_key1 <-
#     Sys.getenv("chatgpt_api_key")
#
#   openai_key2 <-
#     input$translation_model_ai_key
#
#   gemini_key1 <-
#     Sys.getenv("gemini_api_key")
#
#   gemini_key2 <-
#     input$translation_model_ai_key
#
#   if (input$translation_model == "chatgpt") {
#     translation_model_ai_key <-
#       openai_key1
#
#     if (openai_key1 != "") {
#       translation_model_ai_key(openai_key1)
#     } else{
#       if (openai_key2 != "") {
#         translation_model_ai_key(openai_key2)
#       } else{
#         translation_model_ai_key("")
#       }
#     }
#
#     if (translation_model_ai_key() == "") {
#       showModal(
#         modalDialog(
#           title = "Warning",
#           "No OpenAI Key provided. No translation will be generated.",
#           easyClose = TRUE,
#           footer = modalButton("Close")
#         )
#       )
#     } else{
#       mapa::set_chatgpt_api_key(api_key = translation_model_ai_key())
#     }
#
#   } else{
#     if (gemini_key1 != "") {
#       translation_model_ai_key(gemini_key1)
#     } else{
#       if (gemini_key2 != "") {
#         translation_model_ai_key(gemini_key2)
#       } else{
#         translation_model_ai_key("")
#       }
#     }
#
#     if (translation_model_ai_key() == "") {
#       showModal(
#         modalDialog(
#           title = "Warning",
#           "No Gemini Key provided. No translation will be generated.",
#           easyClose = TRUE,
#           footer = modalButton("Close")
#         )
#       )
#     } else{
#       mapa::set_gemini_api_key(api_key = translation_model_ai_key())
#     }
#   }
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
#     # shinyjs::show("loading")
#
#     withProgress(message = 'Analysis in progress...', {
#       tryCatch({
#         enriched_functional_module <-
#           translate_language(
#             text = enriched_functional_module(),
#             engine = input$translation_model,
#             to = input$translation_to
#           )
#       },
#       error = function(e) {
#         showModal(modalDialog(
#           title = "Error",
#           paste("Details:", e$message),
#           easyClose = TRUE,
#           footer = modalButton("Close")
#         ))
#       })
#     })
#
#     enriched_functional_module(enriched_functional_module)
#     enriched_functional_module2(enriched_functional_module)
#
#     # shinyjs::hide("loading")
#
#     ##save code
#     translation_code <-
#       sprintf(
#         '
#         enriched_functional_module <-
#         translate_language(text = enriched_functional_module,
#                            engine = %s,
#                            to = %s)
#         ',
#         paste0('"', input$translation_model, '"'),
#         paste0('"', input$translation_to, '"')
#       )
#
#     translation_code(translation_code)
#   }
# })
#
#
# output$enriched_functional_module_object2 <-
#   renderText({
#     req(enriched_functional_module2())
#     enriched_functional_module <- enriched_functional_module2()
#     captured_output1 <-
#       capture.output(enriched_functional_module,
#                      type = "message")
#     captured_output2 <-
#       capture.output(enriched_functional_module,
#                      type = "output")
#     captured_output <-
#       c(captured_output1,
#         captured_output2)
#     paste(captured_output, collapse = "\n")
#   })
#
#
# ####show code
# observeEvent(input$show_translation_code, {
#   if (is.null(translation_code()) ||
#       length(translation_code()) == 0) {
#     showModal(
#       modalDialog(
#         title = "Warning",
#         "No available code",
#         easyClose = TRUE,
#         footer = modalButton("Close")
#       )
#     )
#   } else{
#     code_content <-
#       translation_code()
#     code_content <-
#       paste(code_content, collapse = "\n")
#     showModal(modalDialog(
#       title = "Code",
#       tags$pre(code_content),
#       easyClose = TRUE,
#       footer = modalButton("Close")
#     ))
#   }
# })
#
#
# output$download_enriched_functional_module_object2 <-
#   shiny::downloadHandler(
#     filename = function() {
#       "enriched_functional_module.rda"
#     },
#     content = function(file) {
#       enriched_functional_module <-
#         enriched_functional_module2()
#       save(enriched_functional_module,
#            file = file)
#     }
#   )
#
# observe({
#   if (is.null(enriched_functional_module2()) ||
#       length(enriched_functional_module2()) == 0) {
#     shinyjs::disable("download_enriched_functional_module_object2")
#   } else {
#     shinyjs::enable("download_enriched_functional_module_object2")
#   }
# })
#
# ###Go to data visualization tab
# ####if there is not enriched_functional_module, show a warning message
# observeEvent(input$go2data_visualization, {
#   # Check if enriched_functional_module is available
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
#     updateTabItems(session = session,
#                    inputId = "tabs",
#                    selected = "data_visualization")
#     if ("enrich_pathway" %in% names(enriched_functional_module()@process_info)) {
#       all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
#     } else {
#       all_choices <- c("NES")
#     }
#     updateSelectInput(
#       session,
#       "x_axis_name",
#       choices = all_choices
#     )
#   }
# })
#
#
# ###if users click skip translation, go to data visualization tab
# observeEvent(input$skip_translation, {
#   # Check if enriched_functional_module is available
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
#     updateTabItems(session = session,
#                    inputId = "tabs",
#                    selected = "data_visualization")
#
#     if ("enrich_pathway" %in% names(enriched_functional_module()@process_info)) {
#       all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
#     } else {
#       all_choices <- c("NES")
#     }
#     updateSelectInput(
#       session,
#       "x_axis_name",
#       choices = all_choices
#     )
#   }
# })
