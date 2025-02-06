server <-
  function(input, output, session) {

    tab_switch <- reactiveVal()
    # Observer to handle tab switching
    observe({
      req(tab_switch())
      updateTabItems(session, "tabs", selected = tab_switch())
    })

    ### Step 1: Upload data ----
    #upload_data_result <- reactive(upload_data_server("upload_data_tab"))
    upload_data_result <- upload_data_server("upload_data_tab", tab_switch)

    ### Step 2 enrich pathways ----
    enriched_pathways_result <- enrich_pathway_server("enrich_pathway_tab", variable_info = upload_data_result, tab_switch)

    ### Step 3 merge pathways ----
    merge_pathways_server("merge_pathways_tab", enriched_pathways = enriched_pathways_result, tab_switch)


    # ###Step 4 merge modules ====
    # # Define enriched_functional_module as a reactive value
    # enriched_functional_module <-
    #   reactiveVal()
    #
    # merge_modules_code <-
    #   reactiveVal()
    #
    # observeEvent(input$submit_merge_modules, {
    #   # Check if enriched_modules is available
    #   if (is.null(enriched_modules()) ||
    #       length(enriched_modules()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched modules data available. Please 'Enrich modules' first.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     # shinyjs::show("loading")
    #     withProgress(message = 'Analysis in progress...', {
    #       tryCatch({
    #         library(clusterProfiler)
    #         library(org.Hs.eg.db)
    #         library(ReactomePA)
    #
    #         result <-
    #           merge_modules(
    #             object = enriched_modules(),
    #             sim.cutoff = input$sim.cutoff.module,
    #             measure_method = input$measure.method.module,
    #             path = "result",
    #             save_to_local = FALSE
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
    #     enriched_functional_module(result)
    #     # shinyjs::hide("loading")
    #
    #     ##save code
    #     merge_modules_code <-
    #       sprintf(
    #         '
    #         enriched_functional_module <-
    #         merge_modules(
    #         object = enriched_modules,
    #         sim.cutoff = %s,
    #         measure_method = %s)
    #         ',
    #         input$sim.cutoff.module,
    #         paste0('"', input$measure.method.module, '"')
    #       )
    #
    #     merge_modules_code(merge_modules_code)
    #   }
    # })
    #
    #
    # output$enriched_functional_module_object <-
    #   renderText({
    #     req(enriched_functional_module())
    #     enriched_functional_module <- enriched_functional_module()
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
    # output$enriched_functional_modules <-
    #   DT::renderDT({
    #     data <- tryCatch(
    #       enriched_functional_modul()@merged_module$functional_module_result,
    #       error = function(e) NULL)
    #
    #     req(data)
    #
    #     DT::datatable(
    #       data,
    #       options = list(pageLength = 10,
    #                      scrollX = TRUE)
    #     )
    #   })
    #
    # output$download_enriched_functional_modules <-
    #   shiny::downloadHandler(
    #     filename = function() {
    #       "merged_modules.csv"
    #     },
    #     content = function(file) {
    #       write.csv(
    #         enriched_functional_module()@merged_module$functional_module_result,
    #         file,
    #         row.names = FALSE
    #       )
    #     }
    #   )
    #
    # observe({
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     shinyjs::disable("download_enriched_functional_modules")
    #   } else {
    #     if (length(enriched_functional_module()@merged_module) == 0) {
    #       shinyjs::disable("download_enriched_functional_modules")
    #     } else{
    #       shinyjs::enable("download_enriched_functional_modules")
    #     }
    #   }
    # })
    #
    # ######data visualization
    # ###define enirched_functional_module_plot
    # enirched_functional_module_plot <-
    #   reactiveVal()
    # observeEvent(input$generate_enirched_functional_module, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched functional modules data available. Please 'Merge modules' first.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     # shinyjs::show("loading")
    #     withProgress(message = 'Analysis in progress...', {
    #       tryCatch(
    #         plot <-
    #           plot_similarity_network(
    #             object = enriched_functional_module(),
    #             level = "functional_module",
    #             degree_cutoff = input$enirched_functional_moduleplot__degree_cutoff,
    #             text = input$enirched_functional_module_plot_text,
    #             text_all = input$enirched_functional_module_plot_text_all
    #           ),
    #         error = function(e) {
    #           showModal(
    #             modalDialog(
    #               title = "Error",
    #               "Please check your input parameters.",
    #               easyClose = TRUE,
    #               footer = modalButton("Close")
    #             )
    #           )
    #         }
    #       )
    #     })
    #
    #     enirched_functional_module_plot(plot)
    #     # shinyjs::hide("loading")
    #   }
    # })
    #
    # output$enirched_functional_module_plot <-
    #   shiny::renderPlot({
    #     req(tryCatch(
    #       enirched_functional_module_plot(),
    #       error = function(e)
    #         NULL
    #     ))
    #   })
    #
    #
    # output$download_enriched_functional_module_object <-
    #   shiny::downloadHandler(
    #     filename = function() {
    #       "enriched_functional_module.rda"
    #     },
    #     content = function(file) {
    #       save(enriched_functional_module(), file = file)
    #     }
    #   )
    #
    # observe({
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     shinyjs::disable("download_enriched_functional_module_object")
    #   } else {
    #     shinyjs::enable("download_enriched_functional_module_object")
    #   }
    # })
    #
    #
    # output$download_enriched_functional_module_object <-
    #   shiny::downloadHandler(
    #     filename = function() {
    #       "enriched_functional_module.rda"
    #     },
    #     content = function(file) {
    #       enriched_functional_module <-
    #         enriched_functional_module()
    #       save(enriched_functional_module, file = file)
    #     }
    #   )
    #
    # observe({
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     shinyjs::disable("download_enriched_functional_module_object")
    #   } else {
    #     shinyjs::enable("download_enriched_functional_module_object")
    #   }
    # })
    #
    #
    # ####show code
    # observeEvent(input$show_merge_modules_code, {
    #   if (is.null(merge_modules_code()) ||
    #       length(merge_modules_code()) == 0) {
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
    #       merge_modules_code()
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
    # ###Go to Translation tab
    # ####if there is not enriched_functional_module, show a warning message
    # observeEvent(input$go2translation, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "Please merge modules first.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     updateTabItems(session = session,
    #                    inputId = "tabs",
    #                    selected = "translation")
    #   }
    # })
    #
    #
    # ###Step 5 Translation ====
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
    #
    # ###Step 6 Data visualization ====
    # # Observe file upload
    # # enriched_functional_module3 <-
    # #   reactiveVal()
    # observe({
    #   req(enriched_functional_module())
    #   if ("enrich_pathway" %in% names(enriched_functional_module()@process_info)) {
    #     all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
    #   } else {
    #     all_choices <- c("NES")
    #   }
    #   updateSelectInput(
    #     session,
    #     "x_axis_name",
    #     choices = all_choices
    #   )
    # })
    #
    # observeEvent(input$upload_enriched_functional_module, {
    #   if (!is.null(input$upload_enriched_functional_module$datapath)) {
    #     message("Loading data")
    #     tempEnv <- new.env()
    #     load(input$upload_enriched_functional_module$datapath,
    #          envir = tempEnv)
    #
    #     names <- ls(tempEnv)
    #
    #     if (length(names) == 1) {
    #       # If enriched_functional_module is another reactiveVal, uncomment the next line
    #       enriched_functional_module(get(names[1], envir = tempEnv))
    #     } else {
    #       message("The .rda file does not contain exactly one object.")
    #       showModal(
    #         modalDialog(
    #           title = "Error",
    #           "The uploaded file should contain exactly one object.",
    #           easyClose = TRUE,
    #           footer = modalButton("Close")
    #         )
    #       )
    #     }
    #   }
    #
    #   if ("enrich_pathway" %in% names(enriched_functional_module()@process_info)) {
    #     all_choices <- c("qscore", "RichFactor", "FoldEnrichment")
    #   } else {
    #     all_choices <- c("NES")
    #   }
    #   updateSelectInput(
    #     session,
    #     "x_axis_name",
    #     choices = all_choices
    #   )
    # })
    #
    # #####barplot =====
    # # Observe generate barplot button click
    # barplot <-
    #   reactiveVal()
    # barplot_code <-
    #   reactiveVal()
    #
    # observeEvent(input$generate_barplot, {
    #   message("generating barplot")
    #   if (is.null(enriched_functional_module())) {
    #     # No enriched functional module available
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched functional module data available. Please complete the previous steps or upload the data",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     # shinyjs::show("loading")
    #     withProgress(message = 'Analysis in progress...', {
    #       tryCatch(
    #         plot <-
    #           plot_pathway_bar(
    #             object = enriched_functional_module(),
    #             top_n = input$barplot_top_n,
    #             x_axis_name = input$x_axis_name,
    #             y_label_width = input$barplot_y_lable_width,
    #             p.adjust.cutoff = input$barplot_p_adjust_cutoff,
    #             count.cutoff = input$barplot_count_cutoff,
    #             level = input$barplot_level,
    #             database = input$barplot_database,
    #             line_type = input$line_type,
    #             database_color = c(
    #               GO = input$barplot_go_color,
    #               KEGG = input$barplot_kegg_color,
    #               Reactome = input$barplot_reactome_color
    #             ),
    #             translation = input$barplot_translation
    #           ),
    #         error = function(e) {
    #           showModal(modalDialog(
    #             title = "Error",
    #             paste("Details:", e$message),
    #             easyClose = TRUE,
    #             footer = modalButton("Close")
    #           ))
    #         }
    #       )
    #     })
    #
    #     # shinyjs::hide("loading")
    #
    #     barplot(plot)
    #
    #     ###save code
    #     data_color <-
    #       paste0("c(", paste(paste(
    #         c("GO", "KEGG", "Reactome"),
    #         c(
    #           paste0('"', input$barplot_go_color, '"'),
    #           paste0('"', input$barplot_kegg_color, '"'),
    #           paste0('"', input$barplot_reactome_color, '"')
    #         ),
    #         sep = " = "
    #       ),
    #       collapse = ", "), ")")
    #
    #     barplot_database <-
    #       paste0("c(", paste(unlist(lapply(paste(input$barplot_database), function(x)
    #         paste0('"', x, '"'))),
    #         collapse = ", "), ")")
    #
    #     barplot_code <-
    #       sprintf(
    #         "
    #         plot_pathway_bar(
    #         object = enriched_functional_module,
    #         top_n = %s,
    #         y_lable_width = %s,
    #         p.adjust.cutoff = %s,
    #         count.cutoff = %s,
    #         level = %s,
    #         database = %s,
    #         line_type = %s,
    #         database_color = %s)
    #         ",
    #         input$barplot_top_n,
    #         input$barplot_y_lable_width,
    #         input$barplot_p_adjust_cutoff,
    #         input$barplot_count_cutoff,
    #         paste0('"', input$barplot_level, '"'),
    #         barplot_database,
    #         paste0('"', input$line_type, '"'),
    #         data_color
    #       )
    #
    #     barplot_code(barplot_code)
    #   }
    # })
    #
    #
    # output$barplot <-
    #   renderPlot({
    #     req(barplot())
    #     barplot()
    #   },
    #   res = 96)
    #
    # # output$barplot <-
    # #   renderPlot({
    # #     req(barplot())
    # #     barplot()
    # #   },
    # #   width = function() {
    # #     input$barplot_width_show
    # #   },
    # #   height = function() {
    # #     input$barplot_height_show
    # #   })
    #
    # output$download_barplot <-
    #   downloadHandler(
    #     filename = function() {
    #       paste0("pathway_barplot.", input$barplot_type)
    #     },
    #     content = function(file) {
    #       ggsave(
    #         file,
    #         plot = barplot(),
    #         width = input$barplot_width,
    #         height = input$barplot_height
    #       )
    #     }
    #   )
    #
    # observe({
    #   if (is.null(barplot()) ||
    #       length(barplot()) == 0) {
    #     shinyjs::disable("download_barplot")
    #   } else {
    #     shinyjs::enable("download_barplot")
    #   }
    # })
    #
    # ####show code
    # observeEvent(input$show_barplot_code, {
    #   if (is.null(barplot_code()) ||
    #       length(barplot_code()) == 0) {
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
    #       barplot_code()
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
    #
    # #####module_similarity_network
    # # Observe generate module_similarity_network button click
    # module_similarity_network <-
    #   reactiveVal()
    #
    # module_similarity_network_code <-
    #   reactiveVal()
    #
    # observeEvent(input$generate_module_similarity_network, {
    #   if (is.null(enriched_functional_module())) {
    #     # No enriched functional module available
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched functional module data available. Please complete the previous steps or upload the data.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     # shinyjs::show("loading")
    #
    #     withProgress(message = 'Analysis in progress...', {
    #       tryCatch(
    #         plot <-
    #           plot_similarity_network(
    #             object = enriched_functional_module(),
    #             level = input$module_similarity_network_level,
    #             database = input$module_similarity_network_database,
    #             degree_cutoff = input$module_similarity_network_degree_cutoff,
    #             text = input$module_similarity_network_text,
    #             text_all = input$module_similarity_network_text_all,
    #             translation = input$module_similarity_network_translation
    #           ),
    #         error = function(e) {
    #           showModal(
    #             modalDialog(
    #               title = "Error",
    #               "Please check your input parameters.",
    #               easyClose = TRUE,
    #               footer = modalButton("Close")
    #             )
    #           )
    #         }
    #       )
    #     })
    #
    #     # shinyjs::hide("loading")
    #
    #     module_similarity_network(plot)
    #
    #     ###save code
    #     module_similarity_network_code <-
    #       sprintf(
    #         '
    #         plot_similarity_network(
    #         object = enriched_functional_module,
    #         level = %s,
    #         database = %s,
    #         degree_cutoff = %s,
    #         text = %s,
    #         text_all = %s)
    #         ',
    #         paste0('"', input$module_similarity_network_level, '"'),
    #         paste0('"', input$module_similarity_network_database, '"'),
    #         input$module_similarity_network_degree_cutoff,
    #         input$module_similarity_network_text,
    #         input$module_similarity_network_text_all
    #       )
    #
    #     module_similarity_network_code(module_similarity_network_code)
    #
    #   }
    # })
    #
    # output$module_similarity_network <-
    #   renderPlot({
    #     req(module_similarity_network())
    #     module_similarity_network()
    #   },
    #   res = 96)
    #
    #
    # ######code for module_similarity_network
    # ####show code
    # observeEvent(input$show_module_similarity_network_code, {
    #   if (is.null(module_similarity_network_code()) ||
    #       length(module_similarity_network_code()) == 0) {
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
    #       module_similarity_network_code()
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
    # output$download_module_similarity_network <-
    #   downloadHandler(
    #     filename = function() {
    #       paste0(
    #         "module_similarity_network_",
    #         ifelse(
    #           input$module_similarity_network_level == "module",
    #           input$module_similarity_network_database,
    #           "functional_module"
    #         ),
    #         ".",
    #         input$module_similarity_network_type
    #       )
    #     },
    #     content = function(file) {
    #       ggsave(
    #         file,
    #         plot = module_similarity_network(),
    #         width = input$module_similarity_network_width,
    #         height = input$module_similarity_network_height
    #       )
    #     }
    #   )
    #
    # observe({
    #   if (is.null(module_similarity_network()) ||
    #       length(module_similarity_network()) == 0) {
    #     shinyjs::disable("download_module_similarity_network")
    #   } else {
    #     shinyjs::enable("download_module_similarity_network")
    #   }
    # })
    #
    #
    #
    #
    # ########module information plot
    # # Update the module ID
    # module_information_module_id <-
    #   reactiveVal()
    #
    # observe({
    #   if (!is.null(enriched_functional_module()) &
    #       length(enriched_functional_module()) != 0) {
    #     ####level is functional module
    #     if (input$module_information_level == "functional_module") {
    #       if (length(enriched_functional_module()@merged_module) > 0) {
    #         module_information_module_id <-
    #           unique(
    #             enriched_functional_module()@merged_module$functional_module_result$module
    #           )
    #         module_information_module_id(module_information_module_id)
    #       }
    #     }
    #
    #     ####level is module
    #     if (input$module_information_level == "module") {
    #       ####database is go
    #       if (input$module_information_database == "go") {
    #         if (length(enriched_functional_module()@merged_pathway_go) > 0) {
    #           module_information_module_id <-
    #             unique(
    #               enriched_functional_module()@merged_pathway_go$module_result$module
    #             )
    #           module_information_module_id(module_information_module_id)
    #         }
    #       }
    #
    #       ####database is kegg
    #       if (input$module_information_database == "kegg") {
    #         if (length(enriched_functional_module()@merged_pathway_kegg) > 0) {
    #           module_information_module_id <-
    #             unique(
    #               enriched_functional_module()@merged_pathway_kegg$module_result$module
    #             )
    #           module_information_module_id(module_information_module_id)
    #         }
    #       }
    #
    #       ####database is reactome
    #       if (input$module_information_database == "reactome") {
    #         if (length(enriched_functional_module()@merged_pathway_reactome) > 0) {
    #           module_information_module_id <-
    #             unique(
    #               enriched_functional_module()@merged_pathway_reactome$module_result$module
    #             )
    #           module_information_module_id(module_information_module_id)
    #         }
    #       }
    #     }
    #
    #     updateSelectInput(
    #       session,
    #       "module_information_module_id",
    #       choices = stringr::str_sort(module_information_module_id(), numeric = TRUE),
    #       selected = stringr::str_sort(module_information_module_id(), numeric = TRUE)[1]
    #     )
    #   }
    # })
    #
    # #####module information plot =====
    # # Observe generate module information button click
    # module_information <-
    #   reactiveVal()
    # module_information1 <-
    #   reactiveVal()
    # module_information2 <-
    #   reactiveVal()
    # module_information3 <-
    #   reactiveVal()
    # module_information_code <-
    #   reactiveVal()
    #
    # ####if the module_information_module_id() is null, then show warning
    # observeEvent(input$generate_module_information, {
    #   if (is.null(enriched_functional_module())) {
    #     # No enriched functional module available
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched functional module data available. Please complete the previous steps or upload the data.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     if (is.null(module_information_module_id())) {
    #       # No enriched functional module available
    #       showModal(
    #         modalDialog(
    #           title = "Warning",
    #           "Select a module ID first",
    #           easyClose = TRUE,
    #           footer = modalButton("Close")
    #         )
    #       )
    #     } else{
    #       # shinyjs::show("loading")
    #       # browser()
    #       withProgress(message = 'Analysis in progress...', {
    #         tryCatch(
    #           plot <-
    #             plot_module_info(
    #               object = enriched_functional_module(),
    #               level = input$module_information_level,
    #               database = input$module_information_database,
    #               module_id = input$module_information_module_id,
    #               translation = input$module_information_translation
    #             ),
    #           error = function(e) {
    #             showModal(
    #               modalDialog(
    #                 title = "Error",
    #                 "Please check your input parameters.",
    #                 easyClose = TRUE,
    #                 footer = modalButton("Close")
    #               )
    #             )
    #           }
    #         )
    #
    #       })
    #
    #       # shinyjs::hide("loading")
    #       if (is(plot, "ggplot")) {
    #         plot_all <-
    #           plot + plot + plot +
    #           patchwork::plot_layout(ncol = 1)
    #
    #         module_information(plot_all)
    #         module_information1(plot)
    #         module_information2(plot)
    #         module_information3(plot)
    #       } else{
    #         plot_all <-
    #           plot[[1]] + plot[[2]] + plot[[3]] +
    #           patchwork::plot_layout(ncol = 1)
    #
    #         module_information(plot_all)
    #         module_information1(plot[[1]])
    #         module_information2(plot[[2]])
    #         module_information3(plot[[3]])
    #       }
    #
    #
    #       ###save code
    #       module_information_code <-
    #         sprintf(
    #           '
    #         plot_module_info(
    #         object = enriched_functional_module,
    #         level = %s,
    #         database = %s,
    #         module_id = %s)
    #         ',
    #         paste0('"', input$module_information_level, '"'),
    #         paste0('"', input$module_information_database, '"'),
    #         paste0('"', input$module_information_module_id, '"')
    #         )
    #
    #       module_information_code(module_information_code)
    #
    #     }
    #   }
    # })
    #
    # output$module_information <-
    #   renderPlot({
    #     req(module_information())
    #     module_information()
    #   }, res = 96)
    #
    # output$module_information1 <-
    #   renderPlot({
    #     req(module_information1())
    #     module_information1()
    #   }, res = 96)
    #
    # output$module_information2 <-
    #   renderPlot({
    #     req(module_information2())
    #     module_information2()
    #   }, res = 96)
    #
    # output$module_information3 <-
    #   renderPlot({
    #     req(module_information3())
    #     module_information3()
    #   }, res = 96)
    #
    # ######code for module_information
    # ####show code
    # observeEvent(input$show_module_information_code, {
    #   if (is.null(module_information_code()) ||
    #       length(module_information_code()) == 0) {
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
    #       module_information_code()
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
    # output$download_module_information <-
    #   downloadHandler(
    #     filename = function() {
    #       paste0(
    #         "module_information_",
    #         input$module_information_module_id,
    #         ".",
    #         input$module_information_type
    #       )
    #     },
    #     content = function(file) {
    #       ggsave(
    #         file,
    #         plot = module_information(),
    #         width = input$module_information_width,
    #         height = input$module_information_height
    #       )
    #     }
    #   )
    #
    # observe({
    #   if (is.null(module_information1()) ||
    #       length(module_information1()) == 0) {
    #     shinyjs::disable("download_module_information")
    #   } else {
    #     shinyjs::enable("download_module_information")
    #   }
    # })
    #
    #
    # #####relationship network plot====
    # # Update the module ID
    # observe({
    #   if (!is.null(enriched_functional_module()) &
    #       length(enriched_functional_module()) != 0) {
    #     ####level is functional module
    #     if (input$relationship_network_level == "functional_module") {
    #       if (length(enriched_functional_module()@merged_module) > 0) {
    #         relationship_network_module_id <-
    #           unique(
    #             enriched_functional_module()@merged_module$functional_module_result$module
    #           )
    #       }
    #     }
    #
    #     ####level is module
    #     if (input$relationship_network_level == "module") {
    #       if (length(enriched_functional_module()@merged_pathway_go) > 0) {
    #         relationship_network_module_id_go <-
    #           unique(enriched_functional_module()@merged_pathway_go$module_result$module)
    #       } else{
    #         relationship_network_module_id_go <- NULL
    #       }
    #
    #       if (length(enriched_functional_module()@merged_pathway_kegg) > 0) {
    #         relationship_network_module_id_kegg <-
    #           unique(
    #             enriched_functional_module()@merged_pathway_kegg$module_result$module
    #           )
    #       } else{
    #         relationship_network_module_id_kegg <- NULL
    #       }
    #
    #       ####database is reactome
    #       if (length(enriched_functional_module()@merged_pathway_reactome) > 0) {
    #         relationship_network_module_id_reactome <-
    #           unique(
    #             enriched_functional_module()@merged_pathway_reactome$module_result$module
    #           )
    #       } else{
    #         relationship_network_module_id_reactome <- NULL
    #       }
    #
    #       relationship_network_module_id <-
    #         c(
    #           relationship_network_module_id_go,
    #           relationship_network_module_id_kegg,
    #           relationship_network_module_id_reactome
    #         )
    #     }
    #
    #     updateSelectInput(
    #       session,
    #       "relationship_network_module_id",
    #       choices = stringr::str_sort(relationship_network_module_id, numeric = TRUE),
    #       selected = stringr::str_sort(relationship_network_module_id, numeric = TRUE)[1]
    #     )
    #   }
    # })
    #
    # # Observe generate relationship network button click
    # relationship_network <-
    #   reactiveVal()
    # relationship_network_code <-
    #   reactiveVal()
    #
    # ####get the filtered enriched_functional_module
    # object <-
    #   reactiveVal()
    #
    # observeEvent(input$generate_relationship_network, {
    #   if (is.null(enriched_functional_module())) {
    #     # No enriched functional module available
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched functional module data available. Please complete the previous steps or upload the data.",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     ####if filtered by functiobal module and modules
    #     object <-
    #       enriched_functional_module()
    #     if (input$relationship_network_filter) {
    #       tryCatch({
    #         object <-
    #           filter_functional_module(
    #             object,
    #             level = input$relationship_network_level,
    #             remain_id = input$relationship_network_module_id
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
    #     }
    #
    #     object(object)
    #
    #     # shinyjs::show("loading")
    #
    #     withProgress(message = 'Analysis in progress...', {
    #       tryCatch(
    #         plot <-
    #           plot_relationship_network(
    #             object = object(),
    #             include_functional_modules = input$relationship_network_include_functional_modules,
    #             include_modules = input$relationship_network_include_modules,
    #             include_pathways = input$relationship_network_include_pathways,
    #             include_molecules = input$relationship_network_include_molecules,
    #             functional_module_text = input$relationship_network_functional_module_text,
    #             module_text = input$relationship_network_module_text,
    #             pathway_text = input$relationship_network_pathway_text,
    #             molecule_text = input$relationship_network_molecule_text,
    #             circular_plot = input$relationship_network_circular_plot,
    #             functional_module_color = input$relationship_network_functional_module_color,
    #             module_color = input$relationship_network_module_color,
    #             pathway_color = input$relationship_network_pathway_color,
    #             molecule_color = input$relationship_network_molecule_color,
    #             functional_module_arrange_position = input$relationship_network_functional_module_arrange_position,
    #             module_arrange_position = input$relationship_network_module_arrange_position,
    #             pathway_arrange_position = input$relationship_network_pathway_arrange_position,
    #             molecule_arrange_position = input$relationship_network_molecule_arrange_position,
    #             functional_module_position_limits = c(
    #               input$relationship_network_functional_module_position_limits[1],
    #               input$relationship_network_functional_module_position_limits[2]
    #             ),
    #             module_position_limits = c(
    #               input$relationship_network_module_position_limits[1],
    #               input$relationship_network_module_position_limits[2]
    #             ),
    #             pathway_position_limits = c(
    #               input$relationship_network_pathway_position_limits[1],
    #               input$relationship_network_pathway_position_limits[2]
    #             ),
    #             molecule_position_limits = c(
    #               input$relationship_network_molecule_position_limits[1],
    #               input$relationship_network_molecule_position_limits[2]
    #             ),
    #             translation = input$relationship_network_translation
    #           ),
    #         error = function(e) {
    #           showModal(
    #             modalDialog(
    #               title = "Error",
    #               "Please check your input parameters.",
    #               easyClose = TRUE,
    #               footer = modalButton("Close")
    #             )
    #           )
    #         }
    #       )
    #     })
    #
    #     # shinyjs::hide("loading")
    #
    #     relationship_network(plot)
    #
    #     ###save code
    #     relationship_network_module_id <-
    #       paste0("c(",
    #              paste0(
    #                paste0('"',
    #                       input$relationship_network_module_id,
    #                       '"'),
    #                collapse = ", "
    #              ),
    #              ")")
    #
    #     relationship_network_code1 <-
    #       sprintf(
    #         '
    #         object <-
    #         filter_functional_module(
    #           object,
    #           level = %s,
    #           remain_id = %s
    #         )
    #         ',
    #         paste0('"', input$relationship_network_level, '"'),
    #         relationship_network_module_id
    #       )
    #
    #     functional_module_position_limits <-
    #       paste0(
    #         "c(",
    #         paste0(
    #           input$relationship_network_functional_module_position_limits,
    #           collapse = ", "
    #         ),
    #         ")"
    #       )
    #
    #     module_position_limits <-
    #       paste0(
    #         "c(",
    #         paste0(
    #           input$relationship_network_module_position_limits,
    #           collapse = ", "
    #         ),
    #         ")"
    #       )
    #
    #     pathway_position_limits <-
    #       paste0(
    #         "c(",
    #         paste0(
    #           input$relationship_network_pathway_position_limits,
    #           collapse = ", "
    #         ),
    #         ")"
    #       )
    #
    #     molecule_position_limits <-
    #       paste0(
    #         "c(",
    #         paste0(
    #           input$relationship_network_molecule_position_limits,
    #           collapse = ", "
    #         ),
    #         ")"
    #       )
    #
    #     relationship_network_code2 <-
    #       sprintf(
    #         '
    #         plot_relationship_network(
    #         object = object,
    #         include_functional_modules = %s,
    #         include_modules = %s,
    #         include_pathways = %s,
    #         include_molecules = %s,
    #         functional_module_text = %s,
    #         module_text = %s,
    #         pathway_text = %s,
    #         molecule_text = %s,
    #         circular_plot = %s,
    #         functional_module_color = %s,
    #         module_color = %s,
    #         pathway_color = %s,
    #         molecule_color = %s,
    #         functional_module_arrange_position = %s,
    #         module_arrange_position = %s,
    #         pathway_arrange_position = %s,
    #         molecule_arrange_position = %s,
    #         functional_module_position_limits = %s,
    #         module_position_limits = %s,
    #         pathway_position_limits = %s,
    #         molecule_position_limits = %s)
    #         ',
    #         input$relationship_network_include_functional_modules,
    #         input$relationship_network_include_modules,
    #         input$relationship_network_include_pathways,
    #         input$relationship_network_include_molecules,
    #         input$relationship_network_functional_module_text,
    #         input$relationship_network_module_text,
    #         input$relationship_network_pathway_text,
    #         input$relationship_network_molecule_text,
    #         input$relationship_network_circular_plot,
    #         paste0(
    #           '"',
    #           input$relationship_network_functional_module_color,
    #           '"'
    #         ),
    #         paste0('"', input$relationship_network_module_color, '"'),
    #         paste0('"', input$relationship_network_pathway_color, '"'),
    #         paste0('"', input$relationship_network_molecule_color, '"'),
    #         input$relationship_network_functional_module_arrange_position,
    #         input$relationship_network_module_arrange_position,
    #         input$relationship_network_pathway_arrange_position,
    #         input$relationship_network_molecule_arrange_position,
    #         functional_module_position_limits,
    #         module_position_limits,
    #         pathway_position_limits,
    #         molecule_position_limits
    #       )
    #
    #     relationship_network_code <-
    #       paste0(relationship_network_code1,
    #              relationship_network_code2,
    #              sep = "\n")
    #
    #     relationship_network_code(relationship_network_code)
    #
    #   }
    # })
    #
    # output$relationship_network <-
    #   renderPlot({
    #     req(relationship_network())
    #     relationship_network()
    #   },
    #   res = 96)
    #
    # ######code for relationship_network
    # ####show code
    # observeEvent(input$show_relationship_network_code, {
    #   if (is.null(relationship_network_code()) ||
    #       length(relationship_network_code()) == 0) {
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
    #       relationship_network_code()
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
    # output$download_relationship_network <-
    #   downloadHandler(
    #     filename = function() {
    #       paste0(
    #         "relationship_network_",
    #         input$relationship_network_level,
    #         ".",
    #         input$relationship_network_type
    #       )
    #     },
    #     content = function(file) {
    #       ggsave(
    #         file,
    #         plot = relationship_network(),
    #         width = input$relationship_network_width,
    #         height = input$relationship_network_height
    #       )
    #     }
    #   )
    #
    # observe({
    #   if (is.null(relationship_network()) ||
    #       length(relationship_network()) == 0) {
    #     shinyjs::disable("download_relationship_network")
    #   } else {
    #     shinyjs::enable("download_relationship_network")
    #   }
    # })
    #
    #
    # ###Go to llm intepretation tab
    # ####if there is not enriched_functional_module, show a warning message
    # observeEvent(input$go2llm_interpretation_1, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched_functional_module available",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     updateTabItems(session = session,
    #                    inputId = "tabs",
    #                    selected = "llm_interpretation")
    #   }
    # })
    #
    # observeEvent(input$go2llm_interpretation_2, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched_functional_module available",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     updateTabItems(session = session,
    #                    inputId = "tabs",
    #                    selected = "llm_interpretation")
    #   }
    # })
    #
    # observeEvent(input$go2llm_interpretation_3, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched_functional_module available",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     updateTabItems(session = session,
    #                    inputId = "tabs",
    #                    selected = "llm_interpretation")
    #   }
    # })
    #
    # observeEvent(input$go2llm_interpretation_4, {
    #   # Check if enriched_functional_module is available
    #   if (is.null(enriched_functional_module()) ||
    #       length(enriched_functional_module()) == 0) {
    #     showModal(
    #       modalDialog(
    #         title = "Warning",
    #         "No enriched_functional_module available",
    #         easyClose = TRUE,
    #         footer = modalButton("Close")
    #       )
    #     )
    #   } else {
    #     updateTabItems(session = session,
    #                    inputId = "tabs",
    #                    selected = "llm_interpretation")
    #   }
    # })
    #
    # ###Step 7 LLM interpretation ====
    # # Define enriched_functional_module as a reactive value
    # llm_interpretation_result <-
    #   reactiveVal("")
    #
    # llm_interpretation_code <-
    #   reactiveVal()
    #
    # openai_key <-
    #   reactiveVal()
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
    #         llm_interpretation_result <-
    #         interpret_pathways(
    #         object = enriched_functional_module,
    #         p.adjust.cutoff = %s,
    #         disease = %s,
    #         count.cutoff = %s,
    #         top_n = %s)
    #         ',
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
    # output$download_llm_interpretation_result <-
    #   shiny::downloadHandler(
    #     filename = function() {
    #       "llm_interpretation_result.md"
    #     },
    #     content = function(file) {
    #       writeLines(llm_interpretation_result(), file)
    #     }
    #   )
    #
    # observe({
    #   if (is.null(llm_interpretation_result()) ||
    #       length(llm_interpretation_result()) == 0) {
    #     shinyjs::disable("download_llm_interpretation_result")
    #   } else {
    #     shinyjs::enable("download_llm_interpretation_result")
    #   }
    # })
    #
    #
    # output$llm_enriched_functional_modules1 <-
    #   DT::renderDT({
    #     data <- tryCatch(
    #       enriched_functional_module()@merged_module$functional_module_result,
    #       error = function(e) NULL)
    #
    #     req(data)
    #
    #     DT::datatable(
    #       data,
    #       options = list(pageLength = 10,
    #                      scrollX = TRUE)
    #     )
    #   })
    #
    # output$llm_enriched_functional_modules2 <-
    #   DT::renderDT({
    #     data <- tryCatch(
    #       enriched_functional_module()@merged_module$result_with_module,
    #       error = function(e) NULL
    #     )
    #
    #     req(data)
    #
    #     DT::datatable(
    #       data,
    #       options = list(pageLength = 10,
    #                      scrollX = TRUE)
    #     )
    #   })
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
    # ####show code
    # observeEvent(input$show_llm_interpretation_code, {
    #   if (is.null(llm_interpretation_code()) ||
    #       length(llm_interpretation_code()) == 0) {
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
    #       llm_interpretation_code()
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
    # ###Go to result tab
    # ####if there is not enriched_functional_module, show a warning message
    # observeEvent(input$go2results, {
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
    #                    selected = "results")
    #   }
    # })
    #
    # ###Step 8 Result and report ====
    # report_code <-
    #   reactiveVal()
    # report_path <-
    #   reactiveVal()
    #
    # observeEvent(input$generate_report, {
    #   # Check if enriched_functional_module and llm_interpretation_result are
    #   #  available
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
    #     if (is.null(llm_interpretation_result()) ||
    #         length(llm_interpretation_result()) == 0) {
    #       showModal(
    #         modalDialog(
    #           title = "Warning",
    #           "No LLM interpreatation result available.",
    #           easyClose = TRUE,
    #           footer = modalButton("Close")
    #         )
    #       )
    #     } else{
    #       # shinyjs::show("loading")
    #
    #       withProgress(message = 'Analysis in progress...', {
    #         tryCatch({
    #           report_path <-
    #             file.path("files",
    #                       paste(sample(
    #                         c(letters, LETTERS, 0:9),
    #                         30, replace = TRUE
    #                       ), collapse = ""))
    #
    #           report_functional_module(
    #             object = enriched_functional_module(),
    #             interpretation_result = llm_interpretation_result(),
    #             path = report_path,
    #             type = "html"
    #           )
    #         },
    #         error = function(e) {
    #           showModal(
    #             modalDialog(
    #               title = "Error",
    #               paste("Details:", e$message),
    #               easyClose = TRUE,
    #               footer = modalButton("Close")
    #             )
    #           )
    #         })
    #       })
    #
    #       report_path(report_path)
    #
    #       # shinyjs::hide("loading")
    #
    #       ##save code
    #       report_code <-
    #         sprintf(
    #           '
    #           report_functional_module(
    #           object = enriched_functional_module,
    #           interpretation_result = interpretation_result,
    #           path = %s,
    #           type = "html")
    #         ',
    #         paste0('"', report_path(), '"')
    #         )
    #       report_code(report_code)
    #     }
    #   }
    # })
    #
    # # Update UI to display HTML
    # output$mapa_report <-
    #   renderUI({
    #     tryCatch(
    #       includeHTML(file.path(
    #         report_path(), "Report/mapa_report.html"
    #       )),
    #       error = function(e) {
    #         NULL
    #       }
    #     )
    #   })
    #
    # ###donload the zip report file
    # output$download_report <-
    #   shiny::downloadHandler(
    #     filename = function() {
    #       "Report.zip"
    #     },
    #     content = function(file) {
    #       zip_path <-
    #         paste0(report_path(), "/Report.zip")
    #       zip(zipfile = zip_path,
    #           files = paste0(report_path(),
    #                          "/Report"))
    #       file.copy(zip_path, file)
    #     }
    #   )
    #
    # observe({
    #   if (is.null(report_path()) ||
    #       length(report_path()) == 0) {
    #     shinyjs::disable("download_report")
    #   } else {
    #     shinyjs::enable("download_report")
    #   }
    # })
    #
    # ###To delete the zip file and folder when the user closes the app
    # session$onSessionEnded(function() {
    #   all_files_folders <-
    #     list.files("files", full.names = TRUE)
    #
    #   folders <-
    #     Filter(function(x) {
    #       file.info(x)$isdir
    #     }, all_files_folders)
    #
    #   regex_pattern <- "^[A-Za-z0-9]{30}$"
    #
    #   report_dirs <-
    #     Filter(function(folder) {
    #       folder_name <- basename(folder)
    #       grepl(regex_pattern, folder_name)
    #     }, folders)
    #
    #   unlink(report_dirs, recursive = TRUE)
    # })
    #
    # ####show code
    # observeEvent(input$show_report_code, {
    #   if (is.null(report_code()) ||
    #       length(report_code()) == 0) {
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
    #       report_code()
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
  }
