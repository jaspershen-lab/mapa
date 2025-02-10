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
    enriched_modules_result <- merge_pathways_server("merge_pathways_tab", enriched_pathways = enriched_pathways_result, tab_switch)


    ### Step 4 merge modules ----
    enriched_functional_module_result <- merge_modules_server("merge_modules_tab", enriched_modules = enriched_modules_result, tab_switch)

    ### Step 5 Translation ----

    ### Step 6 Data visualization ----
    data_visualization_server("data_visualization_tab", enriched_functional_module = enriched_functional_module_result, tab_switch)

    ### Step 7 LLM interpretation ----
    llm_interpretation_result <- llm_interpretation_server("llm_interpretation_tab", enriched_functional_module = enriched_functional_module_result, tab_switch)

    ### Step 8 Result and report ----
    results_server("results_tab", enriched_functional_module = enriched_functional_module_result, llm_interpretation_result = llm_interpretation_result, tab_switch)
  }
