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
    processed_info <- reactiveValues(
      variable_info = NULL,
      query_type = NULL,
      organism = NULL
    )
    upload_data_server("upload_data_tab",
                       processed_info = processed_info,
                       tab_switch)

    ### Step 2 enrich pathways ----
    enriched_pathways <- reactiveValues(
      enriched_pathways_res = NULL,
      query_type = NULL,
      organism = NULL,
      available_db = NULL
    )
    enrich_pathway_server("enrich_pathway_tab",
                          processed_info = processed_info,
                          enriched_pathways = enriched_pathways,
                          tab_switch)

    ### Step 3a merge pathways ----
    enriched_modules <- reactiveVal(NULL)
    merge_pathways_server("merge_pathways_tab",
                          enriched_pathways = enriched_pathways,
                          enriched_modules = enriched_modules,
                          tab_switch)

    ### Step 4a merge modules ----
    enriched_functional_module <- reactiveVal(NULL)
    merge_modules_server("merge_modules_tab",
                          enriched_modules = enriched_modules,
                          enriched_functional_module = enriched_functional_module,
                          tab_switch)

    ### Step 3b-4b embed and cluster pathways
    embed_cluster_pathways_server("embed_cluster_pathways_tab",
                                  enriched_pathways = enriched_pathways,
                                  enriched_functional_module = enriched_functional_module,
                                  tab_switch)
    ### Step 5 Translation ----

    ### Step 6 LLM interpretation ----
    llm_interpretation_server("llm_interpretation_tab",
                              enriched_functional_module = enriched_functional_module,
                              tab_switch)

    ### Step 7 Data visualization ----
    data_visualization_server("data_visualization_tab",
                              enriched_functional_module = enriched_functional_module,
                              tab_switch)

    ### Step 8 Result and report ----
    results_server("results_tab",
                   enriched_functional_module = enriched_functional_module_result,
                   llm_interpretation_result = llm_interpretation_result,
                   tab_switch)
  }
