server <-
  function(input, output, session) {
    ###Step 1: Upload data
    ###read the uploaded data
    variable_info <-
      reactive({
        # When 'Use example' is checked
        if (input$use_example) {
          return(read.csv("files/example.csv")) # replace with the path to your built-in example CSV file
        }
        # When a file is uploaded
        in_file <- input$variable_info
        if (is.null(in_file)) {
          return(NULL)
        }
        # Read the uploaded file based on its type
        if (grepl("\\.csv$", in_file$name)) {
          return(read.csv(in_file$datapath))
        } else if (grepl("\\.(xlsx|xls)$", in_file$name)) {
          return(readxl::read_excel(in_file$datapath))
        }
      })

    output$variable_info <-
      shiny::renderDataTable({
        req(variable_info())
        variable_info()
      },
      options = list(pageLength = 10))

    #####define variable_info_new
    variable_info_new <-
      reactiveVal()

    ###if the user don't upload data and don't use example data,
    ####show a warning message
    observeEvent(input$map_id, {
      if (is.null(input$variable_info) && !input$use_example) {
        showModal(
          modalDialog(
            title = "Warning",
            "No data is uploaded",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }
    })

    ####map the IDs
    upload_data_code <-
      reactiveVal()
    observeEvent(input$map_id, {
      req(variable_info())
      variable_info_old <-
        variable_info()
      if (input$id_type == "ensembl") {
        colnames(variable_info_old) <- c("ensembl")
        library(clusterProfiler)
        library(org.Hs.eg.db)
        other_id <-
          tryCatch(
            clusterProfiler::bitr(
              variable_info_old$ensembl,
              fromType = "ENSEMBL",
              toType = c("UNIPROT", "ENTREZID", "SYMBOL"),
              OrgDb = org.Hs.eg.db
            ),
            error = function(e) {
              data.frame(
                ENSEMBL = variable_info_old$ensembl,
                UNIPROT = NA,
                ENTREZID = NA
              )
            }
          )
        ### remove duplicated rows
        other_id <-
          dplyr::distinct(other_id, ENSEMBL, .keep_all = TRUE)
        variable_info_old <-
          dplyr::left_join(variable_info_old,
                           other_id, by = c("ensembl" = "ENSEMBL"))
        colnames(variable_info_old) <-
          c("ensembl", "uniprot", "entrezid", "symbol")

        ####save code
        upload_data_code <-
          sprintf(
            'library(clusterProfiler)
            library(org.Hs.eg.db)
            other_id <-
            clusterProfiler::bitr(
            variable_info$ensembl,
            fromType = "ENSEMBL",
            toType = c("UNIPROT", "ENTREZID", "SYMBOL"),
            OrgDb = org.Hs.eg.db)
            other_id <-
            dplyr::distinct(other_id, ENSEMBL, .keep_all = TRUE)
            variable_info <-
            dplyr::left_join(variable_info,
            other_id, by = c("ensembl" = "ENSEMBL"))
            colnames(variable_info) <-
            c("ensembl", "uniprot", "entrezid", "symbol")
            '
          )
        # Save to a text file
        # writeLines(upload_data_code, "files/upload_data_code.txt")
      }
      if (input$id_type == "uniprot") {
        colnames(variable_info_old) <- c("uniprot")
        library(clusterProfiler)
        library(org.Hs.eg.db)
        other_id <-
          tryCatch(
            clusterProfiler::bitr(
              variable_info_old$uniprot,
              fromType = "UNIPROT",
              toType = c("ENSEMBL", "ENTREZID", "SYMBOL"),
              OrgDb = org.Hs.eg.db
            ),
            error = function(e) {
              data.frame(
                UNIPROT = variable_info_old$uniprot,
                ENSEMBL = NA,
                ENTREZID = NA
              )
            }
          )
        # remove duplicated rows
        other_id <-
          dplyr::distinct(other_id, UNIPROT, .keep_all = TRUE)
        variable_info_old <-
          dplyr::left_join(variable_info_old,
                           other_id, by = c("uniprot" = "UNIPROT"))
        colnames(variable_info_old) <-
          c("uniprot", "ensembl", "entrezid", "symbol")

        ####save code
        upload_data_code <-
          sprintf(
            'library(clusterProfiler)
            library(org.Hs.eg.db)
            other_id <-
            clusterProfiler::bitr(
            variable_info$uniprot,
            fromType = "UNIPROT",
            toType = c("ENSEMBL", "ENTREZID", "SYMBOL"),
            OrgDb = org.Hs.eg.db)
            other_id <-
            dplyr::distinct(other_id, UNIPROT, .keep_all = TRUE)
            variable_info <-
            dplyr::left_join(variable_info,
            other_id, by = c("uniprot" = "UNIPROT"))
            colnames(variable_info) <-
            c("uniprot", "ensembl", "entrezid", "symbol")
            '
          )
        # Save to a text file
        # writeLines(upload_data_code, "files/upload_data_code.txt")
      }
      if (input$id_type == "entrezid") {
        colnames(variable_info_old) <- c("entrezid")
        library(clusterProfiler)
        library(org.Hs.eg.db)
        other_id <-
          tryCatch(
            clusterProfiler::bitr(
              variable_info_old$entrezid,
              fromType = "ENTREZID",
              toType = c("ENSEMBL", "UNIPROT", "SYMBOL"),
              OrgDb = org.Hs.eg.db
            ),
            error = function(e) {
              data.frame(
                ENTREZID = variable_info_old$entrezid,
                ENSEMBL = NA,
                UNIPROT = NA
              )
            }
          )
        # remove duplicated rows
        other_id <-
          dplyr::distinct(other_id, ENTREZID, .keep_all = TRUE)
        variable_info_old <-
          dplyr::left_join(variable_info_old,
                           other_id,
                           by = c("entrezid" = "ENTREZID"))
        colnames(variable_info_old) <-
          c("entrezid", "ensembl", "entrezid", "symbol")

        ####save code
        upload_data_code <-
          sprintf(
            'library(clusterProfiler)
            library(org.Hs.eg.db)
            other_id <-
            clusterProfiler::bitr(
            variable_info$uniprot,
            fromType = "ENTREZID",
            toType = c("ENSEMBL", "UNIPROT", "SYMBOL"),
            OrgDb = org.Hs.eg.db)
            other_id <-
            dplyr::distinct(other_id, ENTREZID, .keep_all = TRUE)
            variable_info <-
            dplyr::left_join(variable_info,
            other_id, by = c("entrezid" = "ENTREZID"))
            colnames(variable_info) <-
            c("entrezid", "ensembl", "entrezid", "symbol")
            '
          )
        # Save to a text file
        # writeLines(upload_data_code, "files/upload_data_code.txt")
      }
      upload_data_code(upload_data_code)
      variable_info_new(variable_info_old)
    })

    output$variable_info <-
      shiny::renderDataTable({
        req(variable_info_new())
        variable_info_new()
      },
      options = list(pageLength = 10))

    ###download the variable_info------------------------------------------
    output$download_variable_info <-
      shiny::downloadHandler(
        filename = function() {
          "variable_info.csv"
        },
        content = function(file) {
          write.csv(variable_info_new(),
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(variable_info_new()) ||
          length(variable_info_new()) == 0) {
        shinyjs::disable("download_variable_info")
      } else {
        shinyjs::enable("download_variable_info")
      }
    })


    ####show code
    observeEvent(input$show_upload_data_code, {
      if (is.null(upload_data_code()) ||
          length(upload_data_code()) == 0) {
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
          upload_data_code()
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

    ####Click next and then go to enrich pathways tab
    observeEvent(input$go2enrich_pathways, {
      # Check if variable_info_new is available
      if (is.null(variable_info_new()) ||
          length(variable_info_new()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please upload data and click 'Submit' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        # Navigate to the enrich_pathways tab
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "enrich_pathways")
      }
    })




    ###--------------------------------------------------------------------
    ###step 2 enrich pathways
    ###when the user click submit_enrich_pathways, begin enrich pathways
    # Define enriched_pathways as a reactive value
    enriched_pathways <-
      reactiveVal()

    enrich_pathways_code <-
      reactiveVal()

    observeEvent(input$submit_enrich_pathways, {
      # Check if variable_info_new is available
      if (is.null(variable_info_new()) ||
          length(variable_info_new()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No data available. Please 'Upload data' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ReactomePA)
        if (length(input$pathway_database) == 0) {
          showModal(
            modalDialog(
              title = "Warning",
              "Please select at least one pathway database.",
              easyClose = TRUE,
              footer = modalButton("Close")
            )
          )
        } else{
          shinyjs::show("loading")
          result <-
            enrich_pathway(
              variable_info_new(),
              database = input$pathway_database,
              save_to_local = FALSE,
              path = "result",
              OrgDb = org.Hs.eg.db,
              organism = input$organism,
              keyType = "ENTREZID",
              use_internal_data = FALSE,
              ont = "ALL",
              pvalueCutoff = input$p_value_cutoff,
              pAdjustMethod = input$p_adjust_method,
              qvalueCutoff = 0.2,
              minGSSize = input$gene_set_size[1],
              maxGSSize = input$gene_set_size[2],
              readable = FALSE,
              pool = FALSE
            )

          # enriched_pathways <-
          enriched_pathways(result)

          shinyjs::hide("loading")

          ###save code
          pathway_database <-
            paste0("c(", paste(unlist(
              lapply(paste(input$pathway_database), function(x)
                paste0('"', x, '"'))
            ),
            collapse = ", "), ")")

          enrich_pathways_code <-
            sprintf(
              '
              enriched_pathways <-
              enrich_pathway(
              variable_info,
              database = %s,
              OrgDb = org.Hs.eg.db,
              organism = %s,
              pvalueCutoff = %s,
              pAdjustMethod = %s,
              minGSSize = %s,
              maxGSSize = %s)
              ',
              pathway_database,
              input$organism,
              input$p_value_cutoff,
              input$p_adjust_method,
              input$gene_set_size[1],
              input$gene_set_size[2]
            )

          enrich_pathways_code(enrich_pathways_code)

          # # Save to a text file
          # writeLines(enrich_pathways_code,
          #            "files/enrich_pathways_code.txt")

        }
      }
    })


    output$enriched_pathways_go <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_pathways()@enrichment_go_result@result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$enriched_pathways_kegg <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_pathways()@enrichment_kegg_result@result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$enriched_pathways_reactome <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_pathways()@enrichment_reactome_result@result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$download_enriched_pathways_go <-
      shiny::downloadHandler(
        filename = function() {
          "enriched_pathways_go.csv"
        },
        content = function(file) {
          write.csv(enriched_pathways()@enrichment_go_result@result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        shinyjs::disable("download_enriched_pathways_go")
      } else {
        if (length(enriched_pathways()@enrichment_go_result) == 0) {
          shinyjs::disable("download_enriched_pathways_go")
        } else{
          shinyjs::enable("download_enriched_pathways_go")
        }
      }
    })

    output$download_enriched_pathways_kegg <-
      shiny::downloadHandler(
        filename = function() {
          "enriched_pathways_kegg.csv"
        },
        content = function(file) {
          write.csv(enriched_pathways()@enrichment_kegg_result@result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        shinyjs::disable("download_enriched_pathways_kegg")
      } else {
        if (length(enriched_pathways()@enrichment_kegg_result) == 0) {
          shinyjs::disable("download_enriched_pathways_kegg")
        } else{
          shinyjs::enable("download_enriched_pathways_kegg")
        }
      }
    })

    output$download_enriched_pathways_reactome <-
      shiny::downloadHandler(
        filename = function() {
          "enriched_pathways_reactome.csv"
        },
        content = function(file) {
          write.csv(enriched_pathways()@enrichment_reactome_result@result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        shinyjs::disable("download_enriched_pathways_reactome")
      } else {
        if (length(enriched_pathways()@enrichment_reactome_result) == 0) {
          shinyjs::disable("download_enriched_pathways_reactome")
        } else{
          shinyjs::enable("download_enriched_pathways_reactome")
        }
      }
    })

    ####show code
    observeEvent(input$show_enrich_pathways_code, {
      if (is.null(enrich_pathways_code()) ||
          length(enrich_pathways_code()) == 0) {
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
          enrich_pathways_code()
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

    ####Go to merge pathways tab
    ####if there is not enriched_pathways,
    ####show a warning message
    observeEvent(input$go2merge_pathways, {
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please enrich pathways first",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "merge_pathways")
      }
    })


    ###--------------------------------------------------------------------
    ###step 3 merge pathways
    # Define enriched_modules as a reactive value
    enriched_modules <-
      reactiveVal()

    merge_pathways_code <-
      reactiveVal()

    observeEvent(input$submit_merge_pathways, {
      # Check if variable_info_new is available
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched pathways data available. Please 'Enrich pathways' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ReactomePA)
        shinyjs::show("loading")
        result <-
          merge_pathways(
            object = enriched_pathways(),
            p.adjust.cutoff.go = input$p.adjust.cutoff.go,
            p.adjust.cutoff.kegg = input$p.adjust.cutoff.kegg,
            p.adjust.cutoff.reactome = input$p.adjust.cutoff.reactome,
            count.cutoff.go = input$count.cutoff.go,
            count.cutoff.kegg = input$count.cutoff.kegg,
            count.cutoff.reactome = input$count.cutoff.reactome,
            sim.cutoff.go = input$sim.cutoff.go,
            sim.cutoff.kegg = input$sim.cutoff.kegg,
            sim.cutoff.reactome = input$sim.cutoff.reactome,
            measure.method.go = input$measure.method.go,
            measure.method.kegg = input$measure.method.kegg,
            measure.method.reactome = input$measure.method.reactome,
            path = "result",
            save_to_local = FALSE
          )

        # enriched_modules <-
        enriched_modules(result)

        shinyjs::hide("loading")

        ###save code
        merge_pathways_code <-
          sprintf(
            '
    enriched_modules <-
    merge_pathways(
    object = enriched_pathways,
    p.adjust.cutoff.go = %s,
    p.adjust.cutoff.kegg = %s,
    p.adjust.cutoff.reactome = %s,
    count.cutoff.go = %s,
    count.cutoff.kegg = %s,
    count.cutoff.reactome = %s,
    sim.cutoff.go = %s,
    sim.cutoff.kegg = %s,
    sim.cutoff.reactome = %s,
    measure.method.go = %s,
    measure.method.kegg = %s,
    measure.method.reactome = %s
  )
    ',
  input$p.adjust.cutoff.go,
  input$p.adjust.cutoff.kegg,
  input$p.adjust.cutoff.reactome,
  input$count.cutoff.go,
  input$count.cutoff.kegg,
  input$count.cutoff.reactome,
  input$sim.cutoff.go,
  input$sim.cutoff.kegg,
  input$sim.cutoff.reactome,
  input$measure.method.go,
  input$measure.method.kegg,
  input$measure.method.reactome
          )
        # Save to a text file
        # writeLines(merge_pathways_code, "files/merge_pathways_code.txt")
        merge_pathways_code(merge_pathways_code)
      }
    })

    output$merged_pathway_go <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_modules()@merged_pathway_go$module_result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$merged_pathway_kegg <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_modules()@merged_pathway_kegg$module_result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$merged_pathway_reactome <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_modules()@merged_pathway_reactome$module_result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$download_merged_pathway_go <-
      shiny::downloadHandler(
        filename = function() {
          "merged_pathway_go.csv"
        },
        content = function(file) {
          write.csv(enriched_modules()@merged_pathway_go$module_result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        shinyjs::disable("download_merged_pathway_go")
      } else {
        if (length(enriched_modules()@merged_pathway_go) == 0) {
          shinyjs::disable("download_merged_pathway_go")
        } else{
          shinyjs::enable("download_merged_pathway_go")
        }
      }
    })

    output$download_merged_pathway_kegg <-
      shiny::downloadHandler(
        filename = function() {
          "merged_pathway_kegg.csv"
        },
        content = function(file) {
          write.csv(enriched_modules()@merged_pathway_kegg$module_result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        shinyjs::disable("download_merged_pathway_kegg")
      } else {
        if (length(enriched_modules()@merged_pathway_kegg) == 0) {
          shinyjs::disable("download_merged_pathway_kegg")
        } else{
          shinyjs::enable("download_merged_pathway_kegg")
        }
      }
    })

    output$download_merged_pathway_reactome <-
      shiny::downloadHandler(
        filename = function() {
          "merged_pathway_reactome.csv"
        },
        content = function(file) {
          write.csv(enriched_modules()@merged_pathway_reactome$module_result,
                    file,
                    row.names = FALSE)
        }
      )

    observe({
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        shinyjs::disable("download_merged_pathway_reactome")
      } else {
        if (length(enriched_modules()@merged_pathway_reactome) == 0) {
          shinyjs::disable("download_merged_pathway_reactome")
        } else{
          shinyjs::enable("download_merged_pathway_reactome")
        }
      }
    })


    ######data visualization
    # GO Plot generation logic
    enirched_module_go_plot <-
      reactiveVal()
    observeEvent(input$generate_enirched_module_plot_go, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_similarity_network(
            object = enriched_modules(),
            level = "module",
            database = "go",
            degree_cutoff = input$enirched_module_plot_degree_cutoff_go,
            text = input$enirched_module_plot_text_go,
            text_all = input$enirched_module_plot_text_all_go
          )
        enirched_module_go_plot(plot)

        shinyjs::hide("loading")
      }
    })

    output$enirched_module_go_plot <-
      shiny::renderPlot({
        req(tryCatch(
          enirched_module_go_plot(),
          error = function(e)
            NULL
        ))
      })

    # kegg Plot generation logic
    enirched_module_kegg_plot <-
      reactiveVal()
    observeEvent(input$generate_enirched_module_plot_kegg, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_similarity_network(
            object = enriched_modules(),
            level = "module",
            database = "kegg",
            degree_cutoff = input$enirched_module_plot_degree_cutoff_kegg,
            text = input$enirched_module_plot_text_kegg,
            text_all = input$enirched_module_plot_text_all_kegg
          )

        enirched_module_kegg_plot(plot)

        shinyjs::hide("loading")
      }
    })

    output$enirched_module_kegg_plot <-
      shiny::renderPlot({
        req(tryCatch(
          enirched_module_kegg_plot(),
          error = function(e)
            NULL
        ))
      })

    # reactome Plot generation logic
    enirched_module_reactome_plot <-
      reactiveVal()
    observeEvent(input$generate_enirched_module_plot_reactome, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_similarity_network(
            object = enriched_modules(),
            level = "module",
            database = "reactome",
            degree_cutoff = input$enirched_module_plot_degree_cutoff_reactome,
            text = input$enirched_module_plot_text_reactome,
            text_all = input$enirched_module_plot_text_all_reactome
          )

        enirched_module_reactome_plot(plot)
        shinyjs::hide("loading")
      }
    })

    output$enirched_module_reactome_plot <-
      shiny::renderPlot({
        req(tryCatch(
          enirched_module_reactome_plot(),
          error = function(e)
            NULL
        ))
      })

    ####show code
    observeEvent(input$show_merge_pathways_code, {
      if (is.null(merge_pathways_code()) ||
          length(merge_pathways_code()) == 0) {
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
          merge_pathways_code()
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

    ###Go to merge modules tab
    ####if there is not enriched_modules, show a warning message
    observeEvent(input$go2merge_modules, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please merge pathways first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "merge_modules")
      }
    })

    ###--------------------------------------------------------------------
    ###Step 4 merge modules
    # Define enriched_functional_module as a reactive value
    enriched_functional_module <-
      reactiveVal()

    merge_modules_code <-
      reactiveVal()

    observeEvent(input$submit_merge_modules, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Enrich modules' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ReactomePA)
        shinyjs::show("loading")
        # Perform analysis with user-provided parameters
        result <-
          merge_modules(
            object = enriched_modules(),
            sim.cutoff = input$sim.cutoff.module,
            measure_method = input$measure.method.module,
            path = "result",
            save_to_local = FALSE
          )

        # enriched_functional_module <-
        enriched_functional_module(result)
        shinyjs::hide("loading")

        ##save code
        merge_modules_code <-
          sprintf(
            '
            enriched_functional_module <-
            merge_modules(
            object = enriched_modules,
            sim.cutoff = %s,
            measure_method = %s)
            ',
            input$sim.cutoff.module,
            input$measure.method.module
          )

        merge_modules_code(merge_modules_code)
      }
    })

    output$enriched_functional_modules <-
      shiny::renderDataTable({
        req(tryCatch(
          enriched_functional_module()@merged_module$functional_module_result,
          error = function(e)
            NULL
        ))
      },
      options = list(pageLength = 10))

    output$download_enriched_functional_modules <-
      shiny::downloadHandler(
        filename = function() {
          "merged_modules.csv"
        },
        content = function(file) {
          write.csv(
            enriched_functional_module()@merged_module$functional_module_result,
            file,
            row.names = FALSE
          )
        }
      )

    observe({
      if (is.null(enriched_functional_module()) ||
          length(enriched_functional_module()) == 0) {
        shinyjs::disable("download_enriched_functional_modules")
      } else {
        if (length(enriched_functional_module()@merged_module) == 0) {
          shinyjs::disable("download_enriched_functional_modules")
        } else{
          shinyjs::enable("download_enriched_functional_modules")
        }
      }
    })

    ######data visualization
    ###define enirched_functional_module_plot
    enirched_functional_module_plot <-
      reactiveVal()
    observeEvent(input$generate_enirched_functional_module, {
      # Check if enriched_functional_module is available
      if (is.null(enriched_functional_module()) ||
          length(enriched_functional_module()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched functional modules data available. Please 'Merge modules' first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_similarity_network(
            object = enriched_functional_module(),
            level = "functional_module",
            degree_cutoff = input$enirched_functional_moduleplot__degree_cutoff,
            text = input$enirched_functional_module_plot_text,
            text_all = input$enirched_functional_module_plot_text_all
          )

        enirched_functional_module_plot(plot)
        shinyjs::hide("loading")
      }
    })

    output$enirched_functional_module_plot <-
      shiny::renderPlot({
        req(tryCatch(
          enirched_functional_module_plot(),
          error = function(e)
            NULL
        ))
      })

    ####show code
    observeEvent(input$show_merge_modules_code, {
      if (is.null(merge_modules_code()) ||
          length(merge_modules_code()) == 0) {
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
          merge_modules_code()
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


    ###Go to data visualization tab
    ####if there is not enriched_functional_module, show a warning message
    observeEvent(input$go2data_visualization, {
      # Check if enriched_functional_module is available
      if (is.null(enriched_functional_module()) ||
          length(enriched_functional_module()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please merge modules first.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "data_visualization")
      }
    })

    ###--------------------------------------------------------------------
    ###Step 5 Data visualization
    # Observe file upload
    # enriched_functional_module <-
    #   reactiveVal()
    observeEvent(input$upload_enriched_functional_module, {
      if (!is.null(input$upload_enriched_functional_module$datapath)) {
        tempEnv <- new.env()
        load(input$upload_enriched_functional_module$datapath,
             envir = tempEnv)

        names <- ls(tempEnv)
        # Handle user response
        enriched_functional_module(get(names[1], envir = tempEnv))
      } else {
        showModal(
          modalDialog(
            title = "Error",
            "The uploaded file should contain exactly one object.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      }
    })

    #####barplot
    # Observe generate barplot button click
    barplot <-
      reactiveVal()
    barplot_code <-
      reactiveVal()

    observeEvent(input$generate_barplot, {
      if (is.null(enriched_functional_module())) {
        # No enriched functional module available
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched functional module data available. Please complete the previous steps or upload the data.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_pathway_bar(
            object = enriched_functional_module(),
            top_n = input$barplot_top_n,
            y_lable_width = input$barplot_y_lable_width,
            p.adjust.cutoff = input$barplot_p_adjust_cutoff,
            count.cutoff = input$barplot_count_cutoff,
            level = input$barplot_level,
            database = input$barplot_database,
            line_type = input$line_type,
            database_color = c(
              GO = input$barplot_go_color,
              KEGG = input$barplot_kegg_color,
              Reactome = input$barplot_reactome_color
            )
          )

        barplot(plot)

        ###save code
        data_color <-
          paste0("c(", paste(paste(
            c("GO", "KEGG", "Reactome"),
            c(
              paste0('"', input$barplot_go_color, '"'),
              paste0('"', input$barplot_kegg_color, '"'),
              paste0('"', input$barplot_reactome_color, '"')
            ),
            sep = " = "
          ),
          collapse = ", "), ")")

        barplot_database <-
          paste0("c(", paste(unlist(lapply(paste(input$barplot_database), function(x)
            paste0('"', x, '"'))),
            collapse = ", "), ")")

        barplot_code <-
          sprintf(
            "plot_pathway_bar(
                  object = enriched_functional_module,
                  top_n = %s,
                  y_lable_width = %s,
                  p.adjust.cutoff = %s,
                  count.cutoff = %s,
                  level = %s,
                  database = %s,
                  line_type = %s,
                  database_color = %s
                  )",
            input$barplot_top_n,
            input$barplot_y_lable_width,
            input$barplot_p_adjust_cutoff,
            input$barplot_count_cutoff,
            input$barplot_level,
            barplot_database,
            input$line_type,
            data_color
          )

        barplot_code(barplot_code)
        shinyjs::hide("loading")
      }
    })

    # output$barplot <-
    #   shiny::renderPlot({
    #     if (is.null(barplot()) ||
    #         length(barplot()) == 0) {
    #       NULL
    #     } else{
    #       barplot()
    #     }
    #   })

    output$barplot <-
      renderPlot({
        req(barplot())  # Request the reactive value
        barplot()       # Display the plot
      })

    output$download_barplot <-
      downloadHandler(
        filename = function() {
          paste0("pathway_barplot.", input$barplot_type)
        },
        content = function(file) {
          ggsave(
            file,
            plot = barplot(),
            width = input$barplot_width,
            height = input$barplot_height
          )
        }
      )

    observe({
      if (is.null(barplot()) ||
          length(barplot()) == 0) {
        shinyjs::disable("download_barplot")
      } else {
        shinyjs::enable("download_barplot")
      }
    })

    ######code for barplot
    ####show code
    observeEvent(input$show_barplot_code, {
      if (is.null(barplot_code()) ||
          length(barplot_code()) == 0) {
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
          barplot_code()
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







    #####module_similarity_network
    # Observe generate module_similarity_network button click
    module_similarity_network <-
      reactiveVal()
    module_similarity_network_code <-
      reactiveVal()

    observeEvent(input$generate_module_similarity_network, {
      if (is.null(enriched_functional_module())) {
        # No enriched functional module available
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched functional module data available. Please complete the previous steps or upload the data.",
            easyClose = TRUE,
            footer = modalButton("Close")
          )
        )
      } else {
        shinyjs::show("loading")
        plot <-
          plot_similarity_network(
            object = enriched_functional_module(),
            level = "module",
            database = input$module_similarity_network_database,
            degree_cutoff = input$module_similarity_network_degree_cutoff,
            text = input$module_similarity_network_text,
            text_all = input$module_similarity_network_text_all
          )

        module_similarity_network(plot)

        ###save code
        module_similarity_network_code <-
          sprintf(
            '
            plot_similarity_network(
            object = enriched_functional_module,
            level = "module",
            database = %s,
            degree_cutoff = %s,
            text = %s,
            text_all = %s
          )',
          input$module_similarity_network_database,
          input$module_similarity_network_degree_cutoff,
          input$module_similarity_network_text,
          input$module_similarity_network_text_all
          )

        module_similarity_network_code(module_similarity_network_code)
        shinyjs::hide("loading")
      }
    })

    output$module_similarity_network <-
      renderPlot({
        req(module_similarity_network())
        module_similarity_network()
      })

    output$download_module_similarity_network <-
      downloadHandler(
        filename = function() {
          paste0("module_similarity_network_",
                 input$module_similarity_network_database,
                 ".",
                 input$module_similarity_network_type)
        },
        content = function(file) {
          ggsave(
            file,
            plot = barplot(),
            width = input$module_similarity_network_width,
            height = input$module_similarity_network_height
          )
        }
      )

    observe({
      if (is.null(module_similarity_network()) ||
          length(module_similarity_network()) == 0) {
        shinyjs::disable("download_module_similarity_network")
      } else {
        shinyjs::enable("download_module_similarity_network")
      }
    })

    ######code for module_similarity_network
    ####show code
    observeEvent(input$show_module_similarity_network_code, {
      if (is.null(module_similarity_network_code()) ||
          length(module_similarity_network_code()) == 0) {
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
          module_similarity_network_code()
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

  }
