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

    output$contents <-
      shiny::renderDataTable({
        variable_info()
      },
      options = list(pageLength = 10))

    ###if the user don't upload data and don't use example data, show a warning message
    observeEvent(input$map_id, {
      if (is.null(input$variable_info) && !input$use_example) {
        showModal(
          modalDialog(
            title = "Warning",
            "No data is uploaded",
            easyClose = TRUE,
            footer = NULL
          )
        )
      }
    })

    variable_info_new <-
      reactiveVal()
    ####map the IDs
    observeEvent(input$map_id, {
      req(variable_info())
      variable_info_old <- variable_info()
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
      }
      variable_info_new(variable_info_old)

      output$contents <-
        shiny::renderDataTable({
          req(variable_info_new()) # Only render when variable_info has been updated
          variable_info_new()
        },
        options = list(pageLength = 10))

      output$download_variable_info <-
        shiny::downloadHandler(
          filename = function() {
            "variable_info.csv"
          },
          content = function(file) {
            write.csv(variable_info_new(), file, row.names = FALSE)
          }
        )

    })

    # observe({
    #   if (!is.null(variable_info_new()) &&
    #       length(variable_info_new()) > 0) {
    #     shinyjs::show("download_variable_info")
    #   } else {
    #     shinyjs::hide("download_variable_info")
    #   }
    # })


    ######code for this step
    output$data_upload_code <- renderUI({
      if (input$show_code_upload_data %% 2 == 1) {
        code_content <-
          readLines("files/code_upload_data.R")
        code_content <-
          paste(code_content, collapse = "\n")
        tags$pre(code_content)
      }
    })

    ####if there is not variable_info_new, show a warning message
    observeEvent(input$go2enrich_pathways, {
      # Check if variable_info_new is available
      if (is.null(variable_info_new()) ||
          length(variable_info_new()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please upload data and click 'Map ID' first.",
            easyClose = TRUE,
            footer = NULL
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

    observeEvent(input$submit_enrich_pathways, {
      # Check if variable_info_new is available
      if (is.null(variable_info_new()) ||
          length(variable_info_new()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No data available. Please 'Upload data' first.",
            easyClose = TRUE,
            footer = NULL
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
              footer = NULL
            )
          )
        } else{
          shinyjs::show("loading")
          # Run the enrich_pathway analysis
          result <-
            enrich_pathway(
              variable_info_new(),
              database = input$pathway_database,
              # Assuming input$database is the user input
              save_to_local = FALSE,
              path = "result",
              OrgDb = org.Hs.eg.db,
              # Make sure this is defined or passed correctly
              organism = input$organism,
              # User input for organism
              keyType = "ENTREZID",
              use_internal_data = FALSE,
              ont = "ALL",
              pvalueCutoff = input$p_value_cutoff,
              # User input for p-value cutoff
              pAdjustMethod = input$p_adjust_method,
              # User input for p-adjust method
              qvalueCutoff = 0.2,
              minGSSize = input$gene_set_size[1],
              # Assuming this is a range input
              maxGSSize = input$gene_set_size[2],
              readable = FALSE,
              pool = FALSE
            )

          output$enriched_pathways_go <-
            shiny::renderDataTable({
              req(tryCatch(
                result@enrichment_go_result@result,
                error = function(e)
                  NULL
              ))
            },
            options = list(pageLength = 10))

          output$enriched_pathways_kegg <-
            shiny::renderDataTable({
              req(tryCatch(
                result@enrichment_kegg_result@result,
                error = function(e)
                  NULL
              ))
            },
            options = list(pageLength = 10))

          output$enriched_pathways_reactome <-
            shiny::renderDataTable({
              req(tryCatch(
                result@enrichment_reactome_result@result,
                error = function(e)
                  NULL
              ))
            },
            options = list(pageLength = 10))

          enriched_pathways <-
            enriched_pathways(result)

          shinyjs::hide("loading")
        }
      }
    })

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

    ######code for this step
    output$enrich_pathways_code <- renderUI({
      if (input$show_code_enrich_pathways %% 2 == 1) {
        code_content <-
          readLines("files/enrich_pathways_code.R")
        code_content <-
          paste(code_content, collapse = "\n")
        tags$pre(code_content)
      }
    })

    ####if there is not enriched_pathways, show a warning message
    observeEvent(input$go2merge_pathways, {
      # Check if enriched_pathways is available
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "Please enrich pathways first.",
            easyClose = TRUE,
            footer = NULL
          )
        )
      } else {
        # Navigate to the enrich_pathways tab
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

    observeEvent(input$submit_merge_pathways, {
      # Check if variable_info_new is available
      if (is.null(enriched_pathways()) ||
          length(enriched_pathways()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched pathways data available. Please 'Enrich pathways' first.",
            easyClose = TRUE,
            footer = NULL
          )
        )
      } else {
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(ReactomePA)
        shinyjs::show("loading")
        # Perform analysis with user-provided parameters
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
        output$merged_pathway_go <-
          shiny::renderDataTable({
            req(tryCatch(
              result@merged_pathway_go$module_result,
              error = function(e)
                NULL
            ))
          },
          options = list(pageLength = 10))
        output$merged_pathway_kegg <-
          shiny::renderDataTable({
            req(tryCatch(
              result@merged_pathway_kegg$module_result,
              error = function(e)
                NULL
            ))
          },
          options = list(pageLength = 10))
        output$merged_pathway_reactome <-
          shiny::renderDataTable({
            req(tryCatch(
              result@merged_pathway_reactome$module_result,
              error = function(e)
                NULL
            ))
          },
          options = list(pageLength = 10))

        enriched_modules <-
          enriched_modules(result)

        shinyjs::hide("loading")

      }
    })


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


    ######code for this step
    output$merge_pathways_code <- renderUI({
      if (input$show_code_merge_pathways %% 2 == 1) {
        code_content <-
          readLines("files/merge_pathways_code.R")
        code_content <-
          paste(code_content, collapse = "\n")
        tags$pre(code_content)
      }
    })


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
            footer = NULL
          )
        )
      } else {
        # Navigate to the merge_modules tab
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "merge_modules")
      }
    })


    ######data visualization
    # GO Plot generation logic
    observeEvent(input$generate_enirched_module_plot_go, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = NULL
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

        output$enirched_module_go_plot <-
          shiny::renderPlot({
            req(tryCatch(
              plot,
              error = function(e)
                NULL
            ))
          })

        shinyjs::hide("loading")

      }
    })


    # kegg Plot generation logic
    observeEvent(input$generate_enirched_module_plot_kegg, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = NULL
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

        output$enirched_module_kegg_plot <-
          shiny::renderPlot({
            req(tryCatch(
              plot,
              error = function(e)
                NULL
            ))
          })

        shinyjs::hide("loading")

      }
    })


    # reactome Plot generation logic
    observeEvent(input$generate_enirched_module_plot_reactome, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Merge pathways' first.",
            easyClose = TRUE,
            footer = NULL
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

        output$enirched_module_reactome_plot <-
          shiny::renderPlot({
            req(tryCatch(
              plot,
              error = function(e)
                NULL
            ))
          })

        shinyjs::hide("loading")
      }
    })







    ###--------------------------------------------------------------------
    ###step 4 merge modules
    # Define enriched_functional_module as a reactive value
    enriched_functional_module <-
      reactiveVal(NULL)

    observeEvent(input$submit_merge_modules, {
      # Check if enriched_modules is available
      if (is.null(enriched_modules()) ||
          length(enriched_modules()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched modules data available. Please 'Enrich modules' first.",
            easyClose = TRUE,
            footer = NULL
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

        output$enriched_functional_modules <-
          shiny::renderDataTable({
            req(tryCatch(
              result@merged_module$functional_module_result,
              error = function(e)
                NULL
            ))
          },
          options = list(pageLength = 10))

        enriched_functional_module <-
          enriched_functional_module(result)

        shinyjs::hide("loading")

      }
    })

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

    ######code for this step
    output$merge_modules_code <- renderUI({
      if (input$show_code_merge_modules %% 2 == 1) {
        code_content <-
          readLines("files/merge_modules_code.R")
        code_content <-
          paste(code_content, collapse = "\n")
        tags$pre(code_content)
      }
    })


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
            footer = NULL
          )
        )
      } else {
        # Navigate to the merge_modules tab
        updateTabItems(session = session,
                       inputId = "tabs",
                       selected = "data_visualization")
      }
    })

    ######data visualization
    observeEvent(input$generate_enirched_functional_module, {
      # Check if enriched_functional_module is available
      if (is.null(enriched_functional_module()) ||
          length(enriched_functional_module()) == 0) {
        showModal(
          modalDialog(
            title = "Warning",
            "No enriched functional modules data available. Please 'Merge modules' first.",
            easyClose = TRUE,
            footer = NULL
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

        output$enirched_functional_module_plot <-
          shiny::renderPlot({
            req(tryCatch(
              plot,
              error = function(e)
                NULL
            ))
          })

        shinyjs::hide("loading")

      }
    })

    ###--------------------------------------------------------------------
    ###step 5 Data visualization
    # Observe file upload
    # enriched_functional_module <-
    #   reactiveVal()
    observeEvent(input$upload_enriched_functional_module, {
      if (!is.null(input$upload_enriched_functional_module$datapath)) {
        tempEnv <- new.env()
        load(input$upload_enriched_functional_module$datapath, envir = tempEnv)

        names <- ls(tempEnv)
          # Handle user response
        enriched_functional_module(get(names[1], envir = tempEnv))
        } else {
          showModal(modalDialog(
            title = "Error",
            "The uploaded file should contain exactly one object.",
            easyClose = TRUE,
            footer = NULL
          ))
        }
    })

    #####barplot
    # Observe generate barplot button click
    observeEvent(input$generate_barplot, {
      if (is.null(enriched_functional_module())) {
        # No enriched functional module available
        showModal(modalDialog(
          title = "Warning",
          "No enriched functional module data available. Please complete the previous steps or upload the data.",
          easyClose = TRUE,
          footer = NULL
        ))
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
            database_color = c(GO = input$barplot_go_color,
                               KEGG = input$barplot_kegg_color,
                               Reactome = input$barplot_reactome_color)
          )

        ###save functions
        data_color <-
        paste0("c(",paste(paste(c("GO", "KEGG", "Reactome"),
                                c(paste0('"',input$barplot_go_color,'"'),
                                  paste0('"',input$barplot_kegg_color,'"'),
                                  paste0('"',input$barplot_reactome_color,'"')),
                                sep = " = "),
                          collapse = ", "),")")
        barplot_database <-
          paste0("c(",paste(unlist(lapply(paste(input$barplot_database), function(x)paste0('"',x,'"'))),
                            collapse = ", "), ")")

        barplot_code <-
          sprintf("plot_pathway_bar(
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
                  data_color)

        # Save to a text file
        writeLines(barplot_code, "files/barplot_code.txt")

        output$barplot <-
          shiny::renderPlot({
            req(plot)
          })

        shinyjs::hide("loading")

        output$download_barplot <-
          downloadHandler(
          filename = function() {
            paste0("pathway_barplot.", input$barplot_type)
          },
          content = function(file) {
            ggsave(file, plot = plot,
                   width = input$barplot_width,
                   height = input$barplot_height)
          }
        )

      }
    })

    ######code for barplot
    output$barplot_code <- renderUI({
      if (input$show_barplot_code %% 2 == 1) {
        code_content <-
          readLines("files/barplot_code.txt")
        code_content <-
          paste(code_content, collapse = "\n")
        tags$pre(code_content)
      }
    })


  }
