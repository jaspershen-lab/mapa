# ####example
# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("enriched_functional_module.rda")
# load("interpretation_result.rda")
#
# report_functional_module(
#   object = met_object,
#   degree_cutoff = 0,
#   path = "demo_data/embedding_output/",
#   type = "html"
# )

#' Generate a Report for a Functional Module
#'
#' This function generates a report for a given object of class 'functional_module'.
#' It creates the report in different formats (HTML, PDF, Word) based on user choice.
#' The function also performs various tasks such as creating directories,
#' fetching report templates, saving parameter and object data,
#' and rendering bar plots for functional modules, modules, and pathways, similarity networks,
#' and interpretation of functional modules.
#'
#' @param object An object of class 'functional_module'. The function will
#'               generate a report based on this object.
#' @param path A character string specifying the directory path where the report
#'             and associated files will be saved. Defaults to the current directory.
#' @param type A character vector specifying the format of the report.
#'             Can be 'html', 'pdf', 'md', 'word', or 'all' to generate the report
#'             in all three formats. Default is 'html'.
#' @param degree_cutoff A numeric value specifying the minimum degree for network
#'                     analysis. Default is 1. Must be 0 or a positive number.
#'
#' @return This function does not return a value but generates a report in the
#'         specified format(s) and saves it to the given path.
#'
#' @details The function first checks for the presence and class of the 'object'
#'          parameter. It then creates the necessary directory structure, fetches
#'          the report template using the 'draft' function from the 'rmarkdown' package,
#'          and saves parameters and object data. It also creates and saves bar plots
#'          for different levels (functional_module, module, pathway) using 'ggplot2'.
#'          Finally, it renders the report in the specified format(s) and cleans up
#'          intermediate files.
#'
#'
#' @importFrom rmarkdown draft render html_document pdf_document word_document
#' @importFrom ggplot2 ggsave
#' @importFrom dplyr bind_rows arrange
#'
#' @export

report_functional_module <-
  function(object,
           path = ".",
           degree_cutoff = 1,
           type = c("html", "pdf", "word", "md", "all")) {

    if (identical(type, "pdf") && Sys.which("pdflatex") == "") {
      stop("PDF output requires a LaTeX distribution. ",
           "Install TinyTeX via install.packages('tinytex') and tinytex::install_tinytex(), ",
           "or re-run with type = 'html' or 'word'.")
    }

    has_kableExtra <- requireNamespace("kableExtra", quietly = TRUE)

    if (!has_kableExtra && identical(type, "pdf")) {
      warning("kableExtra not found; tables will not scale to page width in PDF. ",
              "Run install.packages('kableExtra') for prettier tables.")
    }

    if (missing(object)) {
      stop("object is missing")
    }

    if (!is(object, "functional_module")) {
      stop("object must be functional_module class")
    }

    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    type <- match.arg(type)

    options(warn = -1)

    ###path
    if (length(grep("Report", dir(path))) > 0) {
      idx <-
      max(
        as.numeric(stringr::str_extract(
          grep(pattern = "Report", dir(path), value = TRUE),
          "[0-9]{1,10}"
        )), na.rm = TRUE
      )

      if(is.na(idx)){
        idx <- 0
      }

      if (!is.finite(idx))          # catches -Inf as well as Inf
        idx <- 0

      output_path <-
        file.path(path, paste('Report', idx + 1, sep = "_"))
    } else{
      output_path <- file.path(path, "Report")
    }

    ####get the template
    message("Get report template.")

    rmarkdown::draft(
      file = output_path,
      template = "mapa",
      package = "mapa",
      create_dir = TRUE,
      edit = FALSE
    )

    ## Parameters ====
    message("Saving parameters ...")

    parameters <-
      object@process_info %>%
      lapply(function(x) {
        if (length(x) == 1) {
          mapa_parse_tidymass_parameter(object = x)
        } else{
          x %>%
            lapply(function(y) {
              mapa_parse_tidymass_parameter(object = y)
            }) %>%
            dplyr::bind_rows()
        }
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(time)

    save(parameters, file = file.path(output_path, "parameters.rda"))
    save(object, file = file.path(output_path, "object.rda"))

    ## Barplot to show the top 10 pathways ====
    message("Saving barplots ...")

    if (length(c(object@merged_pathway_go,
                 object@merged_pathway_kegg,
                 object@merged_pathway_reactome,
                 object@merged_pathway_hmdb,
                 object@merged_pathway_metkegg)) == 0) {
      module <- FALSE
    } else {
      module <- TRUE
    }

    if ("llm_interpret_module" %in% names(object@process_info)) {
      llm_text <- TRUE
    } else {
      llm_text <- FALSE
    }

    plot_functional_module <-
      plot_pathway_bar(
        object = object,
        top_n = 10,
        p.adjust.cutoff = 0.05,
        count.cutoff = 5,
        y_label_width = 30,
        level = "functional_module",
        llm_text = llm_text
      )

    if (module) {
      plot_module <-
        plot_pathway_bar(
          object = object,
          top_n = 10,
          p.adjust.cutoff = 0.05,
          count.cutoff = 5,
          y_label_width = 30,
          level = "module"
        )
    }

    plot_pathway <-
      plot_pathway_bar(
        object = object,
        top_n = 10,
        p.adjust.cutoff = 0.05,
        count.cutoff = 5,
        y_label_width = 30,
        level = "pathway"
      )

    tryCatch({
      ggplot2::ggsave(
        filename = file.path(output_path, "plot_functional_module.png"),
        plot = plot_functional_module,
        width = 8,
        height = 6,
        dpi = 600
      )
      message("Plot plot_functional_module.png saved successfully to: ", file.path(output_path, "plot_functional_module.png"), "!")
    }, error = function(e) {
      message("Error saving plot: ", e$message)
    })

    if (!module) {
      message("No module level result was generated - no plot to save.")
    } else {
      tryCatch({
        ggsave(
          filename = file.path(output_path, "plot_module.png"),
          plot = plot_module,
          width = 8,
          height = 6
        )
        message("Plot plot_module.png saved successfully to: ", file.path(output_path, "plot_module.png"), "!")
      }, error = function(e) {
        message("Error saving plot: ", e$message)
      })
    }

    tryCatch({
      ggplot2::ggsave(
        filename = file.path(output_path, "plot_pathway.png"),
        plot = plot_pathway,
        width = 8,
        height = 6
      )
      message("Plot plot_pathway.png saved successfully to: ", file.path(output_path, "plot_pathway.png"), "!")
    }, error = function(e) {
      message("Error saving plot: ", e$message)
    })

    ## Whole module network ====
    message("Saving similarity networks ...")

    if (length(object@merged_pathway_go) != 0) {
      if (sum(object@merged_pathway_go$module_result$module_content_number > degree_cutoff) > 0) {
        similarity_network_go <-
          plot_similarity_network(object = object,
                                  level = "module",
                                  degree_cutoff = degree_cutoff,
                                  database = "go") +
          labs(title = "GO Modules")

        tryCatch({
          ggplot2::ggsave(
            filename = file.path(output_path, "similarity_network_go.png"),
            plot = similarity_network_go,
            width = 8,
            height = 6
          )
          message("Plot similarity_network_go.png saved successfully to: ", file.path(output_path, "similarity_network_go.png"), "!")
        }, error = function(e) {
          message("Error saving plot: ", e$message)
        })
      } else {
        message("GO similarity network plot not generated: No modules with content number > ", degree_cutoff)
      }
    } else {
      message("GO similarity network plot not generated: No GO pathway data available")
    }

    if (length(object@merged_pathway_kegg) != 0) {
      if (sum(object@merged_pathway_kegg$module_result$module_content_number > degree_cutoff) > 0) {
        similarity_network_kegg <-
          plot_similarity_network(object = object,
                                  level = "module",
                                  degree_cutoff = degree_cutoff,
                                  database = "kegg") +
          ggplot2::labs(title = "KEGG Modules")

        tryCatch({
          ggplot2::ggsave(
            filename = file.path(output_path, "similarity_network_kegg.png"),
            plot = similarity_network_kegg,
            width = 8,
            height = 6
          )
          message("Plot similarity_network_kegg.png saved successfully to: ", file.path(output_path, "similarity_network_kegg.png"), "!")
        }, error = function(e) {
          message("Error saving plot: ", e$message)
        })
      } else {
        message("KEGG similarity network plot not generated: No modules with content number > ", degree_cutoff)
      }
    } else {
      message("KEGG similarity network plot not generated: No KEGG pathway data available")
    }

    if (length(object@merged_pathway_reactome) != 0) {
      if (sum(object@merged_pathway_reactome$module_result$module_content_number > degree_cutoff) > 0) {
        similarity_network_reactome <-
          plot_similarity_network(object = object,
                                  level = "module",
                                  degree_cutoff = degree_cutoff,
                                  database = "reactome") +
          ggplot2::labs(title = "Reactome Modules")

        tryCatch({
          ggplot2::ggsave(
            filename = file.path(output_path, "similarity_network_reactome.png"),
            plot = similarity_network_reactome,
            width = 8,
            height = 6
          )
          message("Plot similarity_network_reactome.png saved successfully to: ", file.path(output_path, "similarity_network_reactome.png"), "!")
        }, error = function(e) {
          message("Error saving plot: ", e$message)
        })
      } else {
        message("Reactome similarity network plot not generated: No modules with content number > ", degree_cutoff)
      }
    } else {
      message("Reactome similarity network plot not generated: No Reactome pathway data available")
    }

    if (length(object@merged_pathway_hmdb) != 0) {
      if (sum(object@merged_pathway_hmdb$module_result$module_content_number > degree_cutoff) > 0) {
        similarity_network_hmdb <-
          plot_similarity_network(object = object,
                                  level = "module",
                                  degree_cutoff = degree_cutoff,
                                  database = "hmdb") +
          ggplot2::labs(title = "HMDB Modules")

        tryCatch({
          ggplot2::ggsave(
            filename = file.path(output_path, "similarity_network_hmdb.png"),
            plot = similarity_network_hmdb,
            width = 8,
            height = 6
          )
          message("Plot similarity_network_hmdb.png saved successfully to: ", file.path(output_path, "similarity_network_hmdb.png"), "!")
        }, error = function(e) {
          message("Error saving plot: ", e$message)
        })
      } else {
        message("HMDB similarity network plot not generated: No modules with content number > ", degree_cutoff)
      }
    } else {
      message("HMDB similarity network plot not generated: No HMDB pathway data available")
    }

    if (length(object@merged_pathway_metkegg) != 0) {
      if (sum(object@merged_pathway_metkegg$module_result$module_content_number > degree_cutoff) > 0) {
        similarity_network_metkegg <-
          plot_similarity_network(object = object,
                                  level = "module",
                                  degree_cutoff = degree_cutoff,
                                  database = "metkegg") +
          ggplot2::labs(title = "KEGG Modules")

        tryCatch({
          ggplot2::ggsave(
            filename = file.path(output_path, "similarity_network_metkegg.png"),
            plot = similarity_network_metkegg,
            width = 8,
            height = 6
          )
          message("Plot similarity_network_metkegg.png saved successfully to: ", file.path(output_path, "similarity_network_metkegg.png"), "!")
        }, error = function(e) {
          message("Error saving plot: ", e$message)
        })
      } else {
        message("MetKEGG similarity network plot not generated: No modules with content number > ", degree_cutoff)
      }
    } else {
      message("MetKEGG similarity network plot not generated: No MetKEGG pathway data available")
    }

    if (sum(object@merged_module$functional_module_result$module_content_number > degree_cutoff) > 0) {
      similarity_network_function_module <-
        plot_similarity_network(object = object,
                                level = "functional_module",
                                degree_cutoff = degree_cutoff,
                                llm_text = llm_text) +
        labs(title = "Functional Modules")

      tryCatch({
        ggplot2::ggsave(
          filename = file.path(output_path, "similarity_network_function_module.png"),
          plot = similarity_network_function_module,
          width = 8,
          height = 6
        )
        message("Plot similarity_network_function_module.png saved successfully to: ", file.path(output_path, "similarity_network_function_module.png"), "!")
      }, error = function(e) {
        message("Error saving plot: ", e$message)
      })
    } else {
      message("Functional module similarity network plot not generated: No modules with content number > ", degree_cutoff)
    }

    ## Module network analysis ====
    # functional_module_id <-
    #   object@merged_module$functional_module_result %>%
    #   dplyr::filter(p_adjust < 0.05 & Count >= 5) %>%
    #   dplyr::arrange(p_adjust) %>%
    #   head(10) %>%
    #   pull(module)
    #
    # if (length(functional_module_id) > 0) {
    #   functional_module_length <-
    #     lapply(functional_module_id, function(x) {
    #       sum(object@merged_module$result_with_module$module == x)
    #     }) %>%
    #     unlist()
    #
    #   functional_module_id <-
    #     functional_module_id[functional_module_length > 1]
    # }
    #
    # if (length(functional_module_id) > 0) {
    #   plot <-
    #     plot_module_info(object = object,
    #                      level = "functional_module",
    #                      module_id = functional_module_id[1])
    # }

    ## Interpretation of functional modules ====
    message("Saving the interpretation of functional modules ...")

    if (llm_text) {
      interpretation_result <- extract_llm_module_data(llm_module_interpretation = object@llm_module_interpretation)
    } else {
      interpretation_result <- object@merged_module$functional_module_result
    }

    save(interpretation_result, file = file.path(output_path, "interpretation_result.rda"))

    message("Rendering report ...")

    ##transform rmd to HTML or pdf
    if (type == "html" | type == "all") {
      rmarkdown::render(
        file.path(output_path, "mapa.template.Rmd"),
        output_format = rmarkdown::html_document(),
        params = list(text_data = NULL)
      )

      file.rename(
        from = file.path(output_path, "mapa.template.html"),
        to = file.path(output_path, "mapa_report.html")
      )
    }

    ########render rmarkddown to html or pdf
    if (type == "pdf" | type == "all") {
      rmarkdown::render(
        input = file.path(output_path, "mapa.template.Rmd"),
        output_format = rmarkdown::pdf_document(),
        params = list(text_data = NULL)
      )
      file.rename(
        from = file.path(output_path, "mapa.template.pdf"),
        to = file.path(output_path, "mapa_report.pdf")
      )
    }

    ########render rmarkddown to word or all
    if (type == "word" | type == "all") {
      rmarkdown::render(
        file.path(output_path, "mapa.template.Rmd"),
        output_format = rmarkdown::word_document(),
        params = list(text_data = NULL)
      )
      file.rename(
        from = file.path(output_path, "mapa.template.docx"),
        to = file.path(output_path, "mapa_report.docx")
      )
    }

    ########render rmarkddown to md or all
    if (type == "md" | type == "all") {
      rmarkdown::render(
        file.path(output_path, "mapa.template.Rmd"),
        output_format = rmarkdown::md_document(),
        params = list(text_data =  NULL)
      )
      file.rename(
        from = file.path(output_path, "mapa.template.md"),
        to = file.path(output_path, "mapa_report.md")
      )
    }

    ####remove some files
    message("Remove some files.")
    file = dir(output_path)
    remove_file = grep("png|Rmd|parameters|rda", file, value = TRUE)
    unlink(
      x = file.path(output_path, remove_file),
      recursive = TRUE,
      force = TRUE
    )
}
