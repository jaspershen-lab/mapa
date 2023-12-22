# ####example
# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("enriched_functional_module.rda")
# load("intepretation_result.rda")
#
# report_functional_module(object = object,
#                          intepretation_result = interpretation_result,
#                          path = ".",
#                          type = "all")

#' Generate a Report for a Functional Module
#'
#' This function generates a report for a given object of class 'functional_module'.
#' It creates the report in different formats (HTML, PDF, Word) based on user choice.
#' The function also performs various tasks such as creating directories,
#' fetching report templates, saving parameter and object data,
#' and rendering bar plots for functional modules, modules, and pathways.
#'
#' @param object An object of class 'functional_module'. The function will
#'               generate a report based on this object.
#' @param intepretation_result Placeholder for future functionality or data
#'                             that might be included in the report.
#' @param path A character string specifying the directory path where the report
#'             and associated files will be saved. Defaults to the current directory.
#' @param type A character vector specifying the format of the report.
#'             Can be 'html', 'pdf', 'word', or 'all' to generate the report
#'             in all three formats. Default is 'html'.
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
#' @importFrom massdataset parse_tidymass_parameter
#'
#' @export

report_functional_module <-
  function(object,
           intepretation_result,
           path = ".",
           type = c("html", "pdf", "word", "all")) {
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
      output_path = file.path(path, paste('Report', length(grep(
        "Report", dir(path)
      )) + 1, sep = "_"))
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

    ###parameters
    message("Parameters.")

    parameters <-
      object@process_info %>%
      lapply(function(x) {
        if (length(x) == 1) {
          massdataset::parse_tidymass_parameter(object = x)
        } else{
          x %>%
            lapply(function(y) {
              massdataset::parse_tidymass_parameter(object = y)
            }) %>%
            dplyr::bind_rows()
        }
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::arrange(time)

    save(parameters, file = file.path(output_path, "parameters.rda"))
    save(object, file = file.path(output_path, "object.rda"))

    ####Barplot to show the pathways
    message("Barplot")

    plot_functional_module <-
      plot_pathway_bar(
        object = object,
        top_n = 10,
        p.adjust.cutoff = 0.05,
        count.cutoff = 5,
        y_lable_width = 30,
        level = "functional_module"
      ) +
      labs(title = "Functional module")

    plot_module <-
      plot_pathway_bar(
        object = object,
        top_n = 10,
        p.adjust.cutoff = 0.05,
        count.cutoff = 5,
        y_lable_width = 30,
        level = "module"
      ) +
      labs(title = "Module")

    plot_pathway <-
      plot_pathway_bar(
        object = object,
        top_n = 10,
        p.adjust.cutoff = 0.05,
        count.cutoff = 5,
        y_lable_width = 30,
        level = "pathway"
      ) +
      labs(title = "Pathway")

    ggplot2::ggsave(
      filename = file.path(output_path, "plot_functional_module.png"),
      plot = plot_functional_module,
      width = 8,
      height = 6,
      dpi = 600
    )

    ggplot2::ggsave(
      filename = file.path(output_path, "plot_module.png"),
      plot = plot_module,
      width = 8,
      height = 6
    )

    ggplot2::ggsave(
      filename = file.path(output_path, "plot_pathway.png"),
      plot = plot_pathway,
      width = 8,
      height = 6
    )

    ####Whole module network
    similarity_network_go <-
      plot_similarity_network(object = object,
                              level = "module",
                              database = "go") +
      labs(title = "GO")

    similarity_network_reactome <-
      plot_similarity_network(object = object,
                              level = "module",
                              database = "reactome") +
      labs(title = "Reactome")

    similarity_network_kegg <-
      plot_similarity_network(object = object,
                              level = "module",
                              database = "kegg") +
      labs(title = "KEGG")

    similarity_network_function_module <-
      plot_similarity_network(object = object,
                              level = "functional_module",
                              database = "kegg") +
      labs(title = "Functional modules")

    ggplot2::ggsave(
      filename = file.path(output_path, "similarity_network_function_module.png"),
      plot = similarity_network_function_module,
      width = 8,
      height = 6
    )

    ggplot2::ggsave(
      filename = file.path(output_path, "similarity_network_go.png"),
      plot = similarity_network_go,
      width = 8,
      height = 6
    )

    ggplot2::ggsave(
      filename = file.path(output_path, "similarity_network_reactome.png"),
      plot = similarity_network_reactome,
      width = 8,
      height = 6
    )

    ggplot2::ggsave(
      filename = file.path(output_path, "similarity_network_kegg.png"),
      plot = similarity_network_kegg,
      width = 8,
      height = 6
    )

    ####Module network analysis
    ###
    functional_module_id <-
      enriched_functional_module@merged_module$functional_module_result %>%
      dplyr::filter(p.adjust < 0.05 & Count >= 5) %>%
      dplyr::arrange(p.adjust) %>%
      head(10) %>%
      pull(module)

    if (length(functional_module_id) > 0) {
      functional_module_length <-
        lapply(functional_module_id, function(x) {
          sum(enriched_functional_module@merged_module$result_with_module$module == x)
        }) %>%
        unlist()

      functional_module_id <-
        functional_module_id[functional_module_length > 1]
    }

    if (length(functional_module_id) > 0) {
      plot <-
        plot_module_info(object = enriched_functional_module,
                         level = "functional_module",
                         module_id = functional_module_id[1])
    }

    message("Render report.")

    if (missing(intepretation_result)) {
      intepretation_result <- "> No interpretation result is provided."
    }

    ##transform rmd to HTML or pdf
    if (type == "html" | type == "all") {
      rmarkdown::render(
        file.path(output_path, "mapa.template.Rmd"),
        output_format = rmarkdown::html_document(),
        params = list(text_data = interpretation_result)
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
        params = list(text_data = interpretation_result)
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
        params = list(text_data = interpretation_result)
      )
      file.rename(
        from = file.path(output_path, "mapa.template.pdf"),
        to = file.path(output_path, "mapa_report.pdf")
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
