# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# if (any(dir("result/go/intermediate_data/") == "module_result")) {
#   load("result/go/intermediate_data/module_result")
#   module_result_go <- module_result
#
# } else{
#   module_result_go <- NULL
# }
#
# if (any(dir("result/kegg/intermediate_data/") == "module_result")) {
#   load("result/kegg/intermediate_data/module_result")
#   module_result_kegg <- module_result
#
# } else{
#   module_result_kegg <- NULL
# }
#
# if (any(dir("result/reactome/intermediate_data/") == "module_result")) {
#   load("result/reactome/intermediate_data/module_result")
#   module_result_reactome <- module_result
# } else{
#   module_result_reactome <- NULL
# }
#
# load("demo_data")
#
#
# merge_modules(
#   object = demo_data,
#   module_result_go = module_result_go,
#   module_result_kegg = module_result_kegg,
#   module_result_reactome = module_result_reactome,
#   sim.cutoff = 0.5,
#   measure_method = "jaccard"
# )

#' merge_modules Function
#'
#' This function merges information from three databases (GO, KEGG, Reactome) and performs similarity analysis among pathways/modules across these databases.
#' It filters, rearranges, and aggregates the data before computing the Jaccard Index for similarity.
#' The function also creates network visualizations and saves the results to an Excel file.
#'
#' @param object An mass_dataset object containing gene information.
#' @param module_result_go A data frame containing module results from GO database.
#' @param module_result_kegg A data frame containing module results from KEGG database.
#' @param module_result_reactome A data frame containing module results from Reactome database.
#' @param sim.cutoff A numeric value for similarity cutoff. Default is 0.5.
#' @param measure_method A character vector specifying the method of similarity measurement. Default is "jaccard".
#' @param path A character string specifying the directory to save results. Default is "result".
#'
#' @return Outputs an Excel file and a PDF containing the merged and analyzed data.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @export

merge_modules <-
  function(object,
           module_result_go,
           module_result_kegg,
           module_result_reactome,
           sim.cutoff = 0.5,
           measure_method = c("jaccard"),
           path = "result") {
    path <- file.path(path, "function_modules")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)

    variable_info <-
      massdataset::extract_variable_info(object)

    ######calculate the similarity (jaccard index) between all the pathways
    if (!is.null(module_result_go)) {
      module_result_go <-
        module_result_go %>%
        dplyr::filter(ONTOLOGY != "CC") %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "GO") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }

    if (!is.null(module_result_kegg)) {
      module_result_kegg <-
        module_result_kegg %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "KEGG") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }

    if (!is.null(module_result_reactome)) {
      module_result_reactome <-
        module_result_reactome %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(database = "Reactome") %>%
        dplyr::select(
          module_annotation,
          module,
          Description,
          BgRatio,
          pvalue,
          qvalue,
          p.adjust,
          Count,
          database,
          geneID,
          pathway_id = node
        ) %>%
        dplyr::mutate(Count = as.numeric(Count)) %>%
        dplyr::filter(!is.na(module_annotation))
    }

    message("Calculating the similarity matrix...")
    jaccard_index <-
      get_jaccard_index_for_three_databases(
        module_result_go = module_result_go,
        module_result_kegg = module_result_kegg,
        module_result_reactome = module_result_reactome,
        object = object
      )

    edge_data =
      jaccard_index %>%
      dplyr::filter(value > sim.cutoff) %>%
      dplyr::rename(from = name1, to = name2, sim = value)

    node_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome) %>%
      dplyr::select(module, dplyr::everything()) %>%
      dplyr::rename(node = module)

    graph_data <-
      tidygraph::tbl_graph(nodes = node_data,
                           edges = edge_data,
                           directed = FALSE) %>%
      dplyr::mutate(degree = tidygraph::centrality_degree())

    subnetwork <-
      igraph::cluster_edge_betweenness(graph = graph_data,
                                       weights = abs(edge_attr(graph_data,
                                                               "sim")))
    cluster <-
      paste("Function_module", as.character(igraph::membership(subnetwork)), sep = "_")

    graph_data <-
      graph_data %>%
      tidygraph::mutate(module = cluster)

    ###clustered different GO terms
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(module, p.adjust)

    ##add module content number
    module_content_number <-
      result_with_module %>%
      dplyr::count(module) %>%
      dplyr::rename(module_content_number = n)

    result_with_module <-
      result_with_module %>%
      dplyr::left_join(module_content_number, by = "module")

    dir.create(
      file.path(path, "intermediate_data"),
      showWarnings = FALSE,
      recursive = TRUE
    )
    save(result_with_module,
         file = file.path(path, "intermediate_data/result_with_module"))

    graph_data <-
      graph_data %>%
      activate(what = "nodes") %>%
      dplyr::left_join(module_content_number, by = "module")

    save(graph_data, file = file.path(path, "intermediate_data/graph_data"))

    module_result <-
      result_with_module %>%
      plyr::dlply(.variables = .(module)) %>%
      purrr::map(function(x) {
        # cat(unique(x$module), " ")
        if (nrow(x) == 1) {
          x$module_content <-
            paste(x$node, collapse = ";")
          x <-
            x %>%
            dplyr::select(module, everything()) %>%
            dplyr::distinct(module, .keep_all = TRUE) %>%
            dplyr::select(-node)
          return(x)
        }

        x =
          x %>%
          dplyr::arrange(p.adjust)

        x$module_content <-
          paste(x$node, collapse = ";")

        x$Description <-
          paste(x$Description, collapse = ";")

        x$BgRatio <-
          paste(x$BgRatio, collapse = ";")

        x$pvalue <- min(as.numeric(x$pvalue))
        x$p.adjust <- min(as.numeric(x$p.adjust))
        x$qvalue <- min(as.numeric(x$qvalue))
        x$geneID =
          x$geneID %>%
          stringr::str_split(pattern = "/") %>%
          unlist() %>%
          unique() %>%
          paste(collapse = '/')

        x$Count <-
          length(stringr::str_split(x$geneID[1], pattern = "/")[[1]])

        x <-
          x %>%
          dplyr::select(module, everything()) %>%
          dplyr::distinct(module, .keep_all = TRUE) %>%
          dplyr::select(-node)
        x
      }) %>%
      do.call(rbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(
        module_annotation = case_when(
          module == "Other" ~ Description,
          module != "Other" ~ stringr::str_split(Description, ";")[[1]][1]
        )
      ) %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(module_annotation, everything())

    module_result$module_annotation <-
      stringr::str_split(module_result$Description, ";") %>%
      purrr::map(function(x) {
        x[1]
      }) %>%
      unlist()

    save(module_result,
         file = file.path(path, "intermediate_data/module_result"))

    wb = openxlsx::createWorkbook()
    openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Times New Roma")
    addWorksheet(wb, sheetName = "enriched_pathway_result", gridLines = TRUE)
    addWorksheet(wb, sheetName = "enriched_module_result", gridLines = TRUE)
    freezePane(wb,
               sheet = 1,
               firstRow = TRUE,
               firstCol = TRUE)
    freezePane(wb,
               sheet = 2,
               firstRow = TRUE,
               firstCol = TRUE)
    writeDataTable(
      wb,
      sheet = 1,
      x = result_with_module,
      colNames = TRUE,
      rowNames = FALSE
    )

    writeDataTable(
      wb,
      sheet = 2,
      x = module_result,
      colNames = TRUE,
      rowNames = FALSE
    )

    saveWorkbook(wb,
                 file = file.path(path, "enriched_result.xlsx"),
                 overwrite = TRUE)

    ###plot to show the clusters of GO terms
    cluster_label_module <-
      igraph::as_data_frame(graph_data, what = "vertices") %>%
      dplyr::group_by(module) %>%
      dplyr::filter(p.adjust == min(p.adjust) &
                      Count == max(Count)) %>%
      dplyr::slice_head(n = 1) %>%
      pull(Description)

    cluster_label_all <-
      node_data$Description

    plot <-
      graph_data %>%
      ggraph(layout = 'fr',
             circular = FALSE) +
      geom_edge_link(
        aes(width = sim),
        strength = 1,
        color = "black",
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_point(
        aes(fill = database,
            size = -log(p.adjust, 10)),
        shape = 21,
        alpha = 1,
        show.legend = TRUE
      ) +
      geom_node_text(aes(
        x = x,
        y = y,
        label = ifelse(Description %in% cluster_label_module, Description, NA)
      ),
      size = 3,
      repel = TRUE) +
      guides(fill = guide_legend(ncol = 1)) +
      scale_edge_width_continuous(range = c(0.1, 2)) +
      scale_size_continuous(range = c(1, 7)) +
      scale_fill_manual(values = database_color) +
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA),
      )

    plot

    extrafont::loadfonts()

    ggsave(
      plot,
      filename = file.path(path, "/combined_similarity_plot.pdf"),
      width = 9,
      height = 7
    )

    message("\nDone")

  }
