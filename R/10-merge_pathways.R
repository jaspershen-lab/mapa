# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
# library(ggraph)
# library(igraph)
# library(tidygraph)
# library(tidyverse)
# library(extrafont)
# library(simplifyEnrichment)
# library(plyr)
# library(openxlsx)
# library(ggwordcloud)
# library(patchwork)
# library(pRoloc)
# library(GOSim)
#
# load("result/enrichment_go_result")
# load("result/enrichment_kegg_result")
# load("result/enrichment_reactome_result")
#
# merge_pathways(pathway_result = enrichment_go_result,
#                p.adjust.cutoff = 0.05,
#                count.cutoff = 5,
#                database = "go",
#                sim.cutoff = 0.5,
#                measure_method = "Wang",
#                path = "result")
#
# merge_pathways(pathway_result = enrichment_kegg_result,
#                p.adjust.cutoff = 0.05,
#                count.cutoff = 5,
#                database = "kegg",
#                sim.cutoff = 0.5,
#                path = "result")
#
# merge_pathways(pathway_result = enrichment_reactome_result,
#                p.adjust.cutoff = 0.05,
#                count.cutoff = 5,
#                database = "reactome",
#                sim.cutoff = 0.5,
#                path = "result")

#' Merge Pathways and Analyze Semantic Similarity
#'
#' This function performs a merging of pathways from enrichment analysis, calculates
#' semantic similarity among terms, and identifies modules within the pathways.
#' It also generates various plots and outputs results into organized files.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @param pathway_result A result object from pathway enrichment analysis.
#' @param p.adjust.cutoff A numeric value for the adjusted p-value cutoff (default: 0.05).
#' @param count.cutoff A numeric value for the count cutoff (default: 5).
#' @param database A character vector specifying the database (default: c("go", "kegg", "reactome")).
#' @param sim.cutoff A numeric value for the similarity cutoff (default: 0.5).
#' @param measure_method A character vector specifying the similarity measure method (default: c("Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS", "jaccard")).
#' @param path A character string specifying the path for result output (default: "result").
#' @return NULL. Results are written to files.
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @examples
#' \dontrun{
#' # Assuming `pathway_result` is the result from your pathway enrichment analysis.
#' merge_pathways(pathway_result = pathway_result)
#' }
#' @export

merge_pathways <-
  function(pathway_result,
           p.adjust.cutoff = 0.05,
           count.cutoff = 5,
           database = c("go", "kegg", "reactome"),
           sim.cutoff = 0.5,
           measure_method = c("Wang", "Resnik",
                              "Rel", "Jiang",
                              "Lin", "TCSS",
                              "jaccard"),
           path = "result") {
    measure_method <-
      match.arg(measure_method)
    database <- match.arg(database)
    path <- file.path(path, database)

    if (missing(pathway_result)) {
      stop("pathway_result is required")
    }

    if (!is(pathway_result, "enrichResult")) {
      stop("pathway_result must be result from enrich_pathway function")
    }

    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    dir.create(
      file.path(path, "intermediate_data"),
      showWarnings = FALSE,
      recursive = TRUE
    )

    result <-
      pathway_result@result %>%
      dplyr::filter(p.adjust < p.adjust.cutoff) %>%
      dplyr::filter(Count > count.cutoff) %>%
      dplyr::arrange(p.adjust)

    if (nrow(result) == 0) {
      return(NULL)
    }

    ##get the similartiy matrix
    ######get the similarity between GO terms
    message("Calculating similartiy matrix, it may take a while...")

    if (any(dir(file.path(path, "intermediate_data")) == "sim_matrix")) {
      load(file.path(path, "intermediate_data/sim_matrix"))
    } else{
      if (database == "go") {
        sim_matrix <-
          get_go_result_sim(
            result = result,
            sim.cutoff = sim.cutoff,
            measure_method = measure_method
          )
      }

      if (database == "kegg") {
        sim_matrix <-
          tryCatch(
            sim_matrix <-
              simplifyEnrichment::term_similarity_from_KEGG(term_id = c(result$ID),
                                                            method = "jaccard") %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "name1") %>%
              tidyr::pivot_longer(
                cols = -name1,
                names_to = "name2",
                values_to = "sim"
              ) %>%
              dplyr::filter(name1 != name2),
            error = function(x) {
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }

      if (database == "reactome") {
        sim_matrix <-
          tryCatch(
            sim_matrix <-
              simplifyEnrichment::term_similarity_from_Reactome(term_id = c(result$ID),
                                                                method = "jaccard") %>%
              as.data.frame() %>%
              tibble::rownames_to_column(var = "name1") %>%
              tidyr::pivot_longer(
                cols = -name1,
                names_to = "name2",
                values_to = "sim"
              ) %>%
              dplyr::filter(name1 != name2),
            error = function(x) {
              data.frame(name1 = character(),
                         name2 = character(),
                         sim = numeric())
            }
          )
      }

      save(sim_matrix, file = file.path(path, "intermediate_data/sim_matrix"))
    }

    ####module detection
    message("Identifying modules...")

    if (any(dir(file.path(path, "intermediate_data")) == "graph_data")) {
      load(file.path(path, "intermediate_data/graph_data"))
    } else{
      edge_data <-
        rbind(sim_matrix) %>%
        dplyr::rename(from = name1, to = name2) %>%
        dplyr::filter(sim > sim.cutoff)

      node_data <-
        rbind(result) %>%
        as.data.frame() %>%
        dplyr::select(ID, everything()) %>%
        dplyr::rename(node = ID)

      graph_data <-
        tidygraph::tbl_graph(nodes = node_data,
                             edges = edge_data,
                             directed = FALSE) %>%
        dplyr::mutate(degree = tidygraph::centrality_degree())

      subnetwork <-
        igraph::cluster_edge_betweenness(graph = graph_data,
                                         weights = abs(edge_attr(graph_data,
                                                                 "sim")))

      # save(subnetwork, file = file.path(path, "subnetwork"))
      cluster <-
        paste(database,
              "Module",
              as.character(igraph::membership(subnetwork)),
              sep = "_")

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

      if (database == "go") {
        result_with_module <-
          result_with_module %>%
          dplyr::arrange(ONTOLOGY, module, p.adjust)
      }

      ###add module content number
      module_content_number <-
        result_with_module %>%
        dplyr::count(module) %>%
        dplyr::rename(module_content_number = n)

      result_with_module <-
        result_with_module %>%
        dplyr::left_join(module_content_number, by = "module")

      save(result_with_module,
           file = file.path(path, "intermediate_data/result_with_module"))
      graph_data <-
        graph_data %>%
        activate(what = "nodes") %>%
        dplyr::left_join(module_content_number, by = "module")

      save(graph_data, file = file.path(path, "intermediate_data/graph_data"))
    }

    ###plot to show the clusters of GO terms
    cluster_label_module <-
      igraph::as_data_frame(graph_data, what = "vertices") %>%
      dplyr::group_by(module) %>%
      dplyr::filter(p.adjust == min(p.adjust) &
                      Count == max(Count)) %>%
      dplyr::slice_head(n = 1) %>%
      pull(Description)

    cluster_label_all <-
      igraph::as_data_frame(graph_data, what = "vertices")$Description

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
        aes(fill = module,
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
      ggraph::theme_graph() +
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        legend.position = "right",
        legend.background = element_rect(fill = "transparent", color = NA)
      )

    extrafont::loadfonts()

    ggplot2::ggsave(
      plot,
      filename =
        file.path(path, "similarity_network_plot.pdf"),
      width = 9,
      height = 7
    )

    ###output some files
    result_with_module <-
      igraph::vertex_attr(graph_data) %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(p.adjust = as.numeric(p.adjust)) %>%
      dplyr::arrange(module, p.adjust)

    if (database == "go") {
      result_with_module <-
        result_with_module %>%
        dplyr::arrange(ONTOLOGY, module, p.adjust)
    }

    module_result <-
      result_with_module %>%
      plyr::dlply(.variables = .(module)) %>%
      purrr::map(function(x) {
        # cat(unique(x$module), " ")
        if (nrow(x) == 1) {
          return(x)
        }

        x =
          x %>%
          dplyr::arrange(p.adjust)

        x$node <-
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

        x =
          x %>%
          dplyr::select(module, everything()) %>%
          dplyr::distinct(module, .keep_all = TRUE)

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


    ####output some results
    dir.create(
      file.path(path, "Similarity_plot"),
      recursive = TRUE,
      showWarnings = FALSE
    )

    message("Output module plot...")

    for (temp_cluster in unique(result_with_module$module)) {
      cat(temp_cluster, " ")

      if (sum(result_with_module$module == temp_cluster) == 1) {
        next()
      }

      plot1 <-
        graph_data %>%
        tidygraph::filter(module == temp_cluster) %>%
        ggraph(layout = 'fr',
               circular = FALSE) +
        geom_edge_link(
          aes(width = sim),
          color = "black",
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_point(
          aes(fill = -log(p.adjust, 10),
              size = Count),
          shape = 21,
          alpha = 1,
          show.legend = TRUE
        ) +
        geom_node_text(aes(x = x,
                           y = y,
                           label = Description),
                       check_overlap = TRUE) +
        guides(fill = guide_legend(ncol = 1)) +
        scale_edge_width_continuous(range = c(0.1, 2)) +
        scale_size_continuous(range = c(3, 10)) +
        ggraph::theme_graph() +
        theme(
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          legend.position = "left",
          legend.background = element_rect(fill = "transparent", color = NA)
        )

      plot2 <-
        result_with_module %>%
        dplyr::filter(module == temp_cluster) %>%
        dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::mutate(Description = factor(Description, levels = Description)) %>%
        ggplot(aes(p.adjust, Description)) +
        geom_bar(stat = "identity", fill = "black") +
        geom_text(
          aes(x = 0, Description, label = Description),
          hjust = 0,
          size = 5,
          color = "red"
        ) +
        theme_bw() +
        labs(y = "", x = "-log10(FDR adjusted P value)") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
        theme(
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )

      temp_data =
        result_with_module %>%
        dplyr::filter(module == temp_cluster) %>%
        dplyr::mutate(p.adjust = -log(as.numeric(p.adjust, 10))) %>%
        dplyr::select(Description, p.adjust) %>%
        dplyr::mutate(Description = stringr::str_replace_all(Description, ",", "")) %>%
        plyr::dlply(.variables = .(Description)) %>%
        purrr::map(function(x) {
          data.frame(
            word = stringr::str_split(x$Description, " ")[[1]],
            p.adjust = x$p.adjust
          )
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        plyr::dlply(.variables = .(word)) %>%
        purrr::map(function(x) {
          x$p.adjust <- sum(x$p.adjust)
          x %>%
            dplyr::distinct(word, .keep_all = TRUE)
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        dplyr::filter(!word %in% remove_words)

      plot3 =
        temp_data %>%
        ggplot(aes(label = word, size = p.adjust)) +
        ggwordcloud::geom_text_wordcloud() +
        scale_radius(range = c(5, 15), limits = c(0, NA)) +
        theme_minimal()

      plot <-
        plot1 + plot2 + plot3 + patchwork::plot_layout(nrow = 1)

      ggsave(
        plot,
        filename = file.path(
          path,
          "Similarity_plot",
          paste(temp_cluster, "sim_plot.pdf", sep = "_")
        ),
        width = 21,
        height = 7
      )
    }


    # ##matrix tow show the cluster GO terms
    # result_with_module %>%
    #   dplyr::group_by(ONTOLOGY) %>%
    #   dplyr::summarise(n = n()) %>%
    #   dplyr::mutate(n = n * 10 / max(n) + 2)
    #output the correlation matrix

    # if(database == "go"){
    #   message("Output correlation matrix plot...")
    #   for(ont in c('MF', "BP", "CC")) {
    #     cat(ont, " ")
    #     show_matrix_cluster(
    #       result = result_with_module %>% dplyr::mutate(Direction = "UP"),
    #       ont = ont,
    #       measure = "Wang",
    #       remove_words = remove_words,
    #       margin = 15,
    #       width = 14,
    #       height = 8,
    #       path = path,
    #       top = 15
    #     )
    #   }
    # }

    ###output the cluster annotation for each cluster
    if (database == "go") {
      dir.create(file.path(path, "GO_module_graph"), showWarnings = FALSE)

      unique(module_result$module) %>%
        purrr::map(
          .f = function(x) {
            cat(x, " ")
            number <- module_result %>%
              dplyr::filter(module == x) %>%
              pull(module_content_number) %>%
              as.numeric()
            if (number == 1) {
              return(NULL)
            }

            temp_id <-
              module_result %>%
              dplyr::filter(module == x) %>%
              dplyr::pull(node) %>%
              stringr::str_split(";") %>%
              `[[`(1) %>%
              pRoloc::goIdToTerm(keepNA = FALSE) %>%
              data.frame(id = ., class = "YES") %>%
              tibble::rownames_to_column(var = "name")

            temp_plot =
              GOSim::getGOGraph(term = temp_id$name, prune = Inf) %>%
              igraph::igraph.from.graphNEL() %>%
              tidygraph::as_tbl_graph() %>%
              tidygraph::left_join(temp_id, by = "name") %>%
              dplyr::mutate(class = case_when(is.na(class) ~ "NO",
                                              TRUE ~ class))

            plot =
              temp_plot %>%
              ggraph(layout = 'kk',
                     circular = FALSE) +
              geom_edge_link(
                color = ggsci::pal_aaas()(n = 10)[1],
                alpha = 1,
                arrow = grid::arrow(
                  angle = 10,
                  length = unit(0.2, "inches"),
                  type = "closed"
                ),
                show.legend = FALSE
              ) +
              geom_node_point(
                aes(fill = class),
                shape = 21,
                alpha = 1,
                size = 6,
                show.legend = FALSE
              ) +
              geom_node_text(aes(
                x = x,
                y = y,
                label = ifelse(class == "YES", id, NA)
              ),
              size = 3,
              repel = TRUE) +
              scale_fill_manual(values = c('YES' = "red", 'NO' = "white")) +
              ggraph::theme_graph() +
              theme(
                plot.background = element_rect(fill = "transparent", color = NA),
                panel.background = element_rect(fill = "transparent", color = NA),
                legend.position = "left",
                legend.background = element_rect(fill = "transparent", color = NA)
              )
            # plot
            ggsave(
              plot,
              filename = file.path(
                path,
                "GO_module_graph",
                paste(x, "_GO graph.pdf", sep = "")
              ),
              width = 7,
              height = 7
            )
          }
        )
    }
    message("\nDone")

  }
