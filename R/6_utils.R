# msg <- function(..., startup = FALSE) {
#   if (startup) {
#     if (!isTRUE(getOption("mapa.quiet"))) {
#       packageStartupMessage(text_col(...))
#     }
#   } else {
#     message(text_col(...))
#   }
# }

# text_col <- function(x) {
#   # If RStudio not available, messages already printed in black
#   if (!rstudioapi::isAvailable()) {
#     return(x)
#   }
#
#   if (!rstudioapi::hasFun("getThemeInfo")) {
#     return(x)
#   }
#
#   theme <- rstudioapi::getThemeInfo()
#
#   if (isTRUE(theme$dark))
#     crayon::white(x)
#   else
#     crayon::black(x)
#
# }

#' List all packages in the mapa
#'
#' @param include_self Include mapa in the list?
#' @export
#' @return mapa packages
#' @examples
#' mapa_packages()
mapa_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("mapa")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <-
    vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))

  if (include_self) {
    names <- c(names, "mapa")
  }

  names
}

invert <- function(x) {
  if (length(x) == 0)
    return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(paste0(...),
                crayon::make_style(grDevices::grey(level), grey = TRUE))
}



database_color =
  c(
    GO = "#1F77B4FF",
    KEGG = "#FF7F0EFF",
    Reactome = "#2CA02CFF"
  )

remove_words <-
  c(
    "to",
    "of",
    "in",
    "type",
    "pathway",
    "IX",
    "part",
    "positive",
    "negative",
    "life",
    "control",
    "quality",
    "body",
    "late",
    "cell",
    "species",
    "cells",
    "or",
    "levels",
    "as",
    "on",
    "by",
    "small",
    "other",
    "involved",
    "alpha",
    "specific",
    "number",
    "through",
    "outer",
    "large",
    "rough",
    "early",
    "via",
    "smooth",
    "system",
    "into",
    "entry",
    "and",
    "T",
    "based",
    "within",
    "from",
    "built",
    "mediated",
    "-",
    "_",
    "animal",
    "the",
    "free",
    "a",
    "pool",
    "60S",
    "40S",
    "and",
    "chain",
    "Decay",
    "enhanced",
    "independent",
    "joining",
    "4",
    "2",
    "up",
    "take",
    "release",
    'Like',
    "presentation",
    "Class",
    "I",
    "mediated",
    "exchange",
    "&",
    "events",
    "B",
    "an",
    "",
    "at",
    "B",
    "Base",
    "c",
    "E",
    "during",
    "for",
    "Major",
    "NOTCH",
    "Of",
    "Opening",
    "Pathway",
    "processing",
    "free",
    letters,
    LETTERS,
    "family",
    "them",
    "ii",
    "class",
    1:7,
    "group",
    "phase",
    "ar",
    "orc",
    "new",
    "ap",
    "ends",
    "sars-cov-2",
    "upon",
    "ix",
    "major",
    "System",
    "with",
    "affected",
    "along",
    "AP",
    "AR",
    "associated",
    "Associated",
    "association",
    "containing",
    "down",
    "Ends",
    "II",
    "SARS-CoV",
    "ARS-CoV-2-host",
    "positively",
    "Network",
    "virus",
    "regulation",
    "Processing",
    "protein",
    "Protein",
    "Nonsense",
    "Nonsense-Mediated",
    "production",
    "pathways",
    "multiple",
    "scanning",
    "site",
    "The",
    "start",
    "pattern",
    "Processing",
    "Phase",
    "Packaging",
    "human",
    "Human",
    "gene",
    "genome",
    "foam",
    "classical",
    "beta",
    "2A",
    "11",
    "17",
    1:100,
    "5'-3'",
    "A-I",
    "absence",
    "break",
    "end",
    "second",
    "zone",
    "activity",
    "binding",
    "response",
    "receptor",
    "signaling",
    "process"
  )



#' get_jaccard_index_for_three_databases Function
#'
#' This function computes the Jaccard index, a measure of similarity for the genes listed in the results of three different databases (GO, KEGG, Reactome).
#' It first extracts the relevant gene information from the provided object, then computes the Jaccard index between all pairs of modules across the three databases.
#'
#' @param variable_info An variable_info containing gene information.
#' @param module_result_go A data frame containing module results from GO database.
#' @param module_result_kegg A data frame containing module results from KEGG database.
#' @param module_result_reactome A data frame containing module results from Reactome database.
#' @param analysis_type Character, type of analysis to perform: either `"enrich_pathway"` or `"do_gsea"`.
#'
#' @return A data frame with the Jaccard index values between all pairs of modules across the three databases.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @export

get_jaccard_index_for_three_databases <-
  function(variable_info,
           module_result_go,
           module_result_kegg,
           module_result_reactome,
           analysis_type = c("enrich_pathway", "do_gsea")) {
    check_variable_info(variable_info)

    analysis_type <- match.arg(analysis_type)

    met_data <-
      rbind(module_result_go,
            module_result_kegg,
            module_result_reactome)

    if (is.null(met_data)) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        value = numeric()
      ))
    }

    if (nrow(met_data) == 0 | nrow(met_data) == 1) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        value = numeric()
      ))
    }

    if (analysis_type == "do_gsea") {
      met_data <-
        met_data %>%
        dplyr::rename(geneID = core_enrichment) %>%
        dplyr::mutate(geneID = stringr::str_replace(geneID, ";", "/"))
    }

    temp_data <-
      met_data$geneID %>%
      stringr::str_split("/") %>%
      unique() %>%
      purrr::map(function(x) {
        if (stringr::str_detect(x[1], "ENSG")) {
          return(x)
        }

        if (stringr::str_detect(x[1], "[A-Za-z]")) {
          return(variable_info$ensembl[match(x, variable_info$uniprot)])
        }

        return(variable_info$ensembl[match(x, variable_info$entrezid)])

      })

    names(temp_data) =
      met_data$module

    ##calculate jaccard index
    jaccard_index =
      purrr::map(
        1:(length(temp_data) - 1),
        .f = function(idx) {
          purrr::map(
            temp_data[(idx + 1):length(temp_data)],
            .f = function(y) {
              length(intersect(temp_data[[idx]], y)) / length(union(temp_data[[idx]], y))
            }
          ) %>%
            unlist() %>%
            data.frame(value = .) %>%
            dplyr::mutate(name2 = names(temp_data)[(idx + 1):length(temp_data)]) %>%
            # tibble::rownames_to_column(var = "name2") %>%
            data.frame(name1 = names(temp_data)[idx], .) %>%
            dplyr::select(name1, name2, value)
        }
      ) %>%
      do.call(rbind, .) %>%
      as.data.frame()
    return(jaccard_index)
  }



arrange_coords <- function(coords, ratio = 0.95) {
  coords <-
    coords %>%
    plyr::dlply(.variables = .(class)) %>%
    purrr::map(function(x) {
      x <-
        x %>%
        dplyr::arrange(x)
      if (length(unique(coords$class)) == 4) {
        if (x$class[1] == "Molecule") {
          min_times <- 1
          max_times <- 1
        }

        if (x$class[1] == "Pathway") {
          min_times <- 1.1
          max_times <- ratio
        }

        if (x$class[1] == "Module") {
          min_times <- 1.1 ^ 2
          max_times <- ratio ^ 2
        }

        if (x$class[1] == "Functional_module") {
          min_times <- 1.1 ^ 3
          max_times <- ratio ^ 3
        }

      }

      if (length(unique(coords$class)) == 3) {
        if (x$class[1] == "Pathway") {
          min_times <- 1
          max_times <- 1
        }

        if (x$class[1] == "Module") {
          min_times <- 1.1
          max_times <- ratio
        }

        if (x$class[1] == "Functional_module") {
          min_times <- 1.1 ^ 2
          max_times <- ratio ^ 2
        }

      }

      x$x <-
        seq(
          from = max(coords$x) - max_times * max(coords$x),
          to = max_times * max(coords$x),
          length.out = nrow(x)
        )
      x
    }) %>%
    do.call(rbind, .) %>%
    as.data.frame()

  coords <-
    coords %>%
    dplyr::arrange(index)
  coords
}


#' Check Variable Information
#'
#' This function checks if the input data frame `variable_info` has the required columns,
#' and if these columns have any non-NA values.
#'
#' @param variable_info A data frame. The data frame should have at least the columns
#'   "ensembl", "symbol", "uniprot", and "entrezid". Each of these columns should have
#'   at least one non-NA value.
#'
#' @param order_by A character specifying the column to order by.
#' @return This function does not return any value. It stops execution and throws an
#'   error if any of the required conditions are not met.
#'
#' @export
#'
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#'
#' @examples
#' variable_info <- data.frame(
#'   ensembl = c("ENSG000001", NA, NA),
#'   symbol = c("Gene1", "Gene2", "Gene3"),
#'   uniprot = c(NA, "P12345", "Q67890"),
#'   entrezid = c(101, 102, 103)
#' )
#' check_variable_info(variable_info)

check_variable_info <-
  function(variable_info, order_by = NULL) {
    ###check order_by
    if (!is.null(order_by)) {
      if (!is.character(order_by)) {
        stop("order_by should be character")
      }

      if (all(!order_by %in% colnames(variable_info))) {
        stop("order_by should be in the variable_info")
      }

      ####the order by column should be numeric and should no any NA values
      if (any(is.na(variable_info[, order_by, drop = TRUE]))) {
        stop("order_by column should not have any NA values")
      }

    }

    if (all(colnames(variable_info) != "ensembl")) {
      stop("ensembl should be in the variable_info")
    } else{
      if (all(is.na(variable_info$ensembl))) {
        stop("All ensembl column are NA")
      }
    }

    if (all(colnames(variable_info) != "symbol")) {
      stop("symbol should be in the variable_info")
    } else{
      if (all(is.na(variable_info$symbol))) {
        stop("All symbol column are NA")
      }
    }

    if (all(colnames(variable_info) != "uniprot")) {
      stop("uniprot should be in the variable_info")
    } else{
      if (all(is.na(variable_info$uniprot))) {
        stop("All uniprot column are NA")
      }
    }

    if (all(colnames(variable_info) != "entrezid")) {
      stop("entrezid should be in the variable_info")
    } else{
      if (all(is.na(variable_info$entrezid))) {
        stop("All entrezid column are NA")
      }
    }
  }




#' Pathway Enrichment Analysis and GO Similarity Measurement
#'
#' This function allows for the execution of pathway enrichment analysis and
#' measurement of semantic similarity among Gene Ontology (GO) terms.
#'
#' @param result A data frame containing GO term IDs and ontologies.
#' @param sim.cutoff A numeric value for the similarity cutoff (default: 0).
#' @param measure.method A character vector specifying the semantic similarity
#' measure methods for GO database
#'   (default: c("Sim_XGraSM_2013", "Sim_Resnik_1999", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_Lin_1998", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Wang_2007","Sim_GOGO_2018")).
#' @return A data frame containing pairs of GO terms and their similarity values.
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @examples
#' \dontrun{
#' # Assuming `result` is your data frame containing GO term IDs and ontologies.
#' similarity_matrix <- get_go_result_sim(result = result)
#' }
#' @export

get_go_result_sim <-
  function(result,
           sim.cutoff = 0,
           measure.method = c("Sim_XGraSM_2013", "Sim_Resnik_1999", "Sim_Relevance_2006", "Sim_SimIC_2010", "Sim_Lin_1998", "Sim_EISI_2015", "Sim_AIC_2014", "Sim_Wang_2007","Sim_GOGO_2018")) {
    measure.method <-
      match.arg(measure.method)

    if (is.null(result)) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        sim = numeric()
      ))
    }

    if (nrow(result) == 0) {
      return(data.frame(
        name1 = character(),
        name2 = character(),
        sim = numeric()
      ))
    }

    bp_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "BP"],
                                        ont = "BP",
                                        measure = measure.method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)

    name <- apply(bp_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })

    bp_sim_matrix <-
      bp_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)

    mf_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "MF"],
                                        ont = "MF",
                                        measure = measure.method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)

    name <- apply(mf_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })

    mf_sim_matrix <-
      mf_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)

    # cc_sim_matrix <-
    #   simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "CC"],
    #                                     ont = "CC",
    #                                     measure = measure.method) %>%
    #   as.data.frame() %>%
    #   tibble::rownames_to_column(var = "name1") %>%
    #   tidyr::pivot_longer(cols = -name1,
    #                       names_to = "name2",
    #                       values_to = "sim") %>%
    #   dplyr::filter(name1 != name2) %>%
    #   dplyr::filter(sim > sim.cutoff)
#
#     name <- apply(cc_sim_matrix, 1, function(x) {
#       paste(sort(x[1:2]), collapse = "_")
#     })
#
#     cc_sim_matrix <-
#       cc_sim_matrix %>%
#       dplyr::mutate(name = name) %>%
#       dplyr::arrange(name) %>%
#       dplyr::distinct(name, .keep_all = TRUE) %>%
#       dplyr::select(-name)

    sim_matrix <-
      rbind(bp_sim_matrix, mf_sim_matrix) %>%
      as.data.frame()
    sim_matrix
  }


#' Similarity calculation between KEGG pathways
#'
#' @param term_id Character, KEGG pathway IDs
#' @param measure.method Character, method for calculating the semantic similarity
#' for KEGG terms
#'   (default: c("jaccard", "dice", "overlap", "kappa")).
#'
#' @return A symmetric matrix
#'
#' Add citation
#'
#' @examples
term_similarity_KEGG <- function(term_id, measure.method) {

  measure.method <- match.arg(measure.method)

  species <- gsub("^([a-zA-Z]+)(\\d+$)", "\\1", term_id[1])

  kegg_data <- tryCatch(
    expr = {
      getFromNamespace("prepare_KEGG", "clusterProfiler")(species, "KEGG", "kegg")
      },
    error = function(e) {
      getFromNamespace("get_data_from_KEGG_db", "clusterProfiler")(species)
      })

  gl <- kegg_data$PATHID2EXTID[term_id]

  term_similarity_internal(gl = gl,
                           measure.method = measure.method)
}


#' Similarity calculation between Reactome terms
#'
#' @param term_id Character, Reactome term IDs
#' @param measure.method Character, method for calculating the semantic similarity
#' for Reactome terms
#'   (default: c("jaccard", "dice", "overlap", "kappa")).
#'
#' @return A symmetric matrix
#'
#' Add citation
#'
#' @examples
term_similarity_Reactome <- function(term_id, measure.method) {

  measure.method <- match.arg(measure.method)

  all <- as.list(reactome.db::reactomePATHID2EXTID)
  gl <- all[term_id]

  term_similarity_internal(gl = gl,
                           measure.method = measure.method)
}


#' Similarity calculation between pathways.
#'
#' This function allows for the execution of measurement of semantic similarity
#' among KEGG and Reactome terms based on the overlap of genes.
#'
#' @param gl Named list, genes that are in the enriched pathways.
#' @param measure.method Character, method for calculating the semantic similarity
#' for KEGG and Reactome terms
#' @param remove_negative Logical, if TRUE reset the negative similarity values to zero
#'
#' @return A symmetric matrix.
#'
#' Add citation
#'
#'
#' @examples
term_similarity_internal <-
  function(gl,
           measure.method = c("jaccard", "dice", "overlap", "kappa"),
           remove_negative = TRUE) {

    measure.method <- match.arg(measure.method)

    all <- unique(unlist(gl))
    gl <- lapply(gl, function(x) as.numeric(factor(x, levels = all)))
    n <- length(gl)

    pathway_gene_m <- matrix(0, ncol = length(all), nrow = n)
    for(i in seq_len(n)) {
      pathway_gene_m[i, gl[[i]]] = 1
      }
    pathway_gene_m <- as(pathway_gene_m, "sparseMatrix")

    if(measure.method == "kappa") {
      mat <- kappa_dist(pathway_gene_m, remove_negative = remove_negative)
      } else if(measure.method == "overlap") {
        mat <- overlap_dist(pathway_gene_m)
        } else {
          mat <- proxyC::simil(pathway_gene_m, method = method)
          }

    sim_matrix <-  as.matrix(mat)
    diag(sim_matrix) <- 1
    rownames(sim_matrix) = colnames(sim_matrix) = names(gl)

    return(sim_matrix)
}

kappa_dist <- function(m, remove_negative = TRUE) {
  tab <- ncol(m)
  po <- proxyC::simil(m, method = "simple matching")
  m_yes <- Matrix::rowSums(m)
  m_no <- abs(Matrix::rowSums(m - 1))
  pe <- (outer(m_yes, m_yes, FUN = "*") + outer(m_no, m_no, FUN = "*"))/tab^2
  k <- (po - pe)/(1 - pe)
  if(remove_negative) k[k < 0] <- 0
  return(k)
}

overlap_dist <- function(m) {
  n = Matrix::rowSums(m)
  proxyC::simil(m, method = "dice")*outer(n, n, FUN = "+")/2/outer(n, n, pmin)
}
