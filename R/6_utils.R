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

remove_words <- c(
  # === EMPTY STRINGS AND SYMBOLS ===
  "",
  "-",
  "_",
  "&",

  # === NUMBERS AND RANGES ===
  as.character(1:100),

  # === ALPHANUMERIC CODES ===
  "2A",
  # "40S",
  # "60S",
  "5'-3'",
  "A-I",

  # === SINGLE LETTERS ===
  LETTERS,
  letters,
  # === ROMAN NUMERALS ===
  "i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x",
  "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
  "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii", "xix", "xx",
  "XI", "XII", "XIII", "XIV", "XV", "XVI", "XVII", "XVIII", "XIX", "XX",

  # === PREPOSITIONS AND CONJUNCTIONS ===
  "and",
  "as",
  "at",
  "by",
  "during",
  "for",
  "from",
  "in",
  "into",
  "of",
  "on",
  "or",
  "through",
  "to",
  "upon",
  "via",
  "with", "within",
  "between",

  # === ARTICLES AND DETERMINERS ===
  "a",
  "an",
  "the", "The",
  "This", "this", "them",
  "That", "that", "which",
  "such", "other", "both",
  "it", "their", "Their",
  "study", "studies", "article",

  # === DESCRIPTIVE ADJECTIVES ===
  "early", "late",
  "small", "large",
  "major", "Major",
  "multiple",
  # "negative",
  "new",
  "outer",
  # "positive", "positively",
  "rough", "smooth",
  "specific",
  "free",
  "overall", "additionally",
  "particularly",
  "important", "importantly", "critical",

  # === DIRECTIONAL/POSITIONAL WORDS ===
  "along",
  "down",
  "ends", "Ends",
  "up",

  # === BIOLOGICAL TERMS - GENERAL ===
  "animal",
  "body",
  "cell", "cells",
  "gene", "genome",
  "human", "Human",
  "life",
  "species",

  # === BIOLOGICAL TERMS - MOLECULAR ===
  "alpha",
  "beta",
  "chain",
  "protein", "Protein",
  # "receptor",

  # === BIOLOGICAL PROCESSES ===
  "activity",
  # "binding",
  "process", "processing", "Processing",
  # "production",
  # "regulation",
  # "response",
  # "signaling",

  # === BIOLOGICAL SYSTEMS/STRUCTURES ===
  "pathway", "Pathway", "pathways",
  "system", "System", "module",

  # === CLASSIFICATION TERMS ===
  "class", "Class", "classical",
  "family", "group", "type", "role",

  # === BIOLOGICAL STATES/CONDITIONS ===
  "affected",
  "Associated", "associated", "association",
  "base", "build",
  "contain",
  "enhance",
  "independent",
  "involve", "include",
  "mediate", "control",

  # === PHASES/STAGES ===
  "Phase", "phase",
  "start",
  "end",
  "second",

  # === BIOLOGICAL CONCEPTS ===
  "absence",
  "break",
  "Decay",
  "entry",
  "events",
  "exchange",
  "joining",
  "levels",
  "Network", "network",
  "number",
  "Opening",
  "Packaging",
  "part",
  "pattern",
  "pool",
  "presentation",
  "quality",
  "release",
  "scanning",
  "site", "zone",
  "take",

  # === SPECIFIC BIOLOGICAL ENTITIES ===
  "AP", "ap",
  "AR", "ar",
  # "ARS-CoV-2-host",
  "Base", "foam",
  "Nonsense",
  "Nonsense-Mediated",
  "NOTCH",
  "orc",
  # "SARS-CoV",
  # "sars-cov-2",
  # "virus"

  # === non-biological VERBs ===
  "be", "is", "are", "have",
  "can", "reveal", "suggest", "understand", "underscore",
  "state", "show", "highlight", "inform", "illustrate"
)



#' get_jaccard_index_for_diff_databases Function
#'
#' This function computes the Jaccard index, a measure of similarity between gene or metabolite sets
#' from different pathway databases. For genes, it can compare GO, KEGG, and Reactome databases.
#' For metabolites, it can compare HMDB and KEGG databases.
#'
#' @param query_type Character, the category of biological entity to query. Must be either "gene" or "metabolite".
#' @param variable_info An object containing mapping information between different identifiers
#'        (for genes: uniprot, entrezid, ensembl; for metabolites: keggid, hmdbid).
#' @param module_result_go A data frame containing module results from GO term enrichment analysis.
#'        Required if query_type is "gene". Default is NULL.
#' @param module_result_kegg A data frame containing module results from KEGG pathway enrichment analysis
#'        for genes. Required if query_type is "gene". Default is NULL.
#' @param module_result_reactome A data frame containing module results from Reactome pathway enrichment analysis.
#'        Required if query_type is "gene". Default is NULL.
#' @param module_result_hmdb A data frame containing module results from HMDB pathway enrichment analysis.
#'        Required if query_type is "metabolite". Default is NULL.
#' @param module_result_metkegg A data frame containing module results from KEGG pathway enrichment analysis
#'        for metabolites. Required if query_type is "metabolite". Default is NULL.
#' @param analysis_type Character, type of analysis to perform: either "enrich_pathway" or "do_gsea".
#'        Default is "enrich_pathway".
#'
#' @return A data frame with three columns:
#'         \item{name1}{Character, name of the first module}
#'         \item{name2}{Character, name of the second module}
#'         \item{value}{Numeric, Jaccard index value between the two modules}
#'         Returns an empty data frame if insufficient data is provided.
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{Yifei.ge@outlook.com}
#'
#' @importFrom dplyr rename mutate select %>%
#' @importFrom stringr str_replace str_split str_detect
#' @importFrom purrr map
#' @export

get_jaccard_index_for_diff_databases <- function(
    query_type = c("gene", "metabolite"),
    variable_info,
    module_result_go = NULL,
    module_result_kegg = NULL,
    module_result_reactome = NULL,
    module_result_hmdb = NULL,
    module_result_metkegg = NULL,
    analysis_type = c("enrich_pathway", "do_gsea")) {

    analysis_type <- match.arg(analysis_type)

    if (missing(query_type)) {
      stop("query_type is required")
    }
    query_type <- match.arg(query_type)

    if (query_type == "gene") {
      check_variable_info(variable_info, query_type = "gene")
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
        #unique() %>%
        purrr::map(function(x) {
          if (stringr::str_detect(x[1], "ENS")) {
            return(x)
          }

          if (stringr::str_detect(x[1], "[A-Za-z]")) {
            return(variable_info$ensembl[match(x, variable_info$uniprot)])
          }

          return(variable_info$ensembl[match(x, variable_info$entrezid)])

        })

    } else if (query_type == "metabolite") {
      check_variable_info(variable_info, query_type = "metabolite")
      met_data <-
        rbind(module_result_hmdb,
              module_result_metkegg)
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

      temp_data <-
        met_data$mapped_id %>%
        stringr::str_split("/") %>%
        #unique() %>%
        purrr::map(function(x) {
          if (stringr::str_detect(x[1], "HMDB")) {
            return(x)
          }

          if (stringr::str_detect(x[1], "C")) {
            return(variable_info$hmdbid[match(x, variable_info$keggid)])
          }
        })

    }

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
#' This function validates the input data frame `variable_info` to ensure it contains
#' the required columns for the specified query type and that these columns have
#' valid (non-NA) values. It also validates the `order_by` parameter when provided.
#'
#' @param variable_info A data frame containing biological entity information.
#'   \itemize{
#'     \item For \code{query_type = "gene"}: Must contain an "ensembl" column
#'           with at least one non-NA value.
#'     \item For \code{query_type = "metabolite"}: Must contain a "keggid" column
#'           with at least one non-NA value.
#'   }
#' @param query_type Character string specifying the category of biological entity
#'   to query. Accepted values are "gene" or "metabolite".
#' @param order_by Optional character string specifying the column name to order by
#'   during Gene Set Enrichment Analysis (GSEA). If provided:
#'   \itemize{
#'     \item Must be a character string
#'     \item Must correspond to an existing column in \code{variable_info}
#'     \item The specified column must not contain any NA values
#'     \item The column should contain numeric values for proper ordering
#'   }
#'
#' @return This function does not return any value. It performs validation checks
#'   and stops execution with an informative error message if any required
#'   conditions are not met.
#'
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#'
#' @keywords internal

check_variable_info <-
  function(variable_info, query_type, order_by = NULL) {
    if (query_type == "gene") {
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

      # if (all(colnames(variable_info) != "symbol")) {
      #   stop("symbol should be in the variable_info")
      # } else{
      #   if (all(is.na(variable_info$symbol))) {
      #     stop("All symbol column are NA")
      #   }
      # }

      # if (all(colnames(variable_info) != "uniprot")) {
      #   stop("uniprot should be in the variable_info")
      # } else{
      #   if (all(is.na(variable_info$uniprot))) {
      #     stop("All uniprot column are NA")
      #   }
      # }

      # if (all(colnames(variable_info) != "entrezid")) {
      #   stop("entrezid should be in the variable_info")
      # } else{
      #   if (all(is.na(variable_info$entrezid))) {
      #     stop("All entrezid column are NA")
      #   }
      # }

    } else if (query_type == "metabolite") {
      # if (all(colnames(variable_info) != "hmdbid")) {
      #   stop("HMDB ID should be in the variable_info")
      # } else {
      #   if (all(is.na(variable_info$hmdbid))) {
      #     stop("All hmdbid column are NA")
      #   }
      # }

      if (all(colnames(variable_info) != "keggid")) {
        stop("KEGG ID should be in the variable_info")
      } else {
        if (all(is.na(variable_info$keggid))) {
          stop("All keggid column are NA")
        }
      }
    }
  }


#' Calculate similarity among Gene Ontology (GO) terms
#'
#' @description
#' Internal function that computes similarity between GO terms from enrichment results.
#' This function processes GO terms by ontology category (BP, MF, CC) and calculates pairwise similarities.
#'
#' @param result A data frame containing GO term IDs in column 'ID' and ontologies in column 'ONTOLOGY'.
#' @param go.orgdb An organism-specific database for GO annotations (default: NULL).
#' @param measure.method A character vector specifying the semantic similarity
#'                      measure method for GO terms. Default is `"Sim_XGraSM_2013"`.
#'                      See `simona::all_term_sim_methods()` for available measures.
#' @param control.method A list of parameters passed to the specified measure method for GO term semantic similarity.
#'                      For details, see https://jokergoo.github.io/simona/articles/v05_term_similarity.html.
#'
#' @return A data frame containing pairs of GO terms (name1, name2) and their similarity values (sim).
#'         Also includes an attribute "obsolete_terms" listing any obsolete GO terms encountered.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter
#'
#' @keywords internal
#'
#' @author Xiaotao Shen \email{shenxt1990@outlook.com}
#' @author Yifei Ge \email{yifeii.ge@outlook.com}

get_go_result_sim <-
  function(result,
           go.orgdb = NULL,
           measure.method = "Sim_XGraSM_2013",
           control.method = list()) {

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

    if (sum(result$ONTOLOGY == "BP") > 0) {
      bp_sim_matrix <-
        GO_similarity_internal(go_id = result$ID[result$ONTOLOGY == "BP"],
                               ont = "BP",
                               go.orgdb = go.orgdb,
                               measure = measure.method,
                               control.method = control.method)
      bp_obsolete_terms <- attr(bp_sim_matrix, "obsolete_terms")
      bp_sim_df <-
        bp_sim_matrix %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "name1") |>
        tidyr::pivot_longer(cols = -name1,
                            names_to = "name2",
                            values_to = "sim") |>
        dplyr::filter(name1 < name2)
        # dplyr::filter(sim > sim.cutoff)
      message("Completed GO term (BP) similarity calculation.")
    } else {
      bp_sim_df <- NULL
      bp_obsolete_terms <- NULL
    }

    # name <- apply(bp_sim_df, 1, function(x) {
    #   paste(sort(x[1:2]), collapse = "_")
    # })
    #
    # bp_sim_df <-
    #   bp_sim_df %>%
    #   dplyr::mutate(name = name) %>%
    #   dplyr::arrange(name) %>%
    #   dplyr::distinct(name, .keep_all = TRUE) %>%
    #   dplyr::select(-name)

    if (sum(result$ONTOLOGY == "MF") > 0) {
      mf_sim_matrix <-
        GO_similarity_internal(go_id = result$ID[result$ONTOLOGY == "MF"],
                               ont = "MF",
                               go.orgdb = go.orgdb,
                               measure = measure.method,
                               control.method = control.method)
      mf_obsolete_terms <- attr(mf_sim_matrix, "obsolete_terms")
      mf_sim_df <-
        mf_sim_matrix %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "name1") |>
        tidyr::pivot_longer(cols = -name1,
                            names_to = "name2",
                            values_to = "sim") |>
        dplyr::filter(name1 < name2)
        # dplyr::filter(sim > sim.cutoff)
      message("Completed GO term (MF) similarity calculation.")
    } else {
      mf_sim_df <- NULL
      mf_obsolete_terms <- NULL
    }

    # name <- apply(mf_sim_df, 1, function(x) {
    #   paste(sort(x[1:2]), collapse = "_")
    # })
    #
    # mf_sim_df <-
    #   mf_sim_df %>%
    #   dplyr::mutate(name = name) %>%
    #   dplyr::arrange(name) %>%
    #   dplyr::distinct(name, .keep_all = TRUE) %>%
    #   dplyr::select(-name)

    if (sum(result$ONTOLOGY == "CC") > 0) {
      cc_sim_matrix <-
        GO_similarity_internal(go_id = result$ID[result$ONTOLOGY == "CC"],
                               ont = "CC",
                               go.orgdb = go.orgdb,
                               measure = measure.method,
                               control.method = control.method)
      cc_obsolete_terms <- attr(cc_sim_matrix, "obsolete_terms")
      cc_sim_df <-
        cc_sim_matrix |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "name1") |>
        tidyr::pivot_longer(cols = -name1,
                            names_to = "name2",
                            values_to = "sim") |>
        dplyr::filter(name1 < name2)
        # dplyr::filter(sim > sim.cutoff)
      message("Completed GO term (CC) similarity calculation.")
    } else {
      cc_sim_df <- NULL
      cc_obsolete_terms <- NULL
    }

    # name <- apply(cc_sim_df, 1, function(x) {
    #   paste(sort(x[1:2]), collapse = "_")
    # })
    #
    # cc_sim_df <-
    #   cc_sim_df %>%
    #   dplyr::mutate(name = name) %>%
    #   dplyr::arrange(name) %>%
    #   dplyr::distinct(name, .keep_all = TRUE) %>%
    #   dplyr::select(-name)

    sim_matrix <-
      rbind(bp_sim_df, mf_sim_df, cc_sim_df)

    attr(sim_matrix, "obsolete_terms") <- c(bp_obsolete_terms, mf_obsolete_terms, cc_obsolete_terms)

    return(sim_matrix)
  }

#' Compute Similarity Between GO Terms for each subontology
#'
#' Computes pair-wise semantic similarity scores for a set of Gene Ontology
#' (GO) terms with one of the similarity measures provided by **simona**.
#' The underlying ontology DAG and information‐content (IC) values are
#' constructed from the specified *OrgDb* annotation package and **cached** in a
#' hidden environment so repeated calls with the same settings are fast.
#'
#' For a detailed explanation of the available similarity metrics and their
#' tunable parameters, see the simona vignette:
#' <https://jokergoo.github.io/simona/articles/v05_term_similarity.html>.
#'
#' @note This helper is adapted from `GO_similarity()` in the
#'   **simplifyEnrichment** package (v2.0.0) by Zuguang Gu.
#'
#' @references
#'   Gu Z, Huebschmann D (2021) *simplifyEnrichment: an R/Bioconductor package
#'   for Clustering and Visualizing Functional Enrichment Results.*
#'   *Genomics, Proteomics & Bioinformatics.*
#'   doi:10.1016/j.gpb.2022.04.008
#'
#'   Source code:
#'   <https://www.bioconductor.org/packages/release/bioc/src/contrib/simplifyEnrichment_2.0.0.tar.gz>
#'
#' @param go_id          Character vector of GO term identifiers
#' @param ont            Character string specifying the ontology namespace:
#'                       `"BP"`, `"MF"`, or `"CC"`.
#' @param go.orgdb       OrgDb object or character string naming the *OrgDb* annotation package
#'                       used to derive gene–GO mappings.
#' @param measure        Character string giving the similarity measure to use.
#'                       Default is `"Sim_XGraSM_2013"`.  See
#'                       `simona::all_term_sim_methods()` for the full list.
#' @param control.method Named list of additional arguments passed to the
#'                       selected `measure`.  See the simona documentation for
#'                       details.
#'
#' @return A numeric matrix of pair-wise similarity scores with three
#'   attributes:
#'   \describe{
#'     \item{`"measure"`}{the similarity method used}
#'     \item{`"ontology"`}{the ontology namespace (`"GO:BP"`, `"GO:MF"`, or
#'                         `"GO:CC"`)}
#'     \item{`"obsolete_terms"`}{character vector of input GO terms that were
#'                               absent from the database and therefore
#'                               omitted}
#'   }
#'
#' @keywords internal
#' @noRd

env <- new.env()

GO_similarity_internal = function(go_id,
                                  ont = NULL,
                                  go.orgdb = NULL,
                                  measure = "Sim_XGraSM_2013",
                                  control.method = list()) {

  hash <- digest::digest(list(ont = ont, db = go.orgdb))
  if(is.null(env$go[[hash]])) {
    dag <- simona::create_ontology_DAG_from_GO_db(namespace = ont, org_db = go.orgdb, relations = c("part_of", "regulates"))

    ic <- simona::term_IC(dag, method = "IC_annotation")
    all_go_id <- names(ic[!is.na(ic)])

    env$go[[hash]] <- list(dag = dag, all_go_id = all_go_id)
  } else {
    dag <- env$go[[hash]]$dag
    all_go_id <- env$go[[hash]]$all_go_id
  }

  go_removed <- setdiff(go_id, all_go_id)
  if(length(go_removed)) {
    message(paste0(length(go_removed), "/", length(go_id), " GO term", ifelse(length(go_removed) == 1, ' is', 's are'), " removed."))
  }

  go_id <- intersect(go_id, all_go_id)
  go_sim <- simona::term_sim(dag,
                             go_id,
                             method = measure,
                             control = control.method)

  attr(go_sim, "measure") <- measure
  attr(go_sim, "ontology") <- paste0("GO:", ont)
  attr(go_sim, "obsolete_terms") <- go_removed

  return(go_sim)
}


#' Similarity calculation between KEGG pathways
#'
#' @note This function was adapted from term_similarity__from_KEGG() in the simplifyEnrichment package (version 1.14.0) by Zuguang Gu.
#'
#' @references Gu Z, Huebschmann D (2021). “simplifyEnrichment: an R/Bioconductor package for Clustering and Visualizing Functional Enrichment Results.” Genomics, Proteomics & Bioinformatics. doi:10.1016/j.gpb.2022.04.008.
#'
#' @source Source code download link: https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/simplifyEnrichment/simplifyEnrichment_1.14.0.tar.gz
#'
#' @param term_id Character, KEGG pathway IDs
#' @param measure.method Character, method for calculating the similarity
#' for KEGG terms, Choices are "jaccard", "dice", "overlap" and "kappa". Default is "jaccard".
#'
#' @return A symmetric matrix
#'
#' @keywords internal
term_similarity_KEGG <- function(term_id,
                                 measure.method = c("jaccard", "dice", "overlap", "kappa")) {

  measure.method <- match.arg(measure.method)

  species <- gsub("^([a-zA-Z]+)(\\d+$)", "\\1", term_id[1])

  kegg_data <- tryCatch(
    expr = {
      #A function that downloads KEGG data for a given species, KEGG type, and key type
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
#' @note This function was adapted from term_similarity_from_Reactome() in the simplifyEnrichment package (version 1.14.0) by Zuguang Gu.
#'
#' @references Gu Z, Huebschmann D (2021). “simplifyEnrichment: an R/Bioconductor package for Clustering and Visualizing Functional Enrichment Results.” Genomics, Proteomics & Bioinformatics. doi:10.1016/j.gpb.2022.04.008.
#'
#' @source Source code download link: https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/simplifyEnrichment/simplifyEnrichment_1.14.0.tar.gz
#'
#' @param term_id Character, Reactome term IDs
#' @param measure.method Character, method for calculating the similarity
#' for Reactome terms, Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#'
#' @return A symmetric matrix
#'
#' @keywords internal
term_similarity_Reactome <- function(term_id,
                                     measure.method = c("jaccard", "dice", "overlap", "kappa")) {

  measure.method <- match.arg(measure.method)

  all <- as.list(reactome.db::reactomePATHID2EXTID)
  gl <- all[term_id]

  term_similarity_internal(gl = gl,
                           measure.method = measure.method)
}

#' Similarity calculation between enriched HMDB or KEGG pathways from metabolite enrichment analysis
#'
#' @note This function was adapted from the simplifyEnrichment package (version 1.14.0) by Zuguang Gu.
#'
#' @references Gu Z, Huebschmann D (2021). “simplifyEnrichment: an R/Bioconductor package for Clustering and Visualizing Functional Enrichment Results.” Genomics, Proteomics & Bioinformatics. doi:10.1016/j.gpb.2022.04.008.
#'
#' @source Source code download link: https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/simplifyEnrichment/simplifyEnrichment_1.14.0.tar.gz
#'
#' @param enrichment_result Dataframe, result from metabolite enrichment analysis
#' @param measure.method Character, method for calculating the similarity
#' between enriched pathways, Choices are "jaccard", "dice", "overlap", "kappa". Default is "jaccard".
#'
#' @return A symmetric matrix
#'
#' @keywords internal

term_similarity_metabolite <- function(enrichment_result,
                                       measure.method = c("jaccard", "dice", "overlap", "kappa")) {
  measure.method <- match.arg(measure.method)

  ## Get annotated gene list for pathways
  gene_list <- split(enrichment_result$all_id, enrichment_result$pathway_id)
  for (i in 1:length(gene_list)) {
    metabolite_ids <- stringr::str_split(gene_list[[i]], pattern = ";")[[1]]
    gene_list[[i]] <- metabolite_ids
  }

  ## Calculate similarity
  term_similarity_internal(gl = gene_list,
                           measure.method = measure.method)
}

#' Similarity calculation between pathways.
#'
#' @description
#' This function allows for the execution of measurement of semantic similarity
#' among KEGG and Reactome terms based on the overlap of genes.
#'
#' @note This function was adapted from term_similarity() in the simplifyEnrichment package (version 1.14.0) by Zuguang Gu.
#'
#' @references Gu Z, Huebschmann D (2021). “simplifyEnrichment: an R/Bioconductor package for Clustering and Visualizing Functional Enrichment Results.” Genomics, Proteomics & Bioinformatics. doi:10.1016/j.gpb.2022.04.008.
#'
#' @source Source code download link: https://bioconductor.org/packages/3.19/bioc/src/contrib/Archive/simplifyEnrichment/simplifyEnrichment_1.14.0.tar.gz
#'
#' @param gl Named list, genes that are in the enriched pathways.
#' @param measure.method Character, method for calculating the semantic similarity
#' for KEGG and Reactome terms
#' @param remove_negative Logical, if TRUE reset the negative similarity values to zero
#'
#' @return A symmetric matrix.
#'
#' @keywords internal

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
          mat <- proxyC::simil(pathway_gene_m, method = measure.method)
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


#' Retrieve Metabolic Pathway Databases
#'
#' This function retrieves metabolic pathway data from the HMDB and KEGG databases.
#'
#' @param use_internal_data Logical. Indicates whether to use internal (local) data for the KEGG pathway retrieval.
#'   Note that the HMDB pathway retrieval currently does not support online access.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{hmdb_pathway}{A list of HMDB pathways.}
#'   \item{kegg_pathway}{A list of KEGG pathways.}
#' }
#'
#'
#' @examples
#' \dontrun{
#'   # Retrieve pathway databases using internal data for KEGG
#'   pathways <- get_met_pathway_db(use_internal_data = TRUE)
#'
#'   # Access HMDB pathways
#'   hmdb <- pathways$hmdb_pathway
#'
#'   # Access KEGG pathways
#'   kegg <- pathways$kegg_pathway
#' }
#'
#' @keywords internal
get_met_pathway_db <- function(use_internal_data = TRUE) {
  # Get hmdb pathway
  message("Get HMDB pathway ...")
  hmdb_pathway <-
    get_hmdb_pathways(
      ## TO DO: Get updated information.
      #local = use_internal_data,  # Currently this function unable to get online database
      threads = 3
    )

  # Get KEGG pathway
  message("Get KEGG pathway ...")
  kegg_pathway <-
    get_kegg_pathways(
      local = use_internal_data,
      organism = "hsa"
    )

  return(list(
    "hmdb_pathway" = hmdb_pathway,
    "kegg_pathway" = kegg_pathway
  ))
}

#' Get KEGG pathway database
#'
#' Retrieve pathway information from KEGG for a specified organism.
#' This is an internal function used by get_met_pathway_db().
#'
#' @param local Logical. If TRUE (default), uses locally stored data. If FALSE,
#'        retrieves data from KEGG API.
#' @param organism Character. KEGG organism code (default: "hsa" for Homo sapiens)
#' @return A pathway_database class object containing KEGG pathway information
#'
#' @importFrom KEGGREST keggList keggGet
#' @importFrom pbapply pblapply
#' @importFrom purrr map
#' @importFrom crayon yellow
#' @importFrom stringr str_replace_all
#'
#' @source This function is adapted from the metpath package by Xiaotao Shen
#'         \url{https://github.com/tidymass/metpath/blob/main/R/10_kegg_database.R}
#' @references
#' Shen, X., Yan, H., Wang, C. et al. TidyMass an object-oriented reproducible analysis framework
#' for LC–MS data. Nat Commun 13, 4365 (2022).
#'
#' @keywords internal

get_kegg_pathways <- function(local = TRUE,
                             organism = "hsa") {
  # organism = match.arg(organism)
  if (local) {
    if (organism == "hsa") {
      data("kegg_hsa_pathway", package = "metpath", envir = environment())
      message(
        crayon::yellow(
          "This database is downloaded in",
          kegg_hsa_pathway@database_info$version
        )
      )
      # cat("\n")
      return(kegg_hsa_pathway)
    }
  } else{
    message(crayon::yellow("It may take a while...\n"))
    # organism = match.arg(organism)
    pathway_ID <-
      KEGGREST::keggList(database = "pathway", organism = organism) %>%
      names() %>%
      unique() %>%
      stringr::str_replace_all(., "path:", "")

    kegg_pathway_database <-
      pbapply::pblapply(pathway_ID, function(x) {
        KEGGREST::keggGet(dbentries = x)[[1]]
      })

    pathway_id =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        unname(x$ENTRY)
      }) %>%
      unlist()

    pathway_name =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        unname(x$PATHWAY_MAP)
        # stringr::str_split(pattern = " - ") %>%
        # `[[`(1) %>%
        # `[`(1)
      }) %>%
      unlist()

    pathway_name =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        unname(x$PATHWAY_MAP)
      }) %>%
      unlist()

    describtion =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        unname(x$DESCRIPTION)
      })

    pathway_class =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        unname(x$CLASS)
      })

    gene_list =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        gene = x$GENE
        if (is.null(gene)) {
          return(data.frame())
        }
        data.frame(
          KEGG.ID = gene[seq(1, length(gene) - 1, by = 2)],
          Gene.name = gene[seq(2, length(gene), by = 2)],
          stringsAsFactors = FALSE
        )
      })

    compound_list =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          KEGG.ID = names(x$COMPOUND),
          Compound.name = x$COMPOUND,
          stringsAsFactors = FALSE
        )
      })

    reference_list =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        purrr::map(
          x$REFERENCE,
          .f = function(y) {
            y = lapply(y, function(z) {
              if (length(z) > 1) {
                paste(z, collapse = "{}")
              } else{
                z
              }
            })
            y = unlist(y)
            if (any(names(y) == "JOURNAL")) {
              names(y)[names(y) == "JOURNAL"] = "JOURNAL1"
              c(y, JOURNAL2 = "")
            }
          }
        ) %>%
          do.call(rbind, .) %>%
          as.data.frame()
      })

    related_disease =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          Disease.ID = names(x$DISEASE),
          Disease.name = x$DISEASE,
          stringsAsFactors = FALSE
        )
      })


    related_module =
      kegg_pathway_database %>%
      purrr::map(function(x) {
        data.frame(
          Module.ID = names(x$MODULE),
          Module.name = x$MODULE,
          stringsAsFactors = FALSE
        )
      })

    pathway =
      new(
        Class = "pathway_database",
        database_info = list(source = "KEGG", version = as.character(Sys.Date())),
        pathway_id = pathway_id,
        pathway_name = pathway_name,
        describtion = describtion,
        pathway_class = pathway_class,
        gene_list = gene_list,
        compound_list = compound_list,
        protein_list = list(),
        reference_list = reference_list,
        related_disease = related_disease,
        related_module = related_module
      )

    if (length(pathway@gene_list) == 0) {
      pathway@gene_list = vector(mode = "list",
                                 length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }

    if (length(pathway@compound_list) == 0) {
      pathway@compound_list = vector(mode = "list",
                                     length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }

    if (length(pathway@protein_list) == 0) {
      pathway@protein_list = vector(mode = "list",
                                    length = length(pathway@pathway_id)) %>%
        purrr::map(function(x) {
          x = data.frame()
          x
        })
    }

    return(pathway)
  }
}

#' Get HMDB pathway database
#'
#' Retrieve pathway information from the Human Metabolome Database (HMDB/SMPDB).
#' This is an internal function used by get_met_pathway_db().
#'
#' @param threads Integer. Number of threads to use for parallel processing (default: 3)
#' @return A pathway_database class object containing HMDB pathway information
#'
#' @importFrom crayon yellow
#'
#' @source This function is adapted from the metpath package by Xiaotao Shen
#'         \url{https://github.com/tidymass/metpath/blob/main/R/9_hmdb_database.R}
#'
#' @references
#' Shen, X., Yan, H., Wang, C. et al. TidyMass an object-oriented reproducible analysis framework
#' for LC–MS data. Nat Commun 13, 4365 (2022).
#'
#' @keywords internal

get_hmdb_pathways <-
  function(threads = 3) {
    data("hmdb_pathway", package = "metpath", envir = environment())
    message(
      crayon::yellow(
        "This database is downloaded in",
        hmdb_pathway@database_info$version
      )
    )
    # cat("\n")
    return(hmdb_pathway)
  }


#' Unify Various ID Types to One Unified ID Type
#'
#' @description
#' An internal function that converts different types of identifiers into a unified format.
#' For genes, it converts to Ensembl IDs, and for metabolites, it converts to HMDB IDs.
#'
#' @param ids A character vector of IDs to be unified.
#' @param variable_info A data frame containing mapping information between different ID types.
#'   For genes, should contain columns: "ensembl", "uniprot", and "entrezid".
#'   For metabolites, should contain columns: "hmdbid" and "keggid".
#' @param query_type Character string specifying the type of IDs to unify.
#'   Must be either "gene" or "metabolite". Default is "gene".
#'
#' @return A character vector of unified IDs in the standard format (Ensembl IDs for genes,
#'   HMDB IDs for metabolites).
#'
#'
#' @examples
#' \dontrun{
#' # Gene ID unification
#' gene_info <- data.frame(
#'   ensembl = c("ENSG00000139618", "ENSG00000141510"),
#'   uniprot = c("P51587", "P04637"),
#'   entrezid = c("675", "7157")
#' )
#' unify_id_internal(gene_info, "gene", c("P51587", "7157"))
#'
#' # Metabolite ID unification
#' metabolite_info <- data.frame(
#'   hmdbid = c("HMDB0000001", "HMDB0000002"),
#'   keggid = c("C00001", "C00002")
#' )
#' unify_id_internal(metabolite_info, "metabolite", c("HMDB0000001", "C00002"))
#' }
#'
#' @keywords internal

unify_id_internal <- function(ids = NULL,
                              variable_info = NULL,
                              query_type = c("gene", "metabolite")) {

  query_type <- match.arg(query_type)

  if (query_type == "gene") {
    unified_ids <-
      ids %>%
      purrr::map_chr(function(x) {
        if (stringr::str_detect(x, "ENS")) {
          return(x)
        }

        if (x %in% variable_info$symbol) {
          return(variable_info$ensembl[match(x, variable_info$symbol)])
        }

        if (x %in% variable_info$uniprot) {
          return(variable_info$ensembl[match(x, variable_info$uniprot)])
        }

        if (x %in% variable_info$entrezid) {
          return(variable_info$ensembl[match(x, variable_info$entrezid)])
        }
      })
  } else if (query_type == "metabolite") {
    if ("hmdbid" %in% colnames(variable_info)) {
      unified_ids <-
        ids %>%
        purrr::map_chr(function(x) {
          if (stringr::str_detect(x, "HMDB")) {
            return(x)
          }

          if (stringr::str_detect(x, "C")) {
            return(variable_info$hmdbid[match(x, variable_info$keggid)])
          }
        })
    } else {
      unified_ids <- ids
    }
  }

  unified_ids <- unified_ids[!is.na(unified_ids)]

  return(unified_ids)
}

# Function to extract LLM module interpretation data into a data frame
extract_llm_module_data <- function(llm_module_interpretation) {

  # Get the names of all functional modules
  module_names <- names(llm_module_interpretation)

  # Initialize empty vectors for each column
  module <- character()
  llm_interpreted_name <- character()
  llm_summary <- character()
  llm_phenotype_analysis <- character()
  confidence_score <- numeric()

  # Loop through each functional module
  for (i in seq_along(module_names)) {
    module_name <- module_names[i]
    module_data <- llm_module_interpretation[[module_name]]

    # Extract data from generated_name section
    generated_name <- module_data$generated_name

    # Append data to vectors
    module <- c(module, module_name)
    llm_interpreted_name <- c(llm_interpreted_name,
                              ifelse(is.null(generated_name$module_name),
                                     NA_character_,
                                     generated_name$module_name))
    llm_summary <- c(llm_summary,
                     ifelse(is.null(generated_name$summary),
                            NA_character_,
                            generated_name$summary))

    llm_phenotype_analysis <- c(llm_phenotype_analysis,
                                ifelse(is.null(generated_name$phenotype_analysis),
                                       NA_character_,
                                       generated_name$phenotype_analysis))

    confidence_score <- c(confidence_score,
                          ifelse(is.null(generated_name$confidence_score),
                                 NA_real_,
                                 generated_name$confidence_score))
  }

  # Create the data frame
  result_df <- data.frame(
    module = module,
    llm_interpreted_name = llm_interpreted_name,
    llm_summary = llm_summary,
    llm_phenotype_analysis = llm_phenotype_analysis,
    confidence_score = confidence_score,
    stringsAsFactors = FALSE
  )

  return(result_df)
}

## Internal functions to get clustering results
merge_by_binary_cut <- function(sim_matrix,
                                sim.cutoff) {
  requireNamespace("flexclust", quietly = TRUE)
  clusters <- simplifyEnrichment::binary_cut(mat = sim_matrix, cutoff = 1 - sim.cutoff)
  cluster_result <-
    data.frame(node = rownames(sim_matrix),
               module = paste("Functional_module", as.character(clusters), sep = "_"))

  return(cluster_result)
}

merge_by_hierarchical <- function(sim_matrix,
                                  hclust.method,
                                  sim.cutoff) {
  cosine_dist <- 1 - sim_matrix
  ## Convert distance matrix to a 'dist' object
  cosine_dist_obj <- as.dist(cosine_dist)
  ## Perform hierarchical clustering
  hc <- hclust(cosine_dist_obj, method = hclust.method)

  clusters <- cutree(hc, h = 1 - sim.cutoff)
  cluster_result <-
    data.frame(node = hc$labels,
               module = paste("Functional_module", as.character(clusters), sep = "_"))

  return(cluster_result)
}

# merge_by_Girvan_Newman <- function(edge_data,
#                                    node_data,
#                                    sim.cutoff) {
#   ## Filter graph data according to sim.cutoff
#   edge_data <-
#     edge_data |>
#     dplyr::filter(sim > sim.cutoff)
#   ## Create tidygraph object
#   graph_data <-
#     tidygraph::tbl_graph(nodes = node_data,
#                          edges = edge_data,
#                          directed = FALSE,
#                          node_key = "node") |>
#     dplyr::mutate(degree = tidygraph::centrality_degree())
#
#   ## Perform clustering
#   subnetwork <-
#     suppressWarnings(igraph::cluster_edge_betweenness(graph = graph_data, weights = abs(igraph::edge_attr(graph_data, "sim"))))
#   ## Assign functional module label for pathways
#   cluster_result <-
    # data.frame(node = node_data$node,
    #            module = paste("Functional_module", as.character(igraph::membership(subnetwork)), sep = "_"))
#
#   return(cluster_result)
# }

calculate_modularity <- function(sim_matrix, edge_data, sim.cutoff, clusters) {
  tryCatch({
    filtered_edges <- edge_data[edge_data$sim >= sim.cutoff, ]
    # Handle case where all nodes are in the same cluster
    if (length(unique(clusters)) <= 1) {
      return(0)
    }
    graph_obj <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE,
                                               vertices = rownames(sim_matrix))
    igraph::modularity(graph_obj, clusters)
  }, error = function(e) {
    warning(paste("Modularity calculation failed:", e$message))
    return(NA)
  })
}

# calculate_silhouette <- function(sim_matrix, clusters) {
#   dist_matrix <- 1 - sim_matrix
#   dist_obj <- as.dist(dist_matrix)
#   tryCatch({
#     if (length(unique(clusters)) < 2 || length(unique(clusters)) > (nrow(sim_matrix) - 1)) return(NA_real_)
#
#     sil <- cluster::silhouette(clusters, dist_obj)
#     mean(sil[, "sil_width"])
#   }, error = function(e) {
#     warning(paste("Silhouette calculation failed:", e$message))
#     return(NA)
#   })
# }

calculate_silhouette <- function(sim_matrix, clusters) {

  cluster_counts <- table(clusters)
  non_singleton_clusters <- names(cluster_counts[cluster_counts > 1])

  # If there are fewer than two non-singleton clusters, silhouette score is not meaningful
  if (length(non_singleton_clusters) < 2) {
    warning("Cannot calculate silhouette score with fewer than two non-singleton clusters.")
    return(NA_real_)
  }

  indices_to_keep <- which(clusters %in% non_singleton_clusters)

  clusters_filtered <- clusters[indices_to_keep]
  sim_matrix_filtered <- sim_matrix[indices_to_keep, indices_to_keep, drop = FALSE]

  dist_matrix <- 1 - sim_matrix_filtered
  dist_obj <- as.dist(dist_matrix)

  tryCatch({
    sil <- cluster::silhouette(clusters_filtered, dist_obj)
    mean(sil[, "sil_width"])
  }, error = function(e) {
    warning(paste("Silhouette calculation failed after removing singletons:", e$message))
    return(NA_real_)
  })
}


#' Parse tidymass_parameter Object to Data Frame
#'
#' This function takes a `tidymass_parameter` object and converts it into a data frame for easier manipulation and reporting.
#' This is a modified version of the original function that excludes the "by" parameter from the do_gsea() output.
#'
#' @param object A `tidymass_parameter` object to be parsed.
#'
#' @return A data frame containing the package name, function name, parameters (excluding "by" for gsea output), and the time when the function was called.
#'
#' @note This function excludes parameters named "by" from the final output, unlike the original implementation.
#'
#' @references
#' Original function adapted from the massdataset package:
#' Shen, X., Yan, H., Wang, C. et al. TidyMass an object-oriented reproducible analysis framework for LC–MS data. Nat Commun 13, 4365 (2022).
#' Author: Xiaotao Shen \email{shenxt1990@@outlook.com}
#'
#' @export

mapa_parse_tidymass_parameter <-
  function(object) {
    if (!is(object, class2 = "tidymass_parameter")) {
      stop("only support tidymass_parameter class.\n")
    }

    if (is.null(names(object@parameter))) {
      names(object@parameter) = paste("parameter",
                                      seq_along(object@parameter),
                                      sep = "_")
    }

    result <-
      data.frame(
        pacakge_name = object@pacakge_name,
        function_name = object@function_name,
        parameter = purrr::map2(names(object@parameter),
                                object@parameter,
                                function(name, value) {
                                  if (name == "by") {
                                    return(NULL)
                                  }
                                  if (length(value) > 5) {
                                    value = head(value, 5)
                                    value = paste(c(value, "..."), collapse = ',')
                                  } else{
                                    value = paste(value, collapse = ',')
                                  }
                                  paste(name, value, sep = ":")
                                }) %>%
          purrr::compact() %>%
          unlist(),
        time = object@time,
        check.names = FALSE
      )

    return(result)

  }
