#' Pathway Enrichment Analysis and GO Similarity Measurement
#'
#' This function allows for the execution of pathway enrichment analysis and
#' measurement of semantic similarity among Gene Ontology (GO) terms.
#'
#' @param result A data frame containing GO term IDs and ontologies.
#' @param sim.cutoff A numeric value for the similarity cutoff (default: 0).
#' @param measure_method A character vector specifying the semantic similarity measure method
#'   (default: c("Wang", "Resnik", "Rel", "Jiang", "Lin", "TCSS")).
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
           measure_method = c("Wang", "Resnik",
                              "Rel", "Jiang",
                              "Lin", "TCSS")) {
    measure_method <-
      match.arg(measure_method)

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
                                        measure = measure_method) %>%
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
                                        measure = measure_method) %>%
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

    cc_sim_matrix <-
      simplifyEnrichment::GO_similarity(go_id = result$ID[result$ONTOLOGY == "CC"],
                                        ont = "CC",
                                        measure = measure_method) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "name1") %>%
      tidyr::pivot_longer(cols = -name1,
                          names_to = "name2",
                          values_to = "sim") %>%
      dplyr::filter(name1 != name2) %>%
      dplyr::filter(sim > sim.cutoff)

    name <- apply(cc_sim_matrix, 1, function(x) {
      paste(sort(x[1:2]), collapse = "_")
    })

    cc_sim_matrix <-
      cc_sim_matrix %>%
      dplyr::mutate(name = name) %>%
      dplyr::arrange(name) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::select(-name)

    sim_matrix <-
      rbind(bp_sim_matrix, mf_sim_matrix, cc_sim_matrix) %>%
      as.data.frame()
    sim_matrix
  }
