.onAttach <- function(libname, pkgname) {
  needed <- core[!is_attached(core)]
  if (length(needed) == 0)
    return()

  crayon::num_colors(TRUE)
  mapa_attach()

  # if (!"package:conflicted" %in% search()) {
  #   x <- mapa_conflicts()
  #   msg(mapa_conflict_message(x), startup = TRUE)
  # }
  packageStartupMessage(paste0("mapa ", mapa_version, " (", update_date, ')'))
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}

utils::globalVariables(c(
  "llm_module_name", "p_adjust", "mapped_number", "mapped_id", "pathway_name",
  "hmdbid", "keggid", "variable_id", "queryID", "p_value_adjust", "hmdb_pathway",
  "kegg_hsa_pathway", "all_id", "all_number", "mapped_percentage", "p_value",
  "degree", "describtion", "RichFactor", "FoldEnrichment", "zScore", "Var1",
  "Var2", "edges", "nodes", "qscore", "cx", "cy", "x_axis", "y_axis", "text_field",
  "module_content", "module_content_number.y", "label", "Compound.name", "HMDB.ID",
  "pmids"
))
