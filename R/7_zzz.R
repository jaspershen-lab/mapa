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
