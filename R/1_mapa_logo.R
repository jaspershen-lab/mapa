#' Display the mapa Logo
#'
#' This function displays the logo and version information for the mapa package.
#'
#' @return None. This function is called for its side effect of printing the logo and messages.
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @importFrom dplyr filter mutate select everything left_join syms case_when pull
#' @importFrom utils packageDescription write.csv
#' @importFrom cli rule symbol cli_abort
#' @importFrom crayon green blue col_align red black white style make_style num_colors
#' @import openxlsx
#' @importFrom purrr map map2
#' @import ggplot2
#' @importFrom methods .hasSlot new is
#' @importFrom stats p.adjust rgamma sd median time
#' @importFrom utils data str head tail packageVersion write.table read.delim
#' @importFrom utils read.table
#' @importFrom tibble add_column
#' @importFrom magrittr %>%
#' @importFrom tidygraph activate tbl_graph centrality_degree as_tbl_graph
#' @importFrom stringr str_wrap str_split str_replace_all str_sort str_detect
#' @importFrom ggwordcloud geom_text_wordcloud
#' @importFrom grid arrow
#' @importFrom extrafont loadfonts
#' @importFrom igraph bipartite_mapping edge_attr membership upgrade_graph vertex_attr
#' @importFrom igraph cluster_edge_betweenness
#' @importFrom massdataset activate_mass_dataset extract_variable_info
#' @importFrom simplifyEnrichment term_similarity_from_KEGG term_similarity_from_Reactome GO_similarity
#' @importFrom httr POST content_type_json content
#' @import ggraph
#' @import igraph
#' @importFrom plyr dlply .
#' @import patchwork
#' @importFrom ggforce geom_link
#' @importClassesFrom massdataset tidymass_parameter
#' @export
#' @examples
#' mapa_logo()

mapa_logo <- function() {
  message("Thank you for using mapa!")
  message("Version ",
          mapa_version,
          " (",
          update_date,
          ')')
  message("More information: jaspershen.github.io")
  cat(
    c(
      "",
      "",
      "  _ __ ___   __ _ _ __   __ _",
      " | '_ ` _ \\ / _` | '_ \\ / _` |",
      " | | | | | | (_| | |_) | (_| |",
      " |_| |_| |_|\\__,_| .__/ \\__,_|",
      "                 | |",
      "                 |_|"
    ),
    sep = "\n"
  )
}

mapa_version <-
  as.character(utils::packageVersion(pkg = "mapa"))

update_date <- as.character(Sys.time())

#' Retrieve the Version of mapa Package
#'
#' This function retrieves the current version of the mapa package installed in the R environment.
#'
#' @return A character string indicating the version of the mapa package.
#' @author Xiaotao Shen \email{shenxt1990@@outlook.com}
#' @examples
#' get_mapa_version()
#' @export
get_mapa_version = function() {
  return(as.character(utils::packageVersion(pkg = "mapa")))
}

# library(cowsay)
# # https://onlineasciitools.com/convert-text-to-ascii-art
# # writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("logo.txt")
# dput(art)
# mapa_logo <-
#   c("", "", "  _ __ ___   __ _ _ __   __ _", " | '_ ` _ \\ / _` | '_ \\ / _` |",
#     " | | | | | | (_| | |_) | (_| |", " |_| |_| |_|\\__,_| .__/ \\__,_|",
#     "                 | |", "                 |_|")
# cat(mapa_logo, sep = "\n")
