# utils::globalVariables(c("count", "module", "omics_type"))
# load("demo_data/demo_multi-omics/multi_omics_modules.rda")
# plot_multi_omics_module_barplot(multi_omics_modules)

#' Plot Multi-Omics Module Composition Stacked Barplot
#'
#' Creates a horizontal stacked bar chart showing the feature composition of
#' each multi-omics functional module, coloured by omics data type
#' (Transcriptomics, Proteomics, Metabolomics).
#'
#' Genes appearing in both transcriptomics and proteomics datasets
#' (\code{dt_src = "P, T"}) are counted once for each contributing layer.
#'
#' @param object A list produced by \code{\link{merge_multi_omics_nodes}},
#'   containing \code{functional_module_result} and \code{result_with_module}.
#' @param top_n Integer. Maximum number of modules to display, ranked by
#'   \code{module_content_number} (largest first). Default \code{30}.
#' @param colors Named character vector of hex colours for each omics type.
#'   Names must include \code{"Transcriptomics"}, \code{"Proteomics"}, and
#'   \code{"Metabolomics"}. Defaults to
#'   \code{c(Transcriptomics = "#fae69e", Proteomics = "#f2b56f",
#'   Metabolomics = "#71b7ed")}.
#' @param bar_width Numeric in \code{(0, 1]}. Width of each bar relative to
#'   the spacing between bars. Default \code{0.7}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#' \dontrun{
#' load("multi_omics_modules.rda")
#' plot_multi_omics_module_barplot(multi_omics_modules, top_n = 20)
#'
#' # Custom colours and thinner bars
#' plot_multi_omics_module_barplot(
#'   multi_omics_modules,
#'   colors = c(
#'     Transcriptomics = "#1f77b4",
#'     Proteomics      = "#ff7f0e",
#'     Metabolomics    = "#2ca02c"
#'   ),
#'   bar_width = 0.5
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual labs
#'   theme_bw theme element_blank
#' @export
plot_multi_omics_module_barplot <- function(
    object,
    top_n     = 30,
    colors    = c(
      Transcriptomics = "#fae69e",
      Proteomics      = "#f2b56f",
      Metabolomics    = "#71b7ed"
    ),
    bar_width = 0.7) {

  if (!all(c("result_with_module", "functional_module_result") %in%
             names(object))) {
    stop("'object' must contain 'result_with_module' and",
         " 'functional_module_result'.")
  }

  rwm <- object$result_with_module
  rwm <- rwm[
    rwm$node_type %in% c("gene", "metabolite") & !is.na(rwm$module),
    , drop = FALSE
  ]

  # Build long-format data frame: one row per (module, omics_type) contribution.
  # Genes with dt_src "P, T" are expanded into two rows (one per layer).
  rows <- lapply(seq_len(nrow(rwm)), function(i) {
    node <- rwm[i, ]
    if (node$node_type == "metabolite") {
      return(data.frame(
        module     = node$module,
        omics_type = "Metabolomics",
        stringsAsFactors = FALSE
      ))
    }
    dt_src <- tryCatch(node$node_info[[1]]$dt_src, error = function(e) "T")
    if (is.null(dt_src) || is.na(dt_src) || !nzchar(dt_src)) dt_src <- "T"
    sources <- trimws(strsplit(dt_src, ",", fixed = TRUE)[[1]])
    do.call(rbind, lapply(sources, function(s) {
      data.frame(
        module     = node$module,
        omics_type = if (s == "T") "Transcriptomics" else "Proteomics",
        stringsAsFactors = FALSE
      )
    }))
  })
  long_df <- do.call(rbind, rows)

  # Count per module x omics_type
  count_df <- as.data.frame(
    table(module = long_df$module, omics_type = long_df$omics_type),
    stringsAsFactors = FALSE
  )
  count_df <- count_df[count_df$Freq > 0L, , drop = FALSE]
  colnames(count_df)[colnames(count_df) == "Freq"] <- "count"

  # Select top_n modules ordered by module_content_number
  fm <- object$functional_module_result[,
    c("module", "module_content_number"), drop = FALSE
  ]
  fm       <- fm[order(fm$module_content_number, decreasing = TRUE), ]
  top_mods <- head(as.character(fm$module), top_n)

  count_df <- count_df[count_df$module %in% top_mods, , drop = FALSE]
  count_df$module <- factor(count_df$module, levels = rev(top_mods))
  count_df$omics_type <- factor(
    count_df$omics_type,
    levels = c("Transcriptomics", "Proteomics", "Metabolomics")
  )

  ggplot2::ggplot(count_df,
                  ggplot2::aes(x = count, y = module, fill = omics_type)) +
    ggplot2::geom_col(width = bar_width, color = "black", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = colors, name = "Omics Type",
                               drop = FALSE) +
    ggplot2::labs(
      x     = "Number of Features",
      y     = "Module ID",
      title = "Module Composition by Omics Type"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position    = "right",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
}
