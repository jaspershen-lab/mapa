# ####example
# setwd(r4projects::get_project_wd())
# setwd("demo_data/")
#
# load("demo_data/demo_mouse/liver/gsea_res/openai_module_annotation_res.rda")
#
# report_functional_module(
#   object = llm_annotated_modules,
#   type = "html"
# )

#' Generate a Report for a Functional Module
#'
#' This function generates a report for a given object of class 'functional_module'.
#' It creates the report in different formats (HTML, PDF, Word) based on user choice.
#' The function also performs various tasks such as creating directories,
#' fetching report templates, saving parameter and object data,
#' and rendering bar plots for functional modules, modules, and pathways, similarity networks,
#' and interpretation of functional modules.
#'
#' @param object An object of class 'functional_module'. The function will
#'               generate a report based on this object.
#' @param path A character string specifying the directory path where the report
#'             and associated files will be saved. Defaults to the current directory.
#' @param type A character vector specifying the format of the report.
#'             Can be 'html', 'pdf', 'md', 'word', or 'all' to generate the report
#'             in all three formats. Default is 'html'.
#' @param degree_cutoff A numeric value specifying the minimum degree for network
#'                     analysis. Default is 1. Must be 0 or a positive number.
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
#'
#' @export

report_functional_module <- function(object, path = ".", ...) {
  UseMethod("report_functional_module")
}

#' @rdname report_functional_module
#' @param degree_cutoff Minimum node degree for the similarity network plot
#'   (default 1).
#' @export
report_functional_module.functional_module <-
  function(object,
           path = ".",
           degree_cutoff = 1,
           type = c("html", "pdf", "word", "md", "all"),
           ...) {

    if (missing(object)) stop("object is missing")
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    options(warn = -1)

    # ── 1. Output directory ───────────────────────────────────────────────────
    out_dir    <- file.path(path, "MAPA_report")
    tables_dir <- file.path(out_dir, "tables")
    dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
    message("Report directory: ", out_dir)

    proc    <- object@process_info
    has_llm <- "llm_interpret_module" %in% names(proc)

    # ── 2. Detect query type and analysis method ──────────────────────────────
    qt <- tryCatch(
      if ("enrich_pathway" %in% names(proc))
        proc$enrich_pathway@parameter$query_type
      else
        proc$do_gsea@parameter$query_type,
      error = function(e) "gene"
    )
    analysis_type <- if ("enrich_pathway" %in% names(proc)) "ORA" else "GSEA"

    fm  <- object@merged_module$functional_module_result
    rwm <- object@merged_module$result_with_module

    # ── 3. Copy MAPA logo ─────────────────────────────────────────────────────
    logo_src <- system.file(
      "rmarkdown/templates/mapa/skeleton/mapa_logo.png",
      package = "mapa"
    )
    has_logo <- nzchar(logo_src) && file.exists(logo_src)
    if (has_logo)
      file.copy(logo_src, file.path(out_dir, "mapa_logo.png"), overwrite = TRUE)

    # ── 4. Export CSV tables ──────────────────────────────────────────────────
    message("Writing CSV tables ...")

    utils::write.csv(fm,  file.path(tables_dir, "functional_module_table.csv"),
                     row.names = FALSE)
    utils::write.csv(rwm, file.path(tables_dir, "pathway_table.csv"),
                     row.names = FALSE)

    if (has_llm) {
      llm_tbl <- do.call(rbind, lapply(
        names(object@llm_module_interpretation), function(mod) {
          x <- tryCatch(
            object@llm_module_interpretation[[mod]]$generated_name,
            error = function(e) list()
          )
          data.frame(
            module             = mod,
            llm_module_name    = .so_str(tryCatch(x$module_name,        error = function(e) NA)),
            summary            = .so_str(tryCatch(x$summary,             error = function(e) NA)),
            phenotype_analysis = .so_str(tryCatch(x$phenotype_analysis,  error = function(e) NA)),
            confidence_score   = .so_str(tryCatch(x$confidence_score,    error = function(e) NA)),
            stringsAsFactors   = FALSE
          )
        }
      ))
      utils::write.csv(llm_tbl,
                       file.path(tables_dir, "llm_interpretation_table.csv"),
                       row.names = FALSE)
    }

    # ── 5. Build summary and report data ─────────────────────────────────────
    message("Preparing report data ...")

    ep_proc   <- if ("enrich_pathway" %in% names(proc)) proc$enrich_pathway
                 else proc$do_gsea
    databases <- tryCatch(
      paste(ep_proc@parameter$database, collapse = ", "),
      error = function(e) "—"
    )

    summary_df <- data.frame(
      Item = c(
        "Analysis mode", "Query type", "Analysis method", "Database(s)",
        "Enriched pathways", "Functional modules identified",
        "LLM interpretation", "MAPA version"
      ),
      Value = c(
        "Single-Omics",
        if (qt == "gene") "Gene" else "Metabolite",
        analysis_type,
        databases,
        as.character(nrow(rwm)),
        as.character(nrow(fm)),
        if (has_llm) "Yes" else "No",
        as.character(utils::packageVersion("mapa"))
      ),
      stringsAsFactors = FALSE
    )

    n_show      <- min(30L, nrow(fm))
    plot_height <- max(4, n_show * 0.28)

    report_data <- list(
      summary_df    = summary_df,
      param_tables  = .build_so_param_tables(proc, qt, analysis_type),
      has_llm       = has_llm,
      n_show        = n_show,
      degree_cutoff = degree_cutoff,
      has_logo      = has_logo
    )

    rda_path     <- file.path(out_dir, "report_data.rda")
    rda_obj_path <- file.path(out_dir, "report_object.rda")
    save(report_data, file = rda_path)
    save(object,      file = rda_obj_path)

    # ── 6. Write and render Rmd ───────────────────────────────────────────────
    message("Rendering HTML report ...")
    rmd_path <- file.path(out_dir, "mapa_so_report.Rmd")
    writeLines(.build_so_rmd_lines(n_show, plot_height), rmd_path)

    rmarkdown::render(
      input         = rmd_path,
      output_format = rmarkdown::html_document(
        theme          = "flatly",
        toc            = TRUE,
        toc_float      = TRUE,
        self_contained = TRUE,
        df_print       = "paged"
      ),
      output_file = "index.html",
      output_dir  = out_dir,
      quiet       = TRUE
    )

    # ── 7. Clean up temp files ────────────────────────────────────────────────
    unlink(rmd_path)
    unlink(rda_path)
    unlink(rda_obj_path)
    if (has_logo) unlink(file.path(out_dir, "mapa_logo.png"))

    message("Done. Report saved to: ", file.path(out_dir, "index.html"))
    invisible(file.path(out_dir, "index.html"))
}

# ── Internal helpers (single-omics) ──────────────────────────────────────────

.so_str <- function(x) {
  if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
}

.so_safe <- function(params, key) {
  v <- tryCatch(params[[key]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0) "—" else paste(v, collapse = ", ")
}

.build_so_param_tables <- function(proc, qt, analysis_type) {
  out <- list()

  ep <- if ("enrich_pathway" %in% names(proc)) proc$enrich_pathway
        else proc$do_gsea
  if (!is.null(ep)) {
    params <- tryCatch(ep@parameter, error = function(e) list())
    out$enrichment <- data.frame(
      Parameter = c("Method", "Query type", "Database(s)", "P-value cutoff",
                    "P-adjust method", "Min gene set size", "Max gene set size"),
      Value = c(
        analysis_type,
        if (qt == "gene") "Gene" else "Metabolite",
        .so_safe(params, "database"),
        .so_safe(params, "pvalueCutoff"),
        .so_safe(params, "pAdjustMethod"),
        .so_safe(params, "minGSSize"),
        .so_safe(params, "maxGSSize")
      ),
      stringsAsFactors = FALSE
    )
  }

  ll <- proc[["llm_interpret_module"]]
  if (!is.null(ll)) {
    lp <- tryCatch(ll@parameter, error = function(e) list())
    out$llm <- data.frame(
      Parameter = c("LLM model", "Embedding model", "API provider",
                    "Module size cutoff", "Phenotype context"),
      Value = c(
        .so_safe(lp, "llm_model"),
        .so_safe(lp, "embedding_model"),
        .so_safe(lp, "api_provider"),
        .so_safe(lp, "module_content_number_cutoff"),
        .so_safe(lp, "phenotype")
      ),
      stringsAsFactors = FALSE
    )
  }

  out
}

.build_so_rmd_lines <- function(n_show, plot_height) {
  barplot_chunk <- paste0(
    "```{r barplot, fig.width=9, fig.height=", round(plot_height, 2),
    ", error=TRUE}"
  )

  c(
    "---",
    'title: "MAPA Single-Omics Analysis Report"',
    "date: \"`r format(Sys.Date(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    theme: flatly",
    "    toc: true",
    "    toc_float: true",
    "    self_contained: true",
    "    df_print: paged",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
    "load(\"report_data.rda\")",
    "load(\"report_object.rda\")",
    "```",
    "",
    "```{r header-info, results='asis'}",
    "cat('<style>h1.title, .date { display:none; }</style>')",
    "mapa_ver <- report_data$summary_df$Value[",
    "  report_data$summary_df$Item == 'MAPA version']",
    "logo_html <- if (report_data$has_logo) {",
    "  '<img src=\"mapa_logo.png\" alt=\"MAPA\"",
    "    style=\"width:64px; height:64px; object-fit:contain;",
    "    flex-shrink:0; margin-right:14px;\">'",
    "} else ''",
    "cat(paste0(",
    "  '<div style=\"display:flex; align-items:center; margin:0.5rem 0 1rem 0;\">',",
    "  logo_html,",
    "  '<span style=\"font-size:1.9em; font-weight:700; color:#1a3a5c;",
    "    line-height:1.2;\">MAPA Single-Omics Analysis Report</span>',",
    "  '</div>',",
    "  '<div style=\"background:#f0f5fb; border-left:4px solid #2c6fac;",
    "    padding:12px 18px; border-radius:6px; margin-bottom:1.5rem;\">',",
    "  '<p style=\"margin:0; font-size:0.88em; color:#444; line-height:1.8;\">',",
    "  'Generated on ', format(Sys.Date(), '%B %d, %Y'), ' by the ',",
    "  '<strong>MAPA</strong> package (v', mapa_ver, '), developed by the ',",
    "  '<strong><a href=\"https://www.shen-lab.org\" target=\"_blank\"",
    "    style=\"color:#2c6fac; text-decoration:none;\">Shen Lab</a></strong>.',",
    "  ' For documentation, tutorials, and updates, visit the ',",
    "  '<a href=\"https://www.shen-lab.org/mapa-website/\" target=\"_blank\"",
    "    style=\"color:#2c6fac; font-weight:600; text-decoration:none;\">",
    "    MAPA website &#8594;</a>.',",
    "  '</p>',",
    "  '</div>'",
    "))",
    "```",
    "",
    "<hr>",
    "",
    "## 1. Report Structure",
    "",
    "This report was generated by the **MAPA** package.",
    "The output folder contains the following files:",
    "",
    "```",
    "MAPA_report/",
    "├── index.html                        # this report",
    "└── tables/",
    "    ├── functional_module_table.csv",
    "    ├── pathway_table.csv",
    "    └── llm_interpretation_table.csv   # if LLM was run",
    "```",
    "",
    "<hr>",
    "",
    "## 2. Analysis Summary",
    "",
    "```{r summary-table}",
    "knitr::kable(report_data$summary_df,",
    "             col.names = c(\"Item\", \"Value\"), align = c(\"l\", \"l\"))",
    "```",
    "",
    "<hr>",
    "",
    "## 3. Key Parameters",
    "",
    "### 3.1 Pathway Enrichment",
    "",
    "```{r enrich-params}",
    "knitr::kable(report_data$param_tables$enrichment,",
    "             col.names = c(\"Parameter\", \"Value\"),",
    "             align = c(\"l\", \"l\"))",
    "```",
    "",
    "```{r llm-header, results='asis', eval=report_data$has_llm}",
    "cat('\\n### 3.2 LLM Annotation\\n')",
    "```",
    "",
    "```{r llm-params, eval=report_data$has_llm}",
    "knitr::kable(report_data$param_tables$llm,",
    "             col.names = c(\"Parameter\", \"Value\"),",
    "             align = c(\"l\", \"l\"))",
    "```",
    "",
    "<hr>",
    "",
    "## 4. Pathway Barplot",
    "",
    paste0("Top ", n_show, " pathways ranked by adjusted p-value."),
    "",
    barplot_chunk,
    paste0("mapa::plot_pathway_bar(object, level = \"pathway\","),
    paste0("  top_n = ", n_show, "L, p.adjust.cutoff = 0.05,"),
    "  count.cutoff = 5, y_label_width = 30,",
    "  llm_text = report_data$has_llm)",
    "```",
    "",
    "<hr>",
    "",
    "## 5. Functional Module Similarity Network",
    "",
    "```{r simnet, fig.width=7, fig.height=7, error=TRUE}",
    "mapa::plot_similarity_network(object,",
    "  level = \"functional_module\",",
    "  degree_cutoff = report_data$degree_cutoff,",
    "  llm_text = report_data$has_llm)",
    "```",
    "",
    "<hr>",
    "",
    "## 6. Module Composition Table",
    "",
    "```{r module-table}",
    "fm <- object@merged_module$functional_module_result",
    "show_cols <- intersect(",
    "  c('module', 'Description', 'Count', 'pvalue', 'p_adjust', 'llm_module_name'),",
    "  colnames(fm))",
    "knitr::kable(head(fm[, show_cols, drop = FALSE], 50), align = 'l')",
    "```",
    "",
    "```{r llm-table-header, results='asis', eval=report_data$has_llm}",
    "cat('\\n<hr>\\n\\n## 7. LLM Interpretation\\n')",
    "```",
    "",
    "```{r llm-table, eval=report_data$has_llm}",
    "llm_tbl <- do.call(rbind, lapply(",
    "  names(object@llm_module_interpretation), function(mod) {",
    "    x <- tryCatch(object@llm_module_interpretation[[mod]]$generated_name,",
    "                  error = function(e) list())",
    "    data.frame(",
    "      module           = mod,",
    "      llm_module_name  = tryCatch(x$module_name,     error = function(e) NA),",
    "      summary          = tryCatch(x$summary,          error = function(e) NA),",
    "      confidence_score = tryCatch(x$confidence_score, error = function(e) NA),",
    "      stringsAsFactors = FALSE",
    "    )",
    "  }",
    "))",
    "knitr::kable(llm_tbl, align = 'l')",
    "```",
    "",
    "<hr>",
    "",
    "<p style=\"color:#888;font-size:0.85em;\">",
    paste0("Report generated by MAPA ",
           "`r report_data$summary_df$Value[",
           "report_data$summary_df$Item == 'MAPA version']`"),
    "</p>"
  )
}

#' @rdname report_functional_module
#' @export
report_functional_module.list <- function(object, path = ".", ...) {
  if (all(c("graph_data", "functional_module_result", "result_with_module")
          %in% names(object))) {
    report_multi_omics_functional_module(object = object, path = path)
  } else {
    stop("'object' is a list but does not appear to be a supported ",
         "multi-omics functional module result.")
  }
}
