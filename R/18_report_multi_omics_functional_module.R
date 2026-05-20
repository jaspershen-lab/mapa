# load("demo_data/demo_multi-omics/multi_omics_modules.rda")
# report_functional_module(multi_omics_modules)

# ── Multi-Omics Report ────────────────────────────────────────────────────────

#' Generate a Report for Multi-Omics Functional Modules
#'
#' Produces a self-contained HTML report together with supporting CSV tables
#' for a multi-omics analysis result produced by
#' \code{\link{merge_multi_omics_nodes}} (and optionally annotated by
#' \code{\link{llm_interpret_module}}).
#'
#' Output layout:
#' \preformatted{
#' MAPA_multiomics_report/
#' ├── index.html
#' └── tables/
#'     ├── module_component_table.csv
#'     ├── component_information_table.csv
#'     └── llm_interpretation_table.csv   # only when LLM was run
#' }
#'
#' @param object A list produced by \code{merge_multi_omics_nodes()}.
#' @param path Character. Directory in which to create the report folder.
#'   Defaults to the current working directory.
#'
#' @return Invisibly returns the path to the created \code{index.html}.
#'
#' @importFrom utils packageVersion write.csv
#' @export
report_multi_omics_functional_module <- function(object, path = ".") {

  if (!all(c("graph_data", "functional_module_result", "result_with_module")
           %in% names(object))) {
    stop("'object' does not appear to be a multi-omics functional module list.")
  }

  # ── 1. Output directory structure ─────────────────────────────────────────
  out_dir    <- file.path(path, "MAPA_multiomics_report")
  tables_dir <- file.path(out_dir, "tables")
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  message("Report directory: ", out_dir)

  proc    <- attr(object, "process_info")
  has_llm <- "llm_interpret_multi_omics_module" %in% names(proc)
  fm      <- object$functional_module_result
  rwm     <- object$result_with_module

  # ── 2. Copy MAPA logo ─────────────────────────────────────────────────────
  logo_src <- system.file(
    "rmarkdown/templates/mapa/skeleton/mapa_logo.png",
    package = "mapa"
  )
  has_logo <- nzchar(logo_src) && file.exists(logo_src)
  if (has_logo) {
    file.copy(logo_src, file.path(out_dir, "mapa_logo.png"), overwrite = TRUE)
  }

  # ── 3. Export CSV tables ───────────────────────────────────────────────────
  message("Writing CSV tables ...")

  ## Module component table
  mc_tbl <- fm
  if (has_llm && "llm_module_interpretation" %in% names(object)) {
    llm_names <- vapply(
      names(object$llm_module_interpretation), function(mod) {
        v <- tryCatch(
          object$llm_module_interpretation[[mod]]$generated_name$module_name,
          error = function(e) NA_character_
        )
        if (is.null(v)) NA_character_ else as.character(v)
      }, character(1L)
    )
    llm_df <- data.frame(
      module          = names(object$llm_module_interpretation),
      llm_module_name = unname(llm_names),
      stringsAsFactors = FALSE
    )
    mc_tbl <- merge(mc_tbl, llm_df, by = "module", all.x = TRUE)
  }
  utils::write.csv(mc_tbl,
                   file.path(tables_dir, "module_component_table.csv"),
                   row.names = FALSE)

  ## Component information table (flatten node_info to key fields)
  flat_rows <- lapply(seq_len(nrow(rwm)), function(i) {
    base <- data.frame(
      node_id     = rwm$node_id[i],
      node_type   = rwm$node_type[i],
      node        = rwm$node[i],
      degree      = rwm$degree[i],
      module      = rwm$module[i],
      module_size = rwm$module_size[i],
      stringsAsFactors = FALSE
    )
    ni <- rwm$node_info[[i]]
    if (rwm$node_type[i] == "gene") {
      base$symbol   <- .null_na(ni$symbol)
      base$ensembl  <- .null_na(ni$ensembl)
      base$entrezid <- .null_na(ni$entrezid)
      base$uniprot  <- .null_na(ni$uniprot)
      base$dt_src   <- .null_na(ni$dt_src)
    } else if (rwm$node_type[i] == "metabolite") {
      base$kegg_id  <- .null_na(ni$keggid)
      base$cpd_name <- .null_na(ni$cpd_name)
    } else if (rwm$node_type[i] == "pathway") {
      base$pathway_name <- .null_na(ni$pathway_name)
      base$p_adjust     <- .null_na(ni$p_adjust)
    }
    base
  })
  comp_tbl <- dplyr::bind_rows(flat_rows)
  utils::write.csv(comp_tbl,
                   file.path(tables_dir, "component_information_table.csv"),
                   row.names = FALSE)

  ## LLM interpretation table
  if (has_llm && "llm_module_interpretation" %in% names(object)) {
    llm_tbl <- do.call(rbind, lapply(
      names(object$llm_module_interpretation), function(mod) {
        x <- object$llm_module_interpretation[[mod]]$generated_name
        data.frame(
          module             = mod,
          llm_module_name    = .null_na(
            tryCatch(x$module_name, error = function(e) NA)),
          summary            = .null_na(
            tryCatch(x$summary, error = function(e) NA)),
          phenotype_analysis = .null_na(
            tryCatch(x$phenotype_analysis, error = function(e) NA)),
          confidence_score   = .null_na(
            tryCatch(x$confidence_score, error = function(e) NA)),
          stringsAsFactors   = FALSE
        )
      }
    ))
    utils::write.csv(llm_tbl,
                     file.path(tables_dir, "llm_interpretation_table.csv"),
                     row.names = FALSE)
  }

  # ── 4. Pre-process data for the Rmd ───────────────────────────────────────
  message("Preparing report data ...")

  gene_nodes <- rwm[rwm$node_type == "gene", ]
  n_transcr  <- sum(vapply(gene_nodes$node_info, function(x) {
    ds <- x$dt_src
    !is.null(ds) && grepl("T", ds, fixed = TRUE)
  }, logical(1L)))
  n_prot <- sum(vapply(gene_nodes$node_info, function(x) {
    ds <- x$dt_src
    !is.null(ds) && grepl("P", ds, fixed = TRUE)
  }, logical(1L)))

  summary_df <- data.frame(
    Item = c(
      "Analysis mode",
      "Omics layers",
      "Transcriptomics features (genes)",
      "Proteomics features (proteins)",
      "Metabolomics features (metabolites)",
      "Enriched pathways (nodes)",
      "Functional modules identified",
      "LLM interpretation",
      "MAPA version"
    ),
    Value = c(
      "Multi-Omics",
      paste(.detect_omics_layers(proc), collapse = ", "),
      as.character(n_transcr),
      as.character(n_prot),
      as.character(sum(rwm$node_type == "metabolite")),
      as.character(sum(rwm$node_type == "pathway")),
      as.character(nrow(fm)),
      if (has_llm) "Yes" else "No",
      as.character(utils::packageVersion("mapa"))
    ),
    stringsAsFactors = FALSE
  )

  n_show     <- min(30L, nrow(fm))
  plot_height <- max(4, n_show * 0.28)

  report_data <- list(
    summary_df   = summary_df,
    param_tables = .build_mo_param_tables(proc, has_llm),
    data_sources = .build_data_sources_tables(proc),
    has_llm      = has_llm,
    n_modules    = nrow(fm),
    n_show       = n_show,
    has_logo     = has_logo
  )
  rda_path <- file.path(out_dir, "report_data.rda")
  save(report_data, file = rda_path)

  ## Save the object separately for the barplot chunk
  rda_obj_path <- file.path(out_dir, "report_object.rda")
  save(object, file = rda_obj_path)

  # ── 5. Write and render Rmd ────────────────────────────────────────────────
  message("Rendering HTML report ...")
  rmd_path <- file.path(out_dir, "mapa_mo_report.Rmd")
  writeLines(.build_mo_rmd_lines(n_show, plot_height), rmd_path)

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

  # ── 6. Clean up temp files ────────────────────────────────────────────────
  unlink(rmd_path)
  unlink(rda_path)
  unlink(rda_obj_path)
  if (has_logo) unlink(file.path(out_dir, "mapa_logo.png"))

  message("Done. Report saved to: ", file.path(out_dir, "index.html"))
  invisible(file.path(out_dir, "index.html"))
}

# ── Internal helpers ──────────────────────────────────────────────────────────

.null_na <- function(x) {
  if (is.null(x) || length(x) == 0) NA_character_
  else as.character(x[[1]])
}

.safe_param <- function(params, key, default = "—") {
  v <- tryCatch(params[[key]], error = function(e) NULL)
  if (is.null(v) || length(v) == 0) return(default)
  paste(v, collapse = ", ")
}

.detect_omics_layers <- function(proc) {
  layers <- character(0)
  if (!is.null(proc$transcriptome_enrich))
    layers <- c(layers, "Transcriptomics")
  if (!is.null(proc$proteome_enrich))
    layers <- c(layers, "Proteomics")
  if (!is.null(proc$metabolome_enrich))
    layers <- c(layers, "Metabolomics")
  if (length(layers) == 0) layers <- "Unknown"
  layers
}

.extract_enrich_params <- function(info) {
  if (is.null(info)) return(NULL)
  ep     <- info$enrich_pathway
  if (is.null(ep)) return(NULL)
  params <- tryCatch(ep@parameter, error = function(e) NULL)
  if (is.null(params)) return(NULL)
  list(
    database      = .safe_param(params, "database"),
    pvalueCutoff  = .safe_param(params, "pvalueCutoff"),
    pAdjustMethod = .safe_param(params, "pAdjustMethod"),
    minGSSize     = .safe_param(params, "minGSSize"),
    maxGSSize     = .safe_param(params, "maxGSSize")
  )
}

.build_mo_param_tables <- function(proc, has_llm) {
  out <- list()

  # Enrichment: one column per omics layer
  t_ep <- .extract_enrich_params(proc$transcriptome_enrich)
  p_ep <- .extract_enrich_params(proc$proteome_enrich)
  m_ep <- .extract_enrich_params(proc$metabolome_enrich)
  keys   <- c("database", "pvalueCutoff", "pAdjustMethod",
              "minGSSize", "maxGSSize")
  labels <- c("Database(s)", "P-value cutoff", "P-adjust method",
              "Min gene set size", "Max gene set size")
  get_val <- function(ep, k) if (!is.null(ep)) ep[[k]] else "—"
  out$enrichment <- data.frame(
    Parameter       = labels,
    Transcriptomics = vapply(keys, get_val, character(1L), ep = t_ep),
    Proteomics      = vapply(keys, get_val, character(1L), ep = p_ep),
    Metabolomics    = vapply(keys, get_val, character(1L), ep = m_ep),
    stringsAsFactors = FALSE
  )

  # Network construction
  bn <- proc$build_network_tables
  bm <- proc$build_MNetwork
  out$network <- data.frame(
    Parameter = c(
      "Organism (NCBI taxon ID)",
      "STRING PPI score cutoff",
      "TF confidence level",
      "Embedding source",
      "Embedding database version"
    ),
    Value = c(
      .safe_param(bn$parameter, "taxon_id"),
      .safe_param(bn$parameter, "string_score_cutoff"),
      .safe_param(bn$parameter, "tf_confidence_levels"),
      .safe_param(bm$parameter, "embedding_source"),
      .safe_param(bm$parameter, "db_version")
    ),
    stringsAsFactors = FALSE
  )

  # Similarity computation (RWR)
  gs <- proc$get_multi_omics_sim
  out$similarity <- data.frame(
    Parameter = c(
      "Min pathway similarity",
      "Min bipartite weight",
      "RWR restart probability (r)",
      "Decay parameter (η)",
      "Path weight (λ)",
      "Node similarity weight (δ1)",
      "Bipartite weight (δ2)"
    ),
    Value = vapply(
      c("min_path_sim", "min_bipartite_weight", "r",
        "eta", "lambda", "delta1", "delta2"),
      function(k) .safe_param(gs$parameter, k),
      character(1L)
    ),
    stringsAsFactors = FALSE
  )

  # Module identification
  mn <- proc$merge_multi_omics_nodes
  out$modules <- data.frame(
    Parameter = c("Similarity cutoff", "Clustering method"),
    Value     = c(
      .safe_param(mn$parameter, "sim_cutoff"),
      .safe_param(mn$parameter, "cluster_method")
    ),
    stringsAsFactors = FALSE
  )

  # LLM annotation (optional)
  if (has_llm) {
    ll <- proc$llm_interpret_multi_omics_module
    out$llm <- data.frame(
      Parameter = c(
        "LLM model", "Embedding model", "API provider",
        "Module size cutoff", "Phenotype context",
        "PubMed search years", "Max PubMed records",
        "Similarity filter top-N", "LLM re-ranking top-N",
        "Parallel threads"
      ),
      Value = vapply(
        c("llm_model", "embedding_model", "api_provider",
          "module_content_number_cutoff", "phenotype",
          "years", "retmax", "similarity_filter_num",
          "GPT_filter_num", "thread"),
        function(k) .safe_param(ll$parameter, k),
        character(1L)
      ),
      stringsAsFactors = FALSE
    )
  }

  out
}

.build_data_sources_tables <- function(proc) {
  out <- list()
  ds  <- proc$build_network_tables$data_sources

  # TF-target regulation
  tf <- ds$tf_target
  out$tf <- data.frame(
    Field = c("Database", "R package", "Package version",
              "Confidence levels"),
    Value = c(
      .null_na(tf$database),
      .null_na(tf$r_package),
      .null_na(tf$package_version),
      paste(tf$confidence_levels, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )

  # Protein-protein interaction
  ppi <- ds$ppi
  out$ppi <- data.frame(
    Field = c("Database", "Score cutoff",
              "Organism (NCBI taxon ID)"),
    Value = c(
      .null_na(ppi$database),
      .safe_param(proc$build_network_tables$parameter,
                  "string_score_cutoff"),
      .null_na(ppi$taxon_id)
    ),
    stringsAsFactors = FALSE
  )

  # Reaction databases (enzyme-metabolite + metabolite-metabolite)
  .fmt_source <- function(src) {
    if (is.null(src)) return(c(NA_character_, NA_character_))
    db  <- .null_na(src$database)
    acc <- if (!is.null(src$access)) {
      pkg <- .null_na(src$r_package)
      ver <- .null_na(src$package_version)
      paste0(src$access, " (", pkg, " v", ver, ")")
    } else {
      "File download"
    }
    c(db, acc)
  }
  em  <- ds$enzyme_metabolite$sources
  mr  <- ds$metabolite_reaction$sources
  em_r <- .fmt_source(em$reactome)
  em_k <- .fmt_source(em$kegg)
  mr_r <- .fmt_source(mr$reactome)
  mr_k <- .fmt_source(mr$kegg)
  out$reaction <- data.frame(
    "Edge Type" = c(
      "Enzyme-metabolite", "Enzyme-metabolite",
      "Metabolite-metabolite", "Metabolite-metabolite"
    ),
    Source   = c("Reactome", "KEGG", "Reactome", "KEGG"),
    Database = c(em_r[1], em_k[1], mr_r[1], mr_k[1]),
    Access   = c(em_r[2], em_k[2], mr_r[2], mr_k[2]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Pathway databases from embedding_db_source_metadata
  meta <- proc$build_MNetwork$parameter$embedding_db_source_metadata
  if (!is.null(meta) && length(meta) >= 3) {
    db_names    <- meta[[1]]
    db_versions <- meta[[2]]
    db_access   <- meta[[3]]
    n <- min(length(db_names), length(db_versions), length(db_access))
    out$pathway <- data.frame(
      Database = db_names[seq_len(n)],
      Version  = db_versions[seq_len(n)],
      "Access Method" = db_access[seq_len(n)],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out
}

.build_mo_rmd_lines <- function(n_show, plot_height) {
  barplot_chunk <- paste0(
    "```{r barplot, fig.width=9, fig.height=",
    round(plot_height, 2),
    ", error=TRUE}"
  )

  c(
    "---",
    'title: "MAPA Multi-Omics Analysis Report"',
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
    "    line-height:1.2;\">MAPA Multi-Omics Analysis Report</span>',",
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
    "## 1. Report Structure",
    "",
    "This report was generated by the **MAPA** package.",
    "The output folder contains the following files:",
    "",
    "```",
    "MAPA_multiomics_report/",
    "├── index.html                          # this report",
    "└── tables/",
    "    ├── module_component_table.csv",
    "    ├── component_information_table.csv",
    "    └── llm_interpretation_table.csv     # if LLM was run",
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
    "## 3. Data Sources",
    "",
    "### 3.1 TF-Target Regulation Database",
    "",
    "```{r tf-db}",
    "knitr::kable(report_data$data_sources$tf,",
    "             col.names = c(\"Field\", \"Value\"), align = c(\"l\", \"l\"))",
    "```",
    "",
    "### 3.2 Protein-Protein Interaction Database",
    "",
    "```{r ppi-db}",
    "knitr::kable(report_data$data_sources$ppi,",
    "             col.names = c(\"Field\", \"Value\"), align = c(\"l\", \"l\"))",
    "```",
    "",
    "### 3.3 Reaction Database",
    "",
    "```{r reaction-db}",
    "knitr::kable(report_data$data_sources$reaction, align = c(\"l\",\"l\",\"l\",\"l\"))",
    "```",
    "",
    "### 3.4 Pathway Database",
    "",
    "```{r pathway-db}",
    "knitr::kable(report_data$data_sources$pathway, align = c(\"l\",\"l\",\"l\"))",
    "```",
    "",
    "<hr>",
    "",
    "## 4. Key Parameters",
    "",
    "### 4.1 Pathway Enrichment",
    "",
    "```{r enrich-params}",
    "knitr::kable(report_data$param_tables$enrichment,",
    "             align = c(\"l\", \"c\", \"c\", \"c\"))",
    "```",
    "",
    "### 4.2 Network Construction",
    "",
    "```{r network-params}",
    "knitr::kable(report_data$param_tables$network,",
    "             col.names = c(\"Parameter\", \"Value\"),",
    "             align = c(\"l\", \"l\"))",
    "```",
    "",
    "### 4.3 Multi-Omics Similarity (RWR)",
    "",
    "```{r sim-params}",
    "knitr::kable(report_data$param_tables$similarity,",
    "             col.names = c(\"Parameter\", \"Value\"),",
    "             align = c(\"l\", \"l\"))",
    "```",
    "",
    "### 4.4 Module Identification",
    "",
    "```{r mod-params}",
    "knitr::kable(report_data$param_tables$modules,",
    "             col.names = c(\"Parameter\", \"Value\"),",
    "             align = c(\"l\", \"l\"))",
    "```",
    "",
    "```{r llm-header, results='asis', eval=report_data$has_llm}",
    "cat('\\n### 4.5 LLM Annotation\\n')",
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
    "## 5. Module Composition",
    "",
    paste0("Top ", n_show,
           " modules ranked by feature count, coloured by omics data type."),
    "Genes present in both transcriptomics and proteomics datasets",
    "are counted once per contributing layer.",
    "",
    barplot_chunk,
    paste0("mapa::plot_multi_omics_module_barplot(object,",
           " top_n = ", n_show, "L)"),
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
