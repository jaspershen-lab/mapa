enriched_pathways <-
  enrich_pathway(
    variable_info,
    database = input$pathway_database,
    save_to_local = FALSE,
    path = "result",
    OrgDb = org.Hs.eg.db,
    organism = input$organism,
    keyType = "ENTREZID",
    use_internal_data = FALSE,
    ont = "ALL",
    pvalueCutoff = input$p_value_cutoff,
    pAdjustMethod = input$p_adjust_method,
    qvalueCutoff = 0.2,
    minGSSize = input$gene_set_size[1],
    maxGSSize = input$gene_set_size[2],
    readable = FALSE,
    pool = FALSE
  )
