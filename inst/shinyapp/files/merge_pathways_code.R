enriched_modules <-
  merge_pathways(
    object = enriched_pathways,
    p.adjust.cutoff.go = p.adjust.cutoff.go,
    p.adjust.cutoff.kegg = p.adjust.cutoff.kegg,
    p.adjust.cutoff.reactome = p.adjust.cutoff.reactome,
    count.cutoff.go = count.cutoff.go,
    count.cutoff.kegg = count.cutoff.kegg,
    count.cutoff.reactome = count.cutoff.reactome,
    sim.cutoff.go = sim.cutoff.go,
    sim.cutoff.kegg = sim.cutoff.kegg,
    sim.cutoff.reactome = sim.cutoff.reactome,
    measure.method.go = measure.method.go,
    measure.method.kegg = measure.method.kegg,
    measure.method.reactome = measure.method.reactome,
    path = "result",
    save_to_local = FALSE
  )
