## **Task Description**
You are an AI tasked with generating a **multi-omics functional module name, a function summary, and a confidence score** based on the following information.

A *multi-omics functional module* integrates:
- **Genes / proteins** (typically gene symbols with gene names)
- **Metabolites** (e.g., KEGG compound IDs with metabolite names)
- **Enriched pathways / GO terms** (names + short descriptions)
- **Related PubMed literature** (titles/abstracts or extracted text chunks)


## **Missing-Component Rules (IMPORTANT)**
- Each omics layer can be **missing or empty** (e.g., no pathways/GO terms for a module).
- If **PathwayNames / PathwayDescription** are missing, `NULL`, `NA`, or an empty list/vector:
  - Treat the pathway layer as **not provided**.
  - **Do not hallucinate or infer** pathway terms.
  - In the summary, you may add one short clause like: *"No pathway enrichment terms were provided for this module."*
- If a layer is provided but sparse, focus on what is available and keep uncertainty explicit.

Your goal is to produce a **coherent, biology-first interpretation** of what the module represents, emphasizing **cross-omics consistency and mechanistic links**.

---

## **Your Task**
1. **Identify the core biological process** that best explains this module **across omics layers**.
2. **Emphasize cross-omics bridges** (when plausible):
   - **Gene/protein → metabolite**: enzymes, transporters, lipid remodeling, sterol synthesis/transport, redox enzymes, etc.
   - **Metabolite → gene/protein**: ligand–receptor signaling, nuclear receptor signaling, epigenetic regulation, membrane microdomain effects, etc.
   - If the bridging mechanism is not explicit in the inputs, propose a **scientifically plausible hypothesis** grounded in known biology.
3. **Avoid isolated component descriptions**:
   - Do **not** provide a long list of independent gene functions.
   - Do **not** restate pathway names without synthesis.
   - Instead, explain how components *collectively* indicate one or a few tightly related themes.
4. **Name the module concisely**:
   - The name should highlight a **key process / compartment / mechanism** (e.g., “Plasma membrane GPCR–G protein signaling and lipid microdomain remodeling”), not a generic term.
6. Assign a **confidence score (0.00–1.00)** reflecting functional coherence:
   - **High (0.8–1.0):** genes + metabolites + pathways converge on a clear theme, and literature supports it.
   - **Medium (0.4–0.79):** partial convergence or plausible but indirect links.
   - **Low (0.0–0.39):** weak convergence, largely unrelated pathways/components, or insufficient evidence.
7. **Incorporate research insights** from the provided literature (if provided):
   - Prefer mechanisms and biological claims supported by the supplied PubMed-derived texts.
   - If literature is sparse or generic, clearly separate **evidence-backed statements** from **hypotheses**.

---

## **Example Input**
The module contains:

**Genes/Proteins:**
- Gene names: ATP-binding cassette sub-family D member 1, G protein subunit alpha s, histone deacetylase 9, …

**Metabolites:**
- Names: phosphatidylcholine, cholesterol

**Enriched pathways/terms:**
plasma membrane region (A membrane that is a (regional) part of the plasma membrane.)

extrinsic component of plasma membrane (The component of a plasma membrane consisting of gene products and protein complexes that are loosely bound to one of its surfaces, but not integrated into the hydrophobic region.)

**Related articles (PubMed-derived):**
Title: Association of DNA Methylation-Derived C-Reactive Protein Predictors With All-Cause Mortality and Blood Lipids in Adults Aged ≥ 50 Years in the United States. (PubMedID:41765372)
Text: Importantly, concurrent elevations of CRPMort and HsCRP were associated with the highest risks of all-cause mortality and dyslipidemia. GrimAge2-derived CRPMort is a robust predictor of long-term all-cause mortality and may capture chronic inflammation linked to triglyceride-rich lipoproteins beyond HsCRP. Combined assessment of both inflammatory markers may enhance risk stratification and inform aging-related cardiometabolic research.

Title: ...
Text: ...


---

## **Example Output**
{
  "module_name": "Plasma membrane signaling and lipid microdomain remodeling",
  "summary": "This multi-omics module converges on plasma membrane organization and signaling. The enriched membrane-associated terms suggest a focus on proteins that localize to or regulate the cell surface. The gene/protein set can be grouped into (i) receptors and G-protein signaling components that mediate extracellular signal transduction, and (ii) transcriptional/epigenetic regulators that may tune downstream responses. The metabolite layer (phosphatidylcholine and cholesterol) supports a membrane lipid composition theme, consistent with altered membrane microdomains that can modulate receptor signaling, trafficking, and immune or stress responsiveness. Together, the module points to coordinated changes in membrane-associated signaling machinery and lipid environment, providing a mechanistic basis for altered cell–environment communication.",
  "confidence_score": "0.86"
}

---

## **Actual Input for Generation**
The multi-omics module is defined by the following components.

**Genes**
- Gene names (vector): {GeneNames_vec}

**Metabolites**
- Metabolite names (vector): {MetNames_vec}

**Pathways / Terms**
{combined_pathway_text}

**Related articles (PubMed-derived):**
{combined_texts}

---

## **Final Output**

When a component (especially pathways/terms) is not provided, keep the summary fully grounded in the available evidence (genes/metabolites/literature).
Please provide your response in **JSON format**, strictly following this structure:

{
  "module_name": "Your concise module name here",
  "summary": "A detailed, literature-integrated explanation of the module's cross-omics function (integrating **only the layers provided**: genes/proteins and/or metabolites and/or pathways/terms), emphasizing functional convergence and plausible mechanisms.",
  "confidence_score": "A value between 0.00 and 1.00, reflecting overall functional coherence across omics layers."
}

When generating your answer, pay attention to the following:
- **module_name** should be concise and specific (process/compartment/mechanism).
- **summary** should include:
  - The central biological theme(s) supported by **pathways/terms**.
  - Functional grouping of key **genes/proteins** (roles, not exhaustive lists).
  - Interpretation of **metabolites** as functional readouts and how they connect to the gene/pathway theme.
  - Key mechanisms or claims supported by the provided literature (if present), clearly distinguishing evidence vs hypothesis.
- **confidence_score** should reflect:
  - Cross-omics convergence (genes ↔ pathways ↔ metabolites)
  - Strength/clarity of mechanism
  - Support from the provided PubMed-derived texts
