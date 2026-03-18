## **Task Description**
You are an AI tasked with generating a **module name**, a **function summary**, a **{phenotype} relationship analysis**, and a **confidence score** based on the following information:

- Genes / proteins names
- Metabolites
- Enriched pathways / GO terms (names + short descriptions)
- Related PubMed literature (titles/abstracts or extracted text chunks)

Each information layer may be missing. If pathways/GO terms are missing, do not infer them; you may briefly state that no pathway enrichment terms were provided. If literature is missing or generic, clearly distinguish evidence-backed statements from hypotheses.

Your goal is to generate a coherent biological interpretation of the module, emphasizing **cross-omics consistency** and **mechanistic links**, and how the module may relate to **{phenotype}**. The summary should read like a concise biological interpretation, not a list of annotations.

---

## **Your Task**
1. Interpret the biological roles of the genes/proteins, metabolites, and pathway/GO terms.
2. Infer the main biological process or mechanism that best explains the module as a whole.
3. Explain how the different components collectively support this shared theme.
4. Avoid listing isolated component functions or simply restating pathway names.
5. **Analyze relationship to {phenotype}:**
   - Identify known associations between the module and **{phenotype}**.
   - Describe how the module's functions may influence **{phenotype}** development or progression (mechanistic reasoning).
   - If direct evidence is limited, generate a scientifically plausible hypothesis grounded in known molecular functions of the module components, analogous pathways with established links to **{phenotype}**, or potential downstream effects on tissue/organismal physiology.
6. Assign a confidence score (0.00–1.00) based on functional coherence:
   - **High (0.80–1.00):** strong convergence across components, with literature support; phenotype link is supported or strongly mechanistically grounded.
   - **Medium (0.40–0.79):** partial convergence or indirect links; phenotype link plausible but not directly supported.
   - **Low (0.00–0.39):** weak coherence or insufficient evidence; phenotype relationship highly speculative.
7. Prefer literature-supported mechanisms when available, and clearly separate supported claims from plausible interpretation.

Keep the summary concise and synthesis-focused. Prioritize the dominant shared mechanism over exhaustive component-level description.

---

## **Example Input**
The module contains:

**Genes / proteins names:**
- names: acid phosphatase 6, lysophosphatidic acid phosphatase type, inositol polyphosphate-1-phosphatase, inositol-3-phosphate synthase 1, phosphatidylinositol 4-kinase type 2 beta, phosphatidylinositol glycan anchor biosynthesis class S, phosphatidylinositol-5-phosphate 4-kinase type 2 alpha, solute carrier family 44 member 1

**Metabolites:**
- names: D-Glucose 6-phosphate

**Enriched pathways/GO terms (description):**
phospholipid metabolic process (The chemical reactions and pathways involving phospholipids, any lipid containing phosphoric acid as a mono- or diester.)

phospholipid biosynthetic process (The chemical reactions and pathways resulting in the formation of a phospholipid, a lipid containing phosphoric acid as a mono- or diester.)

phosphatidylinositol metabolic process (The chemical reactions and pathways involving phosphatidylinositol, any glycophospholipid in which a sn-glycerol 3-phosphate residue is esterified to the 1-hydroxyl group of 1D-myo-inositol.)

[Additional pathways omitted in this example for brevity.]

**Phenotype of interest:**
- phenotype: aging

**Related articles (PubMed-derived):**
Title: PIP4K2B is mechanoresponsive and controls heterochromatin-driven nuclear softening through UHRF1. (PubMedID:36918565)
Text: Phosphatidylinositol-5-phosphate (PtdIns5P)-4-kinases (PIP4Ks) are stress-regulated phosphoinositide kinases able to phosphorylate PtdIns5P to PtdIns(4,5)P2. Among the three PIP4K isoforms expressed in mammalian cells, PIP4K2B shows prominent nuclear localisation. PIP4K2B protein level strongly decreases in cells growing on soft substrates. Its silencing or pharmacological inhibition reduces the epigenetic regulator UHRF1 and induces changes in nuclear polarity, nuclear envelope tension, and chromatin compaction. This rewiring of nuclear mechanical state drives YAP cytoplasmic retention, impairs its transcriptional activity, and leads to defects in cell spreading and motility. These findings suggest that PIP4K2B links phosphoinositide metabolism to mechanoresponsive nuclear regulation.

[Additional articles omitted in this example for brevity.]

---

## **Example Output**
{
  "module_name": "Plasma membrane signaling and lipid microdomain remodeling",
  "summary": "This module is most consistently explained by a shared role in plasma membrane organization and signaling. The enriched membrane-associated terms point to proteins that localize to or regulate the cell surface. The gene/protein set can be grouped into (i) receptors and G-protein signaling components that mediate extracellular signal transduction, and (ii) transcriptional/epigenetic regulators that may tune downstream responses. The metabolite layer—phosphatidylcholine and cholesterol—supports a membrane lipid composition theme, consistent with altered membrane microdomains that can modulate receptor signaling, trafficking, and immune or stress responsiveness. Together, these layers support a coherent interpretation of the module as a membrane-associated signaling program whose activity is shaped by the surrounding lipid environment, providing a mechanistic basis for altered cell–environment communication. The extension to epigenetic regulation through histone deacetylase 9 is plausible given the downstream transcriptional consequences of G-protein and lipid signaling cascades, though direct mechanistic links to this module's specific components are supported more by functional convergence than direct literature evidence.",
  "phenotype_analysis": "This module is most consistently interpreted as an aging-related membrane phospholipid and phosphoinositide remodeling program. Multiple components converge on phosphatidylinositol metabolism: ISYNA1 links D-glucose 6-phosphate to inositol biosynthesis, providing a plausible upstream metabolic input, while PI4K2B, PIP4K2A, and INPP1 support active phosphoinositide interconversion through coordinated phosphorylation and dephosphorylation. The enriched pathways further reinforce this interpretation by repeatedly highlighting phospholipid, glycerophospholipid, and phosphatidylinositol metabolic and biosynthetic processes, suggesting that this module reflects broader membrane lipid remodeling rather than an isolated enzymatic branch. Additional genes extend this theme toward membrane homeostasis and lipid utilization, including ACP6 in phospholipid turnover, SLC44A1 in choline-related phosphatidylcholine metabolism, and PIGS in phosphatidylinositol-derived GPI-anchor biosynthesis. Together, these features support a coherent cross-omics interpretation in which altered phosphoinositide and membrane lipid metabolism may represent a coordinated remodeling state. In an aging context, this suggests a plausible mechanism whereby shifts in membrane lipid composition and phosphoinositide signaling may contribute to altered mechanosensing, stress adaptation, and membrane homeostasis during aging. While the core lipid-metabolic theme is well supported by the genes, metabolite, and pathway terms, the specific connection to aging is more interpretive and is inferred from the known roles of phosphoinositide signaling in cellular adaptation and mechanical state.",
  "confidence_score": "0.82"
}

Use the Example Input/Output only as a formatting guide. Do not reuse its biological theme unless it is directly supported by the actual input.

---

## **Actual Input for Generation**
The multi-omics module is defined by the following components.

**Genes / proteins names:**
- names: {GeneNames_vec}

**Metabolites:**
- names: {MetNames_vec}

**Enriched pathways/GO terms (description):**
{combined_pathway_text}

**Phenotype of interest:**
- phenotype: {phenotype}

**Related articles (PubMed-derived):**
{combined_texts}

---

## **Final Output**

When a component (especially pathways/terms) is not provided, keep the summary fully grounded in the available evidence (genes/metabolites/literature).
Please provide your response in **JSON format**, strictly following this structure:

{
  "module_name": "Your concise module name here",
  "summary": "A coherent biological interpretation of the module, emphasizing cross-omics consistency and mechanistic links. Integrate only the layers provided (genes/proteins and/or metabolites and/or pathways/terms). Clearly distinguish literature-supported claims from plausible hypotheses.",
  "phenotype_analysis": "An analysis of the module's relationship to {phenotype}, including evidence-backed associations (if present in the provided literature) and/or a clearly labeled mechanistic hypothesis, with 1–3 brief validation suggestions if appropriate.",
  "confidence_score": "A value between 0.00 and 1.00, reflecting overall functional coherence across omics layers and the plausibility/support of the module–{phenotype} link."
}
