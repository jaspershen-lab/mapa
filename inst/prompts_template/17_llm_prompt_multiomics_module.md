## **Task Description**
You are an AI tasked with generating a **module name**, a **function summary**, and a **confidence score** for a multi-omics functional module based on the following information:

- Genes / proteins names
- Metabolites
- Enriched pathways / GO terms (names + short descriptions)
- Related PubMed literature (titles/abstracts or extracted text chunks)

Each information layer may be missing. If pathways/GO terms are missing, do not infer them; you may briefly state that no pathway enrichment terms were provided. If literature is missing or generic, clearly distinguish evidence-backed statements from hypotheses.

Your goal is to generate a coherent biological interpretation of the module, emphasizing **cross-omics consistency** and **mechanistic links**. The summary should read like a concise biological interpretation, not a list of annotations.

---

## **Your Task**
1. Interpret the biological roles of the genes/proteins, metabolites, and pathway/GO terms.
2. Infer the main biological process or mechanism that best explains the module as a whole.
3. Explain how the different components collectively support this shared theme.
4. Avoid listing isolated component functions or simply restating pathway names.
5. Generate a concise module name that reflects a specific biological process or mechanism.
6. Assign a confidence score (0.00–1.00) based on functional coherence:
   - **High (0.8–1.0):** strong convergence across components, with literature support
   - **Medium (0.4–0.79):** partial convergence or indirect links
   - **Low (0.0–0.39):** weak coherence or insufficient evidence
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

**Related articles (PubMed-derived):**
Title: PIP4K2B is mechanoresponsive and controls heterochromatin-driven nuclear softening through UHRF1. (PubMedID:36918565)
Text: Phosphatidylinositol-5-phosphate (PtdIns5P)-4-kinases (PIP4Ks) are stress-regulated phosphoinositide kinases able to phosphorylate PtdIns5P to PtdIns(4,5)P2. Among the three PIP4K isoforms expressed in mammalian cells, PIP4K2B shows prominent nuclear localisation. PIP4K2B protein level strongly decreases in cells growing on soft substrates. Its silencing or pharmacological inhibition reduces the epigenetic regulator UHRF1 and induces changes in nuclear polarity, nuclear envelope tension, and chromatin compaction. This rewiring of nuclear mechanical state drives YAP cytoplasmic retention, impairs its transcriptional activity, and leads to defects in cell spreading and motility. These findings suggest that PIP4K2B links phosphoinositide metabolism to mechanoresponsive nuclear regulation.

[Additional articles omitted in this example for brevity.]

---

## **Example Output**
{
  "module_name": "Phosphatidylinositol and membrane phospholipid metabolism",
  "summary": "This module is most consistently explained by a shared role in phosphatidylinositol-centered membrane lipid metabolism and broader glycerophospholipid biosynthesis/remodeling. Several components converge on the inositol–phosphoinositide axis: ISYNA1 links central carbon metabolism to inositol production, making D-glucose 6-phosphate a plausible upstream precursor input for phosphatidylinositol synthesis, while PI4K2B, PIP4K2A, and INPP1 together indicate active interconversion of phosphoinositide species through phosphorylation and dephosphorylation. This core lipid-signaling branch is reinforced by the enriched pathways, which repeatedly point to phospholipid, glycerophospholipid, and phosphatidylinositol metabolic and biosynthetic processes. Additional components broaden the module toward membrane lipid remodeling and utilization: ACP6 is consistent with lysophosphatidic acid/phospholipid turnover, SLC44A1 supports the choline-related branch of phosphatidylcholine metabolism, and PIGS links phosphatidylinositol-derived lipids to GPI-anchor biosynthesis. Together, these layers support a coherent interpretation of the module as a membrane phospholipid metabolic program with a particularly strong phosphoinositide component. The provided literature most directly supports the phosphoinositide branch, showing that PIP4K2B functions as a stress- and mechanoresponsive phosphoinositide kinase that can couple lipid signaling to nuclear mechanical regulation and YAP-related transcriptional control. This literature-backed mechanism strengthens the interpretation that the module is not only involved in structural membrane lipid metabolism but may also contribute to signaling processes responsive to cellular mechanical state. The broader extension of the module to membrane remodeling and lipid biosynthesis is plausible from the combined gene, metabolite, and pathway evidence, although that part is supported more by functional convergence than by the single literature example alone.",
  "confidence_score": 0.89
}

Use the Example Input/Output only as a formatting and reasoning guide. Do not reuse its biological theme unless it is directly supported by the actual input.

---

## **Actual Input for Generation**
The multi-omics module is defined by the following components.

**Genes / proteins names:**
- names: {GeneNames_vec}

**Metabolites:**
- names: {MetNames_vec}

**Enriched pathways/GO terms (description):**
{combined_pathway_text}

**Related articles (PubMed-derived):**
{combined_texts}

---

## **Final Output**

When a component (especially pathways/terms) is not provided, keep the summary fully grounded in the available evidence (genes/metabolites/literature).
Please provide your response in **JSON format**, strictly following this structure:

{
  "module_name": "Your concise module name here",
  "summary": "A coherent biological interpretation of the module, emphasizing cross-omics consistency and mechanistic links",
  "confidence_score": "A value between 0.00 and 1.00, reflecting overall functional coherence across omics layers."
}
