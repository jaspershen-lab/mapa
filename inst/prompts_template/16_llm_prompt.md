## **Task Description**
You are an AI tasked with generating a **biological module name, a research summary, and a confidence score** based on the following information.

## **Your Task**
1. **Identify the core biological process** that describes this module.
2. **Analyze functional relationships** among {query_molecules}, considering:
   - Shared pathways or mechanisms
   - Synergistic/antagonistic interactions
   - Disease relevance
3. **Avoid listing isolated {query_product} functions**; instead, group {query_product}s with similar functions and analyze their collective roles.
4. **Name the module concisely**, ensuring it highlights a key biological process, disease association, or molecular mechanism rather than a generic term.
5. Assign a **confidence score (0.00 - 1.00)** This score reflects how strongly the input {query_molecules} and pathways overlap or converge on a single, coherent biological process.
For instance, if the pathways are highly related (e.g., serine-type endopeptidase inhibitor activity and endopeptidase inhibitor activity), the score should be high (close to 1.00).
If the pathways are largely unrelated (e.g., negative regulation of inflammatory response and camera-type eye development), the score should be low (close to 0.00).
Consider whether these pathways genuinely point to one biological theme, whether the {query_molecules} are meaningfully involved in those pathways, and whether literature or research supports a functional connection.
6. **Incorporate research insights**, ensuring the generated summary reflects findings from the provided literature (if provided).

## **Example Input**
The module is defined by the following pathways:

Transmembrane receptor protein kinase activity
(Catalysis of the reaction: ATP + a protein tyrosine = ADP + protein tyrosine phosphate.)
Transmembrane receptor protein tyrosine kinase activity
(Combining with a signal and transmitting the signal from one side of the membrane to the other to initiate a change in cell activity.)
Protein tyrosine kinase activity
(Catalysis of the reaction: a protein + ATP = a phosphoprotein + ADP.)
The module also focuses on the following genes and/or proteins:
transforming growth factor beta receptor 2, insulin like growth factor 1 receptor, EPH receptor A2, epidermal growth factor receptor, platelet derived growth factor receptor alpha, erb-b2 receptor tyrosine kinase 2, ret proto-oncogene, fms related receptor tyrosine kinase 4.

Below are related articles:

Title: Inhibition of DDR2 reveals kinase-dependent and kinase-independent functions.
Text: Tyrosine kinase signaling plays a crucial role in fibrosis. DDR2, a receptor tyrosine kinase, regulates fibroblast apoptosis and fibrosis progression. Inhibitors like WRG-28 block DDR2 activation, reducing fibrosis severity.

Title: Oncogene-driven receptor tyrosine kinase expression in cancer.
Text: Overexpression of receptor tyrosine kinases like DDR2, RET, and PDGFRA correlates with poor survival in neuroblastoma. The multikinase inhibitor sitravatinib effectively targets DDR2, RET, and PDGFRA, reducing tumor proliferation and metastasis.

## **Example Output**
{
  "module_name": "Receptor Tyrosine Kinase Signaling in Fibrosis and Cancer",
  "summary": "This module focuses on the role of receptor tyrosine kinases (RTKs) in fibrosis and cancer, particularly emphasizing pathways related to transmembrane receptor protein kinase activity. Key genes involved include DDR2, PDGFRA, RET, EGFR, and EPH receptors, which regulate cellular signaling in fibroblasts, neuroblastoma, and other disease contexts. The studies highlight DDR2’s central role in fibrosis through collagen autocrine signaling and its involvement in cancer progression. Inhibitors like sitravatinib and WRG-28 have shown promise in targeting DDR2 and other RTKs, affecting tumor growth and fibrosis progression. Additionally, DDR2-mediated PIK3C2α signaling has been implicated in fibroblast function and fibrosis, presenting a potential therapeutic target."
  "confidence_score": "0.95"
}

## **Actual Input for Generation**
The biological module is defined by following pathways: {pathway_info}
The module also focuses on following {query_molecules}: {query_molecule_names}
Below are related articles: {combined_texts}

## **Final Output**
Please provide your response in **JSON format**, strictly following this structure:

{
  "module_name": "Your concise module name here",
  "summary": "A detailed, literature-integrated explanation of the module's function, key {query_molecules}, and relevant biological processes.",
  "confidence_score": "A value between 0.00 and 1.00, based on {query_molecules}-pathway association and overall functional coherence."
}
When generating your answer, pay attention to the following:

module_name should be concise and descriptive, highlighting the core biological process or disease relevance.
summary should include:
The central biological function or process.
Key {query_molecules} and their roles within the identified pathways.
Insights from the provided articles or research (if provided).
confidence_score should reflect how well the {query_molecules} and pathways associated with each other and support a coherent biological process, with 1.00 indicating very high confidence and 0.00 indicating minimal or no overlap.
