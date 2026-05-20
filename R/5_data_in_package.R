#' Example Dataset for Over-Representation Analysis
#'
#' Significantly downregulated proteins from the muscle of aging mice
#' (6 vs 30 months, male C57BL/6). Contains 66 proteins with |log2FC| >= 0.5
#' and FDR < 0.05, suitable for pathway enrichment analysis.
#'
#' @format A data frame with 66 rows and 3 columns:
#' \describe{
#'   \item{symbol}{Gene symbols}
#'   \item{log2FC (6 vs 30mo)}{Log2 fold changes (all negative)}
#'   \item{FDR (6 vs 30mo)}{False discovery rates (all < 0.05)}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
#' @examples
#' data(example_ora_data)
#' head(example_ora_data)
#'
"example_ora_data"

#' Example Dataset for Gene Set Enrichment Analysis
#'
#' Complete proteomics dataset from liver of aging mice (6 vs 30 months,
#' male C57BL/6). Contains 5,290 proteins with fold changes and adjusted
#' p-values, suitable for gene set enrichment analysis (GSEA).
#'
#' @format A data frame with 5,290 rows and 3 columns:
#' \describe{
#'   \item{symbol}{Gene symbols}
#'   \item{fc}{Fold changes (6 vs 30 months)}
#'   \item{p_value_adjust}{Adjusted p-values}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
#' @examples
#' data(example_gsea_data)
#' head(example_gsea_data)
#'
"example_gsea_data"

#' Example Metabolomics Dataset for Pathway Enrichment Analysis
#'
#' Significantly altered metabolites identified through untargeted metabolomics
#' analysis. Contains 106 metabolite features with KEGG IDs and statistical
#' measures, suitable for metabolite pathway enrichment analysis.
#'
#' @format A data frame with 106 rows and 4 columns:
#' \describe{
#'   \item{variable_id}{Metabolite feature identifiers in the format
#'                     "M[mass]T[retention_time]_[ionization_mode]"}
#'   \item{keggid}{KEGG compound identifiers (e.g., "C05466");
#'                 contains NA for unidentified metabolites}
#'   \item{fdr}{False discovery rate adjusted p-values from differential
#'              abundance testing}
#'   \item{score}{Abundance fold change scores or effect sizes}
#' }
#'
#' @examples
#' data(example_met_data)
#' head(example_met_data)
#'
"example_met_data"

#' Multi-Omics Demo Dataset: Transcriptomics
#'
#' Differentially expressed genes from human transcriptomics data, provided
#' as a demo for multi-omics pathway enrichment analysis.
#' Contains gene symbols with fold changes and adjusted p-values.
#'
#' @format A data frame with gene symbols (human, symbol IDs) and associated
#' statistical measures from differential expression analysis.
#'
#' @examples
#' data(demo_mo_T_data)
#' head(demo_mo_T_data)
#'
"demo_mo_T_data"

#' Multi-Omics Demo Dataset: Proteomics
#'
#' Differentially expressed proteins from human proteomics data, provided
#' as a demo for multi-omics pathway enrichment analysis.
#' Contains protein symbols with fold changes and adjusted p-values.
#'
#' @format A data frame with protein symbols (human, symbol IDs) and associated
#' statistical measures from differential expression analysis.
#'
#' @examples
#' data(demo_mo_P_data)
#' head(demo_mo_P_data)
#'
"demo_mo_P_data"

#' Multi-Omics Demo Dataset: Metabolomics
#'
#' Significantly altered metabolites from human metabolomics data, provided
#' as a demo for multi-omics pathway enrichment analysis.
#' Contains KEGG compound IDs with statistical measures.
#'
#' @format A data frame with KEGG compound identifiers (human/hsa) and
#' associated statistical measures from differential abundance analysis.
#'
#' @examples
#' data(demo_mo_M_data)
#' head(demo_mo_M_data)
#'
"demo_mo_M_data"
