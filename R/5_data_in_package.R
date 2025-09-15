#' Demo Dataset for Over-Representation Analysis
#'
#' Significantly downregulated proteins from the muscle of aging mice
#' (6 vs 30 months, male C57BL/6). Contains 66 proteins with |log2FC| ≥ 0.5
#' and FDR < 0.05, suitable for pathway enrichment analysis.
#'
#' @format A tibble with 66 rows and 3 columns:
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
#' data(demo_data_ora)
#' head(demo_data_ora)
#'
#' # Extract gene symbols for ORA
#' genes <- demo_data_ora$symbol
#'
"demo_data_ora"

#' Demo Dataset for Gene Set Enrichment Analysis
#'
#' Complete proteomics dataset from liver of aging mice (6 vs 30 months,
#' male C57BL/6). Contains 5,290 proteins with fold changes and adjusted
#' p-values, suitable for gene set enrichment analysis (GSEA).
#'
#' @format A tibble with 5,290 rows and 3 columns:
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
#' data(demo_data_gsea)
#' head(demo_data_gsea)
#'
#' # Create ranked gene list for GSEA
#' ranked_genes <- setNames(demo_data_gsea$fc, demo_data_gsea$symbol)
#' ranked_genes <- sort(ranked_genes, decreasing = TRUE)
#'
"demo_data_gsea"

#' Demo Dataset for Metabolite Enrichment Analysis
#'
#' Significantly altered metabolites identified through untargeted metabolomics
#' analysis. Contains 106 metabolite features with their corresponding KEGG IDs,
#' FDR values, suitable for metabolite pathway enrichment
#' analysis.
#'
#' @format A tibble with 106 rows and 4 columns:
#' \describe{
#'   \item{variable_id}{Unique metabolite feature identifiers}
#'   \item{keggid}{KEGG compound identifiers}
#'   \item{fdr}{False discovery rates from statistical testing}
#'   \item{score}{fold-change scores}
#' }
#'
#'
#' @examples
#' data(demo_data_met)
#' head(demo_data_met)
#'
"demo_data_met"

#' Demo Expression Dataset for Relationship Heatmap Visualization
#'
#' Gene expression data from muscle tissue of aging mice (6 vs 30 months,
#' male C57BL/6). Contains 65 genes with expression values across 8 samples
#' (4 samples per age group). This dataset is designed to work with the
#' analysis results from \code{demo_data_ora} to create integrated relationship
#' heatmaps using \code{plot_relationship_heatmap()}.
#'
#' @format A tibble with 65 rows and 9 columns:
#' \describe{
#'   \item{id}{ENSEMBL gene identifiers}
#'   \item{6mo-1, 6mo-2, 6mo-3, 6mo-4}{Expression values for 6-month-old samples}
#'   \item{30mo-1, 30mo-2, 30mo-3, 30mo-4}{Expression values for 30-month-old samples}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
"ora_expression_dt"


#' Demo Expression Dataset for Relationship Heatmap Visualization
#'
#' Gene expression data from liver tissue of aging mice (6 vs 30 months,
#' male C57BL/6). Contains 5,167 genes with expression values across 8 samples
#' (4 samples per age group). This dataset is designed to work with the
#' analysis results from \code{demo_data_gsea} to create integrated relationship
#' heatmaps using \code{plot_relationship_heatmap()}.
#'
#' @format A tibble with 5,167 rows and 9 columns:
#' \describe{
#'   \item{id}{ENSEMBL gene identifiers}
#'   \item{6mo-1, 6mo-2, 6mo-3, 6mo-4}{Expression values for 6-month-old samples}
#'   \item{30mo-1, 30mo-2, 30mo-3, 30mo-4}{Expression values for 30-month-old samples}
#' }
#'
#' @source
#' Takasugi, M., et al. An atlas of the aging mouse proteome reveals the
#' features of age-related post-transcriptional dysregulation.
#' \emph{Nat Commun} \strong{15}, 8520 (2024).
#' \doi{10.1038/s41467-024-52845-x}
#'
"gsea_expression_dt"


