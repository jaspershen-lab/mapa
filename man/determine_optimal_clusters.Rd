% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_1_determine_optimal_clutsers.R
\name{determine_optimal_clusters}
\alias{determine_optimal_clusters}
\title{Determine Optimal Clustering Parameters for Pathway Enrichment Results}
\usage{
determine_optimal_clusters(object, ...)
}
\arguments{
\item{object}{An object containing pathway enrichment results. This can be a \code{list}
with \code{sim_matrix} and \code{enriched_pathway} elements, or a \code{functional_module} S4 object.}

\item{...}{Additional arguments passed to specific methods.}
}
\value{
A list containing three elements:
\item{cluster_result}{A data frame in long format with the complete
evaluation results, including method, cutoff, metric type (modularity or
silhouette), and the calculated metric values.}
\item{evaluation_plot}{A comprehensive \code{ggplot} object with a heatmap-like
visualization showing performance across methods and cutoffs. This plot
makes optimal parameter selection visually intuitive, with color representing
modularity and point size representing the silhouette score.}
\item{best_combination}{A data frame identifying the best-performing
method and cutoff combination for each metric (modularity and silhouette).
Ties for the top score will result in multiple rows per metric.}
}
\description{
This function systematically evaluates various clustering strategies to identify the
optimal parameters for grouping pathway enrichment results. It iterates through
different clustering methods, including network-based community detection,
distance-based hierarchical clustering, and a recursive divisive algorithm.
Each parameter combination is assessed using modularity and silhouette scores
to quantify cluster quality.
}
\details{
The function handles two types of input for the \code{object} parameter:
\enumerate{
\item A \code{list} that must contain a pre-computed \code{sim_matrix} (a similarity matrix)
and an \code{enriched_pathway} object. This is typically used for results from
biotext embedding.
\item An S4 \code{functional_module} object from the \code{mapa} package, which has been
processed by the \code{mapa::merge_pathways()} function. In this case,
the function will first calculate a Jaccard similarity matrix internally
based on the overlapping genes/metabolites between pathways before
proceeding with the evaluation.
}

Available clustering methods:
\itemize{
\item Network-based: "louvain", "walktrap", "infomap", "edge_betweenness", "fast_greedy",
"label_prop", "leading_eigen", "optimal"
\item Distance-based: "hierarchical" (supports various linkage methods)
\item Divisive: "binary_cut"
}
}
\examples{
\dontrun{
# Assuming 'enrichment_obj' is your S4 enrichment analysis result from mapa
# that has been processed by merge_pathways().
cluster_eval <- determine_optimal_clusters(enrichment_obj)
}

}
