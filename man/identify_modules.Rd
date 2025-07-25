% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10_merge_pathways.R
\name{identify_modules}
\alias{identify_modules}
\title{Identify Modules in Similarity Matrix}
\usage{
identify_modules(
  sim_matrix,
  query_type = c("gene", "metabolite"),
  analysis_type = c("enrich_pathway", "do_gsea"),
  result,
  database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
  sim.cutoff = 0.5,
  save_to_local = TRUE,
  path = ""
)
}
\arguments{
\item{sim_matrix}{A data frame containing the similarity matrix with columns \code{name1}, \code{name2}, and \code{sim} representing the similarity between entities.}

\item{query_type}{Character, the category of biological entity to query ("gene", "metabolite") for merging pathway enrichment result.}

\item{analysis_type}{Character. Type of analysis to perform: either \code{"enrich_pathway"} or \code{"do_gsea"}. Default is \code{"enrich_pathway"}.}

\item{result}{A data frame containing the enrichment analysis results, including columns like \code{ID}, \code{Description}, \code{p_adjust}, and other relevant data.}

\item{database}{Character. The database from which the enrichment results were obtained (\code{go}, \code{kegg}, \code{reactome}, \code{hmdb}, \code{metkegg}).}

\item{sim.cutoff}{Numeric. The similarity cutoff value used to filter the edges in the similarity matrix. Default is \code{0.5}.}

\item{save_to_local}{Logical. Whether to save the resulting data to local files. Default is \code{TRUE}.}

\item{path}{Character. The directory path where intermediate results will be saved, if \code{save_to_local = TRUE}. Default is an empty string (current working directory).}
}
\value{
A list containing:
\item{graph_data}{A tidygraph object representing the network with nodes and edges.}
\item{module_result}{A data frame with the identified modules and their associated information, including pathway descriptions and p-values.}
\item{result_with_module}{A data frame with the original result data enriched with module information.}
}
\description{
This function identifies modules (clusters) in a given similarity matrix based on pathway enrichment or gene set enrichment analysis (GSEA). It constructs a network graph using similarity and result data, applies clustering algorithms, and optionally saves the results to a specified path.
}
\details{
The function first constructs a graph from the similarity matrix and result data, then applies clustering to identify modules. For pathway enrichment, the function organizes and processes the result data for easier interpretation. If \code{save_to_local} is \code{TRUE}, it saves intermediate results in the specified \code{path}.
}
\examples{
\dontrun{
sim_matrix <- data.frame(
  name1 = c("A", "B", "C"),
  name2 = c("B", "C", "A"),
  sim = c(0.6, 0.7, 0.8)
)
result <- data.frame(
  ID = c("P1", "P2", "P3"),
  p_adjust = c(0.01, 0.05, 0.03),
  Description = c("Pathway 1", "Pathway 2", "Pathway 3")
)
modules <- identify_modules(
  sim_matrix,
  "enrich_pathway",
  result,
  sim.cutoff = 0.5,
  save_to_local = FALSE
)
}

}
