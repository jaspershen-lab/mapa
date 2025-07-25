% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/13_data_visualization_plot_module_info.R
\name{plot_module_info}
\alias{plot_module_info}
\title{Plot Module Information for Enrichment Analysis}
\usage{
plot_module_info(
  object,
  level = c("module", "functional_module"),
  database = c("go", "kegg", "reactome", "hmdb", "metkegg"),
  database_color = c(GO = "#eeca40", KEGG = "#fd7541", Reactome = "#23b9c7", HMDB =
    "#7998ad"),
  module_id,
  llm_text = FALSE
)
}
\arguments{
\item{object}{An object of class "functional_module" containing enrichment results.}

\item{level}{A character string specifying the level of biological organization, either "module" or "functional_module".
For results generated by get_bioembedsim() and merge_pathways_bioembedsim(),  please use "functional_module".}

\item{database}{A character string specifying the source database ("go", "kegg", "reactome", "hmdb", or "metkegg"). Only required when level is "module".}

\item{database_color}{A named vector specifying colors for each database. Default colors are provided.}

\item{module_id}{A single identifier specifying the module of interest.}

\item{llm_text}{Logical. Whether to use LLM-generated module names in the visualization. Default is FALSE.}
}
\value{
A list containing three ggplot objects: network plot (\code{network}), bar plot (\code{barplot}), and word cloud (\code{wordcloud}).
}
\description{
This function generates various visualizations including a network plot, bar plot,
and word cloud to provide detailed insights into the enriched functional or biological modules.
}
\examples{
\dontrun{
result <- plot_module_info(object = myObject,
                          level = "module",
                          database = "go",
                          module_id = "M123")
}
}
\author{
Xiaotao Shen \href{mailto:shenxt1990@outlook.com}{shenxt1990@outlook.com}

Yifei Ge \href{mailto:yifeii.ge@outlook.com}{yifeii.ge@outlook.com}
}
