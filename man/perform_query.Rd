% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_llm_module_online_retrieval.R
\name{perform_query}
\alias{perform_query}
\title{Perform PubMed Query with Terms and Pathway Information}
\usage{
perform_query(query_terms, pathway_query, years, retmax, chunk_size)
}
\arguments{
\item{query_terms}{A character vector of query terms (e.g., gene symbols, gene names, metabolite names).}

\item{pathway_query}{A character string representing the pathway part of the query.}

\item{years}{An integer specifying how many years to look back in the search.}

\item{retmax}{An integer specifying the maximum number of results to retrieve.}

\item{chunk_size}{An integer specifying the size of query chunks.}
}
\value{
A character vector of unique PubMed IDs retrieved from the search.
}
\description{
This internal function executes a PubMed search using query terms and pathway information,
with fallback strategies for handling failed queries. It implements a three-level search
strategy: first attempting a full query with all terms, then breaking into chunks if that
fails, and finally trying individual terms if chunk queries also fail.
}
\author{
Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}

Yifei Ge \email{yifeii.ge@outlook.com}
}
\keyword{internal}
