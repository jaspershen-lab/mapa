% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_llm_module_embedding_database.R
\name{query_abstracts_and_titles_by_pubmed_ids}
\alias{query_abstracts_and_titles_by_pubmed_ids}
\title{Retrieve Titles and Abstracts by PubMed IDs}
\usage{
query_abstracts_and_titles_by_pubmed_ids(pmid_list, retries = 2, pause = 1)
}
\arguments{
\item{retries}{An integer specifying the number of retry attempts for failed queries (default is 2).}

\item{pause}{A numeric value specifying the time in seconds to pause between retries (default is 1).}

\item{PID_list}{A character vector containing PubMed IDs to query.}
}
\value{
A named list where each element corresponds to a PubMed ID. Each element is a list containing:
\item{title}{The retrieved title for the PubMed ID. If retrieval fails, an error message is included.}
\item{abstract}{The retrieved abstract for the PubMed ID. If retrieval fails, an error message is included.}
}
\description{
Queries PubMed for titles and abstracts using a list of PubMed IDs (PIDs).
}
\details{
This function performs the following steps for each PubMed ID:
\itemize{
\item Queries PubMed for the title and abstract.
\item Handles errors and retries the query up to \code{retries} times in case of failure.
\item Returns a default error message if all retry attempts fail.
}

The function ensures robust error handling and provides meaningful messages in case of retrieval failures.
}
\examples{
\dontrun{
# Example: Query titles and abstracts for a list of PubMed IDs
PID_list <- c("12345678", "23456789")
results <- query_abstracts_and_titles_by_pubmed_ids(PID_list)
print(results)
}


}
\author{
Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
}
\keyword{internal}
