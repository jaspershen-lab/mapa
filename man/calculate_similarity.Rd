% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_llm_module_RAG_strategy.R
\name{calculate_similarity}
\alias{calculate_similarity}
\title{Calculate Similarity Between Embeddings}
\usage{
calculate_similarity(target_embeddings_list, module_embedding)
}
\arguments{
\item{target_embeddings_list}{A matrix or data frame where each row represents an embedding.}

\item{module_embedding}{A numeric vector representing the module embedding to compare against.}
}
\value{
A numeric vector containing similarity scores for each row in \code{target_embeddings_list}.
}
\description{
Calculates the cosine similarity between a target embedding and a list of embeddings.
}
\examples{
\dontrun{
# Example: Calculate similarity scores
target_embeddings <- matrix(runif(100), nrow = 10) # 10 embeddings with random values
module_embedding <- runif(10) # A random module embedding
similarity_scores <- calculate_similarity(target_embeddings, module_embedding)
}

}
\author{
Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
}
\keyword{internal}
