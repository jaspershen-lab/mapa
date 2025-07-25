% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_llm_module_embedding_database.R
\name{embedding_single_pdf}
\alias{embedding_single_pdf}
\title{Process Single PDF for Embedding}
\usage{
embedding_single_pdf(
  pdf_path,
  embedding_model = "text-embedding-3-small",
  api_key,
  api_provider = "openai"
)
}
\arguments{
\item{pdf_path}{A character string specifying the path to the PDF file.}

\item{embedding_model}{A string specifying the embedding model to use (default is \code{"text-embedding-3-small"}).}

\item{api_key}{A character string containing the API key for the embedding service.}

\item{api_provider}{A string indicating the API provider, either \code{"openai"}, \code{"gemini"}, or \code{"siliconflow"} (default is \code{"openai"}).}
}
\value{
A data frame containing the following columns:
\item{paper_title}{The title of the PDF (derived from the file name).}
\item{chunks}{The extracted text chunks from the PDF.}
\item{embedding_vector}{A list of embedding vectors for each chunk.}
}
\description{
Processes a single PDF file by extracting text, splitting it into manageable chunks, and generating embeddings for each chunk.
}
\details{
This function performs the following:
\enumerate{
\item Extracts text from the PDF using
\code{\link[pdftools]{pdf_text}}.
\item Splits the text into smaller chunks using
\code{\link{split_into_chunks}}.
\item Generates embeddings for each chunk in parallel with
\code{\link[parallel]{parLapply}} (Windows/macOS) or
\code{\link[pbmcapply]{pbmclapply}} (Linux).
}
The resulting data frame contains the title of the PDF, the processed chunks, and their corresponding embedding vectors.
}
\examples{
\dontrun{
# Example: Process a single PDF to generate embeddings
pdf_path <- "example.pdf"
api_key <- "your_api_key"
embedding_data <- embedding_single_pdf(pdf_path, api_key)
print(embedding_data)
}
}
\author{
Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
}
\keyword{internal}
