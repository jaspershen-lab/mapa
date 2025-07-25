% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_llm_module_output_generation.R
\name{single_module_generation}
\alias{single_module_generation}
\title{Single Module Generation}
\usage{
single_module_generation(
  module_related_paper,
  module_info,
  phenotype = NULL,
  model = "gpt-4o-mini-2024-07-18",
  api_key,
  output_prompt = TRUE,
  api_provider = "openai",
  thinkingBudget = 0
)
}
\arguments{
\item{module_related_paper}{A list of related papers with titles and cleaned texts for the module.}

\item{module_info}{A list containing pathway names and gene symbols relevant to the module.}

\item{phenotype}{Character string. Phenotype or disease to focus on. Default is NULL.}

\item{model}{A string specifying the GPT model to use. Default is \code{"gpt-4o-mini-2024-07-18"}.}

\item{api_key}{A string containing the API key required to access the AI API.}

\item{output_prompt}{Logical. Whether to output prompt in final annotation result. Default is TRUE.}

\item{api_provider}{A string indicating the API provider, either \code{"openai"}, \code{"gemini"}, or \code{"siliconflow"} (default is \code{"openai"}).}

\item{thinkingBudget}{An integer for the "thinking budget" parameter specific to the Gemini API (default is \code{0}).}
}
\value{
A list containing two elements: \code{module_name} (the generated biological module name)
and \code{summary} (the research summary).
}
\description{
This internal function generates a biological module name and a research summary
for a specific module based on provided pathway, gene, and related article information.
It interacts with an AI API to perform text analysis and ensure the response is in JSON format.
}
\author{
Feifan Zhang \email{FEIFAN004@e.ntu.edu.sg}
}
\keyword{internal}
