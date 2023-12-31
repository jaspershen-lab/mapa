% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/17-llm_interpretation.R
\name{interpret_pathways}
\alias{interpret_pathways}
\title{Interpret Pathways in Functional Modules}
\usage{
interpret_pathways(
  object,
  disease = "pregnancy",
  p.adjust.cutoff = 0.05,
  count.cutoff = 5,
  top_n = 3
)
}
\arguments{
\item{object}{An object of class "functional_module". This object contains the functional module data that needs to be interpreted.}

\item{disease}{A character string specifying the disease to focus on. Default is "pregnancy".}

\item{p.adjust.cutoff}{A numeric value representing the p-value adjustment cutoff for filtering functional modules. Default is 0.05.}

\item{count.cutoff}{A numeric value representing the count cutoff for filtering functional modules. Default is 5.}

\item{top_n}{An integer indicating the number of top pathways to return. Default is 3.}
}
\value{
A character string or list containing the interpreted pathway information. If no enriched functional modules are found,
a message is returned indicating that no enriched functional modules are found and suggesting different cutoffs.
}
\description{
This function interprets the pathways in functional modules, particularly focused on a specific disease.
It filters and arranges functional modules based on p-value adjustment and count cutoffs,
then retrieves and formats information about the top pathways related to the disease.
}
