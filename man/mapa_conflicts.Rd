% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4-conflicts.R
\docType{methods}
\name{mapa_conflicts}
\alias{mapa_conflicts}
\title{Conflicts between the mapa and other packages}
\usage{
mapa_conflicts()
}
\value{
mapa conflicts
}
\description{
This function lists all the conflicts between packages in the mapa
and other packages that you have loaded.
}
\details{
There are four conflicts that are deliberately ignored: \code{intersect},
\code{union}, \code{setequal}, and \code{setdiff} from dplyr. These functions
make the base equivalents generic, so shouldn't negatively affect any
existing code.
}
\examples{
mapa_conflicts()
}
