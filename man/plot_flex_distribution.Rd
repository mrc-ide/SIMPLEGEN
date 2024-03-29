% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_flex_distribution}
\alias{plot_flex_distribution}
\title{Generic function for plotting flexible distributions}
\usage{
plot_flex_distribution(project, name = "duration_acute")
}
\arguments{
\item{project}{a SIMPLEGEN project, as produced by the
\code{simplegen_project()} function.}

\item{name}{the name of the distribution to be plotted, which must be one of:
\itemize{
  \item \code{"duration_acute"}
  \item \code{"duration_chronic"}
  \item \code{"time_treatment"}
  \item \code{"duration_prophylactic"}
  \item \code{"prob_infection"}
  \item \code{"prob_acute"}
  \item \code{"prob_AC"}
  \item \code{"detectability_microscopy_acute"}
  \item \code{"detectability_microscopy_chronic"}
  \item \code{"detectability_PCR_acute"}
  \item \code{"detectability_PCR_chronic"}
  \item \code{"infectivity_acute"}
  \item \code{"infectivity_chronic"}
}}
}
\description{
Many distributions within the inbuilt SIMPLEGEN transmission
  model can be defined flexibly as either a vector, or a list of vectors.
  This function plots these distributions.
}
