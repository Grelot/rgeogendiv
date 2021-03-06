% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preliminary_analysis.R
\name{preliminary_analysis}
\alias{preliminary_analysis}
\title{Preliminary Analysis}
\usage{
preliminary_analysis(bold_res)
}
\arguments{
\item{bold_res}{a data.frame with same fields than resBold$data$markercode and a supplementary column with DNA sequences as string}
}
\value{
a ggplot object
}
\description{
From bold prepared query, generates a figure as ggplot object
Distribution of markercode
Accumulation of georeferenced GenBank accessions and frequency of submissions with mostly georeferenced data
}
