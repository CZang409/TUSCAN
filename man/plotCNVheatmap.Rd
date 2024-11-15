% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCNVheatmap.R
\name{plotCNVheatmap}
\alias{plotCNVheatmap}
\title{Plot the copy number matrix as a heatmap, with spots as rows and genes as columns, ordered according to chromosome}
\usage{
plotCNVheatmap(
  object,
  color_belt = NULL,
  output_format = c("png", "pdf"),
  title = "heatmap"
)
}
\arguments{
\item{object}{TUSCAN object.}

\item{color_belt}{Optional color palette for plotting.}

\item{title}{Plot title.}
}
\value{
Create a PNG or PDF image file.
}
\description{
Visualize the CNV profile by heatmap.
}
