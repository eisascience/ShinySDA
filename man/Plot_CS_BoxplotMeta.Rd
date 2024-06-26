% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fxs.R
\name{Plot_CS_BoxplotMeta}
\alias{Plot_CS_BoxplotMeta}
\title{Plot_CS_BoxplotMeta}
\usage{
Plot_CS_BoxplotMeta(
  SDAScores = NULL,
  MetaDF,
  MetaSelect = "Population",
  KeepCols = NULL,
  scaleMean0 = T,
  col_vector
)
}
\arguments{
\item{SDAScores}{SDA score matrix}

\item{MetaDF}{Metadata Dataframe}

\item{MetaSelect}{A factored feature of Metadata}

\item{KeepCols}{a numeric vector of features to keep, default NULL = all}

\item{scaleMean0}{Boolean default T to scale each compoent to mean 0}
}
\value{
A pheatmap cor heatmap
}
\description{
This function plots Box plots of the cell scores relative to a selected meta feature
}
