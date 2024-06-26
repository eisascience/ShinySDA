% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fxs.R
\name{.run_tSN_scores}
\alias{.run_tSN_scores}
\title{This function to run tsne}
\usage{
.run_tSN_scores(
  SDAres,
  save_path,
  tsnepp,
  tSNE_n.iter,
  num_threads = 8,
  check_duplicates = F,
  selectComps = NULL
)
}
\arguments{
\item{SDAres}{SDA score matrix where cols are the comps and rows are cells}

\item{save_path}{path to save the tsne rds obj}

\item{tsnepp}{tSNE preplexity}

\item{tSNE_n.iter}{N iter for tSNE}

\item{num_threads}{number of threads}

\item{check_duplicates}{check for dup rows}

\item{selectComps}{select SDA compos to run the tSNE}
}
\value{
A list of two vector text of enrichments, one for pos one for neg scored cells
}
\description{
This function to run tsne
}
