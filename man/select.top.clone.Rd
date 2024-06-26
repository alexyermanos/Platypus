% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Echidna_other_functions.R
\name{select.top.clone}
\alias{select.top.clone}
\title{Get the index of top ranking clones.}
\usage{
select.top.clone(clonotypes, top.n)
}
\arguments{
\item{clonotypes}{The output "clonotypes" dataframe from simulation output.}

\item{top.n}{The top n abundant clones to be selected.}
}
\value{
a vector of indexes of top ranking clones
}
\description{
Get the index of top ranking clones.
}
