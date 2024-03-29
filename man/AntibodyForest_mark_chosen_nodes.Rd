% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AntibodyForest_mark_chosen_nodes.R
\name{AntibodyForest_mark_chosen_nodes}
\alias{AntibodyForest_mark_chosen_nodes}
\title{From PnP Sequences to ABForest plots (targets the plot.ready slot)
Takes identifiers from extracted PnP sequences and marks then the respective nodes in the ABForest plot}
\usage{
AntibodyForest_mark_chosen_nodes(
  AntibodyForest.networks,
  Seq_output,
  Frame.color,
  Plot.Suppress,
  specific.node.colors
)
}
\arguments{
\item{AntibodyForest.networks}{List of samples (sample_ids need to be unique), which is a list of clonotypes/igraph objects (AbForest output containing mixcr alignemnt results in network/igraph objects)}

\item{Seq_output}{Dataframe with extracted sequences. Output of VDJ_assemble_for_PnP_github_MV or correspondingly the AntibodyForest_network_to_pnp function}

\item{Frame.color}{Character. Specifying with which color one wants the nodes to be marked.}

\item{Plot.Suppress}{Boolean. Whether the plots should be plotted on function call.}

\item{specific.node.colors}{Vector of strings. Colors for the legend for the specific isotypes. Normally already defined.}
}
\value{
A data.frame containing PnP sequences of selected sequences of networks.
}
\description{
From PnP Sequences to ABForest plots (targets the plot.ready slot)
Takes identifiers from extracted PnP sequences and marks then the respective nodes in the ABForest plot
}
\examples{
\dontrun{
what
To return results for a non default column
physiologically_relevant <- to_something(VDJ = VDJ_GEX_matrix.output[[1]]
, column.to.plot = "VJ_jgene", normalization factor = 20)
}
}
