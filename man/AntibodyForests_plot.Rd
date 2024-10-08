% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AntibodyForests_plot.R
\name{AntibodyForests_plot}
\alias{AntibodyForests_plot}
\title{Custom plots for trees/networks created with AntibodyForests}
\usage{
AntibodyForests_plot(
  network.list,
  graph.type,
  node.color,
  node.label,
  node.shape,
  node.size,
  max.node.size,
  node.scale.factor,
  edge.length,
  edge.width,
  path.list,
  specific.node.colors,
  specific.node.shapes,
  specific.edge.colors,
  color.by.majority,
  cell.color,
  specific.cell.colors,
  cell.size,
  network.layout,
  save.pdf,
  save.dir,
  show.legend,
  color.gradient
)
}
\arguments{
\item{network.list}{nested list of igraph objects, as obtained from the AntibodyForests function. input[[1]][[2]] represents the igraph object for the first sample, second clonotype, if the normal AntibodyForests parameters are used.}

\item{graph.type}{string - graph to be plotted - 'tree' for the normal tree plots, 'heterogeneous' for the single-cell grahs (call AntibodyForests_heterogeneous before), 'phylo' to plot the phylo objects (call AntibodyForests_phylo before), 'dynamic' for the dynamic/temporal graphs.}

\item{node.color}{string specifying the name of the igraph vertex attribute (or the original vgm[[1]]/VDJ column name used in the node.features parameter of AntibodyForests) to be used for coloring the nodes.
If the node.shape is also 'pie', the resulting nodes will be a pie chart of the per-node values denoted by the node.color parameter (if there are multiple different values per node - e.g., a node denotes a sequence, which can be further traced back to multiple cells with different barcodes and with potentially different transcriptomic clusters - multiple feature values per node).
If the node.color parameter is NULL, the default node color will be '#FFCC00' (yellow) for sequence nodes, lighter gray for intermediate/inferred nodes, and darker gray for germline nodes.}

\item{node.label}{string - 'cells' to label the nodes by the number of cells with that specific sequence, 'rank' for cell count-based ranking of the nodes. If NULL, will not add number labels to the sequence nodes.}

\item{node.shape}{string specifying the the name of the igraph vertex attribute (or the original vgm[[1]]/VDJ column name used in the node.features parameter of AntibodyForests) to be used for node shapes. Shapes will be assigned per unique value from the values of this column/vertex attribute. There is a maximum of 7 unique shapes to be chosen from, therefore node.shapes should be used for features with less than 7 unique values.
If the node.shape parameter is null, the node shape for all nodes will default to 'circle'}

\item{node.size}{string denoting either a specific method of scaling the node sizes of the input graph or a specific vertex attribute of numeric values (added via the node.features parameter of the AntibodyForests from the original vgm[[1]]/VDJ dataframe) to be used for node sizes.
If NULL, then the sizes will be equal to node.scale.factor * 1. If scaleByEigen, node sizes will be scaled by the eigenvector centrality of each node; scaleByCloseness - closeness centrality; scaleByBetweenness - betweenness centrality; scaleByExpansion - by the sequence frequency of each node, as originally calculated by the AntibodyForests function.}

\item{max.node.size}{Maximum size of any given node}

\item{node.scale.factor}{integer to further refine the size of each node. Each vertex size will be multiplied by scaling factor.}

\item{edge.length}{either NULL - edge lengths are constant, scaleByWeight - edge lengths will be scaled by the edge weight attribute (the string distance between pairs of sequences/nodes, as calculated by the AntibodyForests function) - larger weights/string distances = longer edges/nodes further apart, or a specific igraph edge attribute name.}

\item{edge.width}{either NULL - edge widths are constant, scaleByWeight - edge widths will be scaled by the edge weight attribute (the string distance between pairs of sequences/nodes, as calculated by the AntibodyForests function) - larger weights/string distances = thinner edges), or a specific igraph edge attribute name.}

\item{path.list}{named list of igraph paths, as obtained from the AntibodyForests_metrics function.}

\item{specific.node.colors}{named list of colors to be used for the node.color parameter. If NULL, colors will be automatically added for each unique value.
For example, if specific.node.colors=list('yes'='blue', 'no'='red'), then the nodes labeled as 'yes' will be colored blue, the others red.}

\item{specific.node.shapes}{named list of node shapes to be used for the node.shapes parameter. Must be shapes compatible with igraph objects, use setdiff(igraph::shapes(), "") to get a list of possible values.
For example, if specific.node.shapes=list('yes'='circle', 'no'='square'), then the nodes labele as 'yes' will be circles, the others squares.}

\item{specific.edge.colors}{named list of edge colors to be used for the edge.colors parameter. The names should be the path metrics names obtained from the nested paths list, from AntibodyForests_metrics.
For example, if specific.edge.colors=list('longest.path.weighted'='blue', 'shortest.path.unweighted'='red'), the longest weighted paths obtained from AbtibodyForests will be colored blue for each igraph object, the rest will be red.}

\item{color.by.majority}{boolean - if T, will color the entire network (all nodes) by the dominant/most frequent node feature, as specified in the node.color parameter.}

\item{cell.color}{string - cell feature column denoting the cell colors - as denoted by the node.features parameter when calling AntibodyForests_heterogeneous.}

\item{specific.cell.colors}{named list of cell colors and their features (e.g., for Seurat clusters: list(1 = 'red', 2 = 'blue')). Optional (will auto search for unique colors per feature).}

\item{cell.size}{integer denoting the size of the cell nodes.}

\item{network.layout}{either NULL - will default to the automatic igraph::layout_nicely(), 'fr' - igraph::layout_with_fr() - for fully connected graphs/graphs with defined connected components, 'tree' - for tree graphs, as obtained from AntibodyForests with network.algorithm='tree'.}

\item{save.pdf}{boolean - if T, plots will be automatically saved to pdf, in the current working directory; F - normal output of the function (plot-ready igraph object and specific layout). New folders will be created for each sample of the input nested list of igraph objects.}

\item{save.dir}{path to the directory oin which the network PDFs will be saved.}

\item{show.legend}{boolean - whether the legend should be showed in the resulting plots}

\item{color.gradient}{string - defualt: NULL - the feature whose value will be plot as a color gradient}
}
\value{
nested list of plot-ready AntibodyForests objects. Can also save the plots as a PDF file.
}
\description{
AntibodyForests_plot takes the input of AntibodyForests and outputs a list of plot-ready graphs (inside the AntibodyForests object) to be further used with plot(). Plots can also be automicatically saved to pdf via the save.pdf parameter. The resulting igraph object have their node/edge colors/shapes/sizes added following the specific parameters in the AntibodyForests_plot function.
}
\examples{
\dontrun{
AntibodyForests_plot(graphs, node.color='clonotype_id',
node.size='scaleByExpansion', network.layout='tree',
save.pdf=T)
}
}
\seealso{
AntibodyForests, AntibodyForests_metrics
}
