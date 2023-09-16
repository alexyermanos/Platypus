#' S4 class for the AntibodyForests object.
#'
#' @slot tree igraph object of tree.
#' @slot sample_id string - sample ID.
#' @slot clonotype_id string - clonotype ID.
#' @slot plot_ready igraph object of tree with plotting config (node size, color, etc).
#' @slot heterogeneous igraph object of cell-cell and sequence-sequence graphs.
#' @slot reactivity igraph object of antibody-antigen graphs.
#' @slot dynamics list of igraph object - longitudinal graphs for a single lineage.
#' @slot metrics dataframe - igraph metrics for a tree (betweenness, evenness, etc).
#' @slot sequences list - sequences for a given tree.
#' @slot germline_sequence string - germline sequence.
#' @slot barcodes list - cell barcodes corresponding to a given tree.
#' @slot node_features dataframe - node features added when creating the tree.
#' @slot edge_list dataframe - edges and links between them.
#' @slot gex_list dataframe - cell-cell edges.
#' @slot paths igraph object - paths (shortest, longest) in a tree.
#' @slot node_transitions dataframe - feature transitions from node to node.
#' @slot adjacency_matrix Matrix object - adjacency matrix for a tree.
#' @slot phylo phylo object for phylogenetic plots/analyses.
#' @slot feature_names list - node features added when creating the tree.
#' @slot network_algorithm string - algorithm for tree/network inference.
#' @slot inferred igraph object - tree with intermediate nodes.
#' @slot permuted_transitions dataframe - node transitions permuted.
methods::setClass('AntibodyForests',
  slots = c(
    tree = 'ANY', #in
    sample_id = 'ANY', #in
    clonotype_id = 'ANY', #in
    plot_ready = 'ANY', #f call
    heterogeneous = 'ANY', #f call
    reactivity = 'ANY', #f call
    dynamic = 'ANY', #f call
    metrics = 'ANY', #f call
    sequences = 'ANY', #in
    germline_sequence = 'ANY', #to add
    barcodes = 'ANY', #in
    node_features = 'ANY', #in
    edge_list = 'ANY', #no/f call
    gex_list = 'ANY', #f call
    paths = 'ANY', #f call
    node_transitions = 'ANY', #f call
    adjacency_matrix = 'ANY', #no/f call
    phylo = 'ANY', #f call
    feature_names = 'ANY',
    network_algorithm = 'ANY',
    inferred = 'ANY',
    permuted_transitions = 'ANY'
  )
)
