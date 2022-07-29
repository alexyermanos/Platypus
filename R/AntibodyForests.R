#' Infer B cell evolutionary networks and/or sequence similarity networks


#'@description Function to infer immune receptor evolutionary networks (trees) or complex/sequence-similarity networks from a Platypus VDJ object/VGM[[1]] - AntibodyForests objects will be created for each sample and each unique clonotype, if the default parameters are used.
#' Networks can be created in a tree-building manner (minimum spanning tree algorithm with custom tie solving methods), by linking sequences with a minimal string distance between them iteratively (and solving distance ties in a hierarchical way, with multiple resolve.ties parameters/configurations). Nodes in the network represent unique sequences per clonotype, edges are determined by the string distance between the nodes/sequences. Sequence types are dictated by the sequence.tyoe parameter, allowing networks to be built for most of the sequence types available in a Platypus VDJ object.
#' Networks can also be created by pruning edges from a fully connected network obtained from all sequences from a specific clonotype - complex similarity networks. Pruning can either be done via a distance threshold (prunes nodes too far apart), a node degree threshold (to prune nodes with a smaller degree/not well connected), or an expansion threshold (to prune nodes for sequences with low expansion/frequency).
#' Lastly, networks can be created by converting a phylogenetic tree inferred via different methods (neighbour-joining, maximum likelihood, maximum parsimony) into an igraph object,


#' @param VDJ VDJ or vgm[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param sequence.type string denoting the sequence types to create the networks for. 'cdr3.aa' - networks for amino-acid CDR3 sequences, 'cdr3.nt' - networks for nucleotide CDR3 sequences, 'VDJ.VJ.nt.trimmed' - full, trimmed VDJ-VJ sequences, as obtained when setting trin.and.align = T for VDJ_GEX_matrix(), 'VDJ.VJ.nt.raw' - full, raw VDJ-VJ sequences, 'VDJ.VJ.aa.mixcr' and 'VDJ.VJ.nt.mixcr' for the VDJ and VJ chains (nt or aa) as inferred by MIXCR, 'VDJ.aa.mixcr' and 'VDJ.nt.mixcr' for the VDJ chain inferred by MIXCR, 'VJ.aa.mixcr' and 'VJ.nt.mixcr' for the VJ chain inferred by MIXCR, 'VDJ.nt.trimmed' for the trimmed VDJ chain as nucleotides, 'VDJ.nt.raw' for the untrimmed VDJ chain as nucleotides, similarly for the VJ chain ('VJ.nt.trimmed' and 'VJ.nt.raw'), 'VDJ.cdr3s.aa' for the CDRH3 region as amino-acids, VDJ.cdr3s.nt' for the CDRH3 region as nucleotides, similarly for the CDRL3 regions ('VJ.cdr3s.aa', 'VJ.cdr3s.nt'), 'VDJ.aa' and 'VJ.aa' for the full VDJ/VJ sequence as amino-acids.
#' Defaults to 'VDJ.VJ.nt.trimmed'.
#' @param include.germline string or vector of strings, denoting the germline column(s) to be used (in the c('VDJ_germline', 'VJ_germline') order). 'trimmed.ref' - the networks will include a germline node, obtained by pasting the VDJ_trimmed_ref and VJ_trimmed_ref sequences for each clonotype, obtained by calling VDJ_call_MIXCR on VDJ. As reconstructed germlines as usually available for full VDJ.VJ.nt sequences, use this with sequence.type=VDJ.VJ.nt.trimmed. NULL will not include a germline.
#' @param network.algorithm string denoting the algorithm used in constructing the networks. 'tree' - will use a tree evolutionary inference algorithm: nodes denoting unique sequences per clonotype will be linked iteratively, as long as their string distance is the minimum. Use the resolve.ties parameter to further dictate the tree topology (when there are multiple ties in the minimum links).
#' 'prune' will create networks by pruning nodes from a fully connected networks. It must always be followed by a pruning method.
#' For example, 'prune.distance' will prune nodes with a larger string distance than the threshold specified in the pruning.threshold parameter.
#' 'prune.degree' will remove nodes with a lower degree than the threshold specified in pruning.threshold.
#' 'prune.expansion' will remove nodes with a lower sequence frequency/expansion than the threshold specified in pruning.threshold.
#' Multiple pruning methods can be used simultaneously, as long as multiple pruning.threshold is a vector - a threshold for each method.
#' For example, 'prune.distance.degree' with pruning.threshold=c(3,2) will remove edges of nodes with a string distance greater than 3, then will remove nodes with a degree smaller than 2. 'prune.distance.degree.expansion' with pruning.threshold=c(3,2,1) will also remove one-of sequences/nodes (with a single cell).
#' 'phylogenetic.tree.nj' will create phylogenetic (binary) trees using the neighbour-joining algorithm via ape::nj().
#' 'phylogenetic.tree.ml' will create phylogenetic (binary) trees using a maximum-likelihood algorithm from the phangorn package (phangorn::ml()).
#' 'phylogenetic.tree.mp' will create phylogenetic (binary) trees using a maximum-parsimony algorithm from the phangorn package (phangorn::mp()).
#' 'mst' will create undirected trees using the minimum spanning tree algorithm from igraph (without specific tie solving mechanisms).
#' 'global' is a custom option to easily create whole-repertoire/multi-repertoire similarity networks: it defaults to the 'prune.distance' option, while also changing some other parameters to ensure consistency (directed is set to F, include.germline is set to F, network.level is set to 'global')
#' @param directed boolean, whether networks obtained using network.algorithm='tree' should be directed (from the germline to the leaf nodes) or not. T - directed; F - undirected trees.
#' @param distance.calculation string, specifying the method for calculating string distances between the sequences. Must be compatible with the method parameter of the stringdist::stringdistmatrix() function. Will default to 'lv' for Levenshtein distances.
#' @param resolve.ties vector of strings, denoting the manner in which ties should be resolved when assembling trees via network.algorithm='tree'. Ties are defined as edges with the same weights=string distances (as determined by the distance matrix for the fully connected network) between nodes already added in the tree and nodes to be added in the tree.
#' There are multiple default and custom configurations for this parameter:
#' 'first' will pick the first edge from a pool of edges of equal string distance (between the sequences) - these are ordered based on each node's expansion/cell count (therefore 'first' will try to add the most expanded node first);
#' 'random' - resolve ties by picking random tied edges;
#' 'close.germline.edges' and 'far.germline.edges' - will prefer the nodes closer or farther to/from the germline, as a number of edges (unweighted) to be next integrated into the network;
#' 'close.germline.distance' and 'far.germline.distance' - picks nodes closer/farther to/from the germline, determined by the string distance;
#' 'close.germline.weighted' and 'far.germline.weighted' - picks edges with nodes closer/farther to/from the germline, as a weighted path from the germline to the most recent integrated node;
#' 'min.expansion' and 'max.expansion' - will pick the most/least expanded sequences;
#' 'min.degree' and 'max.degree' - picks the nodes with the minimum or maximum degree - a distance threshold must be specified in the pruning.threshold parameter, otherwise all nodes will have the same degree;
#' An additional custom configuration can be used: either min/max/specific feature value, tied to a specific feature column as defined in the node.features parameter, using '-'. For example, 'yes-OVA_binder' will select nodes that are OVA binders when resolving ties; for min and max, the node.feature column should be of numeric class.
#' If a vector is provided, ties will be resolved in a hierarchical manner: for example, if resolve.ties=c('max.expansion', 'close.germline.distance'), it will first try to pick nodes with a max expansion that were not added in the network (with edge ties to those already added), then those closer to the germline (minimum string distance).
#' As these two options do not always fully converge, meaning that there could be also expansion ties and distance ties between the nodes to be added, not just edge ties, a 'first' options is always added at the end of the hierarchical tie resolving algorithm, which always converges/picks a specific edge and resolves a tie. Moreover, the 'first' option is also added when only selecting a single option (still in the form of a vector - for e.g., c('min.expansion') turns automatically into c('min.expansion', 'first')).
#' @param connect.germline.to string defining how the germline should be connected for both the pruning and tree building algorithms. When network.algorithm='tree', two options are available: 'min.adjacent' - will first connect the nodes with the min string distance from the germline, then continue adding nodes and building the tree, and 'threshold.adjacent' - will connect nodes with a string distance value lower than the threshold defined in pruning.threshold.
#' As the pruning algorithm starts by pruning all connections out of the specified boundaries, irrespective if they are germline ones or not, the germline needs to be added at the end if it is removed. If not, then this option is ignored.
#' Thus, there are additional options for including the germline when building a network via the pruning algorithm: 'largest.connected.component' - connects the germline to the largest resulting connected component(s), 'all.connected.components' - to all connected components, 'all.connected.components.non.single' - does not connect the germline to single-node components; 'none' - germline is not connected; 'min.adjacent' - connects to the node(s) with the minimum string distance.
#' @param pruning.threshold vector of max size=3, specifying the thresholds for the pruning algorithm when network.algorithm includes 'prune' (as seen. 'prune' can be followed by either 'expansion', 'degree' or 'distance', or a combination of them - 'prune.distance.degree'). If we have 'prune.degree', we need to first specify a distance threshold (as 'prune.degree'='prune.distance.degree', otherwise the degree is the same for all nodes in a fully connected network).
#' See also network.algorithm.
#' For a direct use, if network.algorithms is set to prune.distance, set pruning.threshold to a single integer denoting the distance between two nodes for which edges will be pruned (equal or more).
#' @param remove.singletons integer - in the case of the pruning network algorithm or 'global' network algorithm, it denotes the minimum connected component node number threshold (for e.g., if remove.singletons = 3, it will remove all nodes from a graph that form a 3-node compnent or less: 2-node and singletons). If NULL, will keep all components (including singletons) in the complex similarity graph.
#' @param keep.largest.cc boolean - if T, will only keep the largest connected component in the similarity network (pruning network algorithm or 'global' option for network.algorithm). If F, will keep all components (including singletons unless removed via remove.singletons).
#' @param VDJ.VJ.1chain boolean, T - excludes cells with an aberrant number of VDJ or VJ chains; F - will be kept in for network inference.
#' @param node.features string or vector of strings, denoting the column name(s) in the original VDJ/VGM[[1]] from which the node features to be extracted. This is done by first pooling the cell cell_barcodes per unique sequence in the VDJ (for each clonotype), then adding the features of those cells.
#' @param filter.specific.features list with two elements - first one denotes the column/feature which you wish to filter on, second denotes the specific feature which HAS TO BE INCLUDED IN THE NETWORK (e.g., list('OVA_binder', 'yes') will result in networks that have at least 1 binder for OVA).
#' @param filter.na.features string or vector of strings, denoting the same column name(s) as specified in node.features. This will remove netowrks were ALL nodes values for the specific feature are equal to NA.
#' @param node.limits list of integers or NULLs. node.limits[[1]] determines the least amount of unique sequences to create a network for, otherwise a network will not be created. If node.limits[[1]] is NULL, then there is no lower bound for the number of sequences in each network. node.limits[[2]] defines the upper bound for the number of sequences - networks with more sequences will have the extra ones removed, keeping only the most abundant sequences/largest sequence frequency.
#' @param cell.limits list of integers or NULLs. cell.limits[[1]] the minimum threshold of cells which should produce a unique sequence (sequences below this threshold are removed from the network). If node.limits[[1]] is NULL, then there is no lower bound for the number of cells per unique sequence. node.limits[[2]] defines the upper bound for the number of cells per sequence -  sequence frequency.
#' @param weighted.edges boolean, T - edge weights will be equal to the string distance between a pair of nodes; F - edge weights = 1
#' @param weighted.germline boolean, T - adds weights to the edges connected to the germline, equal to the string distance between the germline and the specific connected nodes; F - edge weights = 1.
#' @param expand.intermediates boolean. T - will add inferred, intermediate nodes between nodes in the original network, determined by the string distance between a pair of nodes (for e.g., 2 nodes with an edge=string distance matrix of 3 will result in in 5 total nodes - 3 inferred nodes and 2 original ones, edges=1)
#' @param specific.networks #either an integer of max sorted clonotypes to be picked for network inference, 'all' for all clonotypes to be used, or list of specific clonotypes to create networks for.
#' @param network.level string determining the level at which networks should be built - 'intraclonal' will create intraclonal networks = networks for each sample and for each clonotype in the VDJ; 'global.clonotype' will create networks for each unique clonotype, irrespective of sample ids; 'global' will pool all clonotypes from all samples into a single global network;
#' 'forest.per.sample' and 'forest.global' are tree-specific methods, used when obtaining networks via network.algorithm='tree': 'forest.per.sample' will join the intraclonal trees in each sample, 'forest.global' will join ALL trees. Joining is determined by the forest.method parameter.
#' @param forest.method string determining how the trees should be joined if network.level='forest.per.sample' or 'forest.global'. 'single.germline' - trees will all be joined at a single germline (the one from the first clonotype), recalculating the string distances for the new adjacent nodes; 'multiple.germlines' - trees will all be joined in the same network, keeping the original germlines; 'multiple.germlines.joined' - same as 'multiple.germlines', but new edges will be added between the germlines.
#' @param random.seed numeric, seed for the random tie resolving method of resolve.ties.
#' @param parallel boolean with T - a parallelized mclapply will be used for each internal function, to accelerate computation; F - normal lapply will be used. Used best when having a large number of networks/clonotypes per sample.
#' @param as.igraph boolean - if T, the resulting networks will be igraph objects. Otherwise, they are converted to tidygraph tibble objects.

#' @return nested list of AntibodyForests objects for each sample and each clonotype. For example, output[[1]][[2]] denotes the AntibodyForests object of the first sample, second clonotype. If only a single clonotype and sample are available in the VDJ (or if the networks are joined via network.level = 'forest.global'), will output a single AntibodyForests object.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, sequence.type='VDJ.VJ.nt.trimmed', include.germline=T, network.algorithm='tree', resolve.ties=c('close.germline.distance', 'max.expansion'), node.features='OVA_binder', expand.intermediates=T, network.level='intraclonal')
#'}



AntibodyForests <- function(VDJ,
                            sequence.type,
                            include.germline,
                            network.algorithm,
                            directed,
                            distance.calculation,
                            resolve.ties,
                            connect.germline.to,
                            pruning.threshold,
                            remove.singletons,
                            keep.largest.cc,
                            VDJ.VJ.1chain,
                            node.features,
                            filter.na.features,
                            filter.specific.features,
                            node.limits,
                            cell.limits,
                            weighted.edges,
                            weighted.germline,
                            expand.intermediates,
                            specific.networks,
                            network.level,
                            forest.method,
                            random.seed,
                            parallel,
                            as.igraph){


  if(missing(VDJ)) stop('Please input your data as VDJ/df per clonotype list')
  if(missing(sequence.type)) sequence.type <- 'VDJ.VJ.nt.trimmed'
  if(missing(include.germline)) include.germline <- 'trimmed.ref'
  if(missing(network.algorithm)) network.algorithm <- 'tree'
  if(missing(directed) & network.algorithm != 'mst' & !stringr::str_detect(network.algorithm, 'prune')) directed <- T
  if(missing(directed) & (network.algorithm == 'mst' | stringr::str_detect(network.algorithm, 'prune'))) directed <- F
  if(missing(distance.calculation)) distance.calculation <- 'lv'
  if(missing(resolve.ties)) resolve.ties <- 'first'
  if(missing(connect.germline.to)) connect.germline.to <- 'min.adjacent'
  if(missing(pruning.threshold)) pruning.threshold <- 3
  if(missing(remove.singletons)) remove.singletons <- NULL
  if(missing(keep.largest.cc)) keep.largest.cc <- F
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
  if(missing(node.features)) node.features <- NULL
  if(missing(filter.na.features)) filter.na.features <- NULL
  if(missing(filter.specific.features)) filter.specific.features <- NULL
  if(missing(node.limits)) node.limits <- list(NULL, NULL)
  if(missing(cell.limits)) cell.limits <- list(NULL, NULL)
  if(missing(weighted.edges)) weighted.edges <- T
  if(missing(weighted.germline)) weighted.germline <- F
  if(missing(expand.intermediates)) expand.intermediates <- F
  if(missing(specific.networks)) specific.networks <- 'all'
  if(missing(network.level)) network.level <- 'intraclonal'
  if(missing(forest.method)) forest.method <- 'multiple.germlines' #join all at germline(single.germline), multiple.germlines, multiple.germlines.joined
  if(missing(random.seed)) random.seed <- 1
  if(missing(parallel)) parallel <- T
  if(missing(as.igraph)) as.igraph <- T


  #Global repertoire networks option:
  if(network.algorithm == 'global'){
    network.algorithm <- 'prune.distance'
    include.germline <- NULL
    directed <- F
    network.level <- 'global'
    keep.largest.cc <- T
  }

  features_to_select <- c('sample_id', 'clonotype_id')
  features_to_select <- unique(c(features_to_select, node.features))

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

  methods::setMethod(f='show', signature='AntibodyForests',
   definition=function(object) {
    cat('AntibodyForests object', '\n')

    cat(length(which(object@node_features$node_type == 'sequence')), ' sequence nodes across ', length(object@sample_id), ' sample(s) and ', length(object@clonotype_id), ' clonotype(s)', '\n')

    if(any(object@node_features$node_type == 'intermediate')){
      cat(length(which(object@node_features$node_type == 'intermediate')), ' intermediate nodes', '\n')
    }

    if(any(object@node_features$node_type == 'inferred')){
      cat(length(which(object@node_features$node_type == 'inferred')), ' inferred sequence nodes', '\n')
    }

    if(any(object@node_features$node_type == 'bulk')){
      cat(length(which(object@node_features$node_type == 'bulk')), ' bulk sequence nodes', '\n')
    }

    cat(length(which(object@node_features$node_type == 'sequence')) + length(which(object@node_features$node_type == 'intermediate')) + length(which(object@node_features$node_type == 'inferred')) + length(which(object@node_features$node_type == 'bulk')),
        ' total network nodes', '\n'
    )


    cat('Sample id(s): ', paste0(object@sample_id, collapse = ', '), '\n')

    cat('Clonotype id(s): ', paste0(object@clonotype_id, collapse = ', '), '\n')

    networks <- c(object@network_algorithm)

    if(!is.null(object@plot_ready)){
      networks <- c(networks, 'plot_ready')
    }

    if(!is.null(object@phylo)){
      networks <- c(networks, 'phylo')
    }

    if(!is.null(object@heterogeneous)){
      networks <- c(networks, 'heterogeneous')

    }

    if(!is.null(object@reactivity)){
      networks <- c(networks, 'reactivity')

    }

    if(!is.null(object@dynamic)){
      networks <- c(networks, 'dynamic')
    }


    if(!is.null(object@metrics)){
      networks <- c(networks, 'metrics')

    }

    if(!is.null(object@node_transitions)){
      networks <- c(networks, 'node_transitions')

    }

    if(!is.null(object@paths)){
      networks <- c(networks, 'paths')
    }

    cat('Networks/analyses available: ', paste0(networks, collapse = ', '))

    cat('\n')

    cat('\n')

    cat('\n')

   }
  )
  get_sequence_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';'){
    if(split.x==T) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
    if(split.y==T) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

    ccombs <- expand.grid(x,y)
    ccombs<-paste0(ccombs[,1], ccombs[,2])
    ccombs <- paste0(ccombs, collapse=collapse.by)

    return(ccombs)
  }

  extract_MIXCR <- function(VDJ.matrix, chain.to.extract, as.nucleotide, regions.to.extract){
    if(missing(VDJ.matrix)) stop('Input the VDJ dataframe obtained after calling VDJ_call_MIXCR')
    if(missing(chain.to.extract)) chain.to.extract <- 'VDJ.VJ'
    if(missing(as.nucleotide)) as.nucleotide <- T
    if(missing(regions.to.extract)) regions.to.extract <- list('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')


    VDJ.matrix$mixcr_assembled <- rep('', nrow(VDJ.matrix))

    if(!('VDJ_SHM' %in% colnames(VDJ.matrix))) stop('Please use the output of the VDJ_call_MIXCR function')

    if(chain.to.extract!='VDJ.VJ'){
      if(as.nucleotide){
        col_name <- paste0(chain.to.extract, '_nSeq')
      }else{
        col_name <- paste0(chain.to.extract, '_aaSeq')
      }

      for(region in regions.to.extract){
        VDJ.matrix$mixcr_assembled <- paste0(VDJ.matrix$mixcr_assembled, gsub('_', '', VDJ.matrix[, paste0(col_name, region)]))
      }
    }else if(chain.to.extract=='VDJ.VJ'){
      if(as.nucleotide==T){
        col_name <- '_nSeq'
      }else{
        col_name <- '_aaSeq'
      }
      extracted_VDJ <- rep('', nrow(VDJ.matrix))
      extracted_VJ <- rep('', nrow(VDJ.matrix))

      for(region in regions.to.extract){
        extracted_VDJ <- paste0(extracted_VDJ, gsub('_', '', VDJ.matrix[,paste0('VDJ',col_name, region)]))
        extracted_VJ <- paste0(extracted_VJ, gsub('_', '', VDJ.matrix[,paste0('VJ',col_name, region)]))
      }

      VDJ.matrix$mixcr_assembled <- paste0(extracted_VDJ, extracted_VJ)

    }
    return(VDJ.matrix)
  }

  transform_clonotype_to_network_df <- function(clonotype_df){
    if(VDJ.VJ.1chain){
      clonotype_df <- clonotype_df[which(clonotype_df$Nr_of_VDJ_chains==1 & clonotype_df$Nr_of_VJ_chains==1),]
      if(nrow(clonotype_df)==0){
        return(NULL)
      }
    }

   if(sequence.type=='cdr3.aa'){
     combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_cdr3s_aa, clonotype_df$VJ_cdr3s_aa)
     clonotype_df$network_sequences <- combined_sequences

   }else if(sequence.type=='cdr3.nt'){
     combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_cdr3s_nt, clonotype_df$VJ_cdr3s_nt)
     clonotype_df$network_sequences <- combined_sequences

   }else if(sequence.type=='VDJ.VJ.nt.trimmed'){
     if(!('VDJ_trimmed_ref' %in% colnames(clonotype_df))){
       stop('Please use trim.and.align=T when creating your VGM object')
     }
     combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_sequence_nt_trimmed, clonotype_df$VJ_sequence_nt_trimmed)
     clonotype_df$network_sequences <- combined_sequences

   }else if(sequence.type=='VDJ.VJ.nt.raw'){
     combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_sequence_nt_raw, clonotype_df$VJ_sequence_nt_raw)
     clonotype_df$network_sequences <- combined_sequences

   }else if(sequence.type=='VDJ.VJ.nt.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VDJ.VJ')
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VDJ.nt.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VDJ')
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VJ.nt.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VJ')
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VDJ.VJ.aa.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VDJ.VJ', as.nucleotide = F)
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VDJ.aa.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VDJ', as.nucleotide = F)
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VJ.aa.mixcr'){
     clonotype_df <- extract_MIXCR(clonotype_df, chain.to.extract = 'VJ', as.nucleotide = F)
     clonotype_df$network_sequences <- clonotype_df$mixcr_assembled

   }else if(sequence.type=='VDJ.nt.trimmed'){
     network_sequences <- clonotype_df$VDJ_sequence_nt_trimmed
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VDJ.nt.raw'){
     network_sequences <- clonotype_df$VDJ_sequence_nt_raw
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VDJ.cdr3s.aa'){
     network_sequences <- clonotype_df$VDJ_cdr3s_aa
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VDJ.cdr3s.nt'){
     network_sequences <- clonotype_df$VDJ_cdr3s_nt
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VJ.nt.trimmed'){
     network_sequences <- clonotype_df$VJ_sequence_nt_trimmed
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VJ.nt.raw'){
     network_sequences <- clonotype_df$VJ_sequence_nt_raw
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VJ.cdr3s.aa'){
     network_sequences <- clonotype_df$VJ_cdr3s_aa
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VJ.cdr3s.nt'){
     network_sequences <- clonotype_df$VJ_cdr3s_nt
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VDJ.aa'){
     network_sequences <- clonotype_df$VDJ_sequence_aa
     clonotype_df$network_sequences <- network_sequences

   }else if(sequence.type=='VJ.aa'){
     network_sequences <- clonotype_df$VJ_sequence_aa
     clonotype_df$network_sequences <- network_sequences

   }else{
     stop('Sequence type unavailable - please check the function documentation for the supported sequence type; or check your spelling :)')
   }


  all_sequences <- unlist(lapply(clonotype_df$network_sequences, function(x) stringr::str_split(x, ';')))

  unique_sequences <- unlist(unique(all_sequences))
  #cell_number <- unlist(lapply(unique_sequences, function(x) length(all_sequences[all_sequences==x])))
  cell_barcodes <- lapply(unique_sequences, function(x) clonotype_df$barcode[which(stringr::str_detect(clonotype_df$network_sequences, x))])
  cell_number <- unlist(lapply(cell_barcodes, function(x) length(x)))


  network_df <- data.frame(network_sequences = unique_sequences,
                           cell_number = cell_number,
                           cell_barcodes = matrix(cell_barcodes))

  minimum.cells <- cell.limits[[1]]
  maximum.cells <- cell.limits[[2]]

  if(!is.null(minimum.cells)){
    network_df <- network_df[network_df$cell_number >= minimum.cells,]
  }

  if(!is.null(maximum.cells)){
    network_df <- network_df[network_df$cell_number <= minimum.cells,]
  }

  if(nrow(network_df) == 0){
    return(NULL)
  }

  minimum.sequences <- node.limits[[1]]
  maximum.sequences <- node.limits[[2]]

  if(!is.null(minimum.sequences)){
    if(nrow(network_df)<minimum.sequences) {
      return(NULL)
    }
  }

  if(!is.null(maximum.sequences)){
    if(nrow(network_df)>maximum.sequences) network_df <- network_df[1:maximum.sequences,]
  }

  network_df <- network_df[order(network_df$cell_number, decreasing = T),]


  for(feature in features_to_select){
    feature_list <- list()
    feature_counts <- list()
    #clonotype_df[[feature]][is.na(clonotype_df[[feature]]) | is.null(clonotype_df[[feature]]) | clonotype_df[[feature]] == ''] <- 'unknown'

    for(i in 1:length(network_df$network_sequences)){
      seq <- network_df$network_sequences[i]
      #counts <- lapply(network_df$cell_barcodes[[i]], function(x) clonotype_df[feature][which(clonotype_df$barcode==x),]) %>%
      #          unlist() %>%
      #          table() %>%
      #          c()

      counts <- clonotype_df[feature][which(clonotype_df$network_sequences==seq),] %>%
                unlist() %>%
                table() %>%
                c()

      feats <- names(counts)

      feats <- feats[counts != 0]
      counts <- counts[counts != 0]

      if(length(feats) == 0){
        feats <- 'unknown'
      }

      if(length(counts) == 0){
        counts <- 1
      }

      feats[is.na(feats) | is.null(feats) | feats == ''] <- 'unknown'


      feature_list[[i]] <- feats
      feature_counts[[i]] <- counts
    }

    network_df[[feature]] <- feature_list

    network_df[[paste0(feature, '_counts')]] <- feature_counts
  }

  if(!is.null(filter.na.features)){
    for(col in filter.na.features){
      if(all(is.na(unlist(network_df[,col])))){
        return(NULL)
      }
    }
  }

  if(!is.null(filter.specific.features)){
    if( !(filter.specific.features[[2]] %in% unlist(network_df[,filter.specific.features[[1]]]))  ) {
      return(NULL)
    }
  }

  #Also highlight the most_expanded sequences before joining (otherwise harder to distinguish per unique clonotype)
  network_df$most_expanded <- 'no'
  network_df$most_expanded[1] <- 'yes' #As we prev ordered the sequence frequencies, the most_expandeds are always the first row in a our networks dfs
  #To decide if most expanded should be 1 (and pick the first seq) or multiple (and pick via==max())

  #Add this before as it will be used to determine where the sequence and intermediate nodes are in the network dataframe
  network_df$germline <- 'no'
  network_df$sequence_id <- 1:nrow(network_df)

  if(!is.null(include.germline)){
    if(any(include.germline != 'trimmed.ref') & any(include.germline != 'cdr3.nt')){
      if(length(include.germline) == 2){
        VDJ <- clonotype_df[[include.germline[1]]]
        VJ <- clonotype_df[[include.germline[2]]]
        if(VDJ=='' | is.na(VDJ) | is.null(VDJ)){
          return(NULL)
        }

        if(VJ=='' | is.na(VJ) | is.null(VJ)){
          return(NULL)
        }
        germline_seq <- paste0(VDJ, VJ)
      }else{
        germline_seq <- clonotype_df[[include.germline]]
      }
      unique_germline_seq <- unlist(unique(germline_seq))
      germline_seq_frequencies <- unlist(lapply(unique_germline_seq, function(x) length(which(germline_seq==x))))
      germline_seq <- unique_germline_seq[which.max(germline_seq_frequencies)]

      if(germline_seq=='' | is.na(germline_seq) | is.null(germline_seq)){
        return(NULL)
      }

    }else if(include.germline == 'trimmed.ref'){
      if(!('VDJ_trimmed_ref' %in% colnames(clonotype_df))){
        stop('Please use trim.and.align=T when creating your VGM object')
      }

      if(stringr::str_detect(sequence.type, 'VDJ.VJ')){
        VDJ <- clonotype_df$VDJ_trimmed_ref
        VJ <- clonotype_df$VJ_trimmed_ref
        if(VDJ=='' | is.na(VDJ) | is.null(VDJ)){
          return(NULL)
        }

        if(VJ=='' | is.na(VJ) | is.null(VJ)){
          return(NULL)
        }
        germline_seq <- paste0(VDJ,VJ)
      }else if(stringr::str_detect(sequence.type, 'VDJ')){
        germline_seq <- clonotype_df$VDJ_trimmed_ref
      }else{
        germline_seq <- clonotype_df$VJ_trimmed_ref
      }

      #germline_seq <- gsub('-', '', germline_seq)  #see how this affects the networks strcuture (previously was not included)
      unique_germline_seq <- unlist(unique(germline_seq))
      germline_seq_frequencies <- unlist(lapply(unique_germline_seq, function(x) length(which(germline_seq==x))))
      germline_seq <- unique_germline_seq[which.max(germline_seq_frequencies)]

      if(germline_seq=='' | is.na(germline_seq) | is.null(germline_seq)){
        return(NULL)
      }

    }else if(include.germline == 'cdr3.nt'){
      VDJ <- clonotype_df$VDJ_cdr3s_nt
      VJ <- clonotype_df$VJ_cdr3s_nt
      germline_seq <- paste0(VDJ,VJ)
      #germline_seq <- gsub('-', '', germline_seq)  #see how this affects the networks strcuture (previously was not included)
      unique_germline_seq <- unlist(unique(germline_seq))
      germline_seq_frequencies <- unlist(lapply(unique_germline_seq, function(x) length(which(germline_seq==x))))
      germline_seq <- unique_germline_seq[which.max(germline_seq_frequencies)]

      if(germline_seq=='' | is.na(germline_seq) | is.null(germline_seq)){
        return(NULL)
      }
    }

    network_df <- rbind(network_df, rep(NA, ncol(network_df)))
    network_df$germline[nrow(network_df)] <- 'yes'
    network_df$network_sequences[nrow(network_df)] <- germline_seq
    network_df$clonotype_id[nrow(network_df)] <- unlist(unique(network_df$clonotype_id))[1]
    network_df$sample_id[nrow(network_df)] <- unlist(unique(network_df$sample_id))[1]
    network_df$cell_number[nrow(network_df)] <- NA
    network_df$cell_barcodes[nrow(network_df)] <- 'germline'
    network_df$most_expanded[nrow(network_df)] <- 'germline'

    #Add distance from germline before joining networks...more efficient downstream integration into AntibodyForests_metrics (otherwise would have to get unique germlines per clonotypes using the clonotype_ids...but this way it's easier)
    network_df$distance_from_germline <- stringdist::stringdistmatrix(network_df$network_sequences, germline_seq)
  }

  #Add sequence ids/labels
  network_df$sequence_id <- 1:nrow(network_df)

  return(network_df)
 }


 calculate_adjacency_matrix_tree <- function(network_df){

   sequences <- network_df$network_sequences
   distance_matrix <- stringdist::stringdistmatrix(sequences, sequences, method=distance.calculation)
   diag(distance_matrix) <- Inf
   final_adjacency_matrix <- matrix(0, nrow(distance_matrix), ncol(distance_matrix))
   diag(final_adjacency_matrix) <- Inf


   nodes_in_network <- NULL
   nodes_not_in_network <- 1:ncol(distance_matrix)
   all_nodes <- 1:ncol(distance_matrix)


   if(!is.null(include.germline)){
     germline_node <- ncol(distance_matrix)

     if(connect.germline.to=='min.adjacent'){
       adjacent_nodes <- which(distance_matrix[germline_node,]==min(distance_matrix[germline_node,]))

     }else if(connect.germline.to=='threshold.adjacent' & !is.null(pruning.threshold[1])){
       adjacent_nodes <- which(distance_matrix[germline_node,]<=pruning.threshold[1])
     }

     for(node in adjacent_nodes){
       if(!weighted.germline) { edge_weight <- 1 }
       else { edge_weight <-  distance_matrix[node,germline_node] }

       final_adjacency_matrix[germline_node, node] <- edge_weight

       if(!directed){
         final_adjacency_matrix[node, germline_node] <- edge_weight
       }
     }

     nodes_not_in_network <- 1:(ncol(distance_matrix)-1)
     all_nodes <- 1:(ncol(distance_matrix)-1)
     distance_matrix_w_germline <- distance_matrix
     distance_matrix[germline_node,] <- Inf
     distance_matrix[,germline_node] <- Inf
     nodes_in_network <- unique(c(nodes_in_network, adjacent_nodes))
     nodes_not_in_network <- nodes_not_in_network[-nodes_in_network]
   }

   while(length(nodes_not_in_network)>0){
     distance_matrix_copy <-  distance_matrix

     if(length(nodes_not_in_network)!=length(all_nodes)){
       distance_matrix_copy[nodes_not_in_network,] <- Inf
       distance_matrix_copy[,nodes_in_network] <- Inf
     }


     current_edge_weight <- min(distance_matrix_copy)
     current_potential_node_pairs <- which(distance_matrix_copy==min(distance_matrix_copy), arr.ind=T)

     if(weighted.edges){
       edge_weight <- current_edge_weight
     }else{
       edge_weight <- 1
     }

     if(nrow(current_potential_node_pairs)!=1){
       current_potential_node_pairs <- c(as.data.frame(t(current_potential_node_pairs)))
       if(is.null(resolve.ties)){
         just_added <- c()
         for(node_pair in current_potential_node_pairs){
           in_network <- node_pair[which(node_pair %in% nodes_in_network)]
           not_in_network <- node_pair[which(!(node_pair %in% nodes_in_network))]

           final_adjacency_matrix[in_network, not_in_network] <- edge_weight

           if(!directed){
             final_adjacency_matrix[not_in_network, in_network] <- edge_weight
           }
           just_added <- unique(c(just_added, in_network, not_in_network))
         }

         nodes_in_network <- unique(c(nodes_in_network, just_added))
         nodes_not_in_network <- all_nodes[-nodes_in_network]
         next

       }else{
        for(tie_algorithm in resolve.ties){
          if(length(current_potential_node_pairs)==1){
            break
          }

          if(tie_algorithm=='first'){

            current_potential_node_pairs <- current_potential_node_pairs[1]

          }else if(tie_algorithm=='random'){

            sampled_index <- sample(1:length(current_potential_node_pairs), 1)
            current_potential_node_pairs <- current_potential_node_pairs[sampled_index]

          }else if(tie_algorithm=='close.germline.distance' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            avg_distance_per_pair <- unlist(lapply(current_potential_node_pairs, function(x) (distance_matrix_w_germline[x[1], germline_node] + distance_matrix_w_germline[x[2], germline_node])/2 ))
            min_index <- which(avg_distance_per_pair==min(avg_distance_per_pair))
            current_potential_node_pairs <- current_potential_node_pairs[min_index]


          }else if(tie_algorithm=='close.germline.edges' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            current_graph <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=T, diag=F)
            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(x%in%nodes_in_network)])
            unique_nodes <- unname(unlist(unique_nodes))

            distances_from_germline <- igraph::distances(current_graph, v=unique_nodes, to=germline_node, weights=NA, algorithm='unweighted')

            min_index <- which(distances_from_germline==min(distances_from_germline))
            current_potential_node_pairs <- current_potential_node_pairs[min_index]

          }else if(tie_algorithm=='close.germline.weighted' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            current_graph <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=T, diag=F)
            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(x%in%nodes_in_network)])
            unique_nodes <- unname(unlist(unique_nodes))

            distances_from_germline <- igraph::distances(current_graph, v=unique_nodes, to=germline_node, weights=NA, algorithm='dijkstra')

            min_index <- which(distances_from_germline==min(distances_from_germline))
            current_potential_node_pairs <- current_potential_node_pairs[min_index]


          }else if(tie_algorithm=='far.germline.distance' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            avg_distance_per_pair <- unlist(lapply(current_potential_node_pairs, function(x) (distance_matrix_w_germline[x[1], germline_node] + distance_matrix_w_germline[x[2], germline_node])/2 ))
            max_index <- which(avg_distance_per_pair==max(avg_distance_per_pair))
            current_potential_node_pairs <- current_potential_node_pairs[max_index]


          }else if(tie_algorithm=='far.germline.edges' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            current_graph <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=T, diag=F)
            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(x%in%nodes_in_network)])
            unique_nodes <- unname(unlist(unique_nodes))

            distances_from_germline <- igraph::distances(current_graph, v=unique_nodes, to=germline_node, weights=NA, algorithm='unweighted')

            max_index <- which(distances_from_germline==max(distances_from_germline))
            current_potential_node_pairs <- current_potential_node_pairs[max_index]

          }else if (tie_algorithm=='far.germline.weighted' & !is.null(include.germline)){

            germline_node <- ncol(distance_matrix)
            current_graph <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=T, diag=F)
            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(x%in%nodes_in_network)])
            unique_nodes <- unname(unlist(unique_nodes))

            distances_from_germline <- igraph::distances(current_graph, v=unique_nodes, to=germline_node, weights=NA, algorithm='dijkstra')

            max_index <- which(distances_from_germline==max(distances_from_germline))
            current_potential_node_pairs <- current_potential_node_pairs[max_index]

          }else if(tie_algorithm=='max.degree'){

            if(is.null(pruning.threshold)){
              stop('Please input a valid pruning distance threshold to resolve ties; otherwise, all nodes will have the same degree')
            }

            pruned_distance_matrix <- distance_matrix_copy
            pruned_distance_matrix[pruned_distance_matrix!=Inf & pruned_distance_matrix>pruning.threshold] <- 0

            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(!(x%in%nodes_in_network))])
            unique_nodes <- unname(unlist(unique_nodes))

            degree_per_node <- unlist(lapply(unique_nodes, function(x) length(which(pruned_distance_matrix[,x]!= Inf & pruned_distance_matrix[,x]!=0))))
            max_index <- which(degree_per_node==max(degree_per_node))
            current_potential_node_pairs <- current_potential_node_pairs[max_index]


          }else if(tie_algorithm=='min.degree'){

            pruned_distance_matrix <- distance_matrix_copy
            pruned_distance_matrix[pruned_distance_matrix!=Inf & pruned_distance_matrix>pruning.threshold] <- 0

            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(!(x%in%nodes_in_network))])
            unique_nodes <- unname(unlist(unique_nodes))

            degree_per_node <- unlist(lapply(unique_nodes, function(x) length(which(pruned_distance_matrix[,x]!= Inf & pruned_distance_matrix[,x]!=0))))

            min_index <- which(degree_per_node==min(degree_per_node))

            current_potential_node_pairs <- current_potential_node_pairs[min_index]


          }else if(tie_algorithm=='max.expansion'){

            expansion_list <- network_df$cell_number

            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(!(x%in%nodes_in_network))])
            unique_nodes <- unname(unlist(unique_nodes))

            expansion_list <- expansion_list[unique_nodes]

            max_index <- which(expansion_list==max(expansion_list))

            current_potential_node_pairs <- current_potential_node_pairs[max_index]


          }else if(tie_algorithm=='min.expansion'){

            expansion_list <- network_df$cell_number

            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(!(x%in%nodes_in_network))])
            unique_nodes <- unname(unlist(unique_nodes))

            expansion_list <- expansion_list[unique_nodes]

            min_index <- which(expansion_list==min(expansion_list))

            current_potential_node_pairs <- current_potential_node_pairs[min_index]


          }else if(stringr::str_detect(tie_algorithm, '-')){
            tie_algorithm <- unlist(stringr::str_split(tie_algorithm, '-'))
            if(!(tie_algorithm[2] %in% colnames(network_df))){
              stop('Please input a valid column name for the custom tie-breaking algorithm')
            }

            tie_values <- network_df[,tie_algorithm[2]]
            unique_nodes <- lapply(current_potential_node_pairs, function(x) x[which(!(x%in%nodes_in_network))])
            unique_nodes <- unname(unlist(unique_nodes))
            tie_values <- tie_values[unique_nodes]

            if(tie_algorithm[1]=='max' & is.numeric(tie_values)){
              max_index <- which(tie_values==max(tie_values))
              current_potential_node_pairs <- current_potential_node_pairs[max_index]

            }else if(tie_algorithm[1]=='min' & is.numeric(tie_values)){
              min_index <- which(tie_values==min(tie_values))
              current_potential_node_pairs <- current_potential_node_pairs[min_index]

            }else{

              index <- which(tie_values==tie_algorithm[1])

              if(length(index)==0){
                index <- 1:length(current_potential_node_pairs)
              }

              current_potential_node_pairs <- current_potential_node_pairs[index]

            }
          }else{
            stop('Not implemented yet')
          }
         }
        }
      }

       node_pair <- unname(unlist(current_potential_node_pairs))

       in_network <- node_pair[which(node_pair %in% nodes_in_network)]
       not_in_network <- node_pair[which(!(node_pair %in% nodes_in_network))]

       if(length(in_network)==0){
         in_network <- node_pair[1]
         not_in_network <- node_pair[2]
       }

       final_adjacency_matrix[in_network, not_in_network] <- edge_weight

       if(!directed){
         final_adjacency_matrix[not_in_network, in_network] <- edge_weight
       }

       nodes_in_network <- unique(c(nodes_in_network, in_network, not_in_network))
       nodes_not_in_network <- all_nodes[-nodes_in_network]
    }

   diag(final_adjacency_matrix) <- 0
   distance_matrix <- NULL
   distance_matrix_copy <- NULL

   return(final_adjacency_matrix)
 }


 calculate_adjacency_matrix_prune <- function(network_df){

    sequences <- network_df$network_sequences
    distance_matrix <- stringdist::stringdistmatrix(sequences, sequences, method=distance.calculation)
    diag(distance_matrix) <- Inf
    final_adjacency_matrix <- distance_matrix
    adjacency_matrix_copy <- final_adjacency_matrix

    pruning.methods <- unlist(stringr::str_split(network.algorithm, '\\.'))[-1]

    for(i in 1:length(pruning.methods)){
      if(pruning.methods[i]=='distance'){
        final_adjacency_matrix[final_adjacency_matrix>pruning.threshold[i]] <- 0
      }else if(pruning.methods[i]=='degree'){
        adjacency_matrix_copy[adjacency_matrix_copy>pruning.threshold[i-1]] <- 0

        for(j in 1:ncol(adjacency_matrix_copy)){
          if(length(adjacency_matrix[j,][adjacency_matrix[j,]!=0])<pruning.threshold[i]){
            final_adjacency_matrix[j,] <- 0
            final_adjacency_matrix[,j] <- 0
          }
        }

        adjacency_matrix_copy <- NULL

      }else if(pruning.methods[i]=='expansion'){
        sequence_counts <- network_df$cell_number
        pruned_nodes <- which(sequence_counts < pruning.threshold[3])
        for(node in pruned_nodes){
          final_adjacency_matrix[node,] <- 0
          final_adjacency_matrix[,node] <- 0
        }
      }
    }

    if(!weighted.edges){
      final_adjacency_matrix[final_adjacency_matrix!=0] <- 1
    }

    #Post-processing - remove singleton nodes
    if(!is.null(remove.singletons)){
      g <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=TRUE, diag=F)
      connected_components <- igraph::components(g)
      components_to_remove <- which(connected_components$csize <= remove.singletons)
      vertices_to_remove <- which(connected_components$membership %in% components_to_remove)
      g <- igraph::delete_vertices(g, vertices_to_remove)
      final_adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(g, attr = 'weight'))
    }

    if(keep.largest.cc){
      g <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix, mode='undirected', weighted=TRUE, diag=F)
      connected_components <- igraph::components(g)
      max_component <- max(connected_components$csize)
      components_to_remove <- which(connected_components$csize < max_component)
      vertices_to_remove <- which(connected_components$membership %in% components_to_remove)
      g <- igraph::delete_vertices(g, vertices_to_remove)
      final_adjacency_matrix <- as.matrix(igraph::as_adjacency_matrix(g, attr = 'weight'))
    }


    if(!is.null(include.germline)){
      germline_node <- which(network_df$germline == 'yes')
      if(sum(final_adjacency_matrix[germline_node,]) == 0){
        graph_no_germline <- igraph::graph_from_adjacency_matrix(final_adjacency_matrix[-germline_node, -germline_node], mode='undirected', weighted=T, diag=F)
        connected_components <- igraph::components(graph_no_germline)

        if(connect.germline.to=='min.adjacent'){
          adjacent_nodes <- which(distance_matrix[germline_node,]==min(distance_matrix[germline_node,]))

          for(node in adjacent_nodes){
            if(!weighted.germline) { edge_weight <- 1 }
            else { edge_weight <-  distance_matrix[node,germline_node] }

            final_adjacency_matrix[node, germline_node] <- edge_weight
            final_adjacency_matrix[germline_node, node] <- edge_weight
          }

        }else if(connect.germline.to!='min.adjacent'){
          if(connect.germline.to=='largest.connected.component'){
            largest_components <- which(connected_components$csize==max(connected_components$csize))

            for(component in largest_components){
              nodes_in_component <- which(connected_components$membership==component)
              adjacent_nodes <- which(distance_matrix[germline_node,nodes_in_component]==min(distance_matrix[germline_node,nodes_in_component]))

              for(node in adjacent_nodes){
                if(weighted.germline){
                  edge_weight <- distance_matrix[node, germline_node]
                }else{
                  edge_weight <- 1
                }

                final_adjacency_matrix[node,germline_node] <- edge_weight
                final_adjacency_matrix[germline_node,node] <- edge_weight
              }
            }
          }else if(connect.germline.to=='all.connected.components'){
            all_components <- unique(connected_components$membership)

            for(component in all_components){
              nodes_in_component <- which(connected_components$membership==component)
              adjacent_nodes <- which(distance_matrix[germline_node,nodes_in_component]==min(distance_matrix[germline_node,nodes_in_component]))

              for(node in adjacent_nodes){
                if(weighted.germline){
                  edge_weight <- distance_matrix[node, germline_node]
                }else{
                  edge_weight <- 1
                }

                final_adjacency_matrix[node,germline_node] <- edge_weight
                final_adjacency_matrix[germline_node,node] <- edge_weight
              }
            }
          }else if(connect.germline.to=='all.connected.components.non.single'){
            all_components <- which(connected_components$csize!=1)
            for(component in all_components){
              nodes_in_component <- which(connected_components$membership==component)
              adjacent_nodes <- which(distance_matrix[germline_node,nodes_in_component]==min(distance_matrix[germline_node,nodes_in_component]))

              for(node in adjacent_nodes){
                if(weighted.germline){
                  edge_weight <- distance_matrix[node, germline_node]
                }else{
                  edge_weight <- 1
                }

                final_adjacency_matrix[node,germline_node] <- edge_weight
                final_adjacency_matrix[germline_node,node] <- edge_weight
              }
            }
          }else if(connect.germline.to=='none'){
            final_adjacency_matrix <- final_adjacency_matrix
          }
        }
      }
    }

    return(final_adjacency_matrix)
 }


 create_phylo_trees <- function(network_df){
   requireNamespace('phangorn')
   requireNamespace('seqinr')

   if(nrow(network_df) < 3){
     stop('Ensure you have at least 3 sequences (incl. germline) to create phylogenetic trees - use the node.limits parameter')
   }

   phylogenetic_method <- unlist(stringr::str_split(network.algorithm, '\\.'))
   sequences <- network_df$network_sequences
   distance_matrix <- stringdist::stringdistmatrix(sequences, sequences, method=distance.calculation)
   phylo_tree <- ape::nj(distance_matrix)

   #Roots at the germline node
   if(!is.null(include.germline)){
     phylo_tree <- phytools::reroot(phylo_tree, node.number=which(network_df$germline=='yes'))
   }

   #Algorithms for mp and ml tree reconstruction - must include the mafft directory to perform multiple sequence alignment
   if(phylogenetic_method[3]!='nj'){

     if(Sys.which('mafft')!=''){

       dir <- './tempdir'
       if(!dir.exists(dir)) dir.create(dir)

       #Unique network id to avoid clashes during parallel computation
       network_id <- paste0(unique(network_df$sample_id), '-', unique(network_df$clonotype_id))
       fasta_file <- paste0(dir, '/', network_id, '_', 'tempseq.fasta')
       fasta_aligned <- paste0(dir, '/', network_id, '_', 'tempseq_aligned.fasta')

       seqinr::write.fasta(as.list(sequences), names=1:nrow(network_df), file.out=fasta_file)

       system(paste0('mafft ', fasta_file, ' > ', fasta_aligned))

       phylo_data <- phangorn::read.phyDat(fasta_aligned, format='fasta')

       system(paste0('rm -r ', fasta_file))
       system(paste0('rm -r ', fasta_aligned))


     }else{
       alignment <- network_df$network_sequences %>%
                    matrix() %>%
                    ape::as.dnabin() %>%
                    ape::clustal()

       phylo_data <- phangorn::read.phyDat(alignment, format='fasta')

     }

     if(phylogenetic_method[3]=='mp'){
       phylo_tree <- phangorn::optim.parsimony(phylo_tree, phylo_data)

     }else if(phylogenetic_method[3]=='ml'){
       fit <- phangorn::pml(phylo_tree, data=phylo_data)
       fit <- phangorn::optim.pml(fit, TRUE)
       phylo_tree <- fit$tree

     }
   }

   germline <- which(network_df$germline == 'yes')

   if(length(germline) > 0){
     phylo_tree <- phytools::reroot(phylo_tree, node.number = germline)
     phylo_tree <- ape::reorder.phylo(phylo_tree, 'postorder')

     edges <- phylo_tree$edge
     mrca_to_germline <- which(edges[, 2] == germline)
     mrca <- edges[mrca_to_germline, 1]

     mrca_to_interm <- which(edges[, 1] == mrca & edges[, 2] != germline)

     phylo_tree$edge <- phylo_tree$edge[-mrca_to_germline,]
     phylo_tree$edge.length <- phylo_tree$edge.length[-mrca_to_germline]
     phylo_tree$edge[mrca_to_interm][1] <- germline

     phylo_tree$edge[,1][phylo_tree$edge[,1] > mrca] <- phylo_tree$edge[,1][phylo_tree$edge[,1] > mrca] - 1
     phylo_tree$edge[,2][phylo_tree$edge[,2] > mrca] <- phylo_tree$edge[,2][phylo_tree$edge[,2] > mrca] - 1
     phylo_tree$Nnode <- phylo_tree$Nnode - 1
   }

   g <- igraph::graph_from_edgelist(phylo_tree$edge, directed = directed)

   igraph::E(g)$weight <- phylo_tree$edge.length

   clonotype_id <- unique(network_df$clonotype_id)
   sample_id <- unique(network_df$sample_id)

   for(i in 1:phylo_tree$Nnode){
     network_df <- rbind(network_df, NA)
     network_df$sequence_id[nrow(network_df)] <- 'intermediate'
     network_df$clonotype_id[nrow(network_df)] <- clonotype_id
     network_df$sample_id[nrow(network_df)] <- sample_id
   }

   adjacency_matrix <- igraph::as_adjacency_matrix(g, sparse = F, attr = 'weight')

   return(list(adjacency_matrix=adjacency_matrix, network_df=network_df))
 }

 create_mst_trees <- function(network_df){

   sequences <- network_df$network_sequences
   distance_matrix <- stringdist::stringdistmatrix(sequences, sequences, method=distance.calculation) #biggest bottleneck currently (for v large graphs)
   g <- igraph::graph_from_adjacency_matrix(distance_matrix, mode='undirected', weighted=T, diag=F)
   g <- igraph::mst(g)

   adjacency_matrix <- igraph::as_adjacency_matrix(g, attr = 'weight')
   return(list(adjacency_matrix=adjacency_matrix, network_df=network_df))

 }

 expand_intermediates <- function(adjacency.network.list){
    network_df <- adjacency.network.list$network_df
    adjacency_matrix <- adjacency.network.list$adjacency_matrix

    if(!directed){
      node_number <- sum(adjacency_matrix)/2 + nrow(adjacency_matrix) - length(which(adjacency_matrix!=0))/2
    }else{
      node_number <- sum(adjacency_matrix) + nrow(adjacency_matrix) - length(which(adjacency_matrix!=0))
    }

    added_intermediates <- node_number - nrow(adjacency_matrix)

    final_adjacency_matrix <- matrix(0, node_number, node_number)
    final_adjacency_matrix[1:nrow(adjacency_matrix), 1:nrow(adjacency_matrix)] <- adjacency_matrix

    new_node <- nrow(adjacency_matrix) + 1
    for(i in 1:nrow(adjacency_matrix)){
      neighbours_to_expand <- which(adjacency_matrix[i, ]>1)
      if(length(neighbours_to_expand)==0){
        next
      }

      for(neighbour in neighbours_to_expand){
        edge_weight <- adjacency_matrix[i, neighbour]
        prev_added_node <- i

        while(edge_weight>1){
          final_adjacency_matrix[prev_added_node, new_node] <- 1
          if(!directed){
            final_adjacency_matrix[new_node, prev_added_node] <- 1
          }
          prev_added_node <- new_node
          new_node <- new_node + 1
          edge_weight <- edge_weight - 1
        }

       final_adjacency_matrix[i, neighbour] <- 0
       final_adjacency_matrix[neighbour, i] <- 0
       adjacency_matrix[i, neighbour] <- 0
       adjacency_matrix[neighbour, i] <- 0
       final_adjacency_matrix[new_node-1, neighbour] <- 1

       if(!directed){
         final_adjacency_matrix[neighbour, new_node-1] <- 1
       }
      }
    }

    clonotype_id <- unique(network_df$clonotype_id)
    sample_id <- unique(network_df$sample_id)

    #Add the intermediate rows into the network_df
    for(i in 1:added_intermediates){
      network_df <- rbind(network_df, NA)
      network_df$sequence_id[nrow(network_df)] <- 'intermediate'
      network_df$clonotype_id[nrow(network_df)] <- clonotype_id
      network_df$sample_id[nrow(network_df)] <- sample_id
    }


    return(list(adjacency_matrix=final_adjacency_matrix, network_df=network_df))
 }

 join_trees <- function(adjacency.network.list){ #Tried to parallelize it...but no parallel methods on array chunks (also unsure how to recombine the chunks afterwards)
   adjacency.matrix.list <- lapply(adjacency.network.list, function(x) x$adjacency_matrix)
   network.dfs.list <- lapply(adjacency.network.list, function(x) x$network_df)

   if(forest.method=='multiple.germlines' | forest.method=='multiple.germlines.joined'){
     output_df <- do.call('rbind', network.dfs.list)
     nodes_per_network <- unname(lapply(adjacency.matrix.list, function(x) nrow(x)))
     output_matrix_nodes <- sum(unlist(nodes_per_network))
     output_matrix <- matrix(0, output_matrix_nodes, output_matrix_nodes)
     #output_df$sequence_id <- 1:nrow(output_df)

     germline_nodes <- which(output_df$germline=='yes')
     m1 <- adjacency.matrix.list[[1]]
     prev_added_node <- 1
     current_node <- nrow(m1)
     output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m1[1:nrow(m1), 1:nrow(m1)]
     prev_added_node <- current_node + 1
     m1 <- NULL

     for(m in adjacency.matrix.list[-1]){
       current_node <- prev_added_node + nrow(m) - 1
       output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m
       prev_added_node <- current_node + 1
       m<-NULL
     }

     if(forest.method=='multiple.germlines.joined'){
       germline_nodes <- which(output_df$germline=='yes')
       germlines <- output_df$network_sequences[germline_nodes]

       distance_matrix <- stringdist::stringdistmatrix(germlines, germlines, method='lv')

       for(i in 1:(length(germline_nodes)-1)){
         if(!weighted.germline){
           edge_weight <- distance_matrix[i,i+1]
         }else{
           edge_weight <- 1
         }

         output_matrix[germline_nodes[i], germline_nodes[i+1]] <- edge_weight
         if(!directed){
           output_matrix[germline_nodes[i+1], germline_nodes[i]] <- edge_weight
         }
       }
     }

   }else if(forest.method=='single.germline'){
     output_df <- do.call('rbind', network.dfs.list)

     nodes_per_network <- lapply(adjacency.matrix.list, function(x) nrow(x))
     germlines_per_network <- lapply(network.dfs.list, function(x) any(x$germline=='yes'))

     output_matrix_nodes <- sum(unlist(nodes_per_network)) - length(which(germlines_per_network==T)) + 1
     output_matrix <- matrix(0, output_matrix_nodes, output_matrix_nodes)

     germline_nodes_new <- which(output_df$germline=='yes')
     germline_nodes_original <- as.integer(output_df$sequence_id[germline_nodes_new])

     germline_row <- output_df[germline_nodes_new[1], ]
     rownames(germline_row) <- NULL
     output_df <- output_df[which(output_df$germline=='no' | is.na(output_df$germline)),]
     rownames(output_df) <- NULL
     output_df <- rbind(output_df, germline_row)

     new_germline_node <- nrow(output_matrix)
     new_germline_sequence <- output_df$network_sequences[new_germline_node]

     nodes_connected_to_germline <- mapply(function(x,y) which(x[y,]!=0), adjacency.matrix.list, germline_nodes_original)
     sequences_connected_to_germline <- mapply(function(x,y) x$network_sequences[y], network.dfs.list, nodes_connected_to_germline)
     new_edge_weights <- lapply(sequences_connected_to_germline, function(x) stringdist::stringdistmatrix(unlist(x), new_germline_sequence))

     nodes_connected_to_germline
     m1 <- adjacency.matrix.list[[1]]
     prev_added_node <- 1
     current_node <- nrow(m1)
     current_node <- current_node - 1
     output_matrix[prev_added_node:current_node, prev_added_node:current_node] <- m1[1:(nrow(m1)-1), 1:(nrow(m1)-1)]
     prev_added_node <- current_node + 1
     m1 <- NULL
     matrices <- adjacency.matrix.list[-1]
     output_matrix[new_germline_node,nodes_connected_to_germline[[1]]] <- new_edge_weights[[1]]

     if(!directed){
       output_matrix[nodes_connected_to_germline[[1]], new_germline_node] <- new_edge_weights[[1]]
     }
     #matrices[[i]][1:(nrow(matrices[[i]])-1), 1:(nrow(matrices[[i]])-1)]
     for(i in 1:length(matrices)){
       current_node <- prev_added_node + nrow(matrices[[i]]) - 1
       current_node <- current_node - 1
       output_matrix[prev_added_node:current_node, prev_added_node:current_node] <-  matrices[[i]][1:(nrow(matrices[[i]])-1), 1:(nrow(matrices[[i]])-1)]

       output_matrix[new_germline_node, unlist(nodes_connected_to_germline[i+1]) + prev_added_node - 1] <- new_edge_weights[[i+1]]
       if(!directed){
         output_matrix[unlist(nodes_connected_to_germline[i+1]) + prev_added_node - 1, new_germline_node] <- new_edge_weights[[i+1]]
       }
       prev_added_node <- current_node + 1
     }
   }

   return(list(adjacency_matrix = output_matrix, network_df = output_df))
 }

 create_network <- function(adjacency.network.list){
   adjacency_matrix <- adjacency.network.list$adjacency_matrix
   network_df <- adjacency.network.list$network_df
   if(!directed){
     g<-igraph::graph_from_adjacency_matrix(adjacency_matrix, mode='undirected', weighted=T, diag=F)
   }else{
     g<-igraph::graph_from_adjacency_matrix(adjacency_matrix, mode='directed', weighted=T, diag=F)
   }

   degrees <- igraph::degree(g, mode='all')

   g<-igraph::set_vertex_attr(g, name='label', index=1:nrow(network_df), value=as.factor(1:nrow(network_df)))
   g<-igraph::set_vertex_attr(g, name='network_sequences', index=which(network_df$germline=='no'), value=unlist(network_df$network_sequences[which(network_df$germline=='no')]))
   g<-igraph::set_vertex_attr(g, name='cell_number', index=1:nrow(network_df), value=unlist(network_df$cell_number))
   g<-igraph::set_vertex_attr(g, name='node_type', index=which(network_df$germline=='no'), value='sequence')
   g<-igraph::set_vertex_attr(g, name='most_expanded', index=1:nrow(network_df), value=unlist(network_df$most_expanded))
   g<-igraph::set_vertex_attr(g, name='hub', index=which.max(degrees), value='yes')
   g<-igraph::set_vertex_attr(g, name='hub', index=which(is.na(igraph::V(g)$hub)), value='no')


   #Also add cell numbers (will see if adding cell_barcodes is necessary)
   g<-igraph::set_vertex_attr(g, name='clonotype_id', index=1:nrow(network_df), value=unlist(network_df$clonotype_id))
   g<-igraph::set_vertex_attr(g, name='sample_id', index=1:nrow(network_df), value=unlist(network_df$sample_id))
   g<-igraph::set_vertex_attr(g, name='cell_barcodes', index=1:nrow(network_df), value=network_df$cell_barcodes)

   if(!is.null(include.germline)){
     g<-igraph::set_vertex_attr(g, name='node_type', index=which(network_df$germline=='yes'), value='germline')
     g<-igraph::set_vertex_attr(g, name='distance_from_germline', index=1:nrow(network_df), value=unlist(network_df$distance_from_germline))
     g<-igraph::set_vertex_attr(g, name='network_sequences', index=which(network_df$germline=='yes'), value=unlist(network_df$network_sequences[which(network_df$germline=='yes')]))
     g<-igraph::set_vertex_attr(g, name='cell_number', index=which(network_df$germline=='yes'), value=1)


   }

   if(any(network_df$sequence_id == 'intermediate')){
     g<-igraph::set_vertex_attr(g, name='node_type', index=which(network_df$sequence_id=='intermediate'), value='intermediate')
     g<-igraph::set_vertex_attr(g, name='label', index=which(igraph::V(g)$node_type=='intermediate'), value=NA)
     g<-igraph::set_vertex_attr(g, name='cell_number', index=which(igraph::V(g)$node_type=='intermediate'), value=1)

   }

   #Add leaf annotation
   if(network.algorithm=='tree'){
     g<-igraph::set_vertex_attr(g, name='leaf', index=which(degrees==1), value='leaf')
     g<-igraph::set_vertex_attr(g, name='leaf', index=which(igraph::V(g)$node_type=='germline'), value='germline')
     g<-igraph::set_vertex_attr(g, name='leaf', index=which(is.na(igraph::V(g)$leaf)), value='internal')
   }



   for(feature in features_to_select){
     if(!(feature %in% c('sample_id', 'clonotype_id'))){
       for(i in 1:nrow(network_df)){
         g<-igraph::set_vertex_attr(g, name=feature, index=i, value=network_df[,feature][i])
         g<-igraph::set_vertex_attr(g, name=paste0(feature, '_counts'), index=i, value=network_df[,paste0(feature, '_counts')][i])
       }

       if(!is.null(include.germline)){
         germline_node <- which(igraph::V(g)$node_type=='germline')
         g<-igraph::set_vertex_attr(g, name=paste0(feature), index=germline_node, value='germline')
         g<-igraph::set_vertex_attr(g, name=paste0(feature,'_counts'), index=germline_node, value=1)

       }

       if(expand.intermediates){
         g<-igraph::set_vertex_attr(g, name=paste0(feature), index=which(igraph::V(g)$node_type=='intermediate'), value='intermediate')
         g<-igraph::set_vertex_attr(g, name=paste0(feature, '_counts'), index=which(igraph::V(g)$node_type=='intermediate'), value=1)
       }
     }
   }


   return(g)
 }

 instantiate_tree_object <- function(g){

   node_df <- igraph::as_data_frame(g, what = 'vertices')
   germline_sequence <- node_df$network_sequences[node_df$node_type == 'germline']
   node_df_subset <- node_df[node_df$node_type == 'sequence',]
   #edge_df <- igraph::as_data_frame(g, what = 'edges')

   labels <- node_df$labels

   network_sequences <- node_df_subset$network_sequences
   names(network_sequences) <- labels

   cell_barcodes <- node_df_subset$cell_barcodes
   names(cell_barcodes) <- labels

   sample_id <- unlist(unique(node_df$sample_id))
   sample_id <- sample_id[order(nchar(sample_id), sample_id)]
   clonotype_id <- unlist(unique(node_df$clonotype_id))
   clonotype_id <- clonotype_id[order(nchar(clonotype_id), clonotype_id)]

   ####optional
   node_df$network_sequences <- NULL
   node_df$cell_barcodes <- NULL
   node_df$sample_id <- NULL
   node_df$clonotype_id <- NULL
   #node_df$node_type <- NULL

   list_column_names <- names(sapply(node_df, class)[which(sapply(node_df, class)=='list')])
   list_column_names <- list_column_names[stringr::str_detect(list_column_names, '_counts')]
   list_column_names <- stringr::str_remove(list_column_names, pattern='_counts')

   for(col in list_column_names){
     node_df <- node_df[,!(names(node_df) %in% paste0(col, '_counts'))]
   }

   if(!as.igraph){
     requireNamespace('tidygraph')
     g <- tidygraph::as_tbl_graph(g)
   }

   final_object <- methods::new('AntibodyForests',
                       tree = g,
                       sample_id = sample_id,
                       clonotype_id = clonotype_id,
                       sequences = network_sequences,
                       barcodes = cell_barcodes,
                       node_features = node_df,
                       germline_sequence = germline_sequence,
                       #edge_list = edge_df,
                       #adjacency_matrix = igraph::as_adjacency_matrix(g),
                       network_algorithm = network.algorithm,
                       feature_names = unlist(node.features)
                      )


   return(final_object)

 }

 set.seed(random.seed)

 if(class(VDJ)=='data.frame'){
   VDJ.GEX.matrix <- list()
   VDJ.GEX.matrix[[1]] <- VDJ
   VDJ <- NULL

 }else if(class(VDJ)=='list'){
    VDJ.GEX.matrix <- list()
    for(i in 1:length(VDJ)){
      VDJ.GEX.matrix[[i]] <- do.call('rbind', VDJ[[i]])
    }
    VDJ.GEX.matrix[[1]] <- do.call('rbind', VDJ.GEX.matrix)
 }

  sample_dfs <- list()
  output_list <- list()
  global_adjacency_network_list <- list()
  sample_names <- c()

  if(network.level!='global.clonotype' & network.level!='global'){
   repertoire.number <- unique(VDJ.GEX.matrix[[1]]$sample_id)

   for(i in 1:length(repertoire.number)){
     sample_dfs[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==repertoire.number[i]),]
     sample_dfs[[i]] <- sample_dfs[[i]]
   }

  }else{
   sample_dfs[[1]] <-VDJ.GEX.matrix[[1]]
   repertoire.number <- 'global'
  }


  if(!is.null(resolve.ties)){
    if(resolve.ties[length(resolve.ties)]!='random' & resolve.ties[length(resolve.ties)]!='first'){
      resolve.ties <- c(resolve.ties, 'first')
    }
  }

  for(i in 1:length(sample_dfs)){

   sample_dfs[[i]] <- sample_dfs[[i]][order(sample_dfs[[i]]$clonotype_frequency, decreasing=T), ]

   if(!(is.numeric(specific.networks))){
    if(specific.networks[1]!='all'){
     unique_clonotypes <- unique(specific.networks)
     sample_dfs[[i]] <- sample_dfs[[i]][which(sample_dfs[[i]]$clonotype_id %in% unique_clonotypes),]
    }
   }

   if(network.level!='global'){
     clonotype_dfs <-  split(sample_dfs[[i]], factor(sample_dfs[[i]]$clonotype_id, levels=unique(sample_dfs[[i]]$clonotype_id)))
   }else{
     clonotype_dfs <- sample_dfs[[i]]
     if(is.numeric(specific.networks)){
       selected_clonotypes <- unique(clonotype_dfs$clonotype_id)

       if(length(selected_clonotypes) > specific.networks){
         selected_clonotypes <- selected_clonotypes[1:specific.networks]
         clonotype_dfs <- clonotype_dfs[which(clonotype_dfs$clonotype_id %in% selected_clonotypes),]
       }
     }
   }
   if(class(clonotype_dfs)!='list'){
     temp_list <- list()
     temp_list[[1]] <- clonotype_dfs
     clonotype_dfs <- temp_list
     temp_list <- NULL
   }

   if(parallel){
     requireNamespace('parallel')
     cores <- parallel::detectCores()

     network_dfs <- parallel::mclapply(clonotype_dfs, transform_clonotype_to_network_df, mc.cores=cores)
     network_dfs <- network_dfs[!sapply(network_dfs, is.null)]

     if(is.null(network_dfs) | length(network_dfs)==0){
       next
     }

     if(is.numeric(specific.networks)){
       if(specific.networks < length(network_dfs)){
         network_dfs <- network_dfs[1:specific.networks]
       }
     }

     clonotype_names <- lapply(network_dfs, function(x) paste0(unique(unlist(x$clonotype_id))[order(nchar(unlist(unique(x$clonotype_id))),unlist(unique(x$clonotype_id)))], collapse=';'))
     sample_names[i] <- paste0(unlist(unique(lapply(network_dfs, function(x) unique(x$sample_id)))), collapse=';')

     #return(network_dfs)
     if(network.algorithm=='tree'){
       adjacency_matrices <- parallel::mclapply(network_dfs, calculate_adjacency_matrix_tree, mc.cores=cores)
       adjacency_network_list <- mapply(function(x,y) list(adjacency_matrix=x, network_df=y), adjacency_matrices, network_dfs, SIMPLIFY=F)
       adjacency_matrices <- NULL
       network_dfs <- NULL

     }else if(stringr::str_detect(network.algorithm, 'prune')){
       directed <- F
       adjacency_matrices <- parallel::mclapply(network_dfs, calculate_adjacency_matrix_prune, mc.cores=cores)
       adjacency_network_list <- mapply(function(x,y) list(adjacency_matrix=x, network_df=y), adjacency_matrices, network_dfs, SIMPLIFY=F)
       adjacency_matrices <- NULL
       network_dfs <- NULL

     }else if(stringr::str_detect(network.algorithm, 'phylogenetic')){
       adjacency_network_list <- parallel::mclapply(network_dfs, create_phylo_trees, mc.cores=cores)
       unlink('./tempdir', recursive=T)
       network_dfs <- NULL

     }else if(network.algorithm == 'mst'){
       directed <- F
       adjacency_network_list <- parallel::mclapply(network_dfs, create_mst_trees, mc.cores=cores)
       network_dfs <- NULL

     }else{
       stop('Unknown network.algorithm!')
     }

     if(expand.intermediates & !(stringr::str_detect(network.algorithm, 'phylogenetic'))) {
       adjacency_network_list <- parallel::mclapply(adjacency_network_list, expand_intermediates, mc.cores=cores)
     }

     if(network.level=='forest.per.sample'){
       adjacency_network_list <- list(join_trees(adjacency_network_list))
       #clonotype_names <- paste0(clonotype_names, collapse=';')
       clonotype_names <- 'joined.per.sample'
     }

     if(network.level=='forest.global'){
       global_adjacency_network_list[[i]] <- adjacency_network_list
       next
     }

     output_list[[i]] <- parallel::mclapply(adjacency_network_list, create_network, mc.cores=cores)
     output_list[[i]] <- parallel::mclapply(output_list[[i]], instantiate_tree_object, mc.cores=cores)

     #output_list[[i]] <- unname(output_list[[i]])
     #output_list[[i]] <- unlist(output_list[[i]])
   }else{

     network_dfs <- lapply(clonotype_dfs, function(x) transform_clonotype_to_network_df(x))
     network_dfs <- network_dfs[!sapply(network_dfs, is.null)]

     if(is.null(network_dfs) | length(network_dfs)==0){
       next
     }

     if(is.numeric(specific.networks)){
       if(specific.networks < length(network_dfs)){
         network_dfs <- network_dfs[1:specific.networks]
       }
     }

     clonotype_names <- lapply(network_dfs, function(x) paste0(unique(unlist(x$clonotype_id))[order(nchar(unique(unlist(x$clonotype_id))),unique(unlist(x$clonotype_id)))], collapse=';'))
     sample_names[i] <- paste0(unlist(unique(lapply(network_dfs, function(x) unique(x$sample_id)))), collapse=';')

     if(network.algorithm=='tree'){
       adjacency_matrices <- lapply(network_dfs, function(x) calculate_adjacency_matrix_tree(x))
       adjacency_network_list <- mapply(function(x,y) list(adjacency_matrix=x, network_df=y), adjacency_matrices, network_dfs, SIMPLIFY=F)
       adjacency_matrices <- NULL
       network_dfs <- NULL

     }else if(stringr::str_detect(network.algorithm, 'prune')){
       adjacency_matrices <- lapply(network_dfs, function(x) calculate_adjacency_matrix_prune(x))
       adjacency_network_list <- mapply(function(x,y) list(adjacency_matrix=x, network_df=y), adjacency_matrices, network_dfs, SIMPLIFY=F)
       adjacency_matrices <- NULL
       network_dfs <- NULL

     }else if(stringr::str_detect(network.algorithm, 'phylogenetic')){
       adjacency_network_list <- lapply(network_dfs, function(x) create_phylo_trees(x))
       unlink('./tempdir', recursive=T)
       network_dfs <- NULL

     }else if(network.algorithm == 'mst'){
       directed <- F
       adjacency_network_list <- lapply(network_dfs, function(x) create_mst_trees(x))
       network_dfs <- NULL

     }else{
       stop('Unknown network.algorithm!')
     }

     if(expand.intermediates & !(stringr::str_detect(network.algorithm, 'phylogenetic'))){
       adjacency_network_list <- lapply(adjacency_network_list, function(x) expand_intermediates(x))
     }

     if(network.level=='forest.per.sample'){
       adjacency_network_list <- list(join_trees(adjacency_network_list))
       #clonotype_names <- paste0(clonotype_names, collapse=';')
       clonotype_names <- 'joined.per.sample'
     }

     if(network.level=='forest.global'){
       global_adjacency_network_list[[i]] <- adjacency_network_list
       next
     }

     output_list[[i]] <- lapply(adjacency_network_list, function(x) create_network(x))
     output_list[[i]] <- lapply(output_list[[i]], function(x) instantiate_tree_object(x))
   }

   names(output_list[[i]]) <- clonotype_names
   output_list[[i]] <- output_list[[i]][order(nchar(names(output_list[[i]])),names(output_list[[i]]))]

 }

 if(network.level=='forest.global'){
   adjacency_network_list <- purrr::flatten(global_adjacency_network_list)

   adjacency_network_list <- join_trees(adjacency_network_list)
   #names(adjacency_network_list) <- paste0(clonotype_names, collapse='-')

   output_list <- create_network(adjacency_network_list)
   output_list <- list(instantiate_tree_object(output_list))
   names(output_list) <- paste0(sample_names[!is.na(sample_names)], collapse=';')
   output_list <- list(output_list)
   names(output_list) <- paste0(sample_names[!is.na(sample_names)], collapse=';')
   return(output_list)
 }

 names(output_list) <- sample_names
 output_list <- output_list[!sapply(output_list, is.null)]

 if(length(output_list) == 1){
   if(length(output_list[[1]]) == 1){
     output_list <- output_list[[1]][[1]]
   }
 }

 #output_list <- purrr::flatten(output_list)

 return(output_list)
}
