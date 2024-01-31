#' Infer B cell evolutionary networks for all clonotypes in VDJ dataframe as obtained from the minimal_VDJ() function.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description 
#' @param VDJ dataframe - VDJ object as obtained from the minimal_VDJ() function in Platypus.
#' @param sequence.columns string or vector of strings - denotes the sequence column(s) in the VDJ dataframe that contain the sequences that will be used to infer B cell lineage trees. Nodes in the trees will represent unique combinations of the selected sequences.
#' @param germline.columns string or vector of strings - denotes the germline column(s) in the VDJ dataframe that contain the sequences that will be used as starting points of the trees. The columns should be in the same order as in 'sequence.columns' and should have an equal length.
#' @param node.features string or vector of strings - denotes the column name(s) in the VDJ dataframe from which the node features should be extracted (which can, for example, be used for plotting of lineage trees later on).
#' @param distance.metric string or vector of strings - denotes the metric that will be calculated with the 'stringdist()' function to measure distance between sequences (when no substitution model is specified). Possible values: 'osa', 'lv', 'dl', 'hamming', 'lcs', 'qgram', 'cosine', and 'soundex'.  Defaults to 'lv' (Levenshtein distance / edit distance).
#' @param substitution.model string or vector of strings - denotes the nucleotide, codon, or amino acid substitution model that will be used to calculate an alignment score of two sequences. Possible values: 'GY94' and 'HLP17'. 
#' @param network.algorithm string - denotes the network algorithm that will be used to convert the distance/score matrix into a lineage tree. Defaults to 'tree' (mst-like algorithm in which the germline node is positioned at the top of the tree, and nodes with the minimum distance to any existing node in the tree are linked iteratively).
#' @param resolve.ties string - denotes the way ties are handled during the conversion of distance matrices into lineage trees using the 'tree' network algorithm. Defaults to 'max.expansion' (when an unconnected clone/node has the same distance to multiple connected clones/nodes in the graph, the node is added to the most expanded/biggest clone/node).
#' @param parallel bool - if TRUE, the per-clone network inference is executed in parallel (parallelized across samples). Defaults to FALSE.
#' @param num.cores integer - number of cores to be used when parallel = TRUE. Defaults to all available cores - 1 or the number of samples in the VDJ dataframe (depending which number is smaller).
#' @return Returns nested list of AntibodyForests objects for each sample and each clonotype, whereby an AntibodyForests object consists of a list of all unique sequences within one clonotype, a list of distance matrices (one for each selected sequence column), and the lineage tree as an igraph object. For example, output[[1]][[2]] denotes the AntibodyForests object of the first sample, second clonotype. 
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, 
#'                 sequence.columns = c("VDJ_sequence_aa_trimmed","VJ_sequence_aa_trimmed"),
#'                 germline.columns = c("VDJ_germline_aa_trimmed","VJ_germline_aa_trimmed"), 
#'                 node.features = c("VDJ_vgene", "VDJ_dgene", "VDJ_jgene", "VDJ_cgene"),
#'                 distance.metric = "lv",
#'                 network.algorithm = "tree",
#'                 resolve.ties = "max.expansion",
#'                 parallel = TRUE)
#'}


AntibodyForests <- function(VDJ,
                            sequence.columns,
                            germline.columns,
                            node.features,
                            distance.metric,
                            substitution.model,
                            network.algorithm,
                            resolve.ties,
                            parallel = FALSE,
                            num.cores = NULL){
  
  # If 'sequence.columns' is missing or if the columns are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(missing(sequence.columns) | !all(sequence.columns %in% colnames(VDJ))){stop("Error: Please provide valid 'sequence.columns' from the input VDJ dataframe.")}
  # If 'germline.columns' is missing or if the columns are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(missing(germline.columns) | !all(germline.columns %in% colnames(VDJ))){stop("Error: Please provide valid 'germline.columns' from the input VDJ dataframe.")}
  # If the columns in 'germline.columns' do not correspond to the columns in 'sequence.columns', a message is returned and execution is stopped
  if(!all(sequence.columns == gsub("germline", "sequence", germline.columns))){stop("Error: Please make sure that the columns in 'germline.columns' correspond to the columns in 'sequence.columns'.")}
  # If 'node.features' contains columns that are not present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(node.features %in% colnames(VDJ))){stop("Error: Not all columns in 'node.features' could be found in the input VDJ dataframe.")}
  # If 'distance.metric', 'substitution.model', and 'network.algorithm' are all missing, the edit distance / Levenshtein distance will be used to calculate distance matrices and the 'tree' algorithm will be used to convert the distance matrices into lineage trees
  if(missing(distance.metric) && missing(substitution.model) && missing(network.algorithm)){distance.metric <- "lv"; network.algorithm <- "tree"}
  # If 'network.algorithm' is set to 'tree', but t
  if(network.algorithm == "tree" && missing(resolve.ties)){resolve.ties <- "max.expansion"}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  
  calculate_edit_score <- function(seq1, seq2){
    
    # Calculates edit distance between two sequences 
    # Function is based on course 'Algorithms for DNA Sequencing', offered by Coursera (https://www.coursera.org/learn/dna-sequencing)
    # Arguments:
    # - seq1: first input sequence
    # - seq2: second input sequence
    # Author: Valentijn Tromp
    
    # Split sequences into lists with individual characters
    seq1 <- unlist(strsplit(seq1, split=""))
    seq2 <- unlist(strsplit(seq2, split=""))
    
    # Create a matrix to store edit distances
    edit_dist_matrix <- matrix(0, nrow = length(seq1)+1, ncol = length(seq2)+1)
    
    # Set row and column names of the matrix
    rownames(edit_dist_matrix) <- c("ε",seq1)
    colnames(edit_dist_matrix) <- c("ε",seq2)
    
    # Initialize the first row of the matrix with sequential numbers
    for(i in 1:(length(seq1)+1)){
      edit_dist_matrix[i,1] <- i-1
    }
    
    # Initialize th first column of the matrix with sequential numbers
    for(i in 1:(length(seq2)+1)){
      edit_dist_matrix[1,i] <- i-1
    }
    
    # Loop over positions in sequence 1 (+1)
    for(i in 2:(length(seq1)+1)){
      # Loop over positions in sequence 2 (+1)
      for(j in 2:(length(seq2)+1)){
        
        # Calculate the horizontal distance 
        dist_hor <- edit_dist_matrix[i,(j-1)] + 1
        
        # Calculate the vertical distance
        dist_ver <- edit_dist_matrix[(i-1),j] + 1
        
        # If the characters at the current positions are equal, no additional costs is given to the diagonal distance
        if(seq1[i-1] == seq2[j-1]){
          dist_diag <- edit_dist_matrix[(i-1),(j-1)]
        }
        
        # If the characters are not equal, a costs is added to the diagonal distance
        if(seq1[i-1] != seq2[j-1]){
          dist_diag <- edit_dist_matrix[(i-1),(j-1)] + 1
        }
        
        # Update the edit distance matrix with the minimum of the three distances
        edit_dist_matrix[i,j] = min(dist_hor, dist_ver, dist_diag)
      }
    }
    
    # Return edit sequence (value bottom right in the matrix)
    return(edit_dist_matrix[length(seq1)+1, length(seq2)+1])
  }
  
  
  convert_dist2tree <- function(dist.matrices, resolve.ties){
    
    # Converts list of distance matrices to lineage tree using mst-like algorithm
    # The germline node is positioned at the top of the tree, and nodes with the minimum distance to any existing node in the tree are linked iteratively
    # Arguments:
    # - dist.matrices: list of distance matrices (second object in AntibodyForests list)
    # - resolve.ties: string denoting how ties should be handled (options: "max.expansion", "min.germline.distance")
    # Author: Valentijn Tromp
    
    # Sum up the distance matrices in input list 'dist.matrices'
    dist_matrix_sum <- Reduce(`+`, dist.matrices)
    
    # Duplicate 'dist_matrix_sum' for later on
    dist_matrix_sum_backup <- dist_matrix_sum
    
    # Replace diagonal 0's by NA
    diag(dist_matrix_sum) <- NA
    
    # Create matrix to store edges of the tree 
    edge_matrix <- matrix(data = NA, nrow = nrow(dist_matrix_sum)-1, ncol = 2)
    
    # Get minimum distance to the germline node
    min_distance <- min(na.exclude(dist_matrix_sum[,"germline"]))
    
    # Select nodes that have this distance to the germline node
    nodes_to_be_connected <- rownames(dist_matrix_sum)[dist_matrix_sum[, "germline"] == min_distance & !is.na(dist_matrix_sum[, "germline"])]
    
    # Add edge(s) between germline and node(s) with 'min_distance' to the germline to 'edge_matrix'
    edge_matrix[1:length(nodes_to_be_connected),1] <- "germline"
    edge_matrix[1:length(nodes_to_be_connected),2] <- nodes_to_be_connected
    
    # Remove minimum distance to the germline node from the distance matrix
    dist_matrix_sum[,"germline"][dist_matrix_sum[,"germline"] == min_distance] <- NA
    dist_matrix_sum["germline",][dist_matrix_sum["germline",] == min_distance] <- NA
    
    # Iteratively add other nodes to 'edge_matrix' while not all nodes are in 'edge_matrix'
    while(all(rownames(dist_matrix_sum) %in% edge_matrix) == FALSE){
      
      # Select nodes that are already connected (present in 'edge_matrix')
      already_connected_nodes <- unique(na.exclude(c(edge_matrix[,1], edge_matrix[,2])))
      
      # Get minimum distance to one of these nodes
      min_distance <- min(na.exclude(dist_matrix_sum[,already_connected_nodes]))
      
      # Select unlinked nodes that have this distance to any existing node in 'already_connected_nodes'
      nodes_to_be_connected <- unique(unlist(lapply(rownames(dist_matrix_sum), function(x) if(min_distance %in% dist_matrix_sum[x,already_connected_nodes] && !(x %in% already_connected_nodes)){return(x)})))
      
      # Iterate through nodes in 'nodes_to_be_connected'
      for(node in nodes_to_be_connected){
        
        # Select the node(s) in 'already_connected_nodes' that has 'min_distance' to the current node
        node_to_connect_to <- unlist(lapply(already_connected_nodes, function(x) if(dist_matrix_sum[x, node] == min_distance){return(x)}))
        
        # If there are multiple nodes in 'already_connected_nodes' that has 'min_distance' to the current node, the 'resolve.ties' parameter comes into play
        if(length(node_to_connect_to) != 1){
          
          # If the 'resolve.ties' parameter is set to 'max.expansion', the current node is added to the node with the smallest number (assuming numerical labels, where nodes/clones are numbered in descending order based on their size)
          if(resolve.ties == "max.expansion"){
            node_to_connect_to <- paste0("node", min(as.numeric(gsub("node", "", node_to_connect_to))))
          }
          
          # If the 'resolve.ties' parameter is set to 'close.germline.distance', the current node is added to the node that has the smallest string distance to the germline node
          if(resolve.ties == "min.germline.distance"){
            node_to_connect_to <- node_to_connect_to[dist_matrix_sum_backup[node_to_connect_to, "germline"] == min(dist_matrix_sum_backup[node_to_connect_to, "germline"])]
            }
          }
        
        # Add edge between node that this already connected and node that has 'min_distance' to this node to 'edge_matrix'
        edge_matrix[length(already_connected_nodes)+which(nodes_to_be_connected == node)-1,] <- c(node_to_connect_to, node)
        
        # Remove this minimum distance from the distance matrix
        dist_matrix_sum[node_to_connect_to, node] <- NA
        dist_matrix_sum[node, node_to_connect_to] <- NA
      }
    }
    
    # Create igraph object by converting 'edge_matrix' into a graph
    tree <- igraph::graph_from_edgelist(edge_matrix)
    
    # Return tree
    return(tree)
  }
  
  infer_single_network <- function(VDJ,
                                   clone,
                                   sequence.columns,
                                   germline.columns,
                                   node.features,
                                   distance.metric,
                                   substitution.model,
                                   network.algorithm,
                                   resolve.ties){
    
    # Infers network/tree for a single clone within one sample
    # Arguments:
    # - VDJ: VDJ dataframe as obtained from the minimal_VDJ() function in Platypus
    # - clone: string denoting the sample ID and the clonotype ID for which the network will be inferred (for example: "S1_clonotype4" or "S2_clonotype3")
    # - sequence.columns: string or vector of strings denoting the sequence columns in the 'VDJ' dataframe that contain the sequences that will be used to infer B cell lineage trees
    # - germline.columns: string or vector of strings denoting the germline columns in the 'VDJ' dataframe that contain the sequences that will be used as starting points of the lineage trees
    # - node.features: string or vector of strings denoting the additional column name(s) in the VDJ dataframe to be selected
    # - distance.metric: string denoting the metric that will be calculated to measure (string) distances between sequences (if no substitution model is specified)
    # - substitution.model: denotes substitution model that will be used to calculate the alignment score when comparing sequences
    # - network.algorithm: string denoting the network algorithm that will be used to convert the distance/score matrices into lineage trees
    # - resolve.ties: string denoting the way ties are handled when the 'tree' network algorithm is being used 
    
    # 1. Create a nested list 'nodes_list', in which each sublist represents a clone/node in the graph, and each node represent a unique combination of sequences of the columns selected in 'sequence.columns'
    # Each node/sublist contains the following items:
    # - a (separate) string for all columns selected in 'sequence.columns' containing the sequence for that node
    # - barcodes: a list of all barcodes of the cells that have the unique combination of sequences
    # - frequency: an integer that corresponds to the number of cells (equal to the length of the list 'barcodes')
    # - a (separate) list for all columns selected in 'node.features'
    
    # Retrieve sample ID and clonotype ID from 'clone'
    sample <- strsplit(clone, split="_")[[1]][1]
    clonotype <- strsplit(clone, split="_")[[1]][2]

    # Make a subset of the VDJ dataframe with all cells from the specified sample and clone and only the columns 'barcode', the specified 'sequence.columns' and 'germline.columns', and additional columns specified in 'node.features'
    vdj_subset <- VDJ[VDJ$sample_id == sample & VDJ$clonotype_id == clonotype,c("barcode", sequence.columns, germline.columns, node.features)]
    
    # Create vectors specifying
    nt.columns <- colnames(vdj_subset)[stringr::str_detect(colnames(vdj_subset), "_nt")]
    aa.columns <- colnames(vdj_subset)[stringr::str_detect(colnames(vdj_subset), "_aa")]
    
    # Trim of last one or two nucleotides to make the length of the nucleotide sequences in 'vdj_subset' multiple of 3
    vdj_subset[nt.columns] <- lapply(vdj_subset[nt.columns], function(x){
      ifelse(nchar(x) %% 3 != 0, substring(x, 1, nchar(x) - (nchar(x) %% 3)))
    })
    
    # Create a dataframe containing unique combinations of the sequences in 'sequence.columns'
    intraclonal_sequences_df <- unique(vdj_subset[,sequence.columns])
    
    # Create a list of sublists where each sublist represents a unique combination of sequences (and each unique combination of sequences will represent a node in the graph later on)
    nodes_list <- lapply(1:nrow(intraclonal_sequences_df), function(x){
      # Extract sequences for the current row and save in 'node_sublist'
      node_sublist <- lapply(sequence.columns, function(y) intraclonal_sequences_df[x,y])
      # Rename indices of 'node_sublist' to 'sequence.columns'
      names(node_sublist) <- sequence.columns
      # Find matching rows in 'vdj_subset' based on the current combination of sequences
      matching_rows <- sapply(1:nrow(vdj_subset), function(y) all(intraclonal_sequences_df[x,] == vdj_subset[y,sequence.columns])) 
      # Add barcodes and frequency information to 'node_sublist'
      node_sublist["barcodes"] <- list(vdj_subset$barcode[matching_rows])
      node_sublist["frequency"] <- sum(matching_rows)
      # Add 'node.features' to 'node_sublist' as well
      for(i in node.features){node_sublist[i] <- list(vdj_subset[,i][matching_rows])}
      # Return list of sublists 
      return(node_sublist)
    })
    
    # Order the list based on the frequency of the (sub)clones/nodes in descending order
    nodes_list <- nodes_list[order(sapply(nodes_list, function(x) x$frequency), decreasing = TRUE)]
    
    # Rename the sublists to 'node1', 'node2', etc.
    names(nodes_list) <- paste0("node", seq_along(nodes_list))
    
    # Add a 'germline' node to the list, containing the the most abundant germline sequences from the specified 'germline.columns' (most of the time, each germline column contains only one unique sequence)
    nodes_list$germline <- lapply(germline.columns, function(x) names(table(vdj_subset[,x]))[which.max(table(vdj_subset[,x]))])
    names(nodes_list$germline) <- sequence.columns
    
    
    # 2. Calculate a (string) distance matrix for each column in 'sequence.columns' using the specified 'distance.metric', if no 'substitution.model' is provided
    
    # If the provided 'distance.metric' is one of the methods for distance calculation of the 'stringdistmatrix()' function, this function is used to calculate the distance matrix with the specified method
    if(distance.metric %in% c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")){
      
      # Calculate a distance matrix for each 'sequence.column' separately
      string_dist_matrices <- lapply(sequence.columns, function(x){
        
        # Extract sequences from 'nodes_list' of one 'sequence.column' and save in 'seqs' list
        seqs <- sapply(nodes_list, function(y) y[x])
        names(seqs) <- names(nodes_list)
        
        # Calculate distance matrix with 'stringdistmatrix()' function using the 'seqs' list as input 
        dist_matrix <- suppressWarnings(stringdist::stringdistmatrix(seqs, seqs, method = distance.metric, useNames = "names"))
        
        # Return the distance matrix
        return(dist_matrix)
      })
      
      # Rename distance matrices in list 'dist_matrices' according to the sequence columns in 'sequence.columns'
      names(string_dist_matrices) <- sequence.columns
    }
    
    
    # 3. Convert distance matrices to a lineage tree
    
    # If the provided 'network.algorithm' equals to 'tree', the function 'convert_dist2tree()' is used to build the lineage tree in a mst-like manner
    if(network.algorithm == "tree"){
      igraph_tree <- convert_dist2tree(dist.matrices = string_dist_matrices, 
                                       resolve.ties = resolve.ties)
    }
    
    
    # 4. Create and return a list containing the phylogenetic tree, the sequences and features of the nodes in the tree, and the distance matrices used to build the tree
    
    # Return 'nodes_list', 'dist_matrices' and 'tree' as a list
    return(list(nodes = nodes_list, 
                distance.matrices = string_dist_matrices,
                lineage.tree = igraph_tree))
  }
  
  
  # Create list with sample IDs
  sample_list <- unique(VDJ$sample_id)
  
  # Create list with clonotype IDs 
  clone_list <- unique(paste(VDJ$sample_id, VDJ$clonotype_id, sep="_"))
  
  # Define partial function for to be executed for each clone
  partial_function <- function(clone){
    single_network <- infer_single_network(VDJ = VDJ,
                                           clone = clone,
                                           sequence.columns = sequence.columns,
                                           germline.columns = germline.columns,
                                           node.features = node.features,
                                           distance.metric = distance.metric,
                                           substitution.model = substitution.model,
                                           network.algorithm = network.algorithm,
                                           resolve.ties = resolve.ties)
    return(single_network)
  }
  
  # If 'parallel' is set to TRUE, the network inference is parallelized across the clones
  if(parallel){
    
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::mclapply(mc.cores=num.cores, X=clone_list, FUN=partial_function)
      # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
      names(output_list) <- clone_list
    }
    
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::parLapply(cluster, X=clone_list, fun=partial_function)
      # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
      names(output_list) <- clone_list
      # Stop cluster
      parallel::stopCluster(cluster)
    }
    
  }
  
  # If 'parallel' is set to FALSE, the network inference is not parallelized
  if(!parallel){
    
    # Create a list for each sample and store in 'output_list'
    output_list <- lapply(clone_list, partial_function)
    # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
    names(output_list) <- clone_list
    
  }
  
  # Reorganize list
  reorganized_output_list <- lapply(sample_list, function(sample){
    # Select all clonotypes from one sample
    matching_sublists <- grep(paste0("^",sample), names(output_list), value = TRUE)
    # Retrieve the networks from 'output_list' of one sample and store in 'sample_sublist'
    sample_sublist <- output_list[matching_sublists]
    # Rename the networks in 'sample_sublist' to their original clonotype ID
    names(sample_sublist) <- sub(paste0("^", sample, "_"), "", matching_sublists)
    # Return the 'sample_sublist'
    return(sample_sublist)
  })
  # Rename the sublists in 'reorganized_output_list' to their original sample ID
  names(reorganized_output_list) <- sample_list
  
  # Return 'output_list'
  return(reorganized_output_list)
}
