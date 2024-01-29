#' Infer B cell evolutionary networks for all clonotypes in VDJ dataframe as obtained from the minimal_VDJ() function.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description 
#' @param VDJ dataframe - VDJ object as obtained from the minimal_VDJ() function in Platypus.
#' @param sequence.columns string or vector of strings - denotes the sequence column(s) in the VDJ dataframe that contain the sequences that will be used to infer B cell lineage trees. Nodes in the trees will represent unique combinations of the selected sequences.
#' @param germline.columns string or vector of strings - denotes the germline column(s) in the VDJ dataframe that contain the sequences that will be used as starting points of the trees. The columns should be in the same order as in 'sequence.columns' and should have an equal length.
#' @param node.features string or vector of strings - denotes the column name(s) in the VDJ dataframe from which the node features should be extracted (which can, for example, be used for plotting of lineage trees later on).
#' @param distance.metric string or vector of strings - denotes the metric that will be calculated with the 'stringdist()' function to measure distance between sequences (when no substitution model is specified). Possible values: 'osa', 'lv', 'dl', 'hamming', 'lcs', 'qgram', 'cosine', and 'soundex'.  Defaults to 'lv' (Levenshtein distance / edit distance).
#' @param substitution.model string or vector of strings - denotes the nucleotide, codon, or amino acid substitution model that will be used to calculate an alignment score of two sequences. Possible values: 'GY94' and 'HLP17'. 
#' @param organism string - specifies the organism of which the sequences are derived. Possible values: 'human' and 'mouse'.
#' @param network.algorithm - string - denotes the network algorithm that will be used to convert the distance/score matrix into a lineage tree. Defaults to 'mst' (Kruskal's minimum spanning tree algorithm).
#' @param parallel bool - if TRUE, the per-clone network inference is executed in parallel (parallelized across samples). Defaults to FALSE.
#' @param num.cores integer - number of cores to be used when parallel = TRUE. Defaults to all available cores - 1 or the number of samples in the VDJ dataframe (depending which number is smaller).
#' @return Returns nested list of AntibodyForests objects for each sample and each clonotype, whereby an AntibodyForests object consists of a list of all unique sequences, a distance matrix, and a phylogenetic tree (in NEWICK format). For example, output[[1]][[2]] denotes the AntibodyForests object of the first sample, second clonotype. 
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, 
#'                 sequence.columns = c("VDJ_sequence_aa_trimmed","VJ_sequence_aa_trimmed"),
#'                 germline.columns = c("VDJ_germline_aa_trimmed","VJ_sequence_aa_trimmed"), 
#'                 node.features = c("VDJ_cgene),
#'                 distance.metric = "lv",
#'                 network.algorithm = "mst")
#'}


AntibodyForests <- function(VDJ,
                            sequence.columns,
                            germline.columns,
                            node.features,
                            distance.metric,
                            substitution.model,
                            organism,
                            network.algorithm,
                            parallel = FALSE,
                            num.cores = NULL){
  
  # If 'sequence.columns' is missing or if the columns are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(missing(sequence.columns) | !all(node.features %in% colnames(VDJ))){stop("")}
  # If the number of specified sequence columns is not the same as the number of specified germline columns, a message is returned and execution is stopped
  if(length(sequence.columns) != length(germline.columns)){stop("")}
  # If 'node.features' contains columns that are not present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(node.features %in% colnames(VDJ))){stop("")}
  # If 'distance.metric' and 'scoring.matrix' are both missing, the edit distance / Levenshtein distance will be used to calculate a distance matrix 
  if(missing(distance.metric) && missing(scoring.matrix)){distance.metric <- "lv"}
  # If a 'nt.substitution.model' is provided, without specifying the organism, a message is returned and execution is stopped
  if(!missing(substitution.model)){if(substitution.model %in% c("GY94", "HLP17") && missing(organism)){stop("")}}
  # If 'num.cores' is not provided, the number of cores is set to all available cores - 1
  if(missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  
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
  
  
  convert_dist2tree <- function(dist_matrices){
    
    # Converts list of distance matrices to lineage tree using mst-like algorithm
    # Germline node is placed on top of tree, and nodes having a minimum distance to one of the nodes already present in the tree are iteratively linked
    # Arguments:
    # - dist_matrices: list of distance matrices
    # Author: Valentijn Tromp
    
    # Sum up the distance matrices in input list 'dist_matrices'
    dist_matrix_sum <- Reduce(`+`, dist_matrices)
    
    # Create matrix to store edges of the tree
    edge_matrix <- matrix(data = NA, nrow = nrow(dist_matrix_sum), ncol = 2)
    
    #...
    
    # Create igraph object by converting 'edge_matrix' into a graph
    tree <- igraph::graph_from_edgelist(edge_matrix, directed = TRUE)
    
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
                                   organism,
                                   network.algorithm){
    
    # Infers network/tree for a single clone within one sample
    # Arguments:
    # - VDJ: VDJ dataframe as obtained from the minimal_VDJ() function in Platypus
    # - clone: string denoting the clonotype ID for which the network will be inferred
    # - sequence.columns: string or vector of strings denoting the sequence columns in the 'VDJ' dataframe that contain the sequences that will be used to infer B cell lineage trees
    # - germline.columns: string or vector of strings denoting the germline columns in the 'VDJ' dataframe that contain the sequences that will be used as starting points of the trees
    # - node.features: string or vector of strings denoting the additional column name(s) in the VDJ dataframe to be selected
    # - distance.metric: string denoting the metric that will be calculated to measure distance between sequences (if no substitution model is specified)
    # - substitution.model: denotes substitution model that will be used to calculate the alignment score when comparing sequences
    # - organism: specifies the organisms of which the sequences are derived  
    # - network.algorithm: string denoting the network algorithm that will be used to convert the distance/score matrix into a lineage tree
    
    # 1. Create a nested list 'intraclonal_sequences_list', in which each sublist represents a node in the graph, and each node represent a unique combination of sequences selected in 'nt.sequence.columns' and 'aa.sequence.columns'.
    # Each node/sublists contains the following items:
    # - a (separate) character value for all columns in 'nt.sequence.columns' and 'aa.sequence.columns' containing the sequence for that node
    # - barcodes: a list of all barcodes of the cells harbouring the unique combination of sequences
    # - frequency: an integer that corresponds to the number of cells (equal to the length of the list 'barcodes')
    # - a (separate) item for all columns selected in 'node.features'
    
    # Retrieve sample ID and clonotype ID from 'clone'
    sample <- strsplit(clone, split="_")[[1]][1]
    clonotype <- strsplit(clone, split="_")[[1]][2]
    
    # Make a subset of the VDJ dataframe with all cells from the specified sample and clone and only the columns 'barcode', the specified 'sequence.columns' and 'germline.columns', and additional columns specified in 'node.features'
    vdj_subset <- subset(VDJ, sample_id == sample & clonotype_id == clonotype)[,c("barcode", sequence.columns, germline.columns, node.features)]
    
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
    intraclonal_sequences_list <- lapply(1:nrow(intraclonal_sequences_df), function(x){
      # Extract sequences for the current row and save in 'sequences' sublist
      sequences <- lapply(sequence.columns, function(y){intraclonal_sequences_df[x,y]})
      # Rename indices of sequences 
      names(sequences) <- sequence.columns
      # Find matching rows in 'vdj_subset' based on the current combination of sequences
      matching_rows <- sapply(1:nrow(vdj_subset), function(y){all(intraclonal_sequences_df[x,] == vdj_subset[y,sequence.columns])}) 
      # Add barcodes and frequency information to the 'sequences' sublist
      sequences["barcodes"] <- list(vdj_subset$barcode[matching_rows])
      sequences["frequency"] <- sum(matching_rows)
      # Add 'node.features' to the 'sequences' sublist as well
      for(i in node.features){sequences[i] <- list(vdj_subset[,i][matching_rows])}
      # Return list of sublists 
      return(sequences)
    })
    
    # Order the list based on the frequency of the cells/nodes in descending order
    intraclonal_sequences_list <- intraclonal_sequences_list[order(sapply(intraclonal_sequences_list, function(x) x$frequency), decreasing = TRUE)]
    
    # Rename the sublists to 'node1', 'node2', etc.
    names(intraclonal_sequences_list) <- paste0("node", seq_along(intraclonal_sequences_list))
    
    # Add 'germline' to the list, containing the the most abundant germline sequences from the specified 'germline.columns' (most of the time, each germline column contains only one unique sequence)
    intraclonal_sequences_list$germline <- lapply(germline.columns, function(x) names(table(vdj_subset[,x]))[which.max(table(vdj_subset[,x]))])
    names(intraclonal_sequences_list$germline) <- sequence.columns
    
    
    # 2. Calculate a distance matrix or score matrix for each column in 'sequence.columns' using the specified 'distance.metric' or 'substitution.model', respectively
    
    # Calculate a distance matrix for each 'sequence.column' separately
    dist_matrices <- lapply(sequence.columns, function(x){
      
      # Extract sequences from 'intraclonal_sequences_list' of one 'sequence.column'
      seqs <- sapply(intraclonal_sequences_list, function(y) y[x])
      names(seqs) <- names(intraclonal_sequences_list)
      
      # If the provided 'distance.metric' is one of the methods for distance calculation of the 'stringdistmatrix()' function, this function is used to calculate the distance matrix with the specified method
      if(distance.metric %in% c("osa", "lv", "dl", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw", "soundex")){
        dist_matrix <- suppressWarnings(stringdist::stringdistmatrix(seqs, seqs, method = distance.metric, useNames = "names"))
      }

      # Return the distance matrix
      return(dist_matrix)
    })
    
    # Rename distance matrices in list 'dist_matrices' according to sequence columns in 'sequence.columns'
    names(dist_matrices) <- sequence.columns
    
    
    # 3. Convert distance matrices to a lineage tree
    
    # If the provided 'network.algorithm' equals to 'tree', the function 'convert_dist2tree()' is used to build the lineage tree in a mst-like manner
    if(network.algorithm == "tree"){
      tree <- convert_dist2tree(dist_matrices)
    }
    
    
    # 4. Create and return a list containing the phylogenetic tree, the sequences and features of the nodes in the tree, and the distance matrices used to build the tree
    
    # Return 'intraclonal_sequences_list', 'dist_matrices' and 'tree' as a list
    return(list(nodes = intraclonal_sequences_list, 
                distance.matrices = dist_matrices,
                lineage.tree = tree))
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
                                           organism = organism,
                                           network.algorithm = network.algorithm)
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
    # Retrieve the networks from'output_list' of one sample and store in 'sample_sublist'
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