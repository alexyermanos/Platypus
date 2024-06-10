#' Function to infer B cell evolutionary networks for all clonotypes in VDJ dataframe as obtained from the 'minimal_VDJ()' function.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description This function takes a VDJ dataframe and uses the specified sequence columns to build a tree/network for each clonotype and stores them in an AntibodyForests object, together with the sequences and other specified features. These trees/networks provide insights into the evolutionary relationships between B cell sequences from each clonotype. The resulting object of class 'AntibodyForests' can be used for downstream analysis as input for...
#' @param VDJ dataframe - VDJ object as obtained from the minimal_VDJ() function in Platypus, or object of class dataframe that contains the columns specified in 'sequence.columns', 'germline.columns', and 'node.features'.
#' @param sequence.columns string or vector of strings - denotes the sequence column(s) in the VDJ dataframe that contain the sequences that will be used to infer B cell lineage trees. Nodes in the trees will represent unique combinations of the selected sequences. Defaults to 'c("VDJ_sequence_nt_trimmed", "VJ_sequence_nt_trimmed")'.
#' @param germline.columns string or vector of strings - denotes the germline column(s) in the VDJ dataframe that contain the sequences that will be used as starting points of the trees. The columns should be in the same order as in 'sequence.columns'. Defaults to 'c("VDJ_germline_nt_trimmed", "VJ_germline_nt_trimmed")'.
#' @param concatenate.sequences bool - if TRUE, sequences from multiple sequence columns are concatenated into one sequence for single distance matrix calculations / multiple sequence alignments, else, a distance matrix is calculated / multiple sequence alignment is performed for each sequence column separately. Defaults to FALSE. 
#' @param node.features string or vector of strings - denotes the column name(s) in the VDJ dataframe from which the node features should be extracted (which can, for example, be used for plotting of lineage trees later on).
#' @param construction.method string - denotes the approach and algorithm that will be used to convert the distance matrix or multiple sequence alignment into a lineage tree. There are two approaches two construct a lineage tree: a tree can be constructed from a network/graph (phylo.network) or from a phylogenetic tree (phylo.tree). There are three algorithm options that take a pairwise distance matrix as input: 'phylo.network.default', 'phylo.network.mst', and 'phylo.tree.nj'. There are two algorithm options that take a multiple sequence alignment as input: 'phylo.tree.ml', and 'phylo.tree.mp'. Defaults to 'phylo.network.default' (mst-like algorithm).
#' 'phylo.network.default': mst-like tree evolutionary network algorithm in which the germline node is positioned at the top of the tree, and nodes with the minimum distance to any existing node in the tree are linked iteratively.
#' 'phylo.network.mst'    : minimum spanning tree (MST) algorithm from 'ape::mst()' constructs networks with the minimum sum of edge lengths/weights, which involves iteratively adding edges to the network in ascending order of edge weights, while ensuring that no cycles are formed, after which the network is reorganized into a germline-rooted lineage tree.
#' 'phylo.tree.nj'        : neighbor-joining (NJ) algorithm from 'ape::nj()' constructs phylogenetic trees by joining pairs of nodes with the minimum distance, creating a bifurcating tree consisting of internal nodes (representing unrecovered sequences) and terminal nodes (representing the recovered sequences). 
#'  NB: More information about the use of the 'ape' package is provided in the book 'Analysis of Phylogenetic and Evolution with R' by Emmanual Paradis (https://link.springer.com/book/10.1007/978-1-4614-1743-9).
#' 'phylo.tree.mp'        : maximum-parsimony (MP) algorithm from 'phangorn::pratchet()' constructs phylogenetic trees by minimizing the total number of edits required to explain the observed differences among sequences.
#' 'phylo.tree.ml'        : maximum-likelihood (ML) algorithm from 'phangorn::pml_bb()' constructs phylogenetic trees by estimating the tree topology and branch lengths that maximize the likelihood of observing the given sequence data under a specified evolutionary model. 
#'  NB: More information about estimating phylogenetic trees with the 'phangorn' package is provided on the Github page of Klaus Schliep: https://klausvigo.github.io/phangorn/articles/Trees.html. 
#' 'phylo.tree.IgPhyML'   : no trees/network are inferred, but trees are directly imported from
#'  NB: More information about the IgPhyML tool and the phylogenetic substitution model for antibody lineages is provded by Hoehn et al. (2017): https://doi.org/10.1534/genetics.116.196303.
#' @param IgPhyML.output.file string - specifies the path to the IgPhyML output file, from which the trees will be imported (if 'construction.method' is set to 'phylo.tree.IgPhyML').
#' @param string.dist.metric string - denotes the metric that will be calculated with the 'stringdist::stringdistmatrix()' function to measure (string) distance between sequences. Options: 'lv', 'dl', 'osa', 'hamming', 'lcs', 'qgram', 'cosine', 'jaccard', and 'jw'.  Defaults to 'lv' (Levenshtein distance / edit distance).
#' 'lv'       : Levensthein distance (also known as edit distance) equals to the minimum number of single-element edits (insertions, deletions, or substitutions) required to transformer one string into another.
#' 'dl'       : Damerau-Levenshtein distance is similar to the Levenshtein distance, but also allows transpositions of adjacent elements as a single-edit operation.
#' 'osa'      : Optimal String Alignment distance is similar to the Damerau-Levensthein distance, but does not allow to apply multiple transformations on a same substring.
#' 'hamming'  : Hamming distance equals to the number of positions at which the corresponding elements differ between two strings (applicable only to strings of equal length).
#' 'lcs'      : Longest Common Subsequence distance is similar to the Levenshtein distance, but only allowing insertions and deletions as single-edit operations.
#' 'qgram'    : Q-gram distance equal to the number of distinct q-grams that appear in either string but not both, whereby q-grams are all possible substrings of length q in both strings (q defaults to 1).
#' 'cosine'   : cosine distance equals to 1 - cosine similarity (the strings are converted into vectors containing the frequency of all single elements (A and B), whereby the cosine similarity (Sc) equals to the dot product of these vectors divided by the product of the magnitude of these vectors, which can be written in a formula as Sc(A, B) = A . B / (||A|| x ||B||)).
#' 'jaccard'  : Jaccard distance equals to 1 - Jaccard index (the strings are converted into sets of single elements (A and B), whereby the Jaccard index (J) equals to the size of the intersection of the two sets divided by the size of the union of the sets, which can be written in a formula as J(A, B) = |A ∩ B| / |A ∪ B|).
#' 'jw'       : Jaro-Winkler distance equals to 1 - Jaro-Winkler similarity (Jaro-Winkler similary is calculated with the following formulas: Sw = Sj + P * L * (1-Sj) in which Sw is the Jaro-Winkler similary, Sj is the Jaro similarity, P is the scaling factor (defaults to 0), and L is the length of the matching prefix; and Sj = 1/3 * (m/|s1| + m/|s2| + (m-t)/m) in which Sj is the Jaro similarity, m is the number of matching elements, |s1| and |s2|are the lengths of the strings, and t is the number of transpositions).
#' @param dna.model string or vector of strings - specifies the DNA model(s) to be used during distance calculation or maximum likelihood tree inference.
#' When using one of the distance-based construction methods ('phylo.network.default', 'phylo.network.mst', or 'phylo.tree.nj'), an evolutionary model can be used to compute a pairwise distance matrix from DNA sequences using the 'ape::dist.dna()' function. 
#' Available DNA models: 'raw', 'N', 'TS', 'TV', 'JC69', 'K80', 'F81', 'K81', 'F84', 'BH87', 'T92', 'TN93', 'GG95', 'logdet', 'paralin', 'indel', and 'indelblock'.
#' When using the 'phylo.tree.ml' construction method, models are compared with each other with the 'phangorn::modelTest()' function, of which the output will be used as input for the  'phangorn::pml_bb()' function to infer the maximum likelihood tree. The best model according to the BIC (Bayesian information criterion) will be used to infer the tree. Defaults to "all" (when nucleotide sequences are found in the specified 'sequence.columns' and the 'germline.columns').
#' Available DNA models: 'JC', 'F81', 'K80', 'HKY', 'TrNe', 'TrN', 'TPM1', 'K81', 'TPM1u', 'TPM2', 'TPM2u', 'TPM3', 'TPM3u', 'TIM1e', 'TIM1', 'TIM2e', 'TIM2', 'TIM3e', 'TIM3', 'TVMe', 'TVM', 'SYM', and 'GTR'.
#' @param aa.model string or vector of strings - specifies the AA model(s) to be used during distance calculation or maximum likelihood tree inference.
#' When using one of the distance-based construction methods ('phylo.network.default', 'phylo.network.mst', or 'phylo.tree.nj'), an evolutionary model can be used to compute a pairwise distance matrix from AA sequences using the 'phangorn::dist.ml()' function. 
#' Available AA models: 'WAG', 'JTT', 'LG', 'Dayhoff', 'cpREV', 'mtmam', 'mtArt', 'MtZoa', 'mtREV24', 'VT', 'RtREV', 'HIVw', 'HIVb', 'FLU', 'Blosum62', 'Dayhoff_DCMut', and 'JTT_DCMut'.
#' When using the 'phylo.tree.ml' construction method, models are compared with each other with the 'phangorn::modelTest()' function, of which the output will be used as input for the  'phangorn::pml_bb' function to infer the maximum likelihood tree. The best model according to the BIC (Bayesian information criterion) will be used to infer the tree. Defaults to "all" (when protein sequences are found in the specified 'sequence.columns' and the 'germline.columns').
#' Available AA models: 'WAG', 'JTT', 'LG', 'Dayhoff', 'cpREV', 'mtmam', 'mtArt', 'MtZoa', 'mtREV24', 'VT', 'RtREV', 'HIVw', 'HIVb', 'FLU', 'Blosum62', 'Dayhoff_DCMut', and 'JTT-DCMut'.
#' @param codon.model string or vector of strings - specifies the codon substitution models to compare with each other with the 'phangorn::codonTest()' function (only possible when the 'construction.method' paramter is set to 'phylo.tree.ml' and when colums with DNA sequences are selected). The best model according to the BIC (Bayesian information criterion) will be used to infer the tree, and this tree will replace the tree inferred with the best model of the model specified in the 'dna.models' parameter. Defaults to NA.
#' Available codon models: 'M0'.
#' @param resolve.ties string or vector of strings - denotes the way ties are handled during the conversion of the distance matrix into lineage trees by the 'phylo.network.tree' algorithm (in the event where an unlinked node, that is to be linked to the tree next, shares identical distances with multiple previously linked nodes in the lineage tree). Options: 'min.expansion', 'max.expansion', 'min.germline.dist', 'max.germline.dist', 'min.germline.edges', 'max.germline.edges', and 'random'. If a vector is provided, ties will be resolved in a hierarchical manner. Defaults to 'c("max.expansion", "close.germline.dist", "close.germline.edges", "random")'.
#' 'min.expansion'        : the node(s) having the smallest size is/are selected.
#' 'max.expansion'        : the node(s) having the biggest size is/are selected.
#' 'min.germline.dist'    : the node(s) having the smallets string distance to the germline node is/are selected.
#' 'max.germline.dist'    : the node(s) having the biggest string distance to the germline node is/are selected.
#' 'min.germline.edges'   : the node(s) having the lowest possible number of edges to the germline node is/are selected.
#' 'max.germline.edges    : the node(s) having the highest possible number of edges to the germline node is/are selected.
#' 'random'               : a random node is selected.
#' @param remove.internal.nodes string - denotes if and how internal nodes should be removed when the 'construction.method' is set to 'phylo.tree.nj', 'phylo.tree.mp', 'phylo.tree.ml' or 'phylo.tree.IgPhyML'. Options: 'zero.length.edges.only', 'connect.to.parent', 'minimum.length', and 'minimum.cost'. Defaults to 'minimum.cost', when 'construction.method' is set to 'phylo.tree.nj'. Defautls to 'connect.to.parent', when 'construction.method' is set to 'phylo.tree.mp', 'phylo.tree.ml', or 'phylo.tree.IgPhyML'. 
#' 'zero.length.edges.only' : only internal nodes with a distance of zero to a terminal node are removed by replacing it with the terminal node.
#' 'connect.to.parent'      : connects all terminal nodes to the first parental sequence-recovered node upper in the tree, resulting in a germline-directed tree.
#' 'minimum.length'         : iteratively replaces internal nodes with terminal nodes that are linked by an edge that has the minimum length.
#' 'minimum.cost'           : iteratively replaces internal nodes with terminal nodes which results in the minimum increase in the sum of all edges (this increase is referred to as the 'cost').
#' @param include string or vector of strings - specifies the objects to be included in the output object for each clonotype (if created). Options: 'nodes', 'dist.matrices', 'msa', 'phylo', 'igraph', 'igraph.with.inner.nodes', 'metrics', or 'all' to select all objects. Defaults to 'all'.
#' 'nodes'                    : nested list wherein for each node, all information is stored (sequences, barcodes, selected column in 'node.features').
#' 'dist'                     : pairwise string distance matrices calculated using the specified 'string.dist.metric', one for each column selected in 'sequence.columns', or only one if 'concatenate_sequences' is set to TRUE.
#' 'msa'                      : multiple sequence alignments, one for each column selected in 'sequence.columns', or only one if 'concatenate_sequences' is set to TRUE.
#' 'phylo'                    : object of class 'phylo' that is created when 'construction.method' is set to 'phylo.tree.nj', 'phylo.tree.mp', or 'phylo.tree.ml', and when the clonotype contains at least three sequences.
#' 'igraph'                   : object of class 'igraph' that represent the B cell lineage tree, which is used for plotting by the 'plot_lineage_tree()' function.
#' 'igraph.with.inner.nodes'  : object of class 'igraph' that represent the B cell lineage tree before the removal of internal nodes (if 'remove.internal.nodes' is set to 'connect.to.parent' or 'all').
#' 'edges'                    : dataframe with the three columns 'upper.node', 'lower.node', and 'edge.length', whereby each row in the dataframe represent an edge in the lineage tree.  
#' 'edges.with.inner.nodes'   : dataframe with the three columns 'upper.node', 'lower.node', and 'edge.length', whereby each row in the dataframe represent an edge in the lineage tree.  
#' 'metrics'                  : list of tree metrics that can only be calculated during the construction of the lineage tree, which includes a 'tie.resolving' matrix, indicating which options were used to handle ties (when 'construction.method' is set to 'phylo.network.default'), and a 'model' string, indicating which model was used to infer the maximum likelihood tree (if 'construction.method' is set to 'phylo.tree.ml').
#' @param parallel bool - if TRUE, the per-clone network inference is executed in parallel (parallelized across samples). Defaults to FALSE.
#' @param num.cores integer - number of cores to be used when parallel = TRUE. Defaults to all available cores - 1 or the number of samples in the VDJ dataframe (depending which number is smaller).
#' @return An object of class 'AntibodyForests', structured as a nested list where each outer list represents a sample, and each inner list represents a clonotype. Each clonotype list contains the output objects specified in the 'include' parameter. For example, AntibodyForests[[1]][[2]] contains the list of output objects for the first sample and third clonotype (which would be equivalent to something like AntibodyForests$S1$clonotype3).
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, 
#'                 sequence.columns = c("VDJ_sequence_aa_trimmed","VJ_sequence_aa_trimmed"),
#'                 germline.columns = c("VDJ_germline_aa_trimmed","VJ_germline_aa_trimmed"), 
#'                 concatenate.sequences = TRUE,
#'                 node.features = c("VDJ_vgene", "VDJ_dgene", "VDJ_jgene", "isotype"),
#'                 string.dist.metric = "lv",
#'                 construction.method = "phylo.network.default",
#'                 resolve.ties = c("max.expansion", "close.germline.dist", "close.germline.edges", "random"),
#'                 include = c("nodes", "dist", "igraph", "edges", "metrics"),
#'                 parallel = TRUE)
#'}

AntibodyForests <- function(VDJ,
                            sequence.columns,
                            germline.columns,
                            concatenate.sequences,
                            node.features,
                            string.dist.metric,
                            dna.model,
                            aa.model,
                            codon.model,
                            construction.method,
                            IgPhyML.output.file,
                            resolve.ties,
                            remove.internal.nodes,
                            include,
                            parallel,
                            num.cores){
  
  # Store original working directory
  original_working_directory <- getwd()
  
  # If no 'VDJ' dataframe is provided, a message is returned and execution is stopped
  if(missing(VDJ)){stop("ERROR: Please provide a VDJ dataframe as obtained from the 'minimal_VDJ()' function in Platypus, or a similar dataframe containing the specified sequence and germline colums.")}
  
  # If the 'sequence.columns' parameter is not specified, and no IgPhyML output file is provided, the 'VDJ_sequence_nt_trimmed' and 'VJ_sequence_nt_trimmed' columns are selected, and a message is returned
  if(missing(sequence.columns) && missing(IgPhyML.output.file)){sequence.columns <- c("VDJ_sequence_nt_trimmed", "VJ_sequence_nt_trimmed"); message("WARNING: No sequence columns are specified. Defaults to 'VDJ_sequence_nt_trimmed' and 'VJ_sequence_nt_trimmed'.\n")}
  # If the 'sequence.columns' parameter is specified, while an IgPhyML output file is provided, a message is returned
  if(!missing(sequence.columns) && !missing(IgPhyML.output.file)){message("NB: The selected sequence column(s) should exactly match the input provided to the IgPhyML tool.")}
  # If the 'sequence.columns' parameter is not specified, while an IgPhyML output file is provided, the 'VDJ_sequence_nt_trimmed' columns is selected, and a message is returned 
  if(missing(sequence.columns) && !missing(IgPhyML.output.file)){sequence.columns <- c("VDJ_sequence_nt_trimmed"); message("WARNING: No sequence columns are specified. As an IgPhyML output file is provided, only the 'VDJ_sequence_nt_trimmed' column is selected.\n")}
  # If the columns specified in the 'sequence.columns' parameter are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(sequence.columns %in% colnames(VDJ))){stop("ERROR: Please provide valid 'sequence.columns' from the input VDJ dataframe.")}
  # If the 'germline.columns' parameter is not specified, and no IgPhyML output file is provided, the 'VDJ_germline_nt_trimmed' and 'VJ_germline_nt_trimmed' columns are selected, and a message is returned
  if(missing(germline.columns) && missing(IgPhyML.output.file)){germline.columns <- c("VDJ_germline_nt_trimmed", "VJ_germline_nt_trimmed"); message("WARNING: No germline columns are specified. Defaults to 'VDJ_germline_nt_trimmed' and 'VJ_germline_nt_trimmed'.\n")}
  # If the 'germline.columns' parameter is specified, while an IgPhyML output file is provided, a message is returned
  if(!missing(germline.columns) && !missing(IgPhyML.output.file)){message("NB: The selected germline column(s) should exactly match the input provided to the IgPhyML tool.")}
  # If the 'germline.columns' parameter is not specified, while an IgPhyML output file is provided, the 'VDJ_sequence_nt_trimmed' columns is selected, and a message is returned 
  if(missing(germline.columns) && !missing(IgPhyML.output.file)){germline.columns <- c("VDJ_germline_nt_trimmed"); message("WARNING: No germline columns are specified. As an IgPhyML output file is provided, only the 'VDJ_germline_nt_trimmed' column is selected.\n")}
  # If the columns specified in the 'germline.columns' parameter are not all present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(germline.columns %in% colnames(VDJ))){stop("ERROR: Please provide valid 'germline.columns' from the input VDJ dataframe.")}
  # If the columns in 'germline.columns' do not correspond to the columns in 'sequence.columns', a warning is returned
  if(!all(sequence.columns == gsub("germline", "sequence", germline.columns))){message("WARNING: Please double check whether the columns in 'germline.columns' correspond to the columns in 'sequence.columns'.\n")}
  
  # Define the 'sequence_type' (DNA sequence or protein sequence) based on the sequences in the specified 'sequence.columns' and the 'germline.columns'
  if(all(sapply(c(sequence.columns, germline.columns), function(x) all(sapply(1:nrow(VDJ), function(y) grepl(pattern = "^[ACTG-]+$", VDJ[y, x])))))){
    sequence_type <- "DNA"
  } else if(all(sapply(c(sequence.columns, germline.columns), function(x) all(sapply(1:nrow(VDJ), function(y) (grepl(pattern = "^[ARNDCQEGHILKMFPSTWYV*-]+$", VDJ[y, x]) && !grepl(pattern = "^[ACTG]+$", VDJ[y, x]))))))){
    sequence_type <- "AA"
  } else{
    # If the 'sequence_type' cannot be defined, a message is returned and execution is stopped
    if(!exists("sequence_type")){stop("ERROR: Please provide an input dataframe that contains valid DNA or protein sequences in the specified sequence and germline columns. NB: the sequences in the 'sequence.columns' and 'germline.columns' should be of the same type.")}
  }
  
  # If the 'concatenate.sequences' parameter is missing, it is set to FALSE
  if(missing(concatenate.sequences)){concatenate.sequences <- FALSE}
  
  # If no columns are specified in 'node.features', it defaults to an empty vector
  if(missing(node.features)){node.features <- c()}
  # If 'node.features' contains columns that are not present in the 'VDJ' dataframe, a message is returned and execution is stopped
  if(!all(node.features %in% colnames(VDJ))){stop("ERROR: Not all columns in 'node.features' could be found in the input VDJ dataframe.")}
  
  # If the 'construction.method' parameter is missing, and the 'IgPhyML.output.file' parameter is not specified, the 'construction.method' parameter is set to 'phylo.network.default'
  if(missing(construction.method) && missing(IgPhyML.output.file)){construction.method <- "phylo.network.default"}
  # If an IgPhyML output file is provided, the 'construction.method' parameter is set to 'phylo.tree.IgPhyML'
  if(!missing(IgPhyML.output.file) && missing(construction.method)){construction.method <- "phylo.tree.IgPhyML"}
  # If the 'construction.method' parameter is set to an unknown method, a message is returned and execution is stopped
  if(!construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj", "phylo.tree.mp", "phylo.tree.ml", "phylo.tree.IgPhyML")){stop("ERROR: The specified tree construction method is not recognized. Please choose from the following options: 'phylo.network.default', 'phylo.network.mst', 'phylo.tree.nj', 'phylo.tree.mp', or 'phylo.tree.ml'.")}
  # If the 'construction.method' is defined, retrieve the 'approach' and 'algorithm' from it
  approach <- unlist(base::strsplit(construction.method, split = "\\."))[2]; algorithm <- unlist(base::strsplit(construction.method, split = "\\."))[3]
  
  # If the 'construction.method' parameter is set to 'phylo.tree.IgPhyML', but no path to an IgPhyML output file is specified, a message is returned and execution is stopped
  if(construction.method == "phylo.tree.IgPhyML" && missing(IgPhyML.output.file)){stop("ERROR: The 'construction.method' parameter is set to 'phylo.tree.IgPhyML', but no IgPhyML output is provided. Please specify the path to the IgPhyML output file!!!")}
  # Check whether the specified path to the IgPhyML output file exists
  if(construction.method == "phylo.tree.IgPhyML"){if(!file.exists(IgPhyML.output.file)){stop("ERROR: The IgPhyML output file does not exist. Please check the specified path!")}}
  # If an IgPhyML output file is provided, while the 'construction.method' parameter is set to a different method than 'phylo.tree.IgPhyML', a message is returned and execution is stopped
  if(!missing(IgPhyML.output.file) && construction.method != "phylo.tree.IgPhyML"){stop("ERROR: If the path to an IgPhyML output file is specified, the 'construction.method' parameter can only be set to 'phylo.tree.IgPhyML' or be left empty.")}
  
  # If the 'construction.method' parameter is set to 'phylo.network.default', but non-relevant parameters are specified ('codon.model', and/or 'remove.internal.nodes'), a warning message is returned
  if(construction.method == "phylo.network.default" && (!missing(codon.model) | !missing(remove.internal.nodes))){stop("ERROR: If the 'construction.method' parameter is set to 'phylo.network.default', the 'codon.model' and 'remove.internal.nodes' parameters cannot be specified.")}
  # If the 'construction.method' parameter is set to 'phylo.network.mst', but non-relevant parameters are specified ('codon.model', 'resolve.ties', and/or 'remove.internal.nodes'), a warning message is returned
  if(construction.method == "phylo.network.mst" && (!missing(codon.model) | !missing(resolve.ties) | !missing(remove.internal.nodes))){stop("ERROR: If the 'construction.method' parameter is set to 'phylo.network.mst', the 'codon.model', 'resolve.ties', and 'remove.internal.nodes' parameters cannot be specified.")}
  # If the 'construction.method' parameter is set to 'phylo.tree.nj', but non-relevant parameters are specified ('codon.model', and/or 'resolve.ties'), a warning message is returned
  if(construction.method == "phylo.tree.nj" && (!missing(codon.model) | !missing(resolve.ties))){stop("ERROR: If the 'construction.method' parameter is set to 'phylo.tree.nj', the 'codon.model', and 'resolve.ties' cannot be specified.")}
  # If the 'construction.method' parameter is to a distance-based method, but both the 'string.dist.metric' and 'dna.model' or 'aa.model' are specified, a message is returned and execution is stopped
  if(construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.network.nj") && !missing(string.dist.metric) && !(missing(dna.model) | missing(aa.model))){stop("ERROR: Both a string distance metric and a model are specified, please choose one.")}
  # If the 'construction.method' parameter is set to 'phylo.tree.mp', but non-relevant parameters are specified ('string.dist.metric', 'dna.model', 'aa.model', 'codon.model', and/or 'resolve.ties'), a warning message is returned
  if(construction.method == "phylo.tree.mp" && (!missing(string.dist.metric) | !missing(dna.model) | !missing(aa.model) | !missing(codon.model) | !missing(resolve.ties))){stop("ERROR: If the 'construction.method' parameter is set to 'phylo.tree.mp', the 'string.dist.metric, 'dna.model', 'aa.model', 'codon.model', and 'resolve.ties' parameters cannot be specified.")}
  # If the 'construction.method' parameter is set to 'phylo.tree.mp', but non-relevant parameters are specified ('string.dist.metric', and/or 'resolve.ties'), a warning message is returned
  if(construction.method == "phylo.tree.ml" && (!missing(string.dist.metric) | !missing(resolve.ties))){stop("ERROR: If the 'construction.method' parameter is set to 'phylo.tree.ml', the 'string.dist.metric', and 'resolve.ties' parameters cannot be speciied.")}
  # If the 'construction.method' parameter is set to 'phylo.tree.IgPhyML' and if an IgPhyML output file is provided, but other non-relevant parameters are specified ('construction.method', 'string.dist.metric', 'dna.model', 'aa.model', 'codon.model', and/or 'resolve.ties'), a warning message is returned
  if(construction.method == "phylo.tree.IgPhyML" && (!missing(string.dist.metric) | !missing(dna.model) | !missing(aa.model) | !missing(codon.model) | !missing(resolve.ties))){warning("ERROR: If a IgPhyML output file is provided, the 'string.dist.metric', 'dna.model', 'aa.model', 'codon.model', and 'resolve.ties' parameters cannot be specified.")}
  
  # If the 'construction.method' is set to a distance-based construction method, but the 'string.dist.metric', 'dna.model', and  'aa.model' parameters are all missing, the 'string.dist.metric' is set to 'lv', which means that a pairwise Levensthein distance matrix will be computed with the 'stringdist::stringdistmatrix()' function 
  if(construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && missing(string.dist.metric) && missing(dna.model) && missing(aa.model)){string.dist.metric <- "lv"}
  # If the 'construction.method' is set to a distance-based construction method, but both the 'string.dist.metric' parameter and 'dna.model' or 'aa.model' parameter are specified, a message is returned and execution is stopped
  if(construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && !missing(string.dist.metric) && (!missing(dna.model) | !missing(aa.model))){stop("ERROR: Please provide a string distance metric OR an evolutionary model to compute the pairwise distance matrices, but not both!")}
  # If the 'construction.method' is set to 'phylo.tree.ml' and the 'sequence_type' is 'DNA', but the 'dna.model' parameter is missing, it is set to 'all'
  if(construction.method == "phylo.tree.ml" && sequence_type == "DNA" && missing(dna.model)){dna.model <- "all"; aa.model <- NA}
  # If the 'construction.method' is set to 'phylo.tree.ml' and the 'sequence_type' is 'AA', but the 'aa.model' parameter is missing, it is set to 'all'
  if(construction.method == "phylo.tree.ml" && sequence_type == "AA" && missing(aa.model)){aa.model <- "all"; dna.model <- NA}
  
  # If the 'aa.model' parameter is specified, while DNA sequences are found, or if the 'dna.model' parameter is specified, while protein sequence are found, a message is returned and execution is stopped
  if(sequence_type == "DNA" && !missing(aa.model)){if(!is.na(aa.model)){stop("ERROR: The sequences in the specified sequence and germline columns are DNA sequences, whereas an amino acid model is specified.")}}
  if(sequence_type == "AA" && !missing(dna.model)){if(!is.na(dna.model)){stop("ERROR: The sequences in the specified sequence and germline columns are protein sequences, whereas a DNA model is specified.")}}
  # If the 'construction.method' parameter is set to 'phylo.tree.ml', and the 'dna.model' or 'aa.model' parameter is set to all, return a message that the tree inference may take up to several hours
  if(construction.method == "phylo.tree.ml" && (dna.model == "all" | aa.model == "all")){message("WARNING: Comparing all available models and using the best fitting model for parameter optimization and phylogenetic tree inference may take up to several hours, when a large set of clonotypes is given as input. Please be patient!\n")}
  # If the 'dna.model' or 'aa.model' parameter is set to 'all', it set to all available models within the 'phangorn' package for the current 'sequence_type'
  if(construction.method == "phylo.tree.ml" && sequence_type == "DNA" && is.na(aa.model)){if(dna.model ==  "all"){dna.model <- c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR")}}
  if(construction.method == "phylo.tree.ml" && sequence_type == "AA" && is.na(dna.model)){if(aa.model == "all"){aa.model <- c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT", "RtREV", "HIVw", "HIVb", "FLU", "Blosum62", "Dayhoff_DCMut", "JTT-DCMut")}}
  
  # If the 'construction.method' is set to 'phylo.network.default', but 'resolve.ties' is missing, the 'max.expansion', 'min.germline.dist', 'min.germline.edges', and 'random' options are used to handle (all) ties
  if(construction.method == "phylo.network.default" && missing(resolve.ties)){resolve.ties <- c("max.expansion", "min.germline.dist", "min.germline.edges", "random")}
  
  # If the 'construction.method' is set to 'phylo.tree.nj', but the 'remove.internal.nodes' parameter is missing, it is set to 'minimum.cost'
  if(construction.method == "phylo.tree.nj" && missing(remove.internal.nodes)){remove.internal.nodes <- "minimum.cost"}
  # If the 'construction.method' is set to 'phylo.tree.mp', 'phylo.tree.ml', or 'phylo.tree.IgPhyML', but the 'remove.internal.nodes' parameter is missing, it is set to 'connect.to.parent'
  if(construction.method %in% c("phylo.tree.mp", "phylo.tree.ml", "phylo.tree.IgPhyML") && missing(remove.internal.nodes)){remove.internal.nodes <- "connect.to.parent"}
  
  # If the 'string.dist.metric', 'dna.model', 'aa.model', 'codon.model', 'resolve.ties', and 'remove.internal.nodes' parameter are still not defined, these parameters will not be used, and are set to NA
  if(missing(string.dist.metric)){string.dist.metric <- NA}
  if(missing(dna.model)){dna.model <- NA}
  if(missing(aa.model)){aa.model <- NA}
  if(missing(codon.model)){codon.model <- NA}
  if(missing(resolve.ties)){resolve.ties <- NA}
  if(missing(remove.internal.nodes)){remove.internal.nodes <- NA}
  
  # If the 'string.dist.metric' is set to an unknown string distance metrics, a message is returned and execution is stopped.
  if(!is.na(string.dist.metric)){if(!string.dist.metric %in% c("lv", "dl", "osa", "hamming", "lcs", "qgram", "cosine", "jaccard", "jw")){stop("ERROR: The specified string distance metric is not recognized. Please choose from the following options: 'lv', 'dl', 'osa', 'hamming', 'lcs', 'qgram', 'cosine', 'jaccard', or 'jw'.")}}
  # If the 'resolve.ties' parameter contains options that are not recognized, a message is returned and execution is stopped
  if(!all(is.na(resolve.ties))){if(!all(resolve.ties %in% c("min.expansion", "max.expansion", "min.germline.dist", "max.germline.dist", "min.germline.edges", "max.germline.edges", "random"))){stop("ERROR: Not all options specified in 'resolve.ties' are recognized. Please choose from the following options: 'min.expansion', 'max.expansion', 'min.germline.dist', 'max.germline.dist', 'min.germline.edges', 'max.germline.edges', and 'random'.")}}
  # If the 'remove.internal.nodes' parameter is set to an algorithm that is not recognized, a message is returned and execution is stopped
  if(!is.na(remove.internal.nodes)){if(!remove.internal.nodes %in% c("zero.length.edges.only", "connect.to.parent", "minimum.length", "minimum.cost")){stop("ERROR: The specified 'remove.internal.nodes' algorithm' is not recognized. Please choose from the following options: 'zero.length.edges.only', 'connect.to.parent', 'minimum.length', or 'minimum.weight'.")}}
  # If the 'dna.model' or 'aa.model' parameter contains models that are not available for the current 'sequence_type' or for the current 'construction.method', a message is returned and execution is stopped
  if(sequence_type == "DNA" && construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj")){if(!all(is.na(dna.model))){if(!all(dna.model %in% c("raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock"))){stop("ERROR: Not all models specified in the 'dna.model' parameter are recognized as available evolutionary DNA models. Please choose from the following options: 'raw', 'N', 'TS', 'TV', 'JC69', 'K80', 'F81', 'K81', 'F84', 'BH87', 'T92', 'TN93', 'GG95', 'logdet', 'paralin', 'indel', or 'indelblock'.")}}}
  if(sequence_type == "AA" && construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj")){if(!all(is.na(dna.model))){if(!all(aa.model %in% c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT", "RtREV", "HIVw", "HIVb", "FLU", "Blosum62", "Dayhoff_DCMut", "JTT_DCMut"))){stop("ERROR: Not all models specified in the 'aa.model' parameter are recognized as available evolutionary AA models. Please choose from the following options: 'WAG', 'JTT', 'LG', 'Dayhoff', 'cpREV', 'mtmam', 'mtArt', 'MtZoa', 'mtREV24', 'VT', 'RtREV', 'HIVw', 'HIVb', 'FLU', 'Blosum62', 'Dayhoff_DCMut', or 'JTT_DCMut'.")}}}
  if(sequence_type == "DNA" && construction.method == "phylo.tree.ml"){if(!all(is.na(dna.model))){if(!all(dna.model %in% c("JC", "F81", "K80", "HKY", "TrNe", "TrN", "TPM1", "K81", "TPM1u", "TPM2", "TPM2u", "TPM3", "TPM3u", "TIM1e", "TIM1", "TIM2e", "TIM2", "TIM3e", "TIM3", "TVMe", "TVM", "SYM", "GTR"))){stop("ERROR: Not all models specified in the 'dna.model' parameter are recognized as available DNA models for the maximum likelihood construction method. Please choose from the following options: 'JC', 'F81', 'K80', 'HKY', 'TrNe', 'TrN', 'TPM1', 'K81', 'TPM1u', 'TPM2', 'TPM2u', 'TPM3', 'TPM3u', 'TIM1e', 'TIM1', 'TIM2e', 'TIM2', 'TIM3e', 'TIM3', 'TVMe', 'TVM', 'SYM', and/or 'GTR'.")}}}
  if(sequence_type == "AA" && construction.method == "phylo.tree.ml"){if(!all(is.na(aa.model))){if(!all(aa.model %in% c("WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT", "RtREV", "HIVw", "HIVb", "FLU", "Blosum62", "Dayhoff_DCMut", "JTT-DCMut"))){stop("ERROR: Not all models specified in the 'aa.model' parameter are recognized as available AA models for the maximum likilihood construction method. Please choose from the following options: 'WAG', 'JTT', 'LG', 'Dayhoff', 'cpREV', 'mtmam', 'mtArt', 'MtZoa', 'mtREV24', 'VT', 'RtREV', 'HIVw', 'HIVb', 'FLU', 'Blosum62', 'Dayhoff_DCMut', and/or 'JTT-DCMut'.")}}}
  # If the 'codon.model' parameter is missing, it set to NA
  if(missing(codon.model)){codon.model <- NA}
  # If the 'codon.model' parameter contains models that are not available, a message is returned and execution is stopped
  if(!all(is.na(codon.model))){if(!all(codon.model %in% c("M0"))){stop("ERROR: Not all models specified in the 'codon.model' parameter are recognized as available codon models. Please choose from the following options: 'M0', 'M1a', and 'M2a'.")}}
  # If the 'codon.model' parameter is specified, while protein sequences are found in the selected 'sequence.columns' and 'germline.columns', a message is returned and execution is stopped
  if(!all(is.na(codon.model)) && sequence_type == "AA"){stop("ERROR: Currently, a codon model can only be used for tree inference, if the sequences in the selected 'sequence.columns' and 'germline.columns' are DNA sequences.")}
  
  # If the 'include' parameter is missing, it is set to 'all'
  if(missing(include)){include <- "all"}
  # If the 'include' parameter is set to 'all', all the object that are created with the specified 'construction.method' will be included in the AntibodyForests object
  if(all(include == "all") && construction.method == "phylo.network.default"){include <- c("nodes", "dist", "igraph", "edges", "metrics")}
  if(all(include == "all") && construction.method == "phylo.network.mst"){include <- c("nodes", "dist", "igraph", "edges")}
  if(all(include == "all") && construction.method == "phylo.tree.nj"){include <- c("nodes", "dist", "phylo", "igraph", "igraph.with.inner.nodes", "edges", "edges.with.inner.nodes")}
  if(all(include == "all") && construction.method == "phylo.tree.mp"){include <- c("nodes", "msa", "phylo", "igraph", "igraph.with.inner.nodes", "edges", "edges.with.inner.nodes")}
  if(all(include == "all") && construction.method == "phylo.tree.ml"){include <- c("nodes", "msa", "phylo", "igraph", "igraph.with.inner.nodes", "edges", "edges.with.inner.nodes", "metrics")}
  if(all(include == "all") && construction.method == "phylo.tree.IgPhyML"){include <- c("nodes", "igraph", "igraph.with.inner.nodes", "edges", "edges.with.inner.nodes")}
  # If the 'include' parameter contains output objects that are not created/recognized, a message is returned and execution is stopped
  if(!all(include %in% c("nodes", "dist", "msa", "phylo", "igraph", "igraph.with.inner.nodes", "edges", "edges.with.inner.nodes", "metrics"))){stop("ERROR: Not all specified output objects are recognized. Please choose from the following options: 'nodes', 'dist', 'msa', 'phylo', 'igraph', 'igraph.with.inner.nodes', 'edges', 'edges.with.inner.nodes'', and 'metrics'.")}
  
  # If the 'parallel' parameter is missing, it is set to FALSE
  if(missing(parallel)){parallel <- FALSE}
  # If 'parallel' is set to TRUE but 'num.cores' is not specified, the number of cores is set to all available cores - 1
  if(parallel == TRUE && missing(num.cores)){num.cores <- parallel::detectCores() -1}
  
  
  build_lineage_tree <- function(clone,
                                 dist.matrix,
                                 msa,
                                 sequence.type,
                                 approach,
                                 algorithm,
                                 IgPhyML.trees,
                                 dna.model,
                                 aa.model,
                                 codon.model,
                                 resolve.ties,
                                 remove.internal.nodes,
                                 nodes.list){
    
    # Builds a B cell lineage tree using a pairwise distance matrix, a multiple sequence alignment (msa), or an IgPhyML output file as input
    # Arguments:
    # - clone: string specifying the sample and clonotype of the distance matrix or msa in the format 'S1_clonotype1'
    # - dist.matrix: distance matrix containing pairwise (string) distances between all possible pairs of sequences/nodes
    # - msa: multiple sequence alignment in the phyDat format
    # - sequence.type: type of the sequences aligned in the 'msa'
    # - approach: string specifying the approach for constructing the lineage tree ('network' for a network/graph-derived lineage tree or 'tree' for a phylogenetic tree-derived lineage tree)
    # - algorithm: string denoting the exact algorithm used to convert the input (distance matrix or msa) into a lineage tree
    # - IgPhyML.trees: list of objects of class 'igraph', as obtained from the IgPhyML output file using the 'alakazam::readIgphyml()' function
    # - dna.model: string or vector of strings specifying the nucleotide substitution models that are compared when the 'sequence.type' is set to 'DNA', the 'approach' is set to 'tree', and the 'algorithm' is set 'ml' (the best model is used to infer the maximum likelihood phylogenetic tree)
    # - aa.model: string or vector of strings specifying the amino acid substitution models that are compared when the 'sequence.type' is set to 'AA', the approach' is set to 'tree', and the 'algorithm' is set 'ml' (the best model is used to infer the maximum likelihood phylogenetic tree)
    # - codon.model: string or vector strings specifying the codon substitution models that are compared when the 'sequence.type' is set to 'DNA', the 'approach' is set to 'tree', and the 'algorithm' is set 'ml' (the best model is used to infer the maximum likelihood phylogenetic tree)
    # - resolve.ties: string or vector of strings denoting the way ties are handled in the 'phylo.network.default' algorithm
    # - remove.internal.nodes: string denoting if and how internal nodes should be removed from phylogenetic tree-derived lineage trees
    # - nodes.list: nested list containing sequences and selected features in a list per node 
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    
    count_edges_between_two_nodes <- function(node1, node2, edge.matrix){
      
      # Counts edges to go from node1 to node2 from a matrix containing edges of a network/tree
      # Arguments:
      # - node1: the starting node (up in the network/tree)
      # - node2: the ending node (down in the network/tree)
      # - edge_matrix: 2-column matrix containing nodes of a directed network in which each row represent an edge, while the first column contains the upper node and the second column contains the lower node
      # Authors: Valentijn Tromp, Daphne van Ginneken
      
      # Save input edge matrix in 'edge_matrix'
      edge_matrix <- edge.matrix

      # Count number of edges by finding path from 'node2' to 'node1' (thereby moving upwards in the network)
      current_node <- node2
      
      # Initialize counter for the edges
      edge_count <- 0
      
      # Keep counting edges until the 'current_node' is equal to 'node1'
      while(sum(current_node == node1) == 0){

        # Find node to which the 'current_node' is connected (the 'current_node' will be present in second column, while the node to which the 'current_node' is connected will be present in the first column)
        next_node <- edge_matrix[edge_matrix[, 2] %in% current_node, 1]
        
        # Update 'edge_count' and 'current_node'
        edge_count <- edge_count+1
        current_node <- next_node
      }
      
      # Return 'edge_count' 
      return(edge_count)
    }
    
    
    igraph_to_phylo <- function(igraph_object){
      
      # Converts an igraph object into phylo format (based on the 'alakazam::graphToPhylo()' function)
      # Arguments:
      # - igraph_object: igraph object
      # Author: Valentijn Tromp, Daphne van Ginneken
      
      # Convert igraph object into 'edges_df' dataframe, in which each row represent one edge
      edges_df <- igraph::as_data_frame(igraph_object)
      
      # For each node, determine the frequency in the 'edges_df', and store in the 'node_counts' table
      node_counts <- table(c(edges_df$to, edges_df$from))
      
      # Using the 'node_counts' table, determine the terminal nodes ('tips')
      tips <- names(node_counts)[node_counts == 1]
        
      # Using the 'node_counts' table, determine the internal nodes ('nodes')
      nodes <- names(node_counts)[node_counts > 1]
      
      # Create a vector 'germline' containing tip nodes that are also parent nodes (germline nodes)
      germline <- tips[tips %in% edges_df$from]
      # If there is a germline node found:
      if(length(germline) > 0){
        # Create extra internal node 
        new_node <- as.character(length(nodes)+1)
        # Add new node to the 'nodes' vector
        nodes <- c(new_node, nodes)
        # Update the 'from' column in 'edges_df' to replace germline nodes with UCA nodes
        edges_df[edges_df$from == germline, ]$from <- new_node
        # Add edges from UCA nodes to their corresponding germline nodes with weight and label as 0
        edges_df <- rbind(edges_df, data.frame(from = new_node, to = germline, weight = 0, label = 0))
      }
      
      # Create a numeric vector 'tipn' containing indices for tip nodes
      tipn <- 1:length(tips)
      names(tipn) <- tips
      
      # Create a numeric vector 'noden' containing indices for internal nodes
      noden <- (length(tips) + 1):(length(tips) + length(nodes))
      names(noden) <- nodes
      
      # Create a renumbering vector 'renumber' that maps original node names to new numeric indices and update 'from' and 'to' columns in the 'edges_df' dataframe with these new numeric indices
      renumber <- c(tipn, noden)
      edges_df$from <- as.numeric(renumber[edges_df$from])
      edges_df$to <- as.numeric(renumber[edges_df$to])
      
      # Initialize 'phylo' as a list, assign edges, edge lengths, tip labels, node labels, and the number of internal nodes, and finally set class attribute to 'phylo'
      phylo <- list()
      phylo$edge <- matrix(cbind(edges_df$from, edges_df$to), ncol = 2)
      phylo$edge.length <- as.numeric(edges_df$weight)
      phylo$tip.label <- tips
      phylo$node.label <- nodes
      phylo$Nnode <- length(nodes)
      class(phylo) <- "phylo"
      
      # Create a list of nodes with their IDs in 'phylo$nodes'
      nnodes <- length(renumber)
      phylo$nodes <- lapply(1:nnodes, function(x) list(id = names(renumber[renumber == x])))
      
      # Ladderize the phylogenetic tree 'phylo' to orient branches
      phylo <- ape::ladderize(phylo, right = FALSE)
      
      # Return the phylo object
      return(phylo)
    }
    
    
    reorder_edges <- function(edges){
      
      # Reorganizes the edges in a dataframe to enable the construction of a directed graph, ensuring that the 'germline' node is placed in the first row and in the 'upper.node' column, and all its descendants are in subsequent rows. 
      # Arguments:
      # - edges: dataframe that contains two columns ('upper.node' and 'lower.node') that contain names of the nodes of the network/tree, whereby each row represent an edge
      # Author: Valentijn Tromp, Daphne van Ginneken
      
      # Retrieve the names of all the nodes in the input dataframe
      nodes <- unique(c(edges$upper.node, edges$lower.node))
      
      # Create new dataframe 'edges_reorganized' to store reorganized edges and select the edges containing the germline node
      edges_reorganized <- edges[edges$upper.node == "germline" | edges$lower.node == "germline", ]
      
      # Make sure that the germline node is in the 'upper.column' and swap the nodes if necessary
      edges_reorganized[edges_reorganized$lower.node == "germline", ] <- edges_reorganized[edges_reorganized$lower.node == "germline", c("lower.node", "upper.node", colnames(edges_reorganized[3:ncol(edges_reorganized)]))]
      
      # Start reordering the nodes in the 'edges' dataframe from the 'germline'
      current_upper_nodes <- edges_reorganized[edges_reorganized$upper.node == "germline", "lower.node"]
      processed_nodes <- c("germline", current_upper_nodes)
      
      # Keep reordering the nodes in the 'edges' dataframe until all nodes are processed and present in the 'processed_nodes' vector
      while(all(nodes %in% processed_nodes) == FALSE){
        
        # Create empty vector to store nodes that will be connected to the nodes in the 'current_upper_nodes' vector
        processed_lower_nodes <- c()
        
        # Iterate through nodes in the 'current_upper_nodes' vector
        for(upper_node in current_upper_nodes){
          
          # Select the rows/edges from 'edges' dataframe that contain the current 'upper_node'
          selected_edges <- edges[edges$upper.node == upper_node | edges$lower.node == upper_node, ]
          
          # Retrieve all the nodes that are present in this selection of edges
          selected_nodes <- unique(c(selected_edges$upper.node, selected_edges$lower.node))
          
          # Remove the nodes that are already processed, the remaining nodes will be the lower nodes of the current 'upper_node'
          current_lower_nodes <- selected_nodes[!selected_nodes %in% processed_nodes]
          
          # Iterate through nodes in 'current_lower_nodes'
          for(lower_node in current_lower_nodes){
            
            # Append edge to the 'edges_reorganized' dataframe
            edges_reorganized <- rbind(edges_reorganized, c(upper_node, lower_node, edges[(edges$upper.node == upper_node & edges$lower.node == lower_node) | (edges$upper.node == lower_node & edges$lower.node == upper_node), colnames(edges_reorganized[3:ncol(edges_reorganized)])]))
          }
          
          # After iterating trough nodes in 'current_lower_nodes', ...
          processed_lower_nodes <- c(processed_lower_nodes, current_lower_nodes)
        }
        
        # All the nodes in 'processed_lower_nodes' will be the 'current_upper_nodes' in the next iteration
        current_upper_nodes <- processed_lower_nodes
        
        # Update 'processed_nodes' vector by appending the nodes in the 'processed_lower_nodes' vector to the 'processed_nodes' vector
        processed_nodes <- c(processed_nodes, processed_lower_nodes)
      }
      
      # Return reorganized dataframe 
      return(edges_reorganized)
    }
    
    
    # 1. Create necessary objects for tree/network construction
    
    # If a distance matrix is provided as input, create a duplicate 
    if(!missing(dist.matrix)){
      dist_matrix <- dist.matrix
      dist_matrix_backup <- dist.matrix
    }
    
    # Create empty list to store warnings during the tree/network construction
    warnings <- list()
    
    # Create empty list to store metrics that are calculated during the tree/network construction
    metrics <- list()
    
    
    # 2. Use specified approach and algorithm to create 'edges' dataframe with the columns "upper.node", "lower.node", and "edge.length" in which each row represents one edge 
    
    # 2.1 If the approach is set to 'network' with the default algorithm, a mst-like network algorithm is used to construct a germline-conditioned minimum spanning tree
    if(approach == "network" && algorithm == "default"){
      
      # Create list of size/frequencies 
      node_sizes <- lapply(names(nodes.list), function(node) nodes.list[[node]][["size"]])
      names(node_sizes) <- names(nodes.list)
      node_sizes["germline"] <- 0
      
      # Create empty matrix to store counts for tie resolving for each node
      tie_resolving <- matrix(data = 0, nrow = nrow(dist.matrix), ncol = length(resolve.ties), dimnames = list(rownames(dist.matrix), resolve.ties))
      
      # Create empty matrix to store edges of the tree 
      edges <- matrix(nrow = 0, ncol = 3)
      colnames(edges) <- c("upper.node", "lower.node", "edge.length")
      
      # The 'phylo.network.default' construction can only be used when there are more than two sequences
      if(nrow(dist_matrix) > 2){
        
        # Replace diagonal 0's by NA
        diag(dist_matrix) <- NA
        
        # Get minimum distance to the germline node
        min_distance <- min(na.exclude(dist_matrix[, "germline"]))
        
        # Select nodes that have this distance to the germline node
        nodes_to_be_connected <- rownames(dist_matrix)[dist_matrix[, "germline"] == min_distance & !is.na(dist_matrix[, "germline"])]
        
        # Add edge(s) between germline and node(s) with 'min_distance' to the germline to the 'edges' matrix
        for(node in nodes_to_be_connected){
          edges <- rbind(edges, c("germline", node, dist_matrix_backup["germline", node]))
        }
        
        # Remove minimum distance to the germline node from the distance matrix
        dist_matrix[, "germline"][dist_matrix[, "germline"] == min_distance] <- NA
        dist_matrix["germline", ][dist_matrix["germline", ] == min_distance] <- NA
        
        # Iteratively add other nodes to the 'edges' matrix until all nodes are present in the 'edges' matrix
        while(!all(rownames(dist_matrix) %in% edges)){
          
          # Select nodes that are already connected (present in the 'edges' matrix) and store these nodes in the 'already_connected_nodes' vector
          already_connected_nodes <- unique(na.exclude(c(edges[, 1], edges[, 2])))
          
          # Get minimum distance to one of the nodes in the 'already_connected_nodes' vector
          min_distance <- min(na.exclude(dist_matrix[, already_connected_nodes]))
          
          # Select unlinked nodes that have this distance to any existing node in 'already_connected_nodes' and store these nodes in the 'nodes_to_be_connected' vector
          nodes_to_be_connected <- unique(unlist(lapply(rownames(dist_matrix), function(x) if(min_distance %in% dist_matrix[x,already_connected_nodes] && !(x %in% already_connected_nodes)){return(x)})))
          
          # Iterate through nodes in 'nodes_to_be_connected'
          for(node in nodes_to_be_connected){
            
            # Select the node(s) in 'already_connected_nodes' that has/have 'min_distance' to the current node
            node_to_connect_to <- unlist(lapply(already_connected_nodes, function(x) if(dist_matrix[x, node] == min_distance){return(x)}))
            
            # If multiple nodes in 'node_to_connect_to' share 'min_distance' to the current node, the 'resolve.ties' options are hierarchically applied to narrow down the candidates for linking the currently unlinked node in the tree to, aiming to select a single node.
            if(length(node_to_connect_to) != 1){
              
              # Loop over provided 'resolve.ties' options
              for(option in resolve.ties){
                
                # Store number of nodes that have 'min_distance' to the current node in 'start_length'
                start_length <- length(node_to_connect_to)
                
                # If the 'resolve.ties' parameter is set to 'min.expansion', the node(s) having the smallest size is/are selected
                if(option == "min.expansion"){
                  node_to_connect_to <- node_to_connect_to[node_sizes[node_to_connect_to] == min(unlist(node_sizes[node_to_connect_to]))]
                }
                
                # If the 'resolve.ties' parameter is set to 'max.expansion', the node(s) having the biggest size is/are selected
                if(option == "max.expansion"){
                  node_to_connect_to <- node_to_connect_to[node_sizes[node_to_connect_to] == max(unlist(node_sizes[node_to_connect_to]))]
                }
                
                # If the 'resolve.ties' parameter is set to 'min.germline.dist', the node(s) having the smallest string distance to the germline node is/are selected
                if(option == "min.germline.dist"){
                  node_to_connect_to <- node_to_connect_to[dist_matrix_backup[node_to_connect_to, "germline"] == min(dist_matrix_backup[node_to_connect_to, "germline"])]
                }
                
                # If the 'resolve.ties' parameter is set to 'max.germline.dist', the node(s) having the biggest string distance to the germline node is/are selected
                if(option == "max.germline.dist"){
                  node_to_connect_to <- node_to_connect_to[dist_matrix_backup[node_to_connect_to, "germline"] == max(dist_matrix_backup[node_to_connect_to, "germline"])]
                }
                
                # If the 'resolve.ties' parameter is set to 'min.germline.edges', the node(s) having the lowest possible number of edges to the germline node is/are selected
                if(option == "min.germline.edges"){
                  number_of_edges <- sapply(node_to_connect_to, function(node) count_edges_between_two_nodes(node1 = "germline", node2 = node, edge.matrix = edges))
                  node_to_connect_to <- node_to_connect_to[number_of_edges == min(number_of_edges)]
                }
                
                # If the 'resolve.ties' parameter is set to 'max.germline.edges', the node(s) having the highest possible number of edges to the germline node is/are selected
                if(option == "max.germline.edges"){
                  number_of_edges <- sapply(node_to_connect_to, function(node) count_edges_between_two_nodes(node1 = "germline", node2 = node, edge.matrix = edges))
                  node_to_connect_to <- node_to_connect_to[number_of_edges == max(number_of_edges)]
                }
                
                # If the 'resolve.ties' parameter is set to 'random', a random node will be selected
                if(option == "random"){
                  node_to_connect_to <- sample(node_to_connect_to, 1)
                }
                
                # If the number of nodes in 'node_to_connect_to' is reduced, add the number of nodes that are excluded to the 'tie_resolving' matrix
                if(length(node_to_connect_to) < start_length){
                  tie_resolving[node, option] <- start_length-length(node_to_connect_to)
                }
              }
            }
            
            # If there is only one node left in 'node_to_connect_to'
            if(length(node_to_connect_to) == 1){
              
              # Add new edge between this linked node in the tree and the currently unlinked node to the 'edges' matrix
              edges <- rbind(edges, c(node_to_connect_to, node, dist_matrix_backup[node_to_connect_to, node]))
            }
            
            # If there are still multiple nodes in 'node_to_connect_to' that has 'min_distance' to the current node...
            if(length(node_to_connect_to) != 1){
              
              # Append warning message to 'warnings' list
              warnings <- c(warnings, paste0("Not all ties could be resolved in ", strsplit(clone, split="_")[[1]][2], " of ", strsplit(clone, split="_")[[1]][1], "!"))
              print("yes")
              # The currently unlinked node is connected to both linked nodes in the tree, thereby creating a cyclic graph
              for(i in node_to_connect_to){
                edges <- rbind(edges, c(i, node, dist_matrix_backup[i, node]))
              }
            }
            
            # Remove this minimum distance from the distance matrix
            dist_matrix[node_to_connect_to, node] <- NA
            dist_matrix[node, node_to_connect_to] <- NA
          }
        }
        
        # Append the 'tie_resolving' matrix to 'metrics' list
        metrics[["tie.resolving"]] <- tie_resolving
      }
      
      # If there are less than three sequences, the tree consists of a germline node and a single descendant
      if(nrow(dist_matrix) < 3){
        edges <- rbind(edges, c("germline", "node1", dist_matrix_backup["germline", "node1"]))
      }
      
      # Convert 'edges' matrix into dataframe
      edges <- as.data.frame(edges)
      
      # Convert 'edges' matrix into igraph object
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
      
      # Set 'phylo_object', 'edges_with_inner_nodes', and 'igraph_object_with_inner_nodes' to NULL
      phylo_object <- NULL
      edges_with_inner_nodes <- NULL
      igraph_object_with_inner_nodes <- NULL
    }
    
    
    # 2.2 If the approach is set to 'network' and the 'mst' algorithm is specified, a minimum spanning tree is constructed with the 'ape::mst()' function and reorganized into a (directed) lineage tree with the B cell on top 
    if(approach == "network" && algorithm == "mst"){
      
      # Create minimum spanning stree with 'ape::mst()' function
      adjacency_matrix <- ape::mst(dist_matrix)
      
      # Convert 'mst' object into 'igraph' object
      igraph_object <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
      
      # Create matrix containing the edges of tree tree between internal nodes (first column) and terminal nodes (second column)
      edges <- igraph::as_edgelist(igraph_object)
      
      # Create dataframe containing the edges plus a third column containing the length of the edges
      edges <- data.frame(do.call(rbind, lapply(1:nrow(edges), function(x) c(edges[x, 1], edges[x, 2], as.numeric(dist_matrix_backup[edges[x, 1], edges[x, 2]])))))
      
      # Rename column names of 'edges' dataframe
      colnames(edges) <- c("upper.node", "lower.node", "edge.length")
      
      # Reorder the 'edges' dataframe to enable the construction of a directed graph
      edges <- reorder_edges(edges)
      
      # Convert reordered 'edges' dataframe back into an object of class 'igraph'
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
      
      # Set 'phylo_object', 'edges_with_inner_nodes', and 'igraph_object_with_inner_nodes' to NULL
      phylo_object <- NULL
      edges_with_inner_nodes <- NULL
      igraph_object_with_inner_nodes <- NULL
    }
    
    
    # 2.3 If the approach is set to 'tree' and the 'nj' algorithm is specified, a neighbor joining tree is constructed with the 'ape::nj()' function  
    if(approach == "tree" && algorithm == "nj"){
      
      # A neighbor joining tree can only be build if there are 3 or more sequences
      if(nrow(dist_matrix) > 2){
        
        # Create neighbor joining tree with 'ape::nj()' function
        phylo_object <- ape::nj(dist_matrix)
      }
      
      # If there are less than 3 sequences, the tree consists of a germline node and a single descendant ('node1'), with a distance that is equal to the string distance from the 'dist_matrix_backup'
      if(nrow(dist_matrix) < 3){
        edges <- data.frame(upper.node = "germline", lower.node = "node1", edge.length = dist_matrix_backup["germline", "node1"])
        igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
        
        # As no 'phylo_object' can be created, internal nodes are absent, and therefore the 'phylo_object', 'igraph_object_with_inner_nodes', and 'edges_with_inner_nodes' are set to NULL
        phylo_object <- NULL
        igraph_object_with_inner_nodes <- NULL
        edges_with_inner_nodes <- NULL
      }
    }
    
    
    # 2.4 If the approach is set to 'tree' and the 'mp' algorithm is specified, a maximum parsimony tree is constructed with the 'phangorn::pratchet()' function
    if(approach == "tree" && algorithm == "mp"){
      
      # A maximum parsimony tree can only be build if there are 3 or more sequences
      if(length(msa) > 2){
        
        # Create maximum parsimony trees with the 'phangorn::pratchet()' function (and hide in-function printed message)
        base::invisible(utils::capture.output(phylo_object <- phangorn::pratchet(msa)))
        
        # Assign branch lengths to the trees 
        phylo_object <- phangorn::acctran(phylo_object, msa)
        
        # Remove node labels to enable conversion into igraph object
        phylo_object$node.label <- NULL
      }
      
      # If there are less than 3 sequences, the tree consists of a germline node and a single descendant ('node1'), with a distance that is equal to the hamming distance
      if(length(msa) < 3){
        edges <- data.frame(upper.node = "germline", lower.node = "node1", edge.length = stringdist::stringdist(as.character(phangorn::as.MultipleAlignment(msa))[1], as.character(phangorn::as.MultipleAlignment(msa))[2], method = "hamming"))
        igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
        
        # As no 'phylo_object' can be created, internal nodes are absent, and therefore the 'phylo_object', 'igraph_object_with_inner_nodes', and 'edges_with_inner_nodes' are set to NULL
        phylo_object <- NULL
        igraph_object_with_inner_nodes <- NULL
        edges_with_inner_nodes <- NULL
      }
    }
    
    
    # 2.5 If the approach is set to 'tree' and the 'ml' algorithm is specified, a maximum likelihood tree is constructed with the 'phangorn::pml_bb()' function
    if(approach == "tree" && algorithm == "ml"){
      
      # A maximum likelihood tree can only be build if there are 3 or more sequences
      if(length(msa) > 2){
        
        # If there are over 4 nodes, next to the germline, bootstrapping is allowed during maximum likelihood tree inference, and 'rearrangement_method' is set to 'stochastic'
        if(length(msa) > 4){rearrangement_method <- "stochastic"}
        
        # If there are 4 nodes or less, no bootstrapping is conducted, to prevent the creation of trees with less than three edges (such a tree cannot be unrooted)
        if(length(msa) <= 4){rearrangement_method <- "NNI"}
        
        # Create maximum likelihood tree with the 'phangorn::pml_bb()' function (and hide in-function printed message), and store the best model according to the BIC value from the 'phangorn::modelTest()' function in the 'metrics' list
        if(sequence.type == "DNA"){base::invisible(utils::capture.output(ml_tree <- phangorn::pml_bb(phangorn::modelTest(object = msa, model = dna.model), rearrangement = rearrangement_method))); metrics["dna.model"] <- ml_tree[["model"]]}
        if(sequence.type == "AA"){base::invisible(utils::capture.output(ml_tree <- phangorn::pml_bb(phangorn::modelTest(object = msa, model = aa.model), rearrangement = rearrangement_method))); metrics["aa.model"] <- ml_tree[["model"]]}
  
        # Store the tree of class 'phylo' in the 'phylo_object'
        phylo_object <- ml_tree[["tree"]]
        
        # If codon models are specified in the 'codon.model' parameter...
        if(!all(is.na(codon.model)) & !is.null(phylo_object)){
          
          # Convert the multiple sequence alignment into a matrix and remove '-' characters and save in 'codon_msa'
          codon_msa <- do.call(rbind, lapply(names(msa), function(x){strsplit(gsub(pattern = "-", replacement = "", do.call(paste0, as.list(as.character(msa)[x,]))), split = "")[[1]]}))
          rownames(codon_msa) <- names(msa)
          
          # Convert the 'codon_msa' into the phyDat format, which will be the input for the 'phangorn::codonTest() function
          codon_msa <- phangorn::as.phyDat(codon_msa, type = "CODON")
          
          # Estimate the specified codon models
          codon_model_test <- phangorn::codonTest(tree = phylo_object, object = codon_msa, model = codon.model)
          
          # Select the best model according to the BIC value and append to the 'metrics' list
          codon.model <- codon_model_test$summary[codon_model_test$summary$BIC == min(codon_model_test$summary$BIC), "model"]
          metrics["codon.model"] <- codon.model
          
          # Replace the tree of class 'phylo' stored in 'phylo_object' by the tree created using the selected codon model
          phylo_object <- codon_model_test$estimates[[codon.model]]$tree
        }
        
        # Remove node labels to enable conversion into igraph object
        phylo_object$node.label <- NULL
      }
      
      # If there are less than 3 sequences, the tree consists of a germline node and a single descendant ('node1'), with a distance that is equal to the hamming distance normalized by the mean sequence length
      if(length(msa) < 3){
        edges <- data.frame(upper.node = "germline", lower.node = "node1", edge.length = stringdist::stringdist(as.character(phangorn::as.MultipleAlignment(msa))[1], as.character(phangorn::as.MultipleAlignment(msa))[2], method = "hamming")/((nchar(as.character(phangorn::as.MultipleAlignment(msa))[1]) + nchar(as.character(phangorn::as.MultipleAlignment(msa))[2]))/2))
        igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
        
        # As no 'phylo_object' can be created, internal nodes are absent, and therefore the 'phylo_object', 'igraph_object_with_inner_nodes', and 'edges_with_inner_nodes' are set to NULL
        phylo_object <- NULL
        igraph_object_with_inner_nodes <- NULL
        edges_with_inner_nodes <- NULL
      }
    }
    
    
    # 2.6 If a list of IgPhyML trees is provided, no tree/network construction is taking place, and trees are retrieved from the IgPhyML output file and processed into the AntibodyForests object
    if(approach == "tree" && algorithm == "IgPhyML"){
      
      # Retrieve tree for the current clonotype
      IgPhyML_igraph_object <- IgPhyML.trees[["trees"]][[clone]]
      
      # If an igraph object could be found in the 'IgPhyML.trees' list...
      if(!is.null(IgPhyML_igraph_object)){
      
        # Remove prefixes and suffixes from the barcodes (if present)
        igraph::V(IgPhyML_igraph_object)$name <- ifelse(grepl("[ACTG]", igraph::V(IgPhyML_igraph_object)$name) & nchar(gsub("[^ACTG]", "", igraph::V(IgPhyML_igraph_object)$name)) >= 10, gsub("[^ACTG]", "", igraph::V(IgPhyML_igraph_object)$name), igraph::V(IgPhyML_igraph_object)$name)
        
        # Rename the names of the nodes in the igraph object 
        igraph::V(IgPhyML_igraph_object)$name <- sapply(igraph::V(IgPhyML_igraph_object)$name, function(x){
          # If a node labels contains 'GERM', the name of the node is set to 'germline'
          if(grepl(pattern = "GERM", x)){return("germline")}
          # If a node labels refers to a barcode, the 'nodes.list' is used to retrieve the node that contains this barcode/cell
          else if(grepl(pattern = "^[ACGT]+$", x)){for(node in names(nodes.list)){if(x %in% nodes.list[[node]][["barcodes"]]){return(node)}}}
          # If the node name does not contain 'GERM' and does not refer to a barcode, the node is an intermediate node and the name of the node remains the same
          else if(grepl(pattern = "^[0-9]+$", x)){return(as.vector(1:sum(grepl(pattern = "^[0-9]+$", igraph::V(IgPhyML_igraph_object)$name)))[grep(pattern = "^[0-9]+$", igraph::V(IgPhyML_igraph_object)$name, value = TRUE) == x])}
        })
        
        # Check whether all nodes in the 'nodes.list' could be found in the igraph object, and check whether all terminal nodes are found once
        igraph_object_check <- (all(names(nodes.list) %in% igraph::V(IgPhyML_igraph_object)$name) && all(table(igraph::V(IgPhyML_igraph_object)$name)[unique(grep(pattern = "^node|germline", igraph::V(IgPhyML_igraph_object)$name, value = TRUE))] == 1))
      } else{igraph_object_check <- FALSE}
      
      # If no igraph object could be found in the 'IgPhyML.trees' object, append warning message to warnings 'list'
      warnings <- c(warnings, paste(c("The tree of", strsplit(clone, split="_")[[1]][], "of", strsplit(clone, split="_")[[1]][1], "could not be found in the specified IgPhyML output file!"), collapse = ""))
      
      # If the igraph object is checked succesfully, onvert the object of class 'igraph' into an object of class 'phylo' using the 'igraph_to_phylo()' function
      if(igraph_object_check){phylo_object <- igraph_to_phylo(IgPhyML_igraph_object)}
      
      # If not, the 'phylo_object' and all other objects that will be derived from it are set to NULL
      if(!igraph_object_check){phylo_object <- NULL; igraph_object <- NULL; igraph_object_with_inner_nodes <- NULL; edges <- NULL; edges_with_inner_nodes <- NULL}
      }
    
    
    # 3.1 If the approach is set to 'tree', the'phylo_object' will be converted into an object of class 'igraph', along with the creation of a dataframe containing the edges (this dataframe serves as input for the internal nodes removal algorithm)
    if(approach == "tree" && !is.null(phylo_object)){
      
      # The 'phylo_object' (created by the 'ape::nj()', 'phangorn::pratchet()', or 'phangorn::pml_bb()' function) is converted into an object of class 'igraph'
      igraph_object <- ape::as.igraph.phylo(phylo_object)
      
      # Remove 'Node' from names of internal nodes
      igraph::V(igraph_object)$name <- gsub(pattern = "Node", replacement = "", x = igraph::V(igraph_object)$name)
      
      # Create matrix containing the edges of tree tree between internal nodes (first column) and terminal nodes (second column)
      edges <- igraph::as_edgelist(igraph_object)
      
      # Create dataframe containing the edges plus a third column containing the length of the edges
      edges <- data.frame(do.call(rbind, lapply(1:nrow(edges), function(x) c(edges[x, 1], edges[x, 2], as.numeric(abs(phylo_object$edge.length[x]))))))
      
      # Rename column names of 'edges' dataframe
      colnames(edges) <- c("upper.node", "lower.node", "edge.length")
      
      # Reorder the 'edges' dataframe to enable the construction of a directed graph
      edges <- reorder_edges(edges)
      
      # Create a duplicate of 'edges' to return later as 'edges_with_inner_nodes'
      edges_with_inner_nodes <- edges
      
      # Create pairwise distance matrix containing all the distances between the nodes, including the internal nodes
      dist_matrix_tree <- ape::dist.nodes(phylo_object)
      
      # Rename column names and row names of 'dist_matrix_NJ'
      rownames(dist_matrix_tree) <- c(phylo_object$tip.label, 1:phylo_object$Nnode)[as.numeric(rownames(dist_matrix_tree))]
      colnames(dist_matrix_tree) <- c(phylo_object$tip.label, 1:phylo_object$Nnode)[as.numeric(colnames(dist_matrix_tree))]
      
      # Convert 'edges' dataframe (without edges of length zero) back into a directed graph of class 'igraph'
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
      
      # Create a duplicate of 'igraph_object' to return later as 'igraph_object_with_inner_nodes'
      igraph_object_with_inner_nodes <- igraph_object
    }
    
    
    # 3.2 If the approach is set to 'tree', edges with a length of zero (if present) or below 1e-8 are removed by replacing the unrecovered internal node (referred to as the outgoing node) in the 'upper.node' column with the sequence recovered terminal node (referred to as the substitute node) in the 'lower.node' by default
    if(approach == "tree" && !is.null(phylo_object) && remove.internal.nodes %in% c("zero.length.edges.only", "connect.to.parent", "minimum.length", "minimum.weight")){
     
      # Make sure the values in the column 'edge.length' are of class 'numeric'
      edges$edge.length <- as.numeric(edges$edge.length)
      
      # Make a subset of edges with a length of 0 from the 'edges' dataframe
      zero_length_edges <- edges[edges$edge.length < 0.0000001,]
      
      # Remove edges with a length of 0 from 'edges' dataframe
      edges <- edges[!rownames(edges) %in% rownames(zero_length_edges), ]
      
      # If the 'zero_length_edges' is not empty...
      if(nrow(zero_length_edges) != 0){
        
        # Loop over edges in 'zero_length_edges' dataframe
        for(i in 1:nrow(zero_length_edges)){
          
          # Select node names connected by the current edge
          node_pair <- as.character(zero_length_edges[i, c("upper.node", "lower.node")])
          
          # If one of the nodes is a terminal sequence recovered node, this node is selected as the 'substitute_node'
          substitute_node <- grep("node|germline", node_pair, value = TRUE)
          
          # If there is a node in 'substitute_node', the other node in the 'node_pair' vector will be selected as the 'outgoing_node', which will be replaced by the 'substitue_node'
          if(length(substitute_node) == 1){
            outgoing_node <- node_pair[!node_pair %in% substitute_node]
          }
          
          # If both nodes in the 'node_pair' vector are internal sequence unrecovered nodes, the node in the 'lower.node' column is selected as the 'substitude_node', and the node in the 'upper.node' column is selected as the 'outgoing_node'
          if(length(substitute_node) == 0){
            substitute_node <- zero_length_edges[i, "lower.node"]
            outgoing_node <- zero_length_edges[i, "upper.node"]
          }
          
          # Replace 'outgoing_node' with 'substitute.node' in the columns 'upper.node' and 'lower.node' and simultaneously update the 'edge.length' column in the 'edges' dataframe
          for(row in rownames(edges[edges$upper.node == outgoing_node, ])){edges[row, c("upper.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "lower.node"]])} 
          for(row in rownames(edges[edges$lower.node == outgoing_node, ])){edges[row, c("lower.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "upper.node"]])}
          
          # Replace 'outgoing_node' with 'substitute.node' in the columns 'upper.node' and 'lower.node' and simultaneously update the 'edge.length' column in the 'zero_length_edges' dataframe
          for(row in rownames(zero_length_edges[zero_length_edges$upper.node == outgoing_node, ])){zero_length_edges[row, c("upper.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, zero_length_edges[row, "lower.node"]])}
          for(row in rownames(zero_length_edges[zero_length_edges$lower.node == outgoing_node, ])){zero_length_edges[row, c("lower.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, zero_length_edges[row, "upper.node"]])}
        }
      }
      
      # Convert 'edges' dataframe (without edges of length zero) back into an object of class 'igraph'
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
    }
    
    
    # 3.3 If the approach is set to 'tree' and 'remove.internal.nodes' is set to 'connect.to.parent', internal nodes (or unrecovered sequences nodes) are removed from the tree by linking terminal sequence recovered nodes to the first upper parental sequence recovered node (resulting in mostly germline-directed treesS)
    if(approach == "tree" && !is.null(phylo_object) && remove.internal.nodes == "connect.to.parent"){
      
      # Retrieve internal nodes present in 'edges' dataframe (names of terminal nodes start with "node" followed by a number, the germline node is named "germline", and all other nodes are internal nodes)
      internal_nodes <- unique(c(edges$upper.node, edges$lower.node)[!startsWith(c(edges$upper.node, edges$lower.node), "node") & c(edges$upper.node, edges$lower.node) != "germline"])
      
      # Iterate through internal nodes
      for(internal_node in internal_nodes){
        
        # Select the rows from the 'edges' dataframe that contain the current 'internal_node'
        selection_edges_to_be_removed <- edges[edges$upper.node == internal_node | edges$lower.node == internal_node, ]
        
        # Remove the edges in 'selection_edges_to_be_removed' from the 'edges' dataframe
        edges <- edges[!rownames(edges) %in% rownames(selection_edges_to_be_removed), ]
        
        # Select the upper (parental) sequence recovered node from the current 'internal_node' 
        upper_node <- selection_edges_to_be_removed[selection_edges_to_be_removed$lower.node == internal_node, "upper.node"]
        
        # Select the lower (descendant) sequence recovered nodes from the current 'internal_node'
        lower_nodes <- selection_edges_to_be_removed[selection_edges_to_be_removed$upper.node == internal_node, "lower.node"]
        
        # Loop over the nodes in the 'lower_nodes' vector
        for(lower_node in lower_nodes){
          
          # Append the new edge between the sequence recovered 'upper_node' and the current 'lower node' to the 'edges' dataframe, together with the distance between these nodes
          edges <- rbind(edges, c(upper_node, lower_node, dist_matrix_tree[upper_node, lower_node]))
        }
      }
      
      # Convert 'edges' dataframe (without edges of length zero) back into an object of class 'igraph'
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
    }
    
    
    
    # 3.4 If the approach is set to 'tree' and 'remove.internal.nodes' is set to 'all', internal nodes (or unrecovered sequences nodes) are removed from the tree by iteratively replacing internal nodes with terminal nodes with the minimum increase in the sum of all edges (this increase is referred to as the 'cost')
    if(approach == "tree" && !is.null(phylo_object) && remove.internal.nodes %in% c("minimum.length", "minimum.cost")){
      
      # Retrieve internal nodes present in 'edges' dataframe (names of terminal nodes start with "node" followed by a number, the germline node is named "germline", and all other nodes are internal nodes)
      internal_nodes <- unique(c(edges$upper.node, edges$lower.node)[!startsWith(c(edges$upper.node, edges$lower.node), "node") & c(edges$upper.node, edges$lower.node) != "germline"])
      
      # Keep removing nodes until the 'internal_nodes' is empty
      while(length(internal_nodes) != 0){
        
        # Select edges from 'edges' dataframe which link an node from the 'internal_nodes' vector to another node (could be another internal node or a terminal node)
        edges_to_be_removed <- edges[edges$upper.node %in% internal_nodes | edges$lower.node %in% internal_nodes, ]
        
        # Apply function that calculated the cost to remove an edge to each edge/row in 'edges_to_be_removed', in which the cost is defined as the increase in the sum of edge lengths
        costs <- unlist(lapply(1:nrow(edges_to_be_removed), function(row){
          
          # Retrieve nodes that are connected by this edge
          node_pair <- c(edges_to_be_removed$upper.node[row], edges_to_be_removed$lower.node[row])
          
          # If one of the nodes is a terminal (or sequence recovered) node, the internal (or unrecovered sequence) node is replaced by this node
          if(!all(node_pair %in% internal_nodes)){
            outgoing_node <- node_pair[node_pair %in% internal_nodes]
            substitute_node <- node_pair[!node_pair %in% internal_nodes]
          }
          
          # If both nodes are internal (or unrecovered sequence) nodes, the upper internal node is replaced by the lower internal node
          if(all(node_pair %in% internal_nodes)){
            outgoing_node <- node_pair[1]
            substitute_node <- node_pair[2]
          }
          
          # Retrieve nodes that are connected to the outgoing node
          nodes_linked_to_outgoing_node <- unique(c(edges[edges$upper.node == outgoing_node | edges$lower.node == outgoing_node, "upper.node"], c(edges[edges$upper.node == outgoing_node | edges$lower.node == outgoing_node, "lower.node"])))
          nodes_linked_to_outgoing_node <- nodes_linked_to_outgoing_node[nodes_linked_to_outgoing_node != outgoing_node & nodes_linked_to_outgoing_node != substitute_node]
          
          # Retrieve nodes that are connected to the substitute node
          nodes_linked_to_substitute_node <- unique(c(edges[edges$upper.node == substitute_node | edges$lower.node == substitute_node, "upper.node"], c(edges[edges$upper.node == substitute_node | edges$lower.node == substitute_node, "lower.node"])))
          nodes_linked_to_substitute_node <- nodes_linked_to_substitute_node[nodes_linked_to_substitute_node != outgoing_node & nodes_linked_to_substitute_node != substitute_node]
          
          # Calculate the cost for removing the edge (by subtracting the sum of the lengths of the edges currently linked to the 'outgoing_node' and 'substitute_node' from the sum of the lengths of the edges that would be linked to the 'substitute_node' after the 'outgoing_node' is replaced with the 'substitute_node')
          cost <- sum(sapply(nodes_linked_to_outgoing_node , function(node) dist_matrix_tree[node, substitute_node])) - sum(sapply(nodes_linked_to_outgoing_node , function(node) dist_matrix_tree[node, outgoing_node])) - dist_matrix_tree[substitute_node, outgoing_node]
          
          # Return cost
          return(cost)
        }))
        
        # If 'remove.internal.nodes' is set to 'minimum.length', select the edges from the 'edges_to_be_removed' dataframe that have the minimum 'edge.length' and store these edges in 'selection_edges_to_be_removed'
        if(remove.internal.nodes == "minimum.length"){selection_edges_to_be_removed <- edges_to_be_removed[edges_to_be_removed$edge.length == min(edges_to_be_removed$edge.length), ]}
        
        # If 'remove.internal.nodes' is set to 'minimum.cost', select the edges from the 'edges_to_be_removed' dataframe that have the minimum cost to get removed and store these edges in 'selection_edges_to_be_removed'
        if(remove.internal.nodes == "minimum.cost"){selection_edges_to_be_removed <- edges_to_be_removed[(1:nrow(edges_to_be_removed))[costs == min(costs)], ]}
        
        # If there is only one edge with this minimum cost...
        if(nrow(selection_edges_to_be_removed) == 1){
          
          # Remove edge in 'selection_edges_to_be_removed' from 'edges' dataframe
          edges <- edges[!rownames(edges) %in% rownames(selection_edges_to_be_removed), ]
         
          # Retrieve nodes that are connected by this edge
          node_pair <- c(selection_edges_to_be_removed$upper.node, selection_edges_to_be_removed$lower.node)
          
          # If one of the nodes is a terminal (or sequence recovered) node, the internal (or unrecovered sequence) node is replaced by this node
          if(!all(node_pair %in% internal_nodes)){
            outgoing_node <- node_pair[node_pair %in% internal_nodes]
            substitute_node <- node_pair[!node_pair %in% internal_nodes]
          }
          
          # If both nodes are internal (or unrecovered sequence) nodes, the upper internal node is replaced by the lower internal node
          if(all(node_pair %in% internal_nodes)){
            outgoing_node <- node_pair[1]
            substitute_node <- node_pair[2]
          }
          
          # Replace all occurrences of 'outgoing_node' with 'substitute_node' in 'edges' dataframe and simultaneously update 'edge.length'
          for(row in rownames(edges[edges$upper.node == outgoing_node, ])){edges[row, c("upper.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "lower.node"]])} 
          for(row in rownames(edges[edges$lower.node == outgoing_node, ])){edges[row, c("lower.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "upper.node"]])}
        }
        
        # If there are multiple edges with this minimum cost...
        if(nrow(selection_edges_to_be_removed) != 1){
          
          # Create vector 'outgoing_nodes' containing the nodes to be replaced that are present in the 'selection_edges_to_be_removed' dataframe
          outgoing_nodes <- unique(unlist(lapply(1:nrow(selection_edges_to_be_removed), function(x){
            
            # Retrieve the nodes that are connected by this edge
            node_pair <- c(selection_edges_to_be_removed$upper.node[x], selection_edges_to_be_removed$lower.node[x])
            
            # If one of the nodes is a terminal (or sequence recovered) node, the internal (or unrecovered sequence) node is replaced by this node
            if(!all(node_pair %in% internal_nodes)){
              outgoing_node <- node_pair[node_pair %in% internal_nodes]
            }
            
            # If both nodes are internal (or unrecovered sequence) nodes, the upper internal node is replaced by the lower internal node
            if(all(node_pair %in% internal_nodes)){
              outgoing_node <- node_pair[1]
            }
            
            # Return the 'outgoing_node'
            return(outgoing_node)
          })))
          
          
          # Loop over nodes in 'outgoing_nodes'
          for(outgoing_node in outgoing_nodes){
            
            # Create vecctor 'substitute_nodes' containing the nodes that are connected to the current outgoing internal node by the edges in 'selection_edges_to_be_removed' with the minimum cost
            substitute_nodes <- unique(c(selection_edges_to_be_removed[selection_edges_to_be_removed$upper.node == outgoing_node, "lower.node"], selection_edges_to_be_removed[selection_edges_to_be_removed$lower.node == outgoing_node, "upper.node"]))
            
            # If there is only node connected to this outgoing node with the minimum cost...
            if(length(substitute_nodes) == 1){
              
              # Save single nodes from 'substitute_nodes' in 'substitute_node'
              substitute_node <- substitute_nodes
              
              # Remove the edge between the 'outgoing_node' and the 'substitute_node' from the 'edges' dataframe and the 'selection_edges_to_be_removed' dataframe
              edges <- edges[!(edges$upper.node == outgoing_node & edges$lower.node == substitute_node), ]
              selection_edges_to_be_removed <- selection_edges_to_be_removed[!(selection_edges_to_be_removed$upper.node == outgoing_node & selection_edges_to_be_removed$lower.node == substitute_node), ]
              
              # Replace all occurrences of 'outgoing_node' with 'substitute_node' in 'edges' dataframe and simultaneously update 'edge.length' in 'edges' dataframe
              for(row in rownames(edges[edges$upper.node == outgoing_node, ])){edges[row, c("upper.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "lower.node"]])} 
              for(row in rownames(edges[edges$lower.node == outgoing_node, ])){edges[row, c("lower.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, edges[row, "upper.node"]])}
              
              # Replace all occurrences of 'outgoing_node' with 'substitute_node' in 'edges' dataframe and simultaneously update 'edge.length' in 'selection_edges_to_be_removed' dataframe
              for(row in rownames(selection_edges_to_be_removed[selection_edges_to_be_removed$upper.node == outgoing_node, ])){selection_edges_to_be_removed[row, c("upper.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, selection_edges_to_be_removed[row, "lower.node"]])} 
              for(row in rownames(selection_edges_to_be_removed[selection_edges_to_be_removed$lower.node == outgoing_node, ])){selection_edges_to_be_removed[row, c("lower.node", "edge.length")] <- c(substitute_node, dist_matrix_tree[substitute_node, selection_edges_to_be_removed[row, "upper.node"]])}
            }
            
            # If all the edges that link the nodes that are currently linked to this outgoing internal node are present in the 'selection_edges_to_be_removed' dataframe, the lower nodes of this outgoing internal node are coupled to the upper node of this outgoing internal node
            if(length(substitute_nodes) > 1){
              
              # Retrieve the upper node of the outgoing internal node from the 'selection_edges_to_be_removed' dataframe and save in 'upper_node' 
              upper_node <- edges[edges$lower.node == outgoing_node, "upper.node"]
              
              # Retrieve the lower nodes of the outgoing internal node from the 'selection_edges_to_be_removed' dataframe and save in 'lower_nodes'
              lower_nodes <- edges[edges$upper.node == outgoing_node, "lower.node"]
              
              # Remove the edges linked to the 'outgoing_node'
              edges <- edges[!(edges$upper.node == outgoing_node | edges$lower.node == outgoing_node), ]
              
              # Append new edges between upper node and lower nodes to the 'edges' dataframe
              for(lower_node in lower_nodes){edges <- base::rbind(edges, c(upper_node, lower_node, dist_matrix_tree[upper_node, lower_node]))}
              colnames(edges) <- c("upper.node", "lower.node", "edge.length")
            }
          }
        }
        
        # Update 'internal_nodes' vector
        internal_nodes <- unique(c(edges$upper.node, edges$lower.node)[!startsWith(c(edges$upper.node, edges$lower.node), "node") & c(edges$upper.node, edges$lower.node) != "germline"])
      }
      
      # Convert 'edges' dataframe (without edges of length zero) back into an object of class 'igraph'
      igraph_object <- igraph::graph_from_data_frame(edges, directed = TRUE)
    }
    
    
    # 4. Combine the 'phylo_object', the 'igraph_object', the 'metrics' list, and the 'warnings' list in one list
    output_list <- list(phylo = phylo_object, 
                        igraph = igraph_object, 
                        igraph.with.inner.nodes = igraph_object_with_inner_nodes, 
                        edges = edges, 
                        edges.with.inner.nodes = edges_with_inner_nodes,
                        metrics = metrics, 
                        warnings = warnings)
    
    # Return the 'output_list'
    return(output_list)
  }
  
  
  infer_single_network <- function(VDJ,
                                   clone,
                                   sequence.columns,
                                   germline.columns,
                                   sequence.type,
                                   concatenate.sequences,
                                   node.features,
                                   construction.method,
                                   IgPhyML.trees,
                                   string.dist.metric,
                                   dna.model,
                                   aa.model,
                                   codon.model,
                                   resolve.ties,
                                   remove.internal.nodes,
                                   include){
    
    # Infers network/tree for a single clone within one sample
    # Arguments:
    # - VDJ: VDJ dataframe as obtained from the minimal_VDJ() function in Platypus
    # - clone: string denoting the sample ID and the clonotype ID for which the network will be inferred (in the format of "S1_clonotype4" or "S2_clonotype3")
    # - sequence.columns: string or vector of strings denoting the sequence columns in the 'VDJ' dataframe that contain the sequences that will be used to infer B cell lineage trees
    # - germline.columns: string or vector of strings denoting the germline columns in the 'VDJ' dataframe that contain the sequences that will be used as starting points of the lineage trees
    # - sequence.type: type of the sequences aligned in the selected 'sequence.columns' and 'germline.columns'
    # - concatenate.sequences: bool indicating whether sequences from different sequence columns should be concatenated before pairwise distance calculations or multiple string alignment
    # - node.features: string or vector of strings denoting the additional column name(s) in the VDJ dataframe to be selected
    # - construction.method: string denoting the network algorithm that will be used to convert the distance matrices or multiple sequence alignments into lineage trees
    # - string.dist.metric: string denoting the metric that will be calculated to measure (string) distances between sequences
    # - dna.model: denotes the evolutionary model that will be used during the pairwise distance matrix calculation or the nucleotide substitution model(s) that will be used during the likelihood calculations
    # - aa.model:  denotes the evolutionary model that will be used during the pairwise distance matrix calculation or the amino acid substitution model(s) that will be used during the likelihood calculations
    # - codon.model: denotes the codon substitution model that can be used during the likelihood calculations
    # - resolve.ties: string denoting the way ties are handled when the 'network.tree' network algorithm is being used 
    # - remove.internal.nodes: bool indicating whether to remove internal nodes from trees constructed with 'phylo.nj' algorithm
    # - include: vector of strings specifying the output object that should be included in the output AntibodyForests object
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # 1. Create a nested list 'nodes_list', in which each sublist represents a clone/node in the graph, and each node represent a unique combination of sequences of the columns selected in 'sequence.columns'
    # Each node/sublist contains the following items:
    # - a (separate) string for all columns selected in 'sequence.columns' containing the sequence for that node
    # - barcodes: a list of all barcodes of the cells that have the unique combination of sequences
    # - frequency: an integer that corresponds to the number of cells (equal to the length of the list 'barcodes')
    # - a (separate) list for all columns selected in 'node.features'
    
    # Retrieve sample ID and clonotype ID from 'clone'
    sample <- strsplit(clone, split="_")[[1]][1]
    clonotype <- strsplit(clone, split="_")[[1]][2]
    
    # Store 'sequence.columns' and 'germline.columns' in vectors 'sequence_columns' and 'germline_columns', respectively
    sequence_columns <- sequence.columns
    germline_columns <- germline.columns
    
    # Make a subset of the VDJ dataframe with all cells from the specified sample and clone and only the columns 'barcode', the specified 'sequence.columns' and 'germline.columns', and additional columns specified in 'node.features'
    VDJ_subset <- VDJ[VDJ$sample_id == sample & VDJ$clonotype_id == clonotype, c("barcode", sequence_columns, germline_columns, node.features)]
    
    # If the 'sequence.type' is set to 'DNA'...
    if(sequence.type == "DNA"){
      # Trim of last one or two nucleotides to make the length of the nucleotide sequences in 'VDJ_subset' multiple of 3
      VDJ_subset[c(sequence_columns, germline_columns)] <- lapply(VDJ_subset[c(sequence_columns, germline_columns)], function(x){
        ifelse(nchar(x) %% 3 != 0, yes = substring(x, 1, nchar(x) - (nchar(x) %% 3)), no = x)
    })}
    
    # If 'concatenate.sequences' is set to TRUE...
    if(concatenate.sequences){
      
      # Concatenate sequences in 'sequence_columns' and 'germline_columns' and store concatenated sequences in 'concatenated_sequence' and 'concatenated_germline', respectively
      VDJ_subset$concatenated_sequence <- sapply(1:nrow(VDJ_subset), function(x) do.call(paste0, as.list(VDJ_subset[x, sequence_columns])))
      VDJ_subset$concatenated_germline <- sapply(1:nrow(VDJ_subset), function(x) do.call(paste0, as.list(VDJ_subset[x, germline_columns])))
      
      # Update column names in 'sequence_columns' and 'germline_columns' vectors
      sequence_columns <- "concatenated_sequence"
      germline_columns <- "concatenated_germline"
    }
    
    # Create a dataframe containing unique combinations of the sequences in 'sequence_columns'
    intraclonal_sequences_df <- as.data.frame(unique(VDJ_subset[,sequence_columns]))
    colnames(intraclonal_sequences_df) <- sequence_columns
    
    # Create a list of sublists where each sublist represents a unique combination of sequences (and each unique combination of sequences will represent a node in the graph later on)
    nodes_list <- lapply(1:nrow(intraclonal_sequences_df), function(x){
      # Extract sequences for the current row and save in 'node_sublist'
      node_sublist <- lapply(sequence_columns, function(y) intraclonal_sequences_df[x,y])
      # Rename indices of 'node_sublist' to 'sequence_columns'
      names(node_sublist) <- sequence_columns
      # Find matching rows in 'VDJ_subset' based on the current combination of sequences
      matching_rows <- sapply(1:nrow(VDJ_subset), function(y) all(intraclonal_sequences_df[x,] == VDJ_subset[y,sequence_columns])) 
      # Add barcodes and size/frequency information to 'node_sublist'
      node_sublist["barcodes"] <- list(VDJ_subset$barcode[matching_rows])
      node_sublist["size"] <- sum(matching_rows)
      # Add 'node.features' to 'node_sublist' as well
      for(i in node.features){node_sublist[i] <- list(VDJ_subset[,i][matching_rows])}
      # Return list of sublists 
      return(node_sublist)
    })
    
    # Order the list based on the size/frequency of the (sub)clones/nodes in descending order
    nodes_list <- nodes_list[order(sapply(nodes_list, function(x) x$size), decreasing = TRUE)]
    
    # Rename the sublists to 'node1', 'node2', etc.
    names(nodes_list) <- paste0("node", seq_along(nodes_list))
    
    # Add a 'germline' node to the list, containing the the most abundant germline sequences from the specified 'germline_columns' (most of the time, each germline column contains only one unique sequence)
    nodes_list$germline <- lapply(germline_columns, function(x) names(table(VDJ_subset[,x]))[which.max(table(VDJ_subset[,x]))])
    names(nodes_list$germline) <- sequence_columns
    
    
    # 2. Create a list 'dist_msa' with the distance matrices and multiple sequence alignments for the current clonotype, which will contain the input for the 'build_lineage_tree()' function in the next step
    
    # Create a list 'dist_msa' to store the distance matrix and multiple sequence alignment per column specified in the 'sequence.columns' vector
    dist_msa <- lapply(sequence_columns, function(x){
      
      # Extract sequences from 'nodes_list' of one 'sequence.column' and save in 'seqs' list
      seqs <- sapply(nodes_list, function(y) y[x])
      names(seqs) <- names(nodes_list)
      
      # If a pairwise distance matrix needs to be computed, while the 'string.dist.metric' parameter is specified, the distance matrix is computed with the 'stringdist::stringdistmatrix()' function using the 'seqs' list as input
      if((construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && !is.na(string.dist.metric))){dist_matrix <- suppressWarnings(stringdist::stringdistmatrix(seqs, seqs, method = string.dist.metric, useNames = "names"))}
      
      # If a multiple sequence alignment needs to be made, the 'msa::msa()' function is used
      if((construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && (!is.na(dna.model) | !is.na(aa.model))) | construction.method %in% c("phylo.tree.mp", "phylo.tree.ml")){
        
        # If the 'sequence.type' is set to "DNA"...
        if(sequence.type == "DNA"){
          
          # Convert character strings in 'seqs' list into DNAString objects
          seqs <- lapply(seqs, function(x) Biostrings::DNAString(x))
          
          # Convers 'seqs' list into an DNAStringSet object
          seqs <- Biostrings::DNAStringSet(seqs)
        }
        
        # If the 'sequence.type' is set to "AA"...
        if(sequence.type == "AA"){
          
          # Convert character strings in 'seqs' list into AAString objects
          seqs <- lapply(seqs, function(x) Biostrings::AAString(x))
          
          # Convers 'seqs' list into an AAStringSet object
          seqs <- Biostrings::AAStringSet(seqs)
        }
        
        # Make the multiple sequence alignment with the 'msa::msa()' function from the msa package (and hide in-function printed messages)
        base::invisible(utils::capture.output(msa <- msa::msa(seqs)))
      }
      
      # If a pairwise distance matrix needs to be computed, while the 'dna.model' parameter is specified, the distance matrix is computed with the 'ape::dist.dna()' function using the 'msa' as input
      if(construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && !all(is.na(dna.model))){dist_matrix <- ape::dist.dna(ape::as.DNAbin(msa), model = dna.model, as.matrix = TRUE)}
      
      # If a pairwise distance matrix needs to be computed, while the 'aa.model' parameter is specified, the distance matrix is computed with the 'phangorn::dist.ml()' function using the 'seqs' list as input
      if(construction.method %in% c("phylo.network.default", "phylo.network.mst", "phylo.tree.nj") && !all(is.na(aa.model))){dist_matrix <- as.matrix(phangorn::dist.ml(phangorn::as.phyDat(msa, type = "AA"), model = aa.model))}
      
      # If a distance matrix is computed, sort the columns and rows to enable proper summation later on
      if(exists("dist_matrix")){dist_matrix <- dist_matrix[names(nodes_list), names(nodes_list)]}
      
      # If a multiple sequence alignment is created, sort the sequences to enable proper merge later on
      if(exists("msa")){msa@unmasked <- msa@unmasked[names(nodes_list)]}
      
      # If the 'dist_matrix' or 'msa' is not created, these objects are set to NA
      if(!exists("dist_matrix")){dist_matrix <- NA}
      if(!exists("msa")){msa <- NA}
      
      # Return the distance matrix and 
      return(list(dist_matrix = dist_matrix, msa = msa))
    })
    
    # Rename the items in the 'dist_msa' according to the column names in 'sequence_columns'
    names(dist_msa) <- sequence.columns
    
    # Retrieve the distance matrices and the multiple sequence alignments from the 'dist_msa' list and store in separate objects
    dist_matrices <- lapply(sequence.columns, function(x) dist_msa[[x]][["dist_matrix"]]); names(dist_matrices) <- sequence.columns
    multiple_sequence_alignments <- lapply(sequence.columns, function(x) dist_msa[[x]][["msa"]]); names(multiple_sequence_alignments) <- sequence.columns
    
    # Summarize the distance matrices into one matrix, which will be the input for the 'build_lineage_tree()' function next
    if(all(!is.na(dist_matrices))){
      if(length(dist_matrices) == 1){input_dist_matrix <- dist_matrices[[sequence.columns]]}
      if(length(dist_matrices) > 1){input_dist_matrix <- base::Reduce(`+`, dist_matrices)}
      }else{input_dist_matrix <- NA}
    
    # Convert the multiple sequence alignments into one phyDat object, which will be the input for the 'build_lineage_tree()' function next (and hide in-function printed messages)
    if(all(!is.na(multiple_sequence_alignments))){
      if(length(multiple_sequence_alignments) == 1){input_msa <- phangorn::as.phyDat(multiple_sequence_alignments[[sequence.columns]], type = sequence.type)}
      if(length(multiple_sequence_alignments) > 1){
        if(sequence.type == "DNA"){concatenated_seqs <- lapply(names(nodes_list), function(x) Biostrings::DNAString(do.call(paste0, lapply(sequence_columns, function(y) nodes_list[[x]][y])))); names(concatenated_seqs) <- names(nodes_list); base::invisible(utils::capture.output(input_msa <- msa::msa(Biostrings::DNAStringSet(concatenated_seqs))))}
        if(sequence.type == "AA"){concatenated_seqs <- lapply(names(nodes_list), function(x) Biostrings::AAString(do.call(paste0, lapply(sequence_columns, function(y) nodes_list[[x]][y])))); names(concatenated_seqs) <- names(nodes_list); base::invisible(utils::capture.output(input_msa <- msa::msa(Biostrings::AAStringSet(concatenated_seqs))))}
        input_msa <- phangorn::as.phyDat(input_msa, type = sequence.type)}
      }else{input_msa <- NA}
    
    
    # 3. Build lineage tree using the objects 'input_dist_matrix' and 'input_msa' 
    message(paste0("Inferring lineage tree of ", strsplit(clone, split = "_")[[1]][2], " of ", strsplit(clone, split = "_")[[1]][1], "."))
    lineage_tree <- build_lineage_tree(clone = clone,
                                       dist.matrix = input_dist_matrix,
                                       msa = input_msa,
                                       sequence.type = sequence.type,
                                       approach = approach,
                                       algorithm = algorithm,
                                       IgPhyML.trees = IgPhyML_trees,
                                       dna.model = dna.model,
                                       aa.model = aa.model,
                                       codon.model = codon.model,
                                       resolve.ties = resolve.ties,
                                       remove.internal.nodes = remove.internal.nodes,
                                       nodes.list = nodes_list)
    
    
    # 4. Create and return an output list containing the requested objects specified in the 'include' parameter
    
    # Create empty list to store output objects
    output_objects <- list()
    
    # For each possible option, check whether it is present in the 'include' vector, and if so, append it to the 'output_list'
    if("nodes" %in% include){output_objects["nodes"] <- list(nodes_list)}
    if("dist" %in% include){output_objects["distance.matrices"] <- list(dist_matrices)}
    if("msa" %in% include){output_objects["multiple.sequence.alignments"] <- list(multiple_sequence_alignments)}
    if("phylo" %in% include){output_objects["phylo"] <- lineage_tree["phylo"]}
    if("igraph" %in% include){output_objects["igraph"] <- lineage_tree["igraph"]}
    if("igraph.with.inner.nodes" %in% include){output_objects["igraph.with.inner.nodes"] <- lineage_tree["igraph.with.inner.nodes"]}
    if("edges" %in% include){output_objects["edges"] <- lineage_tree["edges"]}
    if("edges.with.inner.nodes" %in% include){output_objects["edges.with.inner.nodes"] <- lineage_tree["edges.with.inner.nodes"]}
    if("metrics" %in% include){output_objects["metrics"] <- lineage_tree["metrics"]}
    
    # Create an output list that contains the output objects and the warnings in two separate sublists
    output_list <- list(output.objects = output_objects, warnings = lineage_tree["warnings"])
    
    # Return the 'output_list'
    return(output_list)
  }
  
  # Create list with sample IDs
  sample_list <- unique(VDJ$sample_id)
  
  # Create list with clonotype IDs 
  clone_list <- unique(paste(VDJ$sample_id, VDJ$clonotype_id, sep="_"))
  
  # If an IgPhyML output file is provided, read in the IgPhyML output file and store in the 'igphyml_trees' list
  if(!missing(IgPhyML.output.file)){IgPhyML_trees <- alakazam::readIgphyml(IgPhyML.output.file)}else{igphyml_trees <- NA}
  
  # Define partial function for to be executed for each clone
  partial_function <- function(clone){
    
    # Create a temporary directory for the current clone
    temp_dir <- tempdir()
    
    # Set this 'temp_dir' as the working directory
    setwd(temp_dir)
    
    # Execute the 'infer_single_network' function for the current clone
    single_network <- infer_single_network(VDJ = VDJ,
                                           clone = clone,
                                           sequence.columns = sequence.columns,
                                           germline.columns = germline.columns,
                                           sequence.type = sequence_type,
                                           concatenate.sequences = concatenate.sequences,
                                           node.features = node.features,
                                           construction.method = construction.method,
                                           IgPhyML.trees = IgPhyML_trees,
                                           string.dist.metric = string.dist.metric,
                                           dna.model = dna.model,
                                           aa.model = aa.model,
                                           codon.model = codon.model,
                                           resolve.ties = resolve.ties,
                                           remove.internal.nodes = remove.internal.nodes,
                                           include = include)
    
    # Return the output object
    return(single_network)
  }
  
  # If 'parallel' is set to TRUE, the network inference is parallelized across the clones
  if(parallel){
    
    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]
    
    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::mclapply(mc.cores = num.cores, X = clone_list, FUN = partial_function)
      # Name items in 'output_list' according to the clonotype IDs in 'clone_list'
      names(output_list) <- clone_list
    }
    
    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      # Infer network for each clone 'clone_list' and store output in the list 'output_list'
      output_list <- parallel::parLapply(cl = cluster, X = clone_list, fun = partial_function)
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
    # Retrieve the output objects from the 'output_list' of one sample and store them in the 'sample_sublist' 
    sample_sublist <- lapply(matching_sublists, function(clone) output_list[[clone]][["output.objects"]])
    # Rename the networks in 'sample_sublist' to their original clonotype ID
    names(sample_sublist) <- sub(pattern = paste0("^", sample, "_"), replacement = "", matching_sublists)
    # Return the 'sample_sublist'
    return(sample_sublist)
  })
  
  # Reset the working directory
  setwd(original_working_directory)
  
  # Rename the sublists in 'reorganized_output_list' to their original sample ID
  names(reorganized_output_list) <- sample_list
  
  # Convert 'reorganized_output_list' of class 'list' into object of class 'AntibodyForests'
  AntibodyForests_object <- base::structure(reorganized_output_list, class = "AntibodyForests")
  
  # Retrieve the warnings from the 'output_list' and print them
  warnings <- unlist(lapply(names(output_list), function(x) output_list[[x]][["warnings"]]))
  for(i in warnings){message(paste(c("WARNING: ", i, "\n"), collapse = ""))}
  
  # Return the 'reorganized_output_list'
  return(AntibodyForests_object)
}