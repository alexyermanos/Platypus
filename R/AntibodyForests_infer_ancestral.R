#' Creates phylogenetic trees, infers ancestral sequences, and converts the resulting trees into igraph objects.

#'@description Phylogenetic trees and ancestral sequence reconstruction is performed using the IQ-TREE software. The IQ-TREE directory is required beforehand.
#' @param trees AntibodyForests object/list of AntibodyForests objects - the resulting sequence similarity or minimum spanning tree networks from the AntibodyForests function.
#' @param alignment.method string - method/software to perform multiple sequence alignment before the ancestral sequence reconstruction step. Options include: 'mafft' (requires the MAFFT software to be locally installed beforehand), 'clustal', 'clustalomega', 'tcoffee', 'muscle', which all require the 'ape' R package.
#' @param iqtree.directory string - path to the IQ-TREE software directory.
#' @param collapse.trees boolean - if T, will collapse the resulting phylogenetic trees if an intermediate daughter sequence/node is the same as its parent.
#' @param parallel boolean - whether to execute the main subroutine in parallel or not. Requires the 'parallel' R package to be installed.
#' @return nested list of AntibodyForests objects or single AntibodyForests object, with a modified tree slot including the phylogenetic tree converted into igraph objects and the reconstructed intermediate/ancestral sequences.
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests_infer_ancestral(trees, alignment.method = 'mafft',
#' igtree.directoty = '/Users/.../Desktop/iqtree-1.6.12-MacOSX')
#'}




AntibodyForests_infer_ancestral <- function(trees,
                                            alignment.method,
                                            iqtree.directory,
                                            collapse.trees,
                                            parallel){


    if(missing(trees)) stop('Please input a nested list of AntibodyForests objects - the output of the AntibodyForests function')
    if(missing(alignment.method)) alignment.method <- 'mafft' #or muscle, clustal, clustalomerga, tcoffee
    if(missing(iqtree.directory)) stop('Please input the full path to the IQ-TREE software')
    if(missing(collapse.trees)) collapse.trees <- F
    if(missing(parallel)) parallel <- F

    dir <- './tempdir'
    if(!dir.exists(dir)) dir.create(dir)

    align_call_parse_iqtree <- function(tree){
      sequences <- tree@sequences
      names(sequences) <- paste0('sequence',1:length(sequences))
      germline <- tree@germline_sequence

      if(!is.null(germline)){
        names(germline) <- 'germline'
      }

      sequences <- c(sequences, germline)
      #names(sequences) <- 1:length(sequences)

      dir_name <- paste0(
                      paste0(unique(tree@sample_id), collapse = ';'),
                      '-',
                      paste0(unique(tree@clonotype_id), collapse = ';')
                    )

      temp_dir <- paste0(dir, '/', dir_name)
      if(!dir.exists(temp_dir)) dir.create(temp_dir)

      fasta_file <- file.path(temp_dir, 'sequences.fasta')
      fasta_aligned <- file.path(temp_dir, 'sequences_aligned.fasta')


      seqinr::write.fasta(as.list(sequences), names=1:length(sequences), file.out = fasta_file)

      if(alignment.method == 'mafft'){
        system(paste0('mafft ', fasta_file, ' > ', fasta_aligned))
        alignment <- ape::read.dna(fasta_aligned, format = 'fasta')


      }else if(alignment.method == 'clustal'){
        alignment <- ape::clustal(fasta_file)

      }else if (alignment.method == 'clustalomega'){
        alignment <- ape::clustalomega(fasta_file)

      }else if(alignment.method == 'tcoffee'){
        alignment <- ape::tcoffee(fasta_file)

      }else if(alignment.method == 'muscle'){
        alignment <- ape::muscle(fasta_file)

      }else{
        stop('Unrecognized alignment method!')
      }


      rownames(alignment) <- paste0(rownames(alignment), '\t')
      phangorn::write.phyDat(alignment, file.path(temp_dir, 'infile.phy'))

      system(paste0('sh -c \'cd ', temp_dir, '; ', iqtree.directory, '/bin/iqtree -s infile.phy -m JC+G -asr\''))
      phylo_tree <- ape::read.tree(file.path(temp_dir, "infile.phy.treefile"))

      inferred_sequences <- utils::read.table(file.path(temp_dir, 'infile.phy.state'),header=TRUE)
      node_ids <- unlist(unique(inferred_sequences$Node))
      node_ids <- node_ids[order(nchar(node_ids), node_ids)]
      inferred_sequences <- lapply(node_ids, function(x) paste0(inferred_sequences$State[inferred_sequences$Node == x], collapse = '' ))
      node_ids <- as.numeric(gsub('Node', '', node_ids))
      node_ids <- length(phylo_tree$tip.label) + node_ids


      #Phylo tree processing step - match node and tip ids with the edge ids
      phylo_tree$node.label <- node_ids
      to_tip_id <- as.numeric(phylo_tree$tip.label)
      names(to_tip_id) <- 1:length(to_tip_id)

      names(node_ids) <- node_ids
      to_tip_id <- c(to_tip_id, node_ids)


      phylo_tree$edge[,1] <- to_tip_id[phylo_tree$edge[,1]]
      phylo_tree$edge[,2] <- to_tip_id[phylo_tree$edge[,2]]

      #Convert to graph
      network_df <- igraph::as_data_frame(tree@tree, what='vertices')
      germline <- which(network_df$node_type == 'germline')

      if(length(germline) > 0){

        phylo_tree <- phytools::reroot(phylo_tree, node.number = germline)

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


      for(i in 1:length(inferred_sequences)){
        network_df <- rbind(network_df, NA)
        network_df$node_type[nrow(network_df)] <- 'inferred'
        network_df$clonotype_id[nrow(network_df)] <- tree@clonotype_id
        network_df$sample_id[nrow(network_df)] <- tree@sample_id
        network_df$label[nrow(network_df)] <- nrow(network_df)
        network_df$cell_number[nrow(network_df)] <- 1
        network_df$network_sequences[nrow(network_df)] <- inferred_sequences[[i]]
      }

      features <- tree@feature_names

      if(length(features)!=0){
        for(feature in features){
          network_df[feature][network_df$node_type == 'inferred',] <- 'inferred'
          network_df[paste0(feature,'_counts')][network_df$node_type == 'inferred',] <- 1
        }
      }

      edgelist <- phylo_tree$edge
      directed <- igraph::is_directed(tree@tree)
      g <- igraph::graph_from_data_frame(edgelist, directed = directed, vertices = network_df)

      igraph::E(g)$weight <- phylo_tree$edge.length


      if(collapse.trees){
        nodes_to_check <- 1:length(igraph::V(g))
        germline_node <- which(igraph::V(g)$node_type == 'germline')

        if(length(germline_node) > 0 ){
          nodes_to_check <- nodes_to_check[-germline_node]
        }

        leaf_nodes <- which(igraph::V(g)$leaf == 'yes')
        while(length(nodes_to_check) > 1){

          current_node <- nodes_to_check[1]
          if(current_node %in% leaf_nodes){
            neighb <- igraph::neighbors(g, current_node, mode = 'all')

            if(stringdist::stringdist(igraph::V(g)$network_sequences[current_node], igraph::V(g)$network_sequences[neighb], method = 'lv') == 0){
              g <- igraph::delete_vertices(g, current_node)
              leaf_nodes <- which(igraph::degree(g) == 1)
            }
          }

          nodes_to_check <- nodes_to_check[-current_node]
        }
      }


      unlink(temp_dir, recursive = T)

      return(g)
    }


    if(inherits(trees, 'list')){
      for(i in 1:length(trees)){
        inferred_list <- vector(mode = "list", length = length(trees[[i]]))

        if(parallel){
          #requireNamespace('parallel')
          cores <- parallel::detectCores()
          inferred_list <- parallel::mclapply(trees[[i]], mc.cores = cores,
                                                   FUN = function(x) {x %>% align_call_parse_iqtree()

                                                                      })

        }else{
          inferred_list <- lapply(trees[[i]], function(x) {x %>% align_call_parse_iqtree()
                                                                })
        }

        for(j in 1:length(trees[[i]])){
          trees[[i]][[j]]@tree <- inferred_list[[j]]
        }

      }

    }else if(inherits(trees, 'AntibodyForests')){
      trees@tree <- trees %>% align_call_parse_iqtree()


    }else{
      stop(paste0('Unrecognized input tree class:  ', class(trees), '. Please ensure the input tree is either an AntibodyForests object or a nested list of AntibodyForests objects (per sample, per clonotype).'))
    }

    unlink(dir, recursive = T)
    return(trees)
}
