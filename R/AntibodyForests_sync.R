#' Function to synchronize the node labels/names of all clonotypes within all samples of two AntibodyForests objects.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description The nodes of each clonotype within each sample of the subject AntibodyForests object will be named according to the names of the nodes of the clonotypes within the samples of the reference AntibodyForests object. The node names present in all the objects within the  Therefore, the sample IDs and clonotype IDs should be the same. Note: if a node in the reference AntibodyForests object is divided over two nodes in the subject AntibodyForests object, the nodes will get a letter as suffix (for example, 'node2' in the reference object would become 'node2A' and 'node2B' in the subject object). Note: if multiple nodes in the reference AntibodyForests object are together in one node in the subject AntibodyForests object, the number of the nodes are pasted together with a '+' (for example, 'node5' and 'node6' in the reference object would become 'node5+6' in the subject object). 
#' @param reference AntibodyForests object - AntibodyForests object as obtained from the 'AntibodyForests()' function in Platypus. This object will be used as a reference.
#' @param subject AntibodyForests object - AntibodyForests object as obtained from the 'AntibodyForests()' function in Platypus. For each clonotype, the names of the nodes will be synced with the names of the nodes in the reference AntibodyForests object, by matching the barcodes.  
#' @return Returns the subject AntibdoyForests object in which all nodes of each clonotypes within all samples are renamed.
#' @examples
#' \dontrun{
#' AntibodyForests_sync(reference = AntibodyForests_object1, 
#'                      subject = AntibodyForests_object2)
#'}

AntibodyForests_sync <- function(reference,
                                 subject){
  
  # Check whether both a reference AntibodyForests object and a subject AntibodyForests object is provided
  if(missing(reference) | missing(subject)){stop("ERROR: Please provide both a reference AntibodyForests object and a subject AntibodyForests object.")}
  
  
  rename_nodes <- function(reference,
                           subject,
                           clone){
    
    # Renames the nodes of a single clone in the 'subject' AntibodyForests object using the 'reference' AntibodyForests object
    # Arguments:
    # - reference: AntibodyForests object that will be used as 'reference' and with which the subject AntibodyForests object will be synchronized
    # - subject: AntibodyForests object of which the nodes of all clonotypes will be renamed
    # - clone: string specifying the sample and clonotype in the format 'S1_clonotype1'
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # Retrieve sample ID and clonotype ID from 'clone'
    sample <- strsplit(clone, split="_")[[1]][1]
    clonotype <- strsplit(clone, split="_")[[1]][2]
    
    # Retrieve the names of the nodes of the clonotype in the reference and subject AntibodyForests object, and store them in the 'reference_names' and 'subject_names' vector, respectively
    reference_names <- names(reference[[sample]][[clonotype]][["nodes"]])
    subject_names <- names(subject[[sample]][[clonotype]][["nodes"]])
    
    # Remove the 'germline' node from both the 'reference_names' and 'subject_names' vector
    reference_names <- reference_names[reference_names != "germline"]
    subject_names <- subject_names[subject_names != "germline"]
    
    # For each node in the reference and subject AntibodyForests object, retrieve the barcodes, and store them in the 'reference_barcodes' and 'subject_barcodes' vector, respectively
    reference_barcodes <- lapply(reference_names, function(x) reference[[sample]][[clonotype]][["nodes"]][[x]][["barcodes"]]); names(reference_barcodes) <- reference_names
    subject_barcodes <- lapply(subject_names, function(x) subject[[sample]][[clonotype]][["nodes"]][[x]][["barcodes"]]); names(subject_barcodes) <- subject_names
    
    # Create a list with the names of the subject nodes as indices and the names of the refrence nodes as items
    reference_subject_matching <- lapply(subject_names, function(x){
      # Get the barcodes of the current subject node
      barcodes <- subject_barcodes[[x]]
      # Retrieve the names of the reference nodes that contain at least one barcode present in the 'barcodes' list
      matching_reference_nodes <- unlist(lapply(reference_names, function(y){
        if(sum(reference_barcodes[[y]] %in% barcodes) >= 1){return(y)}
      }))
      # Return the matching reference nodes
      return(matching_reference_nodes)
    })
    names(reference_subject_matching) <- subject_names
    
    # Create a list 'new_names' in which the indices refer to the subject names and the items refer to reference names (later on, this list will be used to rename the subject nodes)
    new_names <- sapply(subject_names, function(x){
      # Retrieve the the reference name(s) for the current subject name
      new_name <- reference_subject_matching[[x]]
      # If multiple reference names are found, paste them together, separated by a '+'
      if(length(new_name) > 1){new_name <- paste0("node", paste(sort(as.numeric(gsub(pattern = "node", replacement = "", new_name))), collapse = "+"))}
      # Return the reference name
      return(new_name)
    })
    
    # Retrieve the names in the 'new_names' vector that occur more than once
    duplicates <- names(table(new_names))[table(new_names) > 1]
    
    # Rename the duplicate names by pasting a letter (A-Z) to the node name
    for(i in duplicates){new_names[new_names == i] <- paste0(i, LETTERS[1:table(new_names)[[i]]])}
    
    # Append the germline node to the 'new_names' list
    new_names["germline"] <- "germline"
    
    # Rename the names of the nodes in all the present objects in the clonotype list of the subject AntibodyForests object using the 'new_names' vector
    if(!is.null(subject[[sample]][[clonotype]][["nodes"]])){names(subject[[sample]][[clonotype]][["nodes"]]) <- new_names[names(subject[[sample]][[clonotype]][["nodes"]])]}
    if(!is.null(subject[[sample]][[clonotype]][["igraph"]])){igraph::V(subject[[sample]][[clonotype]][["igraph"]])$name[igraph::V(subject[[sample]][[clonotype]][["igraph"]])$name %in% names(new_names)] <- new_names[igraph::V(subject[[sample]][[clonotype]][["igraph"]])$name[igraph::V(subject[[sample]][[clonotype]][["igraph"]])$name %in% names(new_names)]]}
    if(!is.null(subject[[sample]][[clonotype]][["igraph.with.inner.nodes"]])){igraph::V(subject[[sample]][[clonotype]][["igraph.with.inner.nodes"]])$name[igraph::V(subject[[sample]][[clonotype]][["igraph.with.inner.nodes"]])$name %in% names(new_names)] <- new_names[igraph::V(subject[[sample]][[clonotype]][["igraph.with.inner.nodes"]])$name[igraph::V(subject[[sample]][[clonotype]][["igraph.with.inner.nodes"]])$name %in% names(new_names)]]}
    if(!is.null(subject[[sample]][[clonotype]][["edges"]])){subject[[sample]][[clonotype]][["edges"]][["upper.node"]][subject[[sample]][[clonotype]][["edges"]][["upper.node"]] %in% names(new_names)] <- new_names[subject[[sample]][[clonotype]][["edges"]][["upper.node"]][subject[[sample]][[clonotype]][["edges"]][["upper.node"]] %in% names(new_names)]]; subject[[sample]][[clonotype]][["edges"]][["lower.node"]][subject[[sample]][[clonotype]][["edges"]][["lower.node"]] %in% names(new_names)] <- new_names[subject[[sample]][[clonotype]][["edges"]][["lower.node"]][subject[[sample]][[clonotype]][["edges"]][["lower.node"]] %in% names(new_names)]]}
    if(!is.null(subject[[sample]][[clonotype]][["edges.with.inner.nodes"]])){subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["upper.node"]][subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["upper.node"]] %in% names(new_names)] <- new_names[subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["upper.node"]][subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["upper.node"]] %in% names(new_names)]]; subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["lower.node"]][subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["lower.node"]] %in% names(new_names)] <- new_names[subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["lower.node"]][subject[[sample]][[clonotype]][["edges.with.inner.nodes"]][["lower.node"]] %in% names(new_names)]]}
    if(!is.null(subject[[sample]][[clonotype]][["distance.matrices"]])){for(i in names(subject[[sample]][[clonotype]][["distance.matrices"]])){rownames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]])[rownames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]]) %in% names(new_names)] <- new_names[rownames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]])[rownames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]]) %in% names(new_names)]]; colnames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]])[colnames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]]) %in% names(new_names)] <- new_names[colnames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]])[colnames(subject[[sample]][[clonotype]][["distance.matrices"]][[i]]) %in% names(new_names)]]}}
    if(!is.null(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]])){for(i in names(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]])){names(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]][[i]]@unmasked)[names(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]][[i]]@unmasked) %in% names(new_names)] <- new_names[names(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]][[i]]@unmasked)[names(subject[[sample]][[clonotype]][["multiple.sequence.alignments"]][[i]]@unmasked) %in% names(new_names)]]}}
    
    # Return the subject AntibodyForests object in the nodes of all clonotypes are renamed
    return(subject)
  }
  
  
  # Retrieve all clonotype names from both the reference and subject AntibodyForests object, add the sample IDs as prefixes to the clonotype IDs, and store them in the 'reference_clones' and 'subjec_clones' vector, respectively
  reference_clones <- unlist(lapply(names(reference), function(x) paste(x, names(reference[[x]]), sep = "_")))
  subject_clones <- unlist(lapply(names(subject), function(x) paste(x, names(subject[[x]]), sep = "_")))
  
  # Create a vector with clones that are present in both the reference and subject AntibodyForests object
  shared_clones <- base::intersect(reference_clones, subject_clones)
  
  # Iterate through the clones in the 'shared_clones' vector and apply the 'rename_nodes()' function to each clone
  for(clone in shared_clones){subject <- rename_nodes(reference = reference, subject = subject, clone = clone)}
  
  # Return 
  return(subject)
}