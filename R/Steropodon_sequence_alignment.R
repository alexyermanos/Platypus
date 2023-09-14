#' Sequence alignment of Steropodon structures


#' @description Sequence alignment of Steropodon structures. Will add a new 'pdbs' slot in each Steropodon object with the aligned sequences.
#'
#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param alignment.method string - sequence alignment algorithm to be used. Currently, only MAFFT is supported (alignment.method = 'mafft').
#' @param fit.to.template boolean - if TRUE, will fit all structures to the template given in steropodon.template. Else, all structures will be aligned to the first one.
#' @param steropodon.template Steropodon object - template structure for alignment.
#' @param annotate.features string - which features already present in the Steropodon object should be added to the new 'pdbs' object with the alignment. Defaults to 'all' to assign all present annotations in the 'pdbs' object.

#' @return a nested list of Steropodon structures superposes/aligned to the steropodon.template structure or the first structure in the list. Alignments are added in the 'pdbs' slot of the Steropodon object.
#' @export
#' @examples
#' \dontrun{
#' superposed_w_pdbs <-
#' superposed_seq_struct %>%
#' Steropodon_sequence_alignment(structure = 'structure')
#'}

Steropodon_sequence_alignment <- function(steropodon.object,
                                          structure,
                                          alignment.method,
                                          fit.to.template,
                                          steropodon.template,
                                          annotate.features
                                         ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object')
  if(missing(structure)) structure <- 'structure'
  if(missing(alignment.method)) alignment.method <- 'mafft'
  if(missing(fit.to.template)) fit.to.template <- F
  if(missing(steropodon.template)) steropodon.template <- NULL
  if(missing(annotate.features)) annotate.features <- 'all'


  sequence_alignment <- function(structure.list,
                                 alignment.method
                                 ){

    sequences <- lapply(structure.list, function(x) paste0(bio3d::pdbseq(x), collapse = ''))
    seqinr::write.fasta(as.list(sequences), names = names(sequences), file.out='temp.fasta')

    if(alignment.method == 'mafft'){
      system(paste0('mafft ', 'temp.fasta', ' > ', 'aligned.fasta'))
    }else{
      stop('Method not implemented yet!')
    }

    aln <- bio3d::read.fasta('aligned.fasta')
    unlink('temp.fasta')
    unlink('aligned.fasta')

    pdbs <- bio3d::read.fasta.pdb(aln = aln, pdblist = structure.list)

    return(pdbs)
  }

  annotate_aligned_pdbs <- function(steropodon.list,
                                    structure,
                                    columns,
                                    alignment.method
                                   ){

    not_included <- c('type', 'eleno', 'alt', 'resid', 'chain', 'resno', 'insert', 'x', 'y', 'z', 'o')
    structure_list <- lapply(steropodon.list, function(x) select_structure(x, structure = structure))
    atom_dfs <- lapply(structure_list, function(x) x$atom[x$atom$elety == 'CA',])

    if(any(columns == 'all')){
      cols <- colnames(atom_dfs[[1]])
      cols <- cols[!(cols %in% not_included)]
    }else{
      cols <- columns
    }

    pdbs <- select_structure(steropodon.list[[1]], structure = 'pdbs')

    if(is.null(pdbs)){
      pdbs <- sequence_alignment(structure.list = structure_list,
                                 alignment.method = alignment.method
                                )
    }

    resno_matrix <- pdbs$resno
    for(col in cols){
      pdbs[[col]] <- resno_matrix

      for(i in 1:length(atom_dfs)){
        vals <- atom_dfs[[i]][[col]][atom_dfs[[i]]$resno %in% resno_matrix[i,]]
        pdbs[[col]][i,][!is.na(pdbs[[col]][i,])] <- vals
      }
    }

    return(pdbs)
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  structure_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))


  pdbs <- sequence_alignment(structure.list = structure_list,
                             alignment.method = alignment.method
                            )

  steropodon_list <- lapply(steropodon_list, function(x) modify_structure(steropodon.object = x, pdb = pdbs, structure = 'pdbs'))


  if(!is.null(annotate.features)){
    pdbs <- annotate_aligned_pdbs(steropodon.list = steropodon_list,
                                  structure = structure,
                                  columns = annotate.features,
                                  alignment.method = alignment.method
                                 )
  }

  if(fit.to.template){
    if(!is.null(steropodon.template)){
      template_id <- paste0(steropodon.template@structure_id, collapse = '.')
    }else{
      template_id <- paste0(steropodon_list[[1]]@structure_id, collapse = '.')
    }

    #Testing for sparse xyz (e.g., if CDRH3)
    #inds <- bio3d::gap.inspect(pdbs$xyz)$f.inds
    inds <- NULL

    pdbs$xyz <- invisible(bio3d::fit.xyz(fixed = pdbs$xyz[template_id,],
                                         mobile=pdbs,
                                         fixed.inds = inds,
                                         mobile.inds = inds
                                         ))
  }

  steropodon_list <- lapply(steropodon_list, function(x) modify_structure(steropodon.object = x, pdb = pdbs, structure = 'pdbs'))
  steropodon_object <- nest_steropodon(steropodon_list)

  return(steropodon_object)
}
