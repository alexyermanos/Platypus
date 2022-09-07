#DONE
Steropodon_interface <- function(steropodon.object,
                                 distance.threshold,
                                 include.hydrogens
                                ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object/ list of objects')
  if(missing(distance.threshold)) distance.threshold <- 5
  if(missing(include.hydrogens)) include.hydrogens <- T

  get_interface <- function(steropodon.object,
                            distance.threshold,
                            include.hydrogens
                            ){

    pdb <- select_structure(steropodon.object, structure = 'complex')

    if(is.null(pdb)){
      stop('Could not find the antibody-antigen complex structure. Please either call Steropodon model using ColabFold, using a specific antigen, or call Steropodon_dock on your main structure!')
    }

    na_residues <- unique(pdb$atom$resno[is.na(pdb$atom$chain)])
    if(length(na_residues) > 0){
      pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = na_residues))
    }

    pdb$atom$temp_chain <- pdb$atom$chain
    pdb$atom$chain[pdb$atom$chain %in% c('VDJ', 'VJ')] <- 'receptor'
    pdb$atom$chain[!(pdb$atom$chain %in% c('VDJ', 'VJ'))] <- 'ligand'

    receptor_inds <- bio3d::atom.select(pdb, chain = 'receptor')
    ligand_inds <- bio3d::atom.select(pdb, chain = 'ligand')

    epitope <- bio3d::binding.site(pdb,
                                   a.inds = ligand_inds,
                                   b.inds = receptor_inds,
                                   cutoff = distance.threshold,
                                   hydrogens = include.hydrogens)

    paratope <- bio3d::binding.site(pdb,
                                    a.inds = receptor_inds,
                                    b.inds = ligand_inds,
                                    cutoff = distance.threshold,
                                    hydrogens = include.hydrogens)

    epitope_struct <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = epitope$resno))
    paratope_struct <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = paratope$resno))

    steropodon.object@epitope <- epitope_struct
    steropodon.object@paratope <- paratope_struct

    pdb$atom$chain <- pdb$atom$temp_chain
    pdb$atom$temp_chain <- NULL

    steropodon.object@complex <- pdb

    return(steropodon.object)
  }


  if(inherits(steropodon.object, 'list')){
     steropodon_list <- unnest_steropodon(steropodon.object)
  }else if(inherits(steropodon.object, 'Steropodon')){
     steropodon_list <- list()
     steropodon_list[[1]] <- steropodon.object
     names(steropodon_list) <- paste0(steropodon.object@structure_id, collapse = '.')
  }else{
     stop('Unrecognized Steropodon object: please input either a Steropodon object or a nested list of objects, as obtained from Steropodon_model')
  }


  steropodon_list <- lapply(steropodon_list, function(x) get_interface(x, distance.threshold = distance.threshold, include.hydrogens = include.hydrogens))
  steropodon_object <- nest_steropodon(steropodon_list)

  return(steropodon_object)
}
