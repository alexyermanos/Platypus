#DONE
Steropodon_find_core <- function(steropodon.object,
                                 structure,
                                 core.per,
                                 alignment.method,
                                 fit.to.core,
                                 core.volume){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(structure)) structure <- 'structure'
  if(missing(core.per)) core.per <- 'all' #or sample, or clonotype
  if(missing(alignment.method)) alignment.method <- 'mafft'
  if(missing(fit.to.core)) fit.to.core <- T
  if(missing(core.volume)) core.volume <- 0.5


  find_core <- function(structure.list,
                        alignment.method,
                        fit.to.core,
                        core.volume){

    #if(is.null(structure.list[[1]]@pdbs)){
    #
    #}

    gaps <- bio3d::gap.inspect(pdbs$xyz)
    core <- bio3d::core.find(pdbs)
    core.inds <- print(core, vol = core.volume)

    pdbs <- sequence_alignment(structure.list,
                               alignment.method = alignment.method
                               )

    #Will add pdbs and core inds objects to the Steropodon object, as well as trimmed cores (from first struct) and resnos found in the core
    pdbs$core.inds <- core.inds
    pdbs$no.gaps.inds <- gaps$f.inds
    pdbs$core <- core

    if(fit.to.core){
      xyz <- bio3d::pdbfit(pdbs, core.inds)
      pdbs$xyz <- xyz
    }

    ids <- which(pdbs$resno[1,] %in% core.inds$resno)
    for(i in 1:length(structure.list)){
      resno <- pdbs$resno[i, ids]
      structure.list[[i]]$core.resno <- resno
      structure.list[[i]]$atom$core <- NA
      structure.list[[i]]$atom$core[structure.list[[i]]$atom$resno %in% resno] <- 'core'
    }

    resno <- pdbs$resno[1, ids] #or structure.list[[1]]$atom$resno %in% core.inds$resno
    core <- bio3d::trim(structure.list[[1]], inds = bio3d::atom.select(structure.list[[1]], elety = 'CA'))
    core <- bio3d::trim(core, inds = bio3d::atom.select(core, resno = resno))

    return(list(core = core, pdbs = pdbs, structure.list = structure.list))
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  seq_ids <- names(steropodon_list)
  structure_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))

  if(core.per == 'all'){
    core_list <- find_core(structure.list = structure_list,
                           alignment.method = alignment.method,
                           fit.to.core = fit.to.core,
                           core.volume = core.volume)

    core <- core_list$core
    pdbs <- core_list$pdbs
    pdb_list <- core_list$structure.list

    steropodon_list <- mapply(function(x,y) modify_structure(steropodon.object = x, pdb = y, structure = structure), steropodon_list, pdb_list)
    steropodon_list <- lapply(steropodon_list, function(x) modify_structure(steropodon.object = x, pdb = core, structure = 'core'))
    steropodon_list <- lapply(steropodon_list, function(x) modify_structure(steropodon.object = x, pdb = pdbs, structure = 'pdbs'))


  }else if(core.per == 'sample'){
    sample_clonotype_list <- lapply(seq_ids, function(x) unlist(stringr::str_split(x, '\\.')))
    sample_list <- lapply(sample_clonotype_list, function(x) x[[1]])

    unique_samples <- unique(sample_list)
    for(sample in unique_samples){
      inds <- which(sample_list == sample)

      subset_structures <- structure_list[inds]

      core_list <- find_core(structure.list = subset_structures,
                             alignment.method = alignment.method,
                             fit.to.core = fit.to.core,
                             core.volume = core.volume)

      core <- core_list$core
      pdbs <- core_list$pdbs
      pdb_list <- core_list$structure.list

      subset_structures <- mapply(function(x,y) modify_structure(steropodon.object = x, pdb = y, structure = structure), subset_structures, pdb_list)
      subset_structures <- lapply(subset_structures, function(x) modify_structure(steropodon.object = x, pdb = core, structure = 'core'))
      subset_structures <- lapply(subset_structures, function(x) modify_structure(steropodon.object = x, pdb = pdbs, structure = 'pdbs'))

      sample_list[inds] <- subset_structures
    }

  }else if(core.per == 'clonotype'){
    sample_clonotype_list <- lapply(seq_ids, function(x) paste0(unlist(stringr::str_split(x, '\\.'))[1:2], collapse = '.'))

    unique_samples <- unique(sample_clonotype_list)
    for(sample in unique_samples){
      inds <- which(sample_clonotype_list == sample)

      subset_structures <- structure_list[inds]

      core_list <- find_core(structure.list = subset_structures,
                             alignment.method = alignment.method,
                             fit.to.core = fit.to.core,
                             core.volume = core.volume)

      core <- core_list$core
      pdbs <- core_list$pdbs
      pdb_list <- core_list$structure.list

      subset_structures <- mapply(function(x,y) modify_structure(steropodon.object = x, pdb = y, structure = structure), subset_structures, pdb_list)
      subset_structures <- lapply(subset_structures, function(x) modify_structure(steropodon.object = x, pdb = pdbs, structure = 'pdbs'))

      sample_list[inds] <- subset_structures
    }

  }else{
    stop('Please choose whether to find the invariant core per sample, clonotype, or across all structures!')
  }

  names(steropodon_list) <- seq_ids

  steropodon_object <- nest_steropodon(steropodon_list)

  return(steropodon_object)
}
