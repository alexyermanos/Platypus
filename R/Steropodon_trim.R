#DONE
#TD: make it work w aligned pdbs
Steropodon_trim <- function(steropodon.object,
                            structure,
                            grouping,
                            specific.values,
                            combine.values,
                            combine.groupings
                            ){

    if(missing(steropodon.object)) stop('Please input your Steropodon object!')
    if(missing(structure)) structure <- 'structure'
    if(missing(grouping)) grouping <- c('chain', 'region')
    if(missing(specific.values)) specific.values <- c('VDJ_CDR1', 'VDJ_CDR2', 'VDJ_CDR3', 'VJ_CDR1', 'VJ_CDR2', 'VJ_CDR3')
    if(missing(combine.values)) combine.values <- T
    if(missing(combine.groupings)) combine.groupings <- T

    split_structure <- function(pdb,
                                grouping,
                                specific.values,
                                combine.groupings,
                                combine.values
                               ){

      if(missing(pdb)) stop('Please input your PDB structure')
      if(missing(grouping)) grouping <- 'chain'
      if(missing(specific.values)) specific.values <- NULL
      if(missing(combine.groupings)) combine.groupings <- F
      if(missing(combine.values)) combine.values <- F

      is_null_pdb <- function(pdb){
        out <- F
        if(nrow(pdb$atom) == 0){
          out <- T
        }
        out
      }

      if(combine.groupings & (length(grouping) > 1)){
        pdb$atom$combined <- rep('', nrow(pdb$atom))
        for(group in grouping){
          pdb$atom$combined <- paste0(pdb$atom$combined, '_', pdb$atom[[group]])
        }
        pdb$atom <- pdb$atom %>%
                    dplyr::mutate(combined = substring(combined,2))

        grouping <- 'combined'
      }

      if(is.null(specific.values)){
        specific.values <- unique(pdb$atom[[grouping]])
      }

      if(length(specific.values) > 1 & !combine.values){
        out_pdbs <- vector(mode = 'list', length = length(specific.values))
        for(i in 1:length(specific.values)){
          res_idx <- pdb$atom$resno[pdb$atom[[grouping]] %in% specific.values[i]]
          idx <- bio3d::atom.select(pdb, resno = res_idx)
          trimmed_struct <- bio3d::trim(pdb, inds = idx)
          trimmed_struct$atom$combined <- NULL
          out_pdbs[[i]] <- trimmed_struct
        }
        names(out_pdbs) <- specific.values
        out_pdbs <- out_pdbs[!sapply(out_pdbs, is_null_pdb)]

      }else{
        res_idx <- pdb$atom$resno[pdb$atom[[grouping]] %in% specific.values]
        idx <- bio3d::atom.select(pdb, resno = res_idx)
        trimmed_struct <- bio3d::trim(pdb, inds = idx)
        trimmed_struct$atom$combined <- NULL
        out_pdbs <- trimmed_struct

        if(is_null_pdb(out_pdbs)){
          return(NULL)
        }
      }

      return(out_pdbs)
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

    steropodon_list <- lapply(steropodon_list, function(x) {
      pdb <- select_structure(x, structure = structure)
      split_pdb <- split_structure(pdb = pdb,
                                   grouping = grouping,
                                   specific.values = specific.values,
                                   combine.groupings = combine.groupings,
                                   combine.values = combine.values)
      obj <- modify_structure(x, structure = structure, pdb = split_pdb)
      return(obj)
    })

    steropodon_object <- nest_steropodon(steropodon_list)
    return(steropodon_object)
}
