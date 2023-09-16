#' Trims the structures in a Steropodon object (removes chains, framework/variable regions, invariant core regions)

#' @description Trims the structures in a Steropodon object (removes chains, framework/variable regions, invariant core regions).
#' The region to be removed must be included in the 'grouping' parameter (e.g., 'chain') and the specific values/labels to be removed in the specific.values (e.g., 'VDJ', removing the heavy chain in a given structure).
#' The combine.values parameter determines whether the regions inserted in specific.values should all be combined in the same final structure.
#' For example, to obtain a new structure with only VDJ and VJ CDR regions we can set combine.values to T and combine.grouping to T (which will paste together all regions if the grouping parameter is a vector of structural region/feature columns).
#' This is useful for calculating distances only across the CDR regions using Steropodon_distances.
#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param grouping string or vector of strings - the feature name in a Steropodon object that includes the parts we want trimmed (e.g., 'chain' if we want to remove a VDJ or VJ chain).
#' @param specific.values string or vector of strings - the specific regions to kept in the structure, depending on the 'grouping' parameter (e.g., grouping = 'chain' and specific.values = 'VDJ' will remove all keep chains).
#' @param combine.values bool - if TRUE, the regions inserted in specific.values will all be combined in the same final structure.
#' @param combine.groupings bool - if TRUE, will combine the 'grouping' features for a more specific trimming (e.g., 'chain' and 'region' to select exactly which hypervariable regions from any chain should be kept - 'VDJ_CDR1', 'VJ_CDR1' etc.).
#'
#' @return a nested list of Steropodon structures with all regions not included in 'specific.values' removed. For more information, see the Steropodon vignette (https://alexyermanos.github.io/Platypus/articles/Steropodon.html).
#' @export
#' @examples
#' \dontrun{
#'steropodon_igfold$s1$clonotype1$`1` %>%
#'  Steropodon_trim(structure = 'structure',
#'                  grouping = c('chain', 'region'),
#'                  specific.values = c('VDJ_CDR1','VDJ_CDR2',
#'                  VDJ_CDR3','VJ_CDR1','VJ_CDR2','VJ_CDR3'),
#'                  combine.values = T,
#'                  combine.groupings = T) %>%
#'  Steropodon_visualize(structure = 'structure',
#'                       color.by = c('chain', 'region'))
#'}

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

    combined <- NULL

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
