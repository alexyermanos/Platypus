#' Save structures modelled using Steropodon as PDB files


#' @description Function to save all structures in the Steropodon class slot specified in the 'structure' parameter as PDB files.
#'
#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param save.full boolean - if TRUE, will save all structures in a Steropodon nested list.
#' @param save.dir string - path to the directory in which all PDB files will be saved.

#' @return no returns. Saves the Steropodon object(s) as PDB files in the directory specified in save.dir.
#' @export
#' @examples
#' \dontrun{
#' steropodon_igfold %>%
#' Steropodon_save(structure = 'structure',
#'                 save.dir = 'example_saved',
#'                 save.full = T)
#'}



Steropodon_save <- function(steropodon.object,
                            structure,
                            save.full,
                            save.dir
                           ){

   if(missing(steropodon.object)) stop('Please input your Steropodon object/ nested list of objects!')
   if(missing(structure)) structure <- 'structure'
   if(missing(save.full)) save.full <- T
   if(missing(save.dir)) save.dir <- 'steropodon_saved'


   if(!dir.exists(save.dir)){
     dir.create(save.dir)
   }

   save.dir <- paste0(save.dir, '/', structure)
   if(!dir.exists(save.dir)){
     dir.create(save.dir)
   }

   write_steropodon_pdbs(steropodon.object, structure = structure, dir = save.dir)

   if(save.full){
     steropodon_list <- unnest_steropodon(steropodon.object)
     file_names <- lapply(steropodon_list, function(x) paste0(x@structure_id, collapse = '_'))

     for(i in 1:length(steropodon_list)){
       struct <- select_structure(steropodon_list[[i]], structure = structure)
       file_path <- paste0(save.dir, '/', file_names[i], '.csv')
       utils::write.csv(struct$atom, file_path)
     }
   }
}
