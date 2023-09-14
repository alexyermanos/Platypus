#' Performs sequence and structural alignment similar to PyMOL's 'align' command


#' @description Sequence and structural alignment of Steropodon structures, similar to PymMOL's 'align' command.
#' Will first perform sequence alignment, then structural superposition on the aligned coordinates for n cycles ('max.cycles' parameter), removing outliers based on a cutoff ('cutoff' parameter)
#'
#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param steropodon.template Steropodon object - template structure for alignment.
#' @param fit.to.core boolean - if TRUE, will align all structures to their respective invariant core, determined using Steropodon_find_core().
#' @param sequence.structure.superpose bool - if TRUE, will perform a sequence alignment followed by an iterative structural superposition (removing outlier atoms in the fit). This is similar to the 'align' command in PyMOL.
#' @param structure.superpose bool - if TRUE, will perform a single structure superposition/ Kabsch algorithm iteration.
#' @param alignment.method string - sequence alignment method to be used when seq.struct.superpose = TRUE. Currently only MAFFT is implemented (alignment.method = 'mafft').
#' @param max.cycles integer - the maximum number of iterations (superposition followed by outlier rejection) to be done in the sequence alignment and iterative structural superposition algorithm (seq.struct.superpose = TRUE).
#' @param cutoff float - the distance cutoff at which outliers will be rejected in the sequence alignment and iterative structural superposition algorithm (seq.struct.superpose = TRUE).
#' @param parallel bool - if TRUE, will perform the structural superposition in parallel, on all available CPU cores - 1.
#'
#' @return a nested list of Steropodon structures superposes/aligned to the steropodon.template structure or the first structure in the list.
#' @export
#' @examples
#' \dontrun{
#' superposed_seq_struct <-
#' Steropodon_superposition(steropodon_igfold,
#' sequence.structure.superpose = T,
#' structure.superpose = F,
#' max.cycles = 10,
#' cutoff = 0.5)
#'}


Steropodon_superposition <- function(steropodon.object,
                                     structure,
                                     steropodon.template,
                                     fit.to.core,
                                     sequence.structure.superpose,
                                     structure.superpose,
                                     alignment.method,
                                     max.cycles,
                                     cutoff,
                                     parallel
                                    ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(structure)) structure <- 'structure'
  if(missing(steropodon.template)) steropodon.template <- NULL
  if(missing(fit.to.core)) fit.to.core <- F
  if(missing(sequence.structure.superpose)) sequence.structure.superpose <- T
  if(missing(structure.superpose)) structure.superpose <- F
  if(missing(alignment.method)) alignment.method <- 'mafft'
  if(missing(max.cycles)) max.cycles <- 10
  if(missing(cutoff)) cutoff <- 0.5
  if(missing(parallel)) parallel <- T

  switch(Sys.info()[['sysname']],
          Windows= {message("Windows system detected")
            operating.system <- "Windows"},
          Linux  = {message("Linux system detected")
            operating.system <- "Linux"},
          Darwin = {message("MAC system detected")
            operating.system <- "Darwin"})


  if(structure == 'all'){
    structure <- c('structure', 'H', 'L', 'CDRH3', 'CDRL3', 'complex', 'paratope', 'epitope', 'core')
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  seq_ids <- names(steropodon_list)

  if(is.null(steropodon.template)){
    steropodon.template <- steropodon_list[[1]]
  }

  for(struct in structure){

    if(is.null(select_structure(steropodon_list[[1]], structure = struct))){
      next
    }

    if(is.null(select_structure(steropodon.template, structure = struct))){
      stop(paste0('Could not find structure ', struct, ' in your Steropodon superpositions template'))
    }

    pdb_list <- lapply(steropodon_list, function(x) select_structure(x, structure = struct))

    if(fit.to.core){
      template_pdb <- select_structure(steropodon.template, structure = 'core')
      if(is.null(template_pdb)){
        stop('Could not find the core structure in your Steropodon object. Please call Steropodon_find_core before!')
      }
      #Trim to C-ALPHA as the core contains only C-ALPHA atoms
      pdb_list <- lapply(pdb_list, function(pdb) {pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = pdb$core.resno))
                                                  pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, elety = 'CA'))
                                                  return(pdb)}
                                                )
    }else{
      template_pdb <- select_structure(steropodon.template, structure = struct)
    }

    if(sequence.structure.superpose){
      partial_function <- function(x) {sequence_structure_superpose(x,
                                                                    fixed = template_pdb,
                                                                    aln.method = alignment.method,
                                                                    max.cycles = max.cycles,
                                                                    cutoff = cutoff,
                                                                    return.mobile = T
                                                                    )}
      if(parallel){
        cores <- parallel::detectCores() - 1


        if(operating.system %in% c('Linux', 'Darwin')){
          pdb_list <- parallel::mclapply(pdb_list, partial_function, mc.cores = cores)

        }else{
          cl <- parallel::makeCluster(cores)
          pdb_list <- parallel::parLapply(cl, pdb_list, partial_function)
          parallel::stopCluster(cl)
        }

      }else{
        pdb_list <- lapply(pdb_list, function(x) partial_function(x))
      }
    }

    if(structure.superpose){
      partial_function <- function(x) {structure_superpose(mobile = x,
                                                           fixed = template_pdb
                                                          )}
      if(parallel){
        cores <- parallel::detectCores() - 1


        if(operating.system %in% c('Linux', 'Darwin')){
          pdb_list <- parallel::mclapply(pdb_list, partial_function, mc.cores = cores)

        }else{
          cl <- parallel::makeCluster(cores)
          pdb_list <- parallel::parLapply(cl, pdb_list, partial_function)
          parallel::stopCluster(cl)
        }

      }else{
        pdb_list <- lapply(pdb_list, function(x) partial_function(x))
      }
    }

    steropodon_list <- mapply(function(x, y) modify_structure(x, y, structure = struct), steropodon_list, pdb_list)
  }

  names(steropodon_list) <- seq_ids
  steropodon_object <- nest_steropodon(steropodon_list)
  return(steropodon_object)
}
