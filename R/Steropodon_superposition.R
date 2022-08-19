#DONE
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
          pdb_list <- parallel::parlapply(cl, pdb_list, partial_function)
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
          pdb_list <- parallel::parlapply(cl, pdb_list, partial_function)
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
