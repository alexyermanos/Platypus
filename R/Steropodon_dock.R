#' Antibody-antigen dockling for a Steropodon object


#' @description Performs antibody-antigen docking using ZDOCK and saves the resulting complex in the 'complex' slot of all Steropodon structures in the nested list/for a single Steropodon object.
#' The antibody-antigen ccomplex can be used for downstream analysis (e.g., obtaining interface and paratope-epitope structures and sequences using Steropodon_interface).

#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param docking.tool string - the docking tool that should be used. Currently, only ZDOCK is implemented (docking.tool = 'zdock').
#' @param tool.directory string - the path to the docking tool main directory.
#' @param antigen string - PDB ID of the target antigen that should be docked.
#' @param antigen.name string - name of the docked target.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param additional.docking.parameters named list - additional parameters for the docking tool. For ZDOCK, these include zdock.n.prediction (number of dockings to be performed), zdock.fixed.receptor (if the receptor should be fixed or flexible when docking).
#' @param parallel bool - if docking should be performed in parallel (requires multiple cores).


#' @return the Steropodon object or nested list of objects with the 'complex' slot including the docked antibody-antigen structure.
#' @export
#' @examples
#' \dontrun{
#'steropodon_docked <-
#'  steropodon_igfold$s1$clonotype1$`1` %>%
#'  Steropodon_dock(docking.tool = 'zdock',
#'                  tool.directory = '/Users/tudorcotet/Desktop/zdock',
#'                 antigen = '2tnf',
#'                 antigen.name = 'TNFR2',
#'                 structure = 'structure',
#'                 additional.docking.params = list(zdock.n.predictions = 10,
#'                                                   zdock.fixed.receptor = T),
#'                 parallel = F)
#'}



Steropodon_dock <- function(steropodon.object,
                            docking.tool,
                            tool.directory,
                            antigen,
                            antigen.name,
                            structure,
                            additional.docking.params,
                            parallel
                           ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(docking.tool)) docking.tool <- 'zdock'
  if(missing(tool.directory)) stop('Please input the ZDOCK directory path')
  if(missing(antigen)) stop('Please input either a antigen PDB ID or a path to the PDB file')
  if(missing(antigen.name)) stop('Please input the antigen name!')
  if(missing(structure)) structure <- 'structure'
  if(missing(additional.docking.params)) additional.docking.params <- list()
  if(missing(parallel)) parallel <- T

  combined <- NULL

  call_zdock <- function(steropodon.object,
                         antigen,
                         antigen.name,
                         structure,
                         zdock.directory,
                         zdock.n.predictions,
                         zdock.specific.regions.to.mask,
                         zdock.region.columns,
                         zdock.fixed.receptor,
                         parallel
                         ){

    if(missing(steropodon.object)) stop('Please input your Steropodon object!')
    if(missing(antigen)) stop('Please input either a antigen PDB ID or a path to the PDB file')
    if(missing(antigen.name)) stop('Please input the antigen name!')
    if(missing(structure)) structure <- 'structure'
    if(missing(zdock.directory)) stop('Please input the ZDOCK directory path')
    if(missing(zdock.n.predictions)) zdock.n.predictions <- 2000
    if(missing(zdock.specific.regions.to.mask)) zdock.specific.regions.to.mask <- c('VDJ_FR1', 'VDJ_FR2', 'VDJ_FR3', 'VDJ_FR4', 'VJ_FR1', 'VJ_FR2', 'VJ_FR3', 'VJ_FR4')
    if(missing(zdock.region.columns)) zdock.region.columns <- c('chain', 'region')
    if(missing(zdock.fixed.receptor)) zdock.fixed.receptor <- F
    if(missing(parallel)) parallel

    #zdock.directory <- '/Users/tudorcotet/Desktop/zdock'

    call_zdock_lapply <- function(receptor.file,
                                  antigen.file,
                                  steropodon.list,
                                  structure,
                                  zdock.directory,
                                  zdock.n.predictions,
                                  zdock.specific.regions.to.mask,
                                  zdock.region.columns,
                                  zdock.fixed.receptor
                                 ){

      if(zdock.fixed.receptor){
        fixed_receptor <- ' -F'
      }else{
        fixed_receptor <- ''
      }

      system(
        paste0(
          zdock.directory, '/mark_sur ', receptor.file, ' ', receptor.file
        )
      )
      id <- stringr::str_split(receptor.file, '\\.')[[1]][1]
      id_1 <- stringr::str_replace_all(id, '_', '\\.')
      out_file_name <- paste0(id, '_predicted.out')

      if(!is.null(zdock.specific.regions.to.mask)){
        pdb <- select_structure(steropodon.list[[id_1]], structure = structure)

        pdb$atom$combined <- rep('', nrow(pdb$atom))
        for(group in zdock.region.columns){
          pdb$atom$combined <- paste0(pdb$atom$combined, '_', pdb$atom[[group]])
        }
        pdb$atom <- pdb$atom %>%
                    dplyr::mutate(combined = substring(combined,2))

        res_list <- unique(pdb$atom$resno[pdb$atom$combined %in% zdock.specific.regions.to.mask])
        res_list <- res_list[order(res_list)]
        masked_file <- paste0(id, '.txt')

        fileConn <- file(masked_file)
        writeLines(as.character(res_list), fileConn)
        close(fileConn)

        receptor_file_masked <- paste0(id, '_masked.pdb')
        system(
          paste0(
            zdock.directory, '/block.pl ', receptor.file, ' ', masked_file, ' > ', receptor_file_masked
          )
        )

        receptor.file <- receptor_file_masked
      }

      system(
        paste0(
          zdock.directory, '/zdock', ' -R ', receptor.file, ' -L ', antigen.file, ' -o ', out_file_name, ' -N ', zdock.n.predictions, fixed_receptor
        )
      )

      return('DONE')
    }

    switch(Sys.info()[['sysname']],
            Windows= {message("Windows system detected")
              operating.system <- "Windows"},
            Linux  = {message("Linux system detected")
              operating.system <- "Linux"},
            Darwin = {message("MAC system detected")
              operating.system <- "Darwin"})


    temp_dir <- './temp_dir'
    if(!dir.exists(temp_dir)) dir.create(temp_dir)
    out_dir <- './out_dir'
    if(!dir.exists(out_dir)) dir.create(out_dir)


    if(inherits(steropodon.object, 'list')){
      steropodon_list <- unnest_steropodon(steropodon.object)
    }else{
      steropodon_list <- list()
      steropodon_list[[1]] <- steropodon.object
      names(steropodon_list) <- paste0(steropodon.object@structure_id, collapse = '.')
    }

    steropodon_names <- names(steropodon_list)

    current_dir <- getwd()
    chain_dict <- write_steropodon_pdbs(steropodon.object, structure = structure, dir = temp_dir)
    chain_dict <- stats::setNames(names(chain_dict), chain_dict)

    receptor_files <- list.files(temp_dir)

    antigen_file <- paste0(antigen.name, '.pdb')
    antigen_file_path <- paste0(temp_dir, '/', antigen_file)
    antigen_pdb <- bio3d::write.pdb(bio3d::read.pdb(antigen), antigen_file_path)

    setwd(temp_dir)

    system(
      paste0(
        'cp ', zdock.directory, '/uniCHARMM ', '.'
      )
    )

    system(
      paste0(
        zdock.directory, '/mark_sur ', antigen_file, ' ', antigen_file
      )
    )

    partial_function <- function(x) call_zdock_lapply(receptor.file = x,
                                                      antigen.file = antigen_file,
                                                      steropodon.list = steropodon_list,
                                                      structure = structure,
                                                      zdock.directory = zdock.directory,
                                                      zdock.n.predictions = zdock.n.predictions,
                                                      zdock.specific.regions.to.mask = zdock.specific.regions.to.mask,
                                                      zdock.region.columns =  zdock.region.columns,
                                                      zdock.fixed.receptor =  zdock.fixed.receptor
                                                     )

    if(parallel){
       cores <- parallel::detectCores() - 1

       if(operating.system %in% c('Linux', 'Darwin')){
         parallel::mclapply(receptor_files, partial_function, mc.cores = cores)

       }else{
         cl <- parallel::makeCluster(cores)
         parallel::parLapply(cl, receptor_files, partial_function)
         parallel::stopCluster(cl)
       }

    }else{
       lapply(receptor_files, function(x) partial_function(x))
    }


    setwd(current_dir)
    out_files <- list.files(temp_dir)
    out_files <- out_files[stringr::str_detect(out_files, '_predicted.out')]
    setwd(temp_dir)

    system(
      paste0(
        'cp ', zdock.directory, '/create_lig ', '.'
      )
    )

    #Post-processing/ trimming/ annotation/ etc.
    #Maybe parallelize?
    for(file in out_files){
      id <- stringr::str_split(file, '_')[[1]][1:3]
      id_1 <- paste0(id, collapse = '.')
      id <- paste0(id, collapse = '_')

      file_name <- paste0('../out_dir', '/', id, '.pdb')

      system(
        paste0(
          zdock.directory, '/create.pl ', file
        )
      )

      file.rename('complex.1.pdb', file_name)

      pdb <- bio3d::read.pdb(file_name)
      receptor_pdb <- select_structure(steropodon_list[[id_1]], structure = structure)

      pdb$atom$eleno <- 1:nrow(pdb$atom)
      docked_receptor_pdb <- pdb$atom[1:nrow(receptor_pdb$atom),]
      docked_antigen_pdb <- pdb$atom[(nrow(receptor_pdb$atom)+1):nrow(pdb$atom),]
      final_resno <- max(docked_receptor_pdb$resno)

      antigen_resnos <- unique(docked_antigen_pdb$resno)
      antigen_resnos <- antigen_resnos[order(antigen_resnos)]
      resno_dict <- final_resno + 1:length(antigen_resnos)
      names(resno_dict) <- antigen_resnos
      docked_antigen_pdb$resno <- resno_dict[as.character(docked_antigen_pdb$resno)]

      docked_receptor_pdb$chain <- chain_dict[docked_receptor_pdb$chain]
      docked_antigen_pdb$chain <- 'antigen'

      docked_receptor_pdb$region <- receptor_pdb$atom$region
      docked_antigen_pdb$region <- 'antigen'

      missing_features <- setdiff(colnames(receptor_pdb$atom), colnames(docked_receptor_pdb))

      for(feature in missing_features){
        docked_receptor_pdb[[feature]] <- receptor_pdb$atom[[feature]]
        docked_antigen_pdb[[feature]] <- rep(NA, nrow(docked_antigen_pdb))
      }


      pdb$atom <- rbind(docked_receptor_pdb, docked_antigen_pdb)

      antigen_structure <- split_structure(pdb, grouping = c('chain'), specific.values = 'antigen')
      receptor_structure <- split_structure(pdb, grouping = c('chain'), specific.values = c('VDJ','VJ'), combine.values = T, combine.groupings = F)

      #Maybe also modify CDRH3s and other split structures that already exist in the Steropodon object
      steropodon_list[[id_1]]@complex <- pdb
      steropodon_list[[id_1]]@antigen <- antigen_structure
      steropodon_list[[id_1]] <- modify_structure(steropodon_list[[id_1]], structure = structure, pdb = receptor_structure)
    }

    names(steropodon_list) <- steropodon_names
    steropodon_object <- nest_steropodon(steropodon_list)
    return(steropodon_object)
  }


  if(docking.tool == 'zdock'){
    params <- list(steropodon.object = steropodon.object,
                   zdock.directory = tool.directory,
                   antigen = antigen,
                   antigen.name = antigen.name,
                   structure = structure,
                   parallel = parallel
                  )

    steropodon_object <- do.call(call_zdock, c(params, additional.docking.params))


  }else{
    stop('Docking tool not implemented yet!')
  }

  return(steropodon_object)
}
