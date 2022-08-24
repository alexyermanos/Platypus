#DONE
Steropodon_dock <- function(steropodon.object,
                            docking.tool,
                            tool.directory,
                            ligand,
                            ligand.name,
                            structure,
                            additional.docking.params,
                            parallel
                           ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(docking.tool)) docking.tool <- 'zdock'
  if(missing(tool.directory)) stop('Please input the ZDOCK directory path')
  if(missing(ligand)) stop('Please input either a ligand PDB ID or a path to the PDB file')
  if(missing(ligand.name)) stop('Please input the ligand name!')
  if(missing(structure)) structure <- 'structure'
  if(missing(additional.docking.params)) additional.docking.params <- list()
  if(missing(parallel)) parallel <- T

  call_zdock <- function(steropodon.object,
                         ligand,
                         ligand.name,
                         structure,
                         zdock.directory,
                         zdock.n.predictions,
                         zdock.specific.regions.to.mask,
                         zdock.region.columns,
                         zdock.fixed.receptor,
                         parallel
                         ){

    if(missing(steropodon.object)) stop('Please input your Steropodon object!')
    if(missing(ligand)) stop('Please input either a ligand PDB ID or a path to the PDB file')
    if(missing(ligand.name)) stop('Please input the ligand name!')
    if(missing(structure)) structure <- 'structure'
    if(missing(zdock.directory)) stop('Please input the ZDOCK directory path')
    if(missing(zdock.n.predictions)) zdock.n.predictions <- 2000
    if(missing(zdock.specific.regions.to.mask)) zdock.specific.regions.to.mask <- c('VDJ_FR1', 'VDJ_FR2', 'VDJ_FR3', 'VDJ_FR4', 'VJ_FR1', 'VJ_FR2', 'VJ_FR3', 'VJ_FR4')
    if(missing(zdock.region.columns)) zdock.region.columns <- c('chain', 'region')
    if(missing(zdock.fixed.receptor)) zdock.fixed.receptor <- F
    if(missing(parallel)) parallel

    #zdock.directory <- '/Users/tudorcotet/Desktop/zdock'

    call_zdock_lapply <- function(receptor.file,
                                  ligand.file,
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
          zdock.directory, '/zdock', ' -R ', receptor.file, ' -L ', ligand.file, ' -o ', out_file_name, ' -N ', zdock.n.predictions, fixed_receptor
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
    chain_dict <- setNames(names(chain_dict), chain_dict)

    receptor_files <- list.files(temp_dir)

    ligand_file <- paste0(ligand.name, '.pdb')
    ligand_file_path <- paste0(temp_dir, '/', ligand_file)
    ligand_pdb <- bio3d::write.pdb(bio3d::read.pdb(ligand), ligand_file_path)

    setwd(temp_dir)

    system(
      paste0(
        'cp ', zdock.directory, '/uniCHARMM ', '.'
      )
    )

    system(
      paste0(
        zdock.directory, '/mark_sur ', ligand_file, ' ', ligand_file
      )
    )

    partial_function <- function(x) call_zdock_lapply(receptor.file = x,
                                                      ligand.file = ligand_file,
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
         parallel::parlapply(cl, receptor_files, partial_function)
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
      docked_ligand_pdb <- pdb$atom[(nrow(receptor_pdb$atom)+1):nrow(pdb$atom),]
      final_resno <- max(docked_receptor_pdb$resno)

      ligand_resnos <- unique(docked_ligand_pdb$resno)
      ligand_resnos <- ligand_resnos[order(ligand_resnos)]
      resno_dict <- final_resno + 1:length(ligand_resnos)
      names(resno_dict) <- ligand_resnos
      docked_ligand_pdb$resno <- resno_dict[as.character(docked_ligand_pdb$resno)]

      docked_receptor_pdb$chain <- chain_dict[docked_receptor_pdb$chain]
      docked_ligand_pdb$chain <- 'antigen'

      docked_receptor_pdb$region <- receptor_pdb$atom$region
      docked_ligand_pdb$region <- 'antigen'

      missing_features <- setdiff(colnames(receptor_pdb$atom), colnames(docked_receptor_pdb))

      for(feature in missing_features){
        docked_receptor_pdb[[feature]] <- receptor_pdb$atom[[feature]]
        docked_ligand_pdb[[feature]] <- rep(NA, nrow(docked_ligand_pdb))
      }


      pdb$atom <- rbind(docked_receptor_pdb, docked_ligand_pdb)

      ligand_structure <- split_structure(pdb, grouping = c('chain'), specific.values = 'antigen')
      receptor_structure <- split_structure(pdb, grouping = c('chain'), specific.values = c('VDJ','VJ'), combine.values = T, combine.groupings = F)

      #Maybe also modify CDRH3s and other split structures that already exist in the Steropodon object
      steropodon_list[[id_1]]@complex <- pdb
      steropodon_list[[id_1]]@antigen <- ligand_structure
      steropodon_list[[id_1]] <- modify_structure(steropodon_list[[id_1]], structure = structure, pdb = receptor_structure)
    }

    names(steropodon_list) <- steropodon_names
    steropodon_object <- nest_steropodon(steropodon_list)
    return(steropodon_object)
  }


  if(docking.tool == 'zdock'){
    params <- list(steropodon.object = steropodon.object,
                   zdock.directory = tool.directory,
                   ligand = ligand,
                   ligand.name = ligand.name,
                   structure = structure,
                   parallel = parallel
                  )

    steropodon_object <- do.call(call_zdock, c(params, additional.docking.params))


  }else{
    stop('Docking tool not implemented yet!')
  }

  return(steropodon_object)
}
