Steropodon_properties <- function(steropodon.object,
                                  structure,
                                  properties,
                                  propka.directory,
                                  dssp.exefile,
                                  parallel
                                 ){

   if(missing(steropodon.object)) stop('Please input your Steropodon object!')
   if(missing(structure)) structure <- 'structure'
   if(missing(properties)) properties <- c('SASA', 'charge', 'hydrophobicity', 'pKa', 'DSSP')
   if(missing(propka.directory) & ('pKa' %in% properties)) stop('Please input your directory to Propka for pKa calculations')
   if(missing(dssp.exefile) & ('DSSP' %in% properties)) stop('Please input path to the DSPP exe file')
   if(missing(parallel)) parallel <- T

   switch(Sys.info()[['sysname']],
           Windows= {message("Windows system detected")
             operating.system <- "Windows"},
           Linux  = {message("Linux system detected")
             operating.system <- "Linux"},
           Darwin = {message("MAC system detected")
             operating.system <- "Darwin"})


    #TO DO:
    calculate_interatomic <- function(){
    }

    #TO DO:
    calculate_sap <- function(pdb){
    }

    #TO DO:
    calculate_pi <- function(pdb){
    }


    #From Lucas' VDJ_structure_analysis function
    calculate_sasa <- function(pdb){

      SASA <- suppressWarnings(vanddraabe::FreeSASA.diff(pdb$atom))
      pdb$atom$SASA.lost <- SASA$SASA.lost
      pdb$atom$SASA.prot <- SASA$SASA.prot
      pdb$atom$SASA.hetatm <- SASA$SASA.hetatm

      return(pdb)
    }

    #From Lucas' VDJ_structure_analysis function
    calculate_charge <- function(pdb){
      df_charge <- data.frame()

      pdb$atom$charge <- NULL
      for(CHAIN in unique(pdb$atom$chain)){
        df_chain <- pdb$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
        df_chain$charge <- pdb$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::charge()
        df_charge <- rbind(df_charge, df_chain)
      }

      pdb$atom <- dplyr::left_join(pdb$atom, df_charge, by = c("resno","chain"))
      return(pdb)
    }

    #From Lucas' VDJ_structure_analysis function
    calculate_hydrophobicity <- function(pdb){
      df_hydph <- data.frame()

      for(CHAIN in unique(pdb$atom$chain)){
        df_chain <- pdb$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
        df_chain$hydrophobicity <- pdb$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::hydrophobicity()

        df_hydph <- rbind(df_hydph,df_chain)
      }

      pdb$atom <- dplyr::left_join(pdb$atom, df_hydph, by = c("resno","chain"))
      return(pdb)
    }


    calculate_dssp <- function(pdb,
                               dssp.exefile = '/opt/homebrew/Caskroom/miniforge/base/bin/mkdssp'){
      template_pdb_columns <- c('type','eleno','elety','alt','resid','chain','resno','insert','x','y','z','o','b','segid','elesy','charge')

      chains <- unique(temp_pdb$atom$chain)
      chain_dict <- toupper(letters)[1:length(chains)]
      names(chain_dict) <- chains
      pdb$atom$chain <- chain_dict[pdb$atom$chain]

      temp_pdb <- pdb
      temp_pdb$atom <- temp_pdb$atom[,template_pdb_columns]

      out <- bio3d::dssp(temp_pdb, exefile = dssp.exefile)

      dssp_df <- data.frame(out$sse)
      dssp_df$id <- rownames(dssp_df)
      colnames(dssp_df) <- c('sse', 'id')
      t <- dssp_df$id %>% stringr::str_split("_") %>% do.call(rbind, .) %>% data.frame()
      colnames(t) <- c('resno', 'chain', 'a')
      dssp_df <- cbind(t, dssp_df)
      dssp_df$a <- NULL
      dssp_df$id <- NULL
      dssp_df$phi <- out$phi
      dssp_df$psi <- out$psi
      dssp_df$acc <- out$acc
      dssp_df$sse[dssp_df$sse == ' '] <- NA
      dssp_df$resno <- as.integer(dssp_df$resno)

      pdb$atom <- dplyr::left_join(pdb$atom, dssp_df, by = c("resno","chain"))
      chain_dict <- setNames(names(chain_dict), chain_dict)
      pdb$atom$chain <- chain_dict[pdb$atom$chain]

      return(pdb)
    }

    calculate_pka <- function(steropodon.object,
                              structure,
                              propka.directory = '/Users/tudorcotet/propka-3.0'
                              ){

      pdb <- select_structure(steropodon.object, structure = structure)
      out <- write_steropodon_pdbs(steropodon.object, structure = structure, dir = temp.dir)
      file_name <- paste0(steropodon.object@structure_id, collapse = '_')
      propka_name <- paste0(file_name, '.pka')

      file_name <- out$file_list
      chain_dict <- out$chain_dict

      system(
        paste0(
          'python ', propka.directory, '/propka.py ', file_name
        )
      )

      propka_out <- readLines(propka_name)
      start <- grep('SUMMARY OF THIS PREDICTION', propka_out)
      end <- grep('Free energy of', propka_out)
      pka <- propka_out[(start+1):(end-3)]
      pka <- paste(pka, collapse = '\n')

      tc <- textConnection(pka)
      df <- read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
      close(tc)

      df <- df[-1,]
      colnames(df) <- c('resid', 'resno', 'chain', 'pKa', 'pKmodel')
      #df$resid <- NULL
      df$resno <- as.integer(df$resno)
      df$pKa <- as.numeric(df$pKa)
      df$pKmodel <- as.numeric(df$pKmodel)

      chain_dict <- setNames(names(chain_dict), chain_dict)
      df$chain <- chain_dict[df$chain]
      pdb$atom <- dplyr::left_join(pdb$atom, df, by = c('resno', 'chain', 'resid'))

      start <- grep('Free energy of', propka_out)
      end <- grep('The pH of optimum stability', propka_out)

      en <- propka_out[(start+1):(end-2)]
      en <- paste(en, collapse = '\n')

      tc <- textConnection(en)
      df <- read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
      close(tc)
      colnames(df) <- c('pH', 'Free energy of folding (kcal/mol')

      pdb$propka.energy <- df

      start <- grep('Protein charge of folded and unfolded state as a function of pH', propka_out)
      end <- grep('The pI is', propka_out)

      ch <- propka_out[(start+2):(end-1)]

      tc <- textConnection(ch)
      df <- read.table(tc, as.is=T, fill=T, blank.lines.skip=F)
      close(tc)

      colnames(df) <- c('pH', 'unfolded', 'folded')
      pdb$propka.charge <- df

      start <- grep('The pI is', propka_out)
      pI <- propka_out[start]
      pI <- stringr::str_split(pI, ' ')[[1]]
      pdb$propka.pI <- list(folded = pI[5], unfolded = pI[9])

      start <- grep('The pH of optimum stability', propka_out)
      optimum <- propka_out[start]
      optimum <- stringr::str_split(optimum, ' ')[[1]]

      pdb$propka.optimum.stability <- list(pH = optimum[8], free_energy = optimum[16])

      unlink(propka_name)
      unlink(file_name)
      return(pdb)
    }


    calculate_properties <- function(steropodon.object,
                                     structure,
                                     properties,
                                     propka.directory,
                                     dssp.exefile
                                    ){
      pdb <- select_structure(steropodon.object, structure = structure)

      if('SASA' %in% properties){
        pdb <- calculate_sasa(pdb)
      }

      if('charge' %in% properties){
        pdb <- calculate_charge(pdb)
      }

      if('hydrophobicity' %in% properties){
        pdb <- calculate_hydrophobicity(pdb)
      }

      if('pKa' %in% properties){
        pdb <- calculate_pka(steropodon.object,
                             structure = structure,
                             propka.directory = propka.directory)
      }

      if('DSSP' %in% properties){
        pdb <- calculate_dssp(pdb,
                              dssp.exefile = dssp.exefile)
      }

      if('SAP' %in% properties){
        stop('Not implemented yet!')
      }

      if('dI' %in% properties){
        stop('Not implemented yet!')
      }

      if('interatomic' %in% properties){
        stop('Not implemented yet!')
      }

      steropodon_object <- modify_structure(steropodon.object, structure = structure, pdb = pdb)
      return(steropodon_object)
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

   steropodon_list <- unnest_steropodon(steropodon.object)
   ids <- names(steropodon_list)
   partial_function <- function(x) calculate_properties(x,
                                                        properties = properties,
                                                        structure = structure,
                                                        propka.directory = propka.directory,
                                                        dssp.exefile = dssp.exefile)
   if(parallel){
      cores <- parallel::detectCores() - 1

      if(operating.system %in% c('Linux', 'Darwin')){
        steropodon_list <- parallel::mclapply(steropodon_list, partial_function, mc.cores = cores)

      }else{
        cl <- parallel::makeCluster(cores)
        steropodon_list <- parallel::parLapply(cl, steropodon_list, partial_function)
        parallel::stopCluster(cl)
      }

   }else{
      steropodon_list <- lapply(receptor_files, function(x) partial_function(x))
   }

   names(steropodon_list) <- ids
   steropodon_object <- nest_steropodon(steropodon_list)

   return(steropodon_object)
}
