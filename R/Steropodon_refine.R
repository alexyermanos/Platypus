#' @importFrom reticulate %as%


Steropodon_refine <- function(steropodon.object,
                              structure,
                              use.conda,
                              env.name,
                              parallel,
                              refining.tools,
                              additional.tool.params){


   if(missing(steropodon.object)) stop('Please input your Steropodon object!')
   if(missing(structure)) structure <- 'structure'
   if(missing(use.conda)) use.conda <- T
   if(missing(env.name)) env.name <- NULL
   if(missing(parallel)) parallel <- T
   if(missing(refining.tools)) refining.tools <- 'openmm'
   if(missing(additional.tool.params)) additional.tool.params <- list()


   merge_refined_pdb <- function(original.pdb, refined.pdb, chain.dict){
     chain.dict <- setNames(names(chain.dict), chain.dict)
     refined.pdb$atom$chain <- chain.dict[refined.pdb$atom$chain]

     cols <- setdiff(colnames(original.pdb$atom), colnames(refined.pdb$atom))

     for(col in cols){
       atom_subset <- original.pdb$atom[,c('resno', 'chain', col)]
       refined.pdb$atom <- dplyr::left_join(refined.pdb$atom, atom_subset, by = c("resno","chain"))
     }

     alpha_chain_subset <- refined.pdb$atom[refined.pdb$atom$elety == 'CA',]
     refined.pdb$atom$resno <- 0:(nrow(alpha_chain_subset) - 1)
     return(refined.pdb)
   }

   call_openmm <- function(steropodon.object,
                           structure,
                           openmm.stiffness,
                           openmm.tolerance,
                           openmm.use.gpu,
                           use.conda,
                           env.name,
                           parallel){

     if(missing(steropodon.object)) stop('Please input your Steropodon object!')
     if(missing(structure)) structure <- 'structure'
     if(missing(openmm.stiffness)) openmm.stiffness <- 10.
     if(missing(openmm.tolerance)) openmm.tolerance <- 2.39
     if(missing(openmm.use.gpu)) openmm.use.gpu <- F
     if(missing(use.conda)) use.conda <- T
     if(missing(env.name)) env.name <- 'openmm'
     if(missing(parallel)) parallel <- T

     openmm_refine <- function(pdb.file,
                               stiffness = 10.,
                               tolerance = 2.39,
                               use.gpu = F){

       openmm <- reticulate::import('openmm')
       pdbfixer <- reticulate::import('pdbfixer')

       ENERGY = openmm$unit$kilocalories_per_mole
       LENGTH = openmm$unit$angstroms

       tolerance = openmm$unit$quantity$Quantity(value = tolerance, unit = ENERGY)
       stiffness = openmm$unit$quantity$Quantity(value = stiffness, unit = ENERGY)
       stiffness = stiffness$`__div__`(LENGTH$`__pow__`(2))

       fixer = pdbfixer$PDBFixer(pdb.file)

       fixer$findMissingResidues()
       fixer$findMissingAtoms()
       fixer$addMissingAtoms()

       force_field = openmm$app$ForceField("amber14/protein.ff14SB.xml")


       modeller = openmm$app$Modeller(fixer$topology, fixer$positions)
       modeller$addHydrogens(force_field)
       system = force_field$createSystem(modeller$topology)

       force = openmm$CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
       force$addGlobalParameter("k", stiffness)

       for(particle in c("x0", "y0", "z0")){
         force$addPerParticleParameter(particle)
       }


       for(residue in reticulate::iterate(modeller$topology$residues())) {
         for(atom in reticulate::iterate(residue$atoms())) {
           if(atom$name %in% c('N', 'CA', 'C', 'CB')){
             force$addParticle(atom$index, modeller$positions[atom$index])
           }
         }
       }

       system$addForce(force)

       integrator = openmm$LangevinIntegrator(0, 0.01, 1.0)
       if(use.gpu){
         platform = openmm$Platform$getPlatformByName("CUDA")
       }else{
         platform = openmm$Platform$getPlatformByName("CPU")
       }

       simulation = openmm$app$Simulation(modeller$topology, system, integrator, platform)
       simulation$context$setPositions(modeller$positions)
       simulation$minimizeEnergy(tolerance)

       py <- reticulate::import_builtins()

       with(py$open(pdb.file, "w") %as% f,{
         openmm$app$PDBFile$writeFile(
             simulation$topology,
             simulation$context$getState(getPositions = T)$getPositions(),
             f,
             keepIds=T)
       })

       pdb <- bio3d::read.pdb(pdb.file)
       return(pdb)
     }


     if(is.null(env.name)){
       env.name <- 'openmm'
     }

     if(use.conda){
       envs <- reticulate::conda_list()
       if(!(env.name %in% envs$name)){
         reticulate::conda_create(env.name)
       }
       reticulate::use_condaenv(env.name)
       if(!reticulate::py_module_available('openmm')){
         reticulate::conda_install(envname = env.name, packages = c('openmm'))
       }

       if(!reticulate::py_module_available('pdbfixer')){
         reticulate::conda_install(envname = env.name, packages = c('pdbfixer'))
       }
     }else{
       envs <- reticulate::virtualenv_list()
       if(!(env.name %in% envs)){
         reticulate::virtualenv_create(env.name)
       }
       reticulate::use_virtualenv(env.name)

       if(!reticulate::py_module_available('openmm')){
         reticulate::py_install('openmm')
       }
       if(!reticulate::py_module_available('pdbfixer')){
         reticulate::py_install('pdbfixer')
       }
     }

     if(!reticulate::py_module_available('openmm')){
       stop('Could not find/install OpenMM in your environment. Please provide a conda/virtualenv environment with OpenMM already installed!')
     }

     if(!reticulate::py_module_available('pdbfixer')){
       stop('Could not find/install pdbfixer in your environment. Please provide a conda/virtualenv environment with pdbfixer already installed!')
     }

     switch(Sys.info()[['sysname']],
             Windows= {message("Windows system detected")
               operating.system <- "Windows"},
             Linux  = {message("Linux system detected")
               operating.system <- "Linux"},
             Darwin = {message("MAC system detected")
               operating.system <- "Darwin"})

     #Create temporary directory
     temp_dir <- './openmm_temp'
     if(!dir.exists(temp_dir)) dir.create(temp_dir)

     #Preprocessing/ changing resno values for numbering schemes
     steropodon_list <- unnest_steropodon(steropodon.object)
     pdb_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
     pdb_list <- lapply(pdb_list, function(pdb) {non_na_resno <- unique(pdb$atom$resno[!is.na(pdb$atom$chain)])
                                                 pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = non_na_resno))
                                                 vj <- pdb$atom[pdb$atom$chain == 'VJ',]
                                                 vj$resno <- vj$resno - min(vj$resno)
                                                 pdb$atom[pdb$atom$chain == 'VJ',] <- vj
                                                 pdb$atom$resno <- 1 + pdb$atom$resno
                                                 return(pdb)}
                                                 )

     steropodon_list <- mapply(function(x, y) modify_structure(x, structure = structure, pdb = y), steropodon_list, pdb_list)
     steropodon.object <- nest_steropodon(steropodon_list)
     out <- write_steropodon_pdbs(steropodon.object, structure = structure, dir = temp_dir)

     pdb_files <- out$file_list
     chain_dict <- out$chain_dict

     partial_function <- function(x) {try(openmm_refine(pdb.file = x,
                                                    stiffness = openmm.stiffness,
                                                    tolerance = openmm.tolerance,
                                                    use.gpu = openmm.use.gpu)
                                               )
                                       }

     if(parallel){
       cores <- parallel::detectCores() - 1

       if(operating.system %in% c('Linux', 'Darwin')){
         pdb_list_refined <- parallel::mclapply(pdb_files, partial_function, mc.cores = cores)

       }else{
         cl <- parallel::makeCluster(cores)
         pdb_list_refined <- parallel::parlapply(cl, pdb_files, partial_function)
         parallel::stopCluster(cl)
       }

     }else{
       pdb_list_refined <- lapply(pdb_files, function(x) partial_function(x))
     }


     is.error <- function(x) inherits(x, "try-error")
     failed <- lapply(pbb_list_refined, function(x) is.error(x))
     if(length(failed)!=0){
       message(paste0('OpenMM failed on the following structures: ', paste0(names(pdb_list)[failed], collapse = ', ')))
       pdb_list_refined[[failed]] <- pdb_list[[failed]]
       pdb_list_refined[[!failed]] <- mapply(function(x,y) merge_refined_pdb(original.pdb = x, refined.pdb = y, chain.dict = chain_dict), pdb_list[[!failed]], pdb_list_refined[[!failed]])
     }else{
       pdb_list_refined <- mapply(function(x,y) merge_refined_pdb(original.pdb = x, refined.pdb = y, chain.dict = chain_dict), pdb_list, pdb_list_refined)
     }

     steropodon_list <- mapply(function(x, y) modify_structure(x, structure = structure, pdb = y), steropodon_list, pdb_list_refined)
     steropodon.object <- nest_steropodon(steropodon_list)

     unlink(temp_dir, recursive = T)
     return(steropodon.object)
   }

   call_ablooper <- function(steropodon.object,
                             structure,
                             ablooper.chains,
                             ablooper.numbering.scheme,
                             ablooper.refine,
                             use.conda,
                             env.name,
                             parallel
                            ){

     if(missing(steropodon.object)) stop('Pleas input your Steropodon object!')
     if(missing(structure)) structure <- 'structure'
     if(missing(ablooper.chains)) ablooper.chains <- list('A', 'B')
     if(missing(ablooper.numbering.scheme)) ablooper.numbering.scheme <- 'imgt'
     if(missing(ablooper.refine)) ablooper.refine <- F
     if(missing(use.conda)) use.conda <- T
     if(missing(env.name)) env.name <- 'ablooper'
     if(missing(parallel)) parallel <- T

     ablooper_refine <- function(pdb.file,
                                 ablooper.chains = list('A', 'B'),
                                 ablooper.numbering.scheme = 'imgt',
                                 ablooper.refine = T
                                ){
        ablooper <- reticulate::import('ABlooper')
        pred <- ablooper$CDR_Predictor(pdb.file,
                                       chains = reticulate::tuple(list('A', 'B')),
                                       model = reticulate::r_to_py(ablooper.numbering.scheme),
                                       refine = reticulate::r_to_py(ablooper.refine)
                                       )

        pred$write_predictions_in_pdb_format(pdb.file)

        pdb <- bio3d::read.pdb(pdb.file)
        pdb$ablooper.cdr.rmsd <- pred$calculate_BB_rmsd_wrt_input()

        return(pdb)
     }

     if(is.null(env.name)){
      env.name <- 'ablooper'
     }

     if(use.conda){
      envs <- reticulate::conda_list()
      if(!(env.name %in% envs$name)){
        reticulate::conda_create(env.name)
      }
      reticulate::use_condaenv(env.name)
      if(!reticulate::py_module_available('ABlooper')){
        system('pip install ABlooper')
      }

     }else{
      envs <- reticulate::virtualenv_list()
      if(!(env.name %in% envs)){
        reticulate::virtualenv_create(env.name)
      }
      reticulate::use_virtualenv(env.name)
      if(!reticulate::py_module_available('ABlooper')){
        system('pip install ABlooper')
      }
     }

     if(!reticulate::py_module_available('ABlooper')){
      stop('Could not find/install ABlooper in your environment. Please provide a conda/virtualenv environment with ABlooper already installed!')
     }

     switch(Sys.info()[['sysname']],
             Windows= {message("Windows system detected")
               operating.system <- "Windows"},
             Linux  = {message("Linux system detected")
               operating.system <- "Linux"},
             Darwin = {message("MAC system detected")
               operating.system <- "Darwin"})

     #Create temporary directory
     temp_dir <- './ablooper_temp'
     if(!dir.exists(temp_dir)) dir.create(temp_dir)

     #Preprocessing/ changing resno values for numbering schemes
     steropodon_list <- unnest_steropodon(steropodon.object)
     pdb_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
     if(!('VDJ' %in% pdb_list[[1]]$atom$chain) | !('VJ' %in% pdb_list[[1]]$atom$chain)){
       stop('ABlooper can only be used on antibody structures with both the heavy and light chain!')
     }
     pdb_list <- lapply(pdb_list, function(pdb) {non_na_resno <- unique(pdb$atom$resno[!is.na(pdb$atom$chain)])
                                                 pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = non_na_resno))
                                                 vj <- pdb$atom[pdb$atom$chain == 'VJ',]
                                                 vj$resno <- vj$resno - min(vj$resno)
                                                 pdb$atom[pdb$atom$chain == 'VJ',] <- vj
                                                 pdb$atom$resno <- 1 + pdb$atom$resno
                                                 i <<- i + 1
                                                 return(pdb)}
                                                 )

     steropodon_list <- mapply(function(x, y) modify_structure(x, structure = structure, pdb = y), steropodon_list, pdb_list)
     steropodon.object <- nest_steropodon(steropodon_list)
     out <- write_steropodon_pdbs(steropodon.object, structure = structure, dir = temp_dir)

     pdb_files <- out$file_list
     chain_dict <- out$chain_dict

     partial_function <- function(x) {ablooper_refine(pdb.file = x,
                                                      ablooper.refine = ablooper.refine,
                                                      ablooper.numbering.scheme = ablooper.numbering.scheme,
                                                      ablooper.chains = ablooper.chains
                                                      )}

     if(parallel){
       cores <- parallel::detectCores() - 1

       if(operating.system %in% c('Linux', 'Darwin')){
         pdb_list_refined <- parallel::mclapply(pdb_files, partial_function, mc.cores = cores)

       }else{
         cl <- parallel::makeCluster(cores)
         pdb_list_refined <- parallel::parlapply(cl, pdb_files, partial_function)
         parallel::stopCluster(cl)
       }

     }else{
       pdb_list_refined <- lapply(pdb_files, function(x) partial_function(x))
     }


     is.error <- function(x) inherits(x, "try-error")
     failed <- lapply(pbb_list_refined, function(x) is.error(x))
     if(length(failed)!=0){
       message(paste0('OpenMM failed on the following structures: ', paste0(names(pdb_list)[failed], collapse = ', ')))
       pdb_list_refined[[failed]] <- pdb_list[[failed]]
       pdb_list_refined[[!failed]] <- mapply(function(x,y) merge_refined_pdb(original.pdb = x, refined.pdb = y, chain.dict = chain_dict), pdb_list[[!failed]], pdb_list_refined[[!failed]])
     }else{
       pdb_list_refined <- mapply(function(x,y) merge_refined_pdb(original.pdb = x, refined.pdb = y, chain.dict = chain_dict), pdb_list, pdb_list_refined)
     }

     steropodon_list <- mapply(function(x, y) modify_structure(x, structure = structure, pdb = y), steropodon_list, pdb_list_refined)
     steropodon.object <- nest_steropodon(steropodon_list)

     unlink(temp_dir, recursive = T)
     return(steropodon.object)
   }


   call_refining_tools <- function(steropodon.object,
                                   structure,
                                   use.conda,
                                   env.name,
                                   parallel,
                                   refining.tools,
                                   additional.tool.params
                                   ){

     params <- list(steropodon.object = steropodon.object,
                    structure = structure,
                    use.conda = use.conda,
                    env.name = env.name,
                    parallel = parallel)

     if('ablooper' %in% refining.tools){
       steropodon.object <- do.call(call_ablooper, c(params, additional.tool.params))
     }

     if('pyrosetta' %in% refining.tools){
       stop('PyRosetta is not implemented yet!')
     }

     if('openmm' %in% refining.tools){
       steropodon.object <- do.call(call_openmm, c(params, additional.tool.params))
     }

     if('openbabel' %in% refining.tools){
       stop('OpenBabel is not implemented yet!')

     }

     if('amber' %in% refining.tools){
       stop('Amber is not implemented yet!')
     }


     return(steropodon.object)
   }


   steropodon.object <- call_refining_tools(steropodon.object = steropodon.object,
                                            structure = structure,
                                            use.conda = use.conda,
                                            env.name = env.name,
                                            parallel = parallel,
                                            refining.tools = refining.tools,
                                            additional.tool.params = additional.tool.params)

   return(steropodon.object)
}
