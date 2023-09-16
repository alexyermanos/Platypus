#' Folds receptor sequences from a given dataframe using ColabFold
#' @description Folds receptor sequences from a given dataframe using ColabFold
#' @param sequence.df dataframe - sequences and their IDs, as obtained from the Steropodon_model function.
#' @param colabfold.num.recycle See ColabFold for more information.
#' @param colabfold.num.ensemble See ColabFold for more information.
#' @param colabfold.num.models See ColabFold for more information.
#' @param colabfold.msa.mode See ColabFold for more information.
#' @param colabfold.amber See ColabFold for more information.
#' @param colabfold.stop.at.score See ColabFold for more information.
#' @param colabfold.stop.at.score.below See ColabFold for more information.
#' @param env.name  string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.
#' @param gpu boolean - if TRUE, will use GPUs for faster structural inference.

#' @return input dataframe with a new column with the folded structure path. Will be processed by Steropodon_model into a Steropodon object/list of objects.
#' @export
#' @examples
#' \dontrun{
#' folded_df <- call_colabfold(sequence_df,
#' env.name = 'colabfold',
#' use.conda = T,
#' gpu = T)
#'}
call_colabfold <- function(sequence.df,
                           colabfold.num.recycle,
                           colabfold.num.ensemble,
                           colabfold.num.models,
                           colabfold.msa.mode,
                           colabfold.amber,
                           colabfold.stop.at.score,
                           colabfold.stop.at.score.below,
                           colabfold.local.database.directory,
                           env.name,
                           use.conda,
                           gpu
                           ){

   if(missing(sequence.df)) stop('Please input the processed sequence dataframe!')
   if(missing(colabfold.num.recycle)) colabfold.num.recycle <- 5
   if(missing(colabfold.num.ensemble)) colabfold.num.ensemble <- 1
   if(missing(colabfold.num.models)) colabfold.num.models <- 5
   #if(missing(colabfold.msa.mode)) colabfold.msa.mode <- 'MMseqs\ (UniRef+Environmental)'
   if(missing(colabfold.amber)) colabfold.amber <- T
   if(missing(colabfold.stop.at.score)) colabfold.stop.at.score <- 100
   if(missing(colabfold.stop.at.score.below)) colabfold.stop.at.score.below <- 0
   if(missing(colabfold.local.database.directory)) colabfold.local.database.directory <- NULL
   if(missing(env.name)) env.name <- NULL
   if(missing(use.conda)) use.conda <- T
   if(missing(gpu)) gpu <- T

   create_fasta_directory <- NULL

   if(Sys.which('colabfold_batch')==''){
     message('Could not find a local installation of colabfold_batch. Will try to install it in a conda/virtualenv environment. However, it is recommended to install ColabFold manually by following the instructions here: https://github.com/sokrypton/ColabFold, https://github.com/YoshitakaMo/localcolabfold and exporting the colabfold_batch and colabfold_search $PATH variables')

     if(is.null(env.name)){
      env.name <- 'colabfold-conda'
     }

     if(use.conda){
       envs <- reticulate::conda_list()
       if(!(env.name %in% envs$name)){
         reticulate::conda_create(env.name)
       }
       reticulate::use_condaenv(env.name)

       if(!reticulate::py_module_available('kalign')){
         system('conda install -c conda-forge -c bioconda kalign2=2.04')
       }

       if(!reticulate::py_module_available('hhsuite')){
         system('conda install -c conda-forge -c bioconda hhsuite=3.3.0')
       }

     }else{

       envs <- reticulate::virtualenv_list()
       if(!(env.name %in% envs)){
         #reticulate::install_python(version = '3.7.9')
         reticulate::virtualenv_create(env.name)
       }

       reticulate::use_virtualenv(env.name)
     }

     if(!reticulate::py_module_available('colabfold')){
       system('pip install "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold"')
     }

     if(!reticulate::py_module_available('jax')){
       system('pip install -q "jax[cuda]>=0.3.8,<0.4" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html')
     }
   }

   temp_dir <- './tempdir_steropodon_fasta'
   pdb_dir <- './tempdir_steropodon_colabfold_pdbs'

   create_fasta_directory(sequence.df, model = 'colabfold')

   if(length(list.files(temp_dir)) == 0){
     stop('Unable to find fasta files to run the ColabFold model for!')
   }

   #Second check to ensure ColabFold was installed
   if(Sys.which('colabfold_batch')==''){
     'Please install the most recent (local) version of ColabFold by following the instructions here: https://github.com/sokrypton/ColabFold, https://github.com/YoshitakaMo/localcolabfold'
   }

   if(!dir.exists(pdb_dir)) dir.create(pdb_dir)


   #Call on MMSEQ2 server (careful regarding number of sequences/queries)
   if(is.null(colabfold.local.database.directory)){
     model_call <- paste0('colabfold_batch ', temp_dir, ' ', pdb_dir, ' ',
              '--stop-at-score ', colabfold.stop.at.score, ' ',
              '--stop-at-score-below ', colabfold.stop.at.score.below, ' ',
              '--num-ensemble ', colabfold.num.ensemble, ' ',
              '--num-recycle ', colabfold.num.recycle, ' ',
              '--num-models ', colabfold.num.models
            #  '--msa-mode ', colabfold.msa.mode
             )

     if(colabfold.amber){
       model_call <- paste0(model_call, ' --amber')
     }

     if(!gpu){
       model_call <- paste0(model_call, ' --cpu')
     }

     system(model_call)


   #Call on local database
   }else{
     model_call <- paste0('colabfold_batch ', msa_dir, ' ', pdb_dir, ' ',
             '--stop-at-score ', colabfold.stop.at.score, ' ',
             '--stop-at-score-below ', colabfold.stop.at.score.below, ' ',
             '--num-ensemble ', colabfold.num.ensemble, ' ',
             '--num-recycle ', colabfold.num.recycle, ' ',
             '--num-models ', colabfold.num.models, ' ',
             '--msa-mode ', colabfold.msa.mode
             )

     if(colabfold.amber){
       model_call <- paste0(model_call, ' --amber')
     }

     if(!gpu){
       model_call <- paste0(model_call, ' --cpu')
     }

     if(Sys.which('colabfold_search')==''){
       'Unable to find colabfold_search to perform MSA on a local database. Please install ColabFold following the instructions here: https://github.com/sokrypton/ColabFold, https://github.com/YoshitakaMo/localcolabfold'
     }

     msa_dir <- './tempdir_steropodon_msa'
     if(!dir.exists(msa_dir)) dir.create(msa_dir)

     system(
       paste0(
         'colabfold_search ', temp_dir, ' ', colabfold.local.database.directory, ' ', msa_dir
       )
     )

     system(model_call)

     unlink(msa_dir, recursive = T)
   }
   #Process resulting PDB files, save resulting plots, cleanup

   if(colabfold.amber){
     sequence.df$pdb_file <- paste0(sequence.df$fasta_name, '_relaxed')

   }else{
     sequence.df$pdb_file <- paste0(sequence.df$fasta_name, '_unrelaxed')
   }

   all_files <- list.files(pdb_dir)
   sequence.df$pdb_file <- lapply(sequence.df$pdb, function(x) {x <- all_files[stringr::str_detect(all_files, x) & stringr::str_detect(all_files, '.pdb')]
                                                                if(length(x) != 0){
                                                                  x <- paste0(pdb_dir, '/', x)
                                                                }else{
                                                                  x <- NA
                                                                }
                                                                return(x)
                                                                })
   return(list(pdb_dir = pdb_dir, sequence_df = sequence.df))
}


#' Folds receptor sequences from a given dataframe using IgFold
#' @description Folds receptor sequences from a given dataframe using IgFold
#' @param sequence.df dataframe - sequences and their IDs, as obtained from the Steropodon_model function.
#' @param igfold.refine boolean - if TRUE, will refine the modelled structure using OpenMM or PyRosetta.
#' @param igfold.renum boolean - if TRUE, will renumber the modelled structure chains.
#' @param igfold.use.openmm boolean - if TRUE, will use OpenMM for structure relaxation. Else, will use PyRosetta.
#' @param igfold.use.abnum boolean - if TRUE, will use AbNum for structure/region renumbering.
#' @param env.name  string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.

#' @return input dataframe with a new column with the folded structure path. Will be processed by Steropodon_model into a Steropodon object/list of objects.
#' @export
#' @examples
#' \dontrun{
#' folded_df <- call_igfold(sequence_df,
#' env.name = 'igfold',
#' use.conda = T)
#'}
call_igfold <- function(sequence.df,
                        igfold.refine,
                        igfold.renum,
                        igfold.use.openmm,
                        igfold.use.abnum,
                        env.name,
                        use.conda
                        ){

  if(missing(sequence.df)) stop('Please input the processed sequence dataframe!')
  if(missing(igfold.refine)) igfold.refine <- T
  if(missing(igfold.renum)) igfold.renum <- F
  if(missing(igfold.use.openmm)) igfold.use.openmm <- T
  if(missing(igfold.use.abnum)) igfold.use.abnum <- F
  if(missing(env.name)) env.name <- NULL
  if(missing(use.conda)) use.conda <- T

  fasta_sequence <- NULL


  pdb_dir <- './tempdir_steropodon_pdbs_igfold_tnfr2'
  if(!dir.exists(pdb_dir)) dir.create(pdb_dir)
  #Setup the reticulate conda environment
  #options(reticulate.repl.quiet = TRUE)
  if(is.null(env.name)){
    env.name <- 'igfold'
  }

  if(use.conda){

    envs <- reticulate::conda_list()
    if(!(env.name %in% envs$name)){
      reticulate::conda_create(env.name)
    }
    reticulate::use_condaenv(env.name)

    if(!reticulate::py_module_available('anarci')){
      system('conda install -c bioconda anarci')
    }

  }else{

    envs <- reticulate::virtualenv_list()
    if(!(env.name %in% envs)){
      #reticulate::install_python(version = '3.7.9')
      reticulate::virtualenv_create(env.name)
    }

    reticulate::use_virtualenv(env.name)

    #Harder to install ANARCI without a conda environment, will avoid entirely as a temporary solution
    igfold.use.abnum <- T
  }

  if(!reticulate::py_module_available('igfold')){
    system('pip install igfold')
    #message(paste0('Installed IgFold in the ', env.name,' environment'))
  }

  #Fix for pytorch3d not being installed on Mac M1
  if(!reticulate::py_module_available('pytorch3d')){
    system('pip install "git+https://github.com/facebookresearch/pytorch3d.git"')
    #message(paste0('Installed pytorch3d in the ', env.name,' environment'))
  }

  if(!reticulate::py_module_available('igfold')){
    stop('Unable to install IgFold. Please select a conda/virtualenv environment with IgFold already installed. Check https://github.com/Graylab/IgFold for install instructions!')
  }

  #Installing additional packages - currently pyRosetta is not installed due to complicated licensing
  if(!reticulate::py_module_available('openmm')){
    system('pip install openmm')
    #message(paste0('Installed OpenMM in the ', env.name,' environment'))
  }

  if(!reticulate::py_module_available('pdbfixer')){
    system('pip install pdbfixer')
    #message(paste0('Installed PDBfixer in the ', env.name,' environment'))
  }

  sequence.df$pdb_file <- paste0(pdb_dir, '/', sequence.df$fasta_name, '.pdb')

  igfold <- reticulate::import('igfold')
  model <- igfold$IgFoldRunner()

  sequence_df <- sequence.df %>%
                 dplyr::distinct(fasta_sequence, .keep_all = T)

  #Currently not useful to parallelize in R
  for(i in 1:nrow(sequence_df)){
    skip_to_next <- FALSE

    seq_dict <- sequence_df$fasta_sequence[i]

    seq_dict <- stringr::str_replace(seq_dict, '\\*', '')
    if(stringr::str_detect(seq_dict, ':')){
      seq_dict <- as.list(stringr::str_split(seq_dict, ':')[[1]])
      names(seq_dict) <- c('H', 'L')
    }else{
      seq_dict <- list(seq_dict)
      names(seq_dict) <- 'H'
    }
    seq_dict <- reticulate::dict(seq_dict)
    pdb_file <- sequence.df$pdb_file[i]

    tryCatch(model$fold(pdb_file,
                        sequences = seq_dict,
                        do_refine = reticulate::r_to_py(igfold.refine),
                        use_openmm = reticulate::r_to_py(igfold.use.openmm),
                        do_renum = reticulate::r_to_py(igfold.renum),
                        use_abnum = reticulate::r_to_py(igfold.use.abnum)),  error = function(e) { skip_to_next <<- TRUE})

    if(skip_to_next) { next }

    message(paste0('Finished folding sequence ', i, ' out of ', nrow(sequence_df)))
  }

  return(list(pdb_dir = pdb_dir, sequence_df = sequence.df))
}


#' Folds receptor sequences from a given dataframe using OmegaFold
#' @description Folds receptor sequences from a given dataframe using OmegaFold
#' @param sequence.df dataframe - sequences and their IDs, as obtained from the Steropodon_model function.
#' @param model.folder string - path to the OmegaFold model directory.
#' @param omegafold.num.cycle integer - number of recycles during inference.
#' @param omegafold.subbatch.size integer - subbatch size.
#' @param omegafold.allow.tf32 boolean - if TRUE, will default to the TensorFloat-32 (TF32) precision format.
#' @param env.name  string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.
#' @param gpu boolean - if TRUE, will use GPUs for faster structural inference.

#' @return input dataframe with a new column with the folded structure path. Will be processed by Steropodon_model into a Steropodon object/list of objects.
#' @export
#' @examples
#' \dontrun{
#' folded_df <- call_omegafold(sequence_df, model.folder = './omegafold',
#' env.name = 'omegafold',
#' use.conda = T)
#'}
call_omegafold <- function(sequence.df,
                           model.folder,
                           omegafold.num.cycle,
                           omegafold.subbatch.size,
                           omegafold.allow.tf32,
                           env.name,
                           use.conda,
                           gpu,
                           ...){

 if(missing(sequence.df)) stop('Please input the processed sequence dataframe!')
 if(missing(model.folder) | is.null(model.folder)) stop('Please input the OmegaFold directory path!')
 if(missing(omegafold.num.cycle)) omegafold.num.cycle <- 1
 if(missing(omegafold.subbatch.size)) omegafold.subbatch.size <- NULL
 if(missing(omegafold.allow.tf32)) omegafold.allow.tf32 <- T
 if(missing(env.name)) env.name <- NULL
 if(missing(use.conda)) use.conda <- T
 if(missing(gpu)) gpu <- F

 if(omegafold.allow.tf32){
   omegafold.allow.tf32 <- 'True'
 }else{
   omegafold.allow.tf32 <- 'False'
 }

 #Fix for M1 Macbooks for MPS usage
 if(gpu==T & system('uname -m', intern = T) == 'arm64'){
   device <- 'mps'
 }else if(gpu==T){
   device <- 'gpu'
 }else if(gpu==F){
   device <- 'cpu'
 }else{
   device <- 'xla' #Test for TPUs on Colab
 }

 fasta_sequence <- NULL

 pdb_dir <- './tempdir_steropodon_pdbs_omegafold'
 if(!dir.exists(pdb_dir)) dir.create(pdb_dir)
 #Setup the reticulate conda environment
 #options(reticulate.repl.quiet = TRUE)
 if(is.null(env.name)){
   env.name <- 'omegafold'
 }

 if(use.conda){

   envs <- reticulate::conda_list()
   if(!(env.name %in% envs$name)){
     reticulate::conda_create(env.name)
   }
   reticulate::use_condaenv(env.name)

 }else{
   envs <- reticulate::virtualenv_list()
   if(!(env.name %in% envs)){
     reticulate::virtualenv_create(env.name)
   }

   reticulate::use_virtualenv(env.name)
 }


 #Install OmegaFold requirements
 system(paste0('pip install -r ', model.folder, '/requirements.txt'))

 #Reinstall PyTorch and BioPython in case not installed properly
 if(!reticulate::py_module_available('torch')){
   system('pip install torch torchvision torchaudio')
   #message(paste0('Installed PyTorch in the ', env.name,' environment'))
 }

 if(!reticulate::py_module_available('biopython')){
   system('pip install biopython')
   #message(paste0('Installed BioPython in the ', env.name,' environment'))
 }

 sequence.df$pdb_file <- paste0(pdb_dir, '/', sequence.df$fasta_name, '.pdb')

 sequence_df <- sequence.df %>%
               dplyr::distinct(fasta_sequence, .keep_all = T)

 sequence_df$fasta_sequence <- stringr::str_replace_all(sequence_df$fasta_sequence, '\\*', '')
 sequence_df$fasta_sequence <- stringr::str_replace_all(sequence_df$fasta_sequence, ':', '') #Only single-chains are supported

 seq_file <- paste0(pdb_dir,'/all_sequences.fasta')
 seqinr::write.fasta(as.list(sequence_df$fasta_sequence), as.list(sequence_df$fasta_name), seq_file)

 command <- paste0('python ', model.folder, '/main.py', ' --num_cycle ', omegafold.num.cycle, ' --device ', device, ' --allow_tf32 ', omegafold.allow.tf32)
 if(!is.null(omegafold.subbatch.size)){
     command <- paste0(command, ' --subbatch_size ', omegafold.subbatch.size)
 }
 command <- paste0(command, ' ', '"', seq_file, '"', ' ', pdb_dir)
 system(command)

 #Might change to only predict from a single fasta w all sequences (avoids reloading the weights)
 #for(i in 1:nrow(sequence_df)){
 #  skip_to_next <- FALSE
 #  seq_dict <- sequence_df$fasta_sequence[i]

 #  seq_dict <- stringr::str_replace(seq_dict, '\\*', '')
 #  seq_dict <- stringr::str_replace(seq_dict, ':', '')

 #  seq_file <- paste0(pdb_dir, '/', sequence_df$fasta_name[i], '.fasta')

 #  seqinr::write.fasta(seq_dict, sequence_df$fasta_name[i], seq_file)

 #  command <- paste0('python ', model.folder, '/main.py', ' --num_cycle ', omegafold.num.cycle, ' --device ', device, ' --allow_tf32 ', omegafold.allow.tf32)
 #  if(!is.null(omegafold.subbatch.size)){
 #    command <- paste0(command, ' --subbatch_size ', omegafold.subbatch.size)
 #  }

 #  command <- paste0(command, ' ', '"', seq_file, '"', ' ', pdb_dir)

 #  tryCatch(system(command), error = function(e) { skip_to_next <<- TRUE})
 #  if(skip_to_next) { next }

 #  message(paste0('Finished folding sequence ', i, ' out of ', nrow(sequence_df)))
 #}

 return(list(pdb_dir = pdb_dir, sequence_df = sequence.df))
}


#' Folds receptor sequences from a given dataframe using DeepAb
#' @description Folds receptor sequences from a given dataframe using DeepAb
#' @param sequence.df dataframe - sequences and their IDs, as obtained from the Steropodon_model function.
#' @param model.folder string - path to the DeepAb model directory.
#' @param deepab.decoys see DeepAb for more information.
#' @param env.name  string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.
#' @param path.to.pyrosetta string - path to the PyRosetta licence.

#' @return input dataframe with a new column with the folded structure path. Will be processed by Steropodon_model into a Steropodon object/list of objects.
#' @export
#' @examples
#' \dontrun{
#' folded_df <- call_deepab(sequence_df, model.folder = './deepab',
#' env.name = 'deepab',
#' use.conda = T)
#'}
call_deepab <- function(sequence.df,
                        model.folder,
                        deepab.decoys,
                        env.name,
                        use.conda,
                        path.to.pyrosetta,
                        ...
                       ){

   if(missing(sequence.df)) stop('Please input the processed sequence dataframe!')
   if(missing(model.folder) | is.null(model.folder)) stop('Please input the DeepAb folder path!')
   if(missing(deepab.decoys)) deepab.decoys <- 5
   if(missing(env.name)) env.name <- NULL
   if(missing(path.to.pyrosetta)) path.to.pyrosetta <- NULL

   fasta_sequence <- NULL

   pdb_dir <- './tempdir_steropodon_pdbs_deepab'
   if(!dir.exists(pdb_dir)) dir.create(pdb_dir)
   #Setup the reticulate conda environment
   #options(reticulate.repl.quiet = TRUE)

   if(is.null(env.name)){
     env.name <- 'deepab'
   }

   if(use.conda){
     envs <- reticulate::conda_list()

     if(!(env.name %in% envs$name)){
       reticulate::install_python(version = '3.7.9')
       reticulate::use_python_version(version = '3.7.9')
       reticulate::conda_create(env.name, python_version = '3.7.9')
     }

     if(!reticulate::py_module_available('pyrosetta')){
       reticulate::py_install(envname = env.name, packages = 'pyrosetta', python_version = '3.7.9', method = 'conda')
     }
     reticulate::use_condaenv(env.name)

   }else{
     envs <- reticulate::virtualenv_list()
     if(!(env.name %in% envs)){
       reticulate::install_python(version = '3.7.9')
       reticulate::use_python_version(version = '3.7.9')
       reticulate::virtualenv_create(env.name, version = '3.7.9')
     }
     reticulate::use_virtualenv(env.name)

     if(!reticulate::py_module_available('pyrosetta')){
       system(
         paste0(
           'pip install --e ', path.to.pyrosetta, '/setup/'
         )
       )
     }

   }

   system(
     paste0(
       'pip install -r ', model.folder, '/', 'requirements.txt'
     )
   )

   models_dir <- paste0('./trained_models')
   if(!dir.exists(models_dir)) dir.create(models_dir)

   weights_dir <- paste0(models_dir, '/ensemble_abresnet')
   if(!dir.exists(weights_dir)) dir.create(weights_dir)

   if(length(list.files(weights_dir))==0){
     if(!('ensemble_abresnet_v1.tar.gz') %in% list.files()){
       system('wget https://data.graylab.jhu.edu/ensemble_abresnet_v1.tar.gz')
     }
     system('tar -xf ensemble_abresnet_v1.tar.gz')
     system(
       paste0(
       'mv ./ensemble_abresnet ', models_dir
      )
     )

     system('rm -r ensemble_abresnet_v1.tar.gz')
     message('Downloaded DeepAb weights!')
   }

   sequence.df$pdb_file <- paste0(pdb_dir, '/', sequence.df$fasta_name, '.pdb')

   sequence_df <- sequence.df %>%
                  dplyr::distinct(fasta_sequence, .keep_all = T)


   for(i in 1:nrow(sequence_df)){
      skip_to_next <- FALSE
      seq_dict <- sequence_df$fasta_sequence[i]

      seq_dict <- stringr::str_replace(seq_dict, '\\*', '')
      if(stringr::str_detect(seq_dict, ':')){
        seq_dict <- as.list(stringr::str_split(seq_dict, ':')[[1]])
        names(seq_dict) <- c(':H', ':L')
      }else{
        seq_dict <- list(seq_dict)
        names(seq_dict) <- ':H'
      }

      seq_file <- paste0(pdb_dir, '/', sequence_df$fasta_name[i], '.fasta')
      seqinr::write.fasta(seq_dict, names(seq_dict), seq_file)

      command <- paste0('python ', model.folder, '/predict.py', ' ', seq_file, ' --decoys ', deepab.decoys)
      if(length(seq_dict) == 1){
        command <- paste0(command, ' --single_chain')
      }

      tryCatch(system(command), error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { next }

      out_dir <- list.files()[which(stringr::str_detect(list.files(), 'pred_'))]
      file.copy(paste0(out_dir, '/pred.deepab.pdb'), pdb_dir)
      file.rename(paste0(pdb_dir, '/pred.deepab.pdb'), sequence_df$pdb_file[i])
      system(paste0('rm -r ', out_dir))

      message(paste0('Finished folding sequence ', i, ' out of ', nrow(sequence_df)))
   }

   return(list(pdb_dir = pdb_dir, sequence_df = sequence.df))
}


#' Calls all implemented models inside Steropodon_model
#' @description Calls all implemented models inside Steropodon_model
#' @param sequence.df dataframe - sequences and their IDs, as obtained from the Steropodon_model function.
#' @param model string - folding model ('colabfold', 'omegafold', 'igfold', 'deepab').
#' @param model.folder boolean - if TRUE, will refine the modelled structure using OpenMM or PyRosetta.
#' @param env.name  string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.
#' @param gpu boolean - if TRUE, will use GPUs for faster structural inference.
#' @param additional.model.parameters named list - additional parameters specific to the folding model.

#' @return input dataframe with a new column with the folded structure path. Will be processed by Steropodon_model into a Steropodon object/list of objects.
#' @export
#' @examples
#' \dontrun{
#' folded_df <- call_models(sequence_df, model = 'igfold',
#' env.name = 'igfold',
#' use.conda = T)
#'}
call_models <- function(sequence.df, model, model.folder, env.name, use.conda, gpu, additional.model.parameters){
  if(model == 'colabfold'){
    params <- list(sequence.df = sequence.df,
                   env.name = env.name,
                   use.conda = use.conda,
                   gpu = gpu)

    out <- do.call(call_colabfold, c(params, additional.model.parameters))

  }else if(model == 'omegafold'){
    params <- list(sequence.df = sequence.df,
                model.folder = model.folder,
                env.name = env.name,
                use.conda = use.conda,
                gpu = gpu)

    out <- do.call(call_omegafold, c(params, additional.model.parameters))

  }else if(model == 'igfold'){
    params <- list(sequence.df = sequence.df,
                   env.name = env.name,
                   use.conda = use.conda
                   )

    out <- do.call(call_igfold, c(params, additional.model.parameters))


  }else if(model == 'deepab'){
    params <- list(sequence.df = sequence.df,
                   model.folder = model.folder,
                   env.name = env.name,
                   use.conda = use.conda
                   )

    out <- do.call(call_deepab, c(params, additional.model.parameters))


  }else{
    stop ('Model not implemented! Check the list of pLM/MSA-based protein folding models available!')
  }

  message(paste0('Finished folding using ', model))
  return(out)
}
