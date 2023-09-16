#' Main function of the Steropodon module - obtain a nested list of Steropodon objects, containing folded receptor structures


#' @description Main function of the Steropodon module - obtain a nested list of Steropodon objects, containing folded receptor structures from the sequences initially selected.
#' Selection of the unique sequences to be folded can occur either per clonotype (given a threshold), from a vector of sequences, or from a vector of unique cell barcodes.
#' Steropodon_model employs several deep-learning models for protein folding, mostly relying on alignment/MSA-free models for a lighter local setup and, usually, faster inference.
#' For example, IgFold, OmegaFold, and DeepAb make use of protein language model-based (pLM) architectures, completely avoiding the MSA step used in models such as AlphaFold, whereas ColabFold uses a public server for MSA queries.
#' These can be installed directly on a local machine using Reticulate and called with a minimal download requirement.
#' Moreover, by specifying an environment name for either conda or Python’s virtualenv, Steropodon will install all dependencies of a given model when first calling Steropodon_model.
#' If this fails, we recommend following the installation instructions for your desired protein folding model and then using that environment as input to env.name.
#' For ColabFold, we recommend first installing following these instructions https://github.com/YoshitakaMo/localcolabfold, as Steropodon_model requires the colabfold_batch command to be exported to your path variable.
#' Therefore, we have highlighted the main design philosophy of the Steropodon pipeline: it can be directly used with a single call on a local machine - it does not require extensive installations of either libraries or sequence databases.
#' This enables the fast exploratory analysis of a given repertoire’s structural features and, most importantly, allows researchers to employ this pipeline on their local machines.
#' Of course, this comes with a computational cost, and the number of structures that can be folded in a reliable amount of time depends on each individual machine.
#' Steropodon can also be set up on a cluster, while the gpu parameter in Steropodon_model will enable GPU usage and faster inference.
#' Steropodon_model allows the folding of either single chains ('VDJ' or 'VJ' in the sequence.type parameter) or paired chains (sequence.type = 'VDJ.VJ').

#' @param VGM the VGM object, used in most analyses in the Platypus computational immunology ecosystem. Obtained from the VDJ_GEX_matrix function.
#' @param AntibodyForests.object the AntibodyForests object - a lineage of set of clonal lineages with unique sequences to be extracted for structural modelling. More multimodal analyses (structural evolution and clonal networks) will be implemented soon.
#' @param model string - the immune receptor folding model to be used. Current options include: 'igfold', 'omegafold', 'colabfold', and 'deepab'. Other deep learning models for protein structure inference will be added soon.
#' @param model.folder string - the path to the model directory (for weights, inference scripts, etc). See Steropodon_protein_folding for the subroutine for each folding model and if/how these paths are needed (as Steropodon aims to initially install the model in a specified environment and call from Python using Reticulate).
#' @param additional.model.parameters named list - additional parameters/options for the folding models. See Steropodon_protein_folding for more information.
#' @param sequence.type string - the chains that should be folded: 'VDJ.VJ' for paired heavy-light chains, 'VDJ' for heavy, 'VJ' for light.
#' @param sequence.vector vector of strings or NULL - vector of receptor sequences that should be folded (will no longer select per clonotype from the VGM object).
#' @param barcode.vector vector of strings or NULL - vector of cell barcodes (from the GEX/Seurat object) for which unique receptors should be selected for structural modelling.
#' @param antigen.columns vector of strings or NULL - if the folding method selected allows modelling antibody-antigen complexes (e.g., model = 'colabfold'), this should define the VGM column name(s) with the antigen sequences.
#' @param antigen.pdbs vector of strings or NULL - if the folding method selected allows modelling antibody-antigen complexes (e.g., model = 'colabfold'), this should define the PDB IDs for the antigen further used in multimer/complex folding.
#' @param antigen.sequences vector of strings or NULL - if the folding method selected allows modelling antibody-antigen complexes (e.g., model = 'colabfold'), this should define the antigen sequences further used in multimer/complex folding.
#' @param max.clonotypes integer - the maximum number of clonotypes from which unique sequences will be picked for folding, per sample. Clonotype will first be ranked by clonal expansion.
#' @param max.per.clonotype integer - the maximum number of unique sequences per clonotype selected for structural prediction.
#' @param min.clonotype.frequency integer - the minimum clonal expansion for a clonotype to be selected from for folding.
#' @param save.pdbs boolean - if TRUE, all structures will be also saved as PDB file.
#' @param save.rds boolean - if TRUE, will save the Steropodon object as an RDS file.
#' @param save.dir string - path to the directory for saving structures/Steropodon objects.
#' @param env.name string - the conda environment name with the model and all dependencies already installed or to be installed by Steropodon_model.
#' @param use.conda boolean - if TRUE, will use conda for managing the Python environments and installing folding model dependencies. Else, it will use virtualenv.
#' @param gpu boolean - if TRUE, will use GPUs for faster structural inference.

#' @return a nested list of Steropodon objects (per sample, per clonotype), with slots for the modelled structure and several other slots for saving structures from the downstream analyses.
#' @export
#' @examples
#' \dontrun{
#' Steropodon_model(VDJ,
#'                 model = 'igfold',
#'                 sequence.type = 'VDJ.VJ',
#'                 max.clonotypes = 10,
#'                 max.per.clonotype = 1,
#'                 additional.model.parameters = list(igfold.refine = T, igfold.use.openmm = T),
#'                 save.rds = T,
#'                 save.dir = './steropodon_RDS',
#'                 use.conda = T,
#'                 env.name = 'vignette_env'
#'}


Steropodon_model <- function(VGM,
                             AntibodyForests.object,
                             model,
                             model.folder,
                             additional.model.parameters,
                             sequence.type,
                             sequence.vector,
                             barcode.vector,
                             antigen.columns,
                             antigen.pdbs,
                             antigen.sequences,
                             max.clonotypes,
                             max.per.clonotype,
                             min.clonotype.frequency,
                             save.pdbs,
                             save.rds,
                             save.dir,
                             env.name,
                             use.conda,
                             gpu
                             ){

  if(missing(VGM)) stop('Please input your VGM object - the output of the VDJ_GEX_matrix function')
  if(missing(AntibodyForests.object)) AntibodyForests.object <- NULL
  if(missing(model)) model <- 'igfold'
  if(missing(model.folder)) model.folder <- '/Users/tudorcotet/OmegaFold'
  if(missing(additional.model.parameters)) additional.model.parameters <- list()
  if(missing(sequence.type)) sequence.type  <- 'VDJ.VJ'
  if(missing(sequence.vector)) sequence.vector <- NULL
  if(missing(barcode.vector)) barcode.vector <- NULL
  if(missing(antigen.columns)) antigen.columns <- NULL
  if(missing(antigen.pdbs)) antigen.pdbs <- NULL
  if(missing(antigen.sequences)) antigen.sequences <- NULL
  if(missing(max.clonotypes)) max.clonotypes <- NULL
  if(missing(max.per.clonotype)) max.per.clonotype <- NULL
  if(missing(min.clonotype.frequency)) min.clonotype.frequency <- NULL
  if(missing(save.pdbs)) save.pdbs <- F
  if(missing(save.rds)) save.rds <- T
  if(missing(save.dir)) save.dir <- './steropodon_annotated'
  if(missing(env.name)) env.name <- NULL
  if(missing(use.conda)) use.conda <- T
  if(missing(gpu)) gpu <- F

  sample_id <- NULL
  clonotype_id <- NULL
  sequence_column <- NULL
  clonotype_frequency <- NULL
  sequence_frequency <- NULL
  fasta_sequence <- NULL
  n_cells <- NULL


  #Define global variables
  aa_short <- c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'Y')
  aa_abbrev <- c('ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PYL', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'SEC', 'VAL', 'TRP', 'TYR')

  abbrev_to_short <- aa_short
  names(abbrev_to_short) <- aa_abbrev


  short_to_abbrev <- aa_abbrev
  names(short_to_abbrev) <- aa_short

  VDJ_regions <- c('VDJ_aaSeqFR1', 'VDJ_aaSeqCDR1', 'VDJ_aaSeqFR2', 'VDJ_aaSeqCDR2', 'VDJ_aaSeqFR3', 'VDJ_aaSeqCDR3', 'VDJ_aaSeqFR4')
  VJ_regions <- c('VJ_aaSeqFR1', 'VJ_aaSeqCDR1', 'VJ_aaSeqFR2', 'VJ_aaSeqCDR2', 'VJ_aaSeqFR3', 'VJ_aaSeqCDR3', 'VJ_aaSeqFR4')


  extract_MIXCR <- function(VDJ.matrix, chain.to.extract, as.nucleotide, regions.to.extract){
      if(missing(VDJ.matrix)) stop('Input the VDJ dataframe obtained after calling VDJ_call_MIXCR')
      if(missing(chain.to.extract)) chain.to.extract <- 'VDJ.VJ'
      if(missing(as.nucleotide)) as.nucleotide <- F
      if(missing(regions.to.extract)) regions.to.extract <- list('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')


      VDJ.matrix$mixcr_assembled <- rep('', nrow(VDJ.matrix))

      if(chain.to.extract!='VDJ.VJ'){
        if(as.nucleotide){
          col_name <- paste0(chain.to.extract, '_nSeq')
        }else{
          col_name <- paste0(chain.to.extract, '_aaSeq')
        }

        for(region in regions.to.extract){
          VDJ.matrix$mixcr_assembled <- paste0(VDJ.matrix$mixcr_assembled, gsub('_', '', VDJ.matrix[, paste0(col_name, region)]))
        }
      }else if(chain.to.extract=='VDJ.VJ'){
        if(as.nucleotide==T){
          col_name <- '_nSeq'
        }else{
          col_name <- '_aaSeq'
        }
        extracted_VDJ <- rep('', nrow(VDJ.matrix))
        extracted_VJ <- rep('', nrow(VDJ.matrix))

        for(region in regions.to.extract){
          extracted_VDJ <- paste0(extracted_VDJ, gsub('_', '', VDJ.matrix[,paste0('VDJ',col_name, region)]))
          extracted_VJ <- paste0(extracted_VJ, gsub('_', '', VDJ.matrix[,paste0('VJ',col_name, region)]))
        }

        VDJ.matrix$mixcr_assembled <- paste0(extracted_VDJ, ':', extracted_VJ)

      }
      return(VDJ.matrix)
  }

  get_unique_sequences <- function(VDJ,
                                   sequence.type,
                                   sequence.vector,
                                   barcode.vector,
                                   antigen.columns,
                                   antigen.pdbs,
                                   antigen.sequences,
                                   max.clonotypes,
                                   max.per.clonotype,
                                   min.clonotype.frequency){


    #FIRST: remove all aberrants!
    VDJ <- VDJ[VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains == 1,]

    #Also remove sequences with missing regions
    if(sequence.type == 'VDJ' | sequence.type == 'VDJ.VJ'){
      for(i in 1:length(VDJ_regions)){
        VDJ <- VDJ[which(VDJ[[VDJ_regions[i]]]!=''),]
        VDJ[[VDJ_regions[i]]] <- stringr::str_replace_all(VDJ[[VDJ_regions[i]]], '_', '')
      }
    }

    if(sequence.type == 'VJ' | sequence.type == 'VDJ.VJ'){
      for(i in 1:length(VJ_regions)){
        VDJ <- VDJ[which(VDJ[[VJ_regions[i]]]!=''),]
        VDJ[[VJ_regions[i]]] <- stringr::str_replace_all(VDJ[[VJ_regions[i]]], '_', '')
      }
    }

    #Currently only obtaining sequences from MIXCR (vs trim and align/ other methods)
    VDJ <- extract_MIXCR(VDJ, chain.to.extract = sequence.type)

    #Can be modified in the future for other sequence types
    VDJ$sequence_column <- VDJ$mixcr_assembled

    VDJ <- VDJ %>%
           dplyr::group_by(sample_id, clonotype_id) %>%
           dplyr::mutate(clonotype_frequency = dplyr::n()) %>%
           dplyr::ungroup()

    VDJ <- VDJ %>%
           dplyr::group_by(sample_id, clonotype_id, sequence_column) %>%
           dplyr::mutate(sequence_frequency = dplyr::n()) %>%
           dplyr::ungroup()

    if(!is.null(sequence.vector) | !is.null(barcode.vector)){

      if(!is.null(barcode.vector)){
        sequence.vector <- unique(VDJ$sequence_column[VDJ$barcode %in% barcode.vector,])
      }


    }else{
      if(!is.null(min.clonotype.frequency)){
        VDJ <- VDJ[VDJ$clonotype_frequency > min.clonotype.frequency,]
        sequence.vector <- unique(VDJ$sequence_column)
      }

      if(!is.null(max.clonotypes)){
        ids <- unique(VDJ$sample_id)
        per_sample_dfs <- list()

        for(i in 1:length(ids)){
          selected_clonotypes <- VDJ[VDJ$sample_id == ids[i],] %>%
                                 dplyr::distinct(clonotype_id, .keep_all = T) %>%
                                 dplyr::arrange(plyr::desc(clonotype_frequency))

          selected_clonotypes <- selected_clonotypes$clonotype_id[1:max.clonotypes]
          per_sample_dfs[[i]] <- VDJ[VDJ$sample_id == ids[i] & (VDJ$clonotype_id %in% selected_clonotypes),]
        }
        VDJ <- do.call('rbind', per_sample_dfs)

        sequence.vector <- unique(VDJ$sequence_column)
      }

      if(!is.null(max.per.clonotype)){
        temp_VDJ <- VDJ %>%
               dplyr::group_by(sample_id, clonotype_id, sequence_column) %>%
               dplyr::arrange(plyr::desc(sequence_frequency))

        temp_VDJ <- temp_VDJ %>%
                dplyr::group_by(sample_id, clonotype_id) %>%
                dplyr::distinct(sequence_column, .keep_all = T) %>%
                dplyr::slice_head(n = max.per.clonotype)

        sequence.vector <- unique(temp_VDJ$sequence_column)
      }
    }

    final_dfs <- list()
    for(i in 1:length(sequence.vector)){
      sample_ids <- VDJ$sample_id[VDJ$sequence_column == sequence.vector[i]]
      clonotype_ids <- VDJ$clonotype_id[VDJ$sequence_column == sequence.vector[i]]
      barcodes <- VDJ$barcode[VDJ$sequence_column == sequence.vector[i]]
      fasta_sequences <- rep(sequence.vector[i], length(barcodes))
      unique_samples <- paste0(unique(sample_ids), collapse = ";")
      unique_clonotypes <- paste0(unique(clonotype_ids), collapse = ';')
      fasta_names <- paste0(c(unique_samples, unique_clonotypes, i), collapse = '_')

      final_dfs[[i]] <- data.frame(sample_id = unlist(sample_ids),
                                   clonotype_id = unlist(clonotype_ids),
                                   fasta_sequence = unlist(fasta_sequences),
                                   barcodes = unlist(barcodes),
                                   fasta_name = rep(fasta_names, length(barcodes)),
                                   sequence_type = rep(sequence.type, length(barcodes))
                                   )
    }


    final_df <- do.call('rbind', final_dfs)

    #Add VDJ/VJ regions inferred from MIXCR to annotate the structures later
    if(sequence.type == 'VDJ' | sequence.type == 'VDJ.VJ'){
      VDJ_subset <- VDJ %>% dplyr::select(c('barcode', VDJ_regions))
      final_df <- merge(final_df, VDJ_subset, by.x = 'barcodes', by.y = 'barcode', all.x = T)
    }

    if(sequence.type == 'VJ' | sequence.type == 'VDJ.VJ'){
      VDJ_subset <- VDJ %>% dplyr::select(c('barcode', VJ_regions))
      final_df <- merge(final_df, VDJ_subset, by.x = 'barcodes', by.y = 'barcode', all.x = T)
    }

    final_df <- final_df %>%
      dplyr::arrange(nchar(sample_id), sample_id, nchar(clonotype_id), clonotype_id)


    if(!is.null(antigen.columns) | !is.null(antigen.sequences) | !is.null(antigen.pdbs)){

      if(!is.null(antigen.columns)){

        antigen_sequences <- c()
        for(i in 1:length(antigen.columns)){
          antigen_sequences[i] <- unique(VDJ[antigen.columns[i],])[1]
        }
        names(antigen_sequences) <- antigen.columns

      }else if(!is.null(antigen.pdbs)){

        antigen_sequences <- c()
        for(i in 1:length(antigen.pdbs)){
          pdb <- bio3d::read.pdb(antigen.pdbs[i])
          antigen_sequences[i] <- paste0(abbrev_to_short[unname(pdb$seqres)], collapse = '')
        }
        names(antigen_sequences) <- antigen.pdbs

      }else if(!is.null(antigen.sequences)){

        antigen_sequences <- antigen.sequences
        if(is.null(names(antigen_sequences))){
          names(antigen_sequences) <- paste0(rep('antigen', length(antigen_sequences)), '_', 1:length(antigen_sequences))
        }
      }

      final_dfs <- list()
      for(i in 1:length(antigen_sequences)){
        final_df$antigen_sequence <- antigen_sequences[i]
        final_df$antigen_name <- names(antigen_sequences)[i]
        final_df$fasta_sequence <- paste0(final_df$fasta_sequence, ':', antigen_sequences[i])
        final_df$fasta_name <- paste0(final_df$fasta_name, '_', names(antigen_sequences)[i])
        final_dfs[[i]] <- final_df
      }

      final_df <- do.call('rbind', final_dfs)
    }

    return(final_df)
  }

  create_fasta_directory <- function(sequence.df,
                                     model
                                     ){

    if(missing(sequence.df)) stop('Please input your sequence dataframe!')
    if(missing(model)) model <- 'colabfold'

    temp_dir <- './tempdir_steropodon_fasta'
    if(!dir.exists(temp_dir)) dir.create(temp_dir)


    sequence_df <- sequence.df %>%
                   dplyr::distinct(fasta_sequence, .keep_all = T)


    if(model == 'colabfold'){
      sequence_df$file_name <- paste0(temp_dir, '/', sequence_df$fasta_name, '.fasta')

      sequences <- sequence_df$fasta_sequence
      fasta_names <- sequence_df$fasta_name
      file_names <- sequence_df$file_name

      for(i in 1:length(sequences)){
        seqinr::write.fasta(sequences[i], fasta_names[i], file_names[i])
      }

    }
  }

  create_single_fasta <- function(sequence.df
                                  ){

    if(missing(sequence.df)) stop('Please input your sequence dataframe!')

    temp_dir <- './tempdir_steropodon_fasta'
    if(!dir.exists(temp_dir)) dir.create(temp_dir)


    sequence_df <- sequence.df %>%
                   dplyr::distinct(fasta_sequence, .keep_all = T)



   sequences <- sequence_df$fasta_sequence
   fasta_names <- sequence_df$fasta_name

   seqinr::write.fasta(sequences, fasta_names, 'sequences.fasta')
  }

  get_region_idx <- function(seq.df){
    seq_df <- seq.df %>%
              dplyr::distinct(fasta_sequence, .keep_all = T)

    sequence_type <- unique(seq_df$sequence_type)

    region_idx <- c()
    chain_idx <- c()

    if(sequence_type == 'VDJ' | sequence_type == 'VDJ.VJ'){
      regions <- c(seq_df[VDJ_regions][1,])
      names(regions) <- c('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')
      regions <- unlist(lapply(1:length(regions), function(i) rep(names(regions)[i], nchar(regions[i]))))
      chains <- rep('VDJ', length(regions))

      region_idx <- c(region_idx, regions)
      chain_idx <- c(chain_idx, chains)
    }

    if(sequence_type == 'VJ' | sequence_type == 'VDJ.VJ'){
      regions <- c(seq_df[VJ_regions][1,])
      names(regions) <- c('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')
      regions <- unlist(lapply(1:length(regions), function(i) rep(names(regions)[i], nchar(regions[i]))))
      chains <- rep('VJ', length(regions))

      region_idx <- c(region_idx, regions)
      chain_idx <- c(chain_idx, chains)
    }

    names(region_idx) <- 0:(length(region_idx) - 1)
    names(chain_idx) <- 0:(length(chain_idx) - 1)

    return(list(chain_idx = chain_idx, region_idx = region_idx))
  }

  create_steropodon_object <- function(file, seq.df){
    if(!(file %in% seq.df$pdb_file) | !(file.exists(file))){
      return(NULL)
    }


    pdb <- bio3d::read.pdb(file)

    seq_df <- seq.df[seq.df$pdb_file==file,]
    barcodes <- unique(seq_df$barcode)
    sequence_type <- unique(seq_df$sequence_type)

    idx_list <- get_region_idx(seq.df)
    region_idx <- idx_list$region_idx
    chain_idx <- idx_list$chain_idx

    pdb <- annotate_structure(pdb, res.idx = region_idx, feature = 'region')
    pdb <- annotate_structure(pdb, res.idx = chain_idx, feature = 'chain')

    VDJ_chain <- NULL
    VDJ_cdr3 <- NULL
    VJ_chain <- NULL
    VJ_cdr3 <- NULL
    complex <- NULL
    antigen_struct <- NULL

    if('antigen_name' %in% colnames(seq_df)){
      pdb$atom$chain[is.na(pdb$atom$chain)] <- 'antigen'
      pdb$atom$region[is.na(pdb$atom$region)] <- 'antigen'

      antigen_struct <- split_structure(pdb, grouping = c('chain'), specific.values = 'antigen')
      complex <- pdb
      pdb <- split_structure(pdb, grouping = c('chain'), specific.values = c('VDJ','VJ'), combine.values = T, combine.groupings = F)
    }

    if(sequence_type == 'VDJ' | sequence_type == 'VDJ.VJ'){
      VDJ_chain <- split_structure(pdb, grouping = 'chain', specific.values = 'VDJ', combine.values = F, combine.groupings = F)
      VDJ_cdr3 <- split_structure(pdb, grouping = c('chain', 'region'), specific.values = 'VDJ_CDR3', combine.values = F, combine.groupings = T)
    }

    if(sequence_type == 'VJ' | sequence_type == 'VDJ.VJ'){
      VJ_chain <- split_structure(pdb, grouping = 'chain', specific.values = 'VJ', combine.values = F, combine.groupings = F)
      VJ_cdr3 <- split_structure(pdb, grouping = c('chain', 'region'), specific.values = 'VJ_CDR3', combine.values = F, combine.groupings = T)
    }
    obj <- methods::new('Steropodon',
                        structure = pdb,
                        sequence = paste0(bio3d::pdbseq(pdb), collapse = ''), #or add bio3d::pdbseq(pdb)
                        barcodes = barcodes,
                        H = VDJ_chain,
                        L = VJ_chain,
                        CDRH3 = VDJ_cdr3,
                        CDRL3 = VJ_cdr3,
                        antigen = antigen_struct,
                        complex = complex,
                        properties = NULL,
                        paratope = NULL,
                        epitope = NULL,
                        core = NULL,
                        pdbs = NULL,
                        structure_id = NULL
                       )
    return(obj)
  }

  parse_pdb_directory <- function(pdb.dir, sequence.df){
    sample_ids <- unique(sequence.df$sample_id)
    pdb_list <- vector(mode = 'list', length = length(sample_ids))

    for(i in 1:length(sample_ids)){
      clonotype_ids <- unique(sequence.df$clonotype_id[sequence.df$sample_id == sample_ids[i]])
      pdb_list[[i]] <- vector(mode = 'list', length = length(clonotype_ids))


      for(j in 1:length(clonotype_ids)){
        sequence_subset <- sequence.df[sequence.df$sample_id == sample_ids[i] & sequence.df$clonotype_id == clonotype_ids[j],]

        temp_subset <- sequence_subset %>%
                           dplyr::group_by(fasta_sequence) %>%
                           dplyr::mutate(n_cells = dplyr::n()) %>%
                           dplyr::ungroup() %>%
                           dplyr::distinct(fasta_sequence, .keep_all = T) %>%
                           dplyr::arrange(dplyr::desc(n_cells))

        files_to_parse <- temp_subset$pdb_file
        pdb_list[[i]][[j]] <- vector(mode = 'list', length = length(files_to_parse))
        ind <- 1
        #Find a way to parallelize
        while(length(files_to_parse) > 0){
          file <- files_to_parse[1]
          pdb_list[[i]][[j]][[ind]] <- create_steropodon_object(file, sequence_subset)
          files_to_parse <- files_to_parse[files_to_parse != file]
          ind <- ind + 1
        }
        if(length(pdb_list[[i]][[j]]) != 0){
          names(pdb_list[[i]][[j]]) <- 1:length(pdb_list[[i]][[j]])
        }
        pdb_list[[i]][[j]] <- pdb_list[[i]][[j]][!sapply(pdb_list[[i]][[j]], is.null)]
      }

      if(length(pdb_list[[i]]) != 0){
        names(pdb_list[[i]]) <- clonotype_ids
      }
      pdb_list[[i]] <- pdb_list[[i]][!sapply(pdb_list[[i]], function(x) length(x) == 0)]
    }
    if(length(pdb_list) != 0){
      names(pdb_list) <- sample_ids
    }
    pdb_list <- pdb_list[!sapply(pdb_list, function(x) length(x) == 0)]

    unlink(pdb.dir, recursive = T)

    return(pdb_list)
  }

  write_steropodon <- function(steropodon.list){

  }

  if(inherits(VGM, 'list')){
    VDJ <- VGM[[1]]
  }else{
    VDJ <- VGM
  }

  VDJ <- extract_MIXCR(VDJ,
                       chain.to.extract = sequence.type,
                       as.nucleotide = F
                       )

  sequence_df <- get_unique_sequences(VDJ = VDJ,
                                      sequence.type = sequence.type,
                                      sequence.vector = sequence.vector,
                                      barcode.vector = barcode.vector,
                                      antigen.columns = antigen.columns,
                                      antigen.pdbs = antigen.pdbs,
                                      antigen.sequences = antigen.sequences,
                                      max.clonotypes = max.clonotypes,
                                      max.per.clonotype = max.per.clonotype,
                                      min.clonotype.frequency = min.clonotype.frequency)

  out_list <- call_models(sequence.df = sequence_df,
                         model = model,
                         model.folder = model.folder,
                         env.name = env.name,
                         use.conda = use.conda,
                         gpu = gpu,
                         additional.model.parameters)

  steropodon_object <- parse_pdb_directory(pdb.dir = out_list$pdb_dir,
                                           sequence.df = out_list$sequence_df)

  #Add structure IDs
  steropodon_object <- unnest_steropodon(steropodon_object)
  ids <- names(steropodon_object)
  steropodon_object <- lapply(1:length(steropodon_object), function(i) {id <- stringr::str_split(ids[i], '\\.')[[1]]
                                                                        steropodon_object[[i]]@structure_id <- list(sample = id[1], clonotype = id[2], rank = id[3])
                                                                        return(steropodon_object[[i]])
                                                                        })
  names(steropodon_object) <- ids
  steropodon_object <- nest_steropodon(steropodon_object)
  if(save.pdbs){
    if(!dir.exists(save.dir)) dir.create(save.dir)
  }

  if(save.rds){
    if(!dir.exists(save.dir)) dir.create(save.dir)
    saveRDS(steropodon_object, file.path(save.dir, paste0(Sys.time(), '_Steropodon_Object.rds')))
  }


  return(steropodon_object)
}
