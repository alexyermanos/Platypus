#' (Re)clonotype a VDJ object using cellranger's enclone tool
#'@description Calls recon to clonotype a VDJ object given a VDJ.directory (with sample folders which should include the all_contig_annotations.json file) - outputs a new VDJ with updated clonotype_id, clonotype_id_10x, and clonotype_frequency columns
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param VDJ.directory string - directory for the VDJ data, should be the main folder which includes the individual sample folders (each with the all_contig_annotations.json file that is used by enclone)
#' @param global.clonotype bool - if T, will use clonotype definitions irrespective of samples. Must also be T is you wish to merge clonotypes from two specific (which should be specified in the samples.to.combine parameter)
#' @param samples.to.clonotype - vector - lists the samples names which should be clonotyped. The unspecified samples will keep their old clonotype defintions.
#' @param samples.to.combine - vector or list of vectors - lists the samples which you wish to have their clonotypes merged (e.g., c('s1','s2') to only merge the first 2 samples, or list(c('s1','s3'), c('s2', 's4')) to merge the first and third, second and fourth, respectively). global.clonotype must be set to T!
#' @param same.origin bool - if the merged samples come from the same donor, with the same or with different origins. If two datasets come from the same origin, enclone will filter to remove certain artifacts.
#' @param output.format string - 'vgm' to output a VGM-specific VDJ dataframe (all samples in the same dataframe).
#' @param operating.system string - operating system on which enclone will be run. 'Windows' for Windows, 'Linux' for Linux, 'Darwin' for MacOS.
#' @param parallel bool - if T, the program will be executed in parallel, on no. cores = max. available cores - 1.

#' @return Reclonotyped VDJ object using the enclone software and 10x-specific clonotype deifinition.
#' @export


#' @examples
#' \donttest{
#' try({
#' VDJ_call_enclone(vdj, VDJ.directory, samples.to.combine = c('s1', 's2', 's3'), global.clonotype = T)
#' })
#'}

VDJ_call_enclone <- function(VDJ,
                             VDJ.directory,
                             global.clonotype,
                             samples.to.clonotype,
                             samples.to.combine,
                             same.origin,
                             output.format,
                             operating.system,
                             parallel){

  if(missing(VDJ)) stop('Please input your VDJ matrix (vgm[[1]]) dat you wish to reclonotype using 10x enclone')
  if(missing(VDJ.directory)) stop('Please input your VDJ directory')
  if(missing(global.clonotype)) global.clonotype <- FALSE
  if(missing(samples.to.clonotype)) samples.to.clonotype <- 'all'
  if(missing(samples.to.combine) & global.clonotype == TRUE) samples.to.combine <- NULL
  if(missing(same.origin)) same.origin <- FALSE
  if(missing(output.format)) output.format <- 'vgm'
  if(missing(parallel)) parallel <- FALSE
  if(missing(operating.system)){
        switch(Sys.info()[['sysname']],
               Windows= {message("Windows system detected")
                         operating.system <- "Windows"},
               Linux  = {message("Linux system detected")
                        operating.system <- "Linux"},
               Darwin = {message("MAC system detected")
                        operating.system <- "Darwin"})

  }

  #SUBROUTINE 1: strip sample id off barcodes for easier mergings

  strip_barcodes <- function(sample_df){
    strip_single_barcode <- function(barcode){
      stripped_barcode <- unlist(stringr::str_split(barcode, '_'))[2]
      return(stripped_barcode)
    }
    sample_df$stripped_barcodes <- unlist(lapply(sample_df$barcode, function(x) strip_single_barcode(x)))
    return(sample_df)
  }

  #SUBROUTINE 2: calls enclone on sample df, returns the output csv path, output sample df, and sample name(s)

  get_enclone_output <- function(sample_df, output_dir = temp_dir, VDJ_directory = VDJ.directory){
    sample_names <- unlist(unique(sample_df$sample_id))

    #sample_names <- paste0(sample_names, collapse = merge_origin_sep)
    file_name <- paste0(sample_names, collapse = '_')
    temp_file <- paste0(output_dir, '/', file_name, '.csv')


    VDJ_sample_paths <- c()
    for(i in 1:length(sample_names)){
      VDJ_sample_paths[i] <- paste0(VDJ_directory, '/', sample_names[i])
    }

    path_with_origin <- paste0(VDJ_sample_paths, collapse = merge_origin_sep)

    #Get cell type
    cell_type <- unique(sample_df$celltype)
    if(length(cell_type) > 1){
      stop('Enclone only works on data with a single cell type - ensure your (merged) sample contains a single cell type or do not merge samples with different cell types')
    }

    if(cell_type == 'B cell'){
      cell_type <- 'BCR'
    }else if(cell_type == 'T cell'){
      cell_type <- 'TCR'
    }else{
      stop(paste0('Unrecognized cell type ', cell_type))
    }

    console_command <- paste0(os_prefix, 'enclone ', cell_type, '=', '"', path_with_origin, '"', ' POUT=', temp_file, ' PCELL PCOLS=barcode,group_id NOPRINT')
    system(console_command)

    return(list('sample_df'=sample_df, 'sample_names'=sample_names, 'csv_file'=temp_file))
  }


  #SUBROUTINE 3: merge the new clonotypes into the sample dfs, using the stripped_barcodes column. Recalculate clonotype frequencies, reorder sample df, clean-up

  merge_enclone_clonotypes <- function(enclone_out){
    new_clonotypes <- utils::read.csv(enclone_out$csv_file)
    new_clonotypes <- new_clonotypes[c('group_id', 'barcode')]
    names(new_clonotypes)[names(new_clonotypes) == 'group_id'] <- 'new_group_id'
    names(new_clonotypes)[names(new_clonotypes) == 'barcode'] <- 'stripped_barcodes'
    new_clonotypes$stripped_barcodes <- lapply(new_clonotypes$stripped_barcodes, function(x) unlist(stringr::str_split(x, '-'))[1])

    sample_df <- enclone_out$sample_df
    sample_df <- strip_barcodes(sample_df)

    sample_out <- merge(sample_df, new_clonotypes, by='stripped_barcodes', all.x = TRUE)
    sample_out$new_group_id <- unlist(lapply(sample_out$new_group_id, function(x) paste0('clonotype', x)))
    sample_out$clonotype_id <- sample_out$new_group_id
    sample_out$clonotype_id_10x <- sample_out$new_group_id

    sample_out$clonotype_frequency <- unlist(lapply(sample_out$clonotype_id, function(x) length(which(sample_out$clonotype_id == x))))

    sample_out$stripped_barcodes <- NULL
    sample_out$new_group_id <- NULL


    return(sample_out)
  }

  #Get the OS prefix for Windows/Darwin and Linux
  if(operating.system == 'Windows'){
    os_prefix <- 'cmd.exe '
  }else{
    os_prefix <- ''
  }

  #Enclone treats differently samples from the same donor but with the same origin/different origins - visit https://10xgenomics.github.io/enclone/ for more information
  if(same.origin){
    merge_origin_sep <- ','
  }else{
    merge_origin_sep <- ':'
  }

  #Lowercase sample names and file names so they match
  VDJ$sample_id <- unlist(lapply(VDJ$sample_id, function(x) tolower(x)))

  #To keep the same VDJ order in the output
  if(!('X' %in% colnames(VDJ))) {
    VDJ$X <- 1:nrow(VDJ)
  }

  samples.to.clonotype <- unlist(lapply(samples.to.clonotype, function(x) tolower(x)))
  #Create the temp directory to output the parseable enclone files (per cell)
  temp_dir <- './temp_call_enclone'
  if(!dir.exists(temp_dir)) dir.create(temp_dir)

  #Check if enclone - the command line version - is installed on the user's system
  enclone_check <- system(paste0(os_prefix, "enclone"), intern = TRUE)
  if(any(stringr::str_detect(enclone_check, 'command not found'))){
    stop('Please ensure enclone is installed on your system - visit https://10xgenomics.github.io/enclone/ for more information!')
  }

  #Check if the VDJ directory sample names correspond to the VDJ sample names
  repertoire.number <- unique(VDJ$sample_id)
  if(any(samples.to.clonotype=='all')){
    samples_to_clonotype <- repertoire.number
  }else{
    samples_to_clonotype <- samples.to.clonotype
  }

  sample_files <- list.files(VDJ.directory)
  sample_files <- unlist(lapply(sample_files, function(x) tolower(x)))
  #for(sample in samples_to_clonotype){
  #  if(!(tolower(sample) %in% sample_files)){
  #    stop(paste0('Unable to find the VDJ out file for sample ', sample))
  #  }
  #}

  samples_not_clonotyped <- VDJ[!(VDJ$sample_id %in% samples_to_clonotype),]
  VDJ <- VDJ[VDJ$sample_id %in% samples_to_clonotype,]

  #Segment the main VDJ into sample dfs - these will be used to indicate the file directory names for enclone (via sample_dfs[[1]]$sample_id)
  sample_dfs <- list()
  if(!global.clonotype){
    for(i in 1:length(repertoire.number)){
      sample_dfs[[i]] <- VDJ[which(VDJ$sample_id==repertoire.number[i]),]
    }

  }else if(!is.null(samples.to.combine) & global.clonotype){
    remaining_dfs <- 0
    samples.to.combine <- unlist(lapply(samples.to.combine, function(x) tolower(x)))
    if(is.vector(samples.to.combine) & !is.list(samples.to.combine)){
      for(i in 1:length(samples.to.combine)){
        if(!(tolower(samples.to.combine[i]) %in% VDJ$sample_id)){
          stop(paste0('Sample ', samples.to.combine[i], ' not found in your VDJ object. Ensure the spelling is correct!'))
        }
      }
      sample_dfs[[1]] <- VDJ[VDJ$sample_id %in% samples.to.combine,]
      remaining_dfs <- remaining_dfs + 1

    }else if(is.list(samples.to.combine) & is.vector(samples.to.combine)){
      for(i in 1:length(samples.to.combine)){
        samples.to.combine[[i]] <- lapply(samples.to.combine[[i]], function(x) tolower(x))

        for(j in 1:length(samples.to.combine[[i]])){
          if(!(tolower(samples.to.combine[[i]][[j]]) %in% VDJ$sample_id)){
            stop(paste0('Sample ', samples.to.combine[[i]][[j]], ' not found in your VDJ object. Ensure the spelling is correct!'))
          }
        }
        sample_dfs[[i]] <- VDJ[VDJ$sample_id %in% samples.to.combine[[i]],]
        remaining_dfs <- remaining_dfs + 1
      }
    }

    remaining_samples <- repertoire.number[which(!(repertoire.number) %in% unlist(samples.to.combine))]

    if(length(remaining_samples)!=0){
      for(i in 1:length(remaining_samples)){
        if(nrow(VDJ[which(VDJ$sample_id==remaining_samples[i]),])!=0){
          sample_dfs[[i+remaining_dfs]] <- VDJ[which(VDJ$sample_id==remaining_samples[i]),]
        }
      }
    }
  }else{
    sample_dfs[[1]] <- VDJ
  }


  #Lapply the subroutines on the sample_dfs list of sample VDJ dataframes
  if(parallel){
    #requireNamespace('parallel')
    cores <- parallel::detectCores() - 1
    enclone_outs <- parallel::mclapply(sample_dfs, get_enclone_output, mc.cores=cores)
    sample_outs <- parallel::mclapply(sample_dfs, merge_enclone_clonotypes, mc.cores=cores)

  }else{
    enclone_outs <- lapply(sample_dfs, function(x) get_enclone_output(x))
    sample_outs <- lapply(enclone_outs, function(x) merge_enclone_clonotypes(x))
  }

  unlink(temp_dir, recursive = TRUE)

  if(output.format == 'vgm'){
    vdj <- do.call('rbind', sample_outs)
    vdj <- rbind(vdj, samples_not_clonotyped)
    vdj <- vdj[order(vdj$X),]
    vdj$X <- NULL
  }else{
    for(i in 1:length(sample_outs)){
      sample_outs[[i]] <- sample_outs[[i]][order(unlist(sample_outs[[i]]$X)),]
    }
    vdj <- sample_outs
  }

  return(vdj)
}
