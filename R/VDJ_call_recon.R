#' Calls the Kaplinsky/RECON tool

#'@description Calls the Kaplinsky/RECON tool on the VDJ/VDJ.GEX.matrix[[1]] object to infer the parent distribution of clonotypes and estimate their diversity. Outputs either a dataframe of the resulting means and weights of the RECON clonotype parent distribution estimation or a plot of the original clonotype distribution along resampled values from the reconstructed parent distribution.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param recon.directory directory containing recon executable. Defaults to {working directory}/Recon
#' @param sample.id boolean - if F, clonotypes will be considered at a global level, irrespective of samples.
#' @param clone.list list - if empty, RECON will be used to estimate the diversity of all clonotypes, else it will only consider the specified clonotypes.
#' @param max.clones integer or 'all' - maximum number of clones ot be considered for the RECON estimation. If 'all', will consider all clonotypes.
#' @param size.threshold integer - the size threshold parameter for the RECON tool, as specified by the '-t' parameter.
#' @param resample boolean - if T, will also perform and output a resample of the clonotype frequencies/sizes from the inferred parent distribution.
#' @param plot.results boolean - if T, will save a pdf of the resampled clonotype frequencies/sizes from the RECON-inferred distribution along the original frequencies.
#' @param max.clone.size integer - the maximum size of clones to be considered in the resulting plot (maximum number of elements on the x axis).
#' @param reticulate boolean - if T, will create a new environment to install python and run the RECON tool, else, your environment must have a python version compatible with RECON installed.
#' @param output.format string - 'vgm' will append the means and weights of the RECON-inferred distribution to the VDJ/VDJ.GEX.matrix[[1]] object, 'recon' will output a new dataframe of these weights, 'plots' will output the ggplot2 objects (if plot.results=T).
#' @param operating.system string - operating system on which RECON will be run. 'Windows' for Windows, 'Linux' for Linux, 'Darwin' for MacOS.

#' @return The resulting means and weights of the RECON-inferred distribution as a seprate dataframe or appended to the VDJ, or a plot of resampled clonotype sizes from the inferred distribution vs original sizes/frequencies.
#' @export


#' @examples
#' \dontrun{
#' VDJ_call_recon(VDJ, recon.directory='./Recon',
#' max.clones='all', sample.id=T, resample=F,
#' plot.results=T, output.format='vgm',
#' reticulate=T, operating.system='Darwin')
#'}

VDJ_call_recon <- function(VDJ,
                            recon.directory,
                            sample.id,
                            clone.list,
                            max.clones,
                            size.threshold,
                            resample,
                            plot.results,
                            max.clone.size,
                            reticulate,
                            output.format,
                            operating.system
                            ){

 if(missing(VDJ)) stop('Please input your data as a vgm[[1]]/VDJ dataframe')
 if(missing(recon.directory)) recon.directory <- paste0(getwd(),'/Recon')
 if(missing(sample.id)) sample.id <- T
 if(missing(clone.list)) clone.list <- list()
 if(missing(max.clones)) max.clones <- 'all'
 if(missing(size.threshold)) size.threshold <- 30
 if(missing(resample)) resample <- F
 if(missing(plot.results)) plot.results <- F
 if(missing(max.clone.size)) max.clone.size <- 50
 if(missing(output.format)) output.format <- 'vgm'
 if(missing(reticulate)) reticulate <- T



 if(missing(operating.system)){
     switch(Sys.info()[['sysname']],
            Windows= {message("Windows system detected")
                      operating.system <- "Windows"},
            Linux  = {message("Linux system detected")
                     operating.system <- "Linux"},
            Darwin = {message("MAC system detected")
                     operating.system <- "Darwin"})
  }

  if(operating.system=='Linux' | operating.system=='Darwin'){
    os_prefix <- ''
  }else{
    os_prefix <- 'cmd.exe '
  }

  size_threshold <- paste0('-t ', size.threshold)

  if(reticulate==T){
    #require(reticulate)

    version <- '2.7.15'

    #reticulate::install_python(version = version)

    path <- reticulate::install_python(version = version)
    Sys.setenv(RETICULATE_PY = path)
    #Sys.getenv('RETICULATE_PY')
    ###reticulate::use_python_version(version=version)
    reticulate::virtualenv_create(envname = "recon", python_version=version)
    reticulate::virtualenv_install('recon', c('scipy', 'numpy'))
    reticulate::use_virtualenv('recon', required = T)
    reticulate::import('numpy')

  }


  #Preprocess VGM for Recon compatibility (output = txt file first col = species/clone name; col2 = clone frequency)
  VDJ.matrix <- VDJ
  VDJ <- NULL
  sample_dfs <- list()
  output_plots <- list()
  if(sample.id==T){
    repertoire.number <- unique(VDJ.matrix$sample_id)
    for(i in 1:length(repertoire.number)){
      sample_dfs[[i]] <- VDJ.matrix[which(VDJ.matrix$sample_id==repertoire.number[i]),]
    }
  }


  for(i in 1:length(sample_dfs)){
    recon_directory <- paste0(recon.directory, '/recon_v3.0.py')
    if(sample.id==F){
      id <- 'global'
    }else{
      id <- repertoire.number[i]
    }

    out_directory <- paste0(getwd(),'/recon_out')
    temp_directory <- paste0(getwd(), '/recon_temp')

    if(!dir.exists(out_directory)) dir.create(out_directory)
    if(!dir.exists(temp_directory)) dir.create(temp_directory)


    if(length(clone.list)==0){
      #Recalculating clonotype frequency to be compatible with global clonotypes
      all_clones <- sample_dfs[[i]]$clonotype_id
      clone_frequencies <- lapply(all_clones, function(x) length(which(sample_dfs[[i]]$clonotype_id==x)))
      sample_dfs[[i]]$recon_clonotype_frequency <- clone_frequencies

      #Order df based on clonotype frequencies
      sample_dfs[[i]] <- sample_dfs[[i]][order(unlist(sample_dfs[[i]]$recon_clonotype_frequency), decreasing=T),]

      if(max.clones=='all'){
        unique_clones <- unique(sample_dfs[[i]]$clonotype_id)
      }else{
        unique_clones <- unique(sample_dfs[[i]]$clonotype_id)[1:max.clones]
      }

      unique_clone_frequencies <- lapply(unique_clones, function(x) unique(sample_dfs[[i]]$recon_clonotype_frequency[which(sample_dfs[[i]]$clonotype_id==x)]))


    }else{
      all_clones <- c(sample_dfs[[i]]$clonotype_id)
      clone_frequencies <- lapply(all_clones, function(x) if(x %in% clone.list) length(which(sample_dfs[[i]]$clonotype_id==x)) else '')
      sample_dfs[[i]]$recon_clonotype_frequency <- clone_frequencies
      sample_dfs[[i]] <- sample_dfs[[i]][order(unlist(sample_dfs[[i]]$recon_clonotype_frequency), decreasing=T),]

      unique_clones <- unique(sample_dfs[[i]]$clonotype_id[which(sample_dfs[[i]]$recon_clonotype_frequency!='')])
      unique_clone_frequencies <- lapply(unique_clones, function(x) unique(sample_dfs[[i]]$recon_clonotype_frequency[which(sample_dfs[[i]]$clonotype_id==x)]))
    }

    temp_df <- data.frame(Clone=unlist(unique_clones),
                          Frequencies=unlist(unique_clone_frequencies))
    table_file_name <- paste0(temp_directory, '/' , id, '_', 'temp.txt')
    utils::write.table(temp_df, file=table_file_name, sep="\t", row.names=F, col.names=F)
    temp_df <- NULL


    fitfile_file_name <- paste0(temp_directory, '/' , id, '_', 'fitfile.txt')
    output_file_name <- paste0(out_directory, '/' , id, '_', 'output.txt')
    system(paste0(os_prefix, 'python ', recon_directory, ' -R ', size_threshold, ' -o ', fitfile_file_name, ' ', table_file_name, ' > ', output_file_name))

    #######Will not include in the vgm - only a test to check if Recon works correctly
    output_string <-scan(output_file_name, what='', sep='\n')
    observed_species_numbers_keep <- stringr::str_extract(output_string[2], pattern='\\{(.*?)\\}')
    observed_species_numbers <- stringr::str_replace(observed_species_numbers_keep, '\\{', '')
    observed_species_numbers <- stringr::str_replace(observed_species_numbers, '\\}', '')
    observed_species_numbers <- stringr::str_split(observed_species_numbers, ',', simplify=T)
    observed_species_numbers <- stringr::str_split_fixed(t(observed_species_numbers), ':', n=2)
    observed_species_numbers <- data.frame(Counts = as.numeric(observed_species_numbers[,1]), Recon_predicted_species_number = as.numeric(observed_species_numbers[,2]))
    actual_species_numbers <- lapply(observed_species_numbers$Counts, function(x) length(unique(sample_dfs[[i]]$clonotype_id[which(sample_dfs[[i]]$recon_clonotype_frequency==x)])))
    observed_species_numbers$Actual_unique_species_numbers <- actual_species_numbers
    ########

    formatted_string <- stringr::str_split(output_string[2], ']', simplify=T)
    weights <- stringr::str_split(stringr::str_replace(formatted_string[1], '\\(\\[', ''), ',')
    sample_dfs[[i]]$recon_weights <- rep(weights, nrow(sample_dfs[[i]]))

    means <- stringr::str_split(stringr::str_replace(formatted_string[2], ', \\[', ''), ',')
    sample_dfs[[i]]$recon_means <- rep(means, nrow(sample_dfs[[i]]))

    missing_species <- unlist(stringr::str_split(stringr::str_replace(formatted_string[3], stringr::fixed(observed_species_numbers_keep), ''), ','))[2]
    sample_dfs[[i]]$recon_missing_species <- rep(as.numeric(missing_species), nrow(sample_dfs[[i]]))

    total_species <- as.numeric(length(unique(sample_dfs[[i]]$clonotype_id))) + as.numeric(missing_species)
    sample_dfs[[i]]$recon_total_species <- rep(as.numeric(total_species), nrow(sample_dfs[[i]]))

    log_likelihood <- unlist(stringr::str_split(stringr::str_replace(formatted_string[3], stringr::fixed(observed_species_numbers_keep), ''), ','))[5]
    sample_dfs[[i]]$recon_log_likelihood <- rep(as.numeric(log_likelihood), nrow(sample_dfs[[i]]))

    if(plot.results==T){
      plotfile_name <- paste0(out_directory, '/' , id, '_plotfile.txt')
      pdf_name <- paste0(out_directory, '/', id, '_plot.pdf')
      recon_directory <- paste0(recon.directory, '/recon_v2.5.py')
      system(paste0(os_prefix, 'python ', recon_directory, ' -x ', ' -o ', plotfile_name, ' -b ', paste0(recon.directory, '/error_bar_parameters.txt'), ' ', fitfile_file_name, ' > test_plot.txt'))

      df <- utils::read.table(plotfile_name, sep='\t', header=T, fill=T)
      if(!(is.null(max.clone.size))) df <- df[1:max.clone.size-1,]

      df <- tidyr::pivot_longer(df, cols=c('sample', 'fit'), names_to='Origin', values_to='Total number of clones')
      df <-  df[which(df$clone_size!=0),]
      df$`Origin`[which(df$`Origin`=='sample')] <- 'Observed in sample'
      df$`Origin`[which(df$`Origin`=='fit')] <- 'Sampled from inferred distribution'

      names(df)[names(df)=='clone_size'] <- 'Clone size'

      #Global variable definition for CRAN checks
      `Clone size` <- NULL
      `Total number of clones` <- NULL
      Origin <- NULL


       output_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x=`Clone size`, y=`Total number of clones`, shape=`Origin`, color=`Origin`)) + ggplot2::scale_shape_manual(values=c('Observed in sample'=1, 'Sampled from inferred distribution'=4))+
            ggplot2::scale_color_manual(values=c('Observed in sample'='black', 'Sampled from inferred distribution'='red')) + ggplot2::theme_bw() +
            ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), axis.ticks = ggplot2::element_line(colour = "black", size = 3),
            axis.line = ggplot2::element_line(colour = 'black', size = 1), axis.text = ggplot2::element_text(face="bold", size=14), text=ggplot2::element_text(size=14)) +
            ggplot2::scale_x_continuous(limits=c(1,max.clone.size), n.breaks = max.clone.size-2) + ggplot2::scale_y_continuous(limits=c(0,max(df$`Total number of clones`))) + ggplot2::ggtitle(id) + ggplot2::geom_point(size=5, stroke=1.8)


      grDevices::pdf(pdf_name)
      plot(output_plots[[i]])
      grDevices::dev.off()
    }

    if(resample==T){
      resampled_file <- paste0(out_directory, '/', id, '_resample_file.txt')
      system(paste0(os_prefix, 'python ', recon_directory, ' -r ', ' -o ', resampled_file, ' ',  fitfile_file_name, ' > test_resample.txt'))
    }

      unlink(temp_directory, recursive=T)
 }


  if(reticulate==T){
    reticulate::virtualenv_remove('recon')
    base::Sys.unsetenv(path)
  }

  if(output.format=='vgm'){
    VDJ.matrix <- do.call('rbind', sample_dfs)

    return(VDJ.matrix)

  }else if(output.format=='recon'){
    VDJ.matrix <- do.call('rbind', sample_dfs)

    recon_weights <- c()
    recon_means <- c()
    recon_missing_species <- c()
    recon_total_species <- c()
    recon_log_likelihood <- c()

    for(i in 1:length(sample_dfs)){
      recon_weights <- c(recon_weights, sample_dfs[[i]]$recon_weights[1])
      recon_means <- c(recon_means, sample_dfs[[i]]$recon_means[1])
      recon_missing_species <- c(recon_missing_species, unique(sample_dfs[[i]]$recon_missing_species))
      recon_total_species <- c(recon_total_species, unique(sample_dfs[[i]]$recon_total_species))
      recon_log_likelihood<- c(recon_log_likelihood, unique(sample_dfs[[i]]$recon_log_likelihood))

      out.df <- data.frame(sample_id=repertoire.number,
                           recon_weights=matrix(recon_weights),
                           recon_means=matrix(recon_means),
                           recon_missing_species=recon_missing_species,
                           recon_total_species=recon_total_species,
                           recon_log_likelihood=recon_log_likelihood)

      return(out.df)
    }

  }else if(output.format=='plots'){
    return(output_plots)
  }
}
