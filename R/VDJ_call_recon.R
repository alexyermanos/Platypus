#' Calls the Kaplinsky/RECON tool

#'@description Calls the Kaplinsky/RECON tool on the VDJ/VDJ.GEX.matrix[[1]] object to infer the parent distribution of species and estimate their diversity. Outputs either a dataframe of the resulting means and weights of the RECON species parent distribution estimation or a plot of the original species distribution along resampled values from the reconstructed parent distribution.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param recon.directory directory containing recon executable. Defaults to {working directory}/Recon
#' @param feature.columns vector of strings/ string - columns denoting the unique species for RECON - e.g., could be the CDRH3s if feature.columns = 'VDJ_cdr3s_aa'. If more than one column is provided (e.g. c("VDJ_cdr3s_aa","VJ_cdr3s_aa")) these columns will be pasted together.
#' @param grouping.column string - column determining the groups/ samples for the species observations. Defaults to 'sample_id' for per-repertoire analysis
#' @param VDJ.VJ.1chain boolean, defaults to TRUE. Whether to filter out aberrant cells (more than 1 VDJ or VJ chain).
#' @param max.features integer or 'all' - maximum number of features/species ot be considered for the RECON estimation. If 'all', will consider all species.
#' @param size.threshold integer - the size threshold parameter for the RECON tool, as specified by the '-t' parameter.
#' @param resample boolean - if T, will also perform and output a resample of the species frequencies/sizes from the inferred parent distribution.
#' @param max.feature.size integer - the maximum size of species/features to be considered in the resulting plot (maximum number of elements on the x axis).
#' @param reticulate boolean - if T, will create a new environment to install python and run the RECON tool, else, your environment must have a python version compatible with RECON installed.
#' @param operating.system string - operating system on which RECON will be run. 'Windows' for Windows, 'Linux' for Linux, 'Darwin' for MacOS.

#' @return The resulting means and weights of the RECON-inferred distribution as a seprate dataframe or appended to the VDJ, or a plot of resampled species sizes from the inferred distribution vs original sizes/frequencies.
#' @export


#' @examples
#' \dontrun{
#' VDJ_call_RECON(VDJ, recon.directory='./Recon',
#' feature.columns = 'VDJ_cdr3s_aa', grouping.column = 'VDJ_cdr3s_aa')
#'}


VDJ_call_RECON <- function(VDJ,
                           recon.directory,
                           feature.columns,
                           grouping.column,
                           VDJ.VJ.1chain,
                           max.features,
                           size.threshold,
                           resample,
                           max.feature.size,
                           reticulate,
                           operating.system
                           ){

 if(missing(VDJ)) stop('Please input your data as a vgm[[1]]/VDJ dataframe')
 if(missing(recon.directory)) recon.directory <- paste0(getwd(),'/Recon')
 if(missing(feature.columns)) feature.columns <- 'VDJ_cdr3s_aa'
 if(missing(grouping.column)) grouping.column <- 'sample_id'
 if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
 if(missing(max.features)) max.features <- 10
 if(missing(size.threshold)) size.threshold <- 30
 if(missing(resample)) resample <- F
 if(missing(max.feature.size)) max.feature.size <- 50
 if(missing(reticulate)) reticulate <- F


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

  size.threshold <- paste0('-t ', size.threshold)

  if(reticulate==T){
    version = '2.7.15'
    path <- reticulate::install_python(version = version)
    Sys.setenv(RETICULATE_PY = path)
    #Sys.getenv('RETICULATE_PY')
    reticulate::use_python_version(version=version)
    #reticulate::conda_create(envname = "recon_temp", packages = c('python=2.7', 'numpy', 'scipy'))
    #reticulate::use_condaenv('recon_temp')
    reticulate::virtualenv_install('recon', c('scipy', 'numpy'))
    reticulate::use_virtualenv('recon')
    reticulate::import('numpy')
    reticulate::import('scipy')
  }


  #Preprocess VGM for Recon compatibility (output = txt file first col = species/clone name; col2 = clone frequency)
  get_abundances <- function(VDJ, feature.columns, grouping.column, VDJ.VJ.1chain){

    if(length(feature.columns) > 1){
      combine.features <- T
    }else{
      combine.features <- F
    }

    abundance_df <- VDJ_abundances(VDJ,
                                   feature.columns = feature.columns,
                                   proportions = 'absolute',
                                   grouping.column = grouping.column,
                                   max.groups = NULL,
                                   specific.groups = 'none',
                                   sample.column = 'none',
                                   VDJ.VJ.1chain = VDJ.VJ.1chain,
                                   treat.incomplete.groups = 'exclude',
                                   treat.incomplete.features = 'exclude',
                                   combine.features = combine.features,
                                   treat.combined.features = 'exclude',
                                   treat.combined.groups = 'exclude',
                                   specific.feature.colors = NULL,
                                   output.format = 'abundance.df')

    return(abundance_df)
  }

  out_directory <- paste0(getwd(),'/recon_out')
  temp_directory <- paste0(getwd(), '/recon_temp')
  recon_directory <- paste0(recon.directory, '/recon.py')


  abundance_df <- get_abundances(VDJ,
                                 feature.columns = feature.columns,
                                 grouping.column = grouping.column,
                                 VDJ.VJ.1chain = VDJ.VJ.1chain)

  unique_groups <- unique(abundance_df$group)
  output_plots <- list()

  for(i in 1:length(unique_groups)){
    if(!dir.exists(out_directory)) dir.create(out_directory)
    if(!dir.exists(temp_directory)) dir.create(temp_directory)

    abundance_df_subset <- abundance_df[abundance_df$group == unique_groups[i],]

    temp_df <- data.frame(Feature=abundance_df$unique_feature_values,
                          Frequencies=abundance_df$feature_value_counts)

    temp_df <- temp_df %>% dplyr::arrange(desc(Frequencies))
    id <- unique_groups[i]

    table_file_name <- paste0(temp_directory, '/' , id, '_', 'temp.txt')
    utils::write.table(temp_df, file=table_file_name, sep="\t", row.names=F, col.names=F)
    temp_df <- NULL


    fitfile_file_name <- paste0(temp_directory, '/' , id, '_', 'fitfile.txt')
    output_file_name <- paste0(out_directory, '/' , id, '_', 'output.txt')
    system(paste0(os_prefix, 'python ', recon_directory, ' -R ', size.threshold, ' -o ', fitfile_file_name, ' ', table_file_name, ' > ', output_file_name))
    output_string <-scan(output_file_name, what='', sep='\n')

    #######Will not include in the vgm - only a test to check if Recon works correctly
    #observed_species_numbers_keep <- stringr::str_extract(output_string[2], pattern='\\{(.*?)\\}')
    #observed_species_numbers <- stringr::str_replace(observed_species_numbers_keep, '\\{', '')
    #observed_species_numbers <- stringr::str_replace(observed_species_numbers, '\\}', '')
    #observed_species_numbers <- stringr::str_split(observed_species_numbers, ',', simplify=T)
    #observed_species_numbers <- stringr::str_split_fixed(t(observed_species_numbers), ':', n=2)
    #observed_species_numbers <- data.frame(Counts = as.numeric(observed_species_numbers[,1]), Recon_predicted_species_number = as.numeric(observed_species_numbers[,2]))
    #actual_species_numbers <- lapply(observed_species_numbers$Counts, function(x) length(unique(sample_dfs[[i]]$clonotype_id[which(sample_dfs[[i]]$recon_clonotype_frequency==x)])))
    #observed_species_numbers$Actual_unique_species_numbers <- actual_species_numbers
    ########

    #formatted_string <- stringr::str_split(output_string[2], ']', simplify=T)
    #weights <- stringr::str_split(stringr::str_replace(formatted_string[1], '\\(\\[', ''), ',')
    #sample_dfs[[i]]$recon_weights <- rep(weights, nrow(sample_dfs[[i]]))

    #means <- stringr::str_split(stringr::str_replace(formatted_string[2], ', \\[', ''), ',')
    #sample_dfs[[i]]$recon_means <- rep(means, nrow(sample_dfs[[i]]))

    #missing_species <- unlist(stringr::str_split(stringr::str_replace(formatted_string[3], stringr::fixed(observed_species_numbers_keep), ''), ','))[2]
    #sample_dfs[[i]]$recon_missing_species <- rep(as.numeric(missing_species), nrow(sample_dfs[[i]]))

    #total_species <- as.numeric(length(unique(sample_dfs[[i]]$clonotype_id))) + as.numeric(missing_species)
    #sample_dfs[[i]]$recon_total_species <- rep(as.numeric(total_species), nrow(sample_dfs[[i]]))

    #log_likelihood <- unlist(stringr::str_split(stringr::str_replace(formatted_string[3], stringr::fixed(observed_species_numbers_keep), ''), ','))[5]
    #sample_dfs[[i]]$recon_log_likelihood <- rep(as.numeric(log_likelihood), nrow(sample_dfs[[i]]))

    plotfile_name <- paste0(out_directory, '/' , id, '_plotfile.txt')
    pdf_name <- paste0(out_directory, '/', id, '_plot.pdf')

    system(paste0(os_prefix, 'python ', recon_directory, ' -x ', ' -o ', plotfile_name, ' -b ', paste0(recon.directory, '/error_bar_parameters.txt'), ' ', fitfile_file_name, ' > ', plotfile_name))

    df <- utils::read.table(plotfile_name, sep='\t', header=T, fill=T)
    if(!(is.null(max.feature.size))) df <- df[1:max.feature.size-1,]

    df <- tidyr::pivot_longer(df, cols=c('sample', 'fit'), names_to='Origin', values_to='Total number of features')
    df <-  df[which(df$clone_size!=0),]
    df$`Origin`[which(df$`Origin`=='sample')] <- 'Observed in sample/group'
    df$`Origin`[which(df$`Origin`=='fit')] <- 'Sampled from inferred distribution'

    names(df)[names(df)=='clone_size'] <- 'Cells per feature'

    #Global variable definition for CRAN checks
    `Clone size` <- NULL
    `Total number of clones` <- NULL
    Origin <- NULL


    output_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x=`Cells per feature`, y=`Total number of features`, shape=`Origin`, color=`Origin`)) + ggplot2::scale_shape_manual(values=c('Observed in sample'=1, 'Sampled from inferred distribution'=4))+
        ggplot2::scale_color_manual(values=c('Observed in sample'='black', 'Sampled from inferred distribution'='red')) + ggplot2::theme_bw() +
        ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), axis.ticks = ggplot2::element_line(colour = "black", size = 3),
        axis.line = ggplot2::element_line(colour = 'black', size = 1), axis.text = ggplot2::element_text(face="bold", size=14), text=ggplot2::element_text(size=14)) +
        ggplot2::scale_x_continuous(limits=c(1,max.feature.size), n.breaks = max.feature.size-2) + ggplot2::scale_y_continuous(limits=c(0,max(df$`Total number of features`))) + ggplot2::ggtitle(id) + ggplot2::geom_point(size=5, stroke=1.8)

    if(save.pdf){
      grDevices::pdf(pdf_name)
      plot(output_plots[[i]])
      grDevices::dev.off()
    }

    if(resample){
      resampled_file <- paste0(out_directory, '/', id, '_resample_file.txt')
      system(paste0(os_prefix, 'python ', recon_directory, ' -r ', ' -o ', resampled_file, ' ',  fitfile_file_name, ' > test_resample.txt'))
    }

      unlink(temp_directory, recursive=T)
  }


  if(reticulate){
    reticulate::virtualenv_remove('recon')
    base::Sys.unsetenv(path)
  }

  return(output_plots)
}
