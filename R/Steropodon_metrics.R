#' Plot the structure-level physicochemical properties obtained form Steropodon_properties


#' @description Function to visualize the physicochemical properties calculated using Steropoodon_properties.
#'
#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param plot.format string - either 'bar' for bar plots of the values denoted in the 'feature' parameter, grouped by the 'grouping parameter, or line plots of property per residue id if plot.format = 'line'
#' @param plot.level string - the level at which properties should be averages for the bar plots. 'global' for all structures, 'sample' to obtain bar plots for each sample, 'clonotype' for bar plots across all samples and clonotypes.
#' @param feature string - the physicochemical property to be quantified/plotted. Options include 'charge', 'hydrophobicity', 'SASA', 'b-factor'.
#' @param grouping string - grouping factor for the bar plots ('region' for framework/hypervariable regions, 'chain' to group across heavy or light chains).
#' @param proportions bool - if TRUE, will create bar plots of proportions (e.g., proportions of residues per region) instead of absolute values.
#' @param combine.lineplot bool - if TRUE, will combine all line plots (defined by plot.level) into a single plot.


#' @return bar plots if plot.format = 'bar' or line plots if plot.format = 'line' of the properties specified in feature.
#' @export
#' @examples
#' \dontrun{
#'steropodon_properties <-
#'  steropodon_igfold %>%
#'  Steropodon_properties(structure = 'structure',
#'  properties = c('SASA', 'charge', 'hydrophobicity', 'DSSP'),
#'  dssp.exefile = '/opt/homebrew/Caskroom/miniforge/base/bin/mkdssp',
#'  parallel = F)
#'
#'steropodon_properties %>%
#'  Steropodon_metrics(plot.format = 'line',
#'                     feature = 'charge',
#'                     plot.level = 'global')
#'}


Steropodon_metrics <- function(steropodon.object,
                               structure,
                               plot.format,
                               plot.level,
                               feature,
                               grouping,
                               proportions,
                               combine.lineplot
                             ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(structure)) structure <- 'structure'
  if(missing(plot.format)) plot.format <- 'bar'
  if(missing(plot.level)) plot.level <- 'global'
  if(missing(feature)) feature <- 'region'
  if(missing(grouping)) grouping <- 'chain'
  if(missing(proportions)) proportions <- T
  if(missing(combine.lineplot)) combine.lineplot <- F


  create_lineplot <- function(df.list,
                              feature,
                              combine.lineplot
                             ){

    if(combine.lineplot){
      temp <- do.call('rbind', df.list)
      df.list <- list()
      df.list[[1]] <- temp
    }

    out_plots <- vector(mode = 'list', length = length(df.list))

    for(i in 1:length(df.list)){
      df <- df.list[[i]]

      #Will always subset for unique residues
      df <- df[df$elety == 'CA',]
      df <- df[!is.na(df[[feature]]),]

      if(combine.lineplot){
        df$feature <- df[[feature]]

        df <- df %>%
              dplyr::group_by(lineplot_group) %>%
              dplyr::group_by(as.character(resno)) %>%
              dplyr::mutate(mean = mean(feature)) %>%
              dplyr::ungroup()

        df <- df %>%
              dplyr::group_by(lineplot_group) %>%
              dplyr::distinct(resno, .keep_all = T) %>%
              dplyr::ungroup()

        df$feature <- df$mean

        out_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = resno, y = feature, group = lineplot_group)) +
                          ggplot2::geom_line(ggplot2::aes(color = lineplot_group))+
                          ggplot2::labs(title = 'Combined', x = 'Residue number', y = feature, color = 'group') +
                          ggplot2::theme(panel.background = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

      }else{
        df$feature <- df[[feature]]
        df <- df %>%
              dplyr::group_by(as.character(resno)) %>%
              dplyr::mutate(mean = mean(feature)) %>%
              dplyr::ungroup()

        df <- df %>%
              dplyr::distinct(resno, .keep_all = T)

        df$feature <- df$mean


        out_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = resno, y = feature, group = 1)) +
                          ggplot2::geom_line()+
                          ggplot2::labs(title = unique(df$lineplot_group), x = 'Residue number', y = feature) +
                          ggplot2::theme(panel.background = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none")

      }
    }

    return(out_plots)
  }


  create_barplot <- function(df.list,
                             feature,
                             grouping,
                             proportions
                            ){

    out_plots <- vector(mode = 'list', length = length(df.list))

    for(i in 1:length(df.list)){
      df <- df.list[[i]]

      if(length(feature) > 1){
        new_feature <- paste0(feature, collapse = '_')
        df$combined <- rep('', nrow(df))

        for(col in feature){
          df$combined <- paste0(df$combined, '_', pdb$atom[[col]])
        }

        df <- df %>%
          dplyr::mutate(combined = substring(combined, 2))

        colnames(df)[colnames(df) == 'combined'] <- new_feature
        feature <- new_feature
      }

      if(length(grouping) > 1){
        new_feature <- paste0(grouping, collapse = '_')
        df$combined <- rep('', nrow(df))

        for(col in color.by){
          df$combined <- paste0(df$combined, '_', df[[col]])
        }

        df <- df %>%
          dplyr::mutate(combined = substring(combined,2))

        colnames(df)[colnames(df) == 'combined'] <- new_feature
        grouping <- new_feature
      }

      #Will always subset for unique residues
      df <- df[df$elety == 'CA',]

      #Remove NAs
      df <- df[!is.na(df[[grouping]]),]
      df <- df[!is.na(df[[feature]]),]

      if(grouping == 'none'){
        df$none <- 'none'
      }

      if(is.numeric(df[[feature]])){
        df$grouping <- df[[grouping]]
        df$feature <- df[[feature]]
        df <- df %>%
              dplyr::group_by(grouping) %>%
              dplyr::mutate(mean = mean(feature))
        df$feature <- df$mean

        df <- df %>%
              dplyr::distinct(grouping, .keep_all = T)

        df$colors <- grDevices::rainbow(nrow(df))

        title <- paste0('Average values of ', feature, ' per ', grouping, ' for ', names(df.list)[i])
        out_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = grouping, y = feature, fill = colors)) +
                          ggplot2::geom_bar(show.legend = F, stat = "identity", width=0.6, color="black") +
                          ggplot2::labs(title = title, x = grouping, y = feature) +
                          ggplot2::theme(panel.background = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none") +
                          ggplot2::coord_flip()

      }else{
        df$grouping <- df[[grouping]]
        df$feature <- df[[feature]]

        df <- df %>%
              dplyr::group_by(grouping) %>%
              dplyr::group_by(feature) %>%
              dplyr::mutate(counts = n()) %>%
              dplyr::ungroup()

        if(proportions){
          df <- df %>%
                dplyr::group_by(grouping) %>%
                dplyr::mutate(group_counts = n()) %>%
                dplyr::ungroup()

          df$counts <- df$counts / df$group_counts
          x_lab <- 'Number of residues'
          title <- paste0('Number of residues', ' per ', feature, ' for ', names(df.list)[i])

        }else{
          x_lab <- 'Proportion of residues'
          title <- paste0('Proportion of residues', ' per ', feature, ' for ', names(df.list)[i])
        }

        if(grouping != 'none'){
          df <- df %>%
                dplyr::group_by(grouping) %>%
                dplyr::distinct(feature, .keep_all = T)
        }else{
          df <- df %>%
                dplyr::distinct(feature, .keep_all = T)
        }



        out_plots[[i]] <- ggplot2::ggplot(df, ggplot2::aes(x = counts, y = feature, fill = feature)) +
                          ggplot2::geom_bar(show.legend = F, stat = "identity", width=0.6, color="black") +
                          ggplot2::labs(title = title, x = x_lab, y = feature) +
                          ggplot2::theme(panel.background = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none") +
                          ggplot2::coord_flip()

        if(grouping != 'none'){
          out_plots[[i]] <- out_plots[[i]] + ggplot2::facet_grid(~grouping, scales = 'free')
        }
      }
    }

    names(out_plots) <- names(df.list)
    return(out_plots)
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

  atom_df_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure)$atom)
  id_df <- names(atom_df_list)
  id_df <- id_df %>% stringr::str_split("\\.") %>% do.call(rbind, .) %>% data.frame()
  colnames(id_df) <- c('sample', 'clonotype', 'id')

  if(plot.level == 'global'){
    df_list <- list()
    df_list[[1]] <- do.call(rbind, atom_df_list)
    df_list[[1]]$lineplot_group <- 'global'

    names(df_list) <- 'global'

  }else if(plot.level == 'sample'){
    unique_samples <- unique(id_df$sample)
    df_list <- vector(mode = 'list', length = length(unique_samples))

    for(i in 1:length(unique_samples)){
      df_list[[i]] <- do.call(rbind, atom_df_list[which(id_df$sample == unique_samples[i])])
      df_list[[i]]$lineplot_group <- unique_samples[i]
    }
    names(df_list) <- unique_samples

  }else if(plot.level == 'clonotype'){
    sample_clonotypes <- paste0(id_df$sample, '.', id_df$clonotype)
    id_df$sample_clonotypes <- sample_clonotypes
    unique_sample_clonotypes <- unique(sample_clonotypes)

    df_list <- vector(mode = 'list', length = length(unique_sample_clonotypes))

    for(i in 1:length(unique_sample_clonotypes)){
      df_list[[i]] <- do.call(rbind, atom_df_list[which(id_df$sample_clonotypes == unique_sample_clonotypes[i])])
      df_list[[i]]$lineplot_group <- unique_sample_clonotypes[i]
    }
    names(df_list) <- unique_sample_clonotypes

  }else{
    df_list <- atom_df_list
    df_list <- lapply(1:length(df_list), function(i) df_list[[i]]$lineplot_group <- names(df_list)[i])
  }



  if(plot.format == 'bar'){
    out_plots <- create_barplot(df.list = df_list,
                                feature = feature,
                                grouping = grouping,
                                proportions = proportions)

  }else if(plot.format == 'line'){
    out_plots <- create_lineplot(df.list = df_list,
                                 feature = feature,
                                 combine.lineplot = combine.lineplot)

  }else{
    stop('Unrecognized plot format!')
  }


  return(out_plots)
}
