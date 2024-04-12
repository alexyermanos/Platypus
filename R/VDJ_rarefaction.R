#' Plots rarefaction curves for species denoted in the feature.columns parameter across groups determined by grouping.columns

#' @param VDJ VDJ dataframe output from the VDJ_GEX_matrix function.
#' @param feature.columns Character vector. One or more column names from the VDJ to indicate the unique species for the rarefaction (to rarefy across). If more than one column is provided (e.g. c("VDJ_cdr3s_aa","VJ_cdr3s_aa")) these columns will be pasted together.
#' @param grouping.column Character. Column name of a column to group the rarefaction by. This could be "sample_id" for rarefaction curves for each sample.
#' @param VDJ.VJ.1chain Boolean defaults to TRUE. Whether to filter out aberrant cells (more than 1 VDJ or VJ chain).
#' @param rarefaction.type Character. Options for the iNEXT rarefaction - 'sample.size','coverage.based', or 'sample.completeness'.
#' @param hill.numbers Integer/ vector of integers. The Hill numbers to be plotted out (0 - species richness, 1 - Shannon diversity, 2 - Simpson  diversity)
#' @param number.resamples Integer. Number of bootstrap replications.
#' @param sample.sizes Vector if integers. The sample size points at which rarefaction should be performed. Defaults to NULL
#' @param endpoint Integer. The maximum sample size for rarefaction extrapolation. Defaults to NULL = 2 times the sample size for each sample.
#' @return Returns a ggplot with the ordination analysis performer across features, groups, or both
#' @export
#' @examples
#'
#' \donttest{
#' try({
#' plot <- VDJ_diversity(VDJ = Platypus::small_vgm[[1]],
#' ,feature.columns = c("VDJ_cdr3s_aa"), grouping.column = "sample_id")
#'})
#'}



VDJ_rarefaction <- function(VDJ,
                            feature.columns,
                            grouping.column,
                            VDJ.VJ.1chain,
                            rarefaction.type,
                            hill.numbers,
                            number.resamples,
                            sample.sizes,
                            endpoint
                            ){


  if(missing(VDJ)) stop('Please input your data as a VDJ matrix')
  if(missing(feature.columns)) feature.columns <- c('VDJ_cdr3s_aa')
  if(missing(grouping.column)) grouping.column <- 'sample_id' # or =='none'
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- FALSE
  if(missing(rarefaction.type)) rarefaction.type <- 'sample.size' #vs coverage.based and sample.completeness
  if(missing(hill.numbers)) hill.numbers <- c(0,1,2)
  if(missing(number.resamples)) number.resamples <- 50
  if(missing(sample.sizes)) sample.sizes <- NULL
  if(missing(endpoint)) endpoint <- NULL

  rarefaction_type_dict <- c(1:3)
  names(rarefaction_type_dict) <- list('sample.size', 'sample.completeness', 'coverage.based')
  diversity_dict <- list('Species richness', 'Shannon diversity', 'Simpson diversity')


  get_abundances <- function(VDJ, feature.columns, grouping.column, VDJ.VJ.1chain){

    if(length(feature.columns) > 1){
      combine.features <- TRUE
    }else{
      combine.features <- FALSE
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

    unique_groups <- unique(abundance_df$group)
    inext_rep_list <- vector(mode = 'list', length = length(unique_groups))

    for(i in 1:length(unique_groups)){
      abundance_df_group <- abundance_df[which(abundance_df$group==unique_groups[i]),]
      inext_rep_list[[i]] <- abundance_df_group$feature_value_counts
    }

    names(inext_rep_list) <- unique_groups

    return(inext_rep_list)
  }


  plot_rarefaction <- function(inext_rep_list){
    output_plots <- list()
    #requireNamespace('iNEXT')

    inext_object <- iNEXT::iNEXT(inext_rep_list, q=hill.numbers, datatype='abundance', nboot = number.resamples, endpoint = endpoint, size = sample.sizes)

    output_plots[[1]] <- iNEXT::ggiNEXT(inext_object, type=as.numeric(rarefaction_type_dict[rarefaction.type]), facet.var="Order.q") +
                         ggplot2::theme_bw() +
                         ggplot2::theme_classic() +
                         ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
                         ggplot2::labs(title = paste0('Rarefaction analysis of ', paste0(feature.columns, collapse = ' / ')), x = 'Number of sampled cells', y = paste0(paste0(feature.columns, collapse = ' / '), ' diversity')) +
                         ggplot2::guides(color = ggplot2::guide_legend(paste0(grouping.column, collapse = ' / ')),  shape = ggplot2::guide_legend(paste0(grouping.column, collapse = ' / ')),  fill = ggplot2::guide_legend(paste0(grouping.column, collapse = ' / '))) +
                         ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey70", size = 0.5, linetype = 2))



    output_plots[[2]] <- iNEXT::ggiNEXT(inext_object, type=as.numeric(rarefaction_type_dict[rarefaction.type]), facet.var="Assemblage") +
                         ggplot2::theme_bw() +
                         ggplot2::theme_classic() +
                         ggplot2::theme(panel.background = ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
                         ggplot2::labs(title = paste0('Rarefaction analysis of ', paste0(feature.columns, collapse = ' / '), ' per groups'), x = 'Number of sampled cells', y = paste0(paste0(feature.columns, collapse = ' / '), ' diversity')) +
                         ggplot2::guides(color = ggplot2::guide_legend('Hill numbers'), fill = ggplot2::guide_legend('Hill numbers'), shape = ggplot2::guide_legend('Hill numbers')) +
                         ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "grey70", size = 0.5, linetype = 2))



    return(output_plots)
  }


  plots <- VDJ %>% get_abundances(feature.columns = feature.columns, grouping.column = grouping.column, VDJ.VJ.1chain = VDJ.VJ.1chain) %>% plot_rarefaction()
  return(plots)
}
