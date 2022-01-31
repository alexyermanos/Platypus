#' Wrapper function for VDJ_abundances to obtain ranked clonotype barplots

#'@description Wraps the VDJ_abundances function and output a barplot of clonotypes ranked by expansion (x axis) with counts of the specific feature values per clonotype (y axis). For a more in-depth configuration of the barplots (e.g., including clonotypes with missing features, different strategies for NA values, etc.), use the VDJ_abundances function with output.format='plots'.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param features string or vector of strings, denoting the columns of the VDJ/VDJ.GEX.matrix[[1]] object from which to extract the unique feature values.
#' @param count.level string, 'absolute' will return the absolute counts, 'group.level.proportions' will return the counts divided by the total number or elements/values in the specific groups (group level proportions), 'sample.level.proportions' will return the counts divided by the total number of elements in the sample.
#' @param max.clonotypes integer or NULL, the maximum number of clonotypes for which to count features. If NULL, it will count for all clonotypes.
#' @param rank.clonotypes boolean, if T - clonotypes will be ranked and order according to their expansion.
#' @param specific.feat.colors named list (or NULL) of specific colors to be used in the final barplots.
#' @return Either a count dataframe with the following columns: group(=unique group value, e.g., 'clonotype1' if grouping.column='clonotype_id'), sample, group_frequency, unique_feature_values, feature_value_counts, total_feature_names
#'or a barplot of the counts/proportions per feature, per group.
#' @export


#' @examples
#' \dontrun{
#' VDJ_clonal_expansion_abundances(VDJ = small_vgm[[1]],
#' features='VDJ_cgene',count.level='absolute',
#' max.clonotypes=30, rank.clonotypes=T, specific.feat.colors=NULL)
#'}

VDJ_clonal_expansion_abundances <- function(VDJ,
                                  features, #specific feature column from the VDJ
                                  count.level, #either absolute counts (absolute), sample.level.proportions, or group.level.proportions
                                  max.clonotypes,
                                  rank.clonotypes,
                                  specific.feat.colors #list = custom color palette
                                  ){


   if(missing(VDJ)) stop('Please input your data as a VDJ matrix')
   if(missing(features)) features <- 'VDJ_cgene'
   if(missing(count.level)) count.level <- 'sample.level.proportions'
   if(missing(max.clonotypes)) max.clonotypes <- 30
   if(missing(rank.clonotypes)) rank.clonotypes <- T
   if(missing(specific.feat.colors)) specific.feat.colors <- NULL

   VDJ.matrix <- VDJ
   VDJ <- NULL

   #Global variable definition for CRAN checks
   unique_feature_values <- NULL
   feature_value_counts <- NULL
   Ranks <- NULL
   group_frequency <- NULL

   clonotype_df  <- VDJ_abundances(VDJ.matrix,
                                 feature.columns=features,
                                 proportions=count.level,
                                 specific.features=NULL,
                                 grouping.column='clonotype_id',
                                 max.groups=max.clonotypes,
                                 specific.groups='none',
                                 sample.column='sample_id',
                                 VDJ.VJ.1chain=T,
                                 treat.incomplete.groups='exclude',
                                 treat.incomplete.features='exclude',
                                 combine.features=F,
                                 treat.combined.features='exclude',
                                 specific.feature.colors = NULL,
                                 output.format='abundance.df')


   sample_ids <- unique(clonotype_df$sample)
   sample_dfs <- list()
   clonotype_plots <- list()
   for(i in 1:length(sample_ids)){
     sample_dfs[[i]] <- clonotype_df[which(clonotype_df$sample==sample_ids[i]),]
     sample_dfs[[i]] <- sample_dfs[[i]][order(sample_dfs[[i]]$group_frequency, decreasing=T),]

     if(rank.clonotypes){
       ranks <- 1:length(unique(sample_dfs[[i]]$group))
       unique_groups <- unique(sample_dfs[[i]]$group)
       sample_dfs[[i]]$Ranks <- rep(NA, nrow(sample_dfs[[i]]))
       for(j in 1:length(unique_groups)){
         sample_dfs[[i]]$Ranks[which(sample_dfs[[i]]$group==unique_groups[j])] <- ranks[j]
       }
       sample_dfs[[i]]$Ranks <- as.factor(sample_dfs[[i]]$Ranks)
     }

    if(rank.clonotypes){
      clonotype_plots[[i]] <-  ggplot2::ggplot(sample_dfs[[i]], ggplot2::aes(fill=unique_feature_values, y=feature_value_counts, x=Ranks)) +
                        ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank()) +
                        ggplot2::labs(fill=paste0(features, collapse='/'), y='Cells') + ggplot2::ggtitle(paste0(sample_ids[i]))

    }else{
      clonotype_plots[[i]] <-  ggplot2::ggplot(sample_dfs[[i]], ggplot2::aes(fill=unique_feature_values, y=feature_value_counts, x=group_frequency)) +
                        ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank()) +
                        ggplot2::labs(fill=paste0(features, collapse='/'), y='Cells') + ggplot2::ggtitle(paste0(sample_ids[i]))

    }

    if(!is.null(specific.feat.colors)){
      clonotype_plots[[i]] <- clonotype_plots[[i]] + ggplot2::scale_fill_manual(values=specific.feat.colors)
    }

    if(length(unique(sample_dfs[[i]]$feature_name))!=1){
      clonotype_plots[[i]] <- clonotype_plots[[i]] +  ggplot2::facet_wrap(~feature_name, scales = "free_x")
    }

    if(count.level!='absolute'){
      clonotype_plots[[i]] <- clonotype_plots[[i]] + ggplot2::labs(y='Proportions')
    }

   }

  return(clonotype_plots)

}
