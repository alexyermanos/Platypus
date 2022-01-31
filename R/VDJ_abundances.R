#' Calculate abundances/counts of specific features for a VDJ dataframe

#'@description Calculate the absolute counts or proportions of a specific cell-level feature (column in the VDJ/VDJ.GEX.matrix[[1]] object), per an optional specific grouping factor (e.g., clonotype via 'clonotype_id') and an optional sample factor(e.g., 'sample_id'). Outputs either a count dataframe of the specific feature or a ggplot2 barplot.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param feature.columns vector of strings, denoting the columns of the VDJ/VDJ.GEX.matrix[[1]] object from which to extract the unique feature values (for which we will calculate the counts or proportions).
#' @param proportions string, 'absolute' will return the absolute counts, 'group.level.proportions' will return the counts divided by the total number or elements/values in the specific groups (group level proportions), 'sample.level.proportions' will return the counts divided by the total number of elements in the sample.
#' @param specific.features vector of specific feature values (or NULL) for which to calculate counts/proportions, from the specified feature.columns parameter (only works if a single feature column is specified in feature.columns).
#' @param grouping.column string or 'none', represents the column from the VDJ/VDJ.GEX.matrix[[1]] object by which to group counting process. This is usually the 'clonotype_id' column to calculate frequencies at the clonotype level. If 'none', no grouping will be done.
#' For example, if feature.columns='VDJ_cgene' and grouping.column='clonotype_id', we will obtain a count dataframe of the frequencies of each isotype per unique clonotype (per sample if sample.column='sample_id').
#' @param max.groups integer or NULL, the maximum number of groups for which to count features. If NULL, it will count for all groups.
#' @param specific.groups vector of strings (or 'none'), if the counting should be done only for specific groups (e.g., count the frequency of isotype only for clonotypes 1 and 2 if feature.columns='VDJ_cgene', grouping.column='clonotype_id' and specific.groups=c('clonotype1', 'clonotype2'))
#' @param sample.column string, represents the sample column if your VDJ/VDJ.GEX.matrix[[1]] object has multiple samples (usually 'sample_id')
#' @param VDJ.VJ.1chain boolean, if T will remove aberrant cells (more than 1 VDJ of VJ chain), if F it will keep them.
#' @param treat.incomplete.groups string, method of dealing with groups which are missing the features in the feature.columns parameter (e.g., a clonotype which does not have any transcriptomic clusters annotations if feature.columns='transcript_cluster').'exclude' - excludes groups with no cells for the specific features, 'unknown' - sets them as unknown
#' @param treat.incomplete.features string, method of dealing with missing feature values (e.g., a clonotype has several NA values for the 'VDJ_cgene' feature.column - cells with NA values). 'unknown' - counted as unknown, 'exclude' - excludes completely, 'max.global' - replaces value by max value of that feature across the repertoire, 'max.group' - replaced by the max feature value inside that group, 'proportional' - iteratively assigns the missing values to the known groups, keeping the same proportions.
#' @param combine.features boolean - if T and we have two columns in feature.columns, will combine the feature values for each cell in the VDJ object, counting them as a single feature when calculating proportions.
#' @param treat.combined.features string, method of dealing with combined features with missing values. 'exclude' will be treated similarly to excluding incomplete feature values (excluding them completely if a single value is missing from the combination), or 'include' and will be treated as a new feature value.
#' @param specific.feature.colors named list of specific colors to be used in the final barplots, for each unique feature value in the VDJ object's feature.columns values.
#' For example, if we have a feature column of binders with unique values=c('yes', 'no'), specific.feature.colors=list('yes'='blue', 'no'='red') will color them accordingly.
#' @param output.format string, either 'plots' to obtain barplots, 'abundance.df' to obtain the count dataframe, or 'abundance.df.list' to obtain a list of count dataframes, for each sample.
#' @return Either a count dataframe with the following columns: group(=unique group value, e.g., 'clonotype1' if grouping.column='clonotype_id'), sample, group_frequency, unique_feature_values, feature_value_counts, total_feature_names
#'or a barplot of the counts/proportions per feature, per group.
#' @export
#' @examples
#' VDJ_abundances(VDJ = small_vgm[[1]],
#' feature.columns='VDJ_cgene',proportions='absolute',
#' grouping.column='clonotype_id',specific.groups='none',
#' output.format='abundance.df')
#'

VDJ_abundances <- function(VDJ,
                           feature.columns,
                           proportions,
                           specific.features,
                           grouping.column,
                           max.groups,
                           specific.groups,
                           sample.column,
                           VDJ.VJ.1chain,
                           treat.incomplete.groups,
                           treat.incomplete.features,
                           combine.features,
                           treat.combined.features,
                           specific.feature.colors,
                           output.format){

  if(missing(VDJ)) stop('Please input your data as a VDJ matrix')
  if(missing(feature.columns)) feature.columns <- c('VDJ_cgene')
  if(missing(proportions)) proportions <- 'absolute'
  if(missing(specific.features)) specific.features <- NULL
  if(missing(grouping.column)) grouping.column <- 'clonotype_id'
  if(missing(max.groups)) max.groups <- NULL
  if(missing(specific.groups)) specific.groups <- 'none'
  if(missing(sample.column)) sample.column <- 'sample_id'
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
  if(missing(treat.incomplete.groups)) treat.incomplete.groups <- 'exclude' #exclude - excludes groups with no cells for the specific features, #unknown - sets them as unknown
  if(missing(treat.incomplete.features)) treat.incomplete.features <- 'exclude'
  if(missing(combine.features)) combine.features <- F
  if(missing(treat.combined.features) & combine.features==T) treat.combined.features <- 'exclude' #will be treated similarly to incomplete features, or include and will be treated as a new feature
  if(missing(specific.feature.colors)) specific.feature.colors <- NULL
  if(missing(output.format)) output.format <- 'plots' #or abundance.df

  #Global variable definitions for CRAN checks
  unique_feature_values <- NULL
  unique_feature_counts <- NULL
  group <- NULL
  feature_value_counts <- NULL


  ###############################UTILITY FUNCTIONS FOR LAPPLY###########################################
  get_feature_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';', combine.sequences=F){
   if(split.x==T) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
   if(split.y==T) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

   ccombs <- expand.grid(x,y)
   if(!combine.sequences){
     ccombs<-paste0(ccombs[,1], ' ', ccombs[,2])
   }else{
     ccombs<-paste0(ccombs[,1], ccombs[,2])

   }
   ccombs <- paste0(ccombs, collapse=collapse.by)

   return(ccombs)
  }



  get_count_df_for_group <- function(df, group, feature){

   group_frequency <- length(which(df[,grouping.column]==group))
   sample_frequency <- nrow(df)
   group_df <- df[which(df[grouping.column]==group),]

   if(!is.null(specific.features)){
     unique_feature_values <- specific.features
   }else{
     unique_feature_values <- unique(df[feature][which(df[,grouping.column]==group),])
     unique_feature_values <- unlist(lapply(unique_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     unique_feature_values <- unique(unique_feature_values)

   }


   if((all(is.na(unique_feature_values)) | all(is.null(unique_feature_values)) | all(unique_feature_values=='') | group=='' | is.na(group) | is.null(group)) & treat.incomplete.groups=='exclude'){
     return(NULL)
   }else if((all(is.na(unique_feature_values==0)) | all(is.null(unique_feature_values)) | all(unique_feature_values=='') | group=='' | is.na(group) | is.null(group)) & treat.incomplete.groups=='unknown'){
     output_df <- data.frame(group = c(group),
                             sample = c(sample_id[i]),
                             sample_frequency = matrix(sample_frequency),
                             group_frequency = matrix(group_frequency),
                             unique_feature_values = c('unknown'),
                             feature_value_counts = matrix(group_frequency),
                             feature_name = feature)
  }



  if(any(is.na(unique_feature_values)) | any(is.null(unique_feature_values)) | any(unique_feature_values=='')){
    missing_indices <- which(is.na(unique_feature_values) | is.null(unique_feature_values) | unique_feature_values=='')

  if(treat.incomplete.features=='exclude'){
     group_df <- group_df[which(!(is.na(group_df[feature])) & !(is.null(group_df[feature])) & group_df[feature]!=''),]
     group_feature_values <- group_df[,feature]

     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }
     unique_feature_values <- unique_feature_values[which(!is.na(unique_feature_values) & !is.null(unique_feature_values) & unique_feature_values!='')]
     feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))

     group_frequency <- nrow(group_df)

  }else if(treat.incomplete.features=='max.global'){
     df_global <- df[which(!(is.na(df[feature])) & !(is.null(df[feature])) & df[feature]!=''),]
     all_global_features <- df_global[,feature]
     if(any(stringr::str_detect(all_global_features, ';'))){
       all_global_features <- unlist(lapply(all_global_features, function(x) unlist(stringr::str_split(x, ';'))))
     }
     unique_global_features <- unique(all_global_features)

     global_feature_counts <- lapply(unique_global_features, function(x) length(which(all_global_features==x)))
     max_global_feature <- unique_global_features[which(global_feature_counts==max(unlist(global_feature_counts)))][1]

     group_df[feature][which(is.na(group_df[feature]) | is.null(group_df[feature]) | group_df[feature]==''),] <- max_global_feature
     group_feature_values <- group_df[,feature]
     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }

     unique_feature_values[missing_indices] <- max_global_feature
     unique_feature_values <- unique(unique_feature_values)

     feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))

     group_frequency <- nrow(group_df)

   }else if(treat.incomplete.features=='unknown'){
     group_df[feature][which(is.na(group_df[feature]) | is.null(group_df[feature]) | group_df[feature]==''),] <- 'unknown'

     group_feature_values <- group_df[,feature]
     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }

     unique_feature_values[missing_indices] <- 'unknown'
     unique_feature_values <- unique(unique_feature_values)

     feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))

     group_frequency <- nrow(group_df)

   }else if(treat.incomplete.features=='max.group'){

     group_df_non_null <- group_df[which(!(is.na(group_df[feature])) & !(is.null(group_df[feature])) & !(group_df[feature]=='')),]
     group_feature_values <- unlist(group_df_non_null[,feature])

     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }
     non_null_feature_values <- unique(group_feature_values)

     feature_counts <- lapply(non_null_feature_values, function(x) length(which(group_feature_values==x)))

     group_max_value <- non_null_feature_values[which(feature_counts==max(unlist(feature_counts)))][1]

     group_df[feature][which(is.na(group_df[feature]) | is.null(group_df[feature]) | group_df[feature]==''),] <- group_max_value
     group_feature_values <- unlist(group_df[,feature])
     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }
     unique_feature_values[missing_indices] <- group_max_value
     unique_feature_values <- unique(unique_feature_values)

     feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))

     group_frequency <- nrow(group_df)


   }else if(treat.incomplete.features=='proportional'){
     group_feature_values <- group_df[,feature]
     if(any(stringr::str_detect(group_feature_values, ';'))){
       group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
     }
     unique_feature_values <- unique(group_feature_values)
     group_frequency <- nrow(group_df)

     feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))
     feature_counts <- unlist(feature_counts)
     count_initial_proportions <- feature_counts/group_frequency

     while(sum(feature_counts) < group_frequency){
       for(i in 1:length(feature_counts)){
         feature_counts[i] <- feature_counts[i] + 1

         if(sum(feature_counts)==group_frequency){
           break
         }

         if( (feature_counts[i]/sum(feature_counts)) > count_initial_proportions[i]){
           next
         }
       }
     }
   }

   output_df <- data.frame(group = c(rep(group, length(unique_feature_values))),
                          sample = c(rep(sample_id[i], length(unique_feature_values))),
                          sample_frequency = matrix(rep(sample_frequency, length(unique_feature_values))),
                          group_frequency = matrix(rep(group_frequency, length(unique_feature_values))),
                          unique_feature_values = matrix(unique_feature_values),
                          feature_value_counts =  matrix(unlist(feature_counts)),
                          feature_name = c(rep(feature, length(unique_feature_values))))

  }else{

    group_feature_values <- group_df[,feature]
    if(any(stringr::str_detect(group_feature_values, ';'))){
      group_feature_values <- unlist(lapply(group_feature_values, function(x) unlist(stringr::str_split(x, ';'))))
    }

    feature_counts <- lapply(unique_feature_values, function(x) length(which(group_feature_values==x)))

    output_df <- data.frame(group = c(rep(group, length(unique_feature_values))),
                           sample = c(rep(sample_id[i], length(unique_feature_values))),
                           sample_frequency = matrix(rep(sample_frequency, length(unique_feature_values))),
                           group_frequency = matrix(rep(group_frequency, length(unique_feature_values))),
                           unique_feature_values = matrix(unique_feature_values),
                           feature_value_counts =  matrix(unlist(feature_counts)),
                           feature_name = c(rep(feature, length(unique_feature_values))))
  }

  if(proportions=='sample.level.proportions'){
    output_df$feature_value_counts <- output_df$feature_value_counts / output_df$sample_frequency
  }
  if(proportions=='group.level.proportions'){
    output_df$feature_value_counts <- output_df$feature_value_counts / output_df$group_frequency
  }

  return(output_df)
 }
  ###############################UTILITY FUNCTIONS FOR LAPPLY###########################################


  VDJ.matrix <- VDJ
  VDJ <- NULL
  sample_dfs <- list()

  if(('CDR3aa' %in% feature.columns) & !('CDR3aa' %in% colnames(VDJ.matrix))){
    VDJ.matrix$CDR3aa <- mapply(function(x,y) get_feature_combinations(x,y,split.x=T,split.y=T, combine.sequences=T), VDJ.matrix$VDJ_cdr3s_aa, VDJ.matrix$VJ_cdr3s_aa)
  }

  for(i in 1:length(feature.columns)){
   if(!feature.columns[i] %in% names(VDJ.matrix)){
     stop("Please provide valid feature column name(s) contained within VDJ")
   }
  }

  if(grouping.column != "none" & !(grouping.column %in% names(VDJ.matrix))){
   stop("The provided grouping.column was not found in VDJ. Please provide a valid name or 'none' to avoid grouping")
  }

  if(VDJ.VJ.1chain==T){
    VDJ.matrix <- VDJ.matrix[which(VDJ.matrix$Nr_of_VDJ_chains==1 & VDJ.matrix$Nr_of_VDJ_chains==1),]
  }

  if(grouping.column=='none'){
    VDJ.matrix$none <- rep('none', nrow(VDJ.matrix))
  }

  if(sample.column!='none'){
   sample_id <- unique(VDJ.matrix[,sample.column])
   for(i in 1:length(sample_id)){
     sample_dfs[[i]] <- VDJ.matrix[which(VDJ.matrix[,sample.column]==sample_id[i]),]
   }
  }else{
   sample_id <- 'global'
   sample_dfs[[1]] <- VDJ.matrix
  }


  if(length(feature.columns)==2 & combine.features==T){
   for(i in 1:length(sample_dfs)){
     if(treat.combined.features=='exclude'){
       combined_features <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T)} else '', sample_dfs[[i]][,feature.columns[[1]]], sample_dfs[[i]][,feature.columns[[2]]])
     }else{
       combined_features <- mapply(function(x,y) get_feature_combinations(x,y,split.x=T,split.y=T), sample_dfs[[i]][,feature.columns[[1]]], sample_dfs[[i]][,feature.columns[[2]]])
     }
     sample_dfs[[i]]$new_feature <- combined_features
     new_feature <- paste0(feature.columns[[1]], '/', feature.columns[[2]])
     names(sample_dfs[[i]])[names(sample_dfs[[i]])=='new_feature'] <- new_feature
   }
   feature.columns <- paste0(feature.columns[[1]], '/', feature.columns[[2]])
  }

  abundance_df_per_sample <- list()
  for(i in 1:length(sample_dfs)){
     if(specific.groups!='none'){
       unique_groups <- specific.groups
     }else{
       unique_groups <- unique(sample_dfs[[i]][,grouping.column])
     }

     if(!is.null(max.groups)){
       if(max.groups<length(unique_groups)){
         group_frequencies <- lapply(unique_groups, function(x) length(which(sample_dfs[[i]][,grouping.column]==x)))

         group_frequencies <- unlist(group_frequencies)

         sorted_unique_groups <- unlist(unique_groups)[order(group_frequencies, decreasing=T)]
         unique_groups <- sorted_unique_groups[1:max.groups]
       }
     }

     all_feature_dfs <- list()
     for(j in 1:length(feature.columns)){
       single_feature_dfs <- lapply(unique_groups, function(x) get_count_df_for_group(df=sample_dfs[[i]], group=x, feature=feature.columns[[j]]))
       single_feature_dfs <- single_feature_dfs[!sapply(single_feature_dfs,is.null)]
       all_feature_dfs[[j]] <- do.call('rbind', single_feature_dfs)
       #all_feature_dfs[[j]] <- single_feature_dfs
     }

     abundance_df_per_sample[[i]] <- do.call('rbind', all_feature_dfs)
  }

   abundance_df <- do.call('rbind', abundance_df_per_sample)


  if(output.format=='abundance.df'){
   return(abundance_df)
  }else if(output.format=='abundance.df.list'){
   return(abundance_df_per_sample)
 }else if(output.format=='plots'){
   sample_ids <- unique(abundance_df$sample)
   sample_dfs <- list()
   plots <- list()
   for(i in 1:length(sample_ids)){
     sample_dfs[[i]] <- abundance_df[which(abundance_df$sample==sample_ids[i]),]
     #sample_dfs[[i]] <- sample_dfs[[i]][order(sample_dfs[[i]]$group_frequency, decreasing=T),]

     #ranks <- 1:length(unique(sample_dfs[[i]]$group))
     #unique_groups <- unique(sample_dfs[[i]]$group)
     #sample_dfs[[i]]$Ranks <- rep(NA, nrow(sample_dfs[[i]]))
     #for(j in 1:length(unique_groups)){
      # sample_dfs[[i]]$Ranks[which(sample_dfs[[i]]$group==unique_groups[j])] <- ranks[j]
     #}
     #sample_dfs[[i]]$Ranks <- as.factor(sample_dfs[[i]]$Ranks)

     plots[[i]] <-  ggplot2::ggplot(sample_dfs[[i]], ggplot2::aes(fill=unique_feature_values, y=feature_value_counts, x=group)) +
                       ggplot2::geom_bar(stat="identity") + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank()) +
                       ggplot2::labs(fill=paste0(feature.columns, collapse='/'), y='Cells', x='Group') + ggplot2::ggtitle(paste0(sample_ids[i]))
    if(length(feature.columns)!=1){
      plots[[i]] <- plots[[i]] + ggplot2::facet_wrap(~feature_name, scales = "free_x")
    }

    if(!is.null(specific.feature.colors)){
      plots[[i]] <- plots[[i]] + ggplot2::scale_fill_manual(values=specific.feature.colors)
    }

    if(proportions!='absolute'){
      plots[[i]] <- plots[[i]] + ggplot2::labs(y='Proportions')
    }
  }
  return(plots)
 }
}
