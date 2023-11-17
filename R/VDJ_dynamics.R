#' Tracks a specific VDJ column across multiple samples/timepoints.

#'@description Track a VDJ column across multiple samples or timepoints. Tracking consists of creating a per sample/timepoint dataframe of unique values for the VDJ column and their respective counts inside that timepoints/repertoire. Also creates alluvial plots to show the temporal dynamics of the tracked elements.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param columns.to.track string or list of strings - VDJ column with values to track (e.g., 'VDJ_cgene' will track the changes in isotype counts/proportions across multiple timepoints, defined by the timepoints.column). If two columns are provided and tracked, then a new values will be created by combining the values from each column.
#' @param starting.point.repertoire string or integer - the repertoire from which to start tracking (1 = will start at the first repertoire, 's3' will start at repertoire 's3').
#' @param track.all.elements boolean - if T (and track.only.common=F), it will track all elements across all repertoires/timepoints.
#' @param track.only.common boolean - if T (and track.all.elements=F), it will only track the common elements across all repertoires/timepoints.
#' @param max.elements.to.track integer or NULL - the maximum number of elements to track (elements are first sorted by frequency/abundance). If NULL, it will track all elements.
#' @param specific.elements.to.track vector of strings or NULL - specific elements we want tracked. If NULL, all elements will be tracked.
#' @param additional.grouping.column string or 'none' - VDJ column for calculating the frequency/counts of elements on a per-group level. If output.format='plot', each unique group will have its own bar plot of timepoints/repertoires (x axis) and feature counts (y axis). If NULL, no additional grouping will be done.
#' @param max.additional.groups integer or NULL - the maximum number of additional groups to consider (groups are first ordered by their frequency = total number of cells in that group in the VDJ matrix). If NULL, all groups will be considered.
#' @param specific.additional.groups vector of strings or NULL - specific grouping factors we want to consider. If NULL, all grouping factors will be considered.
#' @param timepoints.column string - VDJ column with either timepoints or repertoires across which we want to track our elements (usually 'sample_id').
#' @param proportions.level string - 'absolute.counts' for absolute counts, 'group' for per group proportions, 'repertoire' for per repertoire/timepoint proportions.
#' @param output.format string - 'plot' for alluvial barplots, 'df' for count/proportions dataframes of the tracked elements.
#' @param ignore.legend boolean - if T, the legend will not be included in the resulting ggplot object.


#' @return Either a count dataframe of the tracked elements across multiple timepoints/repertoires, or alluvial barplot.
#' @export

#' @examples
#' VDJ_dynamics(VDJ = small_vgm[[1]],columns.to.track='clonotype_id', starting.point.repertoire=1,
#' max.elements.to.track=10, timepoints.column='sample_id',
#' output.format='plot')
#'

VDJ_dynamics <- function(VDJ,
                      columns.to.track,
                      starting.point.repertoire,
                      track.all.elements,
                      track.only.common,
                      max.elements.to.track,
                      specific.elements.to.track,
                      additional.grouping.column,
                      max.additional.groups,
                      specific.additional.groups,
                      timepoints.column,
                      proportions.level,
                      output.format,
                      ignore.legend){

  #global variable definitions for CRAN checks
  unique_feature_values <- NULL

  if(missing(VDJ)) stop('Please input your data as a VDJ object')
  if(missing(columns.to.track)) columns.to.track <- 'VDJ_cdr3s_aa'
  if(missing(starting.point.repertoire)) starting.point.repertoire <- 1 #look for the n most abundant max elements to track from the first repertoire
  if(missing(track.all.elements)) track.all.elements <- F
  if(missing(track.only.common)) track.only.common <- F
  if(missing(max.elements.to.track)) max.elements.to.track <- 20
  if(missing(specific.elements.to.track)) specific.elements.to.track <- NULL
  if(missing(additional.grouping.column)) additional.grouping.column <- NULL
  if(missing(max.additional.groups)) max.additional.groups <- NULL
  if(missing(specific.additional.groups)) specific.additional.groups <- NULL
  if(missing(timepoints.column)) timepoints.column <- 'sample_id'
  if(missing(proportions.level)) proportions.level <- 'repertoire' # or absolute.counts or proportions.group (if additional.grouping.column)
  if(missing(output.format)) output.format <- 'plot'
  if(missing(ignore.legend)) ignore.legend <- F
  
  if(track.only.common){
    starting.point.repertoire <- NULL
  }

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

  if(('CDR3aa' %in% columns.to.track) & !('CDR3aa' %in% colnames(VDJ))){
    VDJ$CDR3aa <- mapply(function(x,y) get_feature_combinations(x,y,split.x=T,split.y=T, combine.sequences=T), VDJ$VDJ_cdr3s_aa, VDJ$VJ_cdr3s_aa)
  }

  for(i in 1:length(columns.to.track)){
   if(!columns.to.track[i] %in% names(VDJ)){
     stop("Please provide valid feature column name(s) contained within VDJ")
   }
  }


  repertoire_numbers <- unique(VDJ[,timepoints.column])

  unique_groups <- list()
  if(!is.null(additional.grouping.column)){
    if(!is.null(specific.additional.groups)){
      unique_groups <- specific.additional.groups
    }else{
      unique_groups <- unlist(unique(VDJ[,additional.grouping.column]))
      unique_group_frequencies <- unlist(lapply(unique_groups, function(x) length(which(VDJ[,additional.grouping.column]==x))))
      unique_groups <- unique_groups[order(unique_group_frequencies, decreasing=T)]
    }
    if(!is.null(max.additional.groups)){
      if(max.additional.groups<length(unique_groups)){
        unique_groups <- unique_groups[1:max.additional.groups]
      }
    }
    unique_groups <- unique_groups[which(!is.null(unique_groups) & !is.na(unique_groups) & unique_groups!='')]
  }

  if(length(unique_groups)==0){
    additional.grouping.column <- 'none'
    unique_groups <- 'none'
  }

  tracked_dfs <- list()

  for(i in 1:length(unique_groups)) {

    if(!is.null(specific.elements.to.track)){
      if(inherits(specific.elements.to.track,'data.frame')) {unique_elements_to_track <- unique(paste0(specific.elements.to.track[,1], ' ', specific.elements.to.track[,2]))}
      else {unique_elements_to_track <- unique(specific.elements.to.track)}
    }

    if(!is.null(starting.point.repertoire)){
      if(unique_groups[i]!='none'){ starting_rep <- VDJ[which(VDJ[additional.grouping.column]==unique_groups[i]), ]
      }else{ starting_rep <- VDJ }

      if(inherits(starting.point.repertoire,'numeric')){
        starting_rep <- starting_rep[which(starting_rep[timepoints.column]==repertoire_numbers[starting.point.repertoire]),]
      }else{
        starting_rep <- starting_rep[which(starting_rep[timepoints.column]==starting.point.repertoire),]
      }

      if(length(columns.to.track)==2){
       combined_features <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T)} else '', starting_rep[,columns.to.track[[1]]], starting_rep[,columns.to.track[[2]]])
       unique_elements_to_track <- unlist(unique(combined_features))
       unique_elements_frequencies <- unlist(lapply(unique_elements_to_track, function(x) length(which(combined_features==x))))
       unique_elements_to_track <- unique_elements_to_track[order(unique_elements_frequencies, decreasing=T)]

     }else{
       unique_elements_to_track <- unlist(unique(starting_rep[,columns.to.track]))
       unique_elements_frequencies <- unlist(lapply(unique_elements_to_track, function(x) length(which(starting_rep[, columns.to.track]==x))))
       unique_elements_to_track <- unique_elements_to_track[order(unique_elements_frequencies, decreasing=T)]
     }
   }else if(track.all.elements==T){

     if(length(columns.to.track)==2){
      combined_features <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T)} else '', VDJ[,columns.to.track[[1]]], VDJ[,columns.to.track[[2]]])
      unique_elements_to_track <- unlist(unique(combined_features))
      unique_elements_frequencies <- unlist(lapply(unique_elements_to_track, function(x) length(which(combined_features==x))))
      unique_elements_to_track <- unique_elements_to_track[order(unique_elements_frequencies, decreasing=T)]

    }else{
      unique_elements_to_track <- unlist(unique(VDJ[,columns.to.track]))
      unique_elements_frequencies <- unlist(lapply(unique_elements_to_track, function(x) length(which(VDJ[, columns.to.track]==x))))
      unique_elements_to_track <- unique_elements_to_track[order(unique_elements_frequencies, decreasing=T)]
    }

   }else if(track.only.common==T){
     unique_per_rep <- list()
     for(j in 1:length(repertoire_numbers)) {
       VDJ_subset <- VDJ[which(VDJ[timepoints.column]==repertoire_numbers[j]),]

       if(length(columns.to.track)==2){
        combined_features <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T)} else '', VDJ_subset[,columns.to.track[[1]]], VDJ_subset[,columns.to.track[[2]]])
        unique_elements_to_track <- unlist(unique(combined_features))

       }else{
        unique_elements_to_track <- unlist(unique(VDJ_subset[,columns.to.track]))
      }
       unique_per_rep[[j]] <- unique_elements_to_track
     }
     unique_elements_to_track <- Reduce(intersect, unique_per_rep)
   }

   if(!is.null(max.elements.to.track)){
     if(max.elements.to.track<length(unique_elements_to_track)){
       unique_elements_to_track <- unique_elements_to_track[1:(max.elements.to.track+1)]
     }
   }

   unique_elements_to_track <- unlist(unique_elements_to_track)
   if(length(columns.to.track)!=2){
     tracked_dfs[[i]] <- VDJ_abundances(VDJ, feature.columns=as.list(columns.to.track), proportions='absolute', specific.features=unique_elements_to_track,
                                        grouping.column = additional.grouping.column, max.groups=NULL, specific.groups=unique_groups[i], sample.column=timepoints.column, treat.incomplete.groups='exclude', treat.incomplete.features='exclude', treat.combined.groups = 'exclude', output.format='abundance.df')

   }else{
     tracked_dfs[[i]] <- VDJ_abundances(VDJ, feature.columns=as.list(columns.to.track), proportions='absolute', specific.features=unique_elements_to_track, combine.features=T,
                                        grouping.column = additional.grouping.column, max.groups=NULL, specific.groups=unique_groups[i], sample.column=timepoints.column, treat.incomplete.groups='exclude', treat.incomplete.features='exclude', treat.combined.features='exclude', treat.combined.groups = 'exclude', output.format='abundance.df')

   }

    if(proportions.level=='repertoire'){
      tracked_dfs[[i]]$proportions <- tracked_dfs[[i]]$feature_value_counts / tracked_dfs[[i]]$sample_frequency
    }else if(proportions.level=='group' & additional.grouping.column!='none'){
      tracked_dfs[[i]]$proportions <- tracked_dfs[[i]]$feature_value_counts / tracked_dfs[[i]]$group_frequency
    }else{
      tracked_dfs[[i]]$proportions <- tracked_dfs[[i]]$feature_value_counts
    }
  }
  tracked_df <- do.call('rbind', tracked_dfs)


  output_plot <- ggplot2::ggplot(tracked_df, ggplot2::aes(x = sample, y = proportions,
    fill = unique_feature_values, stratum = unique_feature_values, alluvium = unique_feature_values, label = unique_feature_values)) +
    ggalluvial::geom_alluvium() + ggalluvial::geom_stratum() + cowplot::theme_cowplot() + ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank()) +
    ggplot2::labs(fill=paste0(columns.to.track, collapse='/'))

  if(additional.grouping.column!='none'){
    output_plot <- output_plot + ggplot2::facet_wrap(~group, scales = "free_x")
  }
  if(proportions.level=='absolute.counts'){
    output_plot <- output_plot + ggplot2::labs(y='absolute counts')
  }
  if(ignore.legend==T){
    output_plot <- output_plot + ggplot2::theme(legend.position = "none")
  }

  if(output.format=='plot'){
    return(output_plot)
  }else if(output.format=='df'){
    return(tracked_df)
  }else{
    stop('Wrong output.format value!')
  }
}
