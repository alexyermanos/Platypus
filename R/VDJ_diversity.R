#' Calculates and plots common diversity and overlap measures for repertoires and alike. Require the vegan package
#' @param VDJ.matrix VDJ dataframe output from either the VDJ_analyse (platypus.version = "v2") or from the VDJ_GEX_matrix function (platypus.version = "v3")
#' @param feature.columns Character vector. One or more column names from the VDJ.matrix of which diversity or overlap metrics are calculated. if more than one column is provided (e.g. c("VDJ_cdr3s_aa","VJ_cdr3s_aa")) these columns will be pasted together before metric calculation. Defaults to "CDRH3_aa" if platypus.version == "v2" and "VDJ_cdr3s_aa" if platypus.version == "v3".
#' @param grouping.column Character. Column name of a column to group metrics by. This could be "sample_id" to calculate the metric for each sample. This column is required if metric = "simpson". If so, the simpson overlap index will be calculated pairwise for all combinations of elements in the grouping.column. Defaults to "none".
#' @param metric Character. Diversity or overlap metric to calculate. Can be c("richness", "remyparker", "simpson", "ginisimpson", "shannon", "shannonevenness", "jaccard"). Defaults to "shannon". If jaccard is selected, a heatmap with the pairwise comparisons between all groups is returned. If any of the others is selected, a dotplot is returned
#' @param pvalues.label.size Numeric. Only used if overlap indices are calculated. Defaults to 4. Is passed on to ggplot theme
#' @param axis.label.size Numeric. Only used if overlap indices are calculated. Defaults to 12. Is passed on to ggplot theme
#' @param platypus.version Version of platypus to use. Defaults to "v2". If an output of the VDJ_analyze function is supplied, set to "v2". If an output of the VDJ_GEX_matrix function is supplied set to "v3"
#' @return Returns a ggplot with the calculated metric for each group (if provided). Data is accessible via ggplot.output$data 
#' @export
#' @examples
#' \dontrun{
#' cluster.distribution.per.sample <- GEX_cluster_membership_per_sample(GEX.output=automate_GEX_out[[i]])
#'}
VDJ_diversity <- function(VDJ.matrix,
                          feature.columns,
                          grouping.column,
                          metric,
                          pvalues.label.size,
                          axis.label.size,
                          platypus.version){
  
  #start with vegan package
  require(vegan)
  
  if(missing(grouping.column)) grouping.column <- "none"
  if(missing(platypus.version)) platypus.version <- "v2"
  if(missing(metric)) metric <- "shannon"
  if(missing(feature.columns)){
    if(platypus.version == "v2"){feature.columns <- "CDR3H_aa"}
    else if(platypus.version == "v3"){feature.columns <- "VDJ_cdr3s_aa"}
  }
  if(missing(pvalues.label.size)) pvalues.label.size <- 4
  if(missing(axis.label.size)) axis.label.size <- 12
  
  for(i in 1:length(feature.columns)){
    if(!feature.columns[i] %in% names(VDJ.matrix)){stop("Please provide valid feature column name(s) contained within VDJ.matrix")} 
  }
  if(grouping.column != "none" & !grouping.column %in% names(VDJ.matrix)){stop("The provided grouping.column was not found in VDJ.matrix. Please provide a valid name or 'none' to avoid grouping")}
  
  #remove any rows that do not contain an entry for a given feature
  to_remove <- c()
  for(n in 1:nrow(VDJ.matrix)){
    if("" %in% VDJ.matrix[n,c(feature.columns)]){
      to_remove <- c(to_remove, n)}
  }
  VDJ.matrix <- VDJ.matrix[-to_remove,]
  
  
  #get basic dataframe with pasted features and group
  if(grouping.column == "none"){
    grouping <- data.frame("group" = rep(1, nrow(VDJ.matrix)))
  } else {
    grouping <- data.frame("group" = VDJ.matrix[, grouping.column])
  }
  
  if(length(feature.columns) > 1){
    grouping$pasted <- do.call(paste, c(VDJ.matrix[, c(feature.columns)], sep="/"))
  } else {
    grouping$pasted <- VDJ.matrix[, c(feature.columns)]
  }

  group.names <- unique(grouping$group)


  #loop over group.names => get index for every group => return table to use for geom_point()
  #to check if a correct metric was selected or whether to go for overlap metrics
  if(metric %in% c("richness", "remyparker", "simpson", "ginisimpson", "shannon", "shannonevenness")){
  
  out <- c()
  for(i in 1:length(group.names)){
  
  #get frequences
  freq <- table(subset(grouping, group == group.names[i])$pasted)
    
  if(metric == "richness"){
  #species richness  "richness"
  out <- c(out,exp(vegan::renyi(freq, scales = 0, hill= F)))
  title_out <- "Species richness"
  } else if (metric == "remyparker"){
  #remy parker "remyparker"
  out <- c(out,max(freq)/sum(freq))
  title_out <- "Remy-Parker index"
  } else if (metric == "simpson"){
  #simpson "simpson"
  out <- c(out,-diversity(freq, "simpson") + 1) #read gini simpson to understand
  title_out <- "Simpson index"
  } else if (metric == "ginisimpson"){ 
  #gini simpson "ginisimpson"
  out <- c(out,diversity(freq, "simpson")) #returns 1-simpson dominance which is the gini simpson index
  #this also works: 1/exp(vegan::renyi(freq, scales = 2, hill= F))
  title_out <- "Gini-Simpson index"
  } else if (metric == "shannon"){ 
  #shannon "shannon"
  out <- c(out,vegan::renyi(freq, scales = 1, hill= F))
  title_out <- "Shannon diversity"
  } else if (metric == "shannonevenness"){ 
  #shannon evenness "shannoneveness"
  out <- c(out,exp(diversity(freq, "shannon")) /length(freq))
  title_out <- "Shannon evenness"
  }
  }
  
  out_df <- data.frame("groups" = group.names, "metric" = out, colors = rainbow(length(group.names)))

  #plot
  plot_out <- ggplot(out_df, aes(x = groups, y = metric, col = colors)) + geom_point(size = 3,show.legend = F) + labs(title = title_out, x = "", y = title_out) + theme(panel.background = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")
  
  return(plot_out)
  
  } else if(metric %in% c("jaccard")){#else if the metric was not in the list of covered above
  
    combs <- as.data.frame(t(combn(group.names, m = 2,simplify = TRUE)))#get combinations to test
    combs[,1] <- ordered(as.factor(combs[,1]), levels = rev(group.names))
    combs[,2] <- ordered(as.factor(combs[,2]), levels = group.names)
    
    
  if(metric == "jaccard"){ #redundant at the moment but not if we add more indices here
   if(length(group.names) == 1){stop("Number of groups must be > 1 to calculate jaccard index")}
    
    combs$metric <- NA
    for(i in 1:nrow(combs)){
      inters <- length(intersect(unique(subset(grouping,group == combs[i,1])$pasted), unique(subset(grouping,group == combs[i,2])$pasted)))
      uni <- length(unique(subset(grouping,group == combs[i,1])$pasted)) + length(unique(subset(grouping,group == combs[i,2])$pasted)) - inters
      combs$metric[i] <- round(inters / uni, 2)
    }
    title_out <- "Jaccard index"
  }  
    
    plot_out <- ggplot(combs, aes(x = combs[,2], y = combs[,1],fill=metric)) + geom_tile() +geom_text(aes(label=metric), size = pvalues.label.size)+ scale_fill_gradient2(low="navy", mid="white", high="red", limits=range(combs$metric)) + theme(panel.background = element_blank(),axis.text = element_text(size = 30), axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks = element_blank(), text = element_text(size=30), legend.key = element_rect(colour = "white"), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25), plot.subtitle = element_text(size = 15),axis.text.x = element_text(angle = 60,vjust = 1, hjust=1, size = axis.label.size),axis.text.y = element_text(size = axis.label.size)) + labs(title = title_out, x = "", y = "", fill = "")
    
    return(plot_out)
  } else {stop("Please input a metric from the available selection listed in the doc")}
}


