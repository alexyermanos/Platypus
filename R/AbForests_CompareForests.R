#' Comparison of distinct B cell repertoires
#'
#' @description CompareForests takes the output of AntibodyForest for 2 distinct repertoires and performs a comparison of these 2 repertoires.
#' @param list1  a list of lists. Each sublist contains an igraph object with the networks of the evolved B clonal lineages in tree format, their legend and the number of randomly considered nodes per network for Repertoire of 1 (Output of AntibodyForest).
#' E.g list1[[1][[1]] is an igraph object, containing the first network of the evolved B clonal lineage in tree format.
#' list1[[1]][[2]] contains the legend parameters of the first network of the evolved B clonal lineage.
#' list1[[1]][[3]] is the number of randomly considered nodes for the first network of the evolved B clonal lineage.
#' @param list2 a list of lists. Each sublist contains an igraph object with the networks of the evolved B clonal lineages in tree format, their legend and the number of randomly considered nodes per network for Repertoire of 2 (Output of AntibodyForest).
#' E.g list2[[1][[1]] is an igraph object, containing the first network of the evolved B clonal lineage in tree format.
#' list2[[1]][[2]] contains the legend parameters of the first network of the evolved B clonal lineage.
#' list2[[1]][[3]] is the number of randomly considered nodes for the first network of the evolved B clonal lineage.
#' @param DAG  a logical variable, when TRUE a directed acyclic graph is produced.
#' @param clonal_frequency  a logical variable, TRUE if labeling of vertices is based on clonal frequency and FALSE otherwise.
#' @param scaleByClonalFreq  logical variable with TRUE if vertex size is scaled by the number of unique sequences per clone and FALSE otherwise.
#' @param weight  logical variable. When its value is FALSE, then the weights of outgoing edges from Germline node are set to 1. When its value is TRUE, the weights are set to the difference between the number of mutations among sequences in germline and connected nodes(value in the corresponding distance matrix) and the absolute value of the difference between the sequence lengths of germline and corresponding connected nodes. In both cases, weights of remaining edges are extracted from the distance matrix.
#' Outgoing edges from Germline represent the number of mutations of sequences having as common ancestor the Germline.
#' @param tie_flag  a string, with options 'rand', 'full', 'close_to_germ', 'far_from_germ', 'close_path_to_germ', 'far_path_from_germ','most_expanded' and 'least_expanded' for removing edges when equal distance (tie) in distance matrix.
#' 'rand' means random pruning in one of nodes, 'full' means keeping all nodes, close_to_germ means pruning of node(s) farthest from germline (based on number of intermediate nodes), 'far_from_germ' means pruning of node(s) closest to germline (based on number of intermediate nodes),
#' 'close_path_to_germ' means pruning of node(s) farthest from germline (based on edge path length), 'far_path_from_germ' meams pruning of node(s) closest to germline (based on edge path length),'most_expanded' means pruning of node(s) with the lowest B cell count(clonal frequency) and least_expanded, which means pruning of node(s) with the hightest B cell count(clonal frequency). In cases of subsequent ties, a random node is selected.
#' @param opt a string with options "isotype" and "cluster". The option "isotype" is utilized when the user desires to do an isotype analysis, while the selection of "cluster" denotes that an analysis based on transcriptome is requested.
#' @return  combined_df. A data.frame that summarizes metrics for both repertoires. In particular, each row represents a single network and networks of both repertoires are combined row wise.
#' Columns of combined_df are:
#' Column1: Weighted.Longest.path.from.germline.
#' Column2: Length.of.weighted.longest.shortest.path.from.germline.
#' Column3: Unweighted.Longest.path.from.germline.
#' Column4: Length.of.unweighted.longest.shortest.path.from.germline.
#' Column5: Average.number.of.daughter.cells.
#' Column6: Std.number.of.daughter.cells.
#' Column7: Min.number.of.daughter.cells.
#' Column8: Max.number.of.daughter.cells.
#' Column9: Weighted.vertex.degree.
#' Column10: Average.number.of.clusters/isotypes.
#' Column11: Isotypes/Clusters.info.
#' Column12: vertex.betweenness.centrality.
#' Column13: edge.betweenness.centrality.
#' Column14: closeness.centrality.of.vertices.
#' Column15: global.clustering.coefficient.
#' Column16: average.clustering.coefficient.
#' Column17: Mean.clonal.expansion.
#' If the labeling or scaling of nodes in graph is based on clonal frequency (arguments: clonal_frequency==TRUE or scaleByClonalFreq==TRUE), then combined_df contains also:
#' Column18: Ratio.Number.of.edges.from.germline.to.each.node.with.clonal.frequency.
#' Column19: Mean.Ratio.Number.of.edges.from.germline.to.each.node.with.clonal.frequency.
#' Column20: Mean.number.of.edges.from.germline.
#' Column21: Ratio.Total.path.length.from.germline.to.each.node.with.clonal.frequency.
#' Column22: Mean.Ratio.Total.path.length.from.germline.to.each.node.with.clonal.frequency.
#' Column23: Mean.Total.path.length.from.germline.
#' Column24: Repertoire.id.
#' Column25: Number.of.sequences.
#' @return  isotype_info_rep1  A data.frame. It summarizes isotype/cluster info for repertoire 1.
#' @return  isotype_info_rep2  A data.frame. It summarizes isotype/cluster info for repertoire 2.
#' @export
#' @seealso AntibodyForest, ForestMetrics
#' @examples
#' \dontrun{
#' CompareForests(list1,list2,DAG=TRUE,
#' clonal_frequency=TRUE,scaleByClonalFreq=TRUE,weight=TRUE,
#' tie_flag='close_to_germ',opt="cluster")
#'}

AbForests_CompareForests<-function(list1,list2,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,opt){

  metrics_list1<-AbForests_ForestMetrics(list1,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,opt)
  metrics_list2<-AbForests_ForestMetrics(list2,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,opt)
  df1<-.CREATE_DF(list1,metrics_list1,1)
  df2<-.CREATE_DF(list2,metrics_list2,2)
  df_total<-rbind(df1, df2)

  iso_info_list1<-sapply(metrics_list1,function(x) x[16])
  res_list1 <- lapply(iso_info_list1, function(x) x %>% dplyr::group_by(x$Parent,x$Child) %>% dplyr::slice(which.max(x$Freq)))
  merged.data.frame_list1 <- Reduce(function(...) merge(..., all=T), res_list1)
  merged.data.frame_list1<-unique(merged.data.frame_list1[ , -which(names(merged.data.frame_list1) %in% c('Freq'))])
  row.names(merged.data.frame_list1) <- NULL
  merged.data.frame_list1<-merged.data.frame_list1[1:2]
  iso_info_lis2<-sapply(metrics_list2,function(x) x[16])
  res_list2 <- lapply(iso_info_lis2, function(x) x %>% dplyr::group_by(x$Parent, x$Child) %>% dplyr::slice(which.max(x$Freq)))
  merged.data.frame_list2 <- Reduce(function(...) merge(..., all=T), res_list2)
  merged.data.frame_list2<-unique(merged.data.frame_list2[ , -which(names(merged.data.frame_list2) %in% c('Freq'))])
  row.names(merged.data.frame_list2) <- NULL
  merged.data.frame_list2<-merged.data.frame_list2[1:2]


  return(list(combined_df=df_total,isotype_info_rep1=merged.data.frame_list1,isotype_info_rep2=merged.data.frame_list2))
}

