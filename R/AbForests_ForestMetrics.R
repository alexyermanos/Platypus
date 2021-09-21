#' Calculate metrics for networks

#' @description ForestMetrics takes the output of AntibodyForest and calculates metrics for each of the networks.
#' @param graphs A list of lists. Each sublist contains an igraph object with the networks of the evolved B clonal lineages in tree format, their legend and the number of randomly considered nodes per network(Output of AntibodyForest function).
#' E.g graphs[[1][[1]] is an igraph object, containing the first network of the evolved B clonal lineage in tree format.
#' graphs[[1]][[2]] contains the legend parameters of the first network of the evolved B clonal lineage.
#' graphs[[1]][[3]] is the number of randomly considered nodes for the first network of the evolved B clonal lineage.
#' @param DAG  a logical variable, when TRUE a directed acyclic graph is produced.
#' @param clonal_frequency  a logical variable, TRUE if labeling of vertices is based on clonal frequency and FALSE otherwise.
#' @param scaleByClonalFreq  logical variable with TRUE if vertex size is scaled by the number of unique sequences per clone and FALSE otherwise.
#' @param weight  logical variable. When its value is FALSE, then the weights of outgoing edges from Germline node are set to 1. When its value is TRUE, the weights are set to the difference between the number of mutations among sequences in germline and connected nodes(value in the corresponding distance matrix) and the absolute value of the difference between the sequence lengths of germline and corresponding connected nodes. In both cases, weights of remaining edges are extracted from the distance matrix.
#' Outgoing edges from Germline represent the number of mutations of sequences having as common ancestor the Germline.
#' @param tie_flag  a string, with options 'rand', 'full', 'close_to_germ', 'far_from_germ', 'close_path_to_germ', 'far_path_from_germ','most_expanded' and 'least_expanded' for removing edges when equal distance (tie) in distance matrix.
#' 'rand' means random pruning in one of nodes, 'full' means keeping all nodes, close_to_germ means pruning of node(s) farthest from germline (based on number of intermediate nodes), 'far_from_germ' means pruning of node(s) closest to germline (based on number of intermediate nodes),
#' 'close_path_to_germ' means pruning of node(s) farthest from germline (based on edge path length), 'far_path_from_germ' meams pruning of node(s) closest to germline (based on edge path length),'most_expanded' means pruning of node(s) with the lowest B cell count(clonal frequency) and least_expanded, which means pruning of node(s) with the hightest B cell count(clonal frequency). In cases of subsequent ties, a random node is selected.
#' @param opt a string with options "isotype" and "cluster". The option "isotype" is utilized when the user desires to do an isotype analysis, while the selection of "cluster" denotes that an analysis based on transcriptome is requested.
#' @return  metrics.  A list of lists. Each list contains various metrics for the quantification of networks.
#' E.g metrics[[1][[1]] is the weighted Longest path from germline for the first network.
#' metrics[[1]][[2]] is the length of weighted longest shortest path from germline for the first network.
#' metrics[[1]][[3]] is the unweighted Longest path from germline for the first network.
#' metrics[[1]][[4]] is the length of unweighted longest shortest path from germline for the first network.
#' metrics[[1]][[5]] is the weighted shortest path network for the first network.
#' metrics[[1]][[6]] is the uweighted shortest path network for the first network.
#' metrics[[1]][[7]] is the average number of daughter cells for the first network.
#' metrics[[1]][[8]] is the std number of daughter cells for the first network.
#' metrics[[1]][[9]] is the min number of daughter cells for the first network.
#' metrics[[1]][[10]] is the max number of daughter cells for the first network.
#' metrics[[1]][[11]] is a ggplot object that contains the plot of Degree Distribution of daughter cells for the first network.
#' metrics[[1]][[12]] is the weighted vertex degree for the first network.
#' metrics[[1]][[13]] is a ggplot object that contains the plot of unweighted Degree Distribution of daughter cells for the first network.
#' metrics[[1]][[14]] is the average number of isotypes for the first network.
#' metrics[[1]][[15]] is a ggplot object that contains the plot of Distribution of isotypes for the first network.
#' metrics[[1]][[16]] is the Isotypes/Clusters info data.frame with columns Parent, Child and Parent-Child, which contains the type of isotypes/clusters for each pair of nodes (Parent-Child relationship in tree) found in the first network.
#' metrics[[1]][[17]] is a ggplot object that contains the plot of Isotype/Cluster Directionality for the first network. In particular, the frequency of all types of isotypes/clusters for each pair of nodes in the tree is depicted.
#' metrics[[1]][[18]] is the vertex betweenness centrality for the first network. It is defined by the number of geodesics (shortest paths) going through a vertex according to igraph documentation.
#' metrics[[1]][[19]] is the edge betweenness centrality for the first network. It is defined by the number of geodesics (shortest paths) going through an edge according to igraph documentation.
#' metrics[[1]][[20]] is the closeness centrality of vertices for the first network. The closeness centrality of a vertex is defined by the inverse of the average length of the shortest paths to/from all the other vertices in the graph according to igraph documentation.
#' metrics[[1]][[21]] is a ggplot object that contains the plot of Path length from Germline vs Node Degree for the first network.
#' metrics[[1]][[22]] is a ggplot object that contains the plot of Number of edges from Germline vs Node Degree for the first network.
#' metrics[[1]][[23]] is an igraph object that contains the Isotype/Cluster transition network for the first network.
#' metrics[[1]][[24]] is the global clustering coefficient for the first network.
#' metrics[[1]][[25]] is the average clustering coefficient for the first network.
#' metrics[[1]][[26]] is the mean clonal expansion for the first network, calculated as the mean of clonal frequencies of all vertices in the network.
#' If the labeling or scaling of nodes in graph is based on clonal frequency (arguments: clonal_frequency==TRUE or scaleByClonalFreq==TRUE), then
#' metrics[[1]][[27]] is the ratio: Number of edges from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[28]] is the mean ratio: Number of edges from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[29]] is a ggplot object that contains the ratio of Number of edges from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[30]] is the mean number of edges from germline for the first network.
#' metrics[[1]][[31]] is the ratio: Total path length from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[32]] is the mean ratio: Total path length from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[33]] is a ggplot object that contains the ratio of Total path length from germline to each node with clonal frequency for the first network.
#' metrics[[1]][[34]] is the mean Total path length from germline for the first network.
#' metrics[[2]][[1]] is the weighted Longest path from germline for the second network.
#' @export
#' @seealso AntibodyForest
#' @examples
#' \dontrun{
#' ForestMetrics(graphs,DAG=TRUE,clonal_frequency=TRUE,scaleByClonalFreq=TRUE,
#' weight=TRUE,tie_flag='close_to_germ',opt="cluster")
#'}

AbForests_ForestMetrics<-function(graphs,DAG,clonal_frequency,scaleByClonalFreq,weight,tie_flag,opt){
  metrics<-list()
  lapply(graphs, function(k){
    tryCatch(
      expr = {
        if(DAG==TRUE & tie_flag!='full'){
          sol<-.LONGEST_PATH_FROM_GERMLINE(k[[1]],"weighted",igraph::E(k[[1]])$weight)
          list_g<-.COLOR_PATH(k[[1]],sol$shortest_paths,k[[1]]$all_colors,opt,k[[1]]$progr_flag)
          sol_un<-.LONGEST_PATH_FROM_GERMLINE(k[[1]],"unweighted",igraph::E(k[[1]])$weight)
          list_g_un<-.COLOR_PATH(k[[1]],sol_un$shortest_paths,k[[1]]$object_colors,opt,k[[1]]$progr_flag)
          list_avg_daughters<-.AVERAGE_DAUGHTERS(k[[1]],igraph::V(k$network)$order)
          list_strengths<-.GRAPH_STRENGTH_DISTRIBUTION(k[[1]],igraph::E(k[[1]])$weight)
          isotypes_avg<-.AVERAGE_ISOTYPES(k[[1]],igraph::V(k[[1]]),k[[1]]$object_colors,opt)
          isotypes<-.ISOTYPE_POSITION(k[[1]],igraph::V(k$network)$order,k$network$edges_list,k[[1]]$object_colors,opt)
          transitions_isotypes<-.ISOTYPE_TRANSITIONS(isotypes$combs)
          p_EdgesFromGerm<-.NODE_DEGREE_EDGES_FROM_GERM(k[[1]],igraph::E(k[[1]])$weight)
          p_path_from_germ<-.NODE_DEGREE_PATH_FROM_GERM(k[[1]],igraph::E(k[[1]])$weight)
         centrality_metrics<-.CENTRALITY_METRICS(k[[1]],igraph::V(k$network)$order,igraph::E(k[[1]]))
          clust_coef<-.CLUSTER_COEFFICIENTS(k[[1]])
          mean_clonal_expansion<-mean(unlist(igraph::V(k$network)$clone_freq))
          if (clonal_frequency==TRUE | scaleByClonalFreq==TRUE){
            ratio_EdgesFromGerm<-.CLONAL_FREQUENCY_EDGES_FROM_GERM(k[[1]],igraph::E(k[[1]])$weight)
            mean_ratio_EdgesFromGerm<-mean(ratio_EdgesFromGerm$ratio)
            ratio_PathFromGerm<-.CLONAL_FREQUENCY_PATH_FROM_GERM(k[[1]],igraph::E(k[[1]])$weight)
            mean_ratio_PathFromGerm<-mean(ratio_PathFromGerm$ratio)
            temp<-list( "Ratio Number of edges from germline to each node with clonal frequency"=ratio_EdgesFromGerm$ratio,"Mean Ratio Number of edges from germline to each node with clonal frequency"=mean_ratio_EdgesFromGerm,
                        "Plot : Number of edges from germline to each node/clonal frequency"=ratio_EdgesFromGerm$p,"Mean number of edges from germline"=ratio_EdgesFromGerm$mean_no_edges_from_germ,"Ratio Total path length from germline to each node with clonal frequency"=ratio_PathFromGerm$ratio,"Mean Ratio Total path length from germline to each node with clonal frequency "=mean_ratio_PathFromGerm,
                        "Plot: Total path length from germline to each node/clonal frequency"=ratio_PathFromGerm$p,"Mean Total path length from germline"=ratio_PathFromGerm$mean_score_idx)
          }
          if(opt=="isotype"){
            output<-list("Weighted Longest path from germline" = sol$shortest_paths, "Length of weighted longest shortest path from germline" =sol$score,"Unweighted Longest path from germline" = sol_un$shortest_paths,
                         "Length of unweighted longest shortest path from germline" =sol_un$score,"Weighted shortest path network"=list_g,"Uweighted shortest path network "=list_g_un,
                         "Average number of daughter cells"=list_avg_daughters$avg_daughters,"Std number of daughter cells"=list_avg_daughters$std_daughters,
                         "Min number of daughter cells"=list_avg_daughters$min_daughters,"Max number of daughter cells"=list_avg_daughters$max_daughters,"Plot: Degree Distribution of daughter cells"=list_avg_daughters$p,
                         "Weighted vertex degree"=list_strengths$edge_strengths,"Plot: Weighted Degree Distribution of daughter cells"=list_strengths$p,"Average number of isotypes"=isotypes_avg$avg_isotypes,"Plot: Distribution of isotypes"=isotypes_avg$p,"Isotypes info"=isotypes$combs,
                         "Plot : Isotype Directionality"=isotypes$p,"vertex betweenness centrality"= centrality_metrics$betweenessCentrality,"edge betweenness centrality"=centrality_metrics$edgeCentrality,
                         "closeness centrality of vertices"=centrality_metrics$clocloseness_metr,"Plot: Path length from Germline vs Node Degree"=p_path_from_germ,"Plot: Number of edges from Germline vs Node Degree"=p_EdgesFromGerm,"Isotype transition network"=transitions_isotypes,
                         "global clustering coefficient"=clust_coef$global_cluster_coef,"average clustering coefficient"=clust_coef$avg_cluster_coef,"Mean clonal expansion"=mean_clonal_expansion)
          }else{
            output<-list("Weighted Longest path from germline" = sol$shortest_paths, "Length of weighted longest shortest path from germline" =sol$score,"Unweighted Longest path from germline" = sol_un$shortest_paths,
                         "Length of unweighted longest shortest path from germline" =sol_un$score,"Weighted shortest path network"=list_g,"Uweighted shortest path network "=list_g_un,
                         "Average number of daughter cells"=list_avg_daughters$avg_daughters,"Std number of daughter cells"=list_avg_daughters$std_daughters,
                         "Min number of daughter cells"=list_avg_daughters$min_daughters,"Max number of daughter cells"=list_avg_daughters$max_daughters,"Plot: Degree Distribution of daughter cells"=list_avg_daughters$p,
                         "Weighted vertex degree"=list_strengths$edge_strengths,"Plot: Weighted Degree Distribution of daughter cells"=list_strengths$p,"Average number of clusters"=isotypes_avg$avg_isotypes,"Plot: Distribution of clusters"=isotypes_avg$p,"Clusters info"=isotypes$combs,
                         "Plot : Cluster Directionality"=isotypes$p,"vertex betweenness centrality"= centrality_metrics$betweenessCentrality,"edge betweenness centrality"=centrality_metrics$edgeCentrality,
                         "closeness centrality of vertices"=centrality_metrics$clocloseness_metr,"Plot: Path length from Germline vs Node Degree"=p_path_from_germ,"Plot: Number of edges from Germline vs Node Degree"=p_EdgesFromGerm,"Cluster transition network"=transitions_isotypes,
                         "global clustering coefficient"=clust_coef$global_cluster_coef,"average clustering coefficient"=clust_coef$avg_cluster_coef,"Mean clonal expansion"=mean_clonal_expansion)
          }

          if (clonal_frequency==TRUE | scaleByClonalFreq==TRUE){
            output<-append(output,temp)
          }
          metrics<-append(metrics,output)
          return(metrics)
        }else{
          return()
        }

      },
      error = function(e){
        message('Error: Metrics cannot be computed as graph is not acyclic!')
        return("Try a different option!")
      }
    )
  }
  )
}
