#' Infer and draw B cell evolutionary networks

#' @description AntibodyForest takes the output of either ConvertStructure or CsvToDf or SubRepertoires or RemoveNets and outputs B cell phylogenetic networks in tree format. There is also the possibility to give the full-length list of clonal lineages, which contains both isotype and transcriptional cluster information, only when no prior data transformation is desired. Each network represents a clonal lineage, referring to the number of B cell receptor sequences originating from an independent V(D)J recombination event. Each vertex represents a unique recovered full-length variable heavy and light chain antibody sequence of a clonal family.
#' Edges separating nodes are drawn given that clonal variants are similarly related according to their Levenshtein distance. Edge weights are extacted from the distance matrix apart from the special case of unmutated germline, in which the weights of outgoing edges from it are either set to 1 or to the difference between the corresponding distance from the matrix and the absolute value of the difference between the sequence lengths of germline and corresponding connected nodes. At tree building, starting from the reference ancestral germline, each node is connected to nodes that can be reached via the minimum distance based on the distance matrix calculation. Therefore, potential edges that go back to previous tree layers along with bidirectional circles are eliminated. Polytomies, displayed by B cell clones producing multiple distinct offsprings, are resolved in case of reaching nodes with equal minimum distance. Indeed, the algorithm removes edges either randomly from the recipient nodes,based on the node closest or farthest from the germline, considering the number of intermediate nodes or edge path length, or the highest/lowest counting of cells on the present node. Additional ties are settled by random edge selection. Consequently, parsimony holds, meaning that each daughter node has only one parent.
#' Distinct tree topologies enable to visually investigate the trade-off between balance and evolution, and further quantify the amount of diversification of the subsequent detected clonal abundant clones during somatic hypermutation and class switching. The minimum decision-based criterion determines the amount of balance presented in the tree, while the maximum decision-based method the amount of evolution presented in the tree.
#' Single color or color distribution on each node demonstrates the proportion of B cells with the specific isotype(s) or transcriptional cluster(s), while setting the size of vertices can be performed based on the number of unique sequences per clone, vertex betweenness and vertex closeness. Scaling of nodes by their relative clonal expansion assists in pinpointing identical antibody sequences across a multitude of B cells. Node labeling can depict clonal frequency.
#' @param full_list a list of clone lineages, represented as data.frames
#' @param csv an indicator variable. TRUE if full_list argument is a list of csv files, FALSE otherwise
#' @param files  a list of data.frames. Each data.frame contains 2 columns, one that describes the sequences and the other which type of information (isotype or cluster) is considered in the analysis. All these cases are determined by the user.
#' @param distance_mat  a custom integer distance matrix, or NULL for using the default distance matrix (calucated based on the levenshtein distance, which counts the number of mutations between sequences).
#' @param clonal_frequency  a logical variable, TRUE if labeling of vertices is based on clonal frequency and FALSE otherwise.
#' @param scaleByClonalFreq  logical variable with TRUE if vertex size is scaled by the number of unique sequences per clone and FALSE otherwise.
#' @param weight  logical variable. When its value is FALSE, then the weights of outgoing edges from Germline node are set to 1. When its value is TRUE, the weights are set to the difference between the number of mutations among sequences in germline and connected nodes(value in the corresponding distance matrix) and the absolute value of the difference between the sequence lengths of germline and corresponding connected nodes. In both cases, weights of remaining edges are extracted from the distance matrix.
#' Outgoing edges from Germline represent the number of mutations of sequences having as common ancestor the Germline.
#' @param tie_flag  a string, with options 'rand', 'full', 'close_to_germ', 'far_from_germ', 'close_path_to_germ', 'far_path_from_germ','most_expanded' and 'least_expanded' for removing edges when equal distance (tie) in distance matrix.
#' 'rand' means random pruning in one of nodes, 'full' means keeping all nodes, close_to_germ means pruning of node(s) farthest from germline (based on number of intermediate nodes), 'far_from_germ' means pruning of node(s) closest to germline (based on number of intermediate nodes),
#' 'close_path_to_germ' means pruning of node(s) farthest from germline (based on edge path length), 'far_path_from_germ' meams pruning of node(s) closest to germline (based on edge path length),'most_expanded' means pruning of node(s) with the lowest B cell count(clonal frequency) and least_expanded, which means pruning of node(s) with the hightest B cell count(clonal frequency). In cases of subsequent ties, a random node is selected.
#' @param scaleBybetweenness  logical variable with TRUE if vertex size is scaled by the vertex betweenness centrality.
#' @param scaleByclocloseness_metr  logical variable with TRUE if vertex size is scaled by closeness centrality of vertices in graph.
#' @param opt a string with options "isotype" and "cluster". The option "isotype" is utilized when the user desires to do an isotype analysis, while the selection of "cluster" denotes that an analysis based on transcriptome is requested.
#' @param random.seed a random seed, specified by the user, when random sampling of sequences happens in each of the cases described in tie_flag argument.
#' @param alg_opt a string denoting the version of the edge selection algorithm used in the construction of networks. Possible choices: "naive", "two-step".
#' @param cdr3 variable with values 0 if the user desires to select full length sequences (only when the input is a list of csv files), 1 for sequences in the CDR3 only (only when the input is a list of csv files) and NULL otherwise.
#' @return  graphs. A list of lists.
#' E.g graphs[[1][[1]] network: an igraph object, containing the first network in tree format.
#' graphs[[1]][[2]] legend: contains the legend parameters of the first network.
#' graphs[[1]][[3]] count.rand: contains the number of randomly considered nodes for the first network.
#' graphs[[1]][[4]] adj.matrix: contains the adjacency matrix for the first network.
#' graphs[[1]][[5]] distance.matrix: contains the distance matrix for the first network.
#' graphs[[1]][[6]] cells.per.network: contains the number of cells for the first network.
#' graphs[[1]][[7]] variants.per.network: contains the number of variants for the first network.
#' graphs[[1]][[8]] variant.sequences: contains the sequences of the variants for the first network.
#' graphs[[1]][[9]] cells.per.variant: contains the number of cells per variant (clonal frequency) for the first network.
#' graphs[[1]][[10]] cell.indicies.per.variant:  the indices of cells per variant for the first network.
#' graphs[[1]][[11]] new.variant.names: contains the names of variants for the first network.
#' graphs[[1]][[12]] germline.index: contains the index of germline sequence for the first network.
#' graphs[[1]][[13]] isotype.per.variant: contains the isotypes corresponding to each variant for the first network.
#' graphs[[1]][[14]] transcriptome.cluster.per.variant: contains the transcriptional clusters corresponding to each variant for the first network.
#' graphs[[1]][[15]] isotype.per.cell: contains the isotype corresponding to each cell for the first network.
#' graphs[[1]][[16]] transcriptome.cluster.per.cell: contains the transcriptional cluster corresponding to each cell for the first network.
#' @export
#' @seealso ConvertStructure, CsvToDf, SubRepertoires, RemoveNets
#' @examples
#' \dontrun{
#' AntibodyForest(full_list,csv=FALSE, files,custom_mat=NULL,clonal_frequency=TRUE,
#' scaleByClonalFreq=TRUE,weight=TRUE,tie_flag='close_to_germ',
#' scaleBybetweenness=FALSE,scaleByclocloseness_metr=FALSE,
#' opt="cluster",seed=165,alg_opt="0",cdr3=NULL)
#'}

AbForests_AntibodyForest<-function(full_list,csv,files,distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleBybetweenness,scaleByclocloseness_metr,opt,random.seed,alg_opt,cdr3){

  graphs<-c()
  i<-0
  if(!is.null(full_list)){
    if (csv==TRUE){
      templ<-AbForests_CsvToDf(full_list)
      df_full<-AbForests_ConvertStructure(templ,opt,cdr3)
      if(opt=="isotype" | opt=="cluster"){
        if(opt=="cluster"){
          no_cluster<-lapply(df_full[[1]], function(z) unique(z$cluster))#length(unique(z[["cluster"]])))
          number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
        }

        isotype.per.cell<-sapply(AbForests_ConvertStructure(templ,"isotype",cdr3),sapply,function(x) x[2])
        transcriptome.cluster.per.cell<-sapply(AbForests_ConvertStructure(templ,"cluster",cdr3),sapply,function(x) x[2])
      }else{
        no_cluster<-lapply(df_full[[1]], function(z) unique(z[[opt]]))#length(unique(z[["cluster"]])))
        number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
      }

      if(.HAS_EMPTY_LIST(df_full)){
        print('No elements exist in list')
        return()
      }
      if(length(df_full)==0){
        return()
      }
      if(opt=="isotype"){
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x$isotype,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"isotype",random.seed,alg_opt,NULL)),recursive = FALSE)
      }else if(opt=="cluster"){
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x$cluster,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"cluster",random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
      }else{
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x[[opt]],parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,opt,random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
      }
      if(opt=="isotype" | opt=="cluster"){
        if(!(length(graphs)==length(isotype.per.cell))){
          isotype.per.cell<-lapply(isotype.per.cell,function(z) unlist(z,recursive = FALSE))
          transcriptome.cluster.per.cell<-lapply(transcriptome.cluster.per.cell,function(z) unlist(z,recursive = FALSE))
        }
        if(length(isotype.per.cell)==length(graphs)){
          for (i in 1:length(isotype.per.cell)){
            graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=isotype.per.cell[[i]]))
            graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=transcriptome.cluster.per.cell[[i]]))
          }
        }
      }
    }else{
      i<-i+1
      df_full<-AbForests_ConvertStructure(full_list,opt,NULL)
      if(opt=="isotype" | opt=="cluster"){
        if(length(df_full)>1){
          isotype.per.cell<-list(sapply(AbForests_ConvertStructure(full_list,"isotype",NULL),sapply, function(x) x[[2]]))
          transcriptome.cluster.per.cell<-list(sapply(AbForests_ConvertStructure(full_list,"cluster",NULL),sapply,function(x) x[[2]]))
        }else{
          isotype.per.cell<-sapply(AbForests_ConvertStructure(full_list,"isotype",NULL),sapply, function(x) x[[2]])
          transcriptome.cluster.per.cell<-sapply(AbForests_ConvertStructure(full_list,"cluster",NULL),sapply,function(x) x[[2]])
        }
        if(opt=="cluster"){
          no_cluster<-lapply(df_full[[1]], function(z) unique(z$cluster))
          number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))

        }
      }else{
        no_cluster<-lapply(df_full[[1]], function(z) unique(z[[opt]]))
        number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
      }
      if(opt=="isotype"){
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x$isotype,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"isotype",random.seed,alg_opt,NULL)),recursive = FALSE)
      }else if (opt=="cluster"){
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x$cluster,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"cluster",random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
      }else{
        graphs<-unlist(lapply(df_full, lapply,function(x) .MAIN(x$Seq,x[[opt]],parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,opt,random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
      }
      if(opt=="isotype" | opt=="cluster"){
        if(!(length(graphs)==length(isotype.per.cell))){
          isotype.per.cell<-lapply(isotype.per.cell,function(z) unlist(z,recursive = FALSE))
          transcriptome.cluster.per.cell<-lapply(transcriptome.cluster.per.cell,function(z) unlist(z,recursive = FALSE))
        }
        if(length(isotype.per.cell)==length(graphs)){
          for (i in 1:length(isotype.per.cell)){
            graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=isotype.per.cell[[i]]))
            graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=transcriptome.cluster.per.cell[[i]]))
          }
        }else{
          for (i in 1:lengths(isotype.per.cell)){
            graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=as.vector(sapply(isotype.per.cell,`[[`, i))))
            graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=as.vector(sapply(transcriptome.cluster.per.cell,`[[`, i))))
          }
        }

      }

      }

  }else{

    if(.HAS_EMPTY_LIST(files)){
      print('No elements exist in list')
      return()
    }
    if(length(files)==0){
      return()
    }
    if(purrr::vec_depth(files)>3){
      files<-do.call(c, files)
    }
    if (any(sapply(files, sapply,is.list))) {
      list_clust_ext<-c()
      if(opt=="isotype" | opt=="cluster"){
        if(opt=="isotype"){
          isotype.per.cell<-list(sapply(files,sapply,function(x) x[[2]]))
          size<-sapply(files,sapply,function(x)lengths(x)[2])
          list_clust<-c("Unknown")
          transcriptome.cluster.per.cell<-mapply(rep, list_clust, size)
          list_clust_ext[[1]]<-transcriptome.cluster.per.cell
          transcriptome.cluster.per.cell<-list_clust_ext
          graphs<-unlist(lapply(files, lapply,function(x) .MAIN(x$Seq,x$isotype,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"isotype",random.seed,alg_opt,NULL)),recursive = FALSE)
        }else{
          transcriptome.cluster.per.cell<-list(sapply(files,sapply,function(x) x[[2]]))
          size<-sapply(files,sapply,function(x)lengths(x)[2])
          list_clust<-c("Unknown")
          isotype.per.cell<-mapply(rep, list_clust, size)
          list_clust_ext[[1]]<-isotype.per.cell
          isotype.per.cell<-list_clust_ext
          if(opt=="cluster"){
            no_cluster<-lapply(df_full[[1]], function(z) unique(z$cluster))#length(unique(z[["cluster"]])))
            number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
          }
          graphs<-unlist(lapply(files, lapply,function(x) .MAIN(x$Seq,x$cluster,parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"cluster",random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
        }
        if(length(isotype.per.cell)==length(graphs)){
          for (i in 1:length(isotype.per.cell)){
            graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=isotype.per.cell[[i]]))
            graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=transcriptome.cluster.per.cell[[i]]))
          }
        }else{
          for (i in 1:lengths(isotype.per.cell)){
            graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=as.vector(sapply(isotype.per.cell,`[[`, i))))
            graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=as.vector(sapply(transcriptome.cluster.per.cell,`[[`, i))))
          }
        }
      }else{
        no_cluster<-lapply(files[[1]], function(z) unique(z[[opt]]))
        number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
        graphs<-unlist(lapply(files, lapply,function(x) .MAIN(x$Seq,x[[opt]],parent.frame()$i[],distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,opt,random.seed,alg_opt,number.of.clusters)),recursive = FALSE)
      }
    }else if (any(sapply(files, class) == "list") ) {
      for (file in files){
        i<-i+1
        if(length(file)==0){
          return()
        }
        list_clust_ext<-c()
        if(opt=="isotype" | opt=="cluster"){
          if(opt=="isotype"){
            if(length(files)>1){
              isotype.per.cell<-list(sapply(files,function(x) x[[2]]))
              size<-sapply(files,lengths)[2,]
              list_clust<-c("Unknown")
              transcriptome.cluster.per.cell<-mapply(rep, list_clust, size)
              list_clust_ext[[1]]<-transcriptome.cluster.per.cell
              transcriptome.cluster.per.cell<-list_clust_ext

            }else{
              isotype.per.cell<-list(list(sapply(files,function(x) x[[2]])))
              size<-sapply(files,lengths)[2,]
              list_clust<-c("Unknown")
              transcriptome.cluster.per.cell<-mapply(rep, list_clust, size)
              list_clust_ext[[1]]<-list(transcriptome.cluster.per.cell)
              transcriptome.cluster.per.cell<-list_clust_ext
            }
            list_net_graph<-.MAIN(file$Seq,file$isotype,i,distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"isotype",random.seed,alg_opt,NULL)
          }else{
            if(length(files)>1){
              transcriptome.cluster.per.cell<-list(sapply(files,function(x) x[[2]]))
              size<-sapply(files,lengths)[2,]
              list_clust<-c("Unknown")
              isotype.per.cell<-mapply(rep, list_clust, size)
              list_clust_ext[[1]]<-isotype.per.cell
              isotype.per.cell<-list_clust_ext

            }else{
              transcriptome.cluster.per.cell<-list(list(sapply(files,function(x) x[[2]])))
              size<-sapply(files,lengths)[2,]
              list_clust<-c("Unknown")
              isotype.per.cell<-mapply(rep, list_clust, size)
              list_clust_ext[[1]]<-list(isotype.per.cell)
              isotype.per.cell<-list_clust_ext
            }
            if(opt=="cluster"){
              no_cluster<-lapply(files, function(z) unique(z$cluster))#length(unique(z[["cluster"]])))
              number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
            }
            list_net_graph<-.MAIN(file$Seq,file$cluster,i,distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,"cluster",random.seed,alg_opt,number.of.clusters)
          }
        }else{
          no_cluster<-lapply(files, function(z) unique(z[[opt]]))
          number.of.clusters<-unique(unlist(lapply(no_cluster, function(z)  subset(z, !grepl(paste(c("germline","Unknown"), collapse= "|"), z)))))
          list_net_graph<-.MAIN(file$Seq,file[[opt]],i,distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,opt,random.seed,alg_opt,number.of.clusters)
        }

        graphs[[i]]<-list_net_graph
      }
      if(opt=="isotype" | opt=="cluster"){
        for (i in 1:length(files)){
          graphs[[i]]<-append(graphs[[i]],list("isotype.per.cell"=as.vector(sapply(isotype.per.cell,`[[`, i))))
          graphs[[i]]<-append(graphs[[i]],list("transcriptome.cluster.per.cell"=as.vector(sapply(transcriptome.cluster.per.cell,`[[`, i))))
        }
      }
    }
  }
  return(graphs)
}
