#' Count the number of unique antibody variants per clonal lineage

#' @description UniqueAntibodyVariants calculates the number of unique antibody sequences, as dictated by the different grouping sequences strategy,for each network in the immune repertoire.
#' @param list a list of data.frames. Each data.frame represents a clone lineage and contains information on the antibody sequences and on the isotype/transcriptional cluster is considered in the analysis based the user's preferences.
#' @param opt a string with options "isotype" and "cluster". The option "isotype" is utilized when the user desires to do an isotype analysis, while the selection of "cluster" denotes that an analysis based on transcriptome is requested.
#' @param distance_mat a custom integer distance matrix, or NULL for using the default distance matrix (calucated based on the levenshtein distance, which counts the number of mutations between sequences). Given the phylogenetic tree, a custom-made distance matrix can be produced by PlyloToMatrix function.
#' @param tie_flag  a string, with options 'rand', 'full', 'close_to_germ', 'far_from_germ', 'close_path_to_germ', 'far_path_from_germ','most_expanded' and 'least_expanded' for removing edges when equal distance (tie) in distance matrix.
#' 'rand' means random pruning in one of nodes, 'full' means keeping all nodes, close_to_germ means pruning of node(s) farthest from germline (based on number of intermediate nodes), 'far_from_germ' means pruning of node(s) closest to germline (based on number of intermediate nodes),
#' 'close_path_to_germ' means pruning of node(s) farthest from germline (based on edge path length), 'far_path_from_germ' meams pruning of node(s) closest to germline (based on edge path length),'most_expanded' means pruning of node(s) with the lowest B cell count(clonal frequency) and least_expanded, which means pruning of node(s) with the hightest B cell count(clonal frequency). In cases of subsequent ties, a random node is selected.
#' @param weight  logical variable. When its value is FALSE, then the weights of outgoing edges from Germline node are set to 1. When its value is TRUE, the weights are set to the difference between the number of mutations among sequences in germline and connected nodes(value in the corresponding distance matrix) and the absolute value of the difference between the sequence lengths of germline and corresponding connected nodes. In both cases, weights of remaining edges are extracted from the distance matrix.
#' @param random.seed a random seed, specified by the user, when random sampling of sequences happens in each of the cases described in tie_flag argument.
#' @param alg_opt a string denoting the version of the edge selection algorithm used in the construction of networks. "0" means the naive version and "1" the advanced one.
#' @param cdr3 variable with values 0 if the user desires to select full length sequences (only when the input is a list of csv files), 1 for sequences in the CDR3 only (only when the input is a list of csv files) and NULL otherwise.
#' @return  uni_seq a vector, same size as list, which contains the number of unique antibody variants for each clonal lineage.
#' @export
#' @examples
#' \dontrun{
#' UniqueAntibodyVariants(list,opt="cluster",
#' distance_mat=NULL,tie_flag=close_to_germ,weight=TRUE,random.seed=165,alg_opt="naive",cdr3=NULL)
#'}


AbForests_UniqueAntibodyVariants<-function(list,opt,distance_mat,tie_flag,weight,random.seed,alg_opt,cdr3){

  if(length(cdr3)>0){
    list<-AbForests_ConvertStructure(list,opt,cdr3)
  }else{
    list<-AbForests_ConvertStructure(list,opt,NULL)
  }
  if (any(sapply(list, sapply,is.list))) {

    uni_seq<-unname(unlist(lapply(list,lapply, function(k) .COUNT_ELEMENTS(k,opt,distance_mat,tie_flag,weight,random.seed,alg_opt))))

  }else{
    uni_seq<-unname(unlist(lapply(list, function(k) .COUNT_ELEMENTS(k,opt,distance_mat,tie_flag,weight,random.seed,alg_opt))))

  }
  return(uni_seq)
}

