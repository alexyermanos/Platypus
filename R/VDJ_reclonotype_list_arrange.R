#' Only Platypus v2 Organizes the top N genes that define each Seurat cluster and converts them into a single dataframe. This can be useful for obtaining insight into cluster-specific phenotypes.
#' @param VDJ_clonotype.output The output object from the VDJ_clonotype function. The column of the merged nucleotide clonotype IDs will be used to rearrange the new object.
#' @param VDJ_analyze.output The output from the initial VDJ_analyze, containing clonotype information based on nucleotide sequence.
#' @param Platypus_list.object The new list object from one of Platypus functions (for example, clonal lineages, VDJ_per_clne, etc) that should be merged based on the VDJ_clonotype output structure. nested list structure, where outer list corresponds to repertoire and the inner list corresponds to clones (on the nucleotide level).
#' @return Returns a dataframe in which the top N genes defining each cluster based on differential expression are selected.
#' @export
#' @examples
#' \dontrun{
#' checking_vdj_reclono <- VDJ_reclonotype_list_arrange(
#' VDJ_clonotype.output = repertoire_reclonotype
#' ,VDJ_analyze.output = repertoire_list
#' ,Platypus_list.object = repertoire_vdj_per_clone)
#'}
VDJ_reclonotype_list_arrange <- function(VDJ_clonotype.output,
                                         VDJ_analyze.output,
                                         Platypus_list.object){
  output_list <- list()
  for(i in 1:length(VDJ_clonotype.output)){
    output_list[[i]] <- list()
    for(j in 1:nrow(VDJ_clonotype.output[[i]])){
      ## string split the nt clonotype ids
      new_clonotype_nt_ids <- stringr::str_split(VDJ_clonotype.output[[i]]$nt_clone_ids[j],pattern = ";")[[1]]
      dfs_to_merge <- match(new_clonotype_nt_ids,VDJ_analyze.output[[i]]$clonotype_id)
      output_list[[i]][[j]] <- do.call("rbind",Platypus_list.object[[i]][dfs_to_merge])

    }
  }
  return(output_list)
}
