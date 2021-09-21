#' Produces a matrix indicating either the number of cells or clones which contain multiple heavy or light chains (or alpha/beta in the case of T cells).
#' @param clonotype.list Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @param clone.level Logical indicating whether the matrix should display information on the clone level. TRUE will result in matrices containing information about the number of chains on the clonal level. FALSE will result in matrices depicting the numnber of cells.
#' @return Returns a list of matrices containing the number of heavy/light chains per either cell or clone depending on the clone.level parameter. This can then be supplied to heatmap functions directly. Each list element corresponds to each of the input list elements of clonotypes.
#' @export
#' @examples
#' \dontrun{
#' example.vdj.analyze <- VDJ_dublets(clonotype.list = "VDJ.analyze.output", clone.level=T)
#'}
#'
VDJ_dublets <- function(clonotype.list,
                        clone.level){
  if(missing(clone.level)) clone.level <- F
  output.matrices <- list()
  for(i in 1:length(clonotype.list)){
    max.hc <- max(clonotype.list[[i]]$HC_count)
    max.lc <- max(clonotype.list[[i]]$LC_count)
    output.matrices[[i]] <- matrix(0,nrow=max.hc+1,ncol=max.lc+1)
    for(j in 1:nrow(output.matrices[[i]])){
      for(k in 1:ncol(output.matrices[[i]])){
        if(clone.level==T) output.matrices[[i]][j,k] <- which(clonotype.list[[i]]$HC_count==j-1 & clonotype.list[[i]]$LC_count==k-1)
        if(clone.level==F) output.matrices[[i]][j,k] <- sum(clonotype.list[[i]]$frequency[which(clonotype.list[[i]]$HC_count==j-1 & clonotype.list[[i]]$LC_count==k-1)])
      }
    }
  }
  return(output.matrices)
}
