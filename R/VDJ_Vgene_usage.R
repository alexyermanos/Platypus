#' Produces a matrix counting the number of occurencies for each IgH and IgK/L Vgene combinations.
#' @param VDJ.clonotype.output Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @return Returns a list of matrices containing the number of Vgene heavy/light chain combinations per repertoire.
#' @seealso
#' @examples
#' \dontrun{
#' example.vdj.vgene_usage <- VDJ_Vgene_usage(VDJ.clonotype.output = VDJ.clonotype.output)
#'}
VDJ_Vgene_usage <- function(VDJ.clonotype.output){

    Vgene_usage_matrix <- list()
    dummy_df <- list()

    for (k in 1:length(VDJ.clonotype.output)){

      #make a new column with both IgH and IgL/K V genes
      VDJ.clonotype.output[[k]]$vgenes <- paste(VDJ.clonotype.output[[k]]$HC_vgene, VDJ.clonotype.output[[k]]$LC_vgene, sep = "_")

      #Create a matrix with rows being heavy chain and columns being light chain v genes
      Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(VDJ.clonotype.output[[k]]$HC_vgene)), ncol = length(unique(VDJ.clonotype.output[[k]]$LC_vgene)))

      #give the row names and col names
      rownames(Vgene_usage_matrix[[k]]) <- unique(VDJ.clonotype.output[[k]]$HC_vgene)
      colnames(Vgene_usage_matrix[[k]]) <- unique(VDJ.clonotype.output[[k]]$LC_vgene)

      #create dummy df which will contain the counts for each combination
      dummy_df[[k]] <- as.data.frame(table(VDJ.clonotype.output[[k]]$vgenes))
      colnames(dummy_df[[k]]) <- c("vgene", "count")

      #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
      for (i in 1:nrow(Vgene_usage_matrix[[k]])){
        for (j in 1:ncol(Vgene_usage_matrix[[k]])){

          if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_df[[k]]$vgene){
            Vgene_usage_matrix[[k]][i,j] <- dummy_df[[k]][which(dummy_df[[k]]$vgene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
          } else {
            Vgene_usage_matrix[[k]][i,j] <- 0
          }

        }
      }
    }
    return(Vgene_usage_matrix)
}
