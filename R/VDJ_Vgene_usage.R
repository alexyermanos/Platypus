#' Produces a matrix counting the number of occurences for each VDJ and VJ Vgene combinations for each list enty in VDJ.clonotype.output or for each sample_id in VDJ.matrix
#' @param VDJ.clonotype.output Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire. IF platypus.version == "v2"
#' @param VDJ.matrix Output VDJ dataframe from VDJ_GEX_matrix IF platypus.version == "v3"
#' @param platypus.version Character. Defaults to "v2". Can be "v2" or "v3" dependent on the input format 
#' @return Returns a list of matrices containing the number of Vgene heavy/light chain combinations per repertoire.
#' @seealso
#' @export
#' @examples
#' \dontrun{
#' example.vdj.vgene_usage <- VDJ_Vgene_usage(VDJ.clonotype.output = VDJ.clonotype.output)
#'}
VDJ_Vgene_usage <- function(VDJ.clonotype.output,
                            VDJ.matrix,
                            platypus.version){
  
    if(missing(platypus.version)) platypus.version <- "v2"

    if(platypus.version == "v2"){
      if(missing(VDJ.clonotype.output)) stop("When using platypus version v2 please provide an input list to VDJ.clonotype.output. If using a VDJ matrix from the function VDJ_GEX_matrix, please switch platypus.version to 'v3'")
    }
    
    if(platypus.version == "v3"){
      if(missing(VDJ.matrix)) stop("When using platypus version v3 please provide an input list to VDJ.matrix. If using the output from VDJ_clonotype, please switch platypus.version to 'v2'")
    }
    
    if(platypus.version == "v2"){ #old 
    
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
    
    
    
    } else if(platypus.version == "v3"){ #new
      
      
      Vgene_usage_matrix <- list()
      dummy_df <- list()
      
      
      #filtering for max 1VDJ 1VJ chain
      VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
      
      
      sample_list <- list()
      for(i in 1:length(unique(VDJ.matrix$sample_id))){
        sample_list[[i]] <- subset(VDJ.matrix, sample_id == unique(VDJ.matrix$sample_id)[i])
      }
      
      for (k in 1:length(sample_list)){
        
        #make a new column with both IgH and IgL/K V genes
        sample_list[[k]]$vgenes <- paste(sample_list[[k]]$VDJ_vgene, sample_list[[k]]$VJ_vgene, sep = "_")
        
        #Create a matrix with rows being heavy chain and columns being light chain v genes
        Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(sample_list[[k]]$VDJ_vgene)), ncol = length(unique(sample_list[[k]]$VJ_vgene)))
        
        #give the row names and col names
        rownames(Vgene_usage_matrix[[k]]) <- unique(sample_list[[k]]$VDJ_vgene)
        colnames(Vgene_usage_matrix[[k]]) <- unique(sample_list[[k]]$VJ_vgene)
        
        #create dummy df which will contain the counts for each combination
        dummy_df[[k]] <- as.data.frame(table(sample_list[[k]]$vgenes))
        colnames(dummy_df[[k]]) <- c("vgene", "count")
        
        #go elementwise in the matrix and count the occurancies of each combination in the sample_list$vgenes column
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
}
