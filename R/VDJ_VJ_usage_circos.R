#'Makes a Circos plot from the VDJ_analyze output. Connects the V gene with the corresponding J gene for each clonotype.
#' @param VDJ.analyze.output The output of the VDJ_analyze function. A list of data frames for each sample containing the clonotype informations.
#' @param A.or.B Determines whether to plot the V J gene pairing of the alpha or beta chain. "A", "B" or "both" as possible inputs. Default: "both".
#' @param label.threshold Minimal amount of clonotypes per gene neccessary to add a gene label to the sector. Default: 0.
#' @param c.threshold Only clonotypes are considered with a frequency higher then c.threshold. Allows to filter for only highly expanded clonotypes.
#' @param cell.level Logical, defines whether weight of connection should be based on number of clonotypes of number of cells. Default: number of clonotypes.
#' @param clonotype.per.gene.threshold How many clonotypes are required to plot a sector for a gene. Filters the rows and colums of the final adjacency matrix.
#' @return Returns list of plots. The first n elements contain the circos plot of the n datasets from the VDJ.analyze function. The n+1 element contains a list of the n adjancey matrices for each dataset.
#' @examples
#' \dontrun{
#'  plots <- VJ_alpha_beta_Vgene_circos(VDJ.analyze.output)
#'}
#' @export

VDJ_VJ_usage_circos <- function(VDJ.analyze.output, A.or.B, label.threshold, cell.level, c.threshold, clonotype.per.gene.threshold){
  if(missing(label.threshold)){label.threshold <- 1}
  if(missing(A.or.B)){A.or.B <- "both"}
  if(missing(cell.level)){cell.level <- F}
  if(missing(c.threshold)){c.threshold <- 0}
  if(missing(clonotype.per.gene.threshold)){clonotype.per.gene.threshold <- 0}

  #filter out clonotypes with less then c.threshold cells
  for(i in 1:length(VDJ.analyze.output)){
    VDJ.analyze.output[[i]] <-VDJ.analyze.output[[i]][which(VDJ.analyze.output[[i]]$frequency >= c.threshold),]
  }


  TRBV <- c()
  TRAV <- c()
  TRBJ <- c()
  TRAJ <- c()

  for (i in 1:length(VDJ.analyze.output)){

    VDJ.analyze.output[[i]]$alpha_VJ_gene <- paste(VDJ.analyze.output[[i]]$HC_vgene, VDJ.analyze.output[[i]]$HC_jgene, sep = "_")
    VDJ.analyze.output[[i]]$beta_VJ_gene <- paste(VDJ.analyze.output[[i]]$LC_vgene, VDJ.analyze.output[[i]]$LC_jgene, sep = "_")

    #get all V and J genes across all samples
    c <- unique(VDJ.analyze.output[[i]]$HC_vgene)
    TRBV <- append(TRBV,c)
    d <- unique(VDJ.analyze.output[[i]]$LC_vgene)
    TRAV <- append(TRAV,d)
    e <- unique(VDJ.analyze.output[[i]]$HC_jgene)
    TRBJ<- append(TRBJ,e)
    f <- unique(VDJ.analyze.output[[i]]$LC_jgene)
    TRAJ <- append(TRAJ,f)
  }
  #make two branches for matrix V vs J genes
  if(A.or.B=="both"){
    TRV <- append(TRAV, TRBV)
    TRJ <- append(TRAJ, TRBJ)
  }else if(A.or.B=="A"){
    TRV <- TRAV
    TRJ <- TRAJ
  }else if(A.or.B=="B"){
    TRV <- TRBV
    TRJ <- TRBJ
  }else{
    print("Please specify A.or.B or leave empty to plot both.")
    return()
  }

  # create matrix Vgenes vs Jgenes
  Vgene_usage_matrix <- list()
  dummy_alpha_df <- list()
  dummy_beta_df <- list()

  for (k in 1:length(VDJ.analyze.output)){
    #Create a matrix with rows being heavy chain and columns being light chain v genes
    Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRV)), ncol = length(unique(TRJ)))

    #give the row names and col names
    rownames(Vgene_usage_matrix[[k]]) <- unique(TRV)
    colnames(Vgene_usage_matrix[[k]]) <- unique(TRJ)
  }

  for (k in 1:length(VDJ.analyze.output)){

    #create dummy df which will contain the counts for each combination
    if(cell.level==F){
      dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.analyze.output[[k]]$alpha_VJ_gene))
      colnames(dummy_alpha_df[[k]]) <- c("vjgene", "count")
      dummy_beta_df[[k]] <- as.data.frame(table(VDJ.analyze.output[[k]]$beta_VJ_gene))
      colnames(dummy_beta_df[[k]]) <- c("vjgene", "count")
    }else{
      dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.analyze.output[[k]]$alpha_VJ_gene))
      colnames(dummy_alpha_df[[k]]) <- c("vjgene", "clonotype_counts")
      for(i in 1:nrow(dummy_alpha_df[[k]])){
        dummy_alpha_df[[k]]$count[i] <- sum(VDJ.analyze.output[[k]][which(VDJ.analyze.output[[k]]$alpha_VJ_gene==dummy_alpha_df[[k]]$vjgene[i]),"frequency"])
      }
      dummy_beta_df[[k]] <- as.data.frame(table(VDJ.analyze.output[[k]]$beta_VJ_gene))
      colnames(dummy_beta_df[[k]]) <- c("vjgene", "clonotype_counts")
      for(i in 1:nrow(dummy_beta_df[[k]])){
        dummy_beta_df[[k]]$count[i] <- sum(VDJ.analyze.output[[k]][which(VDJ.analyze.output[[k]]$beta_VJ_gene==dummy_beta_df[[k]]$vjgene[i]),"frequency"])
      }
    }

    #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
    for (i in 1:nrow(Vgene_usage_matrix[[k]])){
      for (j in 1:ncol(Vgene_usage_matrix[[k]])){
        if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_alpha_df[[k]]$vjgene){
          Vgene_usage_matrix[[k]][i,j] <- dummy_alpha_df[[k]][which(dummy_alpha_df[[k]]$vjgene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
        } else {
          Vgene_usage_matrix[[k]][i,j] <- 0
        }

      }
    }
    #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
    for (i in 1:nrow(Vgene_usage_matrix[[k]])){
      for (j in 1:ncol(Vgene_usage_matrix[[k]])){
        if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_beta_df[[k]]$vjgene){
          Vgene_usage_matrix[[k]][i,j] <- dummy_beta_df[[k]][which(dummy_beta_df[[k]]$vjgene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
        }
        #no else clause. would overwrite alpha chain genes in matrix    }
      }
    }

  }



  # group by TRAV;TRBV; TRAJ and TRBJ
  nm = unique(unlist(dimnames(Vgene_usage_matrix[[1]])))
  group = structure(str_sub(nm, 1,4), names = nm)

  #set colors: same color across all samples
  grid.col <- setNames(rainbow(length(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]])))),sample(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]]))))

  # filter out genes which have fewer clonotypes then clonotype.per.gene.threshold
  for (k in 1:length(Vgene_usage_matrix)){
    Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][which(rowSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold),]
    Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][,which(colSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold)]

  }

  plots <- list()
  for (i in 1:length(Vgene_usage_matrix)){
    plots[[i]] <- VDJ_circos(Vgene_usage_matrix[[i]], group = group, grid.col=grid.col, label.threshold = label.threshold)
  }
  plots[[i+1]] <- Vgene_usage_matrix
  return(plots)
}
