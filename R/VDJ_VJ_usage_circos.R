#'Makes a Circos plot from the VDJ_analyze output. Connects the V gene with the corresponding J gene for each clonotype.
#' @param VDJ.GEX.matrix The output of the VDJ_GEX_integrate function (Platypus platypus.version v2). A list of data frames for each sample containing the clonotype information and cluster membership information. For Platypus platypus.version v3, VDJ_GEX_matrix.output[[1]] has to be supplied.
#' @param A.or.B Determines whether to plot the V J gene pairing of the alpha or beta chain. "A", "B" or "both" as possible inputs. Default: "both".
#' @param label.threshold Minimal amount of clonotypes per gene neccessary to add a gene label to the sector. Default: 0.
#' @param c.threshold Only clonotypes are considered with a frequency higher then c.threshold. Allows to filter for only highly expanded clonotypes.
#' @param cell.level Logical, defines whether weight of connection should be based on number of clonotypes of number of cells. Default: number of clonotypes.
#' @param clonotype.per.gene.threshold How many clonotypes are required to plot a sector for a gene. Filters the rows and colums of the final adjacency matrix.
#' @param c.count Show clonotype or cell count on Circos plot. Default = T.
#' @param platypus.version Which platypus.version of platypus is being used. Default = v2.
#' @param filter1H1L Whether to filter the input VDJ.matrix in "v3" to only include cells with 1 VDJ and 1 VJ chain. Defaults to TRUE
#' @return Returns list of plots. The first n elements contain the circos plot of the n datasets from the VDJ.analyze function. The n+1 element contains a list of the n adjancey matrices for each dataset.
#' @examples
#' \dontrun{
#'  plots <- VJ_alpha_beta_Vgene_circos(VDJ.GEX.matrix)
#'}
#' @export

VDJ_VJ_usage_circos <- function(VDJ.GEX.matrix, A.or.B, label.threshold, cell.level, c.threshold, clonotype.per.gene.threshold, c.count, platypus.version, filter1H1L){
  require(reshape2)
  require(ggplot2)
  require(stringr)
  if(missing(label.threshold)){label.threshold <- 1}
  if(missing(A.or.B)){A.or.B <- "both"}
  if(missing(cell.level)){cell.level <- F}
  if(missing(c.threshold)){c.threshold <- 0}
  if(missing(clonotype.per.gene.threshold)){clonotype.per.gene.threshold <- 0}
  if(missing(c.count)){c.count <- T}
  if(missing(platypus.version)){platypus.version <- "v2"}
  if(missing(filter1H1L)){filter1H1L <- T}

  # If new version with VDJ_GEX_matric output should be used
  if(platypus.version=="v3"){
      print("Reminder: VDJ_VJ_usage_circos() funcion built for new Platypus 3.0.0 is being used. The output VDJ dataframe of the VDJ_GEX_matrix required as input.")
      clonotype <- "clonotype_id_10x"

      if(class(VDJ.GEX.matrix) == "list"){
        stop("Please input the VDJ matrix from the VGM output (usually VDJ_GEX_matrix.output[[1]])")
      }
      #swapping to a list to not change the whole function
      bk <- VDJ.GEX.matrix
      VDJ.GEX.matrix <- list()
      VDJ.GEX.matrix[[1]] <- bk

      #filter for 1H1L
      if(filter1H1L==T){
        VDJ.GEX.matrix[[1]]<-VDJ.GEX.matrix[[1]][which((VDJ.GEX.matrix[[1]]$Nr_of_VDJ_chains==1)&(VDJ.GEX.matrix[[1]]$Nr_of_VJ_chains==1)),]
      }

      #filter out clonotypes with less then c.threshold cells
        VDJ.GEX.matrix[[1]] <-VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$clonotype_frequency >= c.threshold),]

      #replace empty fields "" with "None"
        VDJ.GEX.matrix[[1]]$VJ_vgene[which(VDJ.GEX.matrix[[1]]$VJ_vgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VJ_jgene[which(VDJ.GEX.matrix[[1]]$VJ_jgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VDJ_vgene[which(VDJ.GEX.matrix[[1]]$VDJ_vgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VDJ_jgene[which(VDJ.GEX.matrix[[1]]$VDJ_jgene == "")] <- "None"

      TRBV <- c()
      TRAV <- c()
      TRBJ <- c()
      TRAJ <- c()

      #split up into samples to plot individually
      VDJ.GEX_list <- list()
      for (i in 1:length(unique(VDJ.GEX.matrix[[1]]$sample_id))){
        VDJ.GEX_list[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id== unique(VDJ.GEX.matrix[[1]]$sample_id)[i]),]
      }


      for (i in 1:length(VDJ.GEX_list)){

      print(paste0("Plotting for sample ", unique(VDJ.GEX.matrix[[1]]$sample_id)[i]))

        VDJ.GEX_list[[i]]$alpha_VJ_gene <- paste(VDJ.GEX_list[[i]]$VJ_vgene, VDJ.GEX_list[[i]]$VJ_jgene, sep = "_")
        VDJ.GEX_list[[i]]$beta_VJ_gene <- paste(VDJ.GEX_list[[i]]$VDJ_vgene, VDJ.GEX_list[[i]]$VDJ_jgene, sep = "_")

        #get all V and J genes across all samples
        c <- unique(VDJ.GEX_list[[i]]$VJ_vgene)
        TRBV <- append(TRBV,c)
        d <- unique(VDJ.GEX_list[[i]]$VDJ_vgene)
        TRAV <- append(TRAV,d)
        e <- unique(VDJ.GEX_list[[i]]$VJ_jgene)
        TRBJ<- append(TRBJ,e)
        f <- unique(VDJ.GEX_list[[i]]$VDJ_jgene)
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

      for (k in 1:length(VDJ.GEX_list)){
        #Create a matrix with rows being heavy chain and columns being light chain v genes
        Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRV)), ncol = length(unique(TRJ)))

        #give the row names and col names
        rownames(Vgene_usage_matrix[[k]]) <- unique(TRV)
        colnames(Vgene_usage_matrix[[k]]) <- unique(TRJ)
      }

      for (k in 1:length(VDJ.GEX_list)){

        #create dummy df which will contain the counts for each combination
        if(cell.level==T){
          dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_VJ_gene))
          colnames(dummy_alpha_df[[k]]) <- c("vjgene", "count")
          dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$beta_VJ_gene))
          colnames(dummy_beta_df[[k]]) <- c("vjgene", "count")
        }else{
          print("---")
          print(paste0("Processing sample ", k))
          print("WARNING: If clonotype strategy is not based on unique V or J genes per clonotype, this setting [cell.level=F] might be questionable. One clonotype might then be represented in several Circos connections between V or J genes. The names of genes of simulatneously used chains will be pasted together.")
          print(paste("Chosen clonotype column: ", clonotype))
          print("WARNING: If Circos plotting error occurs: Maybe your `gap.degree` is too large so that there is no space to allocate sectors -> You might want to increase clonotype.per.gene.threshold to reduce number of sectors in your Circos plots")

          dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype]],VDJ.GEX_list[[k]]$alpha_VJ_gene, sep="/and/")))
          colnames(dummy) <- c("pasted")
          dummy$clonotype <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)
          dummy$gene <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)[,2]

          dummy_alpha_df[[k]] <- as.data.frame(table(dummy$gene))
          colnames(dummy_alpha_df[[k]]) <- c("vjgene", "count")


          dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype]],VDJ.GEX_list[[k]]$beta_VJ_gene, sep="/and/")))
          colnames(dummy) <- c("pasted")
          dummy$clonotype <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)
          dummy$gene <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)[,2]

          dummy_beta_df[[k]] <- as.data.frame(table(dummy$gene))
          colnames(dummy_beta_df[[k]]) <- c("vjgene", "count")

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
      group = structure(stringr::str_sub(nm, 1,4), names = nm)
      group = factor(group)
      #set colors: same color across all samples
      grid.col <- stats::setNames(grDevices::rainbow(length(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]])))),sample(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]]))))

      # filter out genes which have fewer clonotypes then clonotype.per.gene.threshold
      for (k in 1:length(Vgene_usage_matrix)){
        Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][which(rowSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold),]
        Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][,which(colSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold)]

      }

  }else if (platypus.version=="v2"){

  # In Case old platypus.version of Platypus is being used: VDJ_analyze output is the input

    print("Reminder: VDJ_VJ_usage_circos() funcion built for Platypus 2.0.5 is being used. Output of VDJ_analyze() required as input. Set [platypus.version = 'v3'] for compatibility with VDJ_GEX_matrix().")
    for (i in 1:length(VDJ.GEX.matrix)) {
      VDJ.GEX.matrix[[i]] <- VDJ.GEX.matrix[[i]][which(VDJ.GEX.matrix[[i]]$frequency >=
                                                                 c.threshold), ]
    }
    TRBV <- c()
    TRAV <- c()
    TRBJ <- c()
    TRAJ <- c()
    for (i in 1:length(VDJ.GEX.matrix)) {
      VDJ.GEX.matrix[[i]]$alpha_VJ_gene <- paste(VDJ.GEX.matrix[[i]]$HC_vgene,
                                                     VDJ.GEX.matrix[[i]]$HC_jgene, sep = "_")
      VDJ.GEX.matrix[[i]]$beta_VJ_gene <- paste(VDJ.GEX.matrix[[i]]$LC_vgene,
                                                    VDJ.GEX.matrix[[i]]$LC_jgene, sep = "_")
      c <- unique(VDJ.GEX.matrix[[i]]$HC_vgene)
      TRBV <- append(TRBV, c)
      d <- unique(VDJ.GEX.matrix[[i]]$LC_vgene)
      TRAV <- append(TRAV, d)
      e <- unique(VDJ.GEX.matrix[[i]]$HC_jgene)
      TRBJ <- append(TRBJ, e)
      f <- unique(VDJ.GEX.matrix[[i]]$LC_jgene)
      TRAJ <- append(TRAJ, f)
    }
    if (A.or.B == "both") {
      TRV <- append(TRAV, TRBV)
      TRJ <- append(TRAJ, TRBJ)
    }
    else if (A.or.B == "A") {
      TRV <- TRAV
      TRJ <- TRAJ
    }
    else if (A.or.B == "B") {
      TRV <- TRBV
      TRJ <- TRBJ
    }
    else {
      print("Please specify A.or.B or leave empty to plot both.")
      return()
    }
    Vgene_usage_matrix <- list()
    dummy_alpha_df <- list()
    dummy_beta_df <- list()
    for (k in 1:length(VDJ.GEX.matrix)) {
      Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRV)),
                                        ncol = length(unique(TRJ)))
      rownames(Vgene_usage_matrix[[k]]) <- unique(TRV)
      colnames(Vgene_usage_matrix[[k]]) <- unique(TRJ)
    }
    for (k in 1:length(VDJ.GEX.matrix)) {
      if (cell.level == F) {
        dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_VJ_gene))
        colnames(dummy_alpha_df[[k]]) <- c("vjgene",
                                           "count")
        dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$beta_VJ_gene))
        colnames(dummy_beta_df[[k]]) <- c("vjgene",
                                          "count")
      }
      else {
        dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_VJ_gene))
        colnames(dummy_alpha_df[[k]]) <- c("vjgene",
                                           "clonotype_counts")
        for (i in 1:nrow(dummy_alpha_df[[k]])) {
          dummy_alpha_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$alpha_VJ_gene ==
                                                                              dummy_alpha_df[[k]]$vjgene[i]), "frequency"])
        }
        dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$beta_VJ_gene))
        colnames(dummy_beta_df[[k]]) <- c("vjgene",
                                          "clonotype_counts")
        for (i in 1:nrow(dummy_beta_df[[k]])) {
          dummy_beta_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$beta_VJ_gene ==
                                                                             dummy_beta_df[[k]]$vjgene[i]), "frequency"])
        }
      }
      for (i in 1:nrow(Vgene_usage_matrix[[k]])) {
        for (j in 1:ncol(Vgene_usage_matrix[[k]])) {
          if (paste0(rownames(Vgene_usage_matrix[[k]])[i],
                     "_", colnames(Vgene_usage_matrix[[k]])[j]) %in%
              dummy_alpha_df[[k]]$vjgene) {
            Vgene_usage_matrix[[k]][i, j] <- dummy_alpha_df[[k]][which(dummy_alpha_df[[k]]$vjgene ==
                                                                         paste0(rownames(Vgene_usage_matrix[[k]])[i],
                                                                                "_", colnames(Vgene_usage_matrix[[k]])[j])),
                                                                 "count"]
          }
          else {
            Vgene_usage_matrix[[k]][i, j] <- 0
          }
        }
      }
      for (i in 1:nrow(Vgene_usage_matrix[[k]])) {
        for (j in 1:ncol(Vgene_usage_matrix[[k]])) {
          if (paste0(rownames(Vgene_usage_matrix[[k]])[i],
                     "_", colnames(Vgene_usage_matrix[[k]])[j]) %in%
              dummy_beta_df[[k]]$vjgene) {
            Vgene_usage_matrix[[k]][i, j] <- dummy_beta_df[[k]][which(dummy_beta_df[[k]]$vjgene ==
                                                                        paste0(rownames(Vgene_usage_matrix[[k]])[i],
                                                                               "_", colnames(Vgene_usage_matrix[[k]])[j])),
                                                                "count"]
          }
        }
      }
    }
    nm = unique(unlist(dimnames(Vgene_usage_matrix[[1]])))
    group = structure(stringr::str_sub(nm, 1, 4), names = nm)
    grid.col <- stats::setNames(grDevices::rainbow(length(union(rownames(Vgene_usage_matrix[[1]]),
                                              colnames(Vgene_usage_matrix[[1]])))), sample(union(rownames(Vgene_usage_matrix[[1]]),
                                                                                                 colnames(Vgene_usage_matrix[[1]]))))
    for (k in 1:length(Vgene_usage_matrix)) {
      Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][which(rowSums(Vgene_usage_matrix[[k]]) >=
                                                                 clonotype.per.gene.threshold), ]
      Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][,
                                                         which(colSums(Vgene_usage_matrix[[k]]) >= clonotype.per.gene.threshold)]
    }


    #--- Forward Vgene_usage_matrix to Plotting function
  }else{
    print("Please specify platypus platypus.version as either v2 or v3.")
  }
  plots <- list()
  for (i in 1:length(Vgene_usage_matrix)){
    plots[[i]] <- VDJ_circos(Vgene_usage_matrix[[i]], group = group, grid.col=grid.col, label.threshold = label.threshold, c.count=c.count)
  }
  plots[[i+1]] <- Vgene_usage_matrix
  return(plots)
}
