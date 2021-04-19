#'Makes a Circos plot from the VDJ_analyze output. Connects the V-alpha with the corresponding V-beta gene for each clonotype.
#' @param VDJ.GEX.matrix The output of the VDJ_GEX_integrate function (Platypus platypus.version v2). A list of data frames for each sample containing the clonotype information and cluster membership information. For Platypus platypus.version v3, the VDJ_GEX_matrix() output has to be supplied.
#' @param V.or.J Determines whether to plot the alpha beta gene pairing of the V or J genes. "V", "J" or "both" as possible inputs. Default: "both".
#' @param label.threshold Minimal amount of clonotypes per gene neccessary to add a gene label to the sector. Default: 0.
#' @param c.threshold Only clonotypes are considered with a frequency higher then c.threshold. Allows to filter for only highly expanded clonotypes.
#' @param cell.level Logical, defines whether weight of connection should be based on number of clonotypes of number of cells. Default: number of clonotypes.
#' @param clonotype.per.gene.threshold How many clonotypes are required to plot a sector for a gene. Filters the rows and colums of the final adjacency matrix.
#' @param B.or.Tcells Specify whether B or T cells are being analyzed ("B" or "T"). If not specified, function attempts to decide based on gene names.
#' @param c.count Show clonotype or cell count on Circos plot. Default = T.
#' @param platypus.version Which platypus.version of platypus is beeing used. Default = v2.
#' @return Returns list of plots. The first n elements contain the circos plot of the n datasets from the VDJ.analyze function. The n+1 element contains a list of the n adjancey matrices for each dataset.
#' @examples
#' \dontrun{
#'  plots <- VJ_alpha_beta_Vgene_circos(vdj_repertoire_tcells)
#'}
#' @export

VJ_alpha_beta_Vgene_circos <- function(VDJ.GEX.matrix, V.or.J, B.or.Tcells, label.threshold, c.threshold, cell.level, clonotype.per.gene.threshold, c.count, platypus.version, filter1H1L){
if(missing(V.or.J)){V.or.J <- "both"}
if(missing(label.threshold)){label.threshold <- 0}
if(missing(c.threshold)){c.threshold <- 0}
if(missing(cell.level)){cell.level <- F}
if(missing(clonotype.per.gene.threshold)){clonotype.per.gene.threshold <- 0}
if(missing(c.count)){c.count <- T}
if(missing(platypus.version)){platypus.version ="v2"}
if(missing(filter1H1L)){filter1H1L <- T}

  
if(platypus.version=="v3"){
  #########################################################################
  print("Reminder: VDJ_VJ_usage_circos() funcion built for new Platypus v3.0.0 is being used. Output of VDJ_GEX_matrix() required as input.")
  clonotype <- "clonotype_id_10x"
  
  if(missing(B.or.Tcells)){
    for(i in 1:nrow(VDJ.GEX.matrix[[1]])){
      if(substr(VDJ.GEX.matrix[[1]]$VDJ_vgene[[i]],start=1, stop = 2)=="IG"){
        B.or.Tcells <- "B"
        print("Bcells used")
        break
      }
      if(substr(VDJ.GEX.matrix[[1]]$VDJ_vgene[[1]],start=1, stop = 2)=="TR"){
        B.or.Tcells <- "T"
        print("Tcells used")
        break
      }
      if(i == nrow(VDJ.GEX.matrix[[1]])){
        print("Please specify whether B or T cells are beeing analyzed (Parameter B.or.Tcells)")
        break
      }
    }
  }
  
  plots <- list()
  
  #filter for 1H1L
  
  if(filter1H1L==T){
    VDJ.GEX.matrix[[1]]<-VDJ.GEX.matrix[[1]][which((VDJ.GEX.matrix[[1]]$Nr_of_VDJ_chains==1)&(VDJ.GEX.matrix[[1]]$Nr_of_VJ_chains==1)),]
  }
  
  #filter out clonotypes with less then c.threshold cells
    VDJ.GEX.matrix[[1]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$clonotype_frequency >= c.threshold),]

  
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
  for (i in 1:length(table(VDJ.GEX.matrix[[1]]$sample_id))){
    VDJ.GEX_list[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==paste0("s",i)),]
  }
  
  
  for (i in 1:length(VDJ.GEX_list)){
    
    VDJ.GEX_list[[i]]$alpha_beta_Vgene <- paste(VDJ.GEX_list[[i]]$VJ_vgene, VDJ.GEX_list[[i]]$VDJ_vgene, sep = "_")
    VDJ.GEX_list[[i]]$alpha_beta_Jgene <- paste(VDJ.GEX_list[[i]]$VJ_jgene, VDJ.GEX_list[[i]]$VDJ_jgene, sep = "_")
    
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
  #make two branches for matrix alpha vs beta genes
  if(V.or.J == "both"){
    TRA <- append(TRAV, TRAJ)
    TRB <- append(TRBV, TRBJ)
  }else if(V.or.J == "V"){
    TRA <- TRAV
    TRB <- TRBV
  }else if(V.or.J == "J"){
    TRA <- TRAJ
    TRB <- TRBJ
  }else{
    print("Please specify V.or.J as 'V' or 'J'. Leave empty to plot both")
  }
  
  
  # create matrix alpha vs beta
  Vgene_usage_matrix <- list()
  dummy_Vgene_df <- list()
  dummy_Jgene_df <- list()
  
  
  for (k in 1:length(VDJ.GEX_list)){
    #Create a matrix with rows being heavy chain and columns being light chain v genes
    Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRA)), ncol = length(unique(TRB)))
    
    #give the row names and col names
    rownames(Vgene_usage_matrix[[k]]) <- unique(TRA)
    colnames(Vgene_usage_matrix[[k]]) <- unique(TRB)
  }
  
  for (k in 1:length(VDJ.GEX_list)){
    
    #create dummy df which will contain the counts for each combination
    
    if(cell.level == T){
      print("---")
      print(paste0("Processing sample ", k))
      dummy_Vgene_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_beta_Vgene))
      colnames(dummy_Vgene_df[[k]]) <- c("gene", "count")
      dummy_Jgene_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_beta_Jgene))
      colnames(dummy_Jgene_df[[k]]) <- c("gene", "count")
    }else{
      print("---")
      print(paste0("Processing sample ", k))
      print("WARNING: If clonotype strategy is not based on unique V or J genes per clonotype, this setting [cell.level=F] might be questionable. One clonotype might then be represented in several Circos connections between V or J genes. The names of genes of simulatneously used chains will be pasted together.")
      print(paste("Chosen clonotype column: ", clonotype))
      print("WARNING: If Circos plotting error occurs: Maybe your `gap.degree` is too large so that there is no space to allocate sectors -> You might want to increase clonotype.per.gene.threshold to reduce number of sectors in your Circos plots")
      
      dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype]],VDJ.GEX_list[[k]]$alpha_beta_Vgene, sep="/and/")))
      colnames(dummy) <- c("pasted")
      dummy$clonotype <- str_split_fixed(dummy$pasted, "/and/", 2)
      dummy$gene <- str_split_fixed(dummy$pasted, "/and/", 2)[,2]
      
      dummy_Vgene_df[[k]] <- as.data.frame(table(dummy$gene))
      colnames(dummy_Vgene_df[[k]]) <- c("gene", "count")
      
      
      dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype]],VDJ.GEX_list[[k]]$alpha_beta_Jgene, sep="/and/")))
      colnames(dummy) <- c("pasted")
      dummy$clonotype <- str_split_fixed(dummy$pasted, "/and/", 2)[,1]
      dummy$gene <- str_split_fixed(dummy$pasted, "/and/", 2)[,2]
      
      dummy_Jgene_df[[k]] <- as.data.frame(table(dummy$gene))
      colnames(dummy_Jgene_df[[k]]) <- c("gene", "count")
      
      # dummy_Vgene_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_beta_Vgene))
      # colnames(dummy_Vgene_df[[k]]) <- c("gene", "clonotype_counts")
      # for(i in 1:nrow(dummy_Vgene_df[[k]])){
      #   dummy_Vgene_df[[k]]$count[i] <- sum(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]]$alpha_beta_Vgene==dummy_Vgene_df[[k]]$gene[i]),"frequency"])
      # }
      # dummy_Jgene_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_beta_Jgene))
      # colnames(dummy_Jgene_df[[k]]) <- c("gene", "clonotype_counts")
      # for(i in 1:nrow(dummy_Jgene_df[[k]])){
      #   dummy_Jgene_df[[k]]$count[i] <- sum(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]]$alpha_beta_Jgene==dummy_Jgene_df[[k]]$gene[i]),"frequency"])
      # }
    }
    
    #go elementwise in the matrix and count the occurancies of each combination
    for (i in 1:nrow(Vgene_usage_matrix[[k]])){
      for (j in 1:ncol(Vgene_usage_matrix[[k]])){
        if (paste0(colnames(Vgene_usage_matrix[[k]])[j], "_",rownames(Vgene_usage_matrix[[k]])[i]) %in% dummy_Vgene_df[[k]]$gene){
          Vgene_usage_matrix[[k]][i,j] <- dummy_Vgene_df[[k]][which(dummy_Vgene_df[[k]]$gene == paste0(colnames(Vgene_usage_matrix[[k]])[j], "_",rownames(Vgene_usage_matrix[[k]])[i])),"count"]
        } else {
          Vgene_usage_matrix[[k]][i,j] <- 0
        }
        
      }
    }
    
    
    #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
    for (i in 1:nrow(Vgene_usage_matrix[[k]])){
      for (j in 1:ncol(Vgene_usage_matrix[[k]])){
        if (paste0(colnames(Vgene_usage_matrix[[k]])[j], "_",rownames(Vgene_usage_matrix[[k]])[i]) %in% dummy_Jgene_df[[k]]$gene){
          Vgene_usage_matrix[[k]][i,j] <- dummy_Jgene_df[[k]][which(dummy_Jgene_df[[k]]$gene == paste0(colnames(Vgene_usage_matrix[[k]])[j], "_",rownames(Vgene_usage_matrix[[k]])[i])),"count"]
        }
        #no else clause. would overwrite alpha chain genes in matrix    }
      }
    }
    
  }
  
  #set colors: same color across all samples
  grid.col <- setNames(rainbow(length(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]])))),sample(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]]))))
  
  for (i in 1:length(Vgene_usage_matrix)){
    # filter out genes which have fewer clonotypes then clonotype.per.gene.threshold
    Vgene_usage_matrix[[i]] <- Vgene_usage_matrix[[i]][which(rowSums(Vgene_usage_matrix[[i]])>=clonotype.per.gene.threshold),]
    Vgene_usage_matrix[[i]] <- Vgene_usage_matrix[[i]][,which(colSums(Vgene_usage_matrix[[i]])>=clonotype.per.gene.threshold)]
    
    nm = unique(unlist(dimnames(Vgene_usage_matrix[[i]])))
    group = structure(str_sub(nm, 1,4), names = nm)
    
    
    
    
    #set levels for grouping.
    #Differentiate whether B or T cells are analyzed
    
    
    #Differentiate whether "None" has to be included as its own group for unpaired clonotypes.
    if(B.or.Tcells == "T"){
      if(V.or.J == "both"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("TRAV", "TRBV", "TRAJ","TRBJ","None")
        }else{
          levels <- c("TRAV", "TRBV", "TRAJ","TRBJ")
        }
      }else if(V.or.J == "V"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("TRAV", "TRBV","None")
        }else{
          levels <- c("TRAV", "TRBV")
        }
      }else if(V.or.J == "J"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("TRAJ","TRBJ","None")
        }else{
          levels <- c("TRAJ", "TRBJ")
        }
      }
    }else if(B.or.Tcells == "B"){
      if(V.or.J == "both"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("IGHV","IGLV","IGKV","IGHJ","IGLJ","IGKJ","None")
        }else{
          levels <- c("IGHV","IGLV","IGKV","IGHJ","IGLJ","IGKJ")
        }
      }else if(V.or.J == "V"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("IGHV","IGLV","IGKV","None")
        }else{
          levels <- c("IGHV","IGLV","IGKV")
        }
      }else if(V.or.J == "J"){
        if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
          levels <- c("IGHJ","IGLJ","IGKJ","None")
        }else{
          levels <- c("IGHJ","IGLJ","IGKJ")
        }
      }
    }else{
      print("Please specify whether B or T cells are analyzed. (Parameter B.or.Tcells)")
    }
    
    # Set grouping factors based on previously defined levels. Order of levels defines order of groups in Circos plot.
    group = factor(group[sample(length(group), length(group))], levels = levels)
    
    plots[[i]] <- VDJ_circos(Vgene_usage_matrix[[i]], group = group, grid.col=grid.col, label.threshold = label.threshold, c.count = c.count)
  }
  plots[[i+1]] <- Vgene_usage_matrix
  
  

  
}else if(platypus.version=="v2"){
  
  
  ###############################################################################
  print("Reminder: VDJ_VJ_usage_circos() funcion built for Platypus 2.0.5 is being used. Output of VDJ_analyze() required as input. Set [platypus.version = new] for compatibility with VDJ_GEX_matrix().")
  
  if(missing(B.or.Tcells)){
    for(i in 1:nrow(VDJ.GEX.matrix[[1]])){
      if(substr(VDJ.GEX.matrix[[1]]$HC_vgene[[i]],start=1, stop = 2)=="IG"){
        B.or.Tcells <- "B"
        break
      }
      if(substr(VDJ.GEX.matrix[[1]]$HC_vgene[[1]],start=1, stop = 2)=="TR"){
        B.or.Tcells <- "T"
        break
      }
      if(i == nrow(VDJ.GEX.matrix[[1]])){
        print("Please specify whether B or T cells are beeing analyzed (Parameter B.or.Tcells)")
        break
      }
    }
  }
    
      plots <- list()
    
      #filter out clonotypes with less then c.threshold cells
      for(i in 1:length(VDJ.GEX.matrix)){
        VDJ.GEX.matrix[[i]] <- VDJ.GEX.matrix[[i]][which(VDJ.GEX.matrix[[i]]$frequency >= c.threshold),]
      }
    
      TRBV <- c()
      TRAV <- c()
      TRBJ <- c()
      TRAJ <- c()
    
      for (i in 1:length(VDJ.GEX.matrix)){
    
        VDJ.GEX.matrix[[i]]$alpha_beta_Vgene <- paste(VDJ.GEX.matrix[[i]]$LC_vgene, VDJ.GEX.matrix[[i]]$HC_vgene, sep = "_")
        VDJ.GEX.matrix[[i]]$alpha_beta_Jgene <- paste(VDJ.GEX.matrix[[i]]$LC_jgene, VDJ.GEX.matrix[[i]]$HC_jgene, sep = "_")
    
        #get all V and J genes across all samples
        c <- unique(VDJ.GEX.matrix[[i]]$HC_vgene)
        TRBV <- append(TRBV,c)
        d <- unique(VDJ.GEX.matrix[[i]]$LC_vgene)
        TRAV <- append(TRAV,d)
        e <- unique(VDJ.GEX.matrix[[i]]$HC_jgene)
        TRBJ<- append(TRBJ,e)
        f <- unique(VDJ.GEX.matrix[[i]]$LC_jgene)
        TRAJ <- append(TRAJ,f)
      }
      #make two branches for matrix alpha vs beta genes
      if(V.or.J == "both"){
        TRA <- append(TRAV, TRAJ)
        TRB <- append(TRBV, TRBJ)
      }else if(V.or.J == "V"){
        TRA <- TRAV
        TRB <- TRBV
      }else if(V.or.J == "J"){
        TRA <- TRAJ
        TRB <- TRBJ
      }else{
        print("Please specify V.or.J as 'V' or 'J'. Leave empty to plot both")
      }
    
    
      # create matrix alpha vs beta
      Vgene_usage_matrix <- list()
      dummy_Vgene_df <- list()
      dummy_Jgene_df <- list()
    
    
      for (k in 1:length(VDJ.GEX.matrix)){
        #Create a matrix with rows being heavy chain and columns being light chain v genes
        Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRA)), ncol = length(unique(TRB)))
    
        #give the row names and col names
        rownames(Vgene_usage_matrix[[k]]) <- unique(TRA)
        colnames(Vgene_usage_matrix[[k]]) <- unique(TRB)
      }
    
      for (k in 1:length(VDJ.GEX.matrix)){
    
    
    
        #create dummy df which will contain the counts for each combination
    
        if(cell.level == F){
          dummy_Vgene_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_beta_Vgene))
          colnames(dummy_Vgene_df[[k]]) <- c("gene", "count")
          dummy_Jgene_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_beta_Jgene))
          colnames(dummy_Jgene_df[[k]]) <- c("gene", "count")
        }else{
          dummy_Vgene_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_beta_Vgene))
          colnames(dummy_Vgene_df[[k]]) <- c("gene", "clonotype_counts")
          for(i in 1:nrow(dummy_Vgene_df[[k]])){
            dummy_Vgene_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$alpha_beta_Vgene==dummy_Vgene_df[[k]]$gene[i]),"frequency"])
          }
          dummy_Jgene_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_beta_Jgene))
          colnames(dummy_Jgene_df[[k]]) <- c("gene", "clonotype_counts")
          for(i in 1:nrow(dummy_Jgene_df[[k]])){
            dummy_Jgene_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$alpha_beta_Jgene==dummy_Jgene_df[[k]]$gene[i]),"frequency"])
          }
        }
    
        #go elementwise in the matrix and count the occurancies of each combination
        for (i in 1:nrow(Vgene_usage_matrix[[k]])){
          for (j in 1:ncol(Vgene_usage_matrix[[k]])){
            if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_Vgene_df[[k]]$gene){
              Vgene_usage_matrix[[k]][i,j] <- dummy_Vgene_df[[k]][which(dummy_Vgene_df[[k]]$gene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
            } else {
              Vgene_usage_matrix[[k]][i,j] <- 0
            }
    
          }
        }
    
    
        #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
        for (i in 1:nrow(Vgene_usage_matrix[[k]])){
          for (j in 1:ncol(Vgene_usage_matrix[[k]])){
            if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_Jgene_df[[k]]$gene){
              Vgene_usage_matrix[[k]][i,j] <- dummy_Jgene_df[[k]][which(dummy_Jgene_df[[k]]$gene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
            }
            #no else clause. would overwrite alpha chain genes in matrix    }
          }
        }
    
      }
    
      #set colors: same color across all samples
      grid.col <- setNames(rainbow(length(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]])))),sample(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]]))))
    
      for (i in 1:length(Vgene_usage_matrix)){
        # filter out genes which have fewer clonotypes then clonotype.per.gene.threshold
        Vgene_usage_matrix[[i]] <- Vgene_usage_matrix[[i]][which(rowSums(Vgene_usage_matrix[[i]])>=clonotype.per.gene.threshold),]
        Vgene_usage_matrix[[i]] <- Vgene_usage_matrix[[i]][,which(colSums(Vgene_usage_matrix[[i]])>=clonotype.per.gene.threshold)]
    
        nm = unique(unlist(dimnames(Vgene_usage_matrix[[i]])))
        group = structure(str_sub(nm, 1,4), names = nm)
    
    
    
    
        #set levels for grouping.
        #Differentiate whether B or T cells are analyzed
    
    
        #Differentiate whether "None" has to be included as its own group for unpaired clonotypes.
        if(B.or.Tcells == "T"){
            if(V.or.J == "both"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("TRAV", "TRBV", "TRAJ","TRBJ","None")
            }else{
              levels <- c("TRAV", "TRBV", "TRAJ","TRBJ")
            }
          }else if(V.or.J == "V"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("TRAV", "TRBV","None")
            }else{
              levels <- c("TRAV", "TRBV")
            }
          }else if(V.or.J == "J"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("TRAJ","TRBJ","None")
            }else{
              levels <- c("TRAJ", "TRBJ")
            }
          }
        }else if(B.or.Tcells == "B"){
          if(V.or.J == "both"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("IGHV","IGLV","IGKV","IGHJ","IGLJ","IGKJ","None")
            }else{
              levels <- c("IGHV","IGLV","IGKV","IGHJ","IGLJ","IGKJ")
            }
          }else if(V.or.J == "V"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("IGHV","IGLV","IGKV","None")
            }else{
              levels <- c("IGHV","IGLV","IGKV")
            }
          }else if(V.or.J == "J"){
            if(ifelse(is.na((colSums(Vgene_usage_matrix[[i]])["None"]>0)),FALSE,(colSums(Vgene_usage_matrix[[i]])["None"]>0))|ifelse(is.na(rowSums(Vgene_usage_matrix[[i]])["None"]>0), FALSE, (rowSums(Vgene_usage_matrix[[i]])["None"]>0))){
              levels <- c("IGHJ","IGLJ","IGKJ","None")
            }else{
              levels <- c("IGHJ","IGLJ","IGKJ")
            }
          }
        }else{
          print("Please specify whether B or T cells are analyzed. (Parameter B.or.Tcells)")
        }
    
        # Set grouping factors based on previously defined levels. Order of levels defines order of groups in Circos plot.
        group = factor(group[sample(length(group), length(group))], levels = levels)
    
        plots[[i]] <- VDJ_circos(Vgene_usage_matrix[[i]], group = group, grid.col=grid.col, label.threshold = label.threshold, c.count = c.count)
      }
      plots[[i+1]] <- Vgene_usage_matrix
}else{
  print("Please specify platypus platypus.version as either v2 or v3.")
  }
  
  return(plots)
}


