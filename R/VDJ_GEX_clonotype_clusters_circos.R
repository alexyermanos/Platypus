#'Makes a Circos plot from the VDJ_GEX_integrate output. Connects the clonotypes with the corresponding clusters.
#' @param VDJ.GEX.matrix The output of the VDJ_GEX_integrate function (Platypus platypus.version v2). A list of data frames for each sample containing the clonotype information and cluster membership information. For Platypus platypus.version v3, the VDJ_GEX_matrix() output has to be supplied.
#' @param topX Plots only the top X most expanded clonotypes. By default all clonotypes are shown.
#' @param label.threshold Minimal amount of clonotypes per gene neccessary to add a gene label to the sector. Default: 0.
#' @param axis Character. Defaults to "max". Passed to VDJ_circos
#' @param n_cluster Integer. No default.
#' @param c.count Show clonotype or cell count on Circos plot. Default = T.
#' @param platypus.version Which platypus.version of platypus is beeing used. Default = v2.
#' @return Returns list of plots. The first n elements contain the circos plot of the n datasets from the VDJ.analyze function. The n+1 element contains a list of the n adjancey matrices for each dataset.
#' @examples
#' \dontrun{
#'  plots <- VDJ_clonotype_clusters_circos(vdj_gex_integrate_test, topX=100, label.threshold=5)
#'}

VDJ_clonotype_clusters_circos <- function(VDJ.GEX.matrix, topX, label.threshold, axis, c.count, n_cluster, platypus.version){
  if(missing(topX)){topX <- "all"}
  if(missing(n_cluster)){print("Please specify cluster number n_cluster")}
  if(missing(label.threshold)){label.threshold <- 1}
  if(missing(axis)){axis <- "max"}
  if(missing(c.count)){c.count <-T}
  if(missing(platypus.version)){platypus.version <- "v2"}


  if(platypus.version == "v3"){
      print("Reminder: VDJ_VJ_usage_circos() funcion built for new Platypus v3.0.0 is being used. Output of VDJ_GEX_matrix() required as input.")
      adj.matrix <- list()
      clonotypes <- c()

      #fill empty entries
      VDJ.GEX.matrix[[1]]$clonotype_id_10x[which(VDJ.GEX.matrix[[1]]$clonotype_id_10x == "")] <- "None"

      #split VDJ.GEX.matrix into samples
      VDJ.GEX_list <- list()
      for (i in 1:length(table(VDJ.GEX.matrix[[1]]$sample_id))){
        VDJ.GEX_list[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==paste0("s",i)),]
      }

      # MAKE Clonotype df
      clonotypes <- list()
      if(topX != "all"){
        for(k in 1:length(VDJ.GEX_list)){
          clonotypes[[k]] <- names(head(sort(table(VDJ.GEX_list[[1]]$clonotype_id_10x),decreasing = T),topX))
        }
      }else{
        for (k in 1:length(VDJ.GEX_list)){
          topX <- length(table(VDJ.GEX_list[[k]]$clonotype_id_10x))
          clonotypes[[k]] <- names(head(sort(table(VDJ.GEX_list[[1]]$clonotype_id_10x),decreasing = T),topX))
        }
      }
      #filter and keep only cells of topX clonotypes
      for (k in 1:length(VDJ.GEX_list)){
        VDJ.GEX_list[[k]] <- VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]]$clonotype_id_10x %in% clonotypes[[k]]),]
      }
      clonotypes_all <- NA

      for (k in 1:length(VDJ.GEX_list)){
        print(paste0("Processing Sample ", k))
        n_cluster <- length(table(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]]$clonotype_id_10x %in% clonotypes[[k]]),"seurat_clusters"]))
        adj.matrix[[k]] <- matrix(nrow =length(clonotypes[[k]]), ncol = n_cluster)
        rownames(adj.matrix[[k]]) <- clonotypes[[k]]
        colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ")
        clonotypes_all <- append(clonotypes_all, clonotypes[[k]], after= length(clonotypes_all))

        for (i in 1:nrow(adj.matrix[[k]])){

           adj.matrix[[k]][i,] <- table(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]]$clonotype_id_10x==rownames(adj.matrix[[k]])[i]),"seurat_clusters"])

           # adj.matrix[[k]][i,j]<-as.numeric(str_split(VDJ.GEX_list[[k]]$cluster_membership_percent, pattern = ",")[[i]])[j]
            # adj.matrix[[k]][i,j] <- adj.matrix[[k]][i,j]*length(str_split(VDJ.GEX_list[[k]]$cell_index[[i]], pattern=";")[[1]])/100 #Put here #cell_index -> has to be splitted by ; to get length()
            # rownames(adj.matrix[[k]]) <- VDJ.GEX_list[[k]]$clonotype_id
            # clonotypes <- append(clonotypes, VDJ.GEX_list[[k]]$clonotype_id, after= length(clonotypes))
            # colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ") #c("cluster 0", "cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10")
        }
        adj.matrix[[k]][is.nan(adj.matrix[[k]])] = 0
      }


      print("Set colours")

      ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
      }
      cluster_col <- ggplotColours(n=n_cluster)

      clonotypes_all <- unique(clonotypes_all)
      cluster_col <- setNames(cluster_col,colnames(adj.matrix[[1]]))
      clonotypes_col <- setNames(rainbow(length(clonotypes_all)),sample(clonotypes_all))
      grid.col <- append(cluster_col, clonotypes_col)

      print("Plotting")
      plot <- list()
      for (i in 1:length(VDJ.GEX_list)){
        print(i)
        print(length(VDJ.GEX_list))
        print(adj.matrix[[i]])
        nm = unique(unlist(dimnames(adj.matrix[[i]])))
        group = structure(gsub('[[:digit:]]+', '', nm), names = nm)
        group = factor(group[sample(length(group), length(group))], levels = c("cluster ", "clonotype"))
        plot[[i]] <- VDJ_circos(adj.matrix[[i]], group = group, grid.col = grid.col, label.threshold = label.threshold, axis = axis, c.count=c.count)
      }

      plot[[i+1]] <- adj.matrix

  }else if(platypus.version == "v2"){
    ########################
    print("Reminder: VDJ_VJ_usage_circos() funcion built for Platypus v2.0.0 is being used. Output of VDJ_analyze() required as input. Set [platypus.version = new] for compatibility with VDJ_GEX_matrix().")
    adj.matrix <- list()
    clonotypes <- c()

    if(topX != "all"){
      for(k in 1:length(VDJ.GEX.matrix)){
        VDJ.GEX.matrix[[k]] <- head(VDJ.GEX.matrix[[k]], topX)
      }
    }
    for (k in 1:length(VDJ.GEX.matrix)){
      print(k)
      n_cluster <-  length(str_split(VDJ.GEX.matrix[[k]]$cluster_membership_percent, pattern = ",")[[1]])
      adj.matrix[[k]] <- matrix(nrow =nrow(VDJ.GEX.matrix[[k]]), ncol = n_cluster)
      for (i in 1:nrow(VDJ.GEX.matrix[[k]])){
        for (j in 1:n_cluster){
          adj.matrix[[k]][i,j]<-as.numeric(str_split(VDJ.GEX.matrix[[k]]$cluster_membership_percent, pattern = ",")[[i]])[j]
          adj.matrix[[k]][i,j] <- adj.matrix[[k]][i,j]*length(str_split(VDJ.GEX.matrix[[k]]$cell_index[[i]], pattern=";")[[1]])/100 #Put here #cell_index -> has to be splitted by ; to get length()
          rownames(adj.matrix[[k]]) <- VDJ.GEX.matrix[[k]]$clonotype_id
          clonotypes <- append(clonotypes, VDJ.GEX.matrix[[k]]$clonotype_id, after= length(clonotypes))
          colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ") #c("cluster 0", "cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10")
        }
      }
      adj.matrix[[k]][is.nan(adj.matrix[[k]])] = 0
    }
    print("Set colours")

    ggplotColours <- function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    cluster_col <- ggplotColours(n=n_cluster)

    clonotypes <- unique(clonotypes)
    cluster_col <- setNames(cluster_col,colnames(adj.matrix[[1]]))
    clonotypes_col <- setNames(rainbow(length(clonotypes)),sample(clonotypes))
    grid.col <- append(cluster_col, clonotypes_col)

    print("Plotting")
    plot <- list()
    for (i in 1:length(VDJ.GEX.matrix)){
      print(i)
      print(length(VDJ.GEX.matrix))
      print(adj.matrix[[i]])
      nm = unique(unlist(dimnames(adj.matrix[[i]])))
      group = structure(gsub('[[:digit:]]+', '', nm), names = nm)
      group = factor(group[sample(length(group), length(group))], levels = c("cluster ", "clonotype"))
      plot[[i]] <- VDJ_circos(adj.matrix[[i]], group = group, grid.col = grid.col, label.threshold = label.threshold, axis = axis, c.count=c.count)
    }
    plot[[i+1]] <- adj.matrix
  }else{
    print("Please specify platypus platypus.version as either v2 or v3.")
  }

  return(plot)
}
