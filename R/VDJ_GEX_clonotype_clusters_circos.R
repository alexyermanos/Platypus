#'Makes a Circos plot from the VDJ_GEX_integrate output. Connects the clonotypes with the corresponding clusters.
#' @param vdj.gex.integrate.output The output of the VDJ_GEX_integrate function. A list of data frames for each sample containing the clonotype information and cluster membership information.
#' @param TopX Plots only the top X most expanded clonotypes. By default all clonotypes are shown.
#' @param label.threshold Minimal amount of clonotypes per gene neccessary to add a gene label to the sector. Default: 0.
#' @param c.count Show clonotype or cell count on Circos plot. Default = T.
#' @return Returns list of plots. The first n elements contain the circos plot of the n datasets from the VDJ.analyze function. The n+1 element contains a list of the n adjancey matrices for each dataset.
#' @examples
#' \dontrun{
#'  plots <- VDJ_clonotype_clusters_circos(vdj_gex_integrate_test, topX=100, label.threshold=5)
#'}

VDJ_clonotype_clusters_circos <- function(vdj.gex.integrate.output, topX, label.threshold, axis, c.count, n_cluster){
  if(missing(topX)){topX <- "all"}
  if(missing(n_cluster)){print("Please specify cluster number n_cluster")}
  if(missing(label.threshold)){label.threshold <- 1}
  if(topX != "all"){
    for(k in 1:length(vdj.gex.integrate.output)){
      vdj.gex.integrate.output[[k]] <- head(vdj.gex.integrate.output[[k]], topX)
    }
  }
  if(missing(axis)){axis <- "max"}
  if(missing(c.count)){c.count <-T}

  adj.matrix <- list()
  clonotypes <- c()
  for (k in 1:length(vdj.gex.integrate.output)){
    print(k)
    adj.matrix[[k]] <- matrix(nrow =nrow(vdj.gex.integrate.output[[k]]), ncol = n_cluster)
    for (i in 1:nrow(vdj.gex.integrate.output[[k]])){
      for (j in 1:n_cluster){
        adj.matrix[[k]][i,j]<-as.numeric(str_split(vdj.gex.integrate.output[[k]]$cluster_membership_percent, pattern = ",")[[i]])[j]
        adj.matrix[[k]][i,j] <- adj.matrix[[k]][i,j]*length(str_split(vdj.gex.integrate.output[[k]]$cell_index[[i]], pattern=";")[[1]])/100 #Put here #cell_index -> has to be splitted by ; to get length()
        rownames(adj.matrix[[k]]) <- vdj.gex.integrate.output[[k]]$clonotype_id
        clonotypes <- append(clonotypes, vdj.gex.integrate.output[[k]]$clonotype_id, after= length(clonotypes))
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
  for (i in 1:length(vdj.gex.integrate.output)){
    print(i)
    print(length(vdj.gex.integrate.output))
    print(adj.matrix[[i]])
    nm = unique(unlist(dimnames(adj.matrix[[i]])))
    group = structure(gsub('[[:digit:]]+', '', nm), names = nm)
    group = factor(group[sample(length(group), length(group))], levels = c("cluster ", "clonotype"))
    plot[[i]] <- VDJ_circos(adj.matrix[[i]], group = group, grid.col = grid.col, label.threshold = label.threshold, axis = axis, c.count=c.count)
  }
  plot[[i+1]] <- adj.matrix
  return(plot)
}
