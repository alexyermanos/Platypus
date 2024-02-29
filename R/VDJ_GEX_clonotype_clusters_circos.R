#'Makes a Circos plot from the VDJ_GEX_integrate output. Connects the clonotypes with the corresponding clusters.
#'

#' @param VGM The output of the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]]) has to be supplied. For Platypus v2: The output of the VDJ_GEX_integrate function (Platypus platypus.version v2). A list of data frames for each sample containing the clonotype information and cluster membership information.
#' @param n_cluster Integer. No default.
#' @param label.threshold Genes are only labeled if the count is larger then the label.threshold. By default all label.threshold = 0 (all genes are labeled).
#' @param c.count.label Boolean, lets the user decide if the gene and count labels should be plotted or not. Default = T.
#' @param c.count.label.size Determines the font size of the gene labels. By default the font size for count labels is 0.6.
#' @param gene.label Boolean, lets the user decide if the gene labels should be plotted or not.
#' @param gene.label.size Determines the font size of the gene labels. By default the labelsize is automatically adjusted to 0.7 for labels with two or less digits, 0.6 for labels between 2 and 6 digits, and 0.4 for all longer labels. A manually defined font size will be the same for all labels!
#' @param platy.theme Allows plotting in the new "pretty" theme or the older "spiky" theme without group labels and radial arrangement of gene.labels. Default = "pretty".
#' @param arr.col Data.frame with three columns where the first two indicate the names of genes, clonotypes or clusters to be connected, and the third corresponds to the color of the arrow. Default set to data.frame(c("dummy.clonotype"), c("dummy.cluster"), c("dummy.color")), so no arrow is drawn.
#' @param arr.direction Either 1 or -1 and determines the direction of the arrow. Default=1.
#' @param c.threshold Only clonotypes are considered with a frequency higher then c.threshold. Allows to filter for only highly expanded clonotypes.
#' @param topX Filters for the top X clonotypes and only plots the respective gene combinations or cluster memberships.
#' @param clonotype.column Which column in VGM contains the clonotyping information? Default="clonotype_id_10X".
#' @param axis Character. Axis scaling. Defaults to "max". Passed to VDJ_circos
#' @param platypus.version Input version to use. Defaults to "v3" for VDJ_GEX_matrix input
#' @return Returns a circos plot and a list object with the following elememts for N samples: [[1 to N]] The first N listelements corresponds to the recorded circos plots for N beeing the number or samples in the VGM. Since Circlize uses the R base plotting funciton, this is not a ggplot object but can still be replotted by calling the first list element. [[N+1]] Adjacency matrix forwarded to VDJ_circos(). This Matrix contains the counts and can be used for manual replotting using VDJ_circos directly. [[N+2]] Contains a named list with colors for each connection drawn and can be used for manual replotting using VDJ_circos directly. [[N+3]] Contains a named list with grouping information and can be used for manual replotting using VDJ_circos directly.
#' @export
#' @examples
#' \donttest{
#' try({
#'  clonotype.clusters <- VDJ_GEX_clonotype_clusters_circos(Platypus::small_vgm[[1]],
#'  n_cluster=8, topX = 20)
#'  clonotype.clusters[[1]]
#'  })
#'}
#'

VDJ_GEX_clonotype_clusters_circos <- function(VGM,
                                          topX,
                                          label.threshold,
                                          axis,
                                          c.threshold,
                                          c.count.label,
                                          c.count.label.size,
                                          n_cluster,
                                          platypus.version,
                                          gene.label,
                                          gene.label.size,
                                          arr.col,
                                          arr.direction,
                                          platy.theme,
                                          clonotype.column){

  if(missing(topX)){topX <- "all"}
  #if(missing(n_cluster)){stop("Please specify cluster number n_cluster")}
  if(missing(label.threshold)){label.threshold <- 0}
  if(missing(axis)){axis <- "max"}
  if(missing(c.count.label)){c.count.label <- TRUE}
  if(missing(c.count.label.size)){c.count.label.size <- 0.6}
  if(missing(platypus.version)){platypus.version <- "v3"}
  if(missing(gene.label)){gene.label <- TRUE}
  if(missing(gene.label.size)){gene.label.size <- "undef"}
  if(missing(arr.col)){arr.col <- data.frame(c("dummy1"), c("dummy2"), c(""))}
  if(missing(arr.direction)){arr.direction <- 1}
  if(missing(c.threshold)){c.threshold <- 0}
  if(missing(platy.theme)){platy.theme <- "pretty"}
  if(missing(clonotype.column)){clonotype.column <- "clonotype_id_10x"}

  #define Variables

  VDJ.GEX.matrix <- NULL
  VDJ.GEX_list <- NULL
  adj.matrix <- NULL
  clonotypes <- NULL
  clonotypes_all <- NULL
  cluster_col <- NULL
  clonotypes_col <- NULL
  grid.col <- NULL
  plot <- NULL
  nm <- NULL
  group <- NULL
  circos.recored <- NULL
  clonotype.freuqency <- NULL


  #naming compatibility
  VDJ.GEX.matrix <- list()
  VDJ.GEX.matrix[[1]] <- VGM
  VGM <- NULL

  if(platypus.version == "v3"){

      adj.matrix <- list()
      clonotypes <- c()

      message(paste("Chosen clonotype column: ", clonotype.column))
      clonotype.frequency <- paste0("clonotype_frequency_", stringr::str_split(clonotype.column, pattern="_")[[1]][3])
      message(paste("Chosen clonotype.frequency column: ", clonotype.frequency))


      #fill empty entries
      VDJ.GEX.matrix[[1]][[clonotype.column]][which(VDJ.GEX.matrix[[1]][[clonotype.column]] == "")] <- "None"


      #filter out clonotypes with less then c.threshold cells
      if(clonotype.frequency %in% colnames(VDJ.GEX.matrix[[1]])){ #check if 10x frequency already exists... will only be created after reclonotyping.
        VDJ.GEX.matrix[[1]] <-VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]][[clonotype.frequency]] >= c.threshold),]
        message(paste("Chosen clonotype.frequency column: ", clonotype.frequency))
      }else{
        VDJ.GEX.matrix[[1]] <-VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]][["clonotype_frequency"]] >= c.threshold),]
        message(paste("Chosen clonotype.frequency column: ", "clonotype_frequency"))
      }

      #split VDJ.GEX.matrix into samples
      VDJ.GEX_list <- list()

      for (i in 1:length(unique(VDJ.GEX.matrix[[1]]$sample_id))){
        VDJ.GEX_list[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==unique(VDJ.GEX.matrix[[1]]$sample_id)[i]),]
      }

      # MAKE Clonotype df
      clonotypes <- list()
      if(topX != "all"){
        for(k in 1:length(VDJ.GEX_list)){
          clonotypes[[k]] <- names(utils::head(sort(table(VDJ.GEX_list[[1]][[clonotype.column]]),decreasing = TRUE),topX))
        }
      }else{
        for (k in 1:length(VDJ.GEX_list)){
          topX <- length(table(VDJ.GEX_list[[k]][[clonotype.column]]))
          clonotypes[[k]] <- names(utils::head(sort(table(VDJ.GEX_list[[1]][[clonotype.column]]),decreasing = TRUE),topX))
        }
      }
      #filter and keep only cells of topX clonotypes
      for (k in 1:length(VDJ.GEX_list)){
        VDJ.GEX_list[[k]] <- VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]][[clonotype.column]] %in% clonotypes[[k]]),]
      }
      clonotypes_all <- NA

      for (k in 1:length(VDJ.GEX_list)){
        n_cluster <- length(table(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]][[clonotype.column]] %in% clonotypes[[k]]),"seurat_clusters"]))
        adj.matrix[[k]] <- matrix(nrow =length(clonotypes[[k]]), ncol = n_cluster)
        rownames(adj.matrix[[k]]) <- clonotypes[[k]]
        colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ")
        clonotypes_all <- append(clonotypes_all, clonotypes[[k]], after= length(clonotypes_all))

        for (i in 1:nrow(adj.matrix[[k]])){

           adj.matrix[[k]][i,] <- table(VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]][[clonotype.column]]==rownames(adj.matrix[[k]])[i]),"seurat_clusters"])

           # adj.matrix[[k]][i,j]<-as.numeric(str_split(VDJ.GEX_list[[k]]$cluster_membership_percent, pattern = ",")[[i]])[j]
            # adj.matrix[[k]][i,j] <- adj.matrix[[k]][i,j]*length(str_split(VDJ.GEX_list[[k]]$cell_index[[i]], pattern=";")[[1]])/100 #Put here #cell_index -> has to be splitted by ; to get length()
            # rownames(adj.matrix[[k]]) <- VDJ.GEX_list[[k]]$clonotype_id
            # clonotypes <- append(clonotypes, VDJ.GEX_list[[k]]$clonotype_id, after= length(clonotypes))
            # colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ") #c("cluster 0", "cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10")
        }
        adj.matrix[[k]][is.nan(adj.matrix[[k]])] = 0
      }

      ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
      }
      cluster_col <- ggplotColours(n=n_cluster)

      clonotypes_all <- unique(clonotypes_all)
      cluster_col <- stats::setNames(cluster_col,colnames(adj.matrix[[1]]))
      clonotypes_col <- stats::setNames(grDevices::rainbow(length(clonotypes_all)),sample(clonotypes_all))
      grid.col <- append(cluster_col, clonotypes_col)

      plot <- list()
      for (i in 1:length(VDJ.GEX_list)){
        nm = unique(unlist(dimnames(adj.matrix[[i]])))
        group = structure(gsub('[[:digit:]]+', '', nm), names = nm)
        group = factor(group[sample(length(group), length(group))], levels = c("cluster ", "clonotype"))

        VDJ_circos(adj.matrix[[i]], group = group, grid.col = grid.col, label.threshold = label.threshold, axis = axis, c.count.label = c.count.label, c.count.label.size = c.count.label.size, gene.label = gene.label, gene.label.size = gene.label.size, arr.col = arr.col, arr.direction = arr.direction, platy.theme = platy.theme)
        circos.recorded <- grDevices::recordPlot()
        plot[[i]] <- circos.recorded
      }

      plot[[i+1]] <- adj.matrix
      plot[[i+2]] <- grid.col
      plot[[i+3]] <- group

  }else if(platypus.version == "v2"){
    ########################

    adj.matrix <- list()
    clonotypes <- c()

    if(topX != "all"){
      for(k in 1:length(VDJ.GEX.matrix)){
        VDJ.GEX.matrix[[k]] <- utils::head(VDJ.GEX.matrix[[k]], topX)
      }
    }
    for (k in 1:length(VDJ.GEX.matrix)){
      n_cluster <-  length(stringr::str_split(VDJ.GEX.matrix[[k]]$cluster_membership_percent, pattern = ",")[[1]])
      adj.matrix[[k]] <- matrix(nrow =nrow(VDJ.GEX.matrix[[k]]), ncol = n_cluster)
      for (i in 1:nrow(VDJ.GEX.matrix[[k]])){
        for (j in 1:n_cluster){
          adj.matrix[[k]][i,j]<-as.numeric(stringr::str_split(VDJ.GEX.matrix[[k]]$cluster_membership_percent, pattern = ",")[[i]])[j]
          adj.matrix[[k]][i,j] <- adj.matrix[[k]][i,j]*length(stringr::str_split(VDJ.GEX.matrix[[k]]$cell_index[[i]], pattern=";")[[1]])/100 #Put here #cell_index -> has to be splitted by ; to get length()
          rownames(adj.matrix[[k]]) <- VDJ.GEX.matrix[[k]]$clonotype_id
          clonotypes <- append(clonotypes, VDJ.GEX.matrix[[k]]$clonotype_id, after= length(clonotypes))
          colnames(adj.matrix[[k]]) <- paste("cluster", 0:(ncol(adj.matrix[[k]])-1), sep = " ") #c("cluster 0", "cluster 1","cluster 2","cluster 3","cluster 4","cluster 5","cluster 6","cluster 7","cluster 8","cluster 9","cluster 10")
        }
      }
      adj.matrix[[k]][is.nan(adj.matrix[[k]])] = 0
    }

    ggplotColours <- function(n = 6, h = c(0, 360) + 15){
      if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
      grDevices::hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
    }
    cluster_col <- ggplotColours(n=n_cluster)

    clonotypes <- unique(clonotypes)
    cluster_col <- stats::setNames(cluster_col,colnames(adj.matrix[[1]]))
    clonotypes_col <- stats::setNames(grDevices::rainbow(length(clonotypes)),sample(clonotypes))
    grid.col <- append(cluster_col, clonotypes_col)

    plot <- list()
    for (i in 1:length(VDJ.GEX.matrix)){
      nm = unique(unlist(dimnames(adj.matrix[[i]])))
      group = structure(gsub('[[:digit:]]+', '', nm), names = nm)
      group = factor(group[sample(length(group), length(group))], levels = c("cluster ", "clonotype"))
      VDJ_circos(adj.matrix[[i]], group = group, grid.col = grid.col, label.threshold = label.threshold, axis = axis, c.count.label=c.count.label, c.count.label.size)
      circos.recorded <- grDevices::recordPlot()
      plot[[i]] <- circos.recorded
    }
    plot[[i+1]] <- adj.matrix
    plot[[i+2]] <- grid.col
    plot[[i+3]] <- group

  }else{
    stop("Please specify platypus.version as either v2 or v3.")
  }

  return(plot)
}
