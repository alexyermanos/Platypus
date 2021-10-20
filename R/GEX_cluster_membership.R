#' Plots the cluster membership for each of the distinct samples in the Seurat object from the automate_GEX function. The distinct samples are determined by "sample_id" field in the Seurat object.
#' @param GEX Output Seurat object containing gene expression data from automate_GEX (platypus.version = "v2") or VDJ_GEX_matrix (platypus.version = "v3", usually VDJ_GEX_matrix.output[[2]])function that contained at least two distinct biological samples. The different biological samples correspond to integer values (v2) or factor values (v3) in the order of the working directories initially supplied to the automate_GEX function.
#' @param by.group Logical indicating whether to look at the cluster distribution per group (using the group_id column). Default is set to FALSE.
#' @param platypus.version Version of platypus to use. Defaults to "v2". If an output of the GEX_automate function is supplied, set to "v2". If an output of the VDJ_GEX_matrix function is supplied set to "v3"
#' @return Returns a ggplot object in which the values on the x axis correspond to each cluster found in the Seurat object. The y axis corresponds to the percentage of cells found in each cluster. The bar and color corresponds to the distinct sample_id.
#' @export
#' @examples
#' #Platypus v2
#' #GEX_cluster_membership(GEX=automate_GEX_out[[2]], platypus.version = "v2")
#' #Platypus v3
#' GEX_cluster_membership(GEX= Platypus::small_vgm[[2]], platypus.version = "v3")

GEX_cluster_membership <- function(GEX,
                                   by.group,
                                   platypus.version){

  Sample <- NULL
  value <- NULL
  L2 <- NULL
  Group <- NULL

  if(missing(by.group)) by.group <- FALSE
  if(missing(platypus.version)) platypus.version <- "v2"

  #platypus version v2
  if(platypus.version == "v2"){

  if(by.group == FALSE){
  unique_samples <- sort(unique(GEX$sample_id),decreasing = F)
  unique_clusters <- sort(unique(GEX$seurat_clusters),decreasing = F)
  cells_per_cluster_per_sample <- list()
  for(i in 1:length(unique_samples)){
    cells_per_cluster_per_sample[[i]] <- list()
    for(j in 1:length(unique_clusters)){
      cells_per_cluster_per_sample[[i]][[j]] <- length(which(GEX$sample_id==unique_samples[i] & GEX$seurat_clusters==unique_clusters[j]))/length(which(GEX$sample_id==unique_samples[i]))
    }
  }
  melting <- reshape2::melt(cells_per_cluster_per_sample)
  melting$L1 <- as.character(melting$L1)
  melting$L2 <- melting$L2 - 1
  colnames(melting) <- c("value", "L2", "Sample")
  output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = Sample, y=value, x=L2,group=Sample)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster membership") + ggplot2::scale_x_continuous(breaks = seq(from = 0, to = length(unique(melting$L2)), by = 1))
  }

  if(by.group == TRUE){
    unique_groups <- sort(unique(GEX$group_id),decreasing = F)
    unique_clusters <- sort(unique(GEX$seurat_clusters),decreasing = F)
    cells_per_cluster_per_group <- list()
    for(i in 1:length(unique_groups)){
      cells_per_cluster_per_group[[i]] <- list()
      for(j in 1:length(unique_clusters)){
        cells_per_cluster_per_group[[i]][[j]] <- length(which(GEX$group_id==unique_groups[i] & GEX$seurat_clusters==unique_clusters[j]))/length(which(GEX$group_id==unique_groups[i]))
      }
    }
    melting <- reshape2::melt(cells_per_cluster_per_group)
    melting$L1 <- as.character(melting$L1)
    melting$L2 <- melting$L2 - 1
    colnames(melting) <- c("value", "L2", "Group")
    output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = Group, y=value, x=L2,group=Group)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster membership") + ggplot2::scale_x_continuous(breaks = seq(from = 0, to = length(unique(melting$L2)), by = 1))
  }

  return(output.plot)

  }
  #Platypus version v3
  else if (platypus.version == "v3"){

    props <- "sample"

    if(by.group == FALSE){
      unique_samples <- unique(GEX$sample_id)
      unique_clusters <- sort(unique(GEX$seurat_clusters),decreasing = F)
      cells_per_cluster_per_sample <- list()
      for(i in 1:length(unique_samples)){
        cells_per_cluster_per_sample[[i]] <- list()
        names(cells_per_cluster_per_sample)[i] <- unique_samples[i]
        for(j in 1:length(unique_clusters)){
          if(props == "sample"){
            cells_per_cluster_per_sample[[i]][[j]] <- length(which(GEX$sample_id==unique_samples[i] & GEX$seurat_clusters==unique_clusters[j]))/length(which(GEX$sample_id==unique_samples[i])) * 100
          } else if(props == "cluster"){
            cells_per_cluster_per_sample[[i]][[j]] <- length(which(GEX$sample_id==unique_samples[i] & GEX$seurat_clusters==unique_clusters[j]))/length(which(GEX$seurat_clusters==unique_clusters[j])) /length(which(GEX$sample_id==unique_samples[i])) * 100
          }
        }
      }
      melting <- reshape2::melt(cells_per_cluster_per_sample)

      melting$L1 <- as.character(melting$L1)
      melting$L2 <- melting$L2 - 1
      colnames(melting) <- c("value", "L2", "Sample")

      #adjust sample names
      for(i in 1:length(unique(melting$Sample))){
        melting$Sample[which(melting$Sample == as.character(i))] <- as.character(unique_samples[i])
      }
      melting$Sample <- ordered(as.factor(melting$Sample), levels = unique_samples)

      output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = Sample, y=value, x=L2,group=Sample)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster membership") + ggplot2::scale_x_continuous(breaks = seq(from = 0, to = length(unique(melting$L2)), by = 1))
    }

    if(by.group == TRUE){
      unique_groups <- unique(GEX$group_id)
      unique_clusters <- sort(unique(GEX$seurat_clusters),decreasing = F)
      cells_per_cluster_per_group <- list()
      for(i in 1:length(unique_groups)){
        cells_per_cluster_per_group[[i]] <- list()
        names(cells_per_cluster_per_group)[i] <- unique_groups[i]
        for(j in 1:length(unique_clusters)){
          cells_per_cluster_per_group[[i]][[j]] <- length(which(GEX$group_id==unique_groups[i] & GEX$seurat_clusters==unique_clusters[j]))/length(which(GEX$group_id==unique_groups[i]))
        }
      }
      melting <- reshape2::melt(cells_per_cluster_per_group)
      melting$L1 <- as.character(melting$L1)
      melting$L2 <- melting$L2 - 1
      colnames(melting) <- c("value", "L2", "Group")

      #adjust group names
      for(i in 1:length(unique(melting$Group))){
        melting$Group[which(melting$Group == as.character(i))] <- as.character(unique_groups[i])
      }
      melting$Group <- ordered(as.factor(melting$Group), levels = unique_groups)

      output.plot <- ggplot2::ggplot(melting, ggplot2::aes(fill = Group, y=value, x=L2,group=Group)) + ggplot2::geom_bar(stat="identity", width=0.6, color="black",position = "dodge") + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::scale_y_continuous(expand = c(0,0)) + ggplot2::ylab("Fraction of cells") + ggplot2::xlab("Cluster membership") + ggplot2::scale_x_continuous(breaks = seq(from = 0, to = length(unique(melting$L2)), by = 1))
    }

    return(output.plot)
  }

}
