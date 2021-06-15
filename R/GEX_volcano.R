#' Plots a volcano plot from the output of the FindMarkers function from the Seurat package or the GEX_cluster_genes function alternatively.
#' @param findmarkers.output output data frame from the FindMarkers function from the Seurat package or GEX_cluster_genes list output.
#' @param cluster.genes.output logical specifying wheter GEX_cluster_genes output is being used.
#' @param condition.1 either character or integer specifying ident.1 that was used in the FindMarkers function from the Seurat package. Should be left empty when using the GEX_cluster_genes output.
#' @param condition.2 either character or integer specifying ident.2 that was used in the FindMarkers function from the Seurat package. Should be left empty when using the GEX_cluster_genes output.
#' @param explicit.title logical specifying whether the title should include logFC information for each condition.
#' @param color.p.threshold numeric specifying the adjusted p-value threshold for geom_points to be colored. Default is set to 0.01.
#' @param color.log.threshold numeric specifying the absolute logFC threshold for geom_points to be colored. Default is set to 0.25.
#' @param label.p.threshold numeric specifying the adjusted p-value threshold for genes to be labeled via geom_text_repel. Default is set to 0.001.
#' @param label.logfc.threshold numeric specifying the absolute logFC threshold for genes to be labeled via geom_text_repel. Default is set to 0.75.
#' @param n.label.up numeric specifying the number of top upregulated genes to be labeled via geom_text_repel. Genes will be ordered by adjusted p-value. Overrides the "label.p.threshold" and "label.logfc.threshold" parameters.
#' @param n.label.down numeric specifying the number of top downregulated genes to be labeled via geom_text_repel. Genes will be ordered by adjusted p-value. Overrides the "label.p.threshold" and "label.logfc.threshold" parameters.
#' @param maximum.overlaps integer specifying removal of labels with too many overlaps. Default is set to Inf.
#' @param plot.adj.pvalue logical specifying whether adjusted p-value should by plotted on the y-axis. 
#' @return Returns a volcano plot from the output of the FindMarkers function from the Seurat package, which is a ggplot object that can be modified or plotted. Infinite p-values are set defined value of the highest -log(p) + 100.
#' @export
#' @examples
#' \dontrun{
#' #using the findmarkers.output
#' GEX_volcano_1 <- GEX_volcano(findmarkers.output = FindMarkers.Output, condition.1 = "cluster1", condition.2 = "cluster2", maximum.overlaps = 20)
#' GEX_volcano_2 <- GEX_volcano(findmarkers.output = FindMarkers.Output, condition.1 = "cluster1", condition.2 = "cluster2", n.label.up = 50, n.label.down = 20)
#' #using the GEX_cluster_genes output
#' GEX_volcano_3 <- GEX_volcano(findmarkers.output = GEX_cluster_genes.Output, cluster.genes.output =T)
#'}
GEX_volcano <- function(findmarkers.output, cluster.genes.output, condition.1, condition.2, explicit.title, RP.MT.filter, color.p.threshold, color.log.threshold, label.p.threshold, label.logfc.threshold, n.label.up, n.label.down ,maximum.overlaps ,plot.adj.pvalue) {
  require(ggrepel)
  
  if(missing(cluster.genes.output)){cluster.genes.output <- F}
  if(missing(condition.1)){condition.1 <- ""}
  if(missing(condition.2)){condition.2 <- ""}
  if(missing(explicit.title)){explicit.title <- T}
  if(missing(RP.MT.filter)){RP.MT.filter <- T}
  if(missing(color.p.threshold)){color.p.threshold <- 0.01}
  if(missing(color.log.threshold)){color.log.threshold <- 0.25}
  if(missing(label.p.threshold)){label.p.threshold <- 0.001}
  if(missing(label.logfc.threshold)){label.logfc.threshold <- 0.75}
  if(missing(n.label.up)){n.label.up <- F}
  if(missing(n.label.down)){n.label.down <- F}
  if(missing(maximum.overlaps)) {maximum.overlaps <- Inf}
  if(missing(plot.adj.pvalue)){plot.adj.pvalue <- F}
  
  if(class(n.label.up) != class(n.label.down)){
    temp <- list(n.label.up, n.label.down)
    n.genes <- as.numeric(Filter(is.numeric, temp))
    n.label.up <- n.genes
    n.label.down <- n.genes
    print("Same number of up- and downregulated genes will be labeled")
  } # Assuming the same number of up- and downregulated genes are to be labeled if only one is specified
  
  if(cluster.genes.output == F) {
    
    if(RP.MT.filter ==T){
      exclude <- c()
      for (i in c("MT-", "RPL", "RPS")) {
        exclude <- c(exclude, stringr::str_which(rownames(findmarkers.output), i))
      }
      findmarkers.output <- findmarkers.output[-exclude,]
    }
    
    findmarkers.output$genes <- rownames(findmarkers.output)
    
    if(plot.adj.pvalue == F) {
      findmarkers.output$minus.log10p <- -log10(findmarkers.output$p_val) 
      findmarkers.output$minus.log10p[which(findmarkers.output$minus.log10p == Inf)] <- findmarkers.output$minus.log10p[which(sort(findmarkers.output$minus.log10p, decreasing = TRUE) != Inf)[1]]+100 #calculating -log10p and setting the ones that are Inf to defined value
      
      if(n.label.up == F & n.label.down == F) {
        output.plot <- ggplot2::ggplot(findmarkers.output, aes(x=avg_log2FC, y=minus.log10p, label = genes)) + ggplot2::geom_point() +
          ggplot2::geom_point(data = subset(findmarkers.output, abs(avg_log2FC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
          ggrepel::geom_text_repel(data =subset(findmarkers.output, abs(avg_log2FC) > label.logfc.threshold & p_val_adj < label.p.threshold),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
          ggplot2::theme_bw() + ggplot2::ylab("-log10(p-value)")
      }
      
      if(class(n.label.up) == "numeric" & class(n.label.down) == "numeric") {
        findmarkers.output <- findmarkers.output[order(findmarkers.output$p_val_adj),]
        posFC_genes <- findmarkers.output$genes[which(findmarkers.output$avg_log2FC > 0)][1:n.label.up]
        negFC_genes <- findmarkers.output$genes[which(findmarkers.output$avg_log2FC < 0)][1:n.label.down]
        label.genes <- c(posFC_genes,negFC_genes)
        
        output.plot <- ggplot2::ggplot(findmarkers.output, aes(x=avg_log2FC, y=minus.log10p, label = genes)) + ggplot2::geom_point() +
          ggplot2::geom_point(data = subset(findmarkers.output, abs(avg_log2FC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
          ggrepel::geom_text_repel(data =subset(findmarkers.output, genes%in%label.genes),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
          ggplot2::theme_bw() + ggplot2::ylab("-log10(p-value)")
      }
      
    }
    if(plot.adj.pvalue == T) {
      findmarkers.output$minus.log10p_adj <- -log10(findmarkers.output$p_val_adj) 
      findmarkers.output$minus.log10p_adj[which(findmarkers.output$minus.log10p_adj == Inf)] <- findmarkers.output$minus.log10p_adj[which(sort(findmarkers.output$minus.log10p_adj, decreasing = TRUE) != Inf)[1]]+100 #calculating -log10p and setting the ones that are Inf to defined value
      
      if(n.label.up == F & n.label.down == F) {
        output.plot <- ggplot2::ggplot(findmarkers.output, aes(x=avg_log2FC, y=minus.log10p_adj, label = genes)) + ggplot2::geom_point() +
          ggplot2::geom_point(data = subset(findmarkers.output, abs(avg_log2FC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
          ggrepel::geom_text_repel(data =subset(findmarkers.output, abs(avg_log2FC) > label.logfc.threshold & p_val_adj < label.p.threshold),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
          ggplot2::theme_bw() + ggplot2::ylab("-log10(adj.p-value)") 
      }
      
      if(class(n.label.up) == "numeric" & class(n.label.down) == "numeric") {
        findmarkers.output <- findmarkers.output[order(findmarkers.output$p_val_adj),]
        posFC_genes <- findmarkers.output$genes[which(findmarkers.output$avg_log2FC > 0)][1:n.label.up]
        negFC_genes <- findmarkers.output$genes[which(findmarkers.output$avg_log2FC < 0)][1:n.label.down]
        label.genes <- c(posFC_genes,negFC_genes)
        
        output.plot <- ggplot2::ggplot(findmarkers.output, aes(x=avg_log2FC, y=minus.log10p_adj, label = genes)) + ggplot2::geom_point() +
          ggplot2::geom_point(data = subset(findmarkers.output, abs(avg_log2FC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
          ggrepel::geom_text_repel(data =subset(findmarkers.output, genes%in%label.genes),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
          ggplot2::theme_bw() + ggplot2::ylab("-log10(adj.p-value)") 
      }
    }
    
    if(explicit.title == F) {
      output.plot <- output.plot + ggplot2::ggtitle(paste(condition.2, "vs.", condition.1))
    }
    if(explicit.title ==T){
      output.plot <- output.plot + ggplot2::ggtitle(paste(condition.2 , "(up=-FC)", "vs.", condition.1, "(up=+FC)"))
    }
  }
  
  if(cluster.genes.output == T){
    output.plot <- list()
    for (i in 1:length(findmarkers.output)) {
      
      if(RP.MT.filter ==T){
        exclude <- c()
        for (j in c("MT-", "RPL", "RPS")) {
          exclude <- c(exclude, stringr::str_which(rownames(findmarkers.output[[i]]), j))
        }
        if(length(exclude) != 0){
        findmarkers.output[[i]] <- findmarkers.output[[i]][-exclude,]
        }
      }
      
      if(plot.adj.pvalue == F) {
        findmarkers.output[[i]]$minus.log10p <- -log10(findmarkers.output[[i]]$p_val) 
        findmarkers.output[[i]]$minus.log10p[which(findmarkers.output[[i]]$minus.log10p == Inf)] <- findmarkers.output[[i]]$minus.log10p[which(sort(findmarkers.output[[i]]$minus.log10p, decreasing = TRUE) != Inf)[1]]+100 #calculating -log10p and setting the ones that are Inf to defined value
        
        if(n.label.up == F & n.label.down == F) {
          cluster_plot <- ggplot2::ggplot(findmarkers.output[[i]], aes(x=avg_logFC, y=minus.log10p, label = SYMBOL)) + ggplot2::geom_point() +
            ggplot2::geom_point(data = subset(findmarkers.output[[i]], abs(avg_logFC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
            ggrepel::geom_text_repel(data =subset(findmarkers.output[[i]], abs(avg_logFC) > label.logfc.threshold & p_val_adj < label.p.threshold),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
            ggplot2::theme_bw() + ggplot2::ylab("-log10(p-value)")
        }
        
        if(class(n.label.up) == "numeric" & class(n.label.down) == "numeric") {
          findmarkers.output[[i]] <- findmarkers.output[[i]][order(findmarkers.output[[i]]$p_val_adj),]
          posFC_genes <- findmarkers.output[[i]]$SYMBOL[which(findmarkers.output[[i]]$avg_logFC > 0)][1:n.label.up]
          negFC_genes <- findmarkers.output[[i]]$SYMBOL[which(findmarkers.output[[i]]$avg_logFC < 0)][1:n.label.down]
          label.genes <- c(posFC_genes,negFC_genes)
          
          cluster_plot <- ggplot2::ggplot(findmarkers.output[[i]], aes(x=avg_logFC, y=minus.log10p, label = SYMBOL)) + ggplot2::geom_point() +
            ggplot2::geom_point(data = subset(findmarkers.output[[i]], abs(avg_logFC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
            ggrepel::geom_text_repel(data =subset(findmarkers.output[[i]], SYMBOL%in%label.genes),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
            ggplot2::theme_bw() + ggplot2::ylab("-log10(p-value)")
        }
      }
      
      if(plot.adj.pvalue == T) {
        findmarkers.output[[i]]$minus.log10p_adj <- -log10(findmarkers.output[[i]]$p_val_adj) 
        findmarkers.output[[i]]$minus.log10p_adj[which(findmarkers.output[[i]]$minus.log10p_adj == Inf)] <- findmarkers.output[[i]]$minus.log10p_adj[which(sort(findmarkers.output[[i]]$minus.log10p_adj, decreasing = TRUE) != Inf)[1]]+100 #calculating -log10p and setting the ones that are Inf to defined value
        
        if(n.label.up == F & n.label.down == F) {
          cluster_plot <- ggplot2::ggplot(findmarkers.output[[i]], aes(x=avg_logFC, y=minus.log10p_adj, label = SYMBOL)) + ggplot2::geom_point() +
            ggplot2::geom_point(data = subset(findmarkers.output[[i]], abs(avg_logFC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
            ggrepel::geom_text_repel(data =subset(findmarkers.output[[i]], abs(avg_logFC) > label.logfc.threshold & p_val_adj < label.p.threshold),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
            ggplot2::theme_bw() + ggplot2::ylab("-log10(adj.p-value)") 
        }
        
        if(class(n.label.up) == "numeric" & class(n.label.down) == "numeric") {
          findmarkers.output[[i]] <- findmarkers.output[[i]][order(findmarkers.output[[i]]$p_val_adj),]
          posFC_genes <- findmarkers.output[[i]]$SYMBOL[which(findmarkers.output[[i]]$avg_logFC > 0)][1:n.label.up]
          negFC_genes <- findmarkers.output[[i]]$SYMBOL[which(findmarkers.output[[i]]$avg_logFC < 0)][1:n.label.down]
          label.genes <- c(posFC_genes,negFC_genes)
          
          cluster_plot <- ggplot2::ggplot(findmarkers.output[[i]], aes(x=avg_logFC, y=minus.log10p_adj, label = SYMBOL)) + ggplot2::geom_point() +
            ggplot2::geom_point(data = subset(findmarkers.output[[i]], abs(avg_logFC) > color.log.threshold & p_val_adj < color.p.threshold), col= "darkred") + 
            ggrepel::geom_text_repel(data =subset(findmarkers.output[[i]], SYMBOL%in%label.genes),  color = 'black', hjust = 0, direction = "y", max.overlaps = maximum.overlaps) +
            ggplot2::theme_bw() + ggplot2::ylab("-log10(adj.p-value)") 
        }
      }
      
      if(explicit.title == F) {
        output.plot[[i]] <- cluster_plot + ggplot2::ggtitle(paste("All clusters vs.", "cluster", i-1))
      }
      if(explicit.title ==T){
        output.plot[[i]] <- cluster_plot + ggplot2::ggtitle(paste("All clusters (up=-FC)", "vs.", "cluster", i-1,"(up=+FC)"))
      } 
    }
  }
  return(output.plot)
}