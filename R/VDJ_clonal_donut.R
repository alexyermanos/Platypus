#' Generate circular plots of clonal expansion per repertoire directly from the VDJ matrix of the VDJ_GEX_matrix function
#' @param VDJ.matrix VDJ dataframe generated using the VDJ_GEX_matrix function. Plots will be made by sample and using the clonal frequencies in the clonotype_frequency column
#' @param counts.to.use How to count clonotypes and cells. If set to "10x" the function uses the clonotype_frequency column derived directly from the cellranger output to calculate expansion. If set to "VGM" the function will base its counts on the number of rows per clonotype in the VDJ.matrix. These two counts may diverge, if cells are filtered out due to overlapping barcodes or if a different clonotyping strategy was applied. Defaults to "VGM" 
#' @param label.size Size of text labels. All parameters below are purely for graphical purposes and optional. If necessary changes should be made in small (0.1) increments. I highly suggest to optimize these ONLY once a format for saving the plot is set. 
#' @param not.expanded.label.vjust Numeric. Regulates the vertical position of the label for non expanded cells 
#' @param not.expanded.label.hjust Numeric. Regulates the horizontal position of the label for non expanded cells
#' @param total.label.vjust Numeric. Regulates the vertical position of the center label
#' @param total.label.hjust Numeric. Regulates the horizontal position of the center label
#' @param expanded.colors Character vector. Colors to use for expanded clones. Should be more than 3 for better visibility. Defaults to a "darkorchid3"-based palette. 
#' @param non.expanded.color Character. Color to use for non expanded clones. Defaults to "black"
#' @param platypus.version Only v3 available. 
#' @return Returns a list of plots. One for each sample in the sample_id column
#' @export
#' @examples
#' \dontrun{
#' 
#' VDJ_clonal_donut(VDJ.matrix = VDJ.matrix.output[[1]])
#'
#'}
VDJ_clonal_donut <- function(VDJ.matrix,
                             counts.to.use,
                             label.size,
                             not.expanded.label.vjust,
                             not.expanded.label.hjust,
                             total.label.vjust,
                             total.label.hjust,
                             expanded.colors,
                             non.expanded.color,
                             platypus.version){

if(missing(counts.to.use)) counts.to.use <- "VGM"
if(missing(label.size)) label.size <- 5
if(missing(not.expanded.label.vjust)) not.expanded.label.vjust <- -0.2
if(missing(not.expanded.label.hjust)) not.expanded.label.hjust <- 1.4
if(missing(total.label.vjust)) total.label.vjust <- 3
if(missing(total.label.hjust)) total.label.hjust <- 0.5
if(missing(expanded.colors)) expanded.colors <- c("darkorchid4","darkorchid1","mediumorchid1","mediumpurple3")
if(missing(non.expanded.color)) non.expanded.color <- "black"

platypus.version = "v3"

VDJ.matrix <- subset(VDJ.matrix, clonotype_id_10x != "") #Filter possible cells with no clonotype. This can cause issues later

if(counts.to.use == "VGM"){

  print("Using counts of entries in the VDJ GEX matrix")
  
clonotypes <- VDJ.matrix %>% group_by(sample_id, clonotype_id_10x) %>% summarise(clonotype_frequency = n())
clonotypes$expanded <- F
clonotypes$expanded[clonotypes$clonotype_frequency > 1] <- T

} else if(counts.to.use == "10X" | counts.to.use == "10x"){
 
  print("Using counts provided by 10X in the clonotype_frequency column")
  
clonotypes <- VDJ.matrix %>% group_by(sample_id, clonotype_id_10x) %>% summarise(clonotype_frequency = as.numeric(clonotype_frequency[1]))
clonotypes$expanded <- F
clonotypes$expanded[clonotypes$clonotype_frequency > 1] <- T

}

plot.list <- list()
for(i in 1:length(unique(clonotypes$sample_id))){
  cur_c <- subset(clonotypes, sample_id == unique(clonotypes$sample_id)[i])
  
  print(head(cur_c))
  
  total_cells <- sum(cur_c$clonotype_frequency)
  expanded_cells <- sum(cur_c$clonotype_frequency[cur_c$expanded == T])
  nonexpanded_cells <- sum(cur_c$clonotype_frequency[cur_c$expanded == F])
  
  cur_c <- cur_c[order(cur_c$clonotype_frequency, decreasing = T),]
  cur_c$clonotype_id_10x[which(cur_c$expanded == F)] <- "1 cell"
  
  print("----------")
  print(paste0("Sample: ", unique(clonotypes$sample_id)[i]))
  
  print(paste0("Clones: Expanded: ", length(which(cur_c$expanded == T)), " / ", round((length(which(cur_c$expanded == T)) / nrow(cur_c)*100),2), "%; 1 cell ", length(which(cur_c$expanded == F)), " / ",round((length(which(cur_c$expanded == F)) / nrow(cur_c)*100),2), "%; total: ", nrow(cur_c)))
  
  print(paste0("Cells: Expanded: ", expanded_cells, " / ", round((expanded_cells / total_cells *100),2), "%; 1 cell ", nonexpanded_cells, " / ",round((nonexpanded_cells / total_cells *100),2), "%; Total: ", total_cells))
  
  
  cur_c$clonotype_id_10x_ord <- ordered(as.factor(cur_c$clonotype_id_10x), levels = unique(cur_c$clonotype_id_10x))
  
  if(nrow(cur_c) > 0){
  c_out <- c()
  for(j in 1:nrow(cur_c)){
    if(is.na(cur_c$expanded[j])){cur_c$expanded[j] <- FALSE}
    if(cur_c$expanded[j]==T){
      c_out <- c(c_out, rep(cur_c$clonotype_id_10x[j], cur_c$clonotype_frequency[j]))
    }
  }
  c_out <- c(c_out, rep("1 cell", length(which(cur_c$expanded == F))))
  c_out <- data.frame("clonotype_id_10x" = c_out, "sample_id" = unique(clonotypes$sample_id)[i])
  
  c_out$clonotype_id_10x <- ordered(as.factor(c_out$clonotype_id_10x), levels = levels(cur_c$clonotype_id_10x_ord))
  
  #generate palette
  pal_cl <- c(rep(expanded.colors,500)[1:(nrow(cur_c[which(cur_c$expanded == T),]))], non.expanded.color)
  
  plot_out <- ggplot(c_out, aes(x=sample_id, fill= clonotype_id_10x))+ geom_bar(width = 1, show.legend = F)+ scale_fill_manual(values = pal_cl)+ coord_polar("y", direction = -1)+ theme(panel.background = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), text = element_text(size=30), legend.key = element_rect(colour = "grey"),legend.key.height = NULL,legend.key.width = NULL, legend.position = "right",legend.direction = "vertical", plot.title = element_text(size = 20,hjust = 0.5)) + labs(x = "", y = "", title = paste0(unique(clonotypes$sample_id)[i])) + geom_text(x = 1, hjust = not.expanded.label.hjust, vjust = not.expanded.label.vjust, y = 1, label = paste0(length(which(cur_c$expanded == F))), color = "white", fontface = "bold", size = label.size) + geom_point(inherit.aes = F, aes(x = 0, y = 1), color = "white", size = 25) + geom_text(x = 1, y = 1, vjust = total.label.vjust, hjust = total.label.hjust, label = paste0(nrow(cur_c)," \n(", total_cells, ")"), color = "black", fontface = "bold", size = label.size) 
  plot.list[[i]] <- plot_out
  }
}
return(plot.list)
}

