#' Generate circular plots of clonal expansion per repertoire directly from the VDJ matrix of the VDJ_GEX_matrix function
#' @param VDJ.matrix VDJ dataframe generated using the VDJ_GEX_matrix function. Plots will be made by sample and using the clonal frequencies in the clonotype_frequency column
#' @param label.size Size of text labels. All parameters below are purely graphic and optional. If necessary changes should be made in small (0.1) increments. I highly suggest to optimize these only once a format for saving the plot set. 
#' @param not.expanded.label.vjust Numeric. Regulates the vertical position of the label for non expanded cells 
#' @param not.expanded.label.hjust Numeric. Regulates the horizontal position of the label for non expanded cells
#' @param total.label.vjust Numeric. Regulates the vertical position of the center label
#' @param total.label.hjust Numeric. Regulates the horizontal position of the center label
#' @param expanded.colors Character vector. Colors to use for expanded clones. Should be more than 3 for better visibility. Defaults to a "darkorchid3"-based palette. 
#' @param non.expanded.color Character. Color to use for non expanded clones.Defaults to "black"
#' @param platypus.version Only v3 available. 
#' @return Returns a list of plots. One for each repertoire
#' @export
#' @examples
#' \dontrun{
#' 
#' VDJ_expansion_pie(VDJ.matrix = VDJ.matrix.output[[1]])
#'
#'}
VDJ_expansion_pie <- function(VDJ.matrix,
                              label.size,
                              not.expanded.label.vjust,
                              not.expanded.label.hjust,
                              total.label.vjust,
                              total.label.hjust,
                              expanded.colors,
                              non.expanded.color,
                              platypus.version){


if(missing(label.size)) label.size <- 5
if(missing(not.expanded.label.vjust)) not.expanded.label.vjust <- -0.2
if(missing(not.expanded.label.hjust)) not.expanded.label.hjust <- 1.4
if(missing(total.label.vjust)) total.label.vjust <- 3
if(missing(total.label.hjust)) total.label.hjust <- 0.7
if(missing(expanded.colors)) expanded.colors <- c("darkorchid4","darkorchid1","mediumorchid1","mediumpurple3")
if(missing(non.expanded.color)) non.expanded.color <- "black"

platypus.version = "v3"

clonotypes <- VDJ.matrix %>% group_by(sample_id, clonotype_id_10x) %>% summarise(freq = as.numeric(clonotype_frequency[1]), expanded = (clonotype_frequency[1] > 1))

cells <- VDJ.matrix %>% group_by(sample_id, barcode) %>% summarise(freq = as.numeric(clonotype_frequency[1]), expanded = (clonotype_frequency[1] > 1))

plot.list <- list()
for(i in 1:length(unique(clonotypes$sample_id))){
  cur_c <- subset(clonotypes, sample_id == unique(clonotypes$sample_id)[i])
  
  cur_c <- cur_c[order(cur_c$freq, decreasing = T),]
  cur_c$clonotype_id_10x[which(cur_c$expanded == F)] <- "1 cell"
  
  print("----------")
  print(paste0("Sample: ", unique(clonotypes$sample_id)[i]))
  
  print(paste0("Clones: Expanded: ", length(which(cur_c$expanded == T)), " / ", (length(which(cur_c$expanded == T)) / nrow(cur_c)*100), "%; 1 cell ", length(which(cur_c$expanded == F)), " / ",(length(which(cur_c$expanded == F)) / nrow(cur_c)*100)))
  
  cur_cells <- subset(cells, sample_id == unique(clonotypes$sample_id)[i])
  
  print(paste0("Cells: Expanded: ", length(which(cur_cells$expanded == T)), " / ", (length(which(cur_cells$expanded == T)) / nrow(cur_cells)*100), "%; 1 cell ", length(which(cur_cells$expanded == F)), " / ",(length(which(cur_cells$expanded == F)) / nrow(cur_cells)*100)))
  
  print("----------")
  
  cur_c$clonotype_id_10x_ord <- ordered(as.factor(cur_c$clonotype_id_10x), levels = unique(cur_c$clonotype_id_10x))
  levels(cur_c$clonotype_id_10x)
  c_out <- c()
  for(j in 1:nrow(cur_c)){
    if(cur_c$expanded[j]==T){
      c_out <- c(c_out, rep(cur_c$clonotype_id_10x[j], cur_c$freq[j]))
    }
  }
  c_out <- c(c_out, rep("1 cell", length(which(cur_c$expanded == F))))
  c_out <- data.frame("clonotype_id_10x" = c_out, "sample_id" = unique(clonotypes$sample_id)[i])
  
  c_out$clonotype_id_10x <- ordered(as.factor(c_out$clonotype_id_10x), levels = levels(cur_c$clonotype_id_10x_ord))
  
  #generate palette
  pal_cl <- c(rep(expanded.colors,500)[1:(nrow(cur_c[which(cur_c$expanded == T),]))], non.expanded.color)
  
  plot_out <- ggplot(c_out, aes(x=sample_id, fill= clonotype_id_10x))+ geom_bar(width = 1, show.legend = F)+ scale_fill_manual(values = pal_cl)+ coord_polar("y", direction = -1)+ theme(panel.background = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), text = element_text(size=30), legend.key = element_rect(colour = "grey"),legend.key.height = NULL,legend.key.width = NULL, legend.position = "right",legend.direction = "vertical", plot.title = element_text(size = 20,hjust = 0.5)) + labs(x = "", y = "", title = paste0(unique(clonotypes$sample_id)[i])) + geom_text(x = 1, hjust = not.expanded.label.hjust, vjust = not.expanded.label.vjust, y = 1, label = paste0(length(which(cur_c$expanded == F))), color = "white", fontface = "bold", size = 5) + geom_point(inherit.aes = F, aes(x = 0, y = 1), color = "white", size = 25) + geom_text(x = 1, y = 1, vjust = total.label.vjust, hjust = total.label.hjust, label = paste0(nrow(cur_c)," \n(", nrow(cur_cells), ")"), color = "black", fontface = "bold", size = 5) 
  plot.list[[i]] <- plot_out
  }
return(plot.list)
}


