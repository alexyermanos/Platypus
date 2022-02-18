#' Generate circular plots of clonal expansion per repertoire directly from the VDJ matrix of the VDJ_GEX_matrix function
#' @param VDJ VDJ dataframe generated using the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]]). Plots will be made by sample and using the clonal frequencies specified by counts.to.use
#' @param counts.to.use How to count clonotypes and cells. A column name of the VDJ matrix containing clonotype IDs. This defaults to "clonotype_id_10x", which reflects clonotypes by Cellranger in an unaltered VGM. To use counts from the VDJ_clonotype_v3 function set this parameter to the relevant column e.g. "clonotype_id_cdr.aa" or   "global_clonotype_id_cdr.aa" are two examples.
#' @param label.size Size of text labels. All parameters below are purely for graphical purposes and optional. If necessary changes should be made in small (0.1) increments. ! It is recommended to optimize these ONLY once a format for saving the plot is set.
#' @param not.expanded.label.vjust Numeric. Regulates the vertical position of the label for non expanded cells
#' @param not.expanded.label.hjust Numeric. Regulates the horizontal position of the label for non expanded cells
#' @param total.label.vjust Numeric. Regulates the vertical position of the center label
#' @param total.label.hjust Numeric. Regulates the horizontal position of the center label
#' @param expanded.colors Character vector. Colors to use for expanded clones. Should be more than 3 for better visibility. Defaults to a "darkorchid3"-based palette.
#' @param non.expanded.color Character. Color to use for non expanded clones. Defaults to "black"
#' @return Returns a list of circular plots showing proportions of expanded clones and non-expanded clones. One plot is generated for each sample in the sample_id column
#' @export
#' @examples
#' VDJ_clonal_donut(VDJ = Platypus::small_vgm[[1]])
#'
VDJ_clonal_donut <- function(VDJ,
                             counts.to.use,
                             label.size,
                             not.expanded.label.vjust,
                             not.expanded.label.hjust,
                             total.label.vjust,
                             total.label.hjust,
                             expanded.colors,
                             non.expanded.color){

  clonotype_id <- NULL
  sample_id <- NULL
  clonotype_frequency <- NULL

if(missing(counts.to.use)) counts.to.use <- "clonotype_id_10x"
if(missing(label.size)) label.size <- 5
if(missing(not.expanded.label.vjust)) not.expanded.label.vjust <- -0.2
if(missing(not.expanded.label.hjust)) not.expanded.label.hjust <- 1.4
if(missing(total.label.vjust)) total.label.vjust <- 3
if(missing(total.label.hjust)) total.label.hjust <- 0.5
if(missing(expanded.colors)) expanded.colors <- c("darkorchid4","darkorchid1","mediumorchid1","mediumpurple3")
if(missing(non.expanded.color)) non.expanded.color <- "black"

platypus.version = "v3"

VDJ <- subset(VDJ, clonotype_id != "") #Filter possible cells with no clonotype. This can cause issues later

if(!counts.to.use %in% names(VDJ)){
  warning("Column name for counts.to.use was not found in VDJ. Defaulting to 'clonotype_id_10x'")
  counts.to.use <- "clonotype_id"
}

if(counts.to.use %in% names(VDJ)){
  message(paste0("Using column ", counts.to.use, " for counting clones"))
  counts.to.use <- "clonotype_id"

  VDJ$for_clonal_donut <- VDJ[,counts.to.use]

  clonotypes <- VDJ %>% dplyr::group_by(sample_id, for_clonal_donut) %>%  dplyr::summarise(clonotype_frequency = dplyr::n())
  clonotypes$expanded <- F
  clonotypes$expanded[clonotypes$clonotype_frequency > 1] <- T
}

plot.list <- list()
for(i in 1:length(unique(clonotypes$sample_id))){
  cur_c <- subset(clonotypes, sample_id == unique(clonotypes$sample_id)[i])

  total_cells <- sum(cur_c$clonotype_frequency)
  expanded_cells <- sum(cur_c$clonotype_frequency[cur_c$expanded == T])
  nonexpanded_cells <- sum(cur_c$clonotype_frequency[cur_c$expanded == F])

  cur_c <- cur_c[order(cur_c$clonotype_frequency, decreasing = T),]
  cur_c$for_clonal_donut[which(cur_c$expanded == F)] <- "1 cell"

  message(paste0("Clones: Expanded: ", length(which(cur_c$expanded == T)), " / ", round((length(which(cur_c$expanded == T)) / nrow(cur_c)*100),2), "%; 1 cell ", length(which(cur_c$expanded == F)), " / ",round((length(which(cur_c$expanded == F)) / nrow(cur_c)*100),2), "%; total: ", nrow(cur_c)))

  message(paste0("Cells: Expanded: ", expanded_cells, " / ", round((expanded_cells / total_cells *100),2), "%; 1 cell ", nonexpanded_cells, " / ",round((nonexpanded_cells / total_cells *100),2), "%; Total: ", total_cells))


  cur_c$clonotype_id_ord <- ordered(as.factor(cur_c$for_clonal_donut), levels = unique(cur_c$for_clonal_donut))

  if(nrow(cur_c) > 0){
  c_out <- c()
  for(j in 1:nrow(cur_c)){
    if(is.na(cur_c$expanded[j])){cur_c$expanded[j] <- FALSE}
    if(cur_c$expanded[j]==T){
      c_out <- c(c_out, rep(cur_c$for_clonal_donut[j], cur_c$clonotype_frequency[j]))
    }
  }
  c_out <- c(c_out, rep("1 cell", length(which(cur_c$expanded == F))))
  c_out <- data.frame("clonotype_id" = c_out, "sample_id" = unique(clonotypes$sample_id)[i])

  c_out$clonotype_id <- ordered(as.factor(c_out$clonotype_id), levels = levels(cur_c$clonotype_id_ord))

  #generate palette
  pal_cl <- c(rep(expanded.colors,500)[1:(nrow(cur_c[which(cur_c$expanded == T),]))], non.expanded.color)

  plot_out <- ggplot2::ggplot(c_out, ggplot2::aes(x=sample_id, fill= clonotype_id))+ ggplot2::geom_bar(width = 1, show.legend = F)+ ggplot2::scale_fill_manual(values = pal_cl)+ ggplot2::coord_polar("y", direction = -1)+ ggplot2::theme(panel.background = ggplot2::element_blank(),axis.text.x = ggplot2::element_blank(),axis.text.y = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), text = ggplot2::element_text(size=30), legend.key = ggplot2::element_rect(colour = "grey"),legend.key.height = NULL,legend.key.width = NULL, legend.position = "right",legend.direction = "vertical", plot.title = ggplot2::element_text(size = 20,hjust = 0.5)) + ggplot2::labs(x = "", y = "", title = paste0(unique(clonotypes$sample_id)[i])) + ggplot2::geom_text(x = 1, hjust = not.expanded.label.hjust, vjust = not.expanded.label.vjust, y = 1, label = paste0(length(which(cur_c$expanded == F))), color = "white", fontface = "bold", size = label.size) + ggplot2::geom_point(inherit.aes = F, ggplot2::aes(x = 0, y = 1), color = "white", size = 25) + ggplot2::geom_text(x = 1, y = 1, vjust = total.label.vjust, hjust = total.label.hjust, label = paste0(nrow(cur_c)," \n(", total_cells, ")"), color = "black", fontface = "bold", size = label.size)
  plot.list[[i]] <- plot_out
  }
}
return(plot.list)
}
