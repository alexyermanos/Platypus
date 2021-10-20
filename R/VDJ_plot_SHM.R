#' Plots for SHM based on MIXCR output generated using the VDJ_call_MIXCR function and appended to the VDJ.GEX.matrix.output
#' @param VDJ.mixcr.matrix Output dataframe from the VDJ_call_MIXCR function or a dataframe generated using the VDJ_GEX_matrix function and supplemented with MIXCR information
#' @param group.by Character. Defaults to "sample_id". Column name of VDJ.matrix to split VDJ.matrix by. For each unique entry in that column a set of plots will be generated. This can be useful to plot SHM by expansion or by transcriptomics-derived clusters
#' @param quantile.label Numeric. Defaults to 0.9. Which points to label in the SHM scatterplot. If set to 0.9, the top 10\% of cells by SHM number will be labelled. If ggrepel throws a warning, concerning overlap it is recommended to attempt to lable less points to avoid cluttering
#' @param point.size Size of points in plots. Passed to geom_jitter()
#' @param mean.line.color Color of mean bar in dotplots. Passed to geom_errorbar()
#' @param stats.to.console Boolean. Defaults to FALSE. Prints basic statistics (AOV \+ post hoc test) to console
#' @param platypus.version Character. Only "v3" available.
#' @return Returns a list of ggplot objects. out\[\[1\]\] is a boxplot comparing SHM by group.by. out\[\[2\]\] to out\[\[n\]\] are plots for each group that visualize VDJ and VJ SHM distribution for each group. Data for any plot can be accessed via out \[\[any\]\]$data
#' @export
#' @examples
#'#Simulating SHM data
#'small_vgm <- Platypus::small_vgm
#'small_vgm[[1]]$VDJ_SHM <- as.integer(rnorm(nrow(small_vgm[[1]]), mean = 5, sd = 3))
#'small_vgm[[1]]$VJ_SHM <- as.integer(rnorm(nrow(small_vgm[[1]]), mean = 5, sd = 3))
#'
#' #Standard plots
#' SHM_plots <- VDJ_plot_SHM(VDJ = small_vgm[[1]]
#' , group.by = "sample_id", quantile.label = 0.9)
#'
#' #Group by transcriptional cluster and label only top 1\%
#' SHM_plots <- VDJ_plot_SHM(VDJ = small_vgm[[1]]
#' , group.by = "seurat_clusters", quantile.label = 0.99)
#'

VDJ_plot_SHM <- function(VDJ.mixcr.matrix,
                     group.by,
                     quantile.label,
                     point.size,
                     mean.line.color,
                     stats.to.console,
                     platypus.version){
  name <- NULL
  group <- NULL
  value <- NULL
  m <- NULL
  VDJ_SHM <- NULL
  VJ_SHM <- NULL
  barcode <- NULL


  platypus.version <- "v3"
  if(missing(quantile.label)) quantile.label <- 0.9
  if(missing(group.by)) group.by <- "sample_id"
  if(missing(point.size)) point.size <- 2
  if(missing(mean.line.color)) mean.line.color <- "black"
  if(missing(stats.to.console)) stats.to.console <- F

  VDJ.matrix <- VDJ.mixcr.matrix

  #get data
  if(!"VDJ_SHM" %in% names(VDJ.matrix) | !"VJ_SHM" %in% names(VDJ.matrix)){
    stop("VDJ_SHM or VJ_SHM column not found in the input dataframe. Please provide a output dataframe of the VDJ_call_MIXCR function")
  }

  if((group.by %in% names(VDJ.matrix)) == F){
    stop(paste0("Specified group.by column " , group.by, " was not found in VDJ.matrix input dataframe. Please provide a valid column name"))
  }

  to_plot <- VDJ.matrix[,c(group.by, "barcode", "VDJ_SHM", "VJ_SHM")]
  names(to_plot)[1] <- "group"

  to_plot$VDJ_SHM[is.na(to_plot$VDJ_SHM)] <- 0
  to_plot$VJ_SHM[is.na(to_plot$VJ_SHM)] <- 0

  to_plot_long <- tidyr::pivot_longer(to_plot, cols = c(3:4))

  #BOXPLOT COMPARING SHM PER GROUP

  means <- to_plot_long %>% dplyr::group_by(name, group) %>% dplyr::summarise(m= mean(value))

  box_plot <- ggplot2::ggplot(to_plot_long, ggplot2::aes(color = group, y= value, x= group)) + ggplot2::geom_jitter(alpha = 0.4, width = 0.35, size = point.size) + ggplot2::geom_errorbar(inherit.aes = F, data = means, ggplot2::aes(ymax = m, ymin = m, x = group), color = mean.line.color, size = 1.6, width = 0.85) + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::ylab("SHM") + ggplot2::xlab("") + ggplot2::ggtitle(label = paste0("SHM per ", group.by))+ ggplot2::theme(strip.background = ggplot2::element_rect(color = "white", fill = "white")) + ggplot2::facet_wrap(~name)

  #SIGNIFICANCE TESTING

  SHM_iso_VDJ <- subset(box_plot$data, name == "VDJ_SHM")
  SHM_iso_VJ <- subset(box_plot$data, name == "VJ_SHM")

  vdjmean <- SHM_iso_VDJ %>% dplyr::group_by(group) %>% dplyr::summarize(m = mean(value))
  vjmean <- SHM_iso_VJ %>% dplyr::group_by(group) %>% dplyr::summarize(m = mean(value))

  #VDJ
  if(stats.to.console) message("VDJ chain SHM data summary")
  if(stats.to.console) message(as.data.frame(dplyr::group_by(SHM_iso_VDJ, group) %>%
  dplyr::summarise(
      count = dplyr::n(),
      mean = mean(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      IQR = stats::IQR(value, na.rm = TRUE)
    )))
  if(stats.to.console) message("\n kruskal.test()")
  if(stats.to.console) message(stats::kruskal.test(value ~ group, data = SHM_iso_VDJ))
  if(stats.to.console) message("\n pairwise.wilcox.test(p.adjust.method = 'BH')")
  if(stats.to.console) suppressWarnings(message(stats::pairwise.wilcox.test(SHM_iso_VDJ$value, SHM_iso_VDJ$group,p.adjust.method = "BH")))

  #VJ
  if(stats.to.console) message("------------")
  if(stats.to.console) message("VJ chain SHM data summary")
  if(stats.to.console) message(as.data.frame(dplyr::group_by(SHM_iso_VJ, group) %>%
          dplyr::summarise(
            count = dplyr::n(),
            mean = mean(value, na.rm = TRUE),
            sd = stats::sd(value, na.rm = TRUE),
            median = stats::median(value, na.rm = TRUE),
            IQR = stats::IQR(value, na.rm = TRUE)
          )))
  if(stats.to.console) message("\n kruskal.test()")
  if(stats.to.console) message(stats::kruskal.test(value ~ group, data = SHM_iso_VJ))
  if(stats.to.console) message("\n pairwise.wilcox.test(p.adjust.method = 'BH')")
  if(stats.to.console) suppressWarnings(message(stats::pairwise.wilcox.test(SHM_iso_VJ$value, SHM_iso_VJ$group,p.adjust.method = "BH")))

  if(stats.to.console) message("\n Please refer to output[[1]] to view boxplot showing group comparions tested here")

  #SCATTERPLOT FOR EVERY GROUP

  out.list <- list()
  out.list[[1]] <- box_plot

  for(j in 1:length(unique(to_plot$group))){

    curr_to_plot <- subset(to_plot, group == unique(to_plot$group)[j])

    qx_HC <- stats::quantile(curr_to_plot$VDJ_SHM, probs = quantile.label)
    qx_LC <- stats::quantile(curr_to_plot$VJ_SHM, probs = quantile.label)

    pos <- ggplot2::position_jitter(width = 0.3, seed = 2)

    out.list[[j+1]] <- ggplot2::ggplot(curr_to_plot, ggplot2::aes(x = VDJ_SHM, y = VJ_SHM, col = VDJ_SHM + VJ_SHM)) + ggplot2::geom_jitter(show.legend = T, size = 3, alpha = 0.8, position = pos) + ggplot2::theme_bw() + ggplot2::theme_classic() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + ggplot2::ylab("VJ SHM") + ggplot2::xlab("VDJ SHM") + ggplot2::ggtitle(label = paste0("SHM in ", unique(to_plot$group)[j]))+ ggrepel::geom_text_repel(inherit.aes = F, data = subset(curr_to_plot, VDJ_SHM > qx_HC | VJ_SHM > qx_LC), ggplot2::aes(x = VDJ_SHM, y = VJ_SHM, label = barcode), color = "black", position = pos) + ggplot2::scale_color_viridis_c(option = "B", end = 0.9)

  }
  return(out.list)
}
