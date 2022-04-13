#'Overlay clones on GEX projection
#'
#' @description Highlights the cells belonging to any number of top clonotypes or of specifically selected clonotypes from one or more samples or groups in a GEX dimensional reduction.
#' @param GEX A single seurat object from VDJ_GEX_matrix, which also includes VDJ information in the metadata (set integrate.VDJ.to.GEX to TRUE in the VDJ_GEX_matrix function) (VDJ_GEX_matrix.output[[2]]) ! Clone ids and frequencies are drawn from the columns "clonotype_id" and "clonotype_frequency"
#' @param reduction Character. Defaults to "umap". Name of the reduction to overlay clones on. Can be "pca", "umap", "tsne"
#' @param n.clones Integer. Defaults to 5. To PLOT TOP N CLONES. Number of Top clones to plot. If either by.sample or by.group is TRUE, n.clones clones from each sample or group will be overlayed
#' @param clones.to.plot Character. Alternative to n.clones. TO PLOT SPECIFIC CLONES. Must reference a column in the GEX@meta.data filled with TRUE and FALSE. Entries with TRUE label are plotted. Such a column may be generated using GEX@metadata$clones_to_plot_column <- GEX@metadata$Some_cell_identifier == "Interesting"
#' @param by.sample Boolean. Defaults to FALSE. Whether to overlay clones by sample. If set to TRUE this will generate a facet_wrap plot with as many facets as samples.
#' @param by.other.group Character string. Defaults to "none". Must be a valid column name of the metadata of the input seurat object. If so, this will generate a facet_wrap plot with as many facets unique entries in the specified column. This may be useful to plot cell type specific clones
#' @param ncol.facet Integer. Defaults to 2. Number of columns in the facet_wrap plot if by.sample or by.group is TRUE
#' @param pt.size Numeric. Defaults to 1. Size of points in DimPlot. Passed to Seurat::DimPlot
#' @param clone.colors Character vector. Defaults to rainbow(n.clones). Colors to use for individual clones. One can provide either a vector of length n.clones or a of length Nr. of samples/groups \* n.clones. In case that a vector of length n.clones is provided and by.group or by.sample is TRUE, colors are repeated for each sample/group
#' @param others.color Character. Color for cells that are not selected i.e. not part of the overlayed clonotypes. Defaults to "grey80". To hide the rest of the umap set to "white"
#' @param split.plot.and.legend Boolean. Defaults to FALSE. Whether to return the plot and the legend separately as a list. This can be useful if legends get large and distort the actual plots. The packages gridExtra and cowplot are required for this. If set to TRUE a list is returned where out[[1]] is the plot which can be printed just by executing out[[1]]; out[[2]] is the legend, which can be printed either using plot(out[[2]]) or grid.arrange(out[[2]])
#' @param platypus.version Character. At the moment this function runs only on the output of the VDJ_GEX_matrix function meaning that it is exclusively part of Platypus "v3". With further updates the functionality will be extended.
#' @return A ggplot object or a list of a ggplot and a gtable legend (if split.plot.and.legend \=\= TRUE). Theme, colors etc. may be changed directly by adding new elements to this output (e.g. out \+ theme_minimal())
#' @export
#' @examples
#'
#' #To return a single plot with top clones across samples
#' overlay_clones_plot <- VDJ_GEX_overlay_clones(
#' GEX = Platypus::small_vgm[[2]], reduction = "umap"
#' ,n.clones = 5, by.sample = FALSE
#' ,by.other.group = "none", pt.size = 1,split.plot.and.legend = FALSE)
#'
#' #To return a facet plot with top clones for each sample
#' overlay_clones_plot <- VDJ_GEX_overlay_clones(
#' GEX = Platypus::small_vgm[[2]], reduction = "umap"
#' ,n.clones = 5, by.sample = TRUE, by.other.group = "none"
#' ,pt.size = 1,ncol.facet = 2, split.plot.and.legend = FALSE)
#'
#' #To return a facet plot and the legend separately with top clones for each group
#' overlay_clones_plot <- VDJ_GEX_overlay_clones(
#' GEX = Platypus::small_vgm[[2]], reduction = "umap"
#' ,n.clones = 5, by.sample = TRUE, by.other.group = "group_id", pt.size = 1
#' ,ncol.facet = 2, split.plot.and.legend = TRUE)
#'
#' #To print both:
#' #overlay_clones_plot[[1]] #Plot
#' #gridExtra::grid.arrange(overlay_clones_plot[[2]]) #Legend
#' #To save, ggsave() is applicable to both
#'
#' #To return a single plot with selected clones
#' #add a clonotype_to_plot column
#' #GEX@meta.data$clonotype_to_plot <- GEX$VJ_vgene == "TRAV5-1"
#' #Column with TRUE for all clones with a particular V gene
#' #overlay_clones_plot <- VDJ_GEX_overlay_clones(GEX = GEX, reduction = "umap"
#' #, clones.to.plot = "clonotype_to_plot", by.sample = TRUE, by.other.group = "none"
#' #, split.plot.and.legend = FALSE, pt.size = 1.5)
#'

VDJ_GEX_overlay_clones <- function(GEX,
                                   reduction,
                                   n.clones,
                                   clones.to.plot,
                                   by.sample,
                                   by.other.group,
                                   ncol.facet,
                                   pt.size,
                                   clone.colors,
                                   others.color,
                                   split.plot.and.legend,
                                   platypus.version){
  sample_id <- NULL
  group_id <- NULL

  #VERSION is set for now:
  platypus.version <- "v3"

  #default is plotting top clones
  if(missing(clones.to.plot)) clones.to.plot <- "none"
  if(missing(n.clones)) n.clones <- 5
  if(n.clones > 10) n.clones <- 10
  if(missing(by.sample)) by.sample <- T
  if(missing(by.other.group)) by.other.group <- "none"
  if(by.sample == T & by.other.group != "none"){by.other.group <- "none"}
  if(missing(reduction)) reduction <- "umap"
  if(missing(ncol.facet)) ncol.facet <- 2
  if(missing(split.plot.and.legend)) split.plot.and.legend <- F
  if(missing(pt.size)) pt.size <- 1
  if(missing(others.color)) others.color <- "grey80"

  if(by.other.group %in% names(GEX@meta.data)){
    by.group <- T
    GEX@meta.data$group_id <- GEX@meta.data[,by.other.group]
  } else{by.group <- F}


  if(!clones.to.plot %in% names(GEX@meta.data) & clones.to.plot != "none"){
    stop("column specified as clones.to.plot was not found in GEX@metadata")
  }

  #make sure that the frequency column is numeric
  GEX@meta.data$clonotype_frequency <- as.numeric(as.character(GEX@meta.data$clonotype_frequency))

  #security check: if clonotype_id is not present in dataframe use clonotype_id_10x
  if(!"clonotype_id" %in% names(GEX@meta.data)){
    GEX@meta.data$clonotype_id <- GEX@meta.data$clonotype_id_10x
    message("Using clonotype_id_10x column")
  }

  #make column sample + clonotype id
  GEX@meta.data$s_cl <- paste0(GEX@meta.data$sample_id, GEX@meta.data$clonotype_id)

  #get subset unique that column + frequency
  cl_unique <- GEX@meta.data[which(duplicated(GEX@meta.data$s_cl) == F),c("s_cl", "sample_id","group_id", "clonotype_id", "clonotype_frequency")]
  if(clones.to.plot != "none"){ #if clones to plot is specified bind that column to cl_unique as well
    cl_unique <- cbind(cl_unique, GEX@meta.data[which(duplicated(GEX@meta.data$s_cl) == F),c(clones.to.plot)])
    names(cl_unique)[ncol(cl_unique)] <- "clonotype_to_plot"
  }
  #order by frequency
  cl_unique <- cl_unique[order(cl_unique$clonotype_frequency, decreasing = T),]

  if(by.sample == F & by.group == F){

    if(clones.to.plot == "none"){
      #get top n.clones
      if(n.clones > nrow(cl_unique)){n_clones_in <- nrow(cl_unique)} else {n_clones_in <- n.clones}
      cl_unique <- cl_unique[1:n_clones_in,]

    } else {
      #Select clones with TRUE label in the clone.to.plot column
      cl_unique <- cl_unique[cl_unique[,"clonotype_to_plot"],]
    }

    #remove any empty columns
    cl_unique <- cl_unique[is.na(cl_unique$clonotype_id) == F,]

    #open new column in original df
    GEX@meta.data$cl_to_plot <- "Not selected"

    #with loop add the clonal rank to the cells of the top n.clones
    #this array is used for ordering later
    track_cl_to_plot <- c()
    for(i in 1:nrow(cl_unique)){
      #paste together clonal rank and frequency
      GEX@meta.data$cl_to_plot[which(GEX@meta.data$s_cl == cl_unique$s_cl[i])] <- paste0(i, " / ",cl_unique$sample_id[i], " / ", cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i])
      track_cl_to_plot <- append(track_cl_to_plot, paste0(i, " / ",cl_unique$sample_id[i], " / ",cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i]))
    }
    #order for plotting
    GEX@meta.data$cl_to_plot <- ordered(as.factor(GEX@meta.data$cl_to_plot), levels = c("Not selected", track_cl_to_plot))

    #Get colors right. Output is a vector of colors shorter by 1 from what is needed. That last color is grey80 and specified in the DimPlot for cells that are not highlighted
    if(missing(clone.colors)) clone.colors <- grDevices::rainbow(length(unique(GEX@meta.data$cl_to_plot))-1)
    if(length(clone.colors) != length(unique(GEX@meta.data$cl_to_plot))-1){stop(paste0("Nr of supplied colors ", length(clone.colors)  ," does not match number of clones to plot ", length(unique(GEX@meta.data$cl_to_plot))-1))}

    #Prep for highlighting cells function

    SeuratObject::Idents(object = GEX) <- "cl_to_plot"
    to_highlight_list <- list()
    for(i in 1:length(track_cl_to_plot)){
      to_highlight_list[[i]] <- SeuratObject::WhichCells(GEX, ident =  track_cl_to_plot[i])
    }
    to_highlight_list <- rev(to_highlight_list) #reversing it to keep plotting order as orders of selected clones

    out.plot <- Seurat::DimPlot(GEX,reduction = reduction, cells.highlight = to_highlight_list, cols.highlight = clone.colors , cols = others.color, shuffle = F, pt.size = pt.size) + ggplot2::labs(col = "Rank / Sample id / Clonotype / Frequency", title = "" ) + ggplot2::scale_color_manual(values = c(others.color, clone.colors), labels = c("Not selected", track_cl_to_plot))

    #dimplot with cols + grey30 to color in the non selected clones
    #out.plot <- Seurat::DimPlot(GEX,reduction = reduction, group.by = "cl_to_plot", cols = c(clone.colors, others.color), shuffle = F, order = c("Not selected","selected"), pt.size = pt.size) + ggplot2::labs(col = "Rank / Sample id / Clonotype / Frequency", title = "")
    #check whether plot and legend should be split
    if(split.plot.and.legend == F){
      return(out.plot)
    } else {
      out.leg <- cowplot::get_legend(out.plot)
      out.plot <- out.plot + ggplot2::theme(legend.position = "none")
      return(list("plot" = out.plot, "legend" = out.leg))
    }
  }

  if(by.sample == T){
    #same routine but now faceted by sample

    if(clones.to.plot == "none"){
      #get top n.clones per sample
      cl_unique_all <- list()
      for(i in 1:length(unique(cl_unique$sample_id))){
        cl_unique_cur <- subset(cl_unique, sample_id == unique(cl_unique$sample_id)[i])

        if(n.clones > nrow(cl_unique_cur)){n_clones_in <- nrow(cl_unique_cur)} else {n_clones_in <- n.clones}
        cl_unique_all[[i]] <- cl_unique_cur[1:n_clones_in,]
      }
      cl_unique <- dplyr::bind_rows(cl_unique_all)

    } else {
      #Select clones with TRUE label
      cl_unique <- cl_unique[cl_unique[,"clonotype_to_plot"],]
    }

    #remove any empty columns
    cl_unique <- cl_unique[is.na(cl_unique$clonotype_id) == F,]

    #open new column in original df
    GEX@meta.data$cl_to_plot <- "Not selected"

    #with loop add the clonal rank to the cells of the top n.clones
    #this array is used for ordering later
    track_cl_to_plot <- c()
    #for ranks
    r <- 1
    for(i in 1:nrow(cl_unique)){
      #paste together clonal rank and frequency
      GEX@meta.data$cl_to_plot[which(GEX@meta.data$s_cl == cl_unique$s_cl[i])] <- paste0(r, " / ",cl_unique$sample_id[i], " / ", cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i])
      track_cl_to_plot <- append(track_cl_to_plot, paste0(r, " / ", cl_unique$sample_id[i], " / ", cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i]))

      if(r == n.clones){r <- 1} else {r <- r + 1}
    }
    #order for plotting
    GEX@meta.data$cl_to_plot <- ordered(as.factor(GEX@meta.data$cl_to_plot), levels = c(track_cl_to_plot, "Not selected"))

    #Here color selection is a little messy. Essentially if there are no colors provided than colors are generated by rainbow and replicated for each sample
    #If colors are present we check whether it is the correct lenght or whether we can replicated the colors by sample to make it the correct length.
    if(clones.to.plot == "none"){
      if(missing(clone.colors)) clone.colors <- rep(grDevices::rainbow(n.clones), length(unique(GEX@meta.data$sample_id)))
      if(length(clone.colors) == n.clones) clone.colors <- rep(clone.colors, length(unique(GEX@meta.data$sample_id)))
      if(length(clone.colors) != n.clones * length(unique(GEX@meta.data$sample_id))){stop("Nr of supplied colors does not match number of clones to plot. Provide either a vector of lenght n.clones or a of length Nr of samples/groups * n.clones")}
    } else{
      #because it is hard to know how many selected clones are present in each sample, if no colors are supplied we define colors individually for each clone. No repping per sample here
      if(missing(clone.colors)) clone.colors <- grDevices::rainbow(length(unique(GEX@meta.data$cl_to_plot))-1)
      if(length(clone.colors) != length(unique(GEX@meta.data$cl_to_plot))-1){stop(paste0("Nr of supplied colors ", length(clone.colors)  ," does not match number of clones to plot ", length(unique(GEX@meta.data$cl_to_plot))-1))}
    }

    #if too many colors were supplied we shorten the vector.
    if(length(unique(GEX$cl_to_plot)) < length(clone.colors)){
      clone.colors <- clone.colors[1:(length(unique(GEX$cl_to_plot))-1)]
    }

    #dimplot with cols + grey30 to color in the non selected clones
    out.plot <- Seurat::DimPlot(GEX,reduction = reduction, group.by = "cl_to_plot", cols = c(clone.colors,others.color), shuffle = F, pt.size = pt.size) + ggplot2::labs(col = "Rank / Sample id / Clonotype / Frequency", title = "") + ggplot2::facet_wrap(~GEX@meta.data$sample_id, ncol = ncol.facet)
    #! Bug in the current version of seurat: shuffle = T does not work with split_by or + facet_wrap(). The bug is reported and should be fixed soon (13.4.21) https://github.com/satijalab/seurat/issues/4300

    if(split.plot.and.legend == F){
      return(out.plot)
    } else {
      out.leg <- cowplot::get_legend(out.plot)
      out.plot <- out.plot + ggplot2::theme(legend.position = "none")
      return(list("plot" = out.plot, "legend" = out.leg))
    }

  }

  if(by.group == T){
    #same routine but now faceted by group

    if(clones.to.plot == "none"){
      #get top n.clones per group
      cl_unique_all <- list()
      for(i in 1:length(unique(cl_unique$group_id))){
        cl_unique_cur <- subset(cl_unique, group_id == unique(cl_unique$group_id)[i])

        if(n.clones > nrow(cl_unique_cur)){n_clones_in <- nrow(cl_unique_cur)} else {n_clones_in <- n.clones}
        cl_unique_all[[i]] <- cl_unique_cur[1:n_clones_in,]
      }
      cl_unique <- dplyr::bind_rows(cl_unique_all)

    } else {
      #Select clones with TRUE label
      cl_unique <- cl_unique[cl_unique[,"clonotype_to_plot"],]
    }

    #remove any empty columns
    cl_unique <- cl_unique[is.na(cl_unique$clonotype_id) == F,]

    #open new column in original df
    GEX@meta.data$cl_to_plot <- "Not selected"

    #with loop add the clonal rank to the cells of the top n.clones
    #this array is used for ordering later
    track_cl_to_plot <- c()
    #for ranks
    r <- 1
    for(i in 1:nrow(cl_unique)){
      #paste together clonal rank and frequency
      GEX@meta.data$cl_to_plot[which(GEX@meta.data$s_cl == cl_unique$s_cl[i])] <- paste0(r, " / ",cl_unique$group_id[i], " / ",cl_unique$sample_id[i], " / ", cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i])
      track_cl_to_plot <- append(track_cl_to_plot, paste0(r, " / ",cl_unique$group_id[i], " / ",cl_unique$sample_id[i], " / ", cl_unique$clonotype_id[i], " / ", cl_unique$clonotype_frequency[i]))

      if(r == n.clones){r <- 1} else {r <- r + 1}
    }
    #order for plotting
    GEX@meta.data$cl_to_plot <- ordered(as.factor(GEX@meta.data$cl_to_plot), levels = c(track_cl_to_plot, "Not selected"))

    if(clones.to.plot == "none"){
      if(missing(clone.colors)) clone.colors <- rep(grDevices::rainbow(n.clones), length(unique(GEX@meta.data$group_id)))
      if(length(clone.colors) == n.clones) clone.colors <- rep(clone.colors, length(unique(GEX@meta.data$group_id)))
      if(length(clone.colors) != n.clones * length(unique(GEX@meta.data$group_id))){stop("Nr of supplied colors does not match number of clones to plot. Provide either a vector of lenght n.clones or a of length Nr of samples/groups * n.clones")}
    } else{
      if(missing(clone.colors)) clone.colors <- grDevices::rainbow(length(unique(GEX@meta.data$cl_to_plot))-1)
      if(length(clone.colors) != length(unique(GEX@meta.data$cl_to_plot))-1){stop(paste0("Nr of supplied colors ", length(clone.colors)  ," does not match number of clones to plot ", length(unique(GEX@meta.data$cl_to_plot))-1))}
    }
    if(length(unique(GEX$cl_to_plot)) < length(clone.colors)){
      clone.colors <- clone.colors[1:(length(unique(GEX$cl_to_plot))-1)]
    }

    #dimplot with cols + grey30 to color in the non selected clones
    out.plot <- Seurat::DimPlot(GEX,reduction = reduction, group.by = "cl_to_plot", cols = c(clone.colors,others.color), shuffle = F, pt.size = pt.size) + ggplot2::labs(col = "Rank / Group id / Sample id / Clonotype / Frequency", title = "") + ggplot2::facet_wrap(~GEX@meta.data$group_id, ncol = ncol.facet)
    #! Bug in the current version of seurat: shuffle = T does not work with split_by or + facet_wrap(). The bug is reported and should be fixed soon (13.4.21) https://github.com/satijalab/seurat/issues/4300

    if(split.plot.and.legend == F){
      return(out.plot)
    } else {
      out.leg <- cowplot::get_legend(out.plot)
      out.plot <- out.plot + ggplot2::theme(legend.position = "none")
      return(list("plot" = out.plot, "legend" = out.leg))
    }
  }
}
