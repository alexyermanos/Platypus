#' Selecion of cells by area
#' @description Allows to select an area on the spatial image and to isolate the cells expressed on this part and repeat this process several times.
#' @param vgm_VDJ Data frame containing all the data on the cell. It must contain the column clonotype_id which describes the number of the clonotype to which the cell belongs. This data frame can be obtained by the assignment functions (VDJ_assignment_random_based, VDJ_assignment_density_based and VDJ_assignment_germline_based).
#' @param alpha Number that give the transparency coefficient (value between 0 and 1). If it is not given it will automatically be 0.
#' @param bcs_merge Data frame containing imagerow, imagecol and barcode of the cells belonging to the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 10.
#' @param images_tibble Tbl-df containing the sample name, grob, height and width of the spatial image. It can also be created by the function scaling_spatial_image_parameter by selecting the output parameter 5.
#' @param sample_names Character vector containing the name of the sample.
#' @param nbpoints Numerical value that limite the maximum number of mouse click for the selection, default = 100.
#' @param title Character vector to name the plot.
#' @param size Number, to define the size of the text, default = 15.
#' @param plotting Character vector to return (TRUE) or not (FALSE) the plot of the selection
#' @return If plotting = TRUE, returns a list containing [[1]] the plot of the selected cells according to their group, [[2]] a data frame that contains all the cells but the selected cells are distinguished. If plotting = FALSE it juste returns the dataframe.
#' @export
#' @examples
#' \dontrun{
#' test<-Spatial_selection_of_cells_on_image(vgm_VDJ = vgm_spatial_simulated$VDJ$B_cells$random_BCR_assignment,
#' images_tibble = scaling_parameters[[5]],bcs_merge = scaling_parameters[[10]],sample_names = sample_names,
#' plotting = TRUE)
#'}

Spatial_selection_of_cells_on_image<-function(vgm_VDJ,alpha,bcs_merge,images_tibble,sample_names,nbpoints,title,size,plotting){
  if(missing(nbpoints)) nbpoints <- 100
  if(missing(alpha)) alpha <- 0
  if(missing(images_tibble)) images_tibble <- scaling_parameters[[5]]
  if(missing(bcs_merge)) bcs_merge <- scaling_parameters[[10]]
  if(missing(title)){
    title=""
  }
  if(missing(size)){
    size = 15
  }
  if(missing(plotting)){
    plotting="FALSE"
  }
  if(missing(sample_names)) stop("Please provide sample_names input for this function")
  if(missing(vgm_VDJ)) stop("Please provide vgm_VDJ input for this function")
  if (!require("dplyr", character.only = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }
  #ggplot-----------------------------------------------------------------
  spatial_plot<-ggplot2::ggplot(vgm_VDJ, ggplot2::aes(x, y)) +
    geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+
    ggplot2::coord_cartesian(expand=FALSE)+
    ggplot2::xlim(0,max(bcs_merge %>%
                 filter(sample ==sample_names[1]) %>%
                 select(width)))+
    ggplot2::ylim(max(bcs_merge %>%
               filter(sample ==sample_names[1]) %>%
               select(height)),0)+
    ggplot2::xlab("") +
    ggplot2::ylab("")+
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=3)))+
    ggplot2::theme_set(ggplot2::theme_bw(base_size = 10))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.ticks = ggplot2::element_blank())
  #common ggplot structure
  ggobj<-ggplot2::ggplot_build(spatial_plot)

  # Extract coordinates --------------------------------------------------
  xr <- ggobj$layout$panel_params[[1]]$x.range
  yr <- ggobj$layout$panel_params[[1]]$y.range

  #Several groups selection-----------------------------------------------------------------------
  continuing_selection<-readline("Do you want to select cell on the plot (Yes = 1, No = 0)?")
  continuing_selection <- as.numeric(unlist(strsplit(continuing_selection, ",")))
  selection_times = 1
  vgm_VDJ$selection_group<-rep(0,length(vgm_VDJ$barcode))
  while(continuing_selection == 1){
    readline ("Select cells on the plot (click on enter to access the plot)")

    # Variable for selected points -----------------------------------------
    selection <- data.frame(x = as.numeric(), y = as.numeric())
    colnames(selection) <- c(ggobj$plot$mapping$x, ggobj$plot$mapping$y)

    # Detect and move to plot area viewport---------------------------------
    plot.new()
    dev.new(width=1,height=1)#Define the dimension of the plot
    suppressWarnings(print(ggobj$plot))
    panels <- unlist(grid::current.vpTree()) %>%
      grep("panel", ., fixed = TRUE, value = TRUE)
    p_n <- length(panels)
    grid::seekViewport(panels, recording=TRUE)
    grid::pushViewport(grid::viewport(width=1, height=1,xscale = xr,yscale = yr))

    # Select point, plot, store and repeat----------------------------------
    for (i in 1:nbpoints){
      tmp <- grid::grid.locator('native')
      if (is.null(tmp)) break
      grid::grid.points(tmp$x,tmp$y, pch = 16, gp=grid::gpar(cex=0.5, col="darkred"))
      selection[i, ] <- as.numeric(tmp)
    }
    grid::grid.polygon(x= grid::unit(selection[,1], "native"), y= grid::unit(selection[,2], "native"), gp=grid::gpar(fill=NA))#to see the selection
    selection<-abs(selection)

    # Selected cells---------------------------------------------------------------------------------------------------------
    nb_cell<-length(vgm_VDJ$x)
    point_in_selection<-list()
    for (a in 1:nb_cell){
      point=select(vgm_VDJ, x, y)
      point=point[a,]
      odd=FALSE
      i=0
      j=nrow(selection)-1
      while (i<nrow(selection)-1) {
        i=i+1
        if(((selection[i,2]>point[2])!=(selection[j,2]>point[2])) && (point[1] < ((selection[j,1] - selection[i,1]) * (point[2] - selection[i,2]) / (selection[j,2] - selection[i,2])) + selection[i,1])){
          odd = !odd
        }
        j=i
      }
      result = odd
      if(result == FALSE){
        value=0
      }else{
        value=1
      }
      vgm_VDJ$point_in_selection[[a]] = value
      remove(point)
      remove(result)
      remove(value)
      remove(odd)
    }
    selected_cell_inside_polygon<-vgm_VDJ %>% filter(point_in_selection == 1)
    for (i in 1:length(selected_cell_inside_polygon$barcode)) {
      for (j in 1:length(vgm_VDJ$barcode)) {
        if (selected_cell_inside_polygon$barcode[[i]]== vgm_VDJ$barcode[[j]]){
          vgm_VDJ$selection_group[[j]]<-selection_times
        } else if (vgm_VDJ$selection_group[[j]]!=0){
          vgm_VDJ$selection_group[[j]]=vgm_VDJ$selection_group[[j]]
        } else if (vgm_VDJ$selection_group[[j]] ==0){
          vgm_VDJ$selection_group[[j]]=0
        }
      }
    }
    selection_times = selection_times+1
    continuing_selection<-readline("Other selection (Yes = 1, No = 0)?")
    continuing_selection <- as.numeric(unlist(strsplit(continuing_selection, ",")))
    remove(selection)
  }
  vgm_VDJ_just_selected_cells<-filter(vgm_VDJ, selection_group!=0)

  #Plot of selected groups-------------------------------------------------------------------
  p<-ggplot2::ggplot(data = vgm_VDJ_just_selected_cells, ggplot2::aes(x=x,y=y,fill=factor(selection_group))) +
    geom_spatial(data=images_tibble[1,], ggplot2::aes(grob=grob), x=0.5, y=0.5)+
    ggplot2::geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    ggplot2::coord_cartesian(expand=FALSE)+
    ggplot2::scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold", "#a65628", "#999999", "black", "grey", "white", "purple"))+
    ggplot2::xlim(0,max(bcs_merge %>%
                 filter(sample ==sample_names[1]) %>%
                 select(width)))+
    ggplot2::ylim(max(bcs_merge %>%
               filter(sample ==sample_names[1]) %>%
               select(height)),0)+
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::ggtitle(sample_names[1], title)+
    ggplot2::theme(axis.text=ggplot2::element_text(size=size),
          axis.title=ggplot2::element_text(size=size))+
    ggplot2::labs(fill = "Group")+
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size=3)))+
    ggplot2::theme_set(ggplot2::theme_bw(base_size = size))+
    ggplot2::theme(legend.key = ggplot2::element_rect(fill = "white"))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank())
  if(plotting=="TRUE"){
    return(list(p,vgm_VDJ))
  }
  if(plotting=="FALSE"){
    return(vgm_VDJ)
  }
}

