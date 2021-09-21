#' !Only for platypus version v2. For platypus v3 refer to: VDJ_GEX_overlay_clones() Visualize selected clonotypes on the tSNE or UMAP projection.
#' @param GEX.list list of Seurat objects, output of the automate_GEX function.
#' @param VDJ.GEX.integrate.list Output of the VDJ_GEX_integrate function.
#' @param highlight.type (Optional) either "None" if representation highlighted by cluster, "clonotype" if want to highlight most expanded clonotypes, or "sample" if several samples are within the same Seurat object. Default is None.
#' @param highlight.number (Optional) an integer or list of integers representing the number of most expanded clonotypes or samples one wants to select eg 4 to highlight the 4th most expanded clonotype or 2:5 to highlight the top 2 to top 5 most expanded clonotype. Only compatible with highlight.type "clonotype" or "sample", will be ignored if type is "None". Default is 1.
#' @param reduction (Optional) Reduction to plot, either "tsne", "umap", or "harmony". Default is "tsne".
#' @return concatenated ggplot2 plot with selected clonotypes highlighted (if None, the coloring is according to the clustering).
#' @export
#' @examples
#' \dontrun{
#' GEX_visualize_clones(GEX.list=automate_GEX.output,
#'  VDJ.per.clone=VDJ_per_clone.output,
#'  highlight.type="clonotype",
#'  highlight.number=1:4,
#'  reduction="umap")
#' }
GEX_visualize_clones <- function(GEX.list, VDJ.GEX.integrate.list, highlight.type, highlight.number, reduction){

  # initialize plot list
  pl_list <- list()

  if(missing(GEX.list)) stop("No provided list of Seurat objects. Please provide the output of the automate_GEX function.")
  if(missing(VDJ.GEX.integrate.list)) stop("No provided per clone VDJ breakdown. Please provide the output of the VDJ_per_clone function.")
  if(missing(highlight.type)) highlight.type <- "None"
  if(missing(highlight.number)) highlight.number <- 1
  if(missing(reduction)) reduction <- "tsne"

  # check if reduction, highlight.type and highlight.number are correctly given
  if(!(highlight.type=="None") && !(highlight.type=="clonotype") &&
     !(highlight.type=="sample")){
    stop("highlight.type can either be clonotype, sample or None")
  }
  if(!(mode(highlight.number)=="numeric")){
    stop("highlight.number can either be an integer or a vector of integers")
  }
  if(!(reduction=="tsne" || reduction=="umap" || reduction=="harmony" )){
    stop("reduction must be tsne or umap")
  }

  if(highlight.type=="None"){
    if(reduction=="tsne"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "tsne", pt.size = 1) +
        ggplot2::ggtitle("tSNE of all repertoires") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    else if(reduction=="umap"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "umap", pt.size = 1) +
        ggplot2::ggtitle("UMAP of all repertoires") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    else if(reduction=="harmony"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "harmony", pt.size = 1) +
        ggplot2::ggtitle("Harmony of all repertoires") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }
  else if(highlight.type=="clonotype"){
    for(i in 1:length(VDJ.GEX.integrate.list)){
      ## now need indices of the cells in the GEX object.
      holding_cells <- list()
      holding_all_cells_integer <- list()
      for(j in 1:length(highlight.number)){
        holding_cells[[j]] <- toString(VDJ.GEX.integrate.list[[i]]$cell_index[highlight.number[j]])
        holding_all_cells_integer[[j]] <- as.integer(stringr::str_split(holding_cells[[j]],pattern = ";")[[1]])

      }

      #holding_all_cells <- toString(VDJ.GEX.integrate.list[[i]]$cell_index[highlight.number])
      #holding_all_cells_integer <- as.integer(str_split(holding_all_cells,pattern = ";")[[1]])
      if(reduction=="tsne"){
        pl_list[[i]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "tsne", pt.size = 1,cells.highlight = (holding_all_cells_integer),cols.highlight = grDevices::rainbow(length(holding_all_cells_integer))) + ggplot2::ggtitle(paste("tSNE highlighting clonotypes from sample ",i,sep=""),)  +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      }
      else if(reduction=="umap"){
        pl_list[[i]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "umap", pt.size = 1,cells.highlight = (holding_all_cells_integer),cols.highlight = grDevices::rainbow(length(holding_all_cells_integer)),order=as.list(c(paste0("Group_",c(max(highlight.number):min(highlight.number))),"Unselected"))) +
          ggplot2::ggtitle(paste("UMAP highlighting clonotypes from sample ",i,sep="")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      }
      else if(reduction=="harmony"){
        pl_list[[i]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "harmony", pt.size = 1,cells.highlight = (holding_all_cells_integer),cols.highlight = grDevices::rainbow(length(holding_all_cells_integer))) +
          ggplot2::ggtitle(paste("Harmony highlighting clonotypes from sample ",i,sep="")) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      }
    }
  }
  else if(highlight.type=="sample"){
    if(reduction=="tsne"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "tsne", pt.size = 1,cells.highlight = which(GEX.list[[1]]$sample_id==i)) +
        ggplot2::ggtitle(paste("tSNE highlighting sample ",i,sep="")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    else if(reduction=="umap"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "umap", pt.size = 1,cells.highlight = which(GEX.list[[1]]$sample_id==i)) +
        ggplot2::ggtitle(paste("UMAP highlighting sample ",i,sep="")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
    else if(reduction=="harmony"){
      pl_list[[1]] <- Seurat::DimPlot(GEX.list[[1]], reduction = "harmony", pt.size = 1,cells.highlight = which(GEX.list[[1]]$sample_id==i)) +
        ggplot2::ggtitle(paste("Harmony highlighting sample ",i,sep="")) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }
  return(pl_list)
}
