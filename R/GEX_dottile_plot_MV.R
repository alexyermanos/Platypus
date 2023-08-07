#'GEX Dottile plots
#'
#'@description Outputs a dotplot for gene expression, where the color of each dot is scaled by the gene expression level and the size is scaled by the \% of cells positive for the gene
#' @param GEX GEX seurat object generated with VDJ_GEX_matrix
#' @param genes Character vector. Genes of those in rownames(GEX) to plot. Can be any number, but more then 30 is discuraged because of cluttering
#' @param group.by Character. Name of a column in GEX@meta.data to split the plot by. If set to \"none\", a plot with a single column will be produced.
#' @param threshold.to.plot Integer 1-100. \% of cells which must be expressing the feature to plot a point. If below, the field will be left empty
#' @param gene_module If missing, the dotplot displays the single genes on the y-axis that are listed in the parameter genes. If gene_module == mathew, the dottile plot shows cell types on the y-axis (Mathew et al., 2021)
#' @param platypus.version This is coded for \"v3\" only, but in practice any Seurat Object can be fed in
#' @return Returns a ggplot object were the dot size indicates the percentage of expressing cells and the dot color indicates the expression level.
#' @export
#' @examples
#' \dontrun{
#' #To return a plot detailing the expression of common genes by seurat cluster
#'GEX_dottile_plot(GEX = Platypus::small_vgm[[2]], genes = c("CD19","CD83"),
#'group.by = "seurat_clusters", threshold.to.plot = 5)
#'
#'
#'# EXAMPLE CODE for getting dottile plot for B cell types
#'
#'
#'dottile <- GEX_dottile_plot_MV(GEX = VGM[[2]], genes = c('CD19'), group.by = 'seurat_clusters', threshold.to.plot = 5, gene_module = 'mathew')
#'}
GEX_dottile_plot_MV <- function(GEX,
                             genes,
                             group.by,
                             threshold.to.plot,
                             gene_module,
                             platypus.version){

  group <- NULL
  name <- NULL
  value <- NULL
  mean_scaled_expression <- NULL
  perc_expressing_cells <- NULL

  platypus.version <- "v3"

  if(missing(GEX)) stop("Please provide a seurat object as input to GEX")
  if(missing(threshold.to.plot)) threshold.to.plot <- 5

  if(missing(genes)){
    if(GEX$celltype[1] == "T cell"){
      genes <- c("CD4","CD8A","CD28","ICOS","CD40LG","PDCD1","LAG3","S1PR1","PTPRC","SELL","TBX21","GATA3","SPI1","IRF4","BCL6","STAT4","STAT6","FOXP3", "MKI67","TCF7","KLRG1","TOX", "GZMB","IFNG","TGFB1","IL10","IL17A","CSF2","IL2RA","IL4RA","IL12RB1","CXCR3","CXCR5","CCR5","CCR8","SELL","ITGA4","ITGB2")
    } else {
      genes <- c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34","CD1D1","CR2","MS4A1","CXCR5","SELL","CD40","CD83","H2-AB1","H2-EB1","CD27","POU2AF1","NT5E","FAS","PDCD1LG2","PRDM1","ITGAM","IL10","IL12A","HAVCR2")
    }
  }
  if (missing(gene_module)) {gene_module <- '' }

  # dotplot with genes on y-axis if gene_module is missing

  if (gene_module == '') {

  #Make unique
  genes <- unique(genes)

  to_del <- c()
  for(i in 1:length(genes)){
    if(!genes[i] %in% rownames(GEX)){
      warning(paste0(genes[i], " not found in seurat object. This gene is skipped"))
      to_del <- c(to_del, i)
    }
  }
  if(length(to_del) > 0){
    genes <- genes[-to_del]
  }

  if(missing(group.by)) group.by <- "none"
  if(!group.by %in% names(GEX@meta.data)){
    warning("group.by column not found. Returning plot with single column")
    group.by <- "singlegroup"
    GEX$singlegroup <- "group"
  }

  to_plot_f <- SeuratObject::FetchData(GEX, vars = c(group.by, genes))

  #pivot
  to_plot_f <- tidyr::pivot_longer(to_plot_f, cols = c(2:ncol(to_plot_f)))
  names(to_plot_f)[1] <- "group"
  #summarize and group
  #for now by cluster
  to_plot_sum_f <- to_plot_f %>% dplyr::group_by(group,name) %>% dplyr::summarise(mean_scaled_expression = mean(value[value > 0]), perc_expressing_cells = (length(value[value > 0])/dplyr::n()*100))

  to_plot_sum_f$name <- ordered(as.factor(to_plot_sum_f$name), levels = rev(genes))

  to_plot_sum_f$perc_expressing_cells[to_plot_sum_f$perc_expressing_cells < threshold.to.plot] <- NA

  plot_out <- ggplot2::ggplot(to_plot_sum_f, ggplot2::aes(x = group, y = name, col = mean_scaled_expression, size = perc_expressing_cells)) + ggplot2::geom_point(show.legend = T) + cowplot::theme_cowplot()  + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1), legend.position = "right",axis.text.x = ggplot2::element_text(angle = 60, vjust = 0.95, hjust=1)) + ggplot2::labs(title = paste0("Expression by ", group.by), x = "", y = "", color = "Scaled expression", size = "% of expressing cells")  + ggplot2::scale_color_viridis_c(option = "B", end = 0.9) + ggplot2::scale_size_binned(range = c(1,9.5))

    return(plot_out)

  }

  else if (gene_module == 'mathew') {

    # setup gene list
    gene_list <- list()
    # NaiveBcell
    #gene_list[[length(gene_list)+1]] <- ("IGHD")
    #names(gene_list)[length(gene_list)] <- "naive"
    # MarginalZone
    gene_list[[length(gene_list)+1]] <- (c("CD1D1", "CD9"))
    names(gene_list)[length(gene_list)] <- "MZ"
    # GerminalCenter
    gene_list[[length(gene_list)+1]] <- (c("AICDA", "FAS"))
    names(gene_list)[length(gene_list)] <- "GC"
    # MemoryCells
    gene_list[[length(gene_list)+1]] <- (c("PML", "PHF21A", "RUNX3", "KLF9", "ZBTB20", "KMT2C", "SKI", "TCF4", "BMI1", "CD80", "IL6RA", "ITGB7", "CCR6", "PECAM1", "NID1", "RGMB", "ACKR3", "CDH17", "ACKR2", "ITGA4", "NAIP2", "BIRC6", "BCL2", "TRAF1", "MALT1", "CCND3", "RB1"))
    names(gene_list)[length(gene_list)] <- "Bmem"
    # PreMemoryCells
    gene_list[[length(gene_list)+1]] <- (c("PBX3", "NCOA1", "ETV6", "EFNB1", "S1PR3", "BMPR1A", "CD86", "CD59A"))
    names(gene_list)[length(gene_list)] <- "PreBmem"
    # naive like
    gene_list[[length(gene_list)+1]] <- (c("CD74", "CCR7", "EBF1", "BTG1", "CD79A", "FCMR", "PAX5", "TCF3", "SELL"))
    names(gene_list)[length(gene_list)] <- "naive_like"
    # DZ
    gene_list[[length(gene_list)+1]] <- (c("PIF1", "OTUB2", "PSRC1", "LRRC49", "CCNB2", "MDN1", "LGR5", "ADHFE1", "PAFAH1B3", "HMMR","TIFA", "RELN", "SELENOP", "CDKN3", "GCNT3", "CDC20", "UBE2C", "SERINC5", "CENPA", "CCNB1", "PLK1", "CENPF", "TGFA", "TPX2", "DEPDC1B", "GPD1L", "POLH", "CEP55", "NEK2", "ZKSCAN16", "LMO4", "GAS2L3", "PTGR1", "ASPM", "RNF125", "CXCR4", "KIF20A", "KIF2C", "RABGAP1L", "BUB1", "NEBL", "LIG4", "NEIL1", "PDE2A", "EHF", "RACGAP1", "CDCA8", "CDC14B", "CDC25C", "LIPC", "MPRIP", "AURKA", "GPSM2", "SPDYA", "AKAP12", "KIF22", "BIRC5", "PFN2", "DDX19A", "LDLRAD4", "GCSAM"))
    names(gene_list)[length(gene_list)] <- "DZ"
    # LZ
    gene_list[[length(gene_list)+1]] <- (c("EGR1", "EGR2", "CCND2", "CD69", "CD86", "CD40", "SLAMF1", "IFI30", "MYC","S1PR3", "IL4I1", "CCR6", "CD83", "GPR183", "NFKBIA", "PLSCR1", "SOCS3", "LIFR", "BCL2A1A", "BCL2A1D", "BCL2A1B", "HHEX", "DOCK10", "HIVEP3", "SAMSN1", "IRF4", "GIMAP4", "PTGER3",  "FSCN1", "MARCKS", "IER2", "SNN", "BHLHE40", "SLPI", "MATR3", "RAPGEF4", "RILPL2", "GADD45G", "CFP", "TYMS", "MDN1"))
    names(gene_list)[length(gene_list)] <- "LZ"
    # PlasmaCells
    gene_list[[length(gene_list)+1]] <- (c("BHLHA15", "SDC1", "PRDM1", "IRF4", "XBP1", "ATF6", "CREB3", "CD28", "LY6C1", "IL13RA1", "ITGB5", "SELPLG", "LGALS1"))
    names(gene_list)[length(gene_list)] <- "PC"
    # LZcells receiving Tfh cells signals to enter/exit the GC
    gene_list[[length(gene_list)+1]] <- (c("ITGB2", "CXCR5", "TGFBR2", "IL21R", "EBF1", "IRF8", "NFKB2", "ICOSL", "SLAMF1", "TNFRSF13C"))
    names(gene_list)[length(gene_list)] <- "TFH_signal"
    # Cell in the lightzone Reentering to DZ
    gene_list[[length(gene_list)+1]] <- (c("AICDA", "CXCR4", "XRCC4", "CD24A", "MKI67", "FOXO1", "TCF3"))
    names(gene_list)[length(gene_list)] <- "DZ_reentry"

    # Add gene modules to GEX object
    Seurat::Idents(GEX) <- "seurat_clusters"
    GEX <- Seurat::AddModuleScore(
      object = GEX,
      features = gene_list,
      name = paste0(names(gene_list), "_module"))
    #Plotting
    plot2 <- Seurat::DotPlot(GEX, features = paste0(names(gene_list), "_module", 1:length(gene_list))) + ggplot2::coord_flip() + ggplot2::RotatedAxis()
    plot_out <- plot2 + ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(angle = 45, hjust = 1, vjust = 0.1), legend.position = "bottom", legend.direction = "vertical")

    return(plot_out)
  }


}

