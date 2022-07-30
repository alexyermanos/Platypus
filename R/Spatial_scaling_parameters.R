#' Scaling modulation for gene expression on spatial image
#' @description Scaling of the spatial parameters to be able to express the gene expression on the spatial image.
#' @param vgm_spatial List containing the output of VDJ_GEX_matrix function from Platypus with at least the gene expression data and the addition of spatial parameters: image, scalefacor, tissue, cluster and matrix.
#' @param GEX.out.directory.list Path to the filtered feature bc matrix data.
#' @param sample_names Character vector containing the name of the sample.
#' @return Returns a list containing all parameters to scale the data on the spatial image. List element [[1]]: images_cl. List element [[2]]: height of the image. List element [[3]]: width of the image. List element [[4]]: grobs. List element [[5]]: images_tibble. List element [[6]]: scales. List element [[7]]: cluster. List element [[8]]: bcs. List element [[9]]: matrix. List element [[10]]: bcs_merge.
#' @export
#' @examples
#' \dontrun{
#' scaling_parameters<-Spatial_scaling_parameters(vgm_spatial = vgm_spatial,
#' GEX.out.directory.list = GEX.out.directory.list,
#' sample_names = sample_names)
#'}
Spatial_scaling_parameters <- function(vgm_spatial,
                                       GEX.out.directory.list,
                                       sample_names){

  if(missing(vgm_spatial)) stop("Please provide vgm_spatial input for this function")
  if(missing(GEX.out.directory.list)) stop("Please provide GEX.out.directory.list input for this function")
  if(missing(sample_names)) stop("Please provide sample_names input for this function")
  #if (!require("hdf5r", character.only = TRUE)) {
  #  install.packages("hdf5r")
  #  library(hdf5r)
  #}
  #if (!require("Matrix", character.only = TRUE)) {
  #  install.packages("Matrix")
  #  library(Matrix)
  #}
  platypus.version <- "v3"

  images_cl <- list()
  for (i in 1:length(sample_names)) {
    images_cl[[i]] <- vgm_spatial$spatial$image[[i]]
  }
  height <- list()
  for (i in 1:length(sample_names)) {
    height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
  }
  height <- dplyr::bind_rows(height)
  width <- list()
  for (i in 1:length(sample_names)) {
    width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
  }
  width <- dplyr::bind_rows(width)
  grobs <- list()
  for (i in 1:length(sample_names)) {
    grobs[[i]] <- grid::rasterGrob(images_cl[[i]], width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
  }
  images_tibble <- tibble::tibble(sample=factor(sample_names), grob=grobs)
  images_tibble$height <- height$height
  images_tibble$width <- width$width

  scales <- list()
  for (i in 1:length(sample_names)) {
    scales[[i]] <- vgm_spatial$spatial$scalefactor[[i]]
  }
  clusters <- list()
  for (i in 1:length(sample_names)) {
    clusters[[i]] <- vgm_spatial$spatial$cluster[[i]]
  }
  bcs <- list()
  for (i in 1:length(sample_names)) {
    bcs[[i]] <- vgm_spatial$spatial$tissue[[i]]
    bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef
    bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
    bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
    bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE)
    bcs[[i]]$height <- height$height[i]
    bcs[[i]]$width <- width$width[i]
  }
  names(bcs) <- sample_names
  matrix <- list()
  for (i in 1:length(sample_names)) {
    matrix[[i]] <- as.data.frame(t(Seurat::Read10X(GEX.out.directory.list[[i]])))
  }
  umi_sum <- list()
  for (i in 1:length(sample_names)) {
    umi_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                               sum_umi = Matrix::rowSums(matrix[[i]]))
  }
  names(umi_sum) <- sample_names
  umi_sum <- dplyr::bind_rows(umi_sum, .id = "sample")
  gene_sum <- list()
  for (i in 1:length(sample_names)) {
    gene_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                                sum_gene = Matrix::rowSums(matrix[[i]] != 0))
  }
  names(gene_sum) <- sample_names
  gene_sum <- dplyr::bind_rows(gene_sum, .id = "sample")
  bcs_merge <- dplyr::bind_rows(bcs, .id = "sample")
  bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))
  bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))
  data<-list()
  data$bcs_merge<-bcs_merge
  data$images_tibble<-images_tibble
  return(list(images_cl,height,width,grobs,images_tibble,scales,clusters,bcs,matrix,bcs_merge))
}
