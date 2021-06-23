#' Automates the transcriptional analysis of the gene expression libraries from cellranger. This function will integrate multiple samples
#' @param GEX.outs.directory.list The path to the output of cellranger vdj runs. Multiple repertoires to be integrated together should be supplied as a character vector in the first element of a list. For example, if two separate VDJ repertoires should be integrated together (e.g. on the same tSNE plot), GEX.outs.directory.list[[1]] <- c("my.VDJ1.path/outs/","my.VDJ2.path/outs/") should be stored as input. If these repertoires should be analyzed separately, >GEX.outs.directory.list[[1]] <-  "my.VDJ1.path/outs/"  >GEX.outs.directory.list[[2]] <-  "my.VDJ2.path/outs/"  should be supplied..
#'  This can be left blank if supplying the clonotypes and all_contig files diretly as input. Multiple analyses can be stored
#' @param GEX.list List containing the output from Seurat Read10x. This must be supplied if GEX.out.directory is not provided.
#' @param integration.method String specifying which data normalization and integration pipeline should be used. Default is "scale.data", which correspondings to the ScaleData function internal to harmony package. 'sct'specifies SCTransform from the Seurat package. "harmony" should be specificied to perform harmony integration. This method requires the harmony package from bioconductor.
#' @param VDJ.gene.filter Logical indicating if variable genes from the b cell receprot and t cell receptor should be removed from the analysis. True is highly recommended to avoid clonal families clustering together.
#' @param mito.filter Numeric specifying which percent of genes are allowed to be composed of mitochondrial genes. This value may require visual inspection and can be specific to each sequencing experiment. Users can visualize the percentage of genes corresponding to mitochondrial genes using the function "investigate_mitochondial_genes".
#' @param norm.scale.factor Scaling factor for the standard Seurat pipeline. Default is set to 10000 as reported in Seurat documentation.
#' @param n.feature.rna Numeric that specifies which cells should be filtered out due to low number of detected genes. Default is set to 0. Seurat standard pipeline uses 2000.
#' @param n.count.rna.min Numeric that specifies which cells should be filtered out due to low RNA count.Default is set to 0. Seurat standard pipeline without VDJ information uses 200.
#' @param n.count.rna.max Numeric that specifies which cells should be filtered out due to high RNA count.Default is set to infinity. Seurat standard pipeline without VDJ information uses 2500.
#' @param n.variable.features Numeric specifying the number of variable features. Default set to 2000 as specified in Seurat standard pipeline.
#' @param cluster.resolution Numeric specifying the resolution that will be supplied to Seurat's FindClusters function. Default is set to 0.5. Increasing this number will increase the number of distinct Seurat clusters. Suggested to examine multiple parameters to ensure gene signatures differentiating clusters remains constant.
#' @param neighbor.dim Numeric vector specifying which dimensions should be supplied in the FindNeighbors function from Seurat. Default input is '1:10'.
#' @param mds.dim Numeric vector specifying which dimensions should be supplied into dimensional reduction techniques in Seurat and Harmony. Default input is '1:10'.
#' @param groups Integer specifying the groups of the different samples. This is needed if there are multiple biological replicates for a given condition sequenced and aligned through cellranger separately.
#' @return Returns a processed Seurat object containing transcriptional information from all samples which can be supplied to the VDJ_GEX_integrate function.
#' @export
#' @examples
#' \dontrun{
#' automate_GEX(out_directory=fullRepertoire.output,
#'  rep.size=3*length(unlist(fullRepertoire.output[[1]])),
#'  distribution="identical",
#'  with.germline="FALSE")
#'}
#'
#'
automate_GEX <- function(GEX.outs.directory.list,
                         GEX.list,integration.method,
                         VDJ.gene.filter,
                         mito.filter,
                         norm.scale.factor,
                         n.feature.rna,
                         n.count.rna.min,
                         n.count.rna.max,
                         n.variable.features,
                         cluster.resolution,
                         neighbor.dim,
                         mds.dim,
                         groups){
  print("This may take longer than other repertoire associated functions. Please see Seurat vignettes for further information")
  if(missing(integration.method)) integration.method <- "scale.data"
  if(integration.method=="harmony") require(harmony)
  if(missing(GEX.outs.directory.list)) print("Missing output directory of cellranger count. Assuming a list of 10x gene expression libraries is supplied as input.")
  if(missing(mito.filter)) mito.filter <- 5
  if(missing(integration.method)) integration.method <- "scale.data"

  if(missing(VDJ.gene.filter)) VDJ.gene.filter <- T
  if(missing(norm.scale.factor)) norm.scale.factor <- 10000
  if(missing(n.count.rna.min)) n.count.rna.min <- 0
  if(missing(n.count.rna.max)) n.count.rna.max <- Inf

  if(missing(n.feature.rna)) n.feature.rna <- 0
  if(missing(n.variable.features)) n.variable.features <- 2000
  if(missing(cluster.resolution)) cluster.resolution <- .5
  if(missing(neighbor.dim)) neighbor.dim <- 1:10
  if(missing(mds.dim)) mds.dim <- 1:10
  if(missing(integration.method)) integration.method <- "scale.data"


  if(missing(GEX.list) & !missing(GEX.outs.directory.list)){
    GEX.list <- list()
    for(i in 1:length(GEX.outs.directory.list)){
      directory_read10x <- Seurat::Read10X(data.dir=GEX.outs.directory.list[[i]])
      rownames(directory_read10x) <- toupper(rownames(directory_read10x))
      GEX.list[[i]] <- Seurat::CreateSeuratObject(directory_read10x)
    }
  }
  for(i in 1:length(GEX.list)){
    if(missing(groups)) groups <- 1:length(GEX.outs.directory.list[[i]])


    if(length(GEX.outs.directory.list[[i]])>1){
      sample1 <- which(grepl(x = names(GEX.list[[i]]$orig.ident),pattern = "\\_")==F)
      GEX.list[[i]]$sample_id <- rep("",ncol(GEX.list[[i]]))
      GEX.list[[i]]$sample_id[sample1] <- i

      sample_rest <- (gsub("\\_.*","",names(GEX.list[[i]]$orig.ident)))
      GEX.list[[i]]$sample_id[(length(sample1)+1):length(sample_rest)] <- as.integer(sample_rest[-c(1:length(sample1))])
      GEX.list[[i]]$group_id <- rep(groups,table(GEX.list[[i]]$sample_id))
    }
    else{
      GEX.list[[i]]$sample_id <- 1
      GEX.list[[i]]$group_id <- 1
    }

    holding_upper_gene_names <- toupper(rownames(GEX.list[[i]]))
    if(VDJ.gene.filter==T){
      antibody_gene_indices <- which(grepl((holding_upper_gene_names),pattern = "IGHA")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHG")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHM")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHD")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHE")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGK")==F &
                                       grepl((holding_upper_gene_names),pattern = "IGHV")==F &
                                       grepl((holding_upper_gene_names),pattern = "JCHAIN")==F&
                                       grepl((holding_upper_gene_names),pattern = "IGL")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRAV")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRAC")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRBC")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRGC")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRDC")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRBD")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRBJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRGV")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRGJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRGJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRDV")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRDD")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRDJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRBV")==F &
                                       grepl((holding_upper_gene_names),pattern = "TRAJ")==F)
      GEX.list[[i]] <- GEX.list[[i]][antibody_gene_indices,]
    }

    GEX.list[[i]][["barcode"]] <- as.character(gsub(".*_","",colnames(GEX.list[[i]])))
    GEX.list[[i]][["barcode"]] <- gsub(GEX.list[[i]]$barcode,pattern = "-1",replacement = "")

    GEX.list[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^MT-") + Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^mt-")
    cell.subset.bool <- (GEX.list[[i]]$percent.mt<20 & GEX.list[[i]]$nFeature_RNA >n.feature.rna & GEX.list[[i]]$nCount_RNA > n.count.rna.min & GEX.list[[i]]$nCount_RNA < n.count.rna.max)
    GEX.list[[i]] <- subset(GEX.list[[i]],cells=which(cell.subset.bool==T))#  & nFeature_RNA>n.feature.rna  & nCount_RNA > n.count.rna)) #subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

    if(integration.method=="sct"){
      GEX.list[[i]] <- Seurat::SCTransform(GEX.list[[i]],vars.to.regress = "percent.mt")
      GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature=Seurat::VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T)
      GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution)

    }
    if(integration.method=="harmony"){
      GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
      all.genes <- rownames(GEX.list[[i]])
      GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = Seurat::VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature=Seurat::VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- harmony::RunHarmony(GEX.list[[i]], "sample_id")
      GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T)
      GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution)

    }
    if(integration.method=="scale.data"){
    GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
    GEX.list[[i]] <- Seurat::FindVariableFeatures(GEX.list[[i]], selection.method = "vst", nfeatures = n.variable.features)
    all.genes <- rownames(GEX.list[[i]])

    GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = all.genes)
    GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]], features = Seurat::VariableFeatures(object = GEX.list[[i]]))
    GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]], dims = neighbor.dim)
    GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]], resolution = cluster.resolution)
    }
    GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim)
    GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim)

  }
  return(GEX.list)
}
