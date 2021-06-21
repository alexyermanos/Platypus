#' Helper function for GEX processing in VDJ_GEX_matrix. DO NOT RUN AS STANDALONE
#' @param ALL_PARAMS Provided by VDJ_GEX_matrix function.
#' @return GEX matrix
#' @export
#' @examples
#' \dontrun{
#' Do not run as standalone!
#' }

GEX_automate_single <- function(GEX.list,
                                GEX.integrate,
                                integration.method,
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
                                group.id){
  
  if(missing(GEX.integrate)) GEX.integrate <- T
  
  if(missing(mito.filter)) mito.filter <- 5
  if(missing(VDJ.gene.filter)) VDJ.gene.filter <- T
  if(missing(norm.scale.factor)) norm.scale.factor <- 10000
  if(missing(n.count.rna.min)) n.count.rna.min <- 0
  if(missing(n.count.rna.max)) n.count.rna.max <- Inf
  if(missing(n.feature.rna)) n.feature.rna <- 0
  if(missing(integration.method)) integration.method <- "scale.data"
  if(integration.method=="harmony") require(harmony)
  if(missing(GEX.integrate)) GEX.integrate <- T
  if(missing(n.variable.features)) n.variable.features <- 2000
  if(missing(cluster.resolution)) cluster.resolution <- .5
  if(missing(neighbor.dim)) neighbor.dim <- 1:10
  if(missing(mds.dim)) mds.dim <- 1:10
  
  for(i in 1:length(GEX.list)){
    #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes)
    GEX.list[[i]][["orig_barcode"]] <- as.character(gsub(".*_","",colnames(GEX.list[[i]])))
    GEX.list[[i]][["orig_barcode"]] <- gsub(GEX.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
  }
  
  if(GEX.integrate == T & length(GEX.list) > 1){ #combine all GEX into one seurat object and add s%number%_ to the FRONT of the barcode
    GEX.merged <- GEX.list[[1]]
    GEX.merged <- RenameCells(GEX.merged, new.names = paste0("s",1,"_",colnames(GEX.merged)))
    GEX.merged@meta.data$sample_id <- paste0("s",1)
    GEX.merged@meta.data$group_id <- group.id[1]
    for(i in 2:length(GEX.list)){
      GEX.list[[i]] <- RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",colnames(GEX.list[[i]])))
      GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
      GEX.list[[i]]@meta.data$group_id <- group.id[i]
      GEX.merged <- merge(GEX.merged, y = GEX.list[[i]], add.cell.ids = c("",""))
    }
    
    GEX.list <- list() #making this into a list item to make the downstream process uniform
    GEX.list[[1]] <- GEX.merged
  } else {
    for(i in 1:length(GEX.list)){ #or do not integrate, but still add the sample identifier to the FRONT of the barcode. This is to make the output uniform and to deal with the possibility of integrating VDJ but not GEX
      GEX.list[[i]] <- RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",colnames(GEX.list[[i]])))
      #add sample and group ID
      GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
      GEX.list[[i]]@meta.data$group_id <- group.id[i]
    }
  }
  
  
  for(i in 1:length(GEX.list)){
    
    
    holding_upper_gene_names <- toupper(rownames(GEX.list[[i]]))
    if(VDJ.gene.filter==T){
      antibody_gene_indices <- which(grepl((holding_upper_gene_names),pattern = "^IGHA")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHG")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHM")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHD")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHE")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGK")==F &
                                       grepl((holding_upper_gene_names),pattern = "^IGHV")==F &
                                       grepl((holding_upper_gene_names),pattern = "^JCHAIN")==F&
                                       grepl((holding_upper_gene_names),pattern = "^IGL")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRAV")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRAC")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRBC")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRGC")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRDC")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRBD")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRBJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRGV")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRDV")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRDD")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRDJ")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TRBV")==F &
                                       grepl((holding_upper_gene_names),pattern = "^TCRG")==F & #New filter for TCR gamma and delta chains
                                       grepl((holding_upper_gene_names),pattern = "^TCRD")==F)
      
      GEX.list[[i]] <- GEX.list[[i]][antibody_gene_indices,]
    }
    
    GEX.list[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^MT-") + Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^mt-")
    cell.subset.bool <- (GEX.list[[i]]$percent.mt< mito.filter & GEX.list[[i]]$nFeature_RNA >n.feature.rna & GEX.list[[i]]$nCount_RNA > n.count.rna.min & GEX.list[[i]]$nCount_RNA < n.count.rna.max)
    GEX.list[[i]] <- subset(GEX.list[[i]],cells=which(cell.subset.bool==T))
    
    if(integration.method=="sct"){
      GEX.list[[i]] <- Seurat::SCTransform(GEX.list[[i]],vars.to.regress = "percent.mt")
      GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature=Seurat::VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T)
      GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution)
      GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim)
      GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,check_duplicates=F)
      
    }
    if(integration.method=="harmony"){
      GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
      GEX.list[[i]] <- Seurat::FindVariableFeatures(GEX.list[[i]], selection.method = "vst", nfeatures = n.variable.features)
      
      all.genes <- rownames(GEX.list[[i]])
      GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature=VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- harmony::RunHarmony(GEX.list[[i]], "sample_id")
      GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T,reduction = "harmony")
      GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution,reduction = "harmony")
      GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim,reduction = "harmony")
      GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,reduction = "harmony",check_duplicates=F)
      
    }
    if(integration.method=="scale.data"){
      GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
      GEX.list[[i]] <- Seurat::FindVariableFeatures(GEX.list[[i]], selection.method = "vst", nfeatures = n.variable.features)
      all.genes <- rownames(GEX.list[[i]])
      
      GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = all.genes)
      GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]], features = Seurat::VariableFeatures(object = GEX.list[[i]]))
      GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]], dims = neighbor.dim)
      GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]], resolution = cluster.resolution)
      GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim)
      
      GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,check_duplicates=F)
    }
    
  }
  
  return(GEX.list)
}
