## ----setup, include=FALSE, eval = FALSE---------------------------------------
#  #to avoid of the dplyr summarise warnings
#  library(dplyr, warn.conflicts = FALSE)
#  options(dplyr.summarise.inform = FALSE)
#  library(stats, warn.conflicts = F)
#  
#  knitr::opts_chunk$set(fig.width=7, fig.height=7)

## ---- eval = FALSE, include=FALSE---------------------------------------------
#  #this is for sourcing all functions under development and is not included in the knitting
#  source_code_dir <- "C:/Users/vickr/Documents/GitHub/Platypus/vignettes"
#  file_path_vec <- list.files(source_code_dir, full.names = T)
#  for(f_path in file_path_vec){
#    print(f_path)
#    tryCatch({source(f_path)}, error = function(e){e})
#  }
#  gc()
#  knitr::opts_chunk$set(fig.width=6, fig.height=6)
#  library(tidyverse)
#  library(Seurat)

## ---- fig.show='hold', message=FALSE, eval = FALSE----------------------------
#  
#  ### Removing any previous versions of the package
#  #First can ensure that there is no previous version installed locally
#  #detach("package:Platypus", unload=TRUE)
#  #remove.packages("Platypus")
#  
#  ### Packages most frequently imported
#  #install.packages("Tidyverse")
#  #install.packages("Biostrings")
#  #install.packages("jsonlite")
#  #install.packages("seqinr")
#  #install.packages("Seurat")
#  
#  ### Downloading and installing Platypus
#  
#  # First we need to download the most recent version from the master branch at https://github.com/alexyermanos/Platypus we can install the package using the following command.
#  # WARNING: This needs to be replaced with your own directory where the downloaded package is found
#  
#  # For MacOS users it may look like this
#  #install.packages("~/Downloads/Platypus_3.1.tar.gz", repos = NULL, type="source")
#  
#  # For windows it will likely look something like this.
#  # WARNING: You will need to replace 'YourPCName' with your user name for the windows account in the directory.
#  #install.packages("C:/Users/YourPCName/Downloads/Platypus_3.1.tar.gz", repos = NULL, type="source")
#  
#  # Now we can load the installed package into the R environment. In case of problems with installing other R packages that are used in Platypus, please see the README file at the https://github.com/alexyermanos/Platypus, where we outline how to install the other R packages for both Windows and MacOS.
#  library(Platypus)
#  
#  # The individual R functions can additionally be found on the github in the Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment in case not all functions are desired. These functions are actively updated and may include more features than the in original tar.gz file.
#  

## ---- fig.show='hold', message=FALSE,results = 'hide', eval = FALSE-----------
#  
#  ### Downloading the test data for VDJ_GEX_matrix
#  # While the Platypus manuscript uses the COVID-19 data, the vignette for Platypus v3 will use the data from murine B cells in the aged CNS, which can be found at the following link https://polybox.ethz.ch/index.php/s/fxQJ3NrRSwiPSSo This small dataset contains VDJ (separate libraries for B) and GEX libraries from the central nervous system of two murine samples. More information can be found https://doi.org/10.1098/rspb.2020.2793
#  
#  # After downloading the zip file named "Platypus_CNS_data.zip", please unzip the file and find the path to the newly formed folder. Typically this will be in the Downloads folder, so the code below should work on MacOS. For windows please uncomment the code and change the user name to match your PC.
#  
#  VDJ.out.directory.list <- list() ### Set directory to the outs folder of cellranger vdj
#  VDJ.out.directory.list[[1]] <- c("~/Downloads/Platypus_CNS_data/VDJ_S1/")
#  VDJ.out.directory.list[[2]] <- c("~/Downloads/Platypus_CNS_data/VDJ_S2/")
#  
#  GEX.out.directory.list <- list() ### Set directory to the outs folder of cellranger count
#  GEX.out.directory.list[[1]] <- c("~/Downloads/Platypus_CNS_data/GEX_S1/")
#  GEX.out.directory.list[[2]] <- c("~/Downloads/Platypus_CNS_data/GEX_S2/")
#  
#  # For windows:
#  #VDJ.out.directory.list[[1]] <- c("C:/Users/YourPCName/Downloads/PlatypusTestData/Patient1_GEX")
#  #VDJ.out.directory.list[[2]] <- c("C:/Users/YourPCName/Downloads/PlatypusTestData/Patient2_GEX")
#  
#  # We will call the output vgm (short for Vdj_Gex_Matrix) - this object can be supplied as input to downstream functions in v3 of the package.
#  vgm <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#                                 GEX.out.directory.list = GEX.out.directory.list,
#                                 GEX.integrate = T,
#                                 VDJ.combine = T,
#                                 integrate.GEX.to.VDJ = T,
#                                 integrate.VDJ.to.GEX = T, #This will adjunct the VDJ information as metadata to the GEX object
#                                 exclude.GEX.not.in.VDJ = F,
#                                 filter.overlapping.barcodes.GEX = T,
#                                 filter.overlapping.barcodes.VDJ = T,
#                                 exclude.on.cell.state.markers = c("CD3E"), #Exclude T cells from this analysis
#                                 get.VDJ.stats = T,
#                                 parallel.processing = "none", #see note at the end of this chunk
#                                 trim.and.align = F, #Do not align BCR sequences to reference
#                                 group.id = c(1,2))
#  
#  ## The output will be a list -
#  # vgm[[1]] corresponds to the VDJ master dataframe
#  # vgm[[2]] corresponds to the GEX in the form of a seurat object
#  # vgm[[3]] corresponds to the output of VDJ_stats subfunction - which provides information about the number of chains, sequencing reads, etc
#  # vgm[[4]] holds the input parameters
#  # vgm[[5]] holds the sessionInfo output at the time of function call
#  
#  ## Setting trim.and.align to TRUE will provide full-length VDJ and VJ sequences but also increase run time significantly. Alignment is done via Biostrings::pairwiseAlignment(). Gap opening and extension costs can be adapted
#  
#  ## This function can similarly be used when only VDJ or GEX data is present. Simply do only provide the path list of either GEX or VDJ
#  # VDJ_comb_gex <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#  #                                GEX.integrate = F,
#  #                                VDJ.combine = T,
#  #                                integrate.GEX.to.VDJ = F,
#  #                                integrate.VDJ.to.GEX = F,
#  #                                exclude.GEX.not.in.VDJ = F,
#  #                                filter.overlapping.barcodes.GEX = F,
#  #                                filter.overlapping.barcodes.VDJ = F,
#  #                                get.VDJ.stats = F,
#  #                                parallel.processing = "none",
#  #                                trim.and.align = T,
#  #                                gap.opening.cost = 10,
#  #                                gap.extension.cost = 4
#  #                                group.id = c(1,2))
#  

## ----eval = FALSE, fig.show='hold', message=FALSE-----------------------------
#  
#  head(colnames(vgm[[1]]))
#  ## By setting integrate.GEX.to.VDJ and integrate.VDJ.to.GEX to T, VDJ and GEX information will be found in vgm[[1]] and vgm[[2]] objects.
#  
#  # For example, the seurat-determined cluster is attached to each cell in the VDJ library by
#  head(vgm[[1]]$seurat_clusters)
#  
#  # which corresponds to cells with the following VDJ_cdr3
#  head(vgm[[1]]$VDJ_cdr3s_aa)
#  # an NA indicates that the cell barcode in the VDJ library was not detected in the GEX object (or was filtered out, depending on mitochondrial gene limits, etc)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  print(vgm[[3]])
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  print(vgm[[4]]) #Runtime params
#  print(vgm[[5]]) #session info
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  # In Platypus version 2, the output from GEX_automate was used as input to other GEX functions. These functions are still compatible with v3 if the vgm[[2]] seurat object is supplied as input.
#  # For example, the following function can be used to calculate the DE genes for each cluster, as before.
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", group.by = "sample_id")
#  
#  Seurat::DimPlot(vgm[[2]],reduction = "pca", group.by = "sample_id")
#  
#  Seurat::DimPlot(vgm[[2]],reduction = "tsne", group.by = "sample_id")
#  
#  #alternatively plot each sample separately
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", split.by = "sample_id")
#  
#  #this also works with any other column of vgm[[2]]@meta.data
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", split.by = "group_id")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  GEX_proportions_barplot(GEX = vgm[[2]], stacked.plot = T, source.group = "sample_id", target.group = "seurat_clusters")
#  #This function is very flexible and can be used to plot proportions of cells from and of any groups. For this use the source.group and target.group parameters to specify metadata columns.

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  Seurat::FeaturePlot(vgm[[2]],reduction = "umap",features = c("CD19","PTPRC", "EBF1", "H2-K1"))
#  
#  #To easily scout through genes in the dataset use:
#  #View(as.data.frame(rownames(vgm[[2]])))
#  

## ---- results='hide', eval = FALSE--------------------------------------------
#  #using defaults
#  vgm[[2]] <- GEX_phenotype(vgm[[2]], default = T)
#  
#  #custom criteria
#  #vgm[[2]] <- GEX_phenotype(vgm[[2]], default = F,cell.state.markers=c("CD8A+;CCL5+;CD44+;IL7R-;CD19-","CD8A+;CCL5-;CD44+;IL7R+;CD19-"),cell.state.names=c("EffectorCD8","MemoryCD8"))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", group.by = "cell.state")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  gene_expression_cluster <- GEX_cluster_genes(vgm[[2]],min.pct = 0.25)
#  
#  length(gene_expression_cluster) # length of 8, corresponding to 8 clusters
#  length(unique(vgm[[2]]$seurat_clusters)) # length of 8
#  
#  print(sapply(gene_expression_cluster,nrow)) #Nr of differentially expressed genes per cluster

## ---- eval = FALSE------------------------------------------------------------
#  head(gene_expression_cluster[[1]])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  agedCNS_heatmap_clusters <- GEX_cluster_genes_heatmap(GEX = vgm[[2]], GEX_cluster_genes.output = gene_expression_cluster,n.genes.per.cluster = 3,max.cell = 30,metric = "avg_logFC", platypus.version = "v3")
#  
#  print(agedCNS_heatmap_clusters)

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  agedCNS_heatmap_volcano <- GEX_volcano(DEGs.input = gene_expression_cluster, input.type = "cluster.genes", RP.MT.filter = T, color.p.threshold = 0.01, n.label.up = 10, n.label.down = 10)
#  
#  print(agedCNS_heatmap_volcano[[1]]) #genes specific to cluster 0

## ----results='hide', eval = FALSE---------------------------------------------
#  
#  ontology_agedCNS <- GEX_GOterm(GEX.cluster.genes.output = gene_expression_cluster, topNgenes = 10, go.plots = F)
#  head(ontology_agedCNS[[1]][[1]]) #Cluster 0
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  top_10_genes_per_cluster <- GEX_topN_DE_genes_per_cluster(GEX_cluster_genes.output = gene_expression_cluster, n.genes = 10, by_FC = T)
#  head(top_10_genes_per_cluster)
#  

## ----GSEA, eval = FALSE, warning=FALSE----------------------------------------
#  
#  gsea_EAE <- GEX_GSEA(GEX.cluster.genes.output = gene_expression_cluster[[1]], MT.Rb.filter = T, path.to.pathways = "~/Downloads/c7.all.v7.4.symbols.gmt")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  DE_genes_samples <- GEX_DEgenes(GEX = vgm[[2]],min.pct = .25, grouping.column = "sample_id",group1 = "s1", group2 = "s2",return.plot = "volcano",up.genes = 5,down.genes = 5,logFC = F)
#  #This function is flexible and takes any column name as grouping.column to allow easy exploration of differences between custom groups
#  
#  head(DE_genes_samples[[1]])
#  DE_genes_samples[[2]]
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  DE_genes_cl1_vs_3 <- GEX_DEgenes(GEX= vgm[[2]],min.pct = .25, grouping.column = "seurat_clusters",group1 = "1", group2 = "3",return.plot = "heatmap",up.genes = 10,down.genes = 10,logFC = F)
#  
#  #head(DE_genes_cl1_vs_3[[1]])
#  DE_genes_cl1_vs_3[[2]]
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #DE_clusters_all <- GEX_pairwise_DEGs(GEX = vgm[[2]], group.by = "seurat_clusters", min.pct = 0.25, RP.MT.filter = T, label.n.top.genes = 10, genes.to.label = c("CD74", "EBF1"), save.plot = F)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  dottile <- GEX_dottile_plot(GEX = vgm[[2]], genes = c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34","CD1D1","CR2","MS4A1","CXCR5","SELL","CD40","CD83","H2-AB1","H2-EB1","CD27","POU2AF1","NT5E","FAS","PDCD1LG2","PRDM1","ITGAM","IL10","IL12A","HAVCR2"), group.by = "seurat_clusters", threshold.to.plot = 5)
#  #threshold.to.plot specifies how many % of cells have to express a gene to show a dot.
#  
#  dottile + ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank(), legend.position = "bottom")
#  #For visualisation purposes in the RMD format
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #overview
#  coexpression_dotmap <- GEX_coexpression_coefficient(GEX = vgm[[2]], genes = c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34"), plot.dotmap = T)
#  
#  coexpression_dotmap + ggplot2::theme(legend.position = "bottom")
#  
#  #detail of two selected genes
#  GEX_scatter_coexpression(GEX = vgm[[2]], gene.1 = "CD19", gene.2 = "EBF1", color.theme = "darkred")
#  
#  #to save use:
#  #ggsave(last_plot(), filename = "Coexpression_scatter.png")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  print(colnames(vgm[[1]]))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  print(unique(vgm[[1]]$VDJ_cgene))
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  cat(" Cell count by number of VDJ chains")
#  print(table(vgm[[1]]$Nr_of_VDJ_chains))
#  cat("\n Cell count by number of VJ chains")
#  print(table(vgm[[1]]$Nr_of_VJ_chains))
#  
#  #Subset the VGM matrix to only include cells with 1 VDJ and 1VJ chain
#  #vgm[[1]] <- subset(vgm[[1]], Nr_of_VJ_chains == 1 & Nr_of_VJ_chains == 1)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  vgm[[1]] <- VDJ_clonotype(VDJ = vgm[[1]], clone.strategy = "cdr3.aa", global.clonotype = F, VDJ.VJ.1chain = F, hierarchical = "single.chains") #Not filtering cells with counts other than 1VDJ 1VJ chain and integrating these cells hierarchically into clonotypes
#  
#  cat(" Nr and distribution of clonotypes using exact CDR3.aa matching \n")
#  print(length(unique(vgm[[1]]$clonotype_id_cdr3.aa)))
#  print(table(vgm[[1]]$clonotype_frequency_cdr3.aa)) #Check distribution of clonotypes with identical CDR3 aa sequences
#  
#  vgm[[1]] <- VDJ_clonotype(VDJ = vgm[[1]], clone.strategy = "hvj.lvj", global.clonotype = F, output.format = "vgm", VDJ.VJ.1chain = F, hierarchical = "single.chains")
#  
#  cat("\n Nr and distribution of clonotypes using germline gene matching \n ")
#  print(length(unique(vgm[[1]]$clonotype_id_hvj.lvj)))
#  print(table(vgm[[1]]$clonotype_frequency_hvj.lvj)) #Check distribution of clonotypes with identical germline genes
#  

## ---- fig.show='hold', message=FALSE,results = 'hide', eval = FALSE-----------
#  
#  vgm <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#                                 GEX.out.directory.list = GEX.out.directory.list,
#                                 GEX.integrate = T,
#                                 VDJ.combine = T,
#                                 integrate.GEX.to.VDJ = T,
#                                 integrate.VDJ.to.GEX = T,
#                                 exclude.GEX.not.in.VDJ = F,
#                                 filter.overlapping.barcodes.GEX = T,
#                                 filter.overlapping.barcodes.VDJ = T,
#                                 exclude.on.cell.state.markers = c("CD3E"),
#                                 get.VDJ.stats = T,
#                                 parallel.processing = "none",
#                                 trim.and.align = T, #Set this to TRUE for full length sequence recovery
#                                 group.id = c(1,2),
#                                 gap.opening.cost = 10,
#                                 gap.extension.cost = 4) #Tweak to optimize alignments if neccessary
#  
#  #saving this for later
#  #saveRDS(vgm, "VDJ_GEX_matrix_agedCNS.rds")

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  print(vgm[[1]][1,])
#  

## ---- eval = FALSE, fig.show='hold',message=FALSE-----------------------------
#  
#  ### WARNING: You will need to download MiXCR and change the mixcr.directory to the location of MiXCR
#  #VDJ_mixcr_out <- VDJ_call_MIXCR(VDJ.matrix = vgm[[1]], mixcr.directory = "~/Downloads/mixcr.jar",species = "mmu", platypus.version = "v3", operating.system = "Darwin", simplify = T)
#  #set simplify to T to append only a selected columns of the MIXCR output to the vgm matrix
#  

## ---- eval = FALSE,message=FALSE----------------------------------------------
#  #check working directory
#  #getwd()
#  
#  ### WARNING: You will need to download MiXCR and have it in your current working directory
#  VDJ_mixcr_out <- VDJ_call_MIXCR(VDJ = vgm[[1]], mixcr.directory = "Is set automatically to current working directory",species = "mmu", platypus.version = "v3", operating.system = "Windows", simplify = T)
#  #set simplify to T to append only a selected columns of the MIXCR output to the vgm matrix
#  #set to False to save as separate object with the complete MIXCR output as in this case
#  

## ---- eval = FALSE, fig.show='hold',message=TRUE------------------------------
#  
#  VDJ_plot_SHM(VDJ = VDJ_mixcr_out, group.by = "sample_id", quantile.label = 0.95)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  donuts <- VDJ_clonal_donut(VDJ = vgm[[1]], expanded.colors = c("grey50", "grey65", "grey80"), non.expanded.color = "black", counts.to.use = "freq_column")
#  #Counts to use = "freq_column" uses the counts in the clonotype_frequency column. Counts to use = "vgm" simply counts the rows of a given clonotype in the VGM table (these counts may differ if cells have been filtered out due to overlapping barcodes or if another clonotyping strategy was used)
#  
#  donuts[[1]]

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #Shannon Evenness for the VDJ chain CDR3
#  diversity_plot <- VDJ_diversity(VDJ = vgm[[1]],feature.columns = c("VDJ_cdr3s_aa"),grouping.column = "sample_id",metric = c("shannonevenness"), platypus.version = "v3", subsample.to.same.n = T)
#  diversity_plot
#  
#  #Gini-Simpson index for pasted VDJ and VJ chain CDR3s
#  diversity_plot <- VDJ_diversity(VDJ = vgm[[1]],feature.columns = c("VDJ_cdr3s_aa", "VJ_cdr3s_aa"),grouping.column = "sample_id",metric = c("ginisimpson"), platypus.version = "v3", subsample.to.same.n = T)
#  diversity_plot
#  
#  #exact values can be retrived as
#  print(head(diversity_plot$data))
#  
#  #Jaccard index between repertoires of the two samples
#  diversity_plot <- VDJ_diversity(VDJ = vgm[[1]],feature.columns = c("VDJ_cgene"),grouping.column = "sample_id",metric = c("jaccard"), platypus.version = "v3", subsample.to.same.n = T)
#  diversity_plot
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  clonal_out <- VDJ_clonal_expansion(VDJ = vgm[[1]],celltype = "Bcells",clones = "30", platypus.version = "v3", group.by = "sample_id", color.by = "isotype", isotypes.to.plot = "all", treat.incomplete.clones = "exclude", treat.incomplete.cells = "proportional")
#  #group by specifies how many separate plots should be generated. If vgm contains global clonotype information this can be set to "none"
#  print(clonal_out[[1]])
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #subsampled_VGM <- dplyr::sample_n(vgm[[1]], 60)
#  
#  agedCNS_igraph <- VDJ_network(VDJ = vgm[[1]],per.sample = T,distance.cutoff = 8, platypus.version = "v3")
#  
#  for(i in 1:2){
#  igraph::plot.igraph(agedCNS_igraph[[1]][[i]],vertex.label=NA,vertex.size=7)
#  }
#  
#  #For a plot including clonotype frequencies
#  agedCNS_igraph <- VDJ_network(VDJ = vgm[[1]],per.sample = F,distance.cutoff = 8, platypus.version = "v3")
#  
#  igraph::plot.igraph(agedCNS_igraph[[4]],vertex.label=NA,vertex.size=3+(.03*as.numeric(agedCNS_igraph[[2]]$clonotype_frequency)),vertex.color= as.factor(vgm[[1]]$sample_id))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #First calculate adjacency matrix for V gene usage
#  agedCNS_Vgene_usage <- VDJ_Vgene_usage(VDJ = vgm[[1]], platypus.version = "v3")
#  
#  #library(pheatmap)
#  pheatmap::pheatmap(agedCNS_Vgene_usage[[1]],show_rownames = F,show_colnames = F)
#  
#  print(class(agedCNS_Vgene_usage[[1]]))
#  print(head(rownames(agedCNS_Vgene_usage[[1]])))
#  print(head(colnames(agedCNS_Vgene_usage[[1]])))
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  agedCNS_Vgene_usage_barplot <- VDJ_Vgene_usage_barplot(VDJ = vgm[[1]], HC.gene.number = 10, LC.Vgene = T, LC.gene.number = 10, platypus.version = "v3")
#  agedCNS_Vgene_usage_barplot[[1]]
#  
#  #VDJ chains
#  agedCNS_Vgene_usage_stackedbarplot <- VDJ_Vgene_usage_stacked_barplot(VDJ = vgm[[1]], LC.Vgene = F,HC.gene.number = 10, Fraction.HC = 1, platypus.version = "v3")
#  agedCNS_Vgene_usage_stackedbarplot
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #VDJ and VJ V and J gene pairing
#  vj_circos_bcells <- VDJ_VJ_usage_circos(VDJ = vgm[[1]], c.threshold = 1,label.threshold=2,cell.level = T, A.or.B = "both", platypus.version = "v3")
#  
#  #VDJ and VJ pairing
#  HL_circos_bcells <- VDJ_alpha_beta_Vgene_circos(VDJ = vgm[[1]], c.threshold = 1,label.threshold=2,cell.level = T, V.or.J= "both", platypus.version = "v3")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  pasted_CDR3s <- paste0(vgm[[1]]$VDJ_cdr3s_aa, vgm[[1]]$VJ_cdr3s_aa)
#  
#  agedCNS_CDR3_logoplot <- VDJ_logoplot_vector(cdr3.vector = pasted_CDR3s, seq_type = "aa", length_cdr3 = "auto")
#  
#  print(agedCNS_CDR3_logoplot)
#  

## ---- eval = FALSE------------------------------------------------------------
#  
#  variants_agedCNS <- VDJ_variants_per_clone(VDJ = vgm[[1]], variants.of = c("VDJ_cdr3s_aa", "VJ_cdr3s_aa"), clonotypes.col = "clonotype_id_10x", split.by = "sample_id", stringDist.method = "levenshtein")
#  
#  head(variants_agedCNS[[1]])
#  
#  #set split.by to "none" if clonotyping was conducted across all samples
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #overlap of VDJ V genes
#  VDJv_overlap <- VDJ_overlap_heatmap(VDJ = vgm[[1]],feature.columns = c("VDJ_vgene") ,grouping.column = "sample_id", axis.label.size = 20, pvalues.label.size = 12, platypus.version = "v3", add.barcode.table = T, plot.type = "ggplot")
#  
#  VDJv_overlap[[2]] #summary of overlap
#  VDJv_overlap[[2]] #Table of overlapping items
#  
#  #overlap of clones
#  VDJv_overlap <- VDJ_overlap_heatmap(VDJ = vgm[[1]],feature.columns = c("VDJ_cdr3s_aa","VJ_cdr3s_aa") ,grouping.column = "sample_id", axis.label.size = 20, pvalues.label.size = 12, platypus.version = "v3", add.barcode.table = T, plot.type = "ggplot")
#  #Pheatmap plots will function only with length(unique(grouping.column)) > 3
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #Nr of cells for which VDJ info is available
#  nrow(vgm[[1]])
#  
#  #Nr of cells for which GEX info is available
#  ncol(vgm[[2]])
#  
#  #VDJ sequences for which GEX is available
#  sum(vgm[[2]]$VDJ_available)
#  
#  #We can also plot this
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", group.by = "VDJ_available", shuffle = T)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  clonal_out <- VDJ_clonal_expansion(VDJ = vgm[[1]],celltype = "Bcells",clones = "30", platypus.version = "v3", group.by = "sample_id", color.by = "seurat_clusters")
#  #group by specifies how many separate plots should be generated. If vgm contains global clonotype information this can be set to "none"
#  print(clonal_out[[1]][[2]])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  clonal_out <- GEX_phenotype_per_clone(GEX = vgm[[2]], GEX.clonotypes = "topclones", GEX.group.by = "seurat_clusters", platypus.version = "v3")
#  #If vgm contains global clonotype information this can be set global.clonotypes to TRUE
#  clonal_out[[1]]

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #here we overlay the top 5 clones
#  plot_out <- VDJ_GEX_overlay_clones(GEX = vgm[[2]], reduction = "umap", n.clones = 5, by.sample = F, ncol.facet = 1, split.plot.and.legend = F, pt.size = 0.5)
#  
#  plot_out[[1]] # the plot
#  
#  #We can also plot any clone of interest by adding a column were TRUE values select which clones are to be plotted
#  interesting_clones <- c("clonotype7", "clonotype42")
#  vgm[[2]]@meta.data$clones_to_plot <- FALSE
#  vgm[[2]]@meta.data$clones_to_plot[which(vgm[[2]]$clonotype_id_10x %in% interesting_clones)] <- TRUE
#  
#  plot_out <- VDJ_GEX_overlay_clones(GEX = vgm[[2]], reduction = "umap", clones.to.plot = "clones_to_plot", by.sample = F, ncol.facet = 1, split.plot.and.legend = T, pt.size = 0.5)
#  
#  plot_out[[1]] # the plot
#  plot(plot_out[[2]]) # the legend. Alternatively use gridExtra::grid.arrange(plot_out[[2]])
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  vgm[[2]]$Expanded <- FALSE
#  vgm[[2]]$Expanded[which(vgm[[2]]$clonotype_frequency > 1)] <- TRUE
#  
#  table(vgm[[2]]$Expanded)
#  
#  #We can use the dottile function to look at selected genes
#  GEX_dottile_plot(GEX = vgm[[2]], genes = c("CD19", "CD74","IL2RA", "CD27","CD80"), group.by = "Expanded", threshold.to.plot = 1) + ggplot2::theme(legend.position = "bottom")

## ---- eval = FALSE------------------------------------------------------------
#  
#  DE_byexpansion <- GEX_DEgenes(GEX = vgm[[2]], min.pct = 0.25, group1 = "TRUE", group2 = "FALSE", grouping.column = "Expanded", return.plot = "volcano", label.n.top.genes = 10, platypus.version = "v3")
#  
#  head(DE_byexpansion[[1]])
#  DE_byexpansion[[2]]

## ---- fig.show='hold'---------------------------------------------------------
sessionInfo()

