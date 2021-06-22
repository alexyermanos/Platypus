## ----setup, include=FALSE-----------------------------------------------------
source_code_dir <- "C:/Users/vickr/Documents/GitHub/platypus_dev_VK" 
file_path_vec <- list.files(source_code_dir, full.names = T)
for(f_path in file_path_vec){
  source(f_path)}
gc()

knitr::opts_chunk$set(fig.width=7, fig.height=7) 

## ---- fig.show='hold', message=FALSE------------------------------------------

### Removing any previous versions of the package
#First can ensure that there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Dependencies 
#install.packages("stringr")


### Downloading and installing Platypus

# First we need to download the most recent version from the master branch at https://github.com/alexyermanos/Platypus we can install the package using the following command. 
# WARNING: This needs to be replaced with your own directory where the downloaded package is found

# For MacOS users it may look like this
#install.packages("~/Downloads/Platypus_3.1.tar.gz", repos = NULL, type="source")

# For windows it will likely look something like this. 
# WARNING: You will need to replace 'YourPCName' with your user name for the windows account in the directory. 
#install.packages("C:/Users/YourPCName/Downloads/Platypus_3.1.tar.gz", repos = NULL, type="source")

# Now we can load the installed package into the R environment. In case of problems with installing other R packages that are used in Platypus, please see the README file at the https://github.com/alexyermanos/Platypus, where we outline how to install the other R packages for both Windows and MacOS.
#library(Platypus)

# The individual R functions can additionally be found on the github in the Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment in case not all functions are desired. These functions are actively updated and may include more features than the in original tar.gz file. 



## ----eval = FALSE, fig.show='hold', message=FALSE,results = 'hide'------------
#  
#  ### Downloading the test data for VDJ_GEX_matrix
#  # While the Platypus manuscript uses the COVID-19 data, the vignette for Platypus v3 will use the data from B cells in the aged CNS, which can be found at the following link https://polybox.ethz.ch/index.php/s/fxQJ3NrRSwiPSSo This small dataset contains VDJ (separate libraries for B) and GEX libraries from the central nervous system of two murine samples. More information can be found https://doi.org/10.1098/rspb.2020.2793
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
#  vgm <- Platypus::VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#                                 GEX.out.directory.list = GEX.out.directory.list,
#                                 GEX.integrate = T,
#                                 VDJ.combine = T,
#                                 integrate.GEX.to.VDJ = T,
#                                 integrate.VDJ.to.GEX = T,
#                                 exclude.GEX.not.in.VDJ = F,
#                                 filter.overlapping.barcodes.GEX = F,
#                                 filter.overlapping.barcodes.VDJ = F,
#                                 exclude.on.cell.state.markers = c("This is not a gene"),
#                                 get.VDJ.stats = T,
#                                 parallel.processing = "none",
#                                 subsample.barcodes = F,
#                                 trim.and.align = F,
#                                 group.id = c(1,2))
#  
#  ## The output will be a list -
#  # vgm[[1]] corresponds to the VDJ master dataframe
#  # vgm[[2]] corresponds to the GEX in the form of a seurat object
#  # vgm[[3]] corresponds to the output of VDJ_stats subfunction - which provides information about the number of chains, sequencing reads, etc
#  # vgm[[4]] holds the input parameters
#  # vgm[[5]] holds the sessionInfo output at the time of function call
#  
#  
#  ## This function can similarly be used when only VDJ or GEX data is present. Simply do only provide the path list of either GEX or VDJ
#  # VDJ_comb_gex <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#  #                                #GEX.out.directory.list = GEX.out.directory.list,
#  #                                GEX.integrate = F,
#  #                                VDJ.combine = T,
#  #                                integrate.GEX.to.VDJ = F,
#  #                                integrate.VDJ.to.GEX = F,
#  #                                exclude.GEX.not.in.VDJ = F,
#  #                                filter.overlapping.barcodes.GEX = F,
#  #                                filter.overlapping.barcodes.VDJ = F,
#  #                                get.VDJ.stats = F,
#  #                                parallel.processing = "none",
#  #                                subsample.barcodes = F,
#  #                                trim.and.align = F,
#  #                                group.id = c(1,2))
#  

## ----eval = FALSE, fig.show='hold', message=FALSE-----------------------------
#  
#  head(colnames(vgm[[1]]))
#  
#  ## By setting integrate.GEX.to.VDJ and integrate.VDJ.to.GEX to T, VDJ and GEX information will be found in vgm[[1]] and vgm[[2]] objects.
#  # For example, the seurat-determined cluster is attached to each cell in the VDJ library by
#  head(vgm[[1]]$seurat_clusters)
#  # which corresponds to cells with the following VDJ_cdr3
#  head(vgm[[1]]$VDJ_cdr3s_aa)
#  # an NA indicates that the cell barcode in the VDJ library was not detected in the GEX object (or was filtered out, depending on mitochondrial gene limits, etc)
#  
#  ## Setting trim.and.align to TRUE will dramatically increase the run time but will also provide full-length VDJ and VJ sequences.
#  
#  ## This function can similarly be used when only VDJ or GEX data is present. Simply leave the list
#  # VDJ_comb_gex <- VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#  #                                #GEX.out.directory.list = GEX.out.directory.list,
#  #                                GEX.integrate = F,
#  #                                VDJ.combine = T,
#  #                                integrate.GEX.to.VDJ = F,
#  #                                integrate.VDJ.to.GEX = F,
#  #                                exclude.GEX.not.in.VDJ = F,
#  #                                filter.overlapping.barcodes.GEX = F,
#  #                                filter.overlapping.barcodes.VDJ = F,
#  #                                get.VDJ.stats = F,
#  #                                parallel.processing = "none",
#  #                                subsample.barcodes = F,
#  #                                trim.and.align = F,
#  #                                group.id = c(1,2))
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
#  GEX_proportions_barplot(GEX.matrix = vgm[[2]], stacked.plot = T)
#  #This function may also be used to plot proportions in other groups. For this use the source.group and target.group parameters.

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  Seurat::FeaturePlot(vgm[[2]],reduction = "umap",features = c("CD19","PTPRC", "EBF1", "H2-K1"))
#  
#  #To easily scout through genes in the dataset use:
#  View(as.data.frame(rownames(vgm[[2]])))
#  

## ---- results='hide', eval = FALSE--------------------------------------------
#  vgm[[2]] <- Platypus::GEX_phenotype(vgm[[2]], default = T)
#  
#  #vgm[[2]] <- Platypus::GEX_phenotype(vgm[[2]], default = F,cell.state.markers=c("CD8A+;CCL5+;CD44+;IL7R-;CD19-","CD8A+;CCL5-;CD44+;IL7R+;CD19-"),cell.state.names=c("EffectorCD8","MemoryCD8"))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  Seurat::DimPlot(vgm[[2]],reduction = "umap", group.by = "cell.state")

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  ## Warning: running this function will take a while
#  gene_expression_cluster <- Platypus::GEX_cluster_genes(vgm[[2]],min.pct = 0.25)
#  
#  length(gene_expression_cluster) ## length of 8, corresponding to 8 clusters
#  length(unique(vgm[[2]]$seurat_clusters)) ## length of 8

## ---- eval = FALSE------------------------------------------------------------
#  head(gene_expression_cluster[[1]])

## ---- eval = FALSE------------------------------------------------------------
#  print(sapply(gene_expression_cluster,nrow))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  agedCNS_heatmap_clusters <- GEX_cluster_genes_heatmap(GEX = vgm[[2]], GEX_cluster_genes.output = gene_expression_cluster,n.genes.per.cluster = 3,max.cell = 30,metric = "top_logFC", platypus.version = "v3")
#  
#  print(agedCNS_heatmap_clusters)

## ----results='hide', eval = FALSE---------------------------------------------
#  ontology_agedCNS <- Platypus::GEX_GOterm(GEX.cluster.genes.output = gene_expression_cluster[1], topNgenes = 10, go.plots = T)
#  head(ontology_agedCNS[[1]])

## ---- eval = FALSE------------------------------------------------------------
#  top_10_genes_per_cluster <- GEX_topN_DE_genes_per_cluster(GEX_cluster_genes.output = gene_expression_cluster, n.genes = 10, by_FC = T)
#  head(top_10_genes_per_cluster)

## ----GSEA, eval = FALSE-------------------------------------------------------
#  gsea_EAE <- GEX_GSEA(GEX.cluster.genes.output = gene_expression_cluster[[1]], MT.Rb.filter = T, path_to_pathways = "~/Downloads/c7.all.v7.4.symbols.gmt")

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  DE_genes_samples <- GEX_DEgenes(GEX = vgm[[2]],min.pct = .25, grouping.column = "sample_id",group1 = "s1", group2 = "s2",return.plot = TRUE,up.genes = 10,down.genes = 10,logFC = F)
#  
#  head(DE_genes_samples[[1]])
#  DE_genes_samples[[2]]

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  DE_genes_cl1_vs_3 <- GEX_DEgenes(GEX= vgm[[2]],min.pct = .25, grouping.column = "seurat_clusters",group1 = "1", group2 = "3",return.plot = TRUE,up.genes = 10,down.genes = 10,logFC = F)
#  
#  head(DE_genes_cl1_vs_3[[1]])
#  DE_genes_cl1_vs_3[[2]]

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  #necessary for plotting
#  library(ggrepel)
#  
#  DE_clusters_all <- GEX_pairwise_DEGs(GEX.matrix = vgm[[2]], group.by = "seurat_clusters", min.pct = 0.25, RP.MT.filter = T, label.n.top.genes = 30, genes.to.label = c("CD74", "EBF1"), save.plot = T)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  GEX_dottile_plot(GEX.matrix = vgm[[2]], genes = c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34","CD1D1","CR2","MS4A1","CXCR5","SELL","CD40","CD83","H2-AB1","H2-EB1","CD27","POU2AF1","NT5E","FAS","PDCD1LG2","PRDM1","ITGAM","IL10","IL12A","HAVCR2"), group.by = "seurat_clusters", threshold.to.plot = 5)
#  #threshold.to.plot specifies how many % of cells have to express a gene to show a dot.
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #overview
#  coex_plot <- GEX_coexpression_coefficient(GEX.matrix = vgm[[2]], genes = c("CD19", "CD74","SDC1", "EBF1","PTPRC","CD93","CD38","CD24A","CD34"))
#  coex_plot[[1]]
#  
#  #detail of two selected genes
#  GEX_scatter_coexpression(GEX.matrix = vgm[[2]], gene.1 = "CD19", gene.2 = "EBF1", color.theme = "darkred")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  print(colnames(vgm[[1]]))

## ---- eval = FALSE------------------------------------------------------------
#  print(unique(vgm[[1]]$VDJ_cgene))

## ---- eval = FALSE------------------------------------------------------------
#  print(table(vgm[[1]]$Nr_of_VDJ_chains))
#  print(table(vgm[[1]]$Nr_of_VJ_chains))
#  
#  #vgm[[1]] <- subset(vgm[[1]], Nr_of_VJ_chains == 1 & Nr_of_VJ_chains == 1)

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  vgm <- Platypus::VDJ_clonotype(platypus.version = "v3", VDJ.GEX.matrix = vgm, clone.strategy = "cdr3.aa", global.clonotype = F, output.format = "vgm")
#  
#  print(length(unique(vgm[[1]]$new_clonal_feature)))
#  print(table(vgm[[1]]$new_clonal_frequency)) #Check distribution of clonotypes with identical CDR3 aa sequences
#  
#  vgm <- Platypus::VDJ_clonotype(platypus.version = "v3", VDJ.GEX.matrix = vgm, clone.strategy = "hvj.lvj", global.clonotype = F, output.format = "vgm")
#  
#  print(length(unique(vgm[[1]]$new_clonal_feature)))
#  print(table(vgm[[1]]$new_clonal_frequency)) #Check distribution of clonotypes with identical germline genes

## ---- fig.show='hold', message=FALSE,results = 'hide', eval = FALSE-----------
#  
#  vgm <- Platypus::VDJ_GEX_matrix(VDJ.out.directory.list = VDJ.out.directory.list,
#                                 GEX.out.directory.list = GEX.out.directory.list,
#                                 GEX.integrate = T,
#                                 VDJ.combine = T,
#                                 integrate.GEX.to.VDJ = T,
#                                 integrate.VDJ.to.GEX = T,
#                                 exclude.GEX.not.in.VDJ = F,
#                                 filter.overlapping.barcodes.GEX = F,
#                                 filter.overlapping.barcodes.VDJ = F,
#                                 exclude.on.cell.state.markers = c("This is not a gene"),
#                                 get.VDJ.stats = T,
#                                 parallel.processing = "none",
#                                 subsample.barcodes = F,
#                                 group.id = c(1,2),
#                                 trim.and.align = T,
#                                 gap.opening.cost = 10,
#                                 gap.extension.cost = 4) #Tweak to optimize alignments if neccessary
#  
#  #saving this for later
#  #saveRDS(vgm, "VDJ_GEX_matrix_agedCNS.rds")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  print(vgm[[1]][1,])
#  

## ---- eval = FALSE, fig.show='hold',message=FALSE-----------------------------
#  
#  ### WARNING: You will need to download MiXCR and change the mixcr.directory to the location of MiXCR
#  #vgm[[1]] <- VDJ_call_MIXCR(VDJ.matrix = vgm[[1]], mixcr.directory = "~/Downloads/mixcr.jar",species = "mmu", platypus.version = "v3", operating.system = "Darwin", simplify = T)
#  #set simplify to T to append only a selected columns of the MIXCR output to the vgm matrix
#  

## ---- eval = FALSE, fig.show='hold',message=FALSE-----------------------------
#  #check working directory
#  #getwd()
#  
#  ### WARNING: You will need to download MiXCR and change the mixcr.directory to the location of MiXCR
#  vgm[[1]] <- VDJ_call_MIXCR(VDJ.matrix = vgm[[1]], mixcr.directory = "Is set automatically",species = "mmu", platypus.version = "v3", operating.system = "Windows", simplify = T)
#  #set simplify to T to append only a selected columns of the MIXCR output to the vgm matrix
#  

## ---- eval = FALSE, fig.show='hold',message=FALSE-----------------------------
#  
#  library(ggrepel)
#  
#  VDJ_plot_SHM(VDJ.matrix = vgm[[1]], group.by = "sample_id", quantile.label = 0.9)
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #FOR MAC / UNIX users
#  #extracted_agedCNS_germline <- Platypus::VDJ_extract_germline(VDJ.per.clone=covid_single_cell,mixcr.directory="~/Downloads/mixcr-3.0.12/mixcr",extract.VDJRegion=T,species = "hsa")
#  
#  #print(colnames(extracted_covid_germline[[1]][[1]])) ## column names of the germlines from the repertoire corresponding to sample 1
#  #print(nrow(extracted_covid_germline[[1]][[1]])) ## germline sequences corresponding to 2298 clones in the first patient
#  #print(nrow(extracted_covid_germline[[2]][[1]])) ## germline sequences corresponding to 2298 clones in the second patient.
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #covid_single_clone_with_JSON <- VDJ_per_cell(clonotype.list = covid_vdj_repertoire_bcells,VDJ.out.directory =VDJ.out.directory.list,JSON = T)
#  
#  
#  # covid_single_cell was the object from VDJ_per_clone
#  
#  #print(covid_single_clone_with_JSON[[1]][[1]]$trimmed_HC_sequence[1])  # trimmmed sequence
#  #print(covid_single_clone_with_JSON[[1]][[1]]$trimmed_HC_germline[1])  # trimmmed germline sequence
#  
#  #stringdist::stringdist(covid_single_clone_with_JSON[[1]][[1]]$sequence_HC[1], VDJ.per.cell[[1]][[1]]$trimmed_ref_HC[1])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  
#  #covid_clonal_lineages <- Platypus::VDJ_clonal_lineages(call_MIXCR.output=covid_vdj_region, VDJ_extract_germline.output=extracted_covid_germline,as.nucleotide=F,with.germline=T)
#  
#  #print(colnames(covid_clonal_lineages[[1]][[1]])) ## dataframe with the columns Seq and Name.
#  
#  #print(covid_clonal_lineages[[1]][[1]]$Seq[1])
#  
#  #print(covid_clonal_lineages[[1]][[1]]$Name[1]) #"clonotype3_1_IGHA1_AACCATGAGTGGAGTC-1"
#  ## Here we can see that the above sequence corresponds to the original clonotype3 but is now the first clonal lineage (as seen by the _1_ before the isotype). Furthermore, the isotype of this cell was of the IGHA1. Lastly, we have the barcode of the cell at the end, allowing us to look back at other cell-specific properties if wanted.
#  
#  
#  #print(covid_clonal_lineages[[1]][[3]]$Name[3]) #"clonotype7_3_IGHA1_AGTGAGGTCGAGAACG-1"
#  ## Again, this was originally the 7th clonotype, and is now the third most clonally expanded lineage ("_3_") of the IGHA1 isotype with given barcode.
#  
#  #print(tail(covid_clonal_lineages[[1]][[3]]$Name))  ## here at the end of the dataframe the user can also see the germline sequence, which has the "Name" of "germline".
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  #covid_trees <- Platypus::VDJ_tree(clonal.lineages = covid_clonal_lineages,with.germline=T,min.sequences = 5,max.sequences = 30,unique.sequences = T)
#  #plot(covid_trees[[1]][[1]])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  donuts <- VDJ_clonal_donut(VDJ.matrix = vgm[[1]], expanded.colors = c("grey50", "grey65", "grey80"), non.expanded.color = "black")
#  
#  donuts

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  library(vegan)
#  
#  #Shannon Evenness for the VDJ chain CDR3
#  diversity_plot <- VDJ_diversity(VDJ.matrix = vgm[[1]],feature.columns = c("VDJ_cdr3s_aa"),grouping.column = "sample_id",metric = c("shannonevenness"), platypus.version = "v3", subsample.to.same.n = T)
#  diversity_plot
#  
#  #Gini-Simpson index for pasted VDJ and VJ chain CDR3s
#  diversity_plot <- VDJ_diversity(VDJ.matrix = vgm[[1]],feature.columns = c("VDJ_cdr3s_aa", "VJ_chain"),grouping.column = "sample_id",metric = c("ginisimpson"), platypus.version = "v3", subsample.to.same.n = T)
#  diversity_plot
#  
#  #exact values can be retrived as
#  print(head(diversity_plot$data))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  clonal_out <- VDJ_clonal_expansion(VDJ.matrix = vgm[[1]],celltype = "Bcells",clones = "30", platypus.version = "v3", group.by = "sample_id", color.by = "isotype")
#  #group by specifies how many separate plots should be generated. If vgm contains global clonotype information this can be set to "none"
#  print(clonal_out[[1]])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  ## Take the output from VDJ_analyze (or subsample as in this case the top 60 clones)
#  
#  agedCNS_igraph <- VDJ_network(VDJ.matrix = vgm[[1]],per.sample = T,distance.cutoff = 8,connected = F, platypus.version = "v3")
#  
#  for(i in 1:2){
#  igraph::plot.igraph(agedCNS_igraph[[i]],vertex.label=NA,vertex.size=7)
#  }
#  
#  #For a plot including clonotype frequencies
#  agedCNS_igraph <- VDJ_network(VDJ.matrix = vgm[[1]],per.sample = F,distance.cutoff = 8,connected = F, platypus.version = "v3")
#  
#  igraph::plot.igraph(agedCNS_igraph[[4]],vertex.label=NA,vertex.size=3+(.03*as.numeric(agedCNS_igraph[[2]]$clonotype_frequency)),vertex.color= as.factor(vgm[[1]]$sample_id))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #First calculate adjacency matrix for V gene usage
#  agedCNS_Vgene_usage <- VDJ_Vgene_usage(VDJ.matrix = vgm[[1]], platypus.version = "v3")
#  library(pheatmap)
#  
#  pheatmap::pheatmap(agedCNS_Vgene_usage[[1]],show_rownames = F,show_colnames = F)
#  
#  print(class(agedCNS_Vgene_usage[[1]]))
#  print(head(rownames(agedCNS_Vgene_usage[[1]])))
#  print(head(colnames(agedCNS_Vgene_usage[[1]])))

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  agedCNS_Vgene_usage_barplot <- VDJ_Vgene_usage_barplot(VDJ.matrix = vgm[[1]], HC.gene.number = 10, LC.Vgene = T, LC.gene.number = 10, platypus.version = "v3")
#  agedCNS_Vgene_usage_barplot[[1]]
#  
#  #VDJ chains
#  agedCNS_Vgene_usage_stackedbarplot <- VDJ_Vgene_usage_stacked_barplot(VDJ.matrix = vgm[[1]], LC.Vgene = F,HC.gene.number = 10, Fraction.HC = 1, platypus.version = "v3")
#  agedCNS_Vgene_usage_stackedbarplot
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  #VDJ and VJ V and J gene pairing
#  vj_circos_bcells <- VDJ_VJ_usage_circos(VDJ.GEX.matrix = vgm[[1]], c.threshold = 1,label.threshold=2,cell.level = T, A.or.B = "both", platypus.version = "v3")
#  
#  #VDJ and VJ pairing
#  HL_circos_bcells <- VDJ_alpha_beta_Vgene_circos(VDJ.GEX.matrix = vgm[[1]], c.threshold = 1,label.threshold=2,cell.level = T, V.or.J= "both", platypus.version = "v3")
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  pasted_CDR3s <- paste0(vgm[[1]]$VDJ_cdr3s_aa, vgm[[1]]$VJ_cdr3s_aa)
#  
#  agedCNS_CDR3_logoplot <- VDJ_logoplot_vector(cdr3.vector = pasted_CDR3s, seq_type = "aa", length_cdr3 = "auto")
#  
#  print(agedCNS_CDR3_logoplot)

## ---- eval = FALSE------------------------------------------------------------
#  variants_agedCNS <- VDJ_variants_per_clone(VDJ.matrix = vgm[[1]], variants.of = c("VDJ_cdr3s_aa", "VJ_cdr3s_aa"), clonotypes.col = "clonotype_id_10x", split.by = "sample_id", stringDist.method = "levenshtein")
#  
#  head(variants_agedCNS[[1]])
#  
#  #set split.by to "none" if clonotyping was conducted across all samples

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  #overlap of VDJ V genes
#  VDJv_overlap <- VDJ_overlap_heatmap(vgm[[1]],feature.columns = c("VDJ_vgene") ,grouping.column = "sample_id", axis.label.size = 20, pvalues.label.size = 12, platypus.version = "v3", add.barcode.table = T)
#  
#  VDJv_overlap[[2]] #summary of overlap
#  VDJv_overlap[[2]] #Table of overlapping items
#  
#  #overlap of clones
#  VDJv_overlap <- VDJ_overlap_heatmap(vgm[[1]],feature.columns = c("VDJ_cdr3s_aa","VJ_cdr3s_aa") ,grouping.column = "sample_id", axis.label.size = 20, pvalues.label.size = 12, platypus.version = "v3", add.barcode.table = T)
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
#  clonal_out <- VDJ_clonal_expansion(VDJ.matrix = vgm[[1]],celltype = "Bcells",clones = "30", platypus.version = "v3", group.by = "sample_id", color.by = "seurat_clusters")
#  #group by specifies how many separate plots should be generated. If vgm contains global clonotype information this can be set to "none"
#  print(clonal_out[[1]][[1]])

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  clonal_out <- GEX_phenotype_per_clone(GEX.matrix = vgm[[2]], GEX.clonotypes = "topclones", GEX.group.by = "seurat_clusters", platypus.version = "v3")
#  #If vgm contains global clonotype information this can be set global.clonotypes to TRUE
#  clonal_out[[1]]

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  library(cowplot)
#  library(gridExtra)
#  
#  #here we overlay the top 5 clones
#  plot_out <- VDJ_GEX_overlay_clones(GEX.output = vgm[[2]], reduction = "umap", n.clones = 5, by.sample = F, ncol.facet = 1, split.plot.and.legend = T, pt.size = 0.5)
#  
#  plot_out[[1]] # the plot
#  grid.arrange(plot_out[[2]]) # the legend
#  
#  #We can also plot any clone of interest by adding a column were TRUE values select which clones are to be plotted
#  interesting_clones <- c("clonotype7", "clonotype42")
#  vgm[[2]]@meta.data$clones_to_plot <- FALSE
#  vgm[[2]]@meta.data$clones_to_plot[which(vgm[[2]]$clonotype_id_10x %in% interesting_clones)] <- TRUE
#  
#  plot_out <- VDJ_GEX_overlay_clones(GEX.output = vgm[[2]], reduction = "umap", clones.to.plot = "clones_to_plot", by.sample = F, ncol.facet = 1, split.plot.and.legend = T, pt.size = 0.5)
#  
#  plot_out[[1]] # the plot
#  grid.arrange(plot_out[[2]]) # the legend
#  

## ---- fig.show='hold', eval = FALSE-------------------------------------------
#  
#  vgm[[2]]$Expanded <- FALSE
#  vgm[[2]]$Expanded[which(vgm[[2]]$clonotype_frequency > 1)] <- TRUE
#  
#  table(vgm[[2]]$Expanded)
#  
#  #We can use the dottile function to look at selected genes
#  GEX_dottile_plot(GEX.matrix = vgm[[2]], genes = c("CD19", "CD74","IL2RA", "CD27","CD80"), group.by = "Expanded", threshold.to.plot = 1)

## ---- eval = FALSE------------------------------------------------------------
#  
#  DE_byexpansion <- GEX_DEgenes(GEX = vgm[[2]], min.pct = 0.25, group1 = "TRUE", group2 = "FALSE", grouping.column = "Expanded", return.plot = T, platypus.version = "v3")
#  
#  head(DE_byexpansion[[1]])
#  DE_byexpansion[[2]]

## ---- fig.show='hold'---------------------------------------------------------
sessionInfo()

