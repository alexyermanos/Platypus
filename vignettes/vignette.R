## ---- fig.show='hold'---------------------------------------------------------

### loading the package
#source("~/Downloads/Platypus-master/Platypus_1.1.tar.gz")
library(Platypus)

dir_to_covid_gex_individual <- list()
dir_to_covid_gex_individual[[1]] <- c("~/Downloads/covid_data/Patient1_GEX/")
dir_to_covid_gex_individual[[2]] <- c("~/Downloads/covid_data/Patient2_GEX/") ## adding the single directory 



# Now we run automate_GEX to perform the standard Seurat pipeline in a single line of code with the ability to change the parameters used by Seurat involving minimum read numbers, mitochondrial gene percentages, cluster resolution, etc. This example will result with a single Seurat GEX object. 


covid_gex_individual <- Platypus::automate_GEX(GEX.outs.directory.list = dir_to_covid_gex_individual,integration.method = "scale.data",mito.filter = 20,cluster.resolution = 0.5,VDJ.gene.filter = T)
#VDJ.gene.filter removes antibody and TCR genes from the data set (e.g. IGHV1-1) to not let clonality impact transcriptional clustering. 
length(covid_gex_individual) ## length of two 
class(covid_gex_individual) ## list class
class(covid_gex_individual[[1]]) ## Seurat object 
ncol(covid_gex_individual[[1]]) ## 6055 cells
ncol(covid_gex_individual[[2]]) ## 6871 cells

## Now can visualize the individual samples using UMAP/tsne plots 
Seurat::DimPlot(covid_gex_individual[[1]],reduction = "umap")



## ---- fig.show='hold'---------------------------------------------------------


Seurat::FeaturePlot(covid_gex_individual[[1]],reduction = "umap",features = c("CD4","CD8A","CD19"))
## Here we can see the B, CD4 and CD8 clusters within the individual patient.


### Integrating two GEX libraries into one Seurat object
dir_to_covid_gex <- list()
dir_to_covid_gex[[1]] <- c("~/PHD/covid_patients/GEX/Both_flowcell_alignments/5921925/filtered_feature_bc_matrix/",
                           "~/PHD/covid_patients/GEX/Both_flowcell_alignments/5922423/filtered_feature_bc_matrix/")

covid_gex <- Platypus::automate_GEX(GEX.outs.directory.list = dir_to_covid_gex[1:1],integration.method = "scale.data",mito.filter = 20,cluster.resolution = 0.5,VDJ.gene.filter = T)


### Now both gene expression samples will be integrated together. The covid_gex[[1]]$sample_id contains the sample origin information based on the order of the input directories. In this case, patient1 is the covid_gex[[1]]$sample_id==1 and patient2 is the covid_gex[[1]]$sample_id==2.

class(covid_gex) ## list
length(covid_gex) ### list of length one
class(covid_gex[[1]]) ## Seurat object

## We can extract which cells come from which sample based on the sample_id in the Seurat object. 
print(table(covid_gex[[1]]$sample_id)) 
## we can see that the first patient had 6871 cells and the second patient had 6055 cells found in their gene expression data based on the current filtering. 


## ---- fig.show='hold'---------------------------------------------------------


## Now can visualize the individual samples using UMAP/tsne plots 
Seurat::DimPlot(covid_gex[[1]],reduction = "umap")

### We observe that under this Seurat pipeline there are 12 distinct clusters. The number of clusters can be changed by altering the cluster.resolution argument in the automate_GEX function. In this current UMAP we do not see which cells are coming from which patient. 

## ---- fig.show='hold'---------------------------------------------------------

Seurat::DimPlot(covid_gex[[1]],reduction = "umap",split.by = "sample_id") # the resulting DimPlot now has separated the cells by the sample_id vector in the Seurat object. We can visually observe that the majority of clusters have cells from both patients, suggesting a similar distribution of transcriptional properties between the two samples. 


## ---- fig.show='hold'---------------------------------------------------------


Seurat::FeaturePlot(covid_gex[[1]],reduction = "umap",features = c("CD4","CD8A","CD19"))
## Again we can see the B, CD4 and CD8 clusters from the cells of both patients. 

## ---- fig.show='hold'---------------------------------------------------------

Seurat::FeaturePlot(covid_gex[[1]],reduction = "umap",features = c("CD4","CD8A","CD19"),split.by = "sample_id")
## Again, splitting by sample_id we can look at specific markers for each patient. We can see that both patients seem to have B cells (CD19+) , CD4 T cells, and CD8 T cells. 





## ---- fig.show='hold'---------------------------------------------------------
## GEX_cluster_genes will return a list containing the genes differentially up or down regulated for each cluster. 
## Warning: running this function will take a while
gene_expression_cluster <- Platypus::GEX_cluster_genes(covid_gex[[1]],min.pct = 0.25) 

length(gene_expression_cluster) ## length of 12, corresponding to 12 clusters
length(unique(covid_gex[[1]]$seurat_clusters)) ## length of 12 


### This gives the genes associated with cluster0 in the previously displayed umap - also corresponding to the cells found at covid_gex[[1]]$seurat_clusters==0 
head(gene_expression_cluster[[1]])
# We see the genes ANXA1, S100A4 etc are highly expressed in cluster0. 

# Here we can quantify the number of differentially expressed genes for all 12 clusters

print(sapply(gene_expression_cluster,nrow))



## ---- fig.show='hold'---------------------------------------------------------

covid_heatmap_clusters <- Platypus::GEX_cluster_genes_heatmap(automate_GEX.output = covid_gex[[1]],
                                                              GEX_cluster_genes.output = gene_expression_cluster,
                                                              n.genes.per.cluster = 5,max.cell = 50,
                                                              metric = "avg_logFC")

print(covid_heatmap_clusters)



## ---- fig.show='hold'---------------------------------------------------------

DE_genes_per_sample <- GEX_DEgenes_persample(automate.GEX=covid_gex[[1]],min.pct = .25,sample1 = "1",sample2 = "2")


head(DE_genes_per_sample) 
nrow(DE_genes_per_sample) ##71 rows 




## ---- fig.show='hold'---------------------------------------------------------
VDJ.out.directory.list <- list()
VDJ.out.directory.list[[1]] <- "~/Downloads/covid_data/Patient1_VDJ_Bcell_out/"
VDJ.out.directory.list[[2]] <- "~/Downloads/covid_data/Patient2_VDJ_Bcell_out/"

covid_vdj_repertoire_bcells <- Platypus::VDJ_analyze(VDJ.out.directory =VDJ.out.directory.list, filter.1HC.1LC = T) 


length(covid_vdj_repertoire_bcells) ## list of length two, first element corresponds to the clones from the first repertoire directory that was set as input. 
print(colnames(covid_vdj_repertoire_bcells[[1]]))

## Can see various information, including which barcodes make up the clonal family, the nt_clone_ids incase the clonotyping method is changed using the VDJ_clonotype function. Furthermore, the majority germline gene per clonal family is extracted from the contigs file. 


## ---- fig.show='hold'---------------------------------------------------------

covid_vdj_aminoacid_clonotype <- Platypus::VDJ_clonotype(clonotype.list=covid_vdj_repertoire_bcells,
                                                         clone.strategy="cdr3.aa")

length(covid_vdj_aminoacid_clonotype) # length of 2

print(colnames(covid_vdj_aminoacid_clonotype[[1]])) 


## ---- fig.show='hold'---------------------------------------------------------

covid_single_cell <- Platypus::VDJ_per_clone(clonotype.list = covid_vdj_repertoire_bcells,VDJ.out.directory =VDJ.out.directory.list)

print(paste("There are",length(covid_single_cell[[1]]),"unique nucleotide B cell clones in patient1"),sep="")

print(paste("There are",nrow(covid_single_cell[[1]][[1]]),"unique B cells in the most abundant clone in patient1"),sep="")

print(paste("There are",nrow(covid_single_cell[[2]][[1]]),"unique B cells in the most abundant clone in patient2"),sep="")


print(colnames(covid_single_cell[[1]][[1]]))



## ---- fig.show='hold'---------------------------------------------------------

covid_vdj_region <- Platypus::call_MIXCR(VDJ.per.clone = covid_single_cell,mixcr.directory = "~/Downloads/mixcr-3.0.12/mixcr",species = "hsa")

print(length(covid_vdj_region[[1]]))
print(colnames(covid_vdj_region[[1]][[1]]))

print(nchar(covid_vdj_region[[1]][[1]]$VDJ.AA.HC[2])) 

print(nchar(covid_vdj_region[[1]][[1]]$full_HC_germline[1]))


## ---- fig.show='hold'---------------------------------------------------------

extracted_covid_germline <- VDJ_extract_germline(VDJ.per.clone=covid_single_cell,mixcr.directory="~/Downloads/mixcr-3.0.12/mixcr",extract.VDJRegion=T,species = "hsa")

print(colnames(extracted_covid_germline[[1]][[1]])) ## column names of the germlines from the repertoire corresponding to patient1. 
print(nrow(extracted_covid_germline[[1]][[1]])) ## germline sequences corresponding to 2298 clones in the first patient 
print(nrow(extracted_covid_germline[[2]][[1]])) ## germline sequences corresponding to 2298 clones in the second patient. 


print(extracted_covid_germline[[1]][[1]]$aaSeqCDR3[1]) ## CAREL_FDYW - we can observe unproductive CDR3s for the germline sequence of the first clonotype 



print((extracted_covid_germline[[1]][[1]]$descrsR1[2]))  ## second row in the dataframe corresponds to clonotype4.
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.HC.LC[2])) ### pasted HC and LC together
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.HC[2])) ### just heavy chain for clonotype4 germline
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.LC[2])) ### just light chain for clonotype4 germline


## ---- fig.show='hold'---------------------------------------------------------


VDJ.per.cell <- VDJ_per_cell(clonotype.list = VDJ_analyze.output, VDJ.out.directory = VDJ.out.directory.list)

print(VDJ.per.cell[[1]][[1]]$sequence_HC)  # trimmmed sequence 
print(VDJ.per.cell[[1]][[1]]$trimmed_ref_HC)  # trimmmed germline sequence 
stringdist::stringdist(VDJ.per.cell[[1]][[1]]$sequence_HC, VDJ.per.cell[[1]][[1]]$trimmed_ref_HC) 


## ---- fig.show='hold'---------------------------------------------------------


covid_clonal_lineages <- VDJ_clonal_lineages(call_MIXCR.output=covid_vdj_region, VDJ_extract_germline.output=extracted_covid_germline,as.nucleotide=F,with.germline=T)

print(colnames(covid_clonal_lineages[[1]][[1]])) ## dataframe with the columns Seq and Name. 

print(covid_clonal_lineages[[1]][[1]]$Seq[1])

print(covid_clonal_lineages[[1]][[1]]$Name[1]) #"clonotype3_1_IGHA1_AACCATGAGTGGAGTC-1"
## Here we can see that the above sequence corresponds to the original clonotype3 but is now the first clonal lineage (as seen by the _1_ before the isotype). Furthermore, the isotype of this cell was of the IGHA1. Lastly we have the barcode of the cell at the end, allowing us to look back at other cell-specific properties if wanted. 


print(covid_clonal_lineages[[1]][[3]]$Name[3]) #"clonotype7_3_IGHA1_AGTGAGGTCGAGAACG-1"
## Again, this was originally the 7th clonotype, and is now the third most clonally expanded lineage ("_3_"). Again of the IGHA1 isotype with the following barcode. 

print(tail(covid_clonal_lineages[[1]][[3]]$Name))  ## here at the end of the dataframe the user can also see the germline sequence, which has the "Name" of "germline". 




## ---- fig.show='hold'---------------------------------------------------------


covid_trees <- VDJ_tree(clonal.lineages = covid_clonal_lineages,with.germline=T,min.sequences = 5,max.sequences = 30,unique.sequences = T)
plot(covid_trees[[1]][[1]])



## ---- fig.show='hold'---------------------------------------------------------

covid_isotypes <- Platypus::VDJ_isotypes_per_clone(VDJ_clonotype_output = covid_vdj_repertoire_bcells, VDJ_per_clone_output = covid_single_cell, clones = 30)

print(covid_isotypes[[1]])


## ---- fig.show='hold'---------------------------------------------------------

print(covid_isotypes[[2]])





## ---- fig.show='hold'---------------------------------------------------------

library(igraph)

## Take the output from VDJ_analyze (or subsample as in this case the top 60 clones)
network_clones_covid <- list()
network_clones_covid[[1]] <- covid_vdj_repertoire_bcells[[1]][1:60,]
network_clones_covid[[2]] <- covid_vdj_repertoire_bcells[[2]][1:60,]


covid_bcell_igraph <- Platypus::VDJ_network(network_clones_covid[1:2],per.mouse = F,distance.cutoff = 8,connected = F)

igraph::plot.igraph(covid_bcell_igraph[[4]],vertex.label=NA,vertex.size=7+(.06*covid_bcell_igraph[[2]]$frequency),vertex.color=factor(covid_bcell_igraph[[2]]$mouse))






## ---- fig.show='hold'---------------------------------------------------------

# 
# covid_vj_gene_usage <- Platypus::VDJ_Vgene_usage(VDJ.clonotype.output = covid_vdj_repertoire_bcells)
# 
# print(class(covid_vj_gene_usage[[1]]))
# print(head(rownames(covid_vj_gene_usage[[1]])))
# print(head(colnames(covid_vj_gene_usage[[1]])))



## ---- fig.show='hold'---------------------------------------------------------

covid_integrating_clonal_level <- Platypus::VDJ_GEX_integrate(GEX.object = covid_gex[[1]], 
                                                 clonotype.list =  covid_vdj_repertoire_bcells,
                                                 VDJ.per.clone = covid_single_cell,
                                                 clonotype.level = TRUE)

print(head(covid_integrating_clonal_level[[1]]$majority_cluster))
print(head(covid_integrating_clonal_level[[1]]$cluster_membership_percent))
print(head(covid_integrating_clonal_level[[1]]$cell_index))



## ---- fig.show='hold'---------------------------------------------------------

covid_integrating_cell_level <- VDJ_GEX_integrate(GEX.object = covid_gex[[1]], 
                                                 clonotype.list =  covid_vdj_repertoire_bcells[1:1],
                                                 VDJ.per.clone = covid_single_cell[1:1],
                                                 clonotype.level = FALSE)


print(head(covid_integrating_cell_level[[1]][[5]]$cluster_membership))
print(head(covid_integrating_cell_level[[1]][[5]]$cell_index))

## ---- fig.show='hold'---------------------------------------------------------

covid_clonotype_clusters_plot <-  Platypus::VDJ_GEX_expansion(GEX.list=covid_gex[[1]],
                                                              VDJ.GEX.integrate.list=covid_integrating_clonal_level,
                                                              highlight.isotype = "None",
                                                              highlight.number=1:20)


print(covid_clonotype_clusters_plot[[1]])



## ---- fig.show='hold'---------------------------------------------------------

covid_top10_umap <- visualize_clones_GEX(GEX.list=covid_gex,
                     VDJ.GEX.integrate.list=covid_integrating_clonal_level[1:1],
                     highlight.type="clonotype",
                     highlight.number=1:10,
                     reduction="umap")
print(covid_top10_umap[[1]])


## ---- fig.show='hold'---------------------------------------------------------

covid_gex_integrate <- clonotype_GEX(GEX.object=covid_gex[[1]],VDJ.per.clone = covid_single_cell)


covid_bcell_gene_heatmap <- GEX_heatmap(covid_gex_integrate,b.or.t = "b",clone.rank.threshold = 20,sample.index = 1)


## ---- fig.show='hold'---------------------------------------------------------

clonal_lineage_integrate <- VDJ_GEX_clonal_lineage_clusters(covid_integrating_cell_level,covid_clonal_lineages[1:1])

print(clonal_lineage_integrate[[1]][[1]]$Name[1])


