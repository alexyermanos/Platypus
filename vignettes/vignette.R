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
dir_to_covid_gex[[1]] <- c("~/Downloads/covid_data/Patient1_GEX/",
                           "~/Downloads/covid_data/Patient2_GEX/")

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




