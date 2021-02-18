## ---- fig.show='hold', message=FALSE------------------------------------------

### Removing any previous versions of the package
#First, ensure there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Downloading and installing Platypus

# Download most recent version from master branch at https://github.com/alexyermanos/Platypus We can install the package using the following command. 
# WARNING: This needs to be replaced with your own directory where the downloaded package is found

# For MacOS users it may look like this:
install.packages("~/Downloads/Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# For windows it will likely look something like this: 
# WARNING: You will need to replace 'YourPCName' with your user name for the windows account in the directory. 
# install.packages("C:\Users\YourPCName\Downloads\Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# Now we can load the installed package into the R environment. In case of problems with installing other R packages that are used in Platypus, please see the README file at the https://github.com/alexyermanos/Platypus, where we outline how to install the other R packages for both Windows and MacOS.
library(Platypus)

# Individual R functions can additionally be found on the github Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment in case not all functions are desired. These functions are actively updated and may include more features than the in original tar.gz file. 


### Downloading the test data
# The COVID-19 data (~136 MB size of the zip file) can be found at the following link https://polybox.ethz.ch/index.php/s/fxQJ3NrRSwiPSSo This dataset contains VDJ (separate libraries for B and T cells) and GEX libraries from two convalescent COVID-19 patients.

# After downloading the zip file named "PlatypusTestData.zip", please unzip the file and find the path to the newly formed folder. Typically this will be in the Downloads folder, so the code below should work on MacOS. For Windows please uncomment the code and change the user name to match your PC.

directory_to_covid_patients_gex <- list()
directory_to_covid_patients_gex[[1]] <- c("~/Downloads/PlatypusTestData/Patient1_GEX/")
directory_to_covid_patients_gex[[2]] <- c("~/Downloads/PlatypusTestData/Patient2_GEX/")

# For Windows: 
#directory_to_covid_patients_gex[[1]] <- c("C:\Users\YourPCName\Downloads\PlatypusTestData\Patient1_GEX")
#directory_to_covid_patients_gex[[2]] <- c("C:\Users\YourPCName\Downloads\PlatypusTestData\Patient2_GEX")



## ---- fig.show='hold', message=FALSE------------------------------------------

covid_gex_patients_not_integrated <- Platypus::GEX_automate(GEX.outs.directory.list = directory_to_covid_patients_gex,integration.method = "scale.data",mito.filter = 20,cluster.resolution = 0.5,VDJ.gene.filter = T)


## ---- fig.show='hold'---------------------------------------------------------
length(covid_gex_patients_not_integrated) ## length of two 
class(covid_gex_patients_not_integrated) ## output is a list
class(covid_gex_patients_not_integrated[[1]]) ## Seurat object 
ncol(covid_gex_patients_not_integrated[[1]]) ## 6871 cells
ncol(covid_gex_patients_not_integrated[[2]]) ## 6055 cells


## ---- fig.show='hold'---------------------------------------------------------
Seurat::DimPlot(covid_gex_patients_not_integrated[[1]],reduction = "umap")

Seurat::DimPlot(covid_gex_patients_not_integrated[[1]],reduction = "pca")

Seurat::DimPlot(covid_gex_patients_not_integrated[[1]],reduction = "tsne")

# UMAP for the second patient
Seurat::DimPlot(covid_gex_patients_not_integrated[[2]],reduction = "umap")


## ---- fig.show='hold',message=FALSE-------------------------------------------
directory_to_covid_patients_integrated <- list()
directory_to_covid_patients_integrated[[1]] <- c("~/Downloads/PlatypusTestData/Patient1_GEX/",
                                                 "~/Downloads/PlatypusTestData/Patient2_GEX/")

## Here we use the previous version of automate_GEX to produce the identical UMAP seen in the manuscript. 
covid_gex <- Platypus::automate_GEX(GEX.outs.directory.list = directory_to_covid_patients_integrated[1:1],integration.method = "scale.data",mito.filter = 20,cluster.resolution = 0.5,VDJ.gene.filter = T)

### In this case, patient1 is the covid_gex[[1]]$sample_id==1 and patient2 is the covid_gex[[1]]$sample_id==2.

class(covid_gex) ## list
length(covid_gex) ### list of length one
class(covid_gex[[1]]) ## Seurat object containing GEX from both samples 


## -----------------------------------------------------------------------------
print(table(covid_gex[[1]]$sample_id)) 

## ---- fig.show='hold'---------------------------------------------------------
Seurat::DimPlot(covid_gex[[1]],reduction = "umap")


## ---- fig.show='hold'---------------------------------------------------------
Seurat::DimPlot(covid_gex[[1]],reduction = "umap",split.by = "sample_id") 

Seurat::DimPlot(covid_gex[[1]],reduction = "umap")

## -----------------------------------------------------------------------------
Platypus::GEX_cluster_membership(automate_GEX.output = covid_gex[[1]])


## ---- fig.show='hold'---------------------------------------------------------
#Seurat::FeaturePlot(covid_gex[[1]],reduction = "umap",features = c("CD4","CD8A","CD19"))
Seurat::FeaturePlot(covid_gex[[1]],reduction = "umap",features = c("CD4"))

## ---- results='hide'----------------------------------------------------------
covid_gex_phenotype <- Platypus::GEX_phenotype(covid_gex[[1]], default = T)

covid_gex[[1]] <- Platypus::GEX_phenotype(covid_gex[[1]], default = F,
                                            cell.state.markers=c("CD8A+;CCL5+;CD44+;IL7R-;CD19-",
                          "CD8A+;CCL5-;CD44+;IL7R+;CD19-"),
                        cell.state.names=c("EffectorCD8",
                        "MemoryCD8"))



## -----------------------------------------------------------------------------
Seurat::DimPlot(covid_gex[[1]],reduction = "umap", group.by = "cell.state") 



## ---- fig.show='hold'---------------------------------------------------------
## Warning: running this function will take a while
gene_expression_cluster <- Platypus::GEX_cluster_genes(covid_gex[[1]],min.pct = 0.25) 

length(gene_expression_cluster) ## length of 12, corresponding to 12 clusters
length(unique(covid_gex[[1]]$seurat_clusters)) ## length of 12 

## -----------------------------------------------------------------------------

head(gene_expression_cluster[[1]])


## -----------------------------------------------------------------------------
print(sapply(gene_expression_cluster,nrow))

## ---- fig.show='hold'---------------------------------------------------------

covid_heatmap_clusters <- Platypus::GEX_cluster_genes_heatmap1(automate_GEX.output = covid_gex[[1]],
                                                              GEX_cluster_genes.output = gene_expression_cluster,
                                                              n.genes.per.cluster = 3,max.cell = 30,
                                                              metric = "top_logFC")

print(covid_heatmap_clusters)

## ----results='hide'-----------------------------------------------------------

ontology_covid <- Platypus::GEX_GOterm(GEX.cluster.genes.output = gene_expression_cluster, topNgenes = 10, go.plots = F)
head(ontology_covid[[1]])

## -----------------------------------------------------------------------------
#top_10_genes_per_cluster <- Platypus::GEX_topN_DE_genes_per_cluster(GEX_cluster_genes.output = gene_expression_cluster, n.genes = 10, by_FC = T)
#head(top_10_genes_per_cluster)


## -----------------------------------------------------------------------------
#gsea_covid <- Platypus::GEX_GSEA(GEX.cluster.genes.output = gene_expression_cluster, MT.Rb.filter = T, path_to_pathways = "~/Downloads/c7.all.v7.2.symbols.gmt")

## ---- fig.show='hold'---------------------------------------------------------
DE_genes_per_sample <- Platypus::GEX_DEgenes_persample(automate.GEX=covid_gex[[1]],min.pct = .25,sample1 = "1",sample2 = "2",return.plot = TRUE,up.genes = 10,down.genes = 10,logFC = F)

head(DE_genes_per_sample) 
nrow(DE_genes_per_sample[[1]]) ##71 rows 


## ---- fig.show='hold'---------------------------------------------------------
#Read in the directories
VDJ.out.directory.list <- list()
VDJ.out.directory.list[[1]] <- "~/Downloads/PlatypusTestData/Patient1_BCR/"
VDJ.out.directory.list[[2]] <- "~/Downloads/PlatypusTestData/Patient2_BCR/"

#Run VDJ_analyze
covid_vdj_repertoire_bcells <- Platypus::VDJ_analyze(VDJ.out.directory =VDJ.out.directory.list, filter.1HC.1LC = T) 



length(covid_vdj_repertoire_bcells) ## list of length two, first element corresponds to the clones from the first repertoire directory that was set as input. 


VDJ.out.directory.list_TCR <- list()
VDJ.out.directory.list_TCR[[1]] <- "~/Downloads/PlatypusTestData/Patient1_TCR/"
VDJ.out.directory.list_TCR[[2]] <- "~/Downloads/PlatypusTestData/Patient2_TCR/"

#Run VDJ_analyze
covid_vdj_repertoire_tcells <- Platypus::VDJ_analyze(VDJ.out.directory =VDJ.out.directory.list_TCR, filter.1HC.1LC = T) 

## -----------------------------------------------------------------------------
print(colnames(covid_vdj_repertoire_bcells[[1]]))

## ---- fig.show='hold'---------------------------------------------------------

covid_vdj_aminoacid_clonotype <- Platypus::VDJ_clonotype(clonotype.list=covid_vdj_repertoire_bcells,
                                                         clone.strategy = "cdr3.aa")

covid_vdj_germline_clonotype <- Platypus::VDJ_clonotype(clonotype.list=covid_vdj_repertoire_bcells,
                                                         clone.strategy = "hvj.lvj")

length(covid_vdj_aminoacid_clonotype) # length of 2
print(colnames(covid_vdj_aminoacid_clonotype[[1]])) 

nrow(covid_vdj_aminoacid_clonotype[[1]]) # 2296 unique clones in patient 1 when clonotyping by amino acid sequence

nrow(covid_vdj_germline_clonotype[[1]]) # 1956 unique clones in patient 1 when clonotyping by germline gene identity


## ---- fig.show='hold'---------------------------------------------------------
covid_single_cell <- Platypus::VDJ_per_clone(clonotype.list = covid_vdj_repertoire_bcells,VDJ.out.directory =VDJ.out.directory.list)


print(paste("There are",length(covid_single_cell[[1]]),"unique nucleotide B cell clones in patient1"),sep="")
print(paste("There are",nrow(covid_single_cell[[1]][[1]]),"unique B cells in the most abundant clone in patient1"),sep="")
print(paste("There are",nrow(covid_single_cell[[2]][[1]]),"unique B cells in the most abundant clone in patient2"),sep="")

print(colnames(covid_single_cell[[1]][[1]]))

covid_single_clone_tcells <- Platypus::VDJ_per_clone(clonotype.list = covid_vdj_repertoire_tcells,VDJ.out.directory =VDJ.out.directory.list_TCR)


## ---- fig.show='hold',message=FALSE-------------------------------------------

### WARNING: You will need to download MiXCR and change the mixcr.directory to the location of MiXCR
covid_vdj_region <- Platypus::call_MIXCR(VDJ.per.clone = covid_single_cell,mixcr.directory = "~/Downloads/mixcr-3.0.12/mixcr",species = "hsa")

print(length(covid_vdj_region[[1]]))
print(colnames(covid_vdj_region[[1]][[1]]))

print(nchar(covid_vdj_region[[1]][[1]]$VDJ.AA.HC[2])) 

print(nchar(covid_vdj_region[[1]][[1]]$full_HC_germline[1]))


## ---- fig.show='hold'---------------------------------------------------------

extracted_covid_germline <- Platypus::VDJ_extract_germline(VDJ.per.clone=covid_single_cell,mixcr.directory="~/Downloads/mixcr-3.0.12/mixcr",extract.VDJRegion=T,species = "hsa")

print(colnames(extracted_covid_germline[[1]][[1]])) ## column names of the germlines from the repertoire corresponding to patient1. 
print(nrow(extracted_covid_germline[[1]][[1]])) ## germline sequences corresponding to 2298 clones in the first patient 
print(nrow(extracted_covid_germline[[2]][[1]])) ## germline sequences corresponding to 2298 clones in the second patient. 


print(extracted_covid_germline[[1]][[1]]$aaSeqCDR3[1]) ## CAREL_FDYW - we can observe unproductive CDR3s for the germline sequence of the first clonotype 



print((extracted_covid_germline[[1]][[1]]$descrsR1[2]))  ## second row in the dataframe corresponds to clonotype4.
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.HC.LC[2])) ### pasted HC and LC together
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.HC[2])) ### just heavy chain for clonotype4 germline
print(nchar(extracted_covid_germline[[1]][[1]]$VDJ.AA.LC[2])) ### just light chain for clonotype4 germline




## ---- fig.show='hold'---------------------------------------------------------

#covid_single_clone_with_JSON <- VDJ_per_cell(clonotype.list = covid_vdj_repertoire_bcells,VDJ.out.directory =VDJ.out.directory.list,JSON = T)


# covid_single_cell was the object from VDJ_per_clone

#print(covid_single_clone_with_JSON[[1]][[1]]$trimmed_HC_sequence[1])  # trimmmed sequence 
#print(covid_single_clone_with_JSON[[1]][[1]]$trimmed_HC_germline[1])  # trimmmed germline sequence 

#stringdist::stringdist(covid_single_clone_with_JSON[[1]][[1]]$sequence_HC[1], VDJ.per.cell[[1]][[1]]$trimmed_ref_HC[1]) 

## ---- fig.show='hold'---------------------------------------------------------


covid_clonal_lineages <- Platypus::VDJ_clonal_lineages(call_MIXCR.output=covid_vdj_region, VDJ_extract_germline.output=extracted_covid_germline,as.nucleotide=F,with.germline=T)

print(colnames(covid_clonal_lineages[[1]][[1]])) ## dataframe with the columns Seq and Name. 

print(covid_clonal_lineages[[1]][[1]]$Seq[1])

print(covid_clonal_lineages[[1]][[1]]$Name[1]) #"clonotype3_1_IGHA1_AACCATGAGTGGAGTC-1"
## Here we can see that the above sequence corresponds to the original clonotype3 but is now the first clonal lineage (as seen by the _1_ before the isotype). Furthermore, the isotype of this cell was of the IGHA1. Lastly, we have the barcode of the cell at the end, allowing us to look back at other cell-specific properties if wanted. 


print(covid_clonal_lineages[[1]][[3]]$Name[3]) #"clonotype7_3_IGHA1_AGTGAGGTCGAGAACG-1"
## Again, this was originally the 7th clonotype, and is now the third most clonally expanded lineage ("_3_") of the IGHA1 isotype with given barcode. 

print(tail(covid_clonal_lineages[[1]][[3]]$Name))  ## here at the end of the dataframe the user can also see the germline sequence, which has the "Name" of "germline". 




## ---- fig.show='hold'---------------------------------------------------------


covid_trees <- Platypus::VDJ_tree(clonal.lineages = covid_clonal_lineages,with.germline=T,min.sequences = 5,max.sequences = 30,unique.sequences = T)
plot(covid_trees[[1]][[1]])



## ---- fig.show='hold'---------------------------------------------------------

covid_isotypes <- Platypus::VDJ_isotypes_per_clone(VDJ_clonotype_output = covid_vdj_repertoire_bcells, VDJ_per_clone_output = covid_single_cell, clones = 30)

covid_isotypes_aa_clonotype <- Platypus::VDJ_isotypes_per_clone(VDJ_clonotype_output = covid_vdj_germline_clonotype, VDJ_per_clone_output = covid_single_cell, clones = 30)

covid_isotypes_germline_clonotype <- Platypus::VDJ_isotypes_per_clone(VDJ_clonotype_output = covid_vdj_aminoacid_clonotype, VDJ_per_clone_output = covid_single_cell, clones = 30)

print(covid_isotypes[[1]])



## ---- fig.show='hold'---------------------------------------------------------
print(covid_isotypes[[2]])

## ---- fig.show='hold'---------------------------------------------------------
## Take the output from VDJ_analyze (or subsample as in this case the top 60 clones)
network_clones_covid <- list()
network_clones_covid[[1]] <- covid_vdj_repertoire_bcells[[1]][1:60,]
network_clones_covid[[2]] <- covid_vdj_repertoire_bcells[[2]][1:60,]
network_clones_covid <- list()
network_clones_covid[[1]] <- covid_vdj_repertoire_tcells[[1]][1:60,]
network_clones_covid[[2]] <- covid_vdj_repertoire_tcells[[2]][1:60,]

covid_bcell_igraph <- Platypus::VDJ_network(network_clones_covid[1:2],per.sample = F,distance.cutoff = 10,connected = F)
covid_bcell_igraph_d14 <- Platypus::VDJ_network(network_clones_covid[1:2],per.sample = F,distance.cutoff = 14,connected = F)


igraph::plot.igraph(covid_bcell_igraph[[4]],vertex.label=NA,vertex.size=7+(.06*covid_bcell_igraph[[2]]$frequency),vertex.color=factor(covid_bcell_igraph[[2]]$mouse))
igraph::plot.igraph(covid_bcell_igraph_d14[[4]],vertex.label=NA,vertex.size=7+(.06*covid_bcell_igraph_d14[[2]]$frequency),vertex.color=factor(covid_bcell_igraph_d14[[2]]$mouse))




## ---- fig.show='hold'---------------------------------------------------------

#First calculate adjacency matrix for V gene usage
covid_Vgene_usage <- Platypus::VDJ_Vgene_usage(VDJ.clonotype.output = covid_vdj_repertoire_bcells)
library(pheatmap)

pheatmap::pheatmap(covid_Vgene_usage[[1]],show_rownames = F,show_colnames = F)


print(class(covid_Vgene_usage[[1]]))
print(head(rownames(covid_Vgene_usage[[1]])))
print(head(colnames(covid_Vgene_usage[[1]])))

## -----------------------------------------------------------------------------
covid_Vgene_usage_barplot <- Platypus::VDJ_Vgene_usage_barplot(covid_vdj_repertoire_bcells, HC.gene.number = 10, LC.Vgene = T, LC.gene.number = 10)
print(covid_Vgene_usage_barplot[[1]])


example.vdj.vgene_usage <- Platypus::VDJ_Vgene_usage_stacked_barplot(clonotype.list = covid_vdj_repertoire_bcells, LC.Vgene = F,HC.gene.number = 10, Fraction.HC = 1)
example.vdj.vgene_usage[[1]]


## -----------------------------------------------------------------------------

vj_circos_bcells <- Platypus::VDJ_VJ_usage_circos(covid_vdj_repertoire_bcells[2:2], c.threshold = 1,label.threshold=50,cell.level = T)

vj_circos_tcells <- Platypus::VDJ_VJ_usage_circos(covid_vdj_repertoire_tcells[1:1], c.threshold = 1,label.threshold=50,cell.level = T)




## -----------------------------------------------------------------------------
covid_CDR3_logoplot <- Platypus::VDJ_logoplot(VDJ.object = covid_vdj_repertoire_bcells, length_cdr3 = 25)

## ---- fig.show='hold'---------------------------------------------------------

covid_integrating_clonal_level <- Platypus::VDJ_GEX_integrate(GEX.object = covid_gex[[1]], 
                                                 clonotype.list =  covid_vdj_repertoire_bcells,
                                                 VDJ.per.clone = covid_single_cell,
                                                 clonotype.level = TRUE)
covid_integrating_clonal_level_tcell <- Platypus::VDJ_GEX_integrate(GEX.object = covid_gex[[1]], 
                                                 clonotype.list =  covid_vdj_repertoire_tcells,
                                                 VDJ.per.clone = covid_single_cell,
                                                 clonotype.level = TRUE)

print(head(covid_integrating_clonal_level[[1]]$majority_cluster))
print(head(covid_integrating_clonal_level[[1]]$cluster_membership_percent))
print(head(covid_integrating_clonal_level[[1]]$cell_index))

## ---- fig.show='hold'---------------------------------------------------------

covid_integrating_cell_level <- Platypus::VDJ_GEX_integrate(GEX.object = covid_gex[[1]], 
                                                 clonotype.list =  covid_vdj_repertoire_bcells[1:1],
                                                 VDJ.per.clone = covid_single_cell[1:1],
                                                 clonotype.level = FALSE)


print(head(covid_integrating_cell_level[[1]][[5]]$cluster_membership))
print(head(covid_integrating_cell_level[[1]][[5]]$cell_index))

## ---- fig.show='hold'---------------------------------------------------------

covid_clonotype_clusters_plot <-  Platypus::VDJ_GEX_expansion(GEX.list=covid_gex[[1]],
                                                              VDJ.GEX.integrate.list = covid_integrating_clonal_level,
                                                              highlight.isotype = "None",
                                                              highlight.number=1:10)


covid_clonotype_clusters_plot_tcells <-  Platypus::VDJ_GEX_expansion(GEX.list=covid_gex[[1]],
                                                              VDJ.GEX.integrate.list = covid_integrating_clonal_level_tcell,
                                                              highlight.isotype = "None",
                                                              highlight.number=1:10)

print(covid_clonotype_clusters_plot_tcells[[2]])




## ---- fig.show='hold'---------------------------------------------------------

covid_top10_umap <- Platypus::GEX_visualize_clones(GEX.list=covid_gex,
                     VDJ.GEX.integrate.list=covid_integrating_clonal_level_tcell[2:2],
                     highlight.type="clonotype",
                     highlight.number=1:10,
                     reduction="umap")
print(covid_top10_umap[[1]]) 


## ---- fig.height=20-----------------------------------------------------------

covid_gex_integrate <- Platypus::GEX_clonotype(GEX.object=covid_gex[[1]],VDJ.per.clone = covid_single_clone_tcells)

covid_gex_integrate_bcells <- Platypus::GEX_clonotype(GEX.object=covid_gex[[1]],VDJ.per.clone = covid_single_cell)


GEX_phenotype_per_clone_plot <- Platypus::GEX_phenotype_per_clone(seurat.object = covid_gex_integrate, clonotype.ids= c(1,2,3,4,5))

print(GEX_phenotype_per_clone_plot)

covid_tcell_gene_heatmap <- Platypus::GEX_heatmap(covid_gex_integrate,b.or.t = "t",clone.rank.threshold = 20,sample.index = 1)


covid_bcell_gene_heatmap <- Platypus::GEX_heatmap(covid_gex_integrate_bcells,b.or.t = "b",clone.rank.threshold = 20,sample.index = 1)




## ---- fig.show='hold'---------------------------------------------------------

clonal_lineage_integrate <- Platypus::VDJ_GEX_clonal_lineage_clusters(covid_integrating_cell_level,covid_clonal_lineages[1:1])

print(clonal_lineage_integrate[[1]][[1]]$Name[1])




## ---- fig.show='hold'---------------------------------------------------------

sessionInfo()



