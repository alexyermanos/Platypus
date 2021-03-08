## ---- fig.show='hold', message=FALSE------------------------------------------

### Removing any previous versions of the package
#First, ensure there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Downloading and installing Platypus

# Download most recent version from master branch at https://github.com/alexyermanos/Platypus We can install the package using the following command. 
# WARNING: This needs to be replaced with your own directory where the downloaded package is found

# For MacOS users it may look like this:
#install.packages("~/Downloads/Platypus_2.0.5.tar.gz", repos = NULL, type="source")

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

covid_heatmap_clusters <- GEX_cluster_genes_heatmap(automate_GEX.output = covid_gex[[1]],
                                                              GEX_cluster_genes.output = gene_expression_cluster,
                                                              n.genes.per.cluster = 3,max.cell = 30,
                                                              metric = "top_logFC")

print(covid_heatmap_clusters)

