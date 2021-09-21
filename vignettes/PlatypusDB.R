## ----setup, eval = FALSE, include=FALSE---------------------------------------
#  gc()
#  knitr::opts_chunk$set(fig.width=7, fig.height=7)

## ---- eval = FALSE, include=FALSE---------------------------------------------
#  #this is for sourcing all functions under development and is not included in the knitting
#  source_code_dir <- "C:/Dokumente usw/Master/Reddy Lab/1_thesis/3_code/4_DB project/PlatypusDB_admin/R"
#  file_path_vec <- list.files(source_code_dir, full.names = T)
#  for(f_path in file_path_vec){
#    print(f_path)
#    tryCatch({source(f_path)}, error = function(e){e})
#  }
#  gc()
#  knitr::opts_chunk$set(fig.width=7, fig.height=7)

## ---- eval = FALSE, fig.show='hold', message=FALSE----------------------------
#  
#  ### Removing any previous versions of the package
#  #First can ensure that there is no previous version installed locally
#  #detach("package:Platypus", unload=TRUE)
#  #remove.packages("Platypus")
#  
#  ### Dependencies
#  #install.packages("stringr")
#  
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
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  projects_metadata <- PlatypusDB_list_projects()
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  names(projects_metadata)
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  #use View(projects_metadata[["kreiner2021a"]]) for a full overview
#  
#  projects_metadata[["kreiner2021a"]][1:2,]
#  
#  projects_metadata[["kreiner2021b"]][1:2,]
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  projects_metadata[["kreiner2021a"]][5:9,]
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  PlatypusDB_fetch(
#      PlatypusDB.links = c("kreiner2021a//ALL"), #specifing no sample id and downloading all project level data available
#      save.to.disk = F, #Whether to save it to a specified path. Necessary if download size exceeds available RAM
#      load.to.enviroment = T, #Whether to load into the enviroment directly
#      load.to.list = F, #Whether to return a list of loaded objects
#      combine.objects = T, #Whether to combine objects appropriately. Needed here to get the full VDJ_GEX_matrix
#      path.to.save = "~/Downloads/PlatypusDB_downloads/") #optional argument
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  plots <- VDJ_clonal_expansion(VDJ = kreiner2021a__VDJGEXmatrix[[1]], clones = 30, color.by = "seurat_clusters", group.by = "sample_id")
#  
#  plots[[1]][[2]]
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  loaded_list <- PlatypusDB_fetch(
#      PlatypusDB.links = c("kreiner2021a//VDJmatrix"),
#      save.to.disk = F,
#      load.to.enviroment = F,
#      load.to.list = T) #Returns a list containing the downloaded objects
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  plots <- VDJ_clonal_donut(VDJ = loaded_list[[1]][[1]])
#  
#  plots[1]
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  data_for_VGM <- PlatypusDB_fetch(
#                  PlatypusDB.links = c("kreiner2021a/s1/ALL"),
#                  save.to.disk = F,
#                  load.to.enviroment = F,
#                  load.to.list = T,
#                  combine.objects = T)
#  
#  
#  s2_VGM <- VDJ_GEX_matrix(Data.in = data_for_VGM)
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  FeaturePlot(s2_VGM[[2]], features = c("PTPRC","CD4"))
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  downloaded_data <- PlatypusDB_fetch(
#                     PlatypusDB.links = c("kreiner2021a/s1/ALL","kreiner2021a/s3/ALL"),
#                     save.to.disk = F,
#                     load.to.enviroment = F,
#                     load.to.list = T,
#                     combine.objects = T)
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  #load in local data
#  local_data <- PlatypusDB_load_from_disk(
#                VDJ.out.directory.list = list("~/Downloads/Local_CNS_data/VDJ_S3"),
#                GEX.out.directory.list = list("~/Downloads/Local_CNS_data/GEX_S3"))
#    #To process data with Feature barcode technology refer to the Platypus Feature Barcode vignette
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  comb_VGM <- VDJ_GEX_matrix(
#              Data.in = list(downloaded_data, local_data),
#              group.id = c("EAE rMOG","EAE MOG35-55", "18m aged CNS"),
#              integration.method = "harmony")
#  
#  DimPlot(comb_VGM[[2]], group.by = "group_id")
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  downloaded_data <- PlatypusDB_fetch(
#                     PlatypusDB.links = c("kreiner2021a/s1/VDJ","kreiner2021a/s3/VDJ"),
#                     save.to.disk = F,
#                     load.to.enviroment = F,
#                     load.to.list = T,
#                     combine.objects = T)
#  
#  
#  #load in local data
#  local_data <- PlatypusDB_load_from_disk(
#                VDJ.out.directory.list = list("~/Downloads/Local_CNS_data/VDJ_S3"))
#  
#  comb_VGM <- VDJ_GEX_matrix(
#          Data.in = list(downloaded_data, local_data),
#          group.id = c("EAE rMOG","EAE MOG35-55", "18m aged CNS"))
#  
#  plots <- VDJ_Vgene_usage_stacked_barplot(VDJ = comb_VGM[[1]], platypus.version = "v3", HC.gene.number = 10)
#  
#  plots
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  
#  public_clones <- PlatypusDB_find_CDR3s(VDJ.cdr3s.aa = "CMRYGNYWYFDVW" , VJ.cdr3s.aa = "CLQHGESPFTF", projects.to.search = "ALL")
#  
#  head(public_clones[[1]]) #subset of VDJ dataframes with query VDJ CDR3s
#  head(public_clones[[2]]) #subset of VDJ dataframes with query VJ CDR3s
#  head(public_clones[[3]]) #subset of VDJ dataframes with query VDJ AND VJ CDR3s
#  nrow(public_clones[[3]]) #Nr of cells in database containing query VDJ and VJ CDR3s
#  

## ---- eval = FALSE, fig.show='hold'-------------------------------------------
#  sessionInfo()

