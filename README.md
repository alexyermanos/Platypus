
# Platypus

Platypus is an R toolkit designed to facilitate the data analysis of single-cell immune repertoire sequencing experiments. 

# System requirements

Platypus has been successfully installed on MacOS X (v10.14.6) and Windows 10 Pro (v1909), and used on R versions 4.0.0 and 3.6.1

# Installation

The package can be installed directly from the tar.gz file on this GitHub. Please see the vignette for examples of how the package can be used. 

```{r}

### Removing any previous versions of the package
# First we will ensure that there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Downloading and installing Platypus


# First we need to download the most recent version from the master branch at https://github.com/alexyermanos/Platypus we can install the package using the following command. 
# WARNING: This needs to be replaced with your own directory where the downloaded package is found

# For MacOS ( users it may look like this
#install.packages("~/Downloads/Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# For windows it will likely look something like this. 
# WARNING: You will need to replace YourPCName with your User name for the windows account in the directory. 
# install.packages("C:\Users\YourPCName\Downloads\Platypus_2.0.4.tar.gz", repos = NULL, type="source")

# Now we can load the installed package into the R environment. 
#library(Platypus)

# The individual R functions can additionally be found on the github in the Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment incase not all functions are desired. Similarly, these functions are actively updated and may include more features than the in original tar.gz file. 


### Downloading the test data
# The COVID-19 data (~136 MB size of the zip file) can be found at the following link https://polybox.ethz.ch/index.php/s/fxQJ3NrRSwiPSSo This dataset contains VDJ (separate libraries for B and T cells) and GEX libraries from two convalescent COVID-19 patients.

# After downloading the zip file named "PlatypusTestData.zip", please unzip the file and find the path to the newly formed folder. Typically this will be in the Downloads folder, so the code below should work on MacOS. For windows please uncomment the code and change the user name to match your PC.

```


Platypus uses a number of different R packages which need prior installation. These can be installed either from CRAN, Bioconductor using "BiocManager" or GitHub using "devtools":


## CRAN

### ape  
### circlize
### do
### dplyr
### ggplot2
### ggseqlogo
### igraph
### jsonlite
### phytools 
### reshape2
### seqinr
### Seurat
### scales
### stringdist
### stringr
### tibble
### tidyverse
### useful 
### utils 

Code to install the packages from CRAN:
```{r}
install.packages("ape")
install.packages("circlize")
install.packages("do")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggseqlogo")
install.packages("igraph")
install.packages("jsonlite")
install.packages("phytools")
install.packages("reshape2")
install.packages("seqinr")
install.packages("Seurat")
install.packages("scales")
install.packages("stringdist")
install.packages("stringr")
install.packages("tibble")
install.packages("tidyverse")
install.packages("useful")
install.packages("utils")
```

## GitHub

### harmony 

Code for installing harmony:
```{r}
install.packages("devtools")
library(devtools)
install_github("immunogenomics/harmony")
```


##Bioconductor

### Biostrings  
### org.Mm.eg.db 
### edgeR
### fgsea
### msa 

Code for installing the packages from Bioconductor:
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("Biostrings")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("edgeR")
BiocManager::install("fgsea")
BiocManager::install("msa")
```


Please post on the github page with any questions or if you would like to contribute.
