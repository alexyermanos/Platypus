
<!-- README.md is generated from README.Rmd. Please edit that file -->

![](https://repository-images.githubusercontent.com/297313954/10e0a180-713e-11eb-9a23-ef93a9d86e8b)

## [==&gt; Platypus and PlatypusDB homepage](https://alexyermanos.github.io/Platypus/index.html "Platypus and PlatypusDB homepage")

# Platypus

Platypus is an R toolkit designed to facilitate the data analysis of
single-cell immune repertoire sequencing experiments. The manuscript
corresponding to Platypus v2 can be found here at Yermanos et al NARGAB
2021 <https://doi.org/10.1093/nargab/lqab023> and the updated manuscript concerning the v3 Platypus ecosystem can be found here: <https://www.biorxiv.org/content/10.1101/2022.09.28.509709v1>

# Ongoing updates in the Platypus pipeline (v3)

Due to the recent changes of the default clonotyping strategy in Cellranger 
(version 5 and version 6) and annotation of the framework and CDR regions 
(version 7) we have rebuild Platypus to revolve around the VDJ\_build function. 
This function creates a dataframe of the repertoire data. The function 
VGM\_build integrates the repertoire and transcriptome information (Seurat object) 
and will serve as the input to all secondary functions in  future iterations of 
the package. The advantage of this is having all repertoire and transcriptome 
information at a per-cell level.

Furthermore we developed PlatypusDB, a publicly available database that
facilitates the download, integration, and analysis of hundred thousands
of single-cells that contain GEX information, VDJ information or both.
With a single line of code, PlatypusDB allows users to either download
and explore available datasets or integrate existing experiments with
their own datasets. Collectively, PlatypusDB serves as a database for
the scientific community interested in exploration of single cell immune
repertoire sequencing experiments from mouse and human.

Stay tuned for updates <https://twitter.com/alexyermanos?lang=en>

## Architecture

![](images/PlatypusV3_abstract.png)

# System requirements

Platypus has been successfully installed on MacOS X (v10.14.6) and
Windows 10 Pro (v1909), and used on R versions 4.4.0, 4.0.0 and 3.6.1

# Installation

Platypus can easily be installed from CRAN or Github. As changes and bugfixes are made regularly, the Github version may be more recent.

Please scroll down for instructions on how to install the necessary dependencies.

``` r

### Removing any previous versions of the package
# First we will ensure that there is no previous version installed locally
#detach("package:Platypus", unload=TRUE)
#remove.packages("Platypus")

### Downloading and installing Platypus from CRAN

install.packages("Platypus")

### Downloading and installing Platypus from Github

# First we need to download the most recent version from the master branch at https://github.com/alexyermanos/Platypus we can install the package using the following command.
# WARNING: This needs to be replaced with your own directory where the downloaded package is found
# For MacOS users it may look like this
install.packages("~/Downloads/Platypus_3.6.0.tar.gz", repos = NULL, type="source")
# For windows it will likely look something like this.
# WARNING: You will need to replace 'YourPCName' with your user name for the windows account in the directory.
install.packages("C:/Users/YourPCName/Downloads/Platypus_3.6.0.tar.gz", repos = NULL, type="source")

# The individual R functions can additionally be found on the github in the Functions branch. Within this branch, there is a folder "R" which contains the individual functions. This can similarly be downloaded and loaded into the R environment incase not all functions are desired. Similarly, these functions are actively updated and may include more features than the in current CRAN version.

```

Platypus uses a number of different R packages, some of which need prior
installation. These can be installed either from CRAN, Bioconductor:

## CRAN

Code to install the packages from CRAN:

``` r
#Essential packages
install.packages("tidyverse")
install.packages("Seurat")
install.packages("utils")

#Optional packages needed for some individual functions
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
install.packages("stringdist")
```

## GitHub

### Harmony (but not required if not using the Harmony integration method)

Code for installing harmony:

``` r
install.packages("devtools")
library(devtools)
install_github("immunogenomics/harmony")
```

\#\#Bioconductor Code for installing the packages from Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("edgeR")
BiocManager::install("fgsea")
BiocManager::install("msa")
BiocManager::install("ggtree")
```

Please post any questions or issues in the respecitve github section:
<https://github.com/alexyermanos/Platypus/issues> Please reach out to us
with any questions or if you would like to contribute.
