
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->
<!-- [![R-CMD-check](https://github.com/alexyermanos/Platypus/workflows/R-CMD-check/badge.svg)](https://github.com/alexyermanos/Platypus/actions) -->

[![license](https://img.shields.io/badge/license-GPL--2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
<!-- badges: end -->


![](https://repository-images.githubusercontent.com/297313954/10e0a180-713e-11eb-9a23-ef93a9d86e8b)

# Platypus 

Platypus is an R toolkit designed to facilitate the data analysis of
single-cell immune repertoire sequencing experiments. The manuscript
corresponding to Platypus v2 can be found here at Yermanos et al NARGAB
2021 <https://doi.org/10.1093/nargab/lqab023>

# Ongoing updates in the Platypus pipeline (v3)

Due to the recent changes of the default clonotyping strategy in
Cellranger (version 5 and version 6) we are currently rebuilding v3 of
Platypus to revolve around the VDJ\_GEX\_matrix function. This function
integrates both repertoire and transcriptome information and will serve
as the input to all secondary functions in future iterations of the
package. The advantage of this is having all repertoire and
transcriptome information at a per-cell level.

The change in clonotyping can be found here -
<https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/clonotyping>

# PlatypusDB

Together with the update to Platypus v3, we also introduce PlatypusDB, a database for , a database with an integrated R component that allows 
the rapid analysis and integration of hundreds of thousands of B and T cells containing both adaptive immune receptor information (VDJ) and single-cell transcriptomes (GEX). 
PlatypusDB both stores raw output files from the commonly used aligner tool cellranger (10x genomics) and also holds the immune-relevant data in the form of an Robject
that can be loaded directly into the R environment without explicitly requiring file download. The user has the ability to
i) download entire experiments from given publications, ii) download individual samples in PlatypusDB, and iii) download and integrate samples on PlatypusDB with their own local samples.
All functions needed are contained within Platypus v3. 

Stay tuned for updates <https://twitter.com/AlexYermanos>

## Architecture

![](docs/images/PlatypusV3_abstract.png)

# System requirements

Platypus has been successfully installed on MacOS X (v10.14.6) and
Windows 10 Pro (v1909), and used on R versions 4.0.0 and 3.6.1

# Installation

The package can be installed directly from the tar.gz file on this
GitHub. Please see the vignette for examples of how the package can be
used.

Please scroll down for instructions on how to install the necessary
dependencies.

``` r
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
```

Platypus uses a number of different R packages, some of which need prior
installation. These can be installed either from CRAN, Bioconductor
using “BiocManager” or GitHub using “devtools”:

## CRAN

Code to install the packages from CRAN:

``` r
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
#install.packages("utils")
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
```

Please post on the github page with any questions or if you would like
to contribute.
