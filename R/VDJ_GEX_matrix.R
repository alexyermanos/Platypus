#' Processes both raw VDJ and GEX Cellranger output to compile a single cell level table containing all available information for each cell.
#' @param VDJ.out.directory.list List containing paths to VDJ output directories from cell ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#'@param GEX.out.directory.list Same as VDJ.out.directory, but for GEX. Order of list items must be the same as for VDJ.
#'@param Seurat.in Alternative to GEX.out.directory.list. List of processed (!) seurat objects. Length of the list can either be 1 if VDJ.integrate = TRUE or equal to the length of VDJ.out.directory.list if VDJ.integrate = FALSE. In metadata the column sample_id and group_id must be present. sample_id must contain ids in the format "s1", "s2" ... "sn" and must be matching the order of VDJ.out.directory.list. No processing (i.e. data normalisation and integration) will be done on these objects.They will be returned as part of the VGM and with additioal VDJ data if integrate.VDJ.to.GEX = T. Filtering parameters such as overlapping barcodes, exclude.GEX.not.in.VDJ and exclude.on.cell.state.markers will be applied to the Seurat.in GEX object(s).
#'@param VDJ.combine Boolean. Defaults to TRUE. Whether to integrate repertoires. A sample identifier will be appended to each barcode both in GEX as well as in VDJ. Recommended for all later functions
#'@param GEX.integrate Boolean. Defaults to TRUE. Whether to integrate GEX data. Default settings use the seurat scale.data option to integrate datasets. Sample identifiers will be appended to each barcode both in GEX and VDJ This is helpful when analysing different samples from the same organ or tissue, while it may be problematic when analysing different tissues.
#'@param integrate.GEX.to.VDJ Boolean. defaults to TRUE. Whether to integrate GEX metadata (not raw counts) into the VDJ output dataframe ! Only possible, if GEX.integrate and VDJ.combine are either both FALSE or both TRUE
#'@param integrate.VDJ.to.GEX Boolean. defaults to TRUE. Whether to integrate VDJ data into GEX seurat object as metadata. ! Only possible, if GEX.integrate and VDJ.combine are either both FALSE or both TRUE
#'@param exclude.GEX.not.in.VDJ Boolean. defaults to FALSE. Whether to delete all GEX cell entries, for which no VDJ information is available. Dependent on data quality and sequencing depth this may reduce the GEX cell count by a significant number
#'@param filter.overlapping.barcodes.GEX Boolean. defaults to TRUE. Whether to remove barcodes which are shared among samples in the GEX analysis. Shared barcodes normally appear at a very low rate.
#'@param filter.overlapping.barcodes.VDJ Boolean. defaults to TRUE. Whether to remove barcodes which are shared among samples in the GEX analysis. Shared barcodes normally appear at a very low rate.
#'@param exclude.on.cell.state.markers Character vector. If no input is provided or input is "none", no cells are excluded. Input format should follow: Character vector containing the gene names for each state. ; is used to use multiple markers within a single gene state. Different vector elements correspond to different states. Example: c("CD4+;CD44-","CD4+;IL7R+;CD44+"). All cells which match any of the given states (in the example case any of the 2) are excluded. This is useful in case different and non lymphocyte cells were co-sequenced. It should give the option to e.g. exclude B cells in the analysis of T cells in a dataset.
#'@param get.VDJ.stats Boolean. defaults to TRUE. Whether to generate general statistics table for VDJ repertoires. This is appended as element [[3]] of the output list.
#'@param parallel.processing Character string. Can be "parlapply" for Windows system, "mclapply" for unix and Mac systems or "none" to use a simple for loop (slow!). Default is "none" for compatibility reasons. For the parlapply option the packages parallel, doParallel and the dependency foreach are required
#'@param numcores Number of cores used for parallel processing. Defaults to number of cores available. If you want to chek how many cores are available use the library Parallel and its command detectCores() (Not setting a limit here when running this function on a cluster may cause a crash)
#'@param trim.and.align Boolean. defaults to FALSE. Whether to trim VJ/VDJ seqs, align them to the 10x reference and trim the reference. This is useful to get full sequences for antibody expression or numbers of somatic hypermutations. !Setting this to TRUE significantly increases computational time
#'@param gap.opening.cost Argument passed to Biostrings::pairwiseAlignment during alignment to reference. Defaults to 10
#'@param gap.extension.cost Argument passed to Biostrings::pairwiseAlignment during alignment to reference. Defaults to 4
#' @param integration.method String specifying which data normalization and integration pipeline should be used. Default is "scale.data", which correspondings to the ScaleData function internal to harmony package. 'sct'specifies SCTransform from the Seurat package. "harmony" should be specificied to perform harmony integration. This method requires the harmony package from bioconductor.
#' @param VDJ.gene.filter Logical indicating if variable genes from the b cell receprot and t cell receptor should be removed from the analysis. True is highly recommended to avoid clonal families clustering together.
#' @param mito.filter Numeric specifying which percent of genes are allowed to be composed of mitochondrial genes. This value may require visual inspection and can be specific to each sequencing experiment. Users can visualize the percentage of genes corresponding to mitochondrial genes using the function "investigate_mitochondial_genes".
#' @param norm.scale.factor Scaling factor for the standard Seurat pipeline. Default is set to 10000 as reported in Seurat documentation.
#' @param n.feature.rna Numeric that specifies which cells should be filtered out due to low number of detected genes. Default is set to 0. Seurat standard pipeline uses 2000.
#' @param n.count.rna.min Numeric that specifies which cells should be filtered out due to low RNA count.Default is set to 0. Seurat standard pipeline without VDJ information uses 200.
#' @param n.count.rna.max Numeric that specifies which cells should be filtered out due to high RNA count.Default is set to infinity. Seurat standard pipeline without VDJ information uses 2500.
#' @param n.variable.features Numeric specifying the number of variable features. Default set to 2000 as specified in Seurat standard pipeline.
#' @param cluster.resolution Numeric specifying the resolution that will be supplied to Seurat's FindClusters function. Default is set to 0.5. Increasing this number will increase the number of distinct Seurat clusters. Suggested to examine multiple parameters to ensure gene signatures differentiating clusters remains constant.
#' @param neighbor.dim Numeric vector specifying which dimensions should be supplied in the FindNeighbors function from Seurat. Default input is '1:10'.
#' @param mds.dim Numeric vector specifying which dimensions should be supplied into dimensional reduction techniques in Seurat and Harmony. Default input is '1:10'.
#'@param group.id vector with integers specifying the group membership. c(1,1,2,2) would specify the first two elements of the input VDJ/GEX lists are in group 1 and the third/fourth input elements will be in group 2.
#'@param subsample.barcodes For development purposes only. If set to TRUE the function will run on 100 cells only to increase speeds of debugging
#'@param Data.in Alternative data input for development purposes only
#' @return Single cell matrix including VDJ and GEX info
#' @import Seurat
#' @import ggplot2
#' @import rmarkdown
#' @export
#' @examples
#' \dontrun{
#'
#' VDJ.out.directory.list <- list()
#' VDJ.out.directory.list[[1]] <- c("~/VDJ/S1/")
#' VDJ.out.directory.list[[2]] <- c("~/VDJ/S2/")
#' GEX.out.directory.list <- list()
#' GEX.out.directory.list[[1]] <- c("~/GEX/S1/")
#' GEX.out.directory.list[[2]] <- c("~/GEX/S2/")
#' VDJ_comb_gex <- VDJ_GEX_matrix(
#' VDJ.out.directory.list = VDJ.out.directory.list
#' ,GEX.out.directory.list = GEX.out.directory.list
#' ,GEX.integrate = T
#' ,VDJ.combine = T
#' ,integrate.GEX.to.VDJ = T
#' ,integrate.VDJ.to.GEX = T
#' ,exclude.GEX.not.in.VDJ = F
#' ,filter.overlapping.barcodes.GEX = F
#' ,filter.overlapping.barcodes.VDJ = F
#' ,get.VDJ.stats = T
#' ,parallel.processing = "none"
#' ,subsample.barcodes = F
#' ,trim.and.align = F
#' ,group.id = c(1,2))
#' }
VDJ_GEX_matrix <- function(VDJ.out.directory.list,
                           GEX.out.directory.list,
                           Seurat.in,
                           VDJ.combine,
                           GEX.integrate,
                           integrate.GEX.to.VDJ,
                           integrate.VDJ.to.GEX,
                           exclude.GEX.not.in.VDJ,
                           filter.overlapping.barcodes.GEX,
                           filter.overlapping.barcodes.VDJ,
                           exclude.on.cell.state.markers,
                           get.VDJ.stats,
                           numcores,
                           trim.and.align,
                           gap.opening.cost,
                           gap.extension.cost,
                           parallel.processing,
                           integration.method,
                           VDJ.gene.filter,
                           mito.filter,
                           norm.scale.factor,
                           n.feature.rna,
                           n.count.rna.min,
                           n.count.rna.max,
                           n.variable.features,
                           cluster.resolution,
                           neighbor.dim,
                           mds.dim,
                           subsample.barcodes,
                           group.id,
                           Data.in){ #LAST ONE only for development to reduce runtime on test runs. Will be removed once the function is stable

  #garbage collector to clear some ram
  gc()

  orig_barcode <- NULL
  match_ex_crit <- NULL
  barcode_VDJ_iteration <- NULL
  GEX_automate_single <- NULL
  VDJ_GEX_stats <- NULL
  do <- NULL
  specificity <- NULL
  affinity <- NULL
  gex.metrics.table <- "none" #error avoidance in case only VDJ is provided and stats are requested

  #Sort out the input situation
  if(missing(Data.in)){ #primary
    Data.in <- "none"
  }

  if(missing(Seurat.in)) Seurat.in <- "none"
  if(class(Seurat.in) == "list"){ #alternative processed seurat input
    for(i in 1:length(Seurat.in)){
      if(!"sample_id" %in% names(Seurat.in[[i]]@meta.data) | !"group_id" %in% names(Seurat.in[[i]]@meta.data)) stop("Seurat.in objects need to contain sample_id and group_id columns")
      if(stringr::str_detect(Seurat.in[[i]]@meta.data$sample_id, "s\\d") == F) stop("Seurat.in objects sample_id column needs to follow sample naming scheme: s1 , s2, ... sn")
    }
    GEX.out.directory.list <- "none"
  }

  if(class(Data.in) == "character"){samples.in <- "none"}#samples in for later
  if(class(Data.in) == "character"){ #if Data.in was provided, this would be false as Data.in would be of class list
    #Worst case: no input was provided at all
    if(missing(VDJ.out.directory.list) & missing(GEX.out.directory.list)){ stop("Please provide data input either as a as a list of local paths to VDJ.out.directory.list and/or GEX.out.directory.list or as list of R objects to Data.in (development only)")

      batches <- "none" #for later to know whether batch numbers should be added as a column. Only if Data.in provided
      #VDJ but not GEX local paths provided
    } else if(missing(VDJ.out.directory.list) == F & missing(GEX.out.directory.list)){
      GEX.out.directory.list <- "none"
      samples.paths.VDJ <- paste0(do.call("c",as.list(VDJ.out.directory.list)), collapse = " ; ")
      if(missing(group.id)) group.id <- 1:length(VDJ.out.directory.list)
      batches <- "none"
      #GEX but not VDJ local paths provided
    } else if(missing(VDJ.out.directory.list) & missing(GEX.out.directory.list) == F){
      VDJ.out.directory.list <- "none"
      samples.paths.GEX <- paste0(do.call("c",as.list(GEX.out.directory.list)), collapse = " ; ")
      if(missing(group.id)) group.id <- 1:length(GEX.out.directory.list)
      batches <- "none"
      #Both local paths provided
    } else if(missing(VDJ.out.directory.list) == F & missing(GEX.out.directory.list) == F)
      if(length(VDJ.out.directory.list) != length(GEX.out.directory.list)){stop("Different number of input paths provided for GEX and VDJ. Please revise input")}
    samples.paths.VDJ <- paste0(do.call("c",as.list(VDJ.out.directory.list)), collapse = " ; ")
    samples.paths.GEX <- paste0(do.call("c",as.list(GEX.out.directory.list)), collapse = " ; ")
    if(missing(group.id)) group.id <- 1:length(GEX.out.directory.list)
    batches <- "none"
  } else if(class(Data.in) == "list"){ #This input is prioritized over local paths input
    GEX.out.directory.list <- "none"
    VDJ.out.directory.list <- "none"

    #figure out the input structure and reorder into a list were elements 1-n are all files for each sample
    samples.in <- list()
    samples.paths.VDJ <- c()
    samples.paths.GEX <- c()
    batches <- c()
    for(i in 1:length(Data.in)){ #First level
      #cat(paste0("i",i))
      if(names(Data.in[[i]])[1] == "VDJ"){ #check if we are already on a sample level
        samples.in[[length(samples.in)+1]] <- Data.in[[i]]
        samples.paths.VDJ <- c(samples.paths.VDJ, Data.in[[i]][[3]])
        samples.paths.GEX <- c(samples.paths.GEX, Data.in[[i]][[4]])
        batches <- c(batches, Data.in[[i]][[5]])
        Data.in[[i]] <- "None" #to limit ram usage
      } else {
        for(j in 1:length(Data.in[[i]])){ #Second level
          #cat(paste0("j",j))
          if(names(Data.in[[i]][[j]])[1] == "VDJ"){ #check if we are a sample level
            samples.in[[length(samples.in)+1]] <- Data.in[[i]][[j]]
            samples.paths.VDJ <- c(samples.paths.VDJ, Data.in[[i]][[j]][[3]])
            samples.paths.GEX <- c(samples.paths.GEX, Data.in[[i]][[j]][[4]])
            batches <- c(batches, Data.in[[i]][[j]][[5]])
            Data.in[[i]][[j]] <- "None" #to limit ram usage
          } else { #Data structure does not match expectations
            stop("Provided datastructure does not match necessary input format")
          }
        }
      }
    }
    samples.paths.GEX <- paste0(samples.paths.GEX, collapse = " ; ")
    samples.paths.VDJ <- paste0(samples.paths.VDJ, collapse = " ; ")

    if(missing(group.id)) group.id <- 1:length(samples.in)
    data.in <- NULL #to save some RAM
  }


  if(missing(subsample.barcodes)) subsample.barcodes = F
  if(missing(group.id)){
    if(missing(GEX.out.directory.list)) group.id <- 1:length(VDJ.out.directory.list)
    if(missing(VDJ.out.directory.list)) group.id <- 1:length(GEX.out.directory.list)
    if(missing(VDJ.out.directory.list) == F & missing(GEX.out.directory.list) == F) group.id <- 1:length(VDJ.out.directory.list)
  }

  if(missing(GEX.out.directory.list)) GEX.out.directory.list <- "none"
  if(missing(VDJ.out.directory.list)) VDJ.out.directory.list <- "none"
  if(missing(VDJ.combine)) VDJ.combine <- T
  if(missing(GEX.integrate)) GEX.integrate <- T
  if(missing(parallel.processing)) parallel.processing <- "none"


  if(missing(integrate.GEX.to.VDJ)) integrate.GEX.to.VDJ <- T
  if(missing(integrate.VDJ.to.GEX)) integrate.VDJ.to.GEX <- T
  if(missing(exclude.GEX.not.in.VDJ)) exclude.GEX.not.in.VDJ <- F
  if(missing(filter.overlapping.barcodes.GEX)) filter.overlapping.barcodes.GEX <- T
  if(missing(filter.overlapping.barcodes.VDJ)) filter.overlapping.barcodes.VDJ <- T
  if(missing(get.VDJ.stats)) get.VDJ.stats <- T

  if(missing(trim.and.align)) trim.and.align <- F
  if(missing(gap.opening.cost)) gap.opening.cost <- 10
  if(missing(gap.extension.cost)) gap.extension.cost <- 4

  if(missing(exclude.on.cell.state.markers)) exclude.on.cell.state.markers <- "none"

  if(parallel.processing == "parlapply" | parallel.processing == "mclapply"){
    if(missing(numcores)) numcores <- parallel::detectCores()
    if(numcores > parallel::detectCores()){numCores <- parallel::detectCores()}
  } else{
    numcores <- 1
  }

  #FOR GEX_automate_single
  if(missing(mito.filter)) mito.filter <- 20
  if(missing(VDJ.gene.filter)) VDJ.gene.filter <- T
  if(missing(norm.scale.factor)) norm.scale.factor <- 10000
  if(missing(n.count.rna.min)) n.count.rna.min <- 0
  if(missing(n.count.rna.max)) n.count.rna.max <- Inf
  if(missing(n.feature.rna)) n.feature.rna <- 0
  if(missing(integration.method)) integration.method <- "scale.data"
  if(missing(GEX.integrate)) GEX.integrate <- T
  if(missing(n.variable.features)) n.variable.features <- 2000
  if(missing(cluster.resolution)) cluster.resolution <- .5
  if(missing(neighbor.dim)) neighbor.dim <- 1:10
  if(missing(mds.dim)) mds.dim <- 1:10
  if(GEX.out.directory.list[[1]] != "none" & VDJ.out.directory.list[[1]] != "none"){
    if(length(VDJ.out.directory.list) != length(GEX.out.directory.list)){stop("Different number of paths supplied for VDJ and GEX")}}

  #save runtime parameters for later
  #params <- do.call("rbind", as.list(environment()))
  params <- c(paste0(do.call("c",as.list(samples.paths.VDJ)), collapse = " / "),
              paste0(do.call("c",as.list(samples.paths.GEX)), collapse = " / "),
              VDJ.combine,
              GEX.integrate,
              integrate.GEX.to.VDJ,
              integrate.VDJ.to.GEX,
              exclude.GEX.not.in.VDJ,
              filter.overlapping.barcodes.GEX,
              filter.overlapping.barcodes.VDJ,
              paste0(exclude.on.cell.state.markers,collapse = ";"),
              get.VDJ.stats,
              numcores,
              trim.and.align,
              gap.opening.cost,
              gap.extension.cost,
              parallel.processing,
              integration.method,
              VDJ.gene.filter,
              mito.filter,
              norm.scale.factor,
              n.feature.rna,
              n.count.rna.min,
              n.count.rna.max,
              n.variable.features,
              cluster.resolution,
              paste0(neighbor.dim, collapse = ";"),
              paste0(mds.dim, collapse = ";"),
              subsample.barcodes,
              paste0(group.id, collapse = ";"))
  names(params) <- c("VDJ.out.directory.list",
                     "GEX.out.directory.list",
                     "VDJ.combine",
                     "GEX.integrate",
                     "integrate.GEX.to.VDJ",
                     "integrate.VDJ.to.GEX",
                     "exclude.GEX.not.in.VDJ",
                     "filter.overlapping.barcodes.GEX",
                     "filter.overlapping.barcodes.VDJ",
                     "exclude.on.cell.state.markers",
                     "get.VDJ.stats",
                     "numcores",
                     "trim.and.align",
                     "gap.opening.cost",
                     "gap.extension.cost",
                     "parallel.processing",
                     "integration.method",
                     "VDJ.gene.filter",
                     "mito.filter",
                     "norm.scale.factor",
                     "n.feature.rna",
                     "n.count.rna.min",
                     "n.count.rna.max",
                     "n.variable.features",
                     "cluster.resolution",
                     "neighbor.dim",
                     "mds.dim",
                     "subsample.barcodes",
                     "group.id")
  gex.list <- list()
  out.list <- list() #open list containing matrices


  ########################################################################################### DEF OF FUNCTIONS

  #Gets statistics on VDJ and GEX
  VDJ_GEX_stats_int <- function(clonotype.list,
                                reference.list,
                                annotations.list,
                                contig.table,
                                vdj.metrics,
                                gex.metrics,
                                samples.paths.VDJ){####START VDJ_GEX_stats

    sample.names <- samples.paths.VDJ
    contig.list <- contig.table

    ### VDJ stats - mainly info coming from the annotated contigs csv
    VDJ.stats.list <- list()
    for(k in 1:length(clonotype.list)){

      cat(paste0("\n Starting with ", k, " of ", length(clonotype.list), "...     "))
      VDJ.stats <- c()

      #gsub to be able to process TCRs as well
      contig.list[[k]]$chain <- gsub(pattern = "TRB", replacement = "IGH", contig.list[[k]]$chain)
      contig.list[[k]]$chain <- gsub(pattern = "TRA", replacement = "IGL", contig.list[[k]]$chain)

      clonotype.list[[k]]$cdr3s_aa <- gsub(pattern = "TRB:", replacement = "IGH:", clonotype.list[[k]]$cdr3s_aa)
      clonotype.list[[k]]$cdr3s_aa <- gsub(pattern = "TRA:", replacement = "IGL:", clonotype.list[[k]]$cdr3s_aa)

      #info on sample
      VDJ.stats[length(VDJ.stats)+1] <- sample.names[[k]]
      names(VDJ.stats)[length(VDJ.stats)] <- "Repertoir path"

      VDJ.stats[length(VDJ.stats)+1] <- sample.names[k]
      names(VDJ.stats)[length(VDJ.stats)] <- "Sample name"

      #Get number of unique barcodes
      VDJ.stats[length(VDJ.stats)+1] <- length(unique(contig.list[[k]]$barcode))
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr unique barcodes"

      #generate lookup table with HC and LC counts and stats per barcode
      cat("\n Getting lookup tables...    ")
      #holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "%",width = 100, style = 3, file = "")
      barcodes <- c()
      nr_HC <- c()
      nr_LC <- c()
      is_cell <-c()
      high_confidence <- c()
      productive <- c()
      full_length <- c()
      nr_bar <- 0
      for(j in unique(contig.list[[k]]$barcode)){
        nr_bar <- nr_bar + 1
        #utils::setTxtProgressBar(value = nr_bar/(length(unique(contig.list[[k]]$barcode)) + nrow(clonotype.list[[k]])),pb = holding_bar)
        barcodes <- append(barcodes, j)
        is_cell <- append(is_cell, min(contig.list[[k]]$is_cell[which(contig.list[[k]]$barcode == j)])) #because most barcodes have two contigs (1HC, 1LC), the min function is used. Normally both contigs have the same "quality stats". In case they do not, the min function always chooses FALSE if present.
        high_confidence <- append(high_confidence, min(contig.list[[k]]$high_confidence[which(contig.list[[k]]$barcode == j)]))
        productive <- append(productive, min(contig.list[[k]]$productive[which(contig.list[[k]]$barcode == j)]))
        full_length <- append(full_length, min(contig.list[[k]]$full_length[which(contig.list[[k]]$barcode == j)]))

        nr_HC <- append(nr_HC,stringr::str_count(paste0(contig.list[[k]]$chain[which(contig.list[[k]]$barcode == j)],collapse = ""), "IGH"))
        nr_LC <- append(nr_LC,stringr::str_count(paste0(contig.list[[k]]$chain[which(contig.list[[k]]$barcode == j)],collapse = ""), "IG(K|L)"))
      }
      lookup_stats <- data.frame(barcodes,nr_HC,nr_LC,is_cell,high_confidence,productive,full_length)
      names(lookup_stats) <- c("barcodes","nr_VDJ","nr_VJ","is_cell","high_confidence","productive","full_length")

      #generate lookup table for clonotypes
      clonotype_ids <- c()
      nr_HC <- c()
      nr_LC <- c()
      for(l in 1:nrow(clonotype.list[[k]])){
        #utils::setTxtProgressBar(value = (nr_bar+l)/(length(unique(contig.list[[k]]$barcode)) + nrow(clonotype.list[[k]])),pb = holding_bar)
        clonotype_ids <- append(clonotype_ids, clonotype.list[[k]]$clonotype_id[l])

        nr_HC <- append(nr_HC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IGH:"))
        nr_LC <- append(nr_LC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IG(K|L):"))
      }
      lookup_stats_clono <- data.frame(clonotype_ids,nr_HC,nr_LC)
      names(lookup_stats_clono) <- c("clonotype_ids","nr_VDJ","nr_VJ")

      #close(holding_bar)

      #number of barcodes with
      #is cell == true
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr barcodes is_cell"

      #number of is.cell with 1 HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 1VJ"

      #number of cells with 1 HC and 0 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 0 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 0VJ"

      #number of cells with 0 HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 0 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 0VDJ 1VJ"

      #number of cells with 2 or more HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more VDJ 1VJ"

      #number of cells with 1 HC and 2 or more LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC > 1 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 2 or more VJ"

      #number of cells with 2 or more HC and 2 or more LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC > 1 & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more VDJ 2 or more VJ"

      #number of cells with
      #full length == true
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$full_length == "true" & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells full_length"

      #number of barcodes with
      #productive == true
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$productive == "true" & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells productive"

      #number of barcodes with
      #high_confidence == true
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$high_confidence == "true" & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells high_confidence"

      #number of cells with
      #all three == true
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$is_cell == "true" & lookup_stats$high_confidence == "true" & lookup_stats$productive == "true" & lookup_stats$full_length == "true" & lookup_stats$is_cell == "true",])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells all true"

      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$is_cell == "true" & lookup_stats$high_confidence == "true" & lookup_stats$productive == "true" & lookup_stats$full_length == "true" & lookup_stats$is_cell == "true" & lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 1,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells all true and 1VDJ 1VJ"

      #number of clonotypes
      VDJ.stats[length(VDJ.stats)+1] <- nrow(clonotype.list[[k]])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes"

      #number of clonotypes with exactly 1HC 1LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC == 1 & lookup_stats_clono$nr_LC == 1,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes 1VDJ 1VJ"

      #number of clonotypes with  < 1HC 1LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC + lookup_stats_clono$nr_LC < 2,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes < 1VDJ 1VJ"

      #number of clonotypes with  > 1HC 1LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC + lookup_stats_clono$nr_LC > 2,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes > 1VDJ 1VJ"

      #More to come

      #percentages
      VDJ.stats.perc <- rep(NA, length(VDJ.stats)-2)
      VDJ.stats.perc[1] <- round(as.numeric(VDJ.stats[3]) / as.numeric(VDJ.stats[3]) *100, digits = 1) #for barcode items
      VDJ.stats.perc[c(2:13)] <- round(as.numeric(VDJ.stats[c(4:15)]) / as.numeric(VDJ.stats[4]) *100, digits = 1) #for barcode and is_cell items
      VDJ.stats.perc[c(14:length(VDJ.stats.perc))] <- round(as.numeric(VDJ.stats[c(16:length(VDJ.stats))]) / as.numeric(VDJ.stats[16])*100,digits = 1)  #for clonotype items

      names(VDJ.stats.perc) <- paste0("% ", names(VDJ.stats[3:length(VDJ.stats)])) #set names to VDJ.stats.perc
      VDJ.stats <- c(VDJ.stats, VDJ.stats.perc) #combine vectors

      VDJ.stats.df <- as.data.frame(t(data.frame(VDJ.stats)))
      names(VDJ.stats.df) <- names(VDJ.stats)
      VDJ.stats.list[[k]] <- VDJ.stats.df

    }
    VDJ.stats.all <- do.call(rbind, VDJ.stats.list)

    cat("\n Getting 10x stats    ")

    VDJ.metrics.list <- vdj.metrics
    GEX.metrics.list <- gex.metrics

    tryCatch({
      #for VDJ
      #check lenght
      #add rep itentifier
      for(ij in 1:length(VDJ.metrics.list)){
        VDJ.metrics.list[[ij]]$rep_id <- ij
        for(ik in 1:ncol(VDJ.metrics.list[[ij]])){
          VDJ.metrics.list[[ij]][,ik] <- as.character(VDJ.metrics.list[[ij]][,ik])
        }
      }

      #### meant to make sure that different samples from different cellranger versions will be bound together
      #get ncols
      ncols <- do.call("c", lapply(VDJ.metrics.list, function(x) ncol(x)))
      #check if length if identical
      if(length(unique(ncols)) == 1){
        VDJ.metrics.all <- dplyr::bind_rows(VDJ.metrics.list)
        #if not merge dataframes sequentially to ensure that all information is maintained
      } else {
        for(m in 1:length(VDJ.metrics.list)){
          if(m == 1){
            ab_1 <- as.data.frame(t(VDJ.metrics.list[[1]]))
            ab_1$idents <- rownames(ab_1)
          } else{
            cur_ab <- as.data.frame(t(VDJ.metrics.list[[m]]))
            cur_ab$idents <- rownames(cur_ab)
            ab_1 <- merge(ab_1, cur_ab, by = "idents", all.x = T, all.y = T)
          }
        }
        VDJ.metrics.all <- as.data.frame(t(ab_1)[2:ncol(ab_1),])
        names(VDJ.metrics.all) <- ab_1$idents
      }

      #for GEX
      #this is a rather inefficient routine to match tables with different columns. This is necessary when outputs from different cellranger versions are combined and the summary metics table is different between samples.
      if(class(GEX.metrics.list) == "list"){

        for(ij in 1:length(GEX.metrics.list)){
          GEX.metrics.list[[ij]]$rep_id <- ij
          for(ik in 1:ncol(GEX.metrics.list[[ij]])){
            GEX.metrics.list[[ij]][,ik] <- as.character(GEX.metrics.list[[ij]][,ik])
          }
        }

        #get ncols
        ncols <- do.call("c", lapply(GEX.metrics.list, function(x) ncol(x)))
        #check if length if identical
        if(length(unique(ncols)) == 1){
          GEX.metrics.all <- dplyr::bind_rows(GEX.metrics.list)
          #if not merge dataframes sequentially to ensure that all information is maintained
        } else {
          for(m in 1:length(GEX.metrics.list)){
            if(m == 1){
              ab_1 <- as.data.frame(t(GEX.metrics.list[[1]]))
              ab_1$idents <- rownames(ab_1)
            } else{
              cur_ab <- as.data.frame(t(GEX.metrics.list[[m]]))
              cur_ab$idents <- rownames(cur_ab)
              ab_1 <- merge(ab_1, cur_ab, by = "idents", all.x = T, all.y = T)
            }
          }
          GEX.metrics.all <- as.data.frame(t(ab_1)[2:ncol(ab_1),])
          names(GEX.metrics.all) <- ab_1$idents
        }
        #bind the two
        VDJ.metrics.all <- cbind(VDJ.metrics.all, GEX.metrics.all)
      }
    }, error = function(e){e
      cat(paste0("\n Adding 10x metrix failed: ", e, "      "))})

    VDJ.stats.all <- cbind(VDJ.stats.all, VDJ.metrics.all)

    cat("\n Done with stats    ")
    return(VDJ.stats.all)
  }

  ########################################################################################### STOP VDJ_GEX_stats

  #Processes and integrates GEX datasets
  GEX_automate_single <- function(GEX.list,
                                  GEX.integrate,
                                  integration.method,
                                  VDJ.gene.filter,
                                  mito.filter,
                                  norm.scale.factor,
                                  n.feature.rna,
                                  n.count.rna.min,
                                  n.count.rna.max,
                                  n.variable.features,
                                  cluster.resolution,
                                  neighbor.dim,
                                  mds.dim,
                                  group.id){

    if(missing(GEX.integrate)) GEX.integrate <- T

    if(missing(mito.filter)) mito.filter <- 5
    if(missing(VDJ.gene.filter)) VDJ.gene.filter <- T
    if(missing(norm.scale.factor)) norm.scale.factor <- 10000
    if(missing(n.count.rna.min)) n.count.rna.min <- 0
    if(missing(n.count.rna.max)) n.count.rna.max <- Inf
    if(missing(n.feature.rna)) n.feature.rna <- 0
    if(missing(integration.method)) integration.method <- "scale.data"
    if(missing(GEX.integrate)) GEX.integrate <- T
    if(missing(n.variable.features)) n.variable.features <- 2000
    if(missing(cluster.resolution)) cluster.resolution <- .5
    if(missing(neighbor.dim)) neighbor.dim <- 1:10
    if(missing(mds.dim)) mds.dim <- 1:10

    for(i in 1:length(GEX.list)){
      #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes)
      GEX.list[[i]][["orig_barcode"]] <- as.character(gsub(".*_","",colnames(GEX.list[[i]])))
      GEX.list[[i]][["orig_barcode"]] <- gsub(GEX.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
    }

    if(GEX.integrate == T & length(GEX.list) > 1){ #combine all GEX into one seurat object and add s%number%_ to the FRONT of the barcode
      GEX.merged <- GEX.list[[1]]
      GEX.merged <- SeuratObject::RenameCells(GEX.merged, new.names = paste0("s",1,"_",colnames(GEX.merged)))
      GEX.merged@meta.data$sample_id <- paste0("s",1)
      GEX.merged@meta.data$group_id <- group.id[1]
      for(i in 2:length(GEX.list)){
        GEX.list[[i]] <- SeuratObject::RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",colnames(GEX.list[[i]])))
        GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
        GEX.list[[i]]@meta.data$group_id <- group.id[i]
        GEX.merged <- merge(GEX.merged, y = GEX.list[[i]], add.cell.ids = c("",""))
      }

      GEX.list <- list() #making this into a list item to make the downstream process uniform
      GEX.list[[1]] <- GEX.merged
    } else {
      for(i in 1:length(GEX.list)){ #or do not integrate, but still add the sample identifier to the FRONT of the barcode. This is to make the output uniform and to deal with the possibility of integrating VDJ but not GEX
        GEX.list[[i]] <- SeuratObject::RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",colnames(GEX.list[[i]])))
        #add sample and group ID
        GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
        GEX.list[[i]]@meta.data$group_id <- group.id[i]
      }
    }

    for(i in 1:length(GEX.list)){

      holding_upper_gene_names <- toupper(rownames(GEX.list[[i]]))
      if(VDJ.gene.filter==T){
        antibody_gene_indices <- which(grepl((holding_upper_gene_names),pattern = "^IGHA")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHG")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHM")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHE")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGK")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^JCHAIN")==F&
                                         grepl((holding_upper_gene_names),pattern = "^IGL")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRAV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRAC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TCRG")==F & #New filter for TCR gamma and delta chains
                                         grepl((holding_upper_gene_names),pattern = "^TCRD")==F)

        GEX.list[[i]] <- GEX.list[[i]][antibody_gene_indices,]
      }

      GEX.list[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^MT-") + Seurat::PercentageFeatureSet(GEX.list[[i]], pattern = "^mt-")
      cell.subset.bool <- (GEX.list[[i]]$percent.mt< mito.filter & GEX.list[[i]]$nFeature_RNA >n.feature.rna & GEX.list[[i]]$nCount_RNA > n.count.rna.min & GEX.list[[i]]$nCount_RNA < n.count.rna.max)
      GEX.list[[i]] <- subset(GEX.list[[i]],cells=which(cell.subset.bool==T))

      if(integration.method=="sct"){
        GEX.list[[i]] <- Seurat::SCTransform(GEX.list[[i]],vars.to.regress = "percent.mt")
        GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature=Seurat::VariableFeatures(object = GEX.list[[i]]))
        GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T)
        GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution)
        GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim)
        GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,check_duplicates=F)

      }
      if(integration.method=="harmony"){
        GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
        GEX.list[[i]] <- Seurat::FindVariableFeatures(GEX.list[[i]], selection.method = "vst", nfeatures = n.variable.features)

        all.genes <- rownames(GEX.list[[i]])
        GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = Seurat::VariableFeatures(object = GEX.list[[i]]))
        GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]],verbose=FALSE,feature= Seurat::VariableFeatures(object = GEX.list[[i]]))
        GEX.list[[i]] <- harmony::RunHarmony(GEX.list[[i]], "sample_id")
        GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]],dims=neighbor.dim,verbose = T,reduction = "harmony")
        GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]],resolution = cluster.resolution,reduction = "harmony")
        GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim,reduction = "harmony")
        GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,reduction = "harmony",check_duplicates=F)

      }
      if(integration.method=="scale.data"){
        GEX.list[[i]] <- Seurat::NormalizeData(GEX.list[[i]], normalization.method = "LogNormalize", scale.factor = norm.scale.factor)
        GEX.list[[i]] <- Seurat::FindVariableFeatures(GEX.list[[i]], selection.method = "vst", nfeatures = n.variable.features)
        all.genes <- rownames(GEX.list[[i]])

        GEX.list[[i]] <- Seurat::ScaleData(GEX.list[[i]], features = all.genes)
        GEX.list[[i]] <- Seurat::RunPCA(GEX.list[[i]], features = Seurat::VariableFeatures(object = GEX.list[[i]]))
        GEX.list[[i]] <- Seurat::FindNeighbors(GEX.list[[i]], dims = neighbor.dim)
        GEX.list[[i]] <- Seurat::FindClusters(GEX.list[[i]], resolution = cluster.resolution)
        GEX.list[[i]] <- Seurat::RunUMAP(GEX.list[[i]], dims = mds.dim)

        GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,check_duplicates=F)
      }

    }

    return(GEX.list)
  }

  ####################################################################################### STOP automate_GEX_single

  #Helper function called in VDJ_GEX_matrix. Do not run as standalone!
  #FUN to call in parlapply mclapply or lapply
  barcode_VDJ_iteration <- function(barcodes, contigs, references, annotations, gap.opening.cost, gap.extension.cost,trim.and.align){


    #Get all the info needed to shrink data usage and search times later in the function
    #Filtering out non productive or non full length contigs from cell. This is neccessary, as a cell labled as productive and full length may still have associated contigs not fullfilling these criteria.
    curr.contigs <- contigs[which(contigs$barcode == barcodes & tolower(contigs$is_cell) == "true" & tolower(contigs$high_confidence) == "true" & tolower(contigs$productive) == "true" & tolower(contigs$full_length) == "true"),]
    if(curr.contigs$raw_clonotype_id[1] != ''){
      curr.references <- references[which(stringr::str_detect(names(references), curr.contigs$raw_clonotype_id[1]))]} else {curr.references <- ""}

    #getting the relevant annotations
    curr.annotations <- annotations[stringr::str_detect(annotations$contig_id, barcodes),]

    #set up data structure
    cols <- c("barcode","sample_id","group_id","clonotype_id_10x","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa", "VJ_cdr3s_aa","VDJ_cdr3s_nt", "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain","VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw","VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_trimmed_ref","VJ_trimmed_ref")
    curr.barcode <- stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 1)), cols)

    #fill in information that do not need processing
    #Contig info on light/a and heavy/b chains is put into different columns (see cols)
    #If only one contig is available, the fields of the other are left blank
    #If more than two contigs of one chain (e.g. 2 TRB) are present, the elements will be pasted separated by a ";" into the relevant fields (in the case of TRB, into the Hb columns)

    #Get number of chains
    HC_count <- sum(stringr::str_count(curr.contigs$chain, pattern = "(TRB|IGH)"))
    LC_count <- sum(stringr::str_count(curr.contigs$chain, pattern = "(TRA|IG(K|L))"))

    #In this case we need to make much less effort with pasting together, so we can save time
    if(HC_count == 1 & LC_count == 1){

      if(which(stringr::str_detect(curr.contigs$chain, "(TRA|IG(K|L))")) == 1){ #make row 1 the heavy chain in case it is not already
        curr.contigs <- curr.contigs[c(2,1),]}

      #fill in the pasted info to curr.barcode directly
      curr.barcode$barcode <- curr.contigs$barcode[1]
      curr.barcode$clonotype_id_10x <- curr.contigs$raw_clonotype_id[1]
      curr.barcode$sample_id <- curr.contigs$raw_clonotype_id[1]
      curr.barcode$group_id <- curr.contigs$raw_clonotype_id[1]
      if(stringr::str_detect(curr.contigs$chain[1], "TR") | stringr::str_detect(curr.contigs$chain[2], "TR")){curr.barcode$celltype <- "T cell"
      } else if(stringr::str_detect(curr.contigs$chain[1], "IG") | stringr::str_detect(curr.contigs$chain[2], "IG")) {curr.barcode$celltype <- "B cell"
      } else {curr.barcode$celltype <- "Unkown"}

      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count

      curr.barcode$VDJ_cdr3s_aa <- curr.contigs$cdr3[1]
      curr.barcode$VJ_cdr3s_aa <- curr.contigs$cdr3[2]
      curr.barcode$VDJ_cdr3s_nt <- curr.contigs$cdr3_nt[1]
      curr.barcode$VJ_cdr3s_nt <- curr.contigs$cdr3_nt[2]
      curr.barcode$VDJ_chain_contig <- curr.contigs$contig_id[1]
      curr.barcode$VJ_chain_contig <- curr.contigs$contig_id[2]
      curr.barcode$VDJ_chain <- curr.contigs$chain[1]
      curr.barcode$VJ_chain <- curr.contigs$chain[2]
      curr.barcode$VDJ_vgene <- curr.contigs$v_gene[1]
      curr.barcode$VJ_vgene <- curr.contigs$v_gene[2]
      curr.barcode$VDJ_dgene <- curr.contigs$d_gene[1]
      curr.barcode$VDJ_jgene <- curr.contigs$j_gene[1]
      curr.barcode$VJ_jgene <- curr.contigs$j_gene[2]
      curr.barcode$VDJ_cgene <- curr.contigs$c_gene[1]
      curr.barcode$VJ_cgene <- curr.contigs$c_gene[2]
      curr.barcode$VDJ_raw_consensus_id <- curr.contigs$raw_consensus_id[1]
      curr.barcode$VJ_raw_consensus_id <- curr.contigs$raw_consensus_id[2]

    } else { # this for cells with abberrant chain numbers

      contigs_pasted <- stats::setNames(data.frame(matrix(ncol = ncol(curr.contigs), nrow = length(unique(curr.contigs$chain)))), names(curr.contigs)) #the dataframe may be one or two rows too long, this will not matter / ROW 1 = Heavy chain information / ROW 2 = Light chain information. This order is maintained even if one of the two chains is not present!

      #Heavy/b chain count
      if(HC_count == 1){
        contigs_pasted[1,] <- curr.contigs[stringr::str_detect(curr.contigs$chain, pattern = "(TRB|IGH)"),]
      } else if(HC_count == 0){
        contigs_pasted[1,] <- ""
      } else if(HC_count > 1){
        for(k in 1:ncol(curr.contigs)){
          contigs_pasted[1,k] <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$chain, pattern = "(TRB|IGH)")), k], collapse = ";")
        }
      }
      ### Order of CDRs with multiple chains is determined here

      #Light/a chain count
      if(LC_count == 1){
        contigs_pasted[2,] <- curr.contigs[stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))"),]
      } else if(LC_count == 0){
        contigs_pasted[2,] <- ""
      } else if(LC_count > 1){
        for(k in 1:ncol(curr.contigs)){
          contigs_pasted[2,k]  <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))")),k],collapse = ";")
        }
      }

      #fill in the pasted info to curr.barcode
      curr.barcode$barcode <- curr.contigs$barcode[1]
      curr.barcode$clonotype_id_10x <- curr.contigs$raw_clonotype_id[1]
      curr.barcode$sample_id <- curr.contigs$raw_clonotype_id[1]
      curr.barcode$group_id <- curr.contigs$raw_clonotype_id[1]
      if(stringr::str_detect(contigs_pasted$chain[1], "TR") | stringr::str_detect(contigs_pasted$chain[2], "TR")){curr.barcode$celltype <- "T cell"
      } else if(stringr::str_detect(contigs_pasted$chain[1], "IG") | stringr::str_detect(contigs_pasted$chain[2], "IG")) {curr.barcode$celltype <- "B cell"
      } else {curr.barcode$celltype <- "Unkown"}


      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count

      curr.barcode$VDJ_cdr3s_aa <- contigs_pasted$cdr3[1]
      curr.barcode$VJ_cdr3s_aa <- contigs_pasted$cdr3[2]
      curr.barcode$VDJ_cdr3s_nt <- contigs_pasted$cdr3_nt[1]
      curr.barcode$VJ_cdr3s_nt <- contigs_pasted$cdr3_nt[2]
      curr.barcode$VDJ_chain_contig <- contigs_pasted$contig_id[1]
      curr.barcode$VJ_chain_contig <- contigs_pasted$contig_id[2]
      curr.barcode$VDJ_chain <- contigs_pasted$chain[1]
      curr.barcode$VJ_chain <- contigs_pasted$chain[2]
      curr.barcode$VDJ_vgene <- contigs_pasted$v_gene[1]
      curr.barcode$VJ_vgene <- contigs_pasted$v_gene[2]
      curr.barcode$VDJ_dgene <- contigs_pasted$d_gene[1]
      curr.barcode$VDJ_jgene <- contigs_pasted$j_gene[1]
      curr.barcode$VJ_jgene <- contigs_pasted$j_gene[2]
      curr.barcode$VDJ_cgene <- contigs_pasted$c_gene[1]
      curr.barcode$VJ_cgene <- contigs_pasted$c_gene[2]
      curr.barcode$VDJ_raw_consensus_id <- stringr::str_split(contigs_pasted$raw_consensus_id[1],";",simplify = T)[1]
      curr.barcode$VJ_raw_consensus_id <- stringr::str_split(contigs_pasted$raw_consensus_id[2],";",simplify = T)[1] #Because we may have more than one consensus ID for light chains, we need to get one of them. Because the consensus ids are always the same for different light or heavy chains of the same cell, we can just take the first element of the str_split

    } #end if HC | LC count > 1

    #Now on to the actual sequences
    reference_HC <- curr.references[names(curr.references) == curr.barcode$VDJ_raw_consensus_id]
    reference_LC <- curr.references[names(curr.references) == curr.barcode$VJ_raw_consensus_id]

    #HEAVY CHAIN / TRB
    #CHECK IF THERE IS 1 2 or 0 chains to process
    if(HC_count == 1){

      #extract match
      HC_contig <- which(curr.annotations$contig_id == curr.barcode$VDJ_chain_contig)
      #get sequence
      curr.barcode$VDJ_sequence_nt_raw <- curr.annotations$sequence[HC_contig]
      if(trim.and.align == T){
        #trim sequence
        curr.barcode$VDJ_sequence_nt_trimmed <- substr(curr.annotations$sequence[HC_contig], as.numeric(curr.annotations$temp_start[HC_contig])+1, as.numeric(curr.annotations$temp_end[HC_contig])-1)
        #translate trimmed sequence
        if(nchar(curr.barcode$VDJ_sequence_nt_trimmed) > 1){
          curr.barcode$VDJ_sequence_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VDJ_sequence_nt_trimmed)))
        } else {to_paste_aa <- ""}
        #align to reference and trim reference
        tryCatch({
          if(nchar(curr.barcode$VDJ_sequence_nt_trimmed) > 1){
            alignments <- Biostrings::pairwiseAlignment(curr.barcode$VDJ_sequence_nt_trimmed, as.character(reference_HC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
            curr.barcode$VDJ_trimmed_ref <- as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))]))
          } else {to_paste_ref_trimmed <-  ""}

        }, error=function(e){
          to_paste_ref_trimmed <- "ALIGNMENT ERROR"
        })
      } else {
        curr.barcode$VDJ_sequence_nt_trimmed <- ""
        curr.barcode$VDJ_sequence_aa <- ""
        curr.barcode$VDJ_trimmed_ref <- ""
      }

    } else if(HC_count == 0){

      curr.barcode$VDJ_sequence_nt_raw <- ""
      curr.barcode$VDJ_sequence_nt_trimmed <- ""
      curr.barcode$VDJ_sequence_aa <- ""
      curr.barcode$VDJ_trimmed_ref <- ""

    } else if(HC_count > 1){ #MORE THAN ONE HC
      #from the annotations extract sequence and paste
      #Heavy/b
      to_paste <- c()
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      #looping contigs in annotation
      for(l in 1:nrow(curr.annotations)){
        #looping over Hb contig ids (as there may be more than 1)
        for(c in 1:length(stringr::str_split(curr.barcode$VDJ_chain_contig, ";",simplify = T))){
          #find a match
          if(curr.annotations$contig_id[l] == stringr::str_split(curr.barcode$VDJ_chain_contig, ";",simplify = T)[c]){
            #get sequence
            to_paste <- append(to_paste, curr.annotations$sequence[l])
            if(trim.and.align == T){
              #trim sequence
              to_paste_trimmed <- append(to_paste_trimmed, substr(curr.annotations$sequence[l], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
              #translate trimmed sequence
              if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
              } else {to_paste_aa <- ""}
              #align to reference and trim reference
              tryCatch({
                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_HC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))])))
                } else {
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
                }
              }, error=function(e){
                to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
              })
            } else {
              to_paste_trimmed <- ""
              to_paste_aa <- ""
              to_paste_ref_trimmed <- ""
            }
          }
        }
      }
      curr.barcode$VDJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
      curr.barcode$VDJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VDJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VDJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }

    #Light/a
    if(LC_count == 1){
      #extract match
      LC_contig <- which(curr.annotations$contig_id == curr.barcode$VJ_chain_contig)
      #get sequence
      curr.barcode$VJ_sequence_nt_raw <- curr.annotations$sequence[LC_contig]
      if(trim.and.align == T){
        #trim sequence
        curr.barcode$VJ_sequence_nt_trimmed <- substr(curr.annotations$sequence[LC_contig], as.numeric(curr.annotations$temp_start[LC_contig])+1, as.numeric(curr.annotations$temp_end[LC_contig])-1)
        #translate trimmed sequence
        if(nchar(curr.barcode$VJ_sequence_nt_trimmed) > 1){
          curr.barcode$VJ_sequence_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VJ_sequence_nt_trimmed)))
        } else {to_paste_aa <- ""}
        #align to reference and trim reference
        tryCatch({
          if(nchar(curr.barcode$VJ_sequence_nt_trimmed) > 1){
            alignments <- Biostrings::pairwiseAlignment(curr.barcode$VJ_sequence_nt_trimmed, as.character(reference_LC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
            curr.barcode$VJ_trimmed_ref <- as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))]))
          } else {to_paste_ref_trimmed <-  ""}

        }, error=function(e){
          to_paste_ref_trimmed <- "ALIGNMENT ERROR"
        })
      } else {
        curr.barcode$VJ_sequence_nt_trimmed <- ""
        curr.barcode$VJ_sequence_aa <- ""
        curr.barcode$VJ_trimmed_ref <- ""
      }

    } else if(LC_count == 0){

      curr.barcode$VJ_sequence_nt_raw <- ""
      curr.barcode$VJ_sequence_nt_trimmed <- ""
      curr.barcode$VJ_sequence_aa <- ""
      curr.barcode$VJ_trimmed_ref <- ""

    } else if(LC_count > 1){ #MORE THAN ONE LC

      to_paste <- c()
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      #looping contigs in annotation
      for(l in 1:nrow(curr.annotations)){
        #looping over Hb contig ids (as there may be more than 1)
        for(c in 1:length(stringr::str_split(curr.barcode$VJ_chain_contig, ";",simplify = T))){
          #find a match
          if(curr.annotations$contig_id[l] == stringr::str_split(curr.barcode$VJ_chain_contig, ";",simplify = T)[c]){
            #get sequence
            to_paste <- append(to_paste, curr.annotations$sequence[l])
            if(trim.and.align == T){
              #trim sequence
              to_paste_trimmed <- append(to_paste_trimmed, substr(curr.annotations$sequence[l], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
              #translate trimmed sequence
              if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
              } else {to_paste_aa <- ""}
              #align to reference and trim reference
              tryCatch({
                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_LC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))])))
                } else {
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
                }
              }, error=function(e){
                to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
              })
            } else {
              to_paste_trimmed <- ""
              to_paste_aa <- ""
              to_paste_ref_trimmed <- ""
            }
          }
        }
      }
      curr.barcode$VJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
      curr.barcode$VJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }

    return(curr.barcode)
  }

  ##################################################################################### STOP barcode_VDJ_iteration


  ##### Start of function
  cat("\n Loaded functions")
  cat("\n Loading in data    ")
  print(Sys.time())
  #Read in raw data for both GEX and VDJ
  vdj.loaded <- F
  if(class(samples.in) == "character" & VDJ.out.directory.list[[1]] != "none"){ #No Data.in => procced with loading by VDJ.out.directory.list

    #Remove possible backslash at the end of the input path
    for(k in 1:length(VDJ.out.directory.list)){
      VDJ.out.directory.list[[k]]<-  gsub("/$", "", VDJ.out.directory.list[[k]])
    }

    vdj_load_error <- tryCatch({suppressWarnings({
      VDJ.out.directory_clonotypes <- paste(VDJ.out.directory.list,"/clonotypes.csv",sep="")
      VDJ.out.directory_reference <- paste(VDJ.out.directory.list,"/concat_ref.fasta",sep="")
      VDJ.out.directory_contigs <- paste(VDJ.out.directory.list,"/filtered_contig_annotations.csv",sep="")
      VDJ.out.directory_annotations <- paste(VDJ.out.directory.list,"/all_contig_annotations.json",sep="")
      VDJ.out.directory_metrics <- paste(VDJ.out.directory.list,"/metrics_summary.csv",sep="")

      clonotype.list <- lapply(VDJ.out.directory_clonotypes, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
      reference.list <- lapply(VDJ.out.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))

      contig.table <- lapply(VDJ.out.directory_contigs, function(x) utils::read.csv(x,sep=",",header=T, )) #better using the table format downstream
      metrics.table <- lapply(VDJ.out.directory_metrics, function(x) utils::read.csv(x,sep=",",header=T, ))

      annotations.list <- lapply(VDJ.out.directory_annotations, function(x) jsonlite::read_json(x))


      #change names so that the barcode function does not have to do that
      for(ijl in 1:length(contig.table)){contig.table[[ijl]]$raw_consensus_id <- gsub("_consensus_","_concat_ref_", contig.table[[ijl]]$raw_consensus_id)}

      ## pulls out the three important features: featureRegions, and of featureRegions. Used for trimming where the V region starts and where the C region ends.
      annotations.table <- list()
      for(i in 1:length(annotations.list)){
        #get annotation table to make parlapply function faster
        annotations.table[[i]] <- do.call(BiocGenerics::rbind,lapply( #get relevant entries out
          annotations.list[[i]]
          , function(y){

            if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")]) == 1){
              temp_start <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")][[1]]$contig_match_start
            } else {
              temp_start <- 10000 #This is to cause substr() in trimming to return an empty string
              #Try substr("ABCDE",100,5) to check
            }
            if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")]) == 1){
              temp_end <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")][[1]]$contig_match_end #!
            } else {
              temp_end <- 0 #This is to cause substr() in trimming to return an empty string
            }
            data.frame("contig_id" = y$contig_name,
                       "sequence" = y$sequence,
                       "temp_start" = temp_start,
                       "temp_end" = temp_end
            )}))### returns this dataframe with these four outputs. if you dont have annotations sufficient for both, then you will just an empty character vector. For high confidence cells we have start and stop.
      }


      #convert everything to character
      for(ijk in 1:length(clonotype.list)){
        for(d in 1:ncol(clonotype.list[[ijk]])){
          clonotype.list[[ijk]][,d] <- as.character(clonotype.list[[ijk]][,d])
        }
        for(d in 1:ncol(annotations.table[[ijk]])){
          annotations.table[[ijk]][,d] <- as.character(annotations.table[[ijk]][,d])
        }
        for(d in 1:ncol(contig.table[[ijk]])){
          contig.table[[ijk]][,d] <- as.character(contig.table[[ijk]][,d])
        }
      }

      vdj.loaded <- T
      cat("\n Loaded VDJ data    ")
      print(Sys.time())
    })}, error = function(e){e
      print(e)})
    if(inherits(vdj_load_error,"error")){
      cat("\n Loading VDJ failed")}

  } else if(class(samples.in) == "list"){ #GET info from samples.in
    #get VDJs

    vdj.loaded <- F
    vdj_load_error <- tryCatch({

      clonotype.list <- lapply(samples.in, function(x){return(x[[1]][[1]])})
      reference.list <- lapply(samples.in, function(x){return(x[[1]][[2]])})
      annotations.list <- lapply(samples.in, function(x){return(x[[1]][[3]])})
      contig.table <- lapply(samples.in, function(x){return(x[[1]][[4]])})
      metrics.table <- lapply(samples.in, function(x){return(x[[1]][[5]])})

      #clear VDJ part in object
      for(i in 1:length(samples.in)){
        samples.in[[i]][[1]] <- "loaded"
      }

      #change names so that the barcode function does not have to do that
      for(ijl in 1:length(contig.table)){contig.table[[ijl]]$raw_consensus_id <- gsub("_consensus_","_concat_ref_", contig.table[[ijl]]$raw_consensus_id)}

      ## pulls out the three important features: featureRegions, and of featureRegions. Used for trimming where the V region starts and where the C region ends.
      annotations.table <- list()
      for(i in 1:length(annotations.list)){
        #get annotation table to make parlapply function faster
        annotations.table[[i]] <- do.call(BiocGenerics::rbind,lapply( #get relevant entries out
          annotations.list[[i]]
          , function(y){

            if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")]) == 1){
              temp_start <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")][[1]]$contig_match_start
            } else {
              temp_start <- 10000 #This is to cause substr() in trimming to return an empty string
              #Try substr("ABCDE",100,5) to check
            }
            if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")]) == 1){
              temp_end <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")][[1]]$contig_match_end #!
            } else {
              temp_end <- 0 #This is to cause substr() in trimming to return an empty string
            }
            data.frame("contig_id" = y$contig_name,
                       "sequence" = y$sequence,
                       "temp_start" = temp_start,
                       "temp_end" = temp_end
            )}))### returns this dataframe with these four outputs. if you dont have annotations sufficient for both, then you will just an empty character vector. For high confidence cells we have start and stop.
      }

      #convert everyhting to character
      for(ijk in 1:length(clonotype.list)){
        for(d in 1:ncol(clonotype.list[[ijk]])){
          clonotype.list[[ijk]][,d] <- as.character(clonotype.list[[ijk]][,d])
        }
        for(d in 1:ncol(annotations.table[[ijk]])){
          annotations.table[[ijk]][,d] <- as.character(annotations.table[[ijk]][,d])
        }
        for(d in 1:ncol(contig.table[[ijk]])){
          contig.table[[ijk]][,d] <- as.character(contig.table[[ijk]][,d])
        }
      }

      vdj.loaded <- T
      cat(paste0(" Loaded VDJ data from Data.in    "))
    }, error = function(e){e
      print(e)})
    if(inherits(vdj_load_error,"error")){
      cat("\n Loading VDJ from Data.in failed")}

  } #end loading from Data.in

  ###################################

  #Load in GEX
  gex.loaded <- F
  if(class(samples.in) == "character" & GEX.out.directory.list[[1]] != "none"){ #No Data.in or input = Seurat.in => procceed with loading by GEX.out.directory.list

    #Remove possible backslash at the end of the input path
    for(k in 1:length(GEX.out.directory.list)){
      GEX.out.directory.list[[k]]<-  gsub("/$", "", GEX.out.directory.list[[k]])
    }

    gex_load_error <- tryCatch({suppressWarnings({
      #add the directory identifier
      GEX.out.directory.list.p <- paste(GEX.out.directory.list,"/filtered_feature_bc_matrix",sep="")
      GEX.out.directory.list.metrics <- paste(GEX.out.directory.list,"/metrics_summary.csv",sep="")
      directory_read10x <- lapply(GEX.out.directory.list.p, function(x) Seurat::Read10X(data.dir=x))
      directory_read10x <- lapply(directory_read10x, function(x){
        rownames(x) <- toupper(rownames(x))
        return(x)})
      gex.list <- lapply(directory_read10x, function(x) Seurat::CreateSeuratObject(x))
      gex.metrics.table <- lapply(GEX.out.directory.list.metrics, function(x) utils::read.csv(x,sep=",",header=T, ))

      for(i in 1:length(gex.list)){
        #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes)
        gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
        gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
      }


      gex.loaded <- T
      cat("\n Loaded GEX data    ")
      print(Sys.time())
    })}, error = function(e){e
      print(e)})
    if(inherits(gex_load_error, "error")){
      cat("\n Loading GEX failed")}

  } else if(class(samples.in) == "list"){ #GET info from samples.in
    gex_load_error <- tryCatch({
      #add the directory identifier

      gex.list <- lapply(samples.in, function(x) Seurat::CreateSeuratObject(x[[2]][[1]]))
      gex.metrics.table <- lapply(samples.in, function(x) return(x[[2]][[2]]))

      #clear GEX part in object
      for(i in 1:length(samples.in)){
        samples.in[[i]][[2]][[1]] <- "loaded"
      }

      for(i in 1:length(gex.list)){
        #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes)
        gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
        gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
      }

      gex.loaded <- T
      cat(paste0(Sys.time()," Loaded GEX data from Data.in"))
    }, error = function(e){e
      print(e)})
    if(inherits(gex_load_error, "error")){
      cat("\n Loading GEX from Data.in failed")}
  }

  ###########################################
  #Call VDJ stats

  #This internal VDJ_GEX_stats function takes already processed input from its wrapper. This shortens loading time, reduces RAM usage and allowes to generate stats when data is provided via Data.in
  stats.done <- F
  if(get.VDJ.stats == T & vdj.loaded == T){
    tryCatch({

      out.stats <- VDJ_GEX_stats_int(clonotype.list = clonotype.list,
                                     reference.list = reference.list,
                                     annotations.list = annotations.list,
                                     contig.table = contig.table,
                                     vdj.metrics = metrics.table,
                                     gex.metrics = gex.metrics.table,
                                     samples.paths.VDJ = stringr::str_split(samples.paths.VDJ, " ; ", simplify = T)[1,]) #needing to split sample paths again as they will end up in separate rows in the stats dataframe
      stats.done <- T
      cat("\n Got VDJ GEX stats    ")
      print(Sys.time())

    }, error = function(e){e
      cat("\n VDJ stats failed: ")
      print(e)})
  }

  #Proceed with selecting the barcodes which will be processed
  #Any barcode that is a cell in the VDJ or in the GEX it will be processed
  #because the filtered GEX is loaded, all barcodes from the GEX will be included in any case.
  #If a barcode is present in GEX but is not a cell in VDJ, it will not be passed through VDJ, to save computational time / This could become a feature to toggle in the future?
  #If a barcode is a cell in VDJ but is not present in the filtered GEX, it will be analyzed for VDJ and GEX rows will be left blank.

  if(class(Seurat.in) == "list"){
    gex.list <- Seurat.in
    Seurat.in <- "loaded"
    vdj.loaded == T
    for(i in 1:length(gex.list)){
      gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
      gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
    }
  }

  if(gex.loaded == T & vdj.loaded == T){ #If both VDJ and GEX are available
    barcodes_GEX <- list()
    barcodes_VDJ <- list()
    for(i in 1:length(gex.list)){
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])

      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true" & tolower(contig.table[[i]]$high_confidence) == "true" & tolower(contig.table[[i]]$productive) == "true" & tolower(contig.table[[i]]$full_length) == "true")])
      #barcodes_VDJ[[i]] <- gsub(unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true" & tolower(contig.table[[i]]$high_confidence) == "true" & tolower(contig.table[[i]]$productive) == "true" & tolower(contig.table[[i]]$full_length) == "true")]),pattern = "-1",replacement = "")

      cat(paste0("\n For sample ", i, ": ", length(barcodes_GEX[[i]])," cell assigned barcodes in GEX, ", length(barcodes_VDJ[[i]]), " cell assigned high confidence barcodes in VDJ. Overlap: ", sum(barcodes_GEX[[i]] %in% barcodes_VDJ[[i]])))

      vdj.gex.available <- colnames(gex.list[[i]]) %in% barcodes_VDJ[[i]]
      gex.list[[i]] <- SeuratObject::AddMetaData(gex.list[[i]], vdj.gex.available, col.name = "VDJ.available")
      #remove all barcodes in GEX which are not present in VDJ (defaults to FALSE)
      if(exclude.GEX.not.in.VDJ == T){
        cat("\n Removing all barcodes from GEX, which are not present in VDJ.")

        vdj.gex.available <- colnames(gex.list[[i]]) %in% barcodes_VDJ[[i]]
        gex.list[[i]]@meta.data$VDJ.available <- vdj.gex.available
        gex.list[[i]] <- subset(gex.list[[i]], cells = colnames(gex.list[[i]])[which(gex.list[[i]]$VDJ.available == T)])

        cat(paste0("Removed ", length(vdj.gex.available)-sum(vdj.gex.available), " GEX entries"))
      }
    }
  }

  if(gex.loaded == T & vdj.loaded == F){
    barcodes_GEX <- list()
    for(i in 1:length(gex.list)){
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])
      cat(paste0("\n For sample ", i, ": ", length(barcodes_GEX[[i]])," cell assigned barcodes in GEX"))
    }}
  if(gex.loaded == F & vdj.loaded == T){
    barcodes_VDJ <- list()
    for(i in 1:length(contig.table)){##barcodes_VDJ holds the unique barcodes
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true" & tolower(contig.table[[i]]$high_confidence) == "true" & tolower(contig.table[[i]]$productive) == "true" & tolower(contig.table[[i]]$full_length) == "true")])
      cat(paste0("\n For sample ", i, ": ", length(barcodes_VDJ[[i]]), " cells assigned with high confidence barcodes in VDJ"))
    }
  }

  #remove sample overlapping barcodes in GEX
  if(filter.overlapping.barcodes.GEX == T & gex.loaded == T){
    if(length(gex.list) > 1){
      barcodes_GEX_c <- do.call("c", lapply(gex.list, function(x) x$orig_barcode))
      unique_barcodes <- names(table(barcodes_GEX_c)[table(barcodes_GEX_c) == 1])
      for(i in 1:length(gex.list)){
        gex.list[[i]] <- subset(gex.list[[i]], subset = orig_barcode %in% unique_barcodes)
      }
      cat(paste0("\n Removed a total of ", length(unique(barcodes_GEX_c)) - length(unique_barcodes), " cells with non unique barcodes in GEX"))
    }
  }

  #remove sample overlapping barcodes in VDJ
  if(filter.overlapping.barcodes.VDJ == T & vdj.loaded == T){
    if(length(barcodes_VDJ) > 1){
      barcodes_VDJ_c <- do.call("c", barcodes_VDJ)
      non_unique_barcodes <- names(table(barcodes_VDJ_c)[table(barcodes_VDJ_c) > 1])
      for(i in 1:length(barcodes_VDJ)){
        barcodes_VDJ[[i]] <- barcodes_VDJ[[i]][which(!barcodes_VDJ[[i]] %in% non_unique_barcodes)]
      }
      cat(paste0("\n Removed a total of ", length(non_unique_barcodes), " cells with non unique barcodes in VDJ"))
    }
  }

  #exclude cells based on marker expression
  #handlers copied from GEX_phenotype.
  if(exclude.on.cell.state.markers[1] != "none" & gex.loaded == T){

    #rename to match GEX_phenotype variables
    cell.state.markers <- exclude.on.cell.state.markers

    Cap<-function(x){
      temp<-c()
      for (i in 1:length(x)){
        s <- strsplit(x, ";")[[i]]
        temp[i]<-paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)), sep="", collapse=";")
      }
      return(temp)
    }

    #check for uppercase genes
    is.hum<-any(useful::find.case(rownames(gex.list[[1]]),case="upper"))

    tryCatch({
      #make user input uppercase in case
      if(is.hum==F){
        cell.state.markers<-Cap(cell.state.markers)
      }
      else if(is.hum==T){
        cell.state.markers <- toupper(cell.state.markers)
      }

      #Replaced with gsub due to package compatibility
      #cell.state.markers<-do:Replace(cell.state.markers,from=";", to="&")
      #cell.state.markers<-do:Replace(cell.state.markers,from="\\+", to=">0")
      #cell.state.markers<-do:Replace(cell.state.markers,from="-", to="==0")

      cell.state.markers<-gsub(pattern =";", replacement ="&",cell.state.markers)
      cell.state.markers<-gsub(pattern ="\\+", replacement =">0",cell.state.markers)
      cell.state.markers<-gsub(pattern ="-", replacement ="==0",cell.state.markers)

      #iterate over seurat objects
      for(j in 1:length(gex.list)){
        cmd<-c()
        #open array containing barcodes which match user input i
        for(i in 1:length(cell.state.markers)){
          barcodes_match_ex_crit <- c()
          #build command
          cmd[i]<-paste0("barcodes_match_ex_crit <- WhichCells(gex.list[[j]], slot = 'counts', expression =", cell.state.markers[i],")")
          is.exist<-tryCatch(expr=length(eval(parse(text=cmd[i]))), error = function(x){
            x<-F
            return(x)})
          #if command worked, subset current seurat object
          if(is.exist!=F){
            cells_unfiltered <- ncol(gex.list[[j]])
            gex.list[[j]]$match_ex_crit <- colnames(gex.list[[j]]) %in% barcodes_match_ex_crit
            gex.list[[j]] <- subset(gex.list[[j]], subset = match_ex_crit == T)
            cat(paste0("\n In GEX sample ", j ," excluded ", cells_unfiltered - ncol(gex.list[[j]])," cells based on ", cell.state.markers[i]))
            #If the Gene was not found in seurat object features
          } else{
            cat(paste0("\n In GEX sample ", j ," failed to exclude cells based on ", cell.state.markers[i], "! Check gene name spelling"))
          }
        }
      }
      #larger error callback should be indepentend of user input.
    }, error = function(e){
      cat("\n Exclusion based on cell markers failed")
      print(e)
    })
  }

  if(subsample.barcodes == T & vdj.loaded == T){
    #FORDEV => shorten barcode list to shorten computational time during DEV
    cat("\n Sampling 50 barcodes from all in VDJ per sample")
    for(i in 1:length(barcodes_VDJ)){
      barcodes_VDJ[[i]] <- sample(barcodes_VDJ[[i]],50)
    }
  }

  #VDJ Processing per cell
  if(vdj.loaded == T){
    VDJ.proc.list <- list()
    for(i in 1:length(contig.table)){

      cat(paste0("\n Starting VDJ barcode iteration ", i , " of ", length(contig.table), "...     "))
      print(Sys.time())
      if(parallel.processing == "parlapply" | parallel.processing == "parLapply"){
        #close any open clusters
        doParallel::stopImplicitCluster()
        #open cluster for parallel computing

        cl <- parallel::makeCluster(numcores)
        cat(paste0("\n Started parlapply cluster with ", numcores, " cores"))

        out.VDJ <- parallel::parLapply(cl, barcodes_VDJ[[i]], barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align)
        #close any open clusters
        doParallel::stopImplicitCluster()
        gc() #garbage collection to reduce ram impact

      } else if(parallel.processing == "mclapply"){
        cat(paste0("\n Started mcapply cluster with ", numcores, " cores"))
        out.VDJ <- parallel::mclapply(X = barcodes_VDJ[[i]], FUN = barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align)
      } else {
        out.VDJ <- lapply(barcodes_VDJ[[i]], barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align)
        gc() #garbage collection to reduce ram impact
      }


      #bind list recieved from parLapply
      VDJ.proc <- dplyr::bind_rows(out.VDJ)
      VDJ.proc[VDJ.proc == ";"] <- "" #fix bug, where if two emtpy strings are concatenated, a ";" is left behind.

      #update barcodes
      VDJ.proc$orig_barcode <- VDJ.proc$barcode
      VDJ.proc$barcode <- paste0("s",i,"_",VDJ.proc$barcode)
      VDJ.proc$sample_id <- paste0("s",i)
      VDJ.proc$group_id <- group.id[i]
      VDJ.proc.list[[i]] <- VDJ.proc

      #add frequency column (e.g. all cells in clonotype 2 will have the same entry, that is the number of cells in clonotype 2)
      VDJ.proc.list[[i]]$clonotype_frequency <- NA
      for(d in 1:nrow(clonotype.list[[i]])){
        VDJ.proc.list[[i]]$clonotype_frequency[which(VDJ.proc.list[[i]]$clonotype_id_10x == clonotype.list[[i]]$clonotype_id[d])] <- clonotype.list[[i]]$frequency[d]
      }
      VDJ.proc.list[[i]]$clonotype_frequency <- as.numeric(as.character(VDJ.proc.list[[i]]$clonotype_frequency))

      #Add further columns to fill in in future updates
      VDJ.proc.list[[i]]$specifity <- NA
      VDJ.proc.list[[i]]$affinity <- NA

      cat(paste0("\n \t Done with ", i , " of ", length(contig.table), "     "))
      print(Sys.time())
    }


    VDJ.proc <- VDJ.proc.list
    if (VDJ.combine == T){ #processing all VDJ files together
      VDJ.proc <- dplyr::bind_rows(VDJ.proc)
    }
  }

  #GEX Processing following the SEURAT pipeline
  if(gex.loaded == T & Seurat.in != "loaded"){ #make sure that Seurat.in is not an input list!
    cat("\n Starting GEX pipeline     ")
    print(Sys.time())
    GEX.proc <- GEX_automate_single(GEX.list = gex.list,
                                    GEX.integrate = GEX.integrate,
                                    integration.method = integration.method,
                                    VDJ.gene.filter = VDJ.gene.filter,
                                    mito.filter = mito.filter,
                                    norm.scale.factor = norm.scale.factor,
                                    n.feature.rna = n.feature.rna,
                                    n.count.rna.min = n.count.rna.min,
                                    n.count.rna.max = n.count.rna.max,
                                    n.variable.features = n.variable.features,
                                    cluster.resolution = cluster.resolution,
                                    neighbor.dim = neighbor.dim,
                                    mds.dim = mds.dim,
                                    group.id = group.id)

    for(i in 1:length(GEX.proc)){
      #remove possible extra underscores before barcode
      GEX.proc[[i]] <- SeuratObject::RenameCells(GEX.proc[[i]], new.names = gsub("^_+","",colnames(GEX.proc[[i]])))

    }
    cat("\n Done with GEX pipeline     ")
    print(Sys.time())
  }

  if(vdj.loaded == T & gex.loaded == T){
    cat("\n Integrating VDJ and GEX...")
    #Combine GEX and VDJ
    #Because a lot of cells are present in only one of the two datasets, we do not add the VDJ table as a metadata item to the Seurat object. This would throw out a lot high confidence VDJ entries. Instead we add metadata columns to both datasets indicating the presence of the cell in the opposing dataset. This should make filtering and joining later easy.
    #Multiple cases:
    if(class(VDJ.proc) == "list" & length(GEX.proc) == 1){ #Multiple VDJ, one GEX
      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], colnames(GEX.proc[[1]]) %in% do.call(c ,lapply(VDJ.proc, function(x) x$barcode)), col.name = "VDJ_available")
      for(i in 1:length(VDJ.proc)){
        VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[1]])
      }

      ########################add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){

        #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
        seur_meta <- SeuratObject::FetchData(GEX.proc[[1]],
                                             vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
        names(seur_meta)[2] <- "orig_barcode_GEX"
        seur_meta$barcode <- rownames(seur_meta)
        #merge to VDJ.proc => into each VDJ we add the relevant info
        for(l in 1:length(VDJ.proc)){
          VDJ.proc[[l]] <- merge(VDJ.proc[[l]], seur_meta, by = "barcode", all.x = T, all.y = F)
        }
      }

      ##############add VDJ info to GEX as metadata columns
      if(integrate.VDJ.to.GEX == T){

        seur_meta <- GEX.proc[[1]]@meta.data
        #set barcodes as columns for merging
        seur_meta$barcode <- rownames(seur_meta)
        #merge
        drop <- c("sample_id", "group_id", "orig.ident", "seurat_clusters", "orig_barcode.x", "orig.barcode.y","orig_barcode", "GEX_available")

        #rbind VDJ proc objects and add them together
        VDJ.proc.all <- dplyr::bind_rows(VDJ.proc)

        seur_meta <- merge(seur_meta, VDJ.proc.all[,which(!names(VDJ.proc.all) %in% drop)], by = "barcode", all.x = T) #these columns are excluded here because they already exist in the seurat object

        VDJ.proc.all <- NULL

        #add rownames to new dataframe, otherwise AddMetaData fails
        rownames(seur_meta) <- seur_meta$barcode
        #add metadata
        GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], seur_meta[,c(10:ncol(seur_meta))], col.name = names(seur_meta)[c(10:ncol(seur_meta))])
      }
      #########################

      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]

    }
    if(class(VDJ.proc) == "data.frame" & length(GEX.proc) == 1){ #one VDJ one GEX ################################################landingsite
      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], colnames(GEX.proc[[1]]) %in% VDJ.proc$barcode, col.name = "VDJ_available")
      VDJ.proc$GEX_available <- VDJ.proc$barcode %in% colnames(GEX.proc[[1]])
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]

      ########################add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){

        #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
        seur_meta <- SeuratObject::FetchData(GEX.proc,
                                             vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
        names(seur_meta)[2] <- "orig_barcode_GEX"
        seur_meta$barcode <- rownames(seur_meta)
        #merge to VDJ.proc
        VDJ.proc <- merge(VDJ.proc, seur_meta, by = "barcode", all.x = T, all.y = F)
      }

      ##############add VDJ info to GEX as metadata columns
      if(integrate.VDJ.to.GEX == T){

        seur_meta <- GEX.proc@meta.data
        #set barcodes as columns for merging
        seur_meta$barcode <- rownames(seur_meta)
        #merge
        drop <- c("sample_id", "group_id", "orig.ident", "seurat_clusters", "orig_barcode.x", "orig.barcode.y","orig_barcode", "GEX_available")
        seur_meta <- merge(seur_meta, VDJ.proc[,which(!names(VDJ.proc) %in% drop)], by = "barcode", all.x = T) #these columns are excluded here because they already exist in the seurat object
        #add rownames to new dataframe, otherwise AddMetaData fails
        rownames(seur_meta) <- seur_meta$barcode
        #add metadata
        GEX.proc <- SeuratObject::AddMetaData(GEX.proc, seur_meta[,c(10:ncol(seur_meta))], col.name = names(seur_meta)[c(10:ncol(seur_meta))])
      }
      #########################

    }
    if(class(VDJ.proc) == "data.frame" & length(GEX.proc) > 1){ #one VDJ multiple GEX (improbable...)
      for(i in 1:length(GEX.proc)){
        GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], colnames(GEX.proc[[i]]) %in% VDJ.proc$barcode, col.name = "VDJ_available")
      }
      VDJ.proc$GEX_available <- VDJ.proc$barcode %in% do.call(c ,lapply(GEX.proc, function(x) colnames(x)))

      ########################add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){
        cat("\n Adding data from multiple GEX objects to one VDJ object")
        #grab metadata from all seurat objects
        seur_meta <- lapply(GEX.proc, function(x){SeuratObject::FetchData(x,
                                                                          vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))})
        seur_meta <- dplyr::bind_rows(seur_meta)
        names(seur_meta)[2] <- "orig_barcode_GEX"
        seur_meta$barcode <- rownames(seur_meta)
        #merge to VDJ.proc
        VDJ.proc <- merge(VDJ.proc, seur_meta, by = "barcode", all.x = T, all.y = F)
      }

      ##############add VDJ info to GEX as metadata columns
      if(integrate.VDJ.to.GEX == T){
        cat("\n Integrating VDJ from a single object to multiple GEX objects is not supported")
      }
      #########################


    }
    if(class(VDJ.proc) == "list" & length(GEX.proc) > 1){ #multiple VDJ multiple GEX
      for(i in 1:length(GEX.proc)){
        GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], colnames(GEX.proc[[i]]) %in% VDJ.proc[[i]]$barcode, col.name = "VDJ_available")
        VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[i]])

        ########################add some GEX columns to VDJ table
        if(integrate.GEX.to.VDJ == T){
          cat("\n Integrating multiple VDJ and GEX in a pariwise way. ! VDJ and GEX input paths or data must be in the same order !")
          #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
          seur_meta <- SeuratObject::FetchData(GEX.proc[[i]],
                                               vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
          names(seur_meta)[2] <- "orig_barcode_GEX"
          seur_meta$barcode <- rownames(seur_meta)
          #merge to VDJ.proc
          VDJ.proc[[i]] <- merge(VDJ.proc[[i]], seur_meta, by = "barcode", all.x = T, all.y = F)
        }

        ##############add VDJ info to GEX as metadata columns
        if(integrate.VDJ.to.GEX == T){
          cat("\n Integrating multiple VDJ and GEX in a pariwise way. ! VDJ and GEX input paths or data must be in the same order !")
          seur_meta <- GEX.proc[[i]]@meta.data
          #set barcodes as columns for merging
          seur_meta$barcode <- rownames(seur_meta)
          #merge
          drop <- c("sample_id", "group_id", "orig.ident", "seurat_clusters", "orig_barcode.x", "orig.barcode.y","orig_barcode", "GEX_available")
          seur_meta <- merge(seur_meta, VDJ.proc[[i]][,which(!names(VDJ.proc[[i]]) %in% drop)], by = "barcode", all.x = T) #these columns are excluded here because they already exist in the seurat object
          #add rownames to new dataframe, otherwise AddMetaData fails
          rownames(seur_meta) <- seur_meta$barcode
          #add metadata
          GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], seur_meta[,c(10:ncol(seur_meta))], col.name = names(seur_meta)[c(10:ncol(seur_meta))])
        }
        #########################
      }
    }
    #Output


    #end multiple cases
    ##############
    out.list <- list("VDJ" = VDJ.proc, "GEX" = GEX.proc)

  }
  if(vdj.loaded == T & gex.loaded == F){
    out.list <- list("VDJ" = VDJ.proc, "GEX" = "none")
  }
  if(vdj.loaded == F & gex.loaded == T){
    if(length(GEX.proc) == 1){
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]
    }
    out.list <- list("VDJ" = "none", "GEX" = GEX.proc)
  }

  if(batches[[1]] != "none"){
    if(class(out.list[[1]]) == "data.frame"){ #check that vdj has been completed
      if(length(batches) == length(unique(out.list[[1]]$sample_id))){
        out.list[[1]]$batch_id <- as.character(out.list[[1]]$sample_id)
        for (i in 1:length(unique(out.list[[1]]$batch_id ))){
          out.list[[1]]$batch_id <- gsub(unique(out.list[[1]]$sample_id)[i], batches[i], out.list[[1]]$batch_id)}
      }
    }
    if(class(out.list[[2]]) != "character" & GEX.integrate == T){ #check that gex has been completed and if there is only a single GEX object
      if(length(batches) == length(unique(out.list[[2]]$sample_id))){
        out.list[[2]]$batch_id <- as.character(out.list[[2]]$sample_id)
        for (i in 1:length(unique(out.list[[2]]$batch_id))){
          out.list[[2]]$batch_id <- gsub(unique(out.list[[2]]$sample_id)[i], batches[i], out.list[[2]]$batch_id)}
      }
    }
  } else if(batches[[1]] == "none"){#append the column anyways to keep formatting consistent
    if(class(out.list[[1]]) == "data.frame"){ #check that vdj has been completed
      out.list[[1]]$batches <- "Unspecified"
    }
    if(class(out.list[[2]]) != "character" & GEX.integrate == T){ #check that gex has been completed
      out.list[[2]]$batches <- "Unspecified"
    }
  }

  #adding VDJ stats
  cat("\n Adding VDJ stats...")
  if(stats.done == T){
    out.list[[3]] <- out.stats
  } else {
    out.list[[3]] <- "VDJ stats failed"
  }

  #adding parameter info for reproducibility
  cat("\n Adding runtime params...")
  out.list[[4]] <- params

  #finally add session info
  out.list[[5]] <- utils::sessionInfo()

  #rename for clarity
  names(out.list) <- c("VDJ", "GEX", "VDJ.GEX.stats", "Running params", "sessionInfo")

  cat("\n Done!   ")
  print(Sys.time())
  # if(class(out.list[[2]])=="Seurat"){
  #   out.list[[2]]$sample_id <- as.integer(gsub(pattern = "s",replacement = "",sub("_.*", "", names(out.list[[2]]$orig_barcode))))
  #   out.list[[2]]$group_id <- rep(group.id,table(out.list[[2]]$sample_id))
  # }

  if(class(out.list[[1]])=="data.frame"){ out.list[[1]]$clonotype_id <- out.list[[1]]$clonotype_id_10x}
  return(out.list)
}
