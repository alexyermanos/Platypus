#' Processes both raw VDJ and GEX Cellranger output to compile a single cell level table containing all available information for each cell. If using Feature Barcodes please note the [FB] paragraph in the description and all "FB." parameters
#' @description
#' # This function is designed as a common input to the Platypus pipeline. Integration of datasets as well as VDJ and GEX information is done here. Please check the Platypus V3 vignette for a detailed walkthrough of the output structure. In short: output[[1]] = VDJ table, output[[2]] = GEX Seurat object and output[[3]] = statistics
#' # [FB] Feature barcode (FB) technology is getting increasingly popular, which is why Platypus V3 fully supports their use as sample delimiters. As of V3, Platpyus does not support Cite-seq data natively, also the VDJ_GEX_matrix function is technically capable of loading a Cite-seq matrix and integrating it with VDJ. For details on how to process sequencing data with FB data and how to supply this information to the VDJ_GEX_matrix function, please consult the dedicated vignette on FB data.
#' @param VDJ.out.directory.list List containing paths to VDJ output directories from cell ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires. ! With version 6. of cellranger and the introduction of the multi function, that allows simultaneous processing of VDJ GEX and feature barcodes, the all_contig_annotations.json file was cut out from the output. For trimming sequences, this function therefore uses the consensus_annotations.csv file
#'@param GEX.out.directory.list List of paths pointing at the outs/ directory of each sample. Order of list items must be the same as for VDJ. ! If using downloaded data or if wishing to use the raw_feature_bc_matrix folder, please add the full path here and set GEX.path.as.is to TRUE
#'@param GEX.path.as.is If TRUE the Seurat function Read10x() is called directly on the GEX.out.directory.list input. If FALSE, "/filtered_feature_bc_matrix" (or sample_feature_bc_matrix in case of Cellranger multi data) is appended to the provided path. Defaults to FALSE. This allows for more flexibility, while still allowing users to simply providing Cellranger output directories.
#'@param FB.out.directory.list [FB] List of paths pointing at the outs/ directory of output from the Cellranger counts function which contain Feature barcode counts. Any input will overwrite potential FB data loaded from the GEX input directories. This may be important, if wanting to input unfiltered FB data that will cover also cells in VDJ not present in GEX.
#'@param Data.in Input for R objects from either the PlatypusDB_load_from_disk or the PlatypusDB_fetch function. If provided, input directories should not be specified. If you wish to integrate local and downloaded data, please load them via load_from_disk and fetch and provide as a list (e.g. Data.in = list(load_from_disk.output, fetch.output))
#'@param Seurat.in Alternative to GEX.out.directory.list. List of processed (!) seurat objects. Length of the list can either be 1 if VDJ.integrate = TRUE or equal to the length of VDJ.out.directory.list if VDJ.integrate = FALSE. In metadata the column sample_id and group_id must be present. sample_id must contain ids in the format "s1", "s2" ... "sn" and must be matching the order of VDJ.out.directory.list. No processing (i.e. data normalisation and integration) will be performed on these objects. They will be returned as part of the VGM and with additional VDJ data if integrate.VDJ.to.GEX = T. Filtering parameters such as overlapping barcodes, exclude.GEX.not.in.VDJ and exclude.on.cell.state.markers will be applied to the Seurat.in GEX object(s).
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
#'@param trim.and.align Boolean. Defaults to FALSE. Whether to trim VJ/VDJ seqs, align them to the 10x reference and trim the reference. This is useful to get full sequences for antibody expression or numbers of somatic hypermutations. !Setting this to TRUE significantly increases computational time
#'@param select.excess.chains.by.umi.count Boolean. Defaults to FALSE. There are several methods of dealing with cells containing reads for more than 1VDJ and 1VJ chain. While many analyses just exclude such cells, the VGM is designed to keep these for downstream evaluation (e.g. in VDJ_clonotype). This option presents an evidenced-based way of selectively keeping only one of the present VDJ  and VJ chains each. E.g. if set to TRUE from all present VDJ chains in a cell only the one with the highest UMI count is kept. (in case of same UMI counts, the chain first in the contig files is kept). Idea source: Zhang W et al. Sci Adv. 2021 (10.1126/sciadv.abf5835)
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
#' @param subsample.barcodes For development purposes only. If set to TRUE the function will run on 100 cells only to increase speeds of debugging
#' @param FB.ratio.threshold Numeric. Defaults to 2 Threshold for assignment of feature barcodes by counts. A feature barcode is assigned to a cell if its counts are >10 and if its counts are FB.ratio.threshold-times higher than the counts of the feature barcode with second most counts.
#'@param group.id vector with integers specifying the group membership. c(1,1,2,2) would specify the first two elements of the input VDJ/GEX lists are in group 1 and the third/fourth input elements will be in group 2.
#' @return Single cell matrix including VDJ and GEX info
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
#' VGM <- VDJ_GEX_matrix(
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
                           GEX.path.as.is,
                           FB.out.directory.list,
                           Data.in,
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
                           select.excess.chains.by.umi.count,
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
                           FB.ratio.threshold,
                           group.id){

  orig_barcode <- NULL
  match_ex_crit <- NULL
  barcode_VDJ_iteration <- NULL
  GEX_automate_single <- NULL
  VDJ_GEX_stats <- NULL
  do <- NULL
  specificity <- NULL
  affinity <- NULL
  gex.metrics.table <- "none" #error avoidance in case only VDJ is provided and stats are requested

  ############################################ Sort out the input situation ####
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
            stop("Provided datastructure does not match required input format")
          }
        }
      }
    }
    samples.paths.GEX <- paste0(samples.paths.GEX, collapse = " ; ")
    samples.paths.VDJ <- paste0(samples.paths.VDJ, collapse = " ; ")

    if(missing(group.id)) group.id <- 1:length(samples.in)
    data.in <- NULL #to save some RAM
  }

  if(missing(GEX.path.as.is)) GEX.path.as.is <- F
  if(missing(subsample.barcodes)) subsample.barcodes = F
  if(missing(group.id)){
    if(missing(GEX.out.directory.list)) group.id <- 1:length(VDJ.out.directory.list)
    if(missing(VDJ.out.directory.list)) group.id <- 1:length(GEX.out.directory.list)
    if(missing(VDJ.out.directory.list) == F & missing(GEX.out.directory.list) == F) group.id <- 1:length(VDJ.out.directory.list)
  }

  if(missing(GEX.out.directory.list)) GEX.out.directory.list <- "none"
  if(missing(VDJ.out.directory.list)) VDJ.out.directory.list <- "none"
  if(missing(FB.out.directory.list)) FB.out.directory.list <- "none"
  if(missing(FB.ratio.threshold)) FB.ratio.threshold <- 2
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
  if(missing(select.excess.chains.by.umi.count)) select.excess.chains.by.umi.count <- F
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

  ############################################ Save runtime parameters ####
  #params <- do.call("rbind", as.list(environment())) #Used this function before: caused crashes when run on MAC
  params <- c(paste0(do.call("c",as.list(samples.paths.VDJ)), collapse = " / "),
              paste0(do.call("c",as.list(samples.paths.GEX)), collapse = " / "),
              paste0(do.call("c",as.list(FB.out.directory.list)), collapse = " / "),
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
              select.excess.chains.by.umi.count,
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
                     "FB.out.directory.list",
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
                     "select.excess.chains.by.umi.count",
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


  ############################################ DEF OF FUNCTIONS ####

  #Assignment of Feature Barcodes to cells (returns a dataframe: rownames= barcodes, column 1 = barcodes, column 2 = feature barcode assignments.) Works identically for VDJ and GEX

  pick_max_feature_barcode <- function(bc_df, #input is a dataframe with the columns beeing the count for each feature barcode and the rows being cells. ! No columns allowed apart from numeric count columns
                                       FB.ratio.threshold){

    if(missing(FB.ratio.threshold)) FB.ratio.threshold <- 2

    #convert intput df to list
    bc_out <- as.list(as.data.frame(t(bc_df)))

    #Choose max barcodes and return the index of the valid maximum BC or 1000 (see below why)
    bc_out <- sapply(bc_out, function(x,FB.ratio.threshold){
      y <- sort(x, decreasing = T)
      if(max(y) > 10){
        if((y[1] / y[2]) > FB.ratio.threshold){
          return(which(x == y[1]))
        } else{
          return(1000)
        }
      } else {
        return(1000)
      }
    },FB.ratio.threshold)

    #from the returned indices pick the right Barcode name from the names of the original dataframe
    bc_out <- sapply(bc_out, function(x, bc_names){
      return(bc_names[x])
    }, names(bc_df))
    #because of the return(1000), cells which could not be assigned will have an NA in the final vector. This is replaced here for clarity
    bc_out[is.na(bc_out)] <- "Not assignable"

    #Open dataframe for results
    bc_match <- data.frame("orig_barcode" = rownames(bc_df), "FB_assignment" = bc_out)
    bc_match[,1] <- gsub("-\\d", "", bc_match[,1]) #remove the -1 at the end of the barcode to allow for merging into GEX and VDJ later
    bc_match[,1] <- gsub(".*_","",bc_match[,1]) #remove any symbols before the actual barcode
    rownames(bc_match) <- rownames(bc_df)
    return(bc_match)
  }

  ############################################ STOP pick_max_feature_barcode ####

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

      cat(paste0("Starting with ", k, " of ", length(clonotype.list), "...\n"))
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
      cat("Getting lookup tables... \n")
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
      names(lookup_stats) <- c("barcodes","nr_HC","nr_LC","is_cell","high_confidence","productive","full_length")

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
      names(lookup_stats_clono) <- c("clonotype_ids","nr_HC","nr_LC")

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

    tryCatch({
      VDJ.metrics.all <- "none" #for error catching later
      if(class(vdj.metrics[[1]]) != "character"){

        cat("Getting 10x stats \n")

        VDJ.metrics.list <- vdj.metrics
        GEX.metrics.list <- gex.metrics

      #for VDJ
      #check length
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
    } #End if vdj.metrics[[1]] != "none"

      #for GEX
      #this is a rather inefficient routine to match tables with different columns. This is necessary when outputs from different cellranger versions are combined and the summary metics table is different between samples.
      GEX.metrics.all <- "none" #for error catching later
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
      }

      if(class(VDJ.metrics.all) != "character" & class(GEX.metrics.all) != "character"){
      #bind the two
      VDJ.metrics.all <- cbind(VDJ.metrics.all, GEX.metrics.all)
      } else if (class(VDJ.metrics.all) == "character" & class(GEX.metrics.all) != "character"){ #got only GEX metrics
        VDJ.metrics.all <- GEX.metrics.all
      } else if (class(VDJ.metrics.all) != "character" & class(GEX.metrics.all) == "character"){ #got only VDJ metrics
        #VDJ.metrics.all <- VDJ.metrics.all no reassignment neccessary
      } else { #got none
        VDJ.metrics.all <- "none"
      }

    }, error = function(e){
      cat(paste0("Adding 10x metrix failed \n"))
      print(e)
      VDJ.metrics.all <- "none"})

    if(class(VDJ.metrics.all) != "character"){ #conditional, only if we got at least one of VDJ and GEX 10x metrics
    VDJ.stats.all <- cbind(VDJ.stats.all, VDJ.metrics.all)
    }

    cat("Done with VDJ_GEX_stats \n")
    return(VDJ.stats.all)
  }

  ############################################ STOP VDJ_GEX_stats ####

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
      GEX.list[[1]] <- "none"
      GEX.merged <- SeuratObject::RenameCells(GEX.merged, new.names = paste0("s",1,"_",colnames(GEX.merged)))
      GEX.merged@meta.data$sample_id <- paste0("s",1)
      GEX.merged@meta.data$group_id <- group.id[1]
      for(i in 2:length(GEX.list)){
        GEX.list[[i]] <- SeuratObject::RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",colnames(GEX.list[[i]])))
        GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
        GEX.list[[i]]@meta.data$group_id <- group.id[i]
        GEX.merged <- merge(GEX.merged, y = GEX.list[[i]], add.cell.ids = c("",""))
        GEX.list[[i]] <- "none"
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

  ############################################ STOP automate_GEX_single ####

  #Helper function called in VDJ_GEX_matrix. Do not run as standalone!
  #FUN to call in parlapply mclapply or lapply
  barcode_VDJ_iteration <- function(barcodes,
                                    contigs,
                                    references,
                                    annotations,
                                    gap.opening.cost,
                                    gap.extension.cost,
                                    trim.and.align,
                                    select.excess.chains.by.umi.count){


    #Get all the info needed to shrink data usage and search times later in the function
    #Filtering out non productive or non full length contigs from cell. This is neccessary, as a cell labled as productive and full length may still have associated contigs not fullfilling these criteria.
    curr.contigs <- contigs[which(contigs$barcode == barcodes & tolower(contigs$is_cell) == "true" & tolower(contigs$high_confidence) == "true" & tolower(contigs$productive) == "true" & tolower(contigs$full_length) == "true"),]

    #Get number of chains
    HC_count <- sum(stringr::str_count(curr.contigs$chain, pattern = "(TRB|IGH)"))
    LC_count <- sum(stringr::str_count(curr.contigs$chain, pattern = "(TRA|IG(K|L))"))

    #In case the cell has more than 1 VJ or VDJ chain, and the user decided to not want to have all these chains included in the later table
    #select one of the excessive VDJ chains by umi count
    if(HC_count > 1 & select.excess.chains.by.umi.count == T){
      curr.contigs$umis_HC[stringr::str_detect(curr.contigs$chain, pattern = "(TRB|IGH)")] <- curr.contigs$umis[stringr::str_detect(curr.contigs$chain, pattern = "(TRB|IGH)")] #getting the counts as a new column entry
      curr.contigs <- curr.contigs[c(which(is.na(curr.contigs$umis_HC) == T), which.max(curr.contigs$umis_HC)),] #getting the VDJ with most umis and all VJs
      HC_count <- 1 #resetting the count to not confuse any loops below
    }
    #select one of the excessive VJ chains by umi count
    if(LC_count > 1 & select.excess.chains.by.umi.count == T){
      curr.contigs$umis_LC[stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))")] <- curr.contigs$umis[stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))")]
      curr.contigs <- curr.contigs[c(which(is.na(curr.contigs$umis_LC) == T), which.max(curr.contigs$umis_LC)),] #getting the VDJ with most umis and all VJs
      LC_count <- 1 #resetting the count to not confuse any loops below
    }

    if(curr.contigs$raw_clonotype_id[1] != ''){ #only getting references if clonotype id is present
      curr.references <- references[which(stringr::str_detect(names(references), curr.contigs$raw_clonotype_id[1]))]
    } else {curr.references <- ""}

    #getting the relevant annotations
    #DEPRECATED after cellranger 6.1 and the multi function, we now use the consensus_annotations table as input. getting the relevant lines is different
    #curr.annotations <- annotations[stringr::str_detect(annotations$contig_id, barcodes),]

    #new annotations read in for cellranger 6.1
    if(all(curr.contigs$raw_consensus_id != "") & trim.and.align == T){ #There are contigs that come without a clonotype id. this protects from erros. Also if trim.and.align == F, we can save time by not doing this lookup
    curr.annotations <- annotations[which(annotations$raw_consensus_id %in% curr.contigs$raw_consensus_id),]
    } else {curr.annotations <- ""}

    #set up data structure
    cols <- c("barcode","sample_id","group_id","clonotype_id_10x","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa", "VJ_cdr3s_aa","VDJ_cdr3s_nt", "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain","VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw","VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_trimmed_ref","VJ_trimmed_ref")
    curr.barcode <- stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 1)), cols)

    #fill in information that do not need processing
    #Contig info on light/a and heavy/b chains is put into different columns (see cols)
    #If only one contig is available, the fields of the other are left blank
    #If more than two contigs of one chain (e.g. 2 TRB) are present, the elements will be pasted separated by a ";" into the relevant fields (in the case of TRB, into the Hb columns)

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
      #adding raw sequences directly
      curr.barcode$VDJ_sequence_nt_raw <- curr.contigs$raw_contig[1]
      curr.barcode$VJ_sequence_nt_raw <- curr.contigs$raw_contig[2]

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
      #adding raw sequences directly
      curr.barcode$VDJ_sequence_nt_raw <- contigs_pasted$raw_contig[1]
      curr.barcode$VJ_sequence_nt_raw <- contigs_pasted$raw_contig[2]

    } #end if HC | LC count > 1

    #Now on to the actual sequences
    reference_HC <- curr.references[names(curr.references) == curr.barcode$VDJ_raw_consensus_id]
    reference_LC <- curr.references[names(curr.references) == curr.barcode$VJ_raw_consensus_id]

    #HEAVY CHAIN / TRB
    #CHECK IF THERE IS 1 2 or 0 chains to process
    if(HC_count == 1){

      if(trim.and.align == T){
        tryCatch({
        #find match in annotations
        HC_contig <- which(curr.annotations$raw_consensus_id == curr.barcode$VDJ_raw_consensus_id)
        #trim sequence
        curr.barcode$VDJ_sequence_nt_trimmed <- substr(curr.barcode$VDJ_sequence_nt_raw, as.numeric(curr.annotations$temp_start[HC_contig])+1, as.numeric(curr.annotations$temp_end[HC_contig])-1)
        #translate trimmed sequence
        if(nchar(curr.barcode$VDJ_sequence_nt_trimmed) > 1){
          curr.barcode$VDJ_sequence_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VDJ_sequence_nt_trimmed)))
        } else {to_paste_aa <- ""}
        #align to reference and trim reference

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
      to_paste <- curr.barcode$VDJ_sequence_nt_raw
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      if(trim.and.align == T){
      tryCatch({
      #looping contigs in annotation
      for(l in 1:nrow(curr.annotations)){
        #looping over Hb contig ids (as there may be more than 1)
        for(c in 1:length(stringr::str_split(curr.barcode$VDJ_raw_consensus_id, ";",simplify = T))){
          #find a match
          if(curr.annotations$raw_consensus_id[l] == stringr::str_split(curr.barcode$VDJ_raw_consensus_id, ";",simplify = T)[c]){
              #trim sequence
              to_paste_trimmed <- append(to_paste_trimmed, substr(stringr::str_split(to_paste, ";",simplify = T)[c], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
              #translate trimmed sequence
              if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
              } else {to_paste_aa <- ""}
              #align to reference and trim reference

                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_HC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))])))
                } else {
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
                }

          }
        }
      }
      }, error=function(e){
        to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
      })
            } else {
              to_paste_trimmed <- ""
              to_paste_aa <- ""
              to_paste_ref_trimmed <- ""
            }

      curr.barcode$VDJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VDJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VDJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }

    #Light/a
    if(LC_count == 1){
      if(trim.and.align == T){
        tryCatch({
         #find match in annotations
        LC_contig <- which(curr.annotations$raw_consensus_id == curr.barcode$VJ_raw_consensus_id)
        #trim sequence
        curr.barcode$VJ_sequence_nt_trimmed <- substr(curr.barcode$VJ_sequence_nt_raw, as.numeric(curr.annotations$temp_start[LC_contig])+1, as.numeric(curr.annotations$temp_end[LC_contig])-1)
        #translate trimmed sequence
        if(nchar(curr.barcode$VJ_sequence_nt_trimmed) > 1){
          curr.barcode$VJ_sequence_aa <- as.character(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VJ_sequence_nt_trimmed)))
        } else {to_paste_aa <- ""}
        #align to reference and trim reference

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

      to_paste <- curr.barcode$VJ_sequence_nt_raw
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      if(trim.and.align == T){
        tryCatch({
        #looping contigs in annotation
      for(l in 1:nrow(curr.annotations)){
        #looping over Hb contig ids (as there may be more than 1)
        for(c in 1:length(stringr::str_split(curr.barcode$VJ_raw_consensus_id, ";",simplify = T))){
          #find a match
          if(curr.annotations$raw_consensus_id[l] == stringr::str_split(curr.barcode$VJ_raw_consensus_id, ";",simplify = T)[c]){
              #trim sequence
              to_paste_trimmed <- append(to_paste_trimmed, substr(stringr::str_split(to_paste, ";",simplify = T)[c], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
              #translate trimmed sequence
              if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
              } else {to_paste_aa <- ""}
              #align to reference and trim reference

                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_LC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))])))
                } else {
                  to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
                }
              }
             }
           }
        }, error=function(e){
          to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
        })
      } else {
        to_paste_trimmed <- ""
        to_paste_aa <- ""
        to_paste_ref_trimmed <- ""
      }
      curr.barcode$VJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }

    return(curr.barcode)
  }

  ############################################ STOP barcode_VDJ_iteration ####

  ##### Start of function
  cat("Loaded functions \n")
  cat("Loading in data \n")
  print(Sys.time())
  ############################################ Load VDJ ####
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

      clonotype.list <- lapply(VDJ.out.directory_clonotypes, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
      reference.list <- lapply(VDJ.out.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
      contig.table <- lapply(VDJ.out.directory_contigs, function(x) utils::read.csv(x,sep=",",header=T)) #better using the table format downstream

      #changes for compatibility with the cellranger multi pipeline
      if(all(file.exists(paste(VDJ.out.directory.list,"/metrics_summary.csv",sep="")))){
        VDJ.out.directory_metrics <- paste(VDJ.out.directory.list,"/metrics_summary.csv",sep="")
        metrics.table <- lapply(VDJ.out.directory_metrics, function(x) utils::read.csv(x,sep=",",header=T))
      } else {
        cat("! metrics_summary.csv file not available in at least one of the VDJ input directories. Loading will be skipped \n")
      }

      #NEW annotations read in. Here we use the consensus_annotations.csv instead of the all_contig_annotations.json, which is not present anymore in the cellranger multi output
      if(all(file.exists(paste(VDJ.out.directory.list,"/consensus_annotations.csv",sep="")))){
        VDJ.out.directory_annotations <- paste(VDJ.out.directory.list,"/consensus_annotations.csv",sep="")
        annotations.list <- lapply(VDJ.out.directory_annotations, function(x) utils::read.csv(x,sep=",",header=T))

        #New processing for this consensus annotation table
        annotations.table <- list()
        for(i in 1:length(annotations.list)){
          #get annotation table to make parlapply function faster
          annotations.table[[i]] <- annotations.list[[i]][,c("consensus_id", "cdr3", "v_start", "j_end")]
          names(annotations.table[[i]]) <- c("raw_consensus_id","sequence","temp_start","temp_end")
          annotations.table[[i]]$raw_consensus_id <- gsub("_consensus","_concat_ref_", annotations.table[[i]]$raw_consensus_id) #This is for conformity with other dataframes and the reference
          #!!! Note the missing _ after _consensus pattern! Interestingly in the 10x v6 output this _ is missing. This will probably get fixed in the next update ! => same thing a bit lower in the Data.in input
        }
      } else { #if file not available
        annotations.list <- "none"
        annotations.table <- "none"
        if(trim.and.align == T){
          cat("! consensus_annotations.csv file was not available in at least one of the VDJ input directories. Trimming and aligning will not be run \n")
          trim.and.align <- F
        }
      }


      #NEW in cellranger 6.1 because the consensus_annotations does not contain the contig anymore
      VDJ.out.directory_raw_contigs <- paste(VDJ.out.directory.list,"/filtered_contig.fasta",sep="")
      raw.contig.table <- lapply(VDJ.out.directory_raw_contigs, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))

      for(ikj in 1:length(raw.contig.table)){
        #coercing this to a dataframe
        raw.contig.table[[ikj]] <- data.frame("contig_id" = names(raw.contig.table[[ikj]]), "raw_contig" = as.character(raw.contig.table[[ikj]]))
        #merge directly with contig dataframe to have annotations and sequence in one place for later
        contig.table[[ikj]] <- merge(contig.table[[ikj]], raw.contig.table[[ikj]], by = "contig_id", all.x = T, all.y = F) #making sure to only merge in raw contig sequences for contigs which are present in the table containing annotations
        if(sum(is.na(contig.table[[ikj]]$raw_contig)) > 0.5*nrow(contig.table[[ikj]])){
          cat("! Merging of raw contigs and filtered_contig_annotations showed unsusually low overlap \n")
        }
      }

      #END NEW STUFF after cellranger 6.1

      #change names so that the barcode function does not have to do that
      for(ijl in 1:length(contig.table)){contig.table[[ijl]]$raw_consensus_id <- gsub("_consensus_","_concat_ref_", contig.table[[ijl]]$raw_consensus_id)}

      #DEPRECATED AFTER CELLRANGER V6.1 AND USAGE OF THE CONSENSUS_ANNOTATIONS.CSV AS SPECIFIED ABOVE!
      # pulls out the three important features: featureRegions, and of featureRegions. Used for trimming where the V region starts and where the C region ends.
      #annotations.table <- list()
      #for(i in 1:length(annotations.list)){
        #get annotation table to make parlapply function faster
      #  annotations.table[[i]] <- do.call(BiocGenerics::rbind,lapply( #get relevant entries out
      #    annotations.list[[i]]
      #    , function(y){
      #
      #      if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")]) == 1){
      #        temp_start <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")][[1]]$contig_match_start
      #      } else {
      #        temp_start <- 10000 #This is to cause substr() in trimming to return an empty string
      #        #Try substr("ABCDE",100,5) to check
      #      }
      #      if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")]) == 1){
      #        temp_end <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")][[1]]$contig_match_end #!
      #      } else {
      #        temp_end <- 0 #This is to cause substr() in trimming to return an empty string
      #      }
      #      data.frame("contig_id" = y$contig_name,
      #                 "sequence" = y$sequence,
      #                 "temp_start" = temp_start,
      #                 "temp_end" = temp_end
      #     )}))### returns this dataframe with these four outputs. if you dont have annotations sufficient for both, then you will just an empty character vector. For high confidence cells we have start and stop.
      #}

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
      cat("Loaded VDJ data \n")
      print(Sys.time())
    })}, error = function(e){e
      cat("Loading VDJ failed \n")
      print(e)})

  } else if(class(samples.in) == "list"){ #GET info from samples.in
    #get VDJs

    vdj.loaded <- F
    vdj_load_error <- tryCatch({

      clonotype.list <- lapply(samples.in, function(x){return(x[[1]][[1]])})
      reference.list <- lapply(samples.in, function(x){return(x[[1]][[2]])})
      annotations.list <- lapply(samples.in, function(x){return(x[[1]][[3]])}) #Note: even after update to cellranger 6.1. this line does not change, given that the Platypus_DB_load_from_disk will be dealing with then new input format
      contig.table <- lapply(samples.in, function(x){return(x[[1]][[4]])})
      metrics.table <- lapply(samples.in, function(x){return(x[[1]][[5]])})

      #clear VDJ part in object
      for(i in 1:length(samples.in)){
        samples.in[[i]][[1]] <- "loaded"
      }

      #change names so that the barcode function does not have to do that
      for(ijl in 1:length(contig.table)){
        contig.table[[ijl]]$raw_consensus_id <- gsub("_consensus_","_concat_ref_", contig.table[[ijl]]$raw_consensus_id)
        }

      #New processing for this consensus annotation table
      annotations.table <- list()
      for(ijl in 1:length(annotations.list)){
        #get annotation table to make parlapply function faster
        annotations.table[[ijl]] <- annotations.list[,c("consensus_id", "cdr3", "v_start", "j_end")]
        names(annotations.table[[ijl]]) <- c("raw_consensus_id","sequence","temp_start","temp_end")
        annotations.table[[ijl]]$raw_consensus_id <- gsub("_consensus","_concat_ref_", annotations.table[[ijl]]$raw_consensus_id) #This is for conformity with other dataframes and the reference
        #!!! Note the missing _ after _consensus pattern! Interestingly in the 10x v6 output this _ is missing. This will probably get fixed in the next update ! => same thing a bit above in the directory input
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
      cat(paste0("Loaded VDJ data from Data.in \n"))
    }, error = function(e){
      cat("Loading VDJ from Data.in failed \n")
      print(e)})
  } #end loading from Data.in

  ############################################ Load GEX ####

  #Load in GEX
  gex.loaded <- F
  FB.loaded <- F
  if(class(samples.in) == "character" & GEX.out.directory.list[[1]] != "none"){ #No Data.in or input = Seurat.in => procceed with loading by GEX.out.directory.list

    #Remove possible backslash at the end of the input path
    for(k in 1:length(GEX.out.directory.list)){
      GEX.out.directory.list[[k]]<-  gsub("/$", "", GEX.out.directory.list[[k]])
    }

    gex_load_error <- tryCatch({suppressWarnings({
      #add the directory identifier

      if(stringr::str_detect(GEX.out.directory.list[[1]], "filtered_feature_bc_matrix") | stringr::str_detect(GEX.out.directory.list[[1]], "raw_feature_bc_matrix") | stringr::str_detect(GEX.out.directory.list[[1]], "sample_feature_bc_matrix")| GEX.path.as.is){
        #Nothing to append
        GEX.out.directory.list.p <- GEX.out.directory.list
      } else{

        if(dir.exists(paste(GEX.out.directory.list[[1]],"/filtered_feature_bc_matrix",sep=""))){ #checking only for first path assuming that all sample inputs are processed with the same cellranger function
          GEX.out.directory.list.p <- paste(GEX.out.directory.list,"/filtered_feature_bc_matrix",sep="")
          cat("Setting GEX directory to provided path/filtered_feature_bc_matrix \n")
        } else if (dir.exists(paste(GEX.out.directory.list[[1]],"/sample_feature_bc_matrix",sep=""))){
          GEX.out.directory.list.p <- paste(GEX.out.directory.list,"/sample_feature_bc_matrix",sep="")
          cat("Setting GEX directory to provided path/sample_feature_bc_matrix \n")
        } else {
          stop("The GEX directory filtered_feature_bc_matrix or sample_feature_bc_matrix was not found at the given path. Please revise GEX input paths")
        }
      }

      directory_read10x <- lapply(GEX.out.directory.list.p, function(x) Seurat::Read10X(data.dir=x))

      #NEW CELLRANGER 6.1. dealing with the possibility of Feature barcode information being included in the GEX directory.
      #=> Procedure: check if there is feature barcode info in the GEX that was just loaded => if not, proceed as normal
      #=> if yes, isolate the matrices for each input sample into two flat lists: GEX.list and FB.list

      FB.list <- list()
      #dealing with possible mixed GEX FB inputs or multiple FB input matrices from the same directory
      for(i in 1:length(directory_read10x)){ #iterating over main list
        if(class(directory_read10x[[i]]) == "list"){ #this returns true only if the GEX directory import contained more than one marices => i.e. there is a GEX and a FB matrix
          cat(paste0("GEX input ", i, " contains multiple count matrices. \n"))
          GEX_ind <- c() #open indices for GEX and FB list elements
          FB_ind <- c()
          for(j in 1:length(directory_read10x[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(directory_read10x[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              cat(paste0("GEX input ", i, " element ", j, " contains > 100 features and will be loaded as GEX \n"))
              GEX_ind <- c(GEX_ind, j)
            } else if(nrow(directory_read10x[[i]][[j]]) < 100){
              cat(paste0("GEX input ", i, " element ", j, " contains < 100 features and will be loaded as FB \n"))
              FB_ind <- c(FB_ind, j)
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }
          if(length(FB_ind) > 0){
            FB.list <- c(FB.list, directory_read10x[[i]][FB_ind]) #add this matrix to the FB list
            directory_read10x[[i]] <- directory_read10x[[i]][-FB_ind] #delete the matrix from the original list, so that only the GEX matrix remains
            #Moreover we want to prevent two user inputs with FB. Here we set FB.loaded to TRUE, so that the FB loading from disk module further down will be skipped.
            FB.loaded <- T
            #at this point we have two lists: the FB.list and the directory_read10x which only contains GEX info, but is still nested. So we have to flatten it
            directory_read10x <- do.call(list, unlist(directory_read10x, recursive=FALSE))
          }
        } else{
          FB.loaded <- F
        }
      }

      #Done => result should be one non-nested list of matrices only containing GEX information and another containing only FB information.

      if(FB.loaded == T){
      #Next we need to check column names / names of feature barcodes. Potentially more than one library of FBs was sequenced giving us more than one matrix but all with the same column names (because same FBs may have been used across samples). We therefore rename the feature barcodes individually to make sure they stay separated
      for(i in 1:length(FB.list)){
        colnames(FB.list[[i]]) <- paste0("i", i, "_", colnames(FB.list[[i]]))
      }
      #We can now convert all to Seurat objects
      FB.list <- lapply(FB.list, function(x) Seurat::CreateSeuratObject(x))
      }
      # Done with FB processing for now.

      #Continuing with GEX processing. Irregardless of FB data presence, the directory_read10x is a flat list of GEX matrices
      directory_read10x <- lapply(directory_read10x, function(x){
        rownames(x) <- toupper(rownames(x))
        return(x)})
      gex.list <- lapply(directory_read10x, function(x) Seurat::CreateSeuratObject(x))
      directory_read10x <- NULL #save some ram

      #load the metrics file conditionally
      GEX.out.directory.list.metrics <- paste(GEX.out.directory.list,"/metrics_summary.csv",sep="")
      if(get.VDJ.stats == T & all(file.exists(GEX.out.directory.list.metrics))){
      gex.metrics.table <- lapply(GEX.out.directory.list.metrics, function(x) utils::read.csv(x,sep=",",header=T, ))
      } else {
        gex.metrics.table <- "none"
      }

      for(i in 1:length(gex.list)){
        #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes => also to integrate FB data later should there be any)
        gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
        gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
      }

      gex.loaded <- T
      if(FB.loaded == F){
        cat("Loaded GEX data \n")
      } else {
        cat("Loaded GEX and FB data \n")
      }
      print(Sys.time())
    })}, error = function(e){
      cat("Loading GEX failed \n")
      print(e)})

  } else if(class(samples.in) == "list"){ #GET info from samples.in
    gex_load_error <- tryCatch({
      #add the directory identifier

      gex.list <- lapply(samples.in, function(x) Seurat::CreateSeuratObject(x[[2]][[1]]))
      if(class(x[[2]][[3]]) != "character"){ #Conditional loading just in case
      FB.list <- lapply(samples.in, function(x) Seurat::CreateSeuratObject(x[[2]][[3]]))
      FB.loaded <- T
      } else {
        FB.loaded <- F
      }
      gex.metrics.table <- lapply(samples.in, function(x) return(x[[2]][[2]]))

      #clear GEX part in object
      for(i in 1:length(samples.in)){
        samples.in[[i]][[2]] <- "loaded"
      }

      for(i in 1:length(gex.list)){
        #Adding column for original barcodes that are not changed upon integration (these are not the colnames, but a metadata column to keep track of original barcodes)
        gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
        gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
      }

      gex.loaded <- T
      cat(paste0(Sys.time()," Loaded GEX data from Data.in \n"))
    }, error = function(e){
      cat("Loading GEX from Data.in failed \n")
      print(e)})
  } else {
    FB.loaded <- F #error catch if no GEX at all was loaded
  }

  ############################################ Load GEX from Seurat object ####

  #check if there is direct seurat input
  if(gex.loaded == F){
  if(class(Seurat.in) == "list"){
    gex.list <- Seurat.in
    Seurat.in <- "loaded"
    for(i in 1:length(gex.list)){
      gex.list[[i]]$orig_barcode <- as.character(gsub(".*_","",colnames(gex.list[[i]])))
      gex.list[[i]]$orig_barcode <- gsub(gex.list[[i]]$orig_barcode,pattern = "-1",replacement = "")
    }
  }
  }

  ############################################ Load Feature barcodes ####

  #! Not checking if FB have already been loaded as part of GEX. If FB directory is supplied, any FB data loaded by GEX is therefore overwritten. This should allow for more flexibility of the user without having to realign the whole GEX data
  if(class(samples.in) == "character" & FB.out.directory.list[[1]] != "none"){ #No Data.in or input => proceed with loading by BC.out.directory.list
    print("HEY")
    #Remove possible backslash at the end of the input path
    for(k in 1:length(FB.out.directory.list)){
      FB.out.directory.list[[k]]<-  gsub("/$", "", FB.out.directory.list[[k]])
    }

    gex_load_error <- tryCatch({suppressWarnings({
      #add the directory identifier

      if(stringr::str_detect(FB.out.directory.list[[1]], "filtered_feature_bc_matrix")){
        FB.out.directory.list.p <- FB.out.directory.list
        cat("\n ! Feature barcode path was specified explicitely to filtered_feature_bc_matrix. For better rates we recommend using the raw_feature_bc_matrix folder content ! \n ")
      } else if (stringr::str_detect(FB.out.directory.list[[1]], "raw_feature_bc_matrix")){
        #Nothing to append
        FB.out.directory.list.p <- FB.out.directory.list
        cat("\n Loading feature barcodes from raw_feature_bc_matrix folder \n")
      } else{
        FB.out.directory.list.p <- paste(FB.out.directory.list,"/raw_feature_bc_matrix",sep="")
        cat("\n Loading feature barcodes from raw_feature_bc_matrix folder \n")
      }
      #Actually loading the data. Critical: if Cellranger 6.1.0 count was run with --libraries input containing both GEX and Feature Barcodes, reading it will result in a list of matrices instead of a single matrix. => the next section deals with this
      directory_read10x <- lapply(FB.out.directory.list.p, function(x) Seurat::Read10X(data.dir=x))

      #dealing with possible mixed GEX FB inputs or multiple FB input matrices from the same directory
      for(i in 1:length(directory_read10x)){ #iterating over main list
        if(class(directory_read10x[[i]]) == "list"){ #this returns true only if the FB directory import contained more than one marices
          cat(paste0("Feature barcode input ", i, " contains multiple count matrices. \n"))
          to_del <- c()
          for(j in 1:length(directory_read10x[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(directory_read10x[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              cat(paste0("Feature barcode input ", i, " element ", j, " contains > 100 features and likely corresponds to GEX data. This matrix will be removed from further FB processing \n"))
              to_del <- c(to_del, j) #Will be deleted later to not mess up the loop
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }
          if(length(to_del) > 0){
            directory_read10x[[i]] <- directory_read10x[[i]][-j] #deleting
          }
        }
      }
      #now flatten the remaining list
      directory_read10x <- do.call(list, unlist(directory_read10x, recursive=FALSE))
      warning(paste0("At least one Feature barcode input contains multiple count matrices. Count matrices containing more than 100 features were filtered out from FB assignment as they most likely correspond to GEX data"))
      #Done => result should be a non-nested list of matrices only containing FB information. With this we can move forward

      #Next we need to check column names / names of feature barcodes. Potentially more than one library of FBs was sequenced giving us more than one matrix but all with the same column names (because same FBs may have been used across samples). We therefore rename the feature barcodes individually to make sure they stay separated
      for(i in 1:length(directory_read10x)){
        rownames(directory_read10x[[i]]) <- paste0("i", i, "_", rownames(directory_read10x[[i]]))
      }

      #We can now convert all to Seurat objects
      FB.list <- lapply(directory_read10x, function(x) Seurat::CreateSeuratObject(x))
      directory_read10x <- NULL

      FB.loaded <- T
      cat("Loaded FB data \n")
      print(Sys.time())
    })}, error = function(e){
      cat("Loading FB failed \n")
      print(e)})

    #TO REVISE, once we know the final format of this r Object
  } else if(class(samples.in) == "list"){ #GET info from samples.in
    #The loading of FBs from data.in input is done during load GEX. In that case FB.loaded was set to T previously, so this part should never be executed
    cat("! Error in FB loading from data.in \n")
  } else {
    FB.loaded <- F
    FB.list <- "none"
  }


  ############################################ Call VDJ stats ####

  #This internal VDJ_GEX_stats function takes already processed input from its wrapper. This shortens loading time, reduces RAM usage and allowes to generate stats when data is provided via Data.in
  stats.done <- F
  if(get.VDJ.stats == T & vdj.loaded == T){
    cat("Getting VDJ GEX stats \n")
    tryCatch({ #error catch (this may be a bit redundant given the error catchers inside the function as well...)
      out.stats <- VDJ_GEX_stats_int(clonotype.list = clonotype.list,
                                     reference.list = reference.list,
                                     annotations.list = annotations.list,
                                     contig.table = contig.table,
                                     vdj.metrics = metrics.table,
                                     gex.metrics = gex.metrics.table,
                                     samples.paths.VDJ = stringr::str_split(samples.paths.VDJ, " ; ", simplify = T)[1,]) #needing to split sample paths again as they will end up in separate rows in the stats dataframe
      stats.done <- T
      cat("Got VDJ GEX stats \n")
      print(Sys.time())

    }, error = function(e){e
      cat("VDJ stats failed: \n")
      print(e)})
  }

  ############################################ Select barcodes for processing ####

  #Proceed with selecting the barcodes which will be processed
  #Any barcode that is a cell in the VDJ or in the GEX it will be processed
  #because the filtered GEX is loaded, all barcodes from the GEX will be included in any case.
  #If a barcode is a cell in VDJ but is not present in the filtered GEX, it will be analyzed for VDJ irregardless

  if(gex.loaded == T & vdj.loaded == T){ #If both VDJ and GEX are available
    barcodes_GEX <- list()
    barcodes_VDJ <- list()
    for(i in 1:length(gex.list)){
      #GEX do not need to be filtered
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])
      #VDJ filtering to make sure => with the new input strategy post Cellranger 6.1. this should not have any effect because of prefiltering by cellranger
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true" & tolower(contig.table[[i]]$high_confidence) == "true" & tolower(contig.table[[i]]$productive) == "true" & tolower(contig.table[[i]]$full_length) == "true")])

      cat(paste0("For sample ", i, ": ", length(barcodes_GEX[[i]])," cell assigned barcodes in GEX, ", length(barcodes_VDJ[[i]]), " cell assigned high confidence barcodes in VDJ. Overlap: ", sum(barcodes_GEX[[i]] %in% barcodes_VDJ[[i]]), " \n"))

      vdj.gex.available <- colnames(gex.list[[i]]) %in% barcodes_VDJ[[i]]
      gex.list[[i]] <- SeuratObject::AddMetaData(gex.list[[i]], vdj.gex.available, col.name = "VDJ.available")

      #remove all barcodes in GEX which are not present in VDJ (defaults to FALSE)
      if(exclude.GEX.not.in.VDJ == T){
        cat("Removing all barcodes from GEX, which are not present in VDJ \n")

        vdj.gex.available <- colnames(gex.list[[i]]) %in% barcodes_VDJ[[i]]
        gex.list[[i]]@meta.data$VDJ.available <- vdj.gex.available
        gex.list[[i]] <- subset(gex.list[[i]], cells = colnames(gex.list[[i]])[which(gex.list[[i]]$VDJ.available == T)])

        cat(paste0("Removed ", length(vdj.gex.available)-sum(vdj.gex.available), " GEX entries \n"))
      }
    }
  }

  #If only GEX is processed
  if(gex.loaded == T & vdj.loaded == F){
    barcodes_GEX <- list()
    for(i in 1:length(gex.list)){
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])
      cat(paste0("For sample ", i, ": ", length(barcodes_GEX[[i]])," cell assigned barcodes in GEX \n"))
    }}
  #If only VDJ is processed
  if(gex.loaded == F & vdj.loaded == T){
    barcodes_VDJ <- list()
    for(i in 1:length(contig.table)){##barcodes_VDJ holds the unique barcodes
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true" & tolower(contig.table[[i]]$high_confidence) == "true" & tolower(contig.table[[i]]$productive) == "true" & tolower(contig.table[[i]]$full_length) == "true")])
      cat(paste0("For sample ", i, ": ", length(barcodes_VDJ[[i]]), " cells assigned with high confidence barcodes in VDJ \n"))
    }
  }

  ############################################ remove sample overlapping barcodes in GEX ####
  if(filter.overlapping.barcodes.GEX == T & gex.loaded == T){
    if(length(gex.list) > 1){
      barcodes_GEX_c <- do.call("c", lapply(gex.list, function(x) x$orig_barcode))
      unique_barcodes <- names(table(barcodes_GEX_c)[table(barcodes_GEX_c) == 1])
      for(i in 1:length(gex.list)){
        gex.list[[i]] <- subset(gex.list[[i]], subset = orig_barcode %in% unique_barcodes)
      }
      cat(paste0("Removed a total of ", length(unique(barcodes_GEX_c)) - length(unique_barcodes), " cells with non unique barcodes in GEX \n"))
    }
  }
  #Sidenote: at this point we do not filter overlapping barcodes in feature barcodes. See further down for this

  ############################################ remove sample overlapping barcodes in VDJ ####
  if(filter.overlapping.barcodes.VDJ == T & vdj.loaded == T){
    if(length(barcodes_VDJ) > 1){
      barcodes_VDJ_c <- do.call("c", barcodes_VDJ)
      non_unique_barcodes <- names(table(barcodes_VDJ_c)[table(barcodes_VDJ_c) > 1])
      for(i in 1:length(barcodes_VDJ)){
        barcodes_VDJ[[i]] <- barcodes_VDJ[[i]][which(!barcodes_VDJ[[i]] %in% non_unique_barcodes)]
      }
      cat(paste0("Removed a total of ", length(non_unique_barcodes), " cells with non unique barcodes in VDJ \n"))
    }
  }

  ############################################ Exclude cells based on marker expression ####
  #handlers copied from GEX_phenotype, Thanks Alex!
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
          cmd[i]<-paste0("barcodes_match_ex_crit <- SeuratObject::WhichCells(gex.list[[j]], slot = 'counts', expression =", cell.state.markers[i],")")
          is.exist<-tryCatch(expr=length(eval(parse(text=cmd[i]))), error = function(x){
            x<-F
            return(x)})
          #if command worked, subset current seurat object
          if(is.exist!=F){
            cells_unfiltered <- ncol(gex.list[[j]])
            gex.list[[j]]$match_ex_crit <- colnames(gex.list[[j]]) %in% barcodes_match_ex_crit
            gex.list[[j]] <- subset(gex.list[[j]], subset = match_ex_crit == T)
            cat(paste0("In GEX sample ", j ," excluded ", cells_unfiltered - ncol(gex.list[[j]])," cells based on ", cell.state.markers[i], " \n"))
            #If the Gene was not found in seurat object features
          } else{
            cat(paste0("In GEX sample ", j ," failed to exclude cells based on: ", cell.state.markers[i], " Please check gene spelling \n"))
          }
        }
      }
      #larger error callback should be independent of user input.
    }, error = function(e){
      cat("Exclusion based on cell markers failed \n")
      print(e)
    })
  }

  ############################################ Subsample VDJ barcodes ####
  if(subsample.barcodes == T & vdj.loaded == T){
    #For development => shorten barcode list to shorten computational time during development
    cat("Sampling 50 barcodes from all in VDJ per sample \n")
    for(i in 1:length(barcodes_VDJ)){
      barcodes_VDJ[[i]] <- sample(barcodes_VDJ[[i]],50)
    }
  }

  ############################################ VDJ Processing per cell ####
  if(vdj.loaded == T){
    VDJ.proc.list <- list()
    for(i in 1:length(contig.table)){

      cat(paste0("Starting VDJ barcode iteration ", i , " of ", length(contig.table), "...\n"))
      print(Sys.time())
      if(parallel.processing == "parlapply" | parallel.processing == "parLapply"){
        #close any open clusters
        doParallel::stopImplicitCluster()
        #open cluster for parallel computing

        cl <- parallel::makeCluster(numcores)
        cat(paste0("Started parlapply cluster with ", numcores, " cores \n"))

        out.VDJ <- parallel::parLapply(cl, barcodes_VDJ[[i]], barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align, select.excess.chains.by.umi.count = select.excess.chains.by.umi.count)
        #close any open clusters
        doParallel::stopImplicitCluster()
        gc() #garbage collection to reduce ram impact

      } else if(parallel.processing == "mclapply"){
        cat(paste0("Started mcapply cluster with ", numcores, " cores \n"))
        out.VDJ <- parallel::mclapply(X = barcodes_VDJ[[i]], FUN = barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align, select.excess.chains.by.umi.count = select.excess.chains.by.umi.count)
      } else { #No parallel computing
        out.VDJ <- lapply(barcodes_VDJ[[i]], barcode_VDJ_iteration, contigs = contig.table[[i]], references = reference.list[[i]], annotations = annotations.table[[i]], gap.extension.cost = gap.extension.cost, gap.opening.cost = gap.opening.cost, trim.and.align = trim.and.align, select.excess.chains.by.umi.count = select.excess.chains.by.umi.count)
        gc() #garbage collection to reduce ram impact
      }

      #bind list recieved from parLapply
      VDJ.proc <- dplyr::bind_rows(out.VDJ)
      VDJ.proc[VDJ.proc == ";"] <- "" #fix bug, where if two emtpy strings are concatenated, a ";" is left behind.

      #update barcodes
      VDJ.proc$orig_barcode <- VDJ.proc$barcode
      VDJ.proc$orig_barcode <- gsub("-\\d", "", VDJ.proc$orig_barcode)
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

      cat(paste0("Done with ", i , " of ", length(contig.table), " \n"))
      print(Sys.time())
    }

    VDJ.proc <- VDJ.proc.list

    #reduce ram impact
    contig.table <- NULL
    reference.list <- NULL
    annotations.table <- NULL

    if (VDJ.combine == T){ #processing all VDJ files together ! THIS RETURNS A DATAFRAME
      VDJ.proc <- dplyr::bind_rows(VDJ.proc)
    }
  } else{
    VDJ.proc <- "none"
  }

  ############################################ GEX Processing ####
  if(gex.loaded == T & Seurat.in != "loaded"){ #make sure that Seurat.in is not an input list!
    cat("Starting GEX pipeline \n")
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
      GEX.proc[[i]]@meta.data$orig_barcode <- as.character(gsub(".*_","",rownames(GEX.proc[[i]]@meta.data)))
      GEX.proc[[i]]@meta.data$orig_barcode <- gsub(GEX.proc[[i]]@meta.data$orig_barcode ,pattern = "-1",replacement = "")
    }
  } else {
    GEX.proc <- "none"
  }
  cat("Done with GEX pipeline \n")
  print(Sys.time())


  ############################################ Feature Barcode assignment to GEX ####
  #Here we merge feature barcodes into the GEX matrices.
  #Because users may choose to not combine individual samples, we need to deal with multiple seurat objects in the GEX.proc and VDJ.proc lists

  if(FB.loaded == T & class(GEX.proc) == "list"){ #ensuring correct input

  cat("Adding Feature barcode information to GEX... \n")

    tryCatch({
    #getting relevant info as a dataframe
    FB.list <- lapply(FB.list, function(x) return(SeuratObject::FetchData(x, rownames(x))))
    #assigning FBs using also the user tunable threshold as an additional argument
    FB.list <- lapply(FB.list, pick_max_feature_barcode, FB.ratio.threshold)
    #This should return a list of dataframes with each 2 columns: "orig_barcode" and "FB_assignment"
    FB.processed <- T #keeping track for VDJ
    print("done processing")

    if(length(GEX.proc) == 1){ #Either only one sample or GEX.integrate == T

      meta_to_merge <- SeuratObject::FetchData(GEX.proc[[1]], "orig_barcode") #getting a reference to merge into
      for(i in 1:length(FB.list)){ #iterating over FB tables
        meta_to_merge <- merge(meta_to_merge, FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F, sort = F) #merging sequentially making sure to not add or remove any rows
      }
      rownames(meta_to_merge) <- colnames(GEX.proc[[1]]) #reconstitute the rownames. Merge deletes those. The sort= F in merge is neccessary to not mess up the order
      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], meta_to_merge[,"FB_assignment"], col.name = "FB_assignment") #add to object. We only add the FB_assignment column. otherwise Seurat throws an error
      #move the column to after sample_id as this column will probably be used a lot
      sample_id_index <- which(names(GEX.proc[[1]]@meta.data) == "sample_id")
      GEX.proc[[1]]@meta.data <- GEX.proc[[1]]@meta.data[,c(1:sample_id_index, ncol(GEX.proc[[1]]@meta.data), (sample_id_index+1):(ncol(GEX.proc[[1]]@meta.data)-1))]

    } else if(length(GEX.proc) > 1){
      if(length(GEX.proc) == length(FB.list)){ #if there is 1-1 correspondence

        for(i in 1:length(FB.list)){ #iterating over FB tables
          meta_to_merge <- SeuratObject::FetchData(GEX.proc[[i]], "orig_barcode") #getting a reference to merge into
          meta_to_merge <- merge(meta_to_merge, FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F, sort = F) #merging sequentially making sure to not add or remove any rows
          rownames(meta_to_merge) <- rownames(GEX.proc[[i]]@meta.data)#reconstitute the rownames. Merge deletes those. The sort= F in merge is neccessary to not mess up the order
          GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], meta_to_merge[,"FB_assignment"], col.name = "FB_assignment") #add to object
          #move the column to after sample_id as this column will probably be used a lot
          sample_id_index <- which(names(GEX.proc[[i]]@meta.data) == "sample_id")
          GEX.proc[[i]]@meta.data <- GEX.proc[[i]]@meta.data[,c(1:sample_id_index, ncol(GEX.proc[[i]]@meta.data), (sample_id_index+1):(ncol(GEX.proc[[i]]@meta.data)-1))]
        }

      } else { #no one to one correspondence (Note: this may be the case if for e.g. 3 GEX samples the FB barcodes were pooled for sequencing or merged in post)
        #here we merge all FB matrices together and then merge that matrix sequentially to all GEX matrices
        if(length(FB.list) > 1){
            FB.list <- lapply(FB.list, function(x){ #remove rownames that would mess up merging if there are any duplicate barcodes
              rownames(x) <- NULL
              return(x)
            })
            FB.list <- dplyr::bind_rows(FB.list) #returns a dataframe
        } else {
            FB.list <- FB.list[[1]] #keep the format consistent and convert to dataframe
          }

        for(i in 1:length(GEX.proc)){ #iterating over FB tables
          meta_to_merge <- SeuratObject::FetchData(GEX.proc[[i]], "orig_barcode") #getting a reference to merge into
          meta_to_merge <- merge(meta_to_merge, FB.list, by = "orig_barcode", all.x = T, all.y = F, sort = F) #merging sequentially making sure to not add or remove any rows
          rownames(meta_to_merge) <- rownames(GEX.proc[[i]]@meta.data)#reconstitute the rownames. Merge deletes those. The sort= F in merge is neccessary to not mess up the order
          GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], meta_to_merge[,"FB_assignment"]) #add to object
          #move the column to after sample_id as this column will probably be used a lot
          sample_id_index <- which(names(GEX.proc[[i]]@meta.data) == "sample_id")
          GEX.proc[[i]]@meta.data <- GEX.proc[[i]]@meta.data[,c(1:sample_id_index, ncol(GEX.proc[[i]]@meta.data), (sample_id_index+1):(ncol(GEX.proc[[i]]@meta.data)-1))]
        }
      }
    }
    meta_to_merge <- NULL
    }, error = function(e){
      cat("Adding Feature barcode information to GEX failed \n")
      print(e)
    })
  } else {
    FB.processed <- F #keeping track for VDJ
  }


  ############################################ Feature Barcode assignment to VDJ ####
  #Here we merge feature barcodes into the VDJ matrices.
  #Because users may choose to not combine individual samples, we need to deal with multiple seurat objects in the GEX.proc and VDJ.proc lists

  if(FB.loaded == T & (class(VDJ.proc) == "list" | class(VDJ.proc) == "data.frame")){ #ensuring correct input

    cat("Adding Feature barcode information to VDJ... \n")

    tryCatch({
    if(FB.processed == F){ #process only if not done before in GEX
    #getting relevant info as a dataframe
    FB.list <- lapply(FB.list, function(x) return(SeuratObject::FetchData(x, vars = rownames(x))))
    #assigning FBs using also the user tunable threshold as an additional argument
    FB.list <- lapply(FB.list, pick_max_feature_barcode, FB.ratio.threshold)
    #This should return a list of dataframes with each 2 columns: "orig_barcode" and "FB_assignment"
    }

    if(class(VDJ.proc) == "data.frame"){ #Either only one sample or VDJ.combine == T

      for(i in 1:length(FB.list)){ #iterating over FB tables
        VDJ.proc <- merge(VDJ.proc, FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F) #merging sequentially making sure to not add or remove any rows
      }
      #move the column to after sample_id as this column will probably be used a lot
      sample_id_index <- which(names(VDJ.proc) == "sample_id")
      VDJ.proc<- VDJ.proc[,c(1:sample_id_index, ncol(VDJ.proc), (sample_id_index+1):(ncol(VDJ.proc)-1))]

    } else if(class(VDJ.proc) == "list"){
      if(length(VDJ.proc) == length(FB.list)){ #if there is 1-1 correspondence

        for(i in 1:length(FB.list)){ #iterating over FB tables
          VDJ.proc[[i]] <- merge(VDJ.proc[[i]], FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F)
          #move the column to after sample_id as this column will probably be used a lot
          sample_id_index <- which(names(VDJ.proc[[i]]) == "sample_id")
          VDJ.proc[[1]]<- VDJ.proc[[1]][,c(1:sample_id_index, ncol(VDJ.proc[[1]]), (sample_id_index+1):(ncol(VDJ.proc[[1]])-1))]
        }

      } else { #no one to one correspondence (Note: this may be the case if for e.g. 3 GEX samples the FB barcodes were pooled for sequencing or merged in post)
        #here we merge all FB matrices together and then merge that matrix sequentially to all GEX matrices
        if(length(FB.list) > 1){
          FB.list.merge <- lapply(FB.list, function(x){ #remove rownames that would mess up merging if there are any duplicate barcodes
            rownames(x) <- NULL
            return(x)
          })
          FB.list.merge <- dplyr::bind_rows(FB.list.merge) #returns a dataframe
        } else {
          FB.list.merge <- FB.list[[1]] #keep the format consistent and convert to dataframe
        }

        for(i in 1:length(VDJ.proc)){ #iterating over FB tables
          VDJ.proc[[i]] <- merge(VDJ.proc[[i]], FB.list.merge, by = "orig_barcode", all.x = T, all.y = F)
          #move the column to after sample_id as this column will probably be used a lot
          sample_id_index <- which(names(VDJ.proc[[i]]) == "sample_id")
          VDJ.proc[[i]]<- VDJ.proc[[i]][,c(1:sample_id_index, ncol(VDJ.proc[[i]]), (sample_id_index+1):(ncol(VDJ.proc[[i]])-1))]
        }
      }
    }
      FB.list <- NULL
      FB.list.merge <- NULL
    }, error = function(e){
      cat("Adding Feature barcode information to VDJ failed \n")
      print(e)
    })
  }

  ############################################ VDJ GEX integration ####

  if(vdj.loaded == T & gex.loaded == T){
    cat("Integrating VDJ and GEX...\n")
    #Combine GEX and VDJ
    #Because a lot of cells are present in only one of the two datasets, we do not add the VDJ table as a metadata item to the Seurat object. This would throw out a lot high confidence VDJ entries. Instead we add metadata columns to both datasets indicating the presence of the cell in the opposing dataset. This should make filtering and joining later easy.
    #Multiple cases:
    if(class(VDJ.proc) == "list" & length(GEX.proc) == 1){ #Multiple VDJ, one GEX
      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], colnames(GEX.proc[[1]]) %in% do.call(c ,lapply(VDJ.proc, function(x) x$barcode)), col.name = "VDJ_available")
      for(i in 1:length(VDJ.proc)){
        VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[1]])
      }

      #add some GEX columns to VDJ table
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

      #add VDJ info to GEX as metadata columns
      if(integrate.VDJ.to.GEX == T){

        seur_meta <- GEX.proc[[1]]@meta.data
        #set barcodes as columns for merging
        seur_meta$barcode <- rownames(seur_meta)
        #merge
        drop <- c("sample_id", "group_id", "orig.ident", "seurat_clusters", "orig_barcode.x", "orig.barcode.y","orig_barcode", "GEX_available", "FB_assignment")

        #rbind VDJ proc objects and add them together
        VDJ.proc.all <- dplyr::bind_rows(VDJ.proc)

        seur_meta <- merge(seur_meta, VDJ.proc.all[,which(!names(VDJ.proc.all) %in% drop)], by = "barcode", all.x = T) #these columns are excluded here because they already exist in the seurat object

        VDJ.proc.all <- NULL

        #add rownames to new dataframe, otherwise AddMetaData fails
        rownames(seur_meta) <- seur_meta$barcode
        #add metadata
        GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], seur_meta[,c(10:ncol(seur_meta))], col.name = names(seur_meta)[c(10:ncol(seur_meta))])
      }

      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]

    } else if (class(VDJ.proc) == "data.frame" & length(GEX.proc) == 1){ #one VDJ one GEX
      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], colnames(GEX.proc[[1]]) %in% VDJ.proc$barcode, col.name = "VDJ_available")
      VDJ.proc$GEX_available <- VDJ.proc$barcode %in% colnames(GEX.proc[[1]])
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]

     #add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){

        #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
        seur_meta <- SeuratObject::FetchData(GEX.proc,
                               vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
        names(seur_meta)[2] <- "orig_barcode_GEX"
        seur_meta$barcode <- rownames(seur_meta)
        #merge to VDJ.proc
        VDJ.proc <- merge(VDJ.proc, seur_meta, by = "barcode", all.x = T, all.y = F)
      }

      #add VDJ info to GEX as metadata columns
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

    } else if(class(VDJ.proc) == "data.frame" & length(GEX.proc) > 1){ #one VDJ multiple GEX (improbable...)
      for(i in 1:length(GEX.proc)){
        GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], colnames(GEX.proc[[i]]) %in% VDJ.proc$barcode, col.name = "VDJ_available")
      }
      VDJ.proc$GEX_available <- VDJ.proc$barcode %in% do.call(c ,lapply(GEX.proc, function(x) colnames(x)))

      #add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){
        cat("Adding data from multiple GEX objects to one VDJ object \n")
        #grab metadata from all seurat objects
        seur_meta <- lapply(GEX.proc, function(x){SeuratObject::FetchData(x,
         vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))})
        seur_meta <- dplyr::bind_rows(seur_meta)
        names(seur_meta)[2] <- "orig_barcode_GEX"
        seur_meta$barcode <- rownames(seur_meta)
        #merge to VDJ.proc
        VDJ.proc <- merge(VDJ.proc, seur_meta, by = "barcode", all.x = T, all.y = F)
      }

      #add VDJ info to GEX as metadata columns
      if(integrate.VDJ.to.GEX == T){
      cat("Integrating VDJ from a single object to multiple GEX objects is not supported \n")
      }

    } else if(class(VDJ.proc) == "list" & length(GEX.proc) > 1){ #multiple VDJ multiple GEX
      for(i in 1:length(GEX.proc)){
        GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], colnames(GEX.proc[[i]]) %in% VDJ.proc[[i]]$barcode, col.name = "VDJ_available")
        VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[i]])

        #add some GEX columns to VDJ table
        if(integrate.GEX.to.VDJ == T){
          cat("Integrating multiple VDJ and GEX in a pariwise way. ! VDJ and GEX input paths or data must be in the same order ! \n")
          #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
          seur_meta <- SeuratObject::FetchData(GEX.proc[[i]],
                                               vars = c("orig.ident","orig_barcode","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
          names(seur_meta)[2] <- "orig_barcode_GEX"
          seur_meta$barcode <- rownames(seur_meta)
          #merge to VDJ.proc
          VDJ.proc[[i]] <- merge(VDJ.proc[[i]], seur_meta, by = "barcode", all.x = T, all.y = F)
        }

        #add VDJ info to GEX as metadata columns
        if(integrate.VDJ.to.GEX == T){
          cat("Integrating multiple VDJ and GEX in a pariwise way. ! VDJ and GEX input paths or data must be in the same order ! \n")
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
      }
    }

    out.list <- list("VDJ" = VDJ.proc, "GEX" = GEX.proc)

  } else if(vdj.loaded == T & gex.loaded == F){
    out.list <- list("VDJ" = VDJ.proc, "GEX" = "none")
  } else if(vdj.loaded == F & gex.loaded == T){
    if(length(GEX.proc) == 1){
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]
    }
    out.list <- list("VDJ" = "none", "GEX" = GEX.proc)
  } else if(vdj.loaded == F & gex.loaded == F){
    stop("Neither VDJ or GEX data loaded. Exiting")
  }

  ############################################ Compile output ####


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
  cat("Adding VDJ stats...\n")
  if(stats.done == T){
    out.list[[3]] <- out.stats
  } else {
    out.list[[3]] <- "VDJ stats failed"
  }

  #adding parameter info for reproducibility
  cat("Adding runtime params...\n")
  out.list[[4]] <- params

  #finally add session info
  out.list[[5]] <- utils::sessionInfo()

  #rename for clarity
  names(out.list) <- c("VDJ", "GEX", "VDJ.GEX.stats", "Running params", "sessionInfo")

  cat("Done!")
  print(Sys.time())
  # if(class(out.list[[2]])=="Seurat"){
  #   out.list[[2]]$sample_id <- as.integer(gsub(pattern = "s",replacement = "",sub("_.*", "", names(out.list[[2]]$orig_barcode))))
  #   out.list[[2]]$group_id <- rep(group.id,table(out.list[[2]]$sample_id))
  # }

  if(class(out.list[[1]])=="data.frame") out.list[[1]]$clonotype_id <- out.list[[1]]$clonotype_id_10x
  return(out.list)
}
