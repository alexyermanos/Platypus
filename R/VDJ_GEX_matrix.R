#'VDJ GEX processing and integration wrapper
#'
#' @description
#' This function is designed as a common input to the Platypus pipeline. Integration of datasets as well as VDJ and GEX information is done here. Please check the Platypus V3 vignette for a detailed walkthrough of the output structure. In short: output[[1]] = VDJ table, output[[2]] = GEX Seurat object and output[[3]] = statistics
#'[FB] Feature barcode (FB) technology is getting increasingly popular, which is why Platypus V3 fully supports their use as sample delimiters. As of V3, Platpyus does not support Cite-seq data natively, also the VDJ_GEX_matrix function is technically capable of loading a Cite-seq matrix and integrating it with VDJ. For details on how to process sequencing data with FB data and how to supply this information to the VDJ_GEX_matrix function, please consult the dedicated vignette on FB data.
#' @param VDJ.out.directory.list List containing paths to VDJ output directories from cell ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires. ! Neccessary files within this folder: filtered_contig_annotations.csv, clonotypes.csv, concat_ref.fasta, all_contig_annotations.csv (only if trim.and.align == T) and metrics_summary.csv (Optional, will be appended to stats table if get.VDJ.stats == T)
#'@param GEX.out.directory.list List containing paths the outs/ directory of each sample or directly the raw or filtered_feature_bc_matrix folder. Order of list items must be the same as for VDJ.
#'@param FB.out.directory.list [FB] List of paths pointing at the outs/ directory of output from the Cellranger counts function which contain Feature barcode counts. ! Single list elements can be a path or "PLACEHOLDER", if the corresponding input in the VDJ or GEX path does not have any adjunct FB data. This is only the case when integrating two datasets of which only one has FB data. See examples for details. Any input will overwrite potential FB data loaded from the GEX input directories. This may be important, if wanting to input unfiltered FB data that will cover also cells in VDJ not present in GEX.
#'@param Data.in Input for R objects from either the PlatypusDB_load_from_disk or the PlatypusDB_fetch function. If provided, input directories should not be specified. If you wish to integrate local and downloaded data, please load them via load_from_disk and fetch and provide as a list (e.g. Data.in = list(load_from_disk.output, fetch.output))
#'@param Seurat.in Alternative to GEX.out.directory.list. A seurat object. VDJ.integrate has to be set to TRUE. In metadata the column of the seurat object, sample_id and group_id must be present. sample_id must contain ids in the format "s1", "s2" ... "sn" and must be matching the order of VDJ.out.directory.list. No processing (i.e. data normalisation and integration) will be performed on these objects. They will be returned as part of the VGM and with additional VDJ data if integrate.VDJ.to.GEX = T. Filtering parameters such as overlapping barcodes, exclude.GEX.not.in.VDJ and exclude.on.cell.state.markers will be applied to the Seurat.in GEX object(s).
#'@param group.id vector with integers specifying the group membership. c(1,1,2,2) would specify the first two elements of the input VDJ/GEX lists are in group 1 and the third/fourth input elements will be in group 2.
#' @param GEX.read.h5 Boolean. defaults to FALSE. Whether to read GEX data from an H5 file. If set to true, please provide the each directory containing a cellranger H5 output file or a direct path to a filtered_feature_bc_matrix.h5 as one GEX.out.directory.list element.
#'@param VDJ.combine Boolean. Defaults to TRUE. Whether to integrate repertoires. A sample identifier will be appended to each barcode both in GEX as well as in VDJ. Recommended for all later functions
#'@param GEX.integrate Boolean. Defaults to TRUE. Whether to integrate GEX data. Default settings use the seurat scale.data option to integrate datasets. Sample identifiers will be appended to each barcode both in GEX and VDJ This is helpful when analysing different samples from the same organ or tissue, while it may be problematic when analysing different tissues.
#'@param integrate.GEX.to.VDJ Boolean. defaults to TRUE. Whether to integrate GEX metadata (not raw counts) into the VDJ output dataframe ! Only possible, if GEX.integrate and VDJ.combine are either both FALSE or both TRUE
#'@param integrate.VDJ.to.GEX Boolean. defaults to TRUE. Whether to integrate VDJ data into GEX seurat object as metadata. ! Only possible, if GEX.integrate and VDJ.combine are either both FALSE or both TRUE
#'@param exclude.GEX.not.in.VDJ Boolean. defaults to FALSE. Whether to delete all GEX cell entries, for which no VDJ information is available. Dependent on data quality and sequencing depth this may reduce the GEX cell count by a significant number
#'@param filter.overlapping.barcodes.GEX Boolean. defaults to TRUE. Whether to remove barcodes which are shared among samples in the GEX analysis. Shared barcodes normally appear at a very low rate.
#'@param filter.overlapping.barcodes.VDJ Boolean. defaults to TRUE. Whether to remove barcodes which are shared among samples in the GEX analysis. Shared barcodes normally appear at a very low rate.
#'@param get.VDJ.stats Boolean. defaults to TRUE. Whether to generate general statistics table for VDJ repertoires. This is appended as element [[3]] of the output list.
#'@param append.raw.reference Boolean. Defaults to TRUE. This appends the raw reference sequence for each contig even if trim.and.align is set to FALSE.
#'@param select.excess.chains.by.umi.count Boolean. Defaults to FALSE. There are several methods of dealing with cells containing reads for more than 1VDJ and 1VJ chain. While many analyses just exclude such cells, the VGM is designed to keep these for downstream evaluation (e.g. in VDJ_clonotype). This option presents an evidenced-based way of selectively keeping or filtering only one of the present VDJ  and VJ chains each. This works in conjunction with the parameter excess.chain.confidence.count.threshold (below) Idea source: Zhang W et al. Sci Adv. 2021 (10.1126/sciadv.abf5835)
#'@param excess.chain.confidence.count.threshold Interger. Defaults to 1000. This sets a umi count threshold for keeping excessive chains in a cell (e.g. T cells with 2 VJ and 1 VDJ chain) and only has an effect if select.excess.chains.by.umi.count is set to TRUE. For a given cell with chains and their UMI counts: VDJ1 = 3, VDJ2 = 7, VJ1 = 6. If count.threshold is kept at default (1000), the VDJ chain with the most UMIs will be kept (VDJ2), while the other is filtered out (VDJ1), leaving the cell as VDJ2, VJ1. If the count.threshold is set to 3, both chains VDJ chains of this cell are kept as their UMI counts are equal or greater to the count.threshold and therefore deemed high confidence chains. In the case of UMI counts being equal for two chains AND below the count.threshold, the first contig entry is kept, while the second is filtered. To avoid filtering excess chains, set select.excess.chains.by.umi.count to FALSE. For further notes on the implication of these please refer to the documentation of the parameter hierarchical in the function VDJ_clonotype_v3.
#'@param trim.and.align Boolean. Defaults to FALSE. Whether to trim VJ/VDJ seqs, align them to the 10x reference and trim the reference. This is useful to get full sequences for antibody expression or numbers of somatic hypermutations. !Setting this to TRUE significantly increases computational time
#'@param parallel.processing Character string. Can be "parlapply" for Windows system, "mclapply" for unix and Mac systems or "none" to use a simple for loop (slow!). Default is "none" for compatibility reasons. For the parlapply option the packages parallel, doParallel and the dependency foreach are required
#'@param numcores Number of cores used for parallel processing. Defaults to number of cores available. If you want to chek how many cores are available use the library Parallel and its command detectCores() (Not setting a limit here when running this function on a cluster may cause a crash)
#'@param gap.opening.cost Argument passed to Biostrings::pairwiseAlignment during alignment to reference. Defaults to 10
#'@param gap.extension.cost Argument passed to Biostrings::pairwiseAlignment during alignment to reference. Defaults to 4
#'@param exclude.on.cell.state.markers Character vector. If no input is provided or input is "none", no cells are excluded. Input format should follow: Character vector containing the gene names for each state. ; is used to use multiple markers within a single gene state. Different vector elements correspond to different states. Example: c("CD4+;CD44-","CD4+;IL7R+;CD44+"). All cells which match any of the given states (in the example case any of the 2) are excluded. This is useful in case different and non lymphocyte cells were co-sequenced. It should give the option to e.g. exclude B cells in the analysis of T cells in a dataset.
#'@param exclude.on.barcodes Character vector. Provide a list of 10x barcodes WITHOUT the terminal id (-1 , -2 etc.) to exclude from GEX and VDJ prior to processing.
#' @param integration.method String specifying which data normalization and integration pipeline should be used. Default is "scale.data", which correspondings to the ScaleData function internal to harmony package. 'anchors' scales data individually and then finds and align cells in similar states as described here: https://satijalab.org/seurat/articles/integration_introduction.html. 'sct'specifies SCTransform from the Seurat package. "harmony" should be specificied to perform harmony integration. This method requires the harmony package from bioconductor.
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
#' @param FB.ratio.threshold Numeric. Defaults to 2 Threshold for assignment of feature barcodes by counts. A feature barcode is assigned to a cell if its counts are >FB.count.threshold and if its counts are FB.ratio.threshold-times higher than the counts of the feature barcode with second most counts.
#' @param FB.count.threshold Numeric. Defaults to 10. For description of Feature Barcode assignment see parameter FB.ratio.threshold above
#' @param FB.exclude.pattern Character (regex compatible). If a feature barcode matches this pattern it will be excluded from the hashing sample assignments. This may be neccessary if CITE-seq barcodes and hashing barcodes are sequenced in the same run.
#'@param verbose if TRUE prints runtime info to console. Defaults to TRUE
#' @return Single cell matrix including VDJ and GEX info. Format is a list with out[[1]] = a VDJ dataframe (or list of dataframes if VDJ.combine == F, not recommended) containing also selected GEX information of integrate.GEX.to.VDJ = T. out[[2]] = GEX Seurat object with the metadata also containing GEX information if integrate.VDJ.to.GEX = T. out[[3]] = Dataframe with statistics on GEX and VDJ. out[[4]] = runtime parameters. out[[5]] = session info
#' @export
#' @examples
#' \dontrun{
#'
#' #FOR EXAMPLES see Platypus vignette at https://alexyermanos.github.io/Platypus/index.html
#'
#'#Run from local directory input. For run from PlatypusDB input see
#'#PlatypusDB vignette
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
#'
#' # With Feature Barcodes
#' ## Option 1: Cellranger multi or Cellranger count with --libraries output
#' VDJ.out.directory.list <- list()
#' VDJ.out.directory.list[[1]] <- "~/VDJ/S1/" #point to outs or per_sample_outs directory content
#' VDJ.out.directory.list[[2]] <- "~/VDJ/S2/"
#' GEX.out.directory.list <- list()
#' GEX.out.directory.list[[1]] <- "~/GEX/S1/"
#' GEX.out.directory.list[[2]] <- "~/GEX/S2/" #These directories contain two matrices (GEX and FB)
#' VGM <- VDJ_GEX_matrix(
#' VDJ.out.directory.list = VDJ.out.directory.list
#' ,GEX.out.directory.list = GEX.out.directory.list,
#' FB.ratio.threshold = 2)
#'
#' ##Option 2: Separate input of FB data from separate Cellranger count run
#' VDJ.out.directory.list <- list()
#' VDJ.out.directory.list[[1]] <- "~/VDJ/S1/"
#' VDJ.out.directory.list[[2]] <- "~/VDJ/S2/"
#' GEX.out.directory.list <- list()
#' GEX.out.directory.list[[1]] <- "~/GEX/S1/"
#' GEX.out.directory.list[[2]] <- "~/GEX/S2/"
#' GEX.out.directory.list <- list()
#' FB.out.directory.list[[1]] <- "~FB/S1/"
#' FB.out.directory.list[[2]] <- "~FB/S1/"
#' VGM <- VDJ_GEX_matrix(
#' VDJ.out.directory.list = VDJ.out.directory.list,
#' GEX.out.directory.list = GEX.out.directory.list,
#' FB.out.directory.list = FB.out.directory.list,
#' FB.ratio.threshold = 2)
#'
#' ##Option 3: FB input for two datasets of which only one contains FB data
#' VDJ.out.directory.list <- list()
#' VDJ.out.directory.list[[1]] <- "~/study1/VDJ/S1/"
#' VDJ.out.directory.list[[2]] <- "~/study2/VDJ/S1/"
#' VDJ.out.directory.list[[3]] <- "~/study2/VDJ/S2/"
#' GEX.out.directory.list <- list()
#' GEX.out.directory.list[[1]] <- "~/study1/GEX/S1/"
#' GEX.out.directory.list[[2]] <- "~/study2/GEX/S1/"
#' GEX.out.directory.list[[2]] <- "~/study2/GEX/S2/"
#' GEX.out.directory.list <- list()
#' FB.out.directory.list[[1]] <- "PLACEHOLDER" #Study 1 does not contain FB data
#' FB.out.directory.list[[2]] <- "~/study2/FB/S1/"
#' FB.out.directory.list[[3]] <- "~/study2/FB/S2/"
#' VGM <- VDJ_GEX_matrix(
#' VDJ.out.directory.list = VDJ.out.directory.list,
#' GEX.out.directory.list = GEX.out.directory.list,
#' FB.out.directory.list = FB.out.directory.list,
#' FB.ratio.threshold = 2)
#'
#' }
VDJ_GEX_matrix <- function(VDJ.out.directory.list,
                           GEX.out.directory.list,
                           FB.out.directory.list,
                           Data.in,
                           Seurat.in,
                           group.id,
                           GEX.read.h5,
                           VDJ.combine,
                           GEX.integrate,
                           integrate.GEX.to.VDJ,
                           integrate.VDJ.to.GEX,
                           exclude.GEX.not.in.VDJ,
                           filter.overlapping.barcodes.GEX,
                           filter.overlapping.barcodes.VDJ,
                           get.VDJ.stats,
                           append.raw.reference,
                           select.excess.chains.by.umi.count,
                           excess.chain.confidence.count.threshold,
                           trim.and.align,
                           parallel.processing,
                           numcores,
                           gap.opening.cost,
                           gap.extension.cost,
                           exclude.on.cell.state.markers,
                           exclude.on.barcodes,
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
                           FB.count.threshold,
                           FB.ratio.threshold,
                           FB.exclude.pattern,
                           subsample.barcodes,
                           verbose){

  #### Variable setup ####

  orig_barcode <- NULL
  match_ex_crit <- NULL
  barcode_VDJ_iteration <- NULL
  GEX_automate_single <- NULL
  VDJ_GEX_stats <- NULL
  do <- NULL
  specificity <- NULL
  affinity <- NULL
  VDJ_raw_ref <- NULL
  VJ_raw_ref <- NULL
  gex.metrics.table <- "none" #error avoidance in case only VDJ is provided and stats are requested
  #Concerning input types
  gex.loaded <- F
  fb.loaded <- F
  vdj.loaded <- F
  seurat.loaded <- F
  data.in.loaded <- F
  #### Function def: pick_max_feature_barcode ####

  #Assignment of Feature Barcodes to cells (returns a dataframe: rownames= barcodes, column 1 = barcodes, column 2 = feature barcode assignments.) Works identically for VDJ and GEX
  pick_max_feature_barcode <- function(bc_df, #input is a dataframe with the columns beeing the count for each feature barcode and the rows being cells. ! No columns allowed apart from numeric count columns
                                       FB.ratio.threshold,
                                       FB.count.threshold){

    #convert intput df to list
    bc_out <- as.list(as.data.frame(t(bc_df)))

    #Choose max barcodes and return the index of the valid maximum BC or 1000 (see below why)
    bc_out <- sapply(bc_out, function(x,FB.ratio.threshold, FB.count.threshold){
      y <- sort(x, decreasing = T)
      if(max(y) > FB.count.threshold){
        if((y[1] / y[2]) > FB.ratio.threshold){
          return(which(x == y[1]))
        } else{
          return(1000)
        }
      } else {
        return(1000)
      }
    },FB.ratio.threshold, FB.count.threshold)

    #from the returned indices pick the right Barcode name from the names of the original dataframe
    bc_out <- sapply(bc_out, function(x, bc_names){
      return(bc_names[x])
    }, names(bc_df))
    #because of the return(1000), cells which could not be assigned will have an NA in the final vector. This is replaced here for clarity
    bc_out[is.na(bc_out)] <- "Not assignable"

    #Open dataframe for results
    bc_match <- data.frame("orig_barcode" = rownames(bc_df), "FB_assignment" = bc_out)
    bc_match[,1] <-  gsub("(^_)|(-\\d+.*$)","",bc_match[,1]) #remove the -1 at the end of the barcode to allow for merging into GEX and VDJ later
    rownames(bc_match) <- rownames(bc_df)
    return(bc_match)
  }

  #### Function def: VDJ_GEX_stats_int ####

  #Gets statistics on VDJ and GEX
  VDJ_GEX_stats_int <- function(clonotype.list,
                                reference.list,
                                annotations.list,
                                contig.table,
                                vdj.metrics,
                                gex.metrics,
                                samples.paths.VDJ,
                                verbose){

    #Getting info from clonotype.list csv
    VDJ.stats.list <- list()
    for(k in 1:length(clonotype.list)){

      if(verbose) message(paste0("Starting with ", k, " of ", length(clonotype.list)))
      VDJ.stats <- c()
      #gsub to unify formatting. B cell notation is picked arbitratily
      contig.table[[k]]$chain <- gsub("(TRB)|(TRG)","IGH", contig.table[[k]]$chain)
      contig.table[[k]]$chain <- gsub("(TRA)|(TRD)","IGL", contig.table[[k]]$chain)
      clonotype.list[[k]]$cdr3s_aa <- gsub("(TRB:)|(TRG:)","IGH:", clonotype.list[[k]]$cdr3s_aa)
      clonotype.list[[k]]$cdr3s_aa <- gsub("(TRA:)|(TRD:)","IGL:", clonotype.list[[k]]$cdr3s_aa)

      #info on sample
      VDJ.stats[length(VDJ.stats)+1] <- samples.paths.VDJ[k]
      names(VDJ.stats)[length(VDJ.stats)] <- "Repertoir path"

      VDJ.stats[length(VDJ.stats)+1] <- paste0("s",k)
      names(VDJ.stats)[length(VDJ.stats)] <- "Sample name"

      #Get number of unique barcodes
      VDJ.stats[length(VDJ.stats)+1] <- length(unique(contig.table[[k]]$barcode))
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr unique barcodes"

      #generate lookup table with HC and LC counts and stats per barcode
      barcodes <- c()
      nr_HC <- c()
      nr_LC <- c()
      is_cell <-c()
      high_confidence <- c()
      productive <- c()
      full_length <- c()
      nr_bar <- 0

      for(j in unique(contig.table[[k]]$barcode)){
        nr_bar <- nr_bar + 1
        barcodes <- append(barcodes, j)
        is_cell <- append(is_cell, min(contig.table[[k]]$is_cell[which(contig.table[[k]]$barcode == j)])) #because most barcodes have two contigs (1HC, 1LC), the min function is used. Normally both contigs have the same "quality stats". In case they do not, the min function always chooses FALSE if present.
        high_confidence <- append(high_confidence, min(contig.table[[k]]$high_confidence[which(contig.table[[k]]$barcode == j)]))
        productive <- append(productive, min(contig.table[[k]]$productive[which(contig.table[[k]]$barcode == j)]))
        full_length <- append(full_length, min(contig.table[[k]]$full_length[which(contig.table[[k]]$barcode == j)]))

        nr_HC <- append(nr_HC,stringr::str_count(paste0(contig.table[[k]]$chain[which(contig.table[[k]]$barcode == j)],collapse = ""), "IGH"))
        nr_LC <- append(nr_LC,stringr::str_count(paste0(contig.table[[k]]$chain[which(contig.table[[k]]$barcode == j)],collapse = ""), "IG(K|L)"))
      }

      #compatibitlity with older cellranger versions using "True" instead of "true"
      is_cell <- tolower(is_cell)
      high_confidence <- tolower(high_confidence)
      productive <- tolower(productive)
      full_length <- tolower(full_length)

      lookup_stats <- data.frame(barcodes,nr_HC,nr_LC,is_cell,high_confidence,productive,full_length)
      names(lookup_stats) <- c("barcodes","nr_HC","nr_LC","is_cell","high_confidence","productive","full_length")

      #generate lookup table for clonotypes
      clonotype_ids <- c()
      nr_HC <- c()
      nr_LC <- c()
      for(l in 1:nrow(clonotype.list[[k]])){
        clonotype_ids <- append(clonotype_ids, clonotype.list[[k]]$clonotype_id[l])
        nr_HC <- append(nr_HC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IGH:"))
        nr_LC <- append(nr_LC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IG(K|L):"))
      }
      lookup_stats_clono <- data.frame(clonotype_ids,nr_HC,nr_LC)
      names(lookup_stats_clono) <- c("clonotype_ids","nr_HC","nr_LC")

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

      #If needed add more stats here...

      #percentages !Column counts are hardcoded! Adjust if adding new stats
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
    VDJ.stats.all <- do.call(rbind, VDJ.stats.list) #Bind dataframes from all samples

    tryCatch({
      VDJ.metrics.all <- "none" #for error catching later
      if(!inherits(vdj.metrics[[1]],"character")){ #10x stats VDJ table provided

        if(verbose) message("Getting 10x stats")

        VDJ.metrics.list <- vdj.metrics

        #for VDJ
        #check length
        #add rep identifier
        for(ij in 1:length(VDJ.metrics.list)){
          VDJ.metrics.list[[ij]]$rep_id <- ij
          for(ik in 1:ncol(VDJ.metrics.list[[ij]])){
            VDJ.metrics.list[[ij]][,ik] <- as.character(VDJ.metrics.list[[ij]][,ik])
          }
        }

        #### To make sure that different samples from different cellranger versions will be bound together
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
              ab_1 <- suppressWarnings(merge(ab_1, cur_ab, by = "idents", all.x = T, all.y = T))
            }
          }
          VDJ.metrics.all <- as.data.frame(t(ab_1)[2:ncol(ab_1),])
          names(VDJ.metrics.all) <- ab_1$idents
        }
      } #End if vdj.metrics[[1]] != "none"

      #for GEX
      #Matching dataframes from potentially different cellranger versions and running modes
      GEX.metrics.list <- gex.metrics
      GEX.metrics.all <- "none" #for error catching later
      if(!inherits(GEX.metrics.list[[1]],"character")){

        #add rep identifier
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
              ab_1 <- suppressWarnings(merge(ab_1, cur_ab, by = "idents", all.x = T, all.y = T))
            }
          }
          GEX.metrics.all <- as.data.frame(t(ab_1)[2:ncol(ab_1),])
          names(GEX.metrics.all) <- ab_1$idents
        }
      }

      if(!inherits(VDJ.metrics.all,"character") & !inherits(GEX.metrics.all, "character")){
        #bind the two
        VDJ.metrics.all <- cbind(VDJ.metrics.all, GEX.metrics.all)
      } else if (inherits(VDJ.metrics.all,"character") & !inherits(GEX.metrics.all,"character")){ #got only GEX metrics
        VDJ.metrics.all <- GEX.metrics.all
      } else if (!inherits(VDJ.metrics.all,"character") & inherits(GEX.metrics.all,"character")){ #got only VDJ metrics
        #VDJ.metrics.all <- VDJ.metrics.all no reassignment neccessary
      } else { #got none
        VDJ.metrics.all <- "none"
      }

    }, error = function(e){
      message(paste0("Adding 10x metrix failed"))
      message(e)
      VDJ.metrics.all <- "none"})

    if(!inherits(VDJ.metrics.all,"character")){ #conditional, only if we got at least one of VDJ and GEX 10x metrics
      VDJ.stats.all <- cbind(VDJ.stats.all, VDJ.metrics.all)
    }
    return(VDJ.stats.all)
  }

  #### Function def: GEX_automate_single ####

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
                                  group.id,
                                  verbose){

    if(integration.method == "scale.data" & length(GEX.list) > 1){
      if(verbose) message("Integrating GEX matrices using the default scale.data function. Other options are 'sct', 'anchors' (recommended in case of batch effects) and 'harmony' (recommended for large datasets) \n")
    }

    if(GEX.integrate == T & length(GEX.list) > 1 & integration.method != "anchors"){ #combine all GEX into one seurat object and add s%number%_ to the FRONT of the barcode
      #In case of the ANCORS integration method, we DO NOT combine datasets here. If integration.method == "anchors" this will execute the else condition below
      GEX.merged <- GEX.list[[1]]
      GEX.list[[1]] <- "none"
      GEX.merged <- SeuratObject::RenameCells(GEX.merged, new.names = paste0("s",1,"_",gsub("(^_)|(-\\d+.*$)","",colnames(GEX.merged))))
      GEX.merged@meta.data$sample_id <- paste0("s",1)
      GEX.merged@meta.data$group_id <- group.id[1]
      for(i in 2:length(GEX.list)){
        GEX.list[[i]] <- SeuratObject::RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",gsub("(^_)|(-\\d+.*$)","",colnames(GEX.list[[i]]))))

        GEX.list[[i]]@meta.data$sample_id <- paste0("s",i)
        GEX.list[[i]]@meta.data$group_id <- group.id[i]
        GEX.merged <- merge(GEX.merged, y = GEX.list[[i]], add.cell.ids = c("",""))
        GEX.list[[i]] <- "none"
      }

      GEX.list <- list() #making this into a list item to make the downstream process uniform
      GEX.list[[1]] <- GEX.merged
      GEX.list[[1]] <- SeuratObject::RenameCells(GEX.list[[1]], new.names = gsub("(^_)","",colnames(GEX.list[[1]])))

    } else {
      for(i in 1:length(GEX.list)){ #or do not integrate, but still add the sample identifier to the FRONT of the barcode. This is to make the output uniform and to deal with the possibility of integrating VDJ but not GEX
        GEX.list[[i]] <- SeuratObject::RenameCells(GEX.list[[i]], new.names = paste0("s",i,"_",gsub("(^_)|(-\\d+.*$)","",colnames(GEX.list[[i]]))))
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
        message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Calculating TSNE embedding"))
        tsne_error <- tryCatch({
          GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[1]], dims = mds.dim,check_duplicates=F)}, error = function(e){
            message(paste0("TSNE failed: ", e))})

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
        message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Calculating TSNE embedding"))
        tsne_error <- tryCatch({
          GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[1]], dims = mds.dim,check_duplicates=F)}, error = function(e){
            message(paste0("TSNE failed: ", e))})

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
        message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Calculating TSNE embedding"))
        tsne_error <- tryCatch({
        GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[i]], dims = mds.dim,check_duplicates=F)}, error = function(e){
          message(paste0("TSNE failed: ", e))})

      }
    } #END LOOP OVER GEX list elements

    if(integration.method=="anchors"){

      #Lapply to normalize data and find variable features for each matrix in the GEX list
      GEX.list <- lapply(GEX.list, function(x, norm.scale.factor, n.variable.features){
        x <- Seurat::NormalizeData(x, scale.factor = norm.scale.factor)
        x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = n.variable.features)
      }, norm.scale.factor, n.variable.features) #these are the additional arguments intput to the lapply call
      #Now select integration features
      int_feat <- Seurat::SelectIntegrationFeatures(GEX.list)
      GEX.merged <- Seurat::FindIntegrationAnchors(GEX.list, anchor.features = int_feat)
      #finally integrating the data and setting the default assay to the integrated one
      GEX.merged <- Seurat::IntegrateData(GEX.merged)
      Seurat::DefaultAssay(GEX.merged) <- "integrated"
      #Matching the formatting of the rest of integration methods
      GEX.list <- list()
      GEX.list[[1]] <- GEX.merged
      GEX.merged <- NULL

      #now calculate embeddings
      all.genes <- rownames(GEX.list[[1]])
      GEX.list[[1]] <- Seurat::ScaleData(GEX.list[[1]], features = all.genes)
      GEX.list[[1]] <- Seurat::RunPCA(GEX.list[[1]], features = Seurat::VariableFeatures(object = GEX.list[[1]]))
      GEX.list[[1]] <- Seurat::FindNeighbors(GEX.list[[1]], dims = neighbor.dim)
      GEX.list[[1]] <- Seurat::FindClusters(GEX.list[[1]], resolution = cluster.resolution)
      GEX.list[[1]] <- Seurat::RunUMAP(GEX.list[[1]], dims = mds.dim)
      message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Calculating TSNE embedding"))
      tsne_error <- tryCatch({
        GEX.list[[i]] <- Seurat::RunTSNE(GEX.list[[1]], dims = mds.dim,check_duplicates=F)}, error = function(e){
          message(paste0("TSNE failed: ", e))})
    }
    return(GEX.list)
  }

  #### Function def: barcode_VDJ_iteration ####

  #Helper function called in VDJ_GEX_matrix. Do not run as standalone!
  #FUN to call in parlapply mclapply or lapply
  barcode_VDJ_iteration <- function(barcodes,
                                    contigs,
                                    clonotypes,
                                    references,
                                    annotations,
                                    gap.opening.cost,
                                    gap.extension.cost,
                                    trim.and.align,
                                    select.excess.chains.by.umi.count,
                                    excess.chain.confidence.count.threshold){

    #Get all the info needed to shrink data usage and search times later in the function
    #Filtering out non productive or non full length contigs from cell. This is neccessary, as a cell labled as productive and full length may still have associated contigs not fullfilling these criteria.
    curr.contigs.HC <- contigs[contigs$barcode == barcodes
                                  & tolower(contigs$is_cell) == "true"
                                  & tolower(contigs$high_confidence) == "true"
                                  & tolower(contigs$productive) == "true"
                                  & tolower(contigs$full_length) == "true"
                                  & stringr::str_detect(contigs$chain, pattern = "(TRG|TRB|IGH)"),]
    curr.contigs.LC <- contigs[contigs$barcode == barcodes
                                  & tolower(contigs$is_cell) == "true"
                                  & tolower(contigs$high_confidence) == "true"
                                  & tolower(contigs$productive) == "true"
                                  & tolower(contigs$full_length) == "true"
                                  & stringr::str_detect(contigs$chain, pattern = "(TRD|TRA|IGK|IGL)"),]

    #Get number of chains
    HC_count <- nrow(curr.contigs.HC)
    LC_count <- nrow(curr.contigs.LC)

    #In case the cell has more than 1 VJ or VDJ chain, and the user decided to not want to have all these chains included in the later table
    #select one of the excessive VDJ chains by umi count
    if(HC_count > 1 & select.excess.chains.by.umi.count == T){
      if(!all(as.numeric(curr.contigs.HC$umis) >= excess.chain.confidence.count.threshold)){ #only filter chains out if not all chains present exceed the set confidence umi count threshold.
        curr.contigs.HC <- curr.contigs.HC[which.max(curr.contigs.HC$umis),] #getting the VDJ with most umis and all VJs
        HC_count <- 1 #resetting the count to not confuse any loops below
      } #no else needed: if all chains exceed the chain.confidence.count.threshold all chains are kept
    }
    #select one of the excessive VJ chains by umi count
    if(LC_count > 1 & select.excess.chains.by.umi.count == T){
      if(!all(as.numeric(curr.contigs.LC$umis) >= excess.chain.confidence.count.threshold)){ #only filter chains out if not all chains present exceed the set confidence umi count threshold.
        curr.contigs.LC <- curr.contigs.LC[which.max(curr.contigs.LC$umis),] #getting the VDJ with most umis and all VJs
        LC_count <- 1 #resetting the count to not confuse any loops below
      } #no else needed: if all chains exceed the chain.confidence.count.threshold all chains are kept
    }

    #only getting references if clonotype id is present for that cell
    if(HC_count > 0) {raw_clonotype_id <- curr.contigs.HC$raw_clonotype_id[1]
    } else {raw_clonotype_id <- curr.contigs.LC$raw_clonotype_id[1]
    }
    if(raw_clonotype_id != ''){
      curr.references <- references[which(stringr::str_detect(names(references), raw_clonotype_id))]
    } else {curr.references <- ""}

    #getting the relevant annotations ! ONLY IF TRIM AND ALIGN IS TRUE (check the VDJ loading part for reference)
    if(trim.and.align == T){
      curr.annotations <- annotations[stringr::str_detect(annotations$contig_id, barcodes),]
    }

    #set up data structure
    cols <- c("barcode","sample_id","group_id","clonotype_id_10x","clonotype_id","clonotype_frequency","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa", "VJ_cdr3s_aa","VDJ_cdr3s_nt", "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain","VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw","VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_raw_ref","VJ_raw_ref","VDJ_trimmed_ref","VJ_trimmed_ref")
    curr.barcode <- stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 1)), cols)

    #fill in information that do not need processing
    #Contig info on light/a and heavy/b chains is put into different columns (see cols)
    #If only one contig is available, the fields of the other are left blank
    #If more than two contigs of one chain (e.g. 2 TRB) are present, the elements will be pasted separated by a ";" into the relevant fields (in the case of TRB, into the Hb columns)

    #In this case we need to make much less effort with pasting together, so we can save time
    if(HC_count == 1 & LC_count == 1){

      #fill in the pasted info to curr.barcode directly
      curr.barcode$barcode <- curr.contigs.HC$barcode[1]
      curr.barcode$clonotype_id_10x <- curr.contigs.HC$raw_clonotype_id[1]
      cl_freq <- clonotypes[(clonotypes[,1] == curr.barcode$clonotype_id_10x),2]
      if(length(cl_freq) > 0){curr.barcode$clonotype_frequency <- clonotypes[(clonotypes[,1] == curr.barcode$clonotype_id_10x),2]} else{curr.barcode$clonotype_frequency <- 0}
      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count
      curr.barcode$VDJ_cdr3s_aa <- curr.contigs.HC$cdr3
      curr.barcode$VJ_cdr3s_aa <- curr.contigs.LC$cdr3
      curr.barcode$VDJ_cdr3s_nt <- curr.contigs.HC$cdr3_nt
      curr.barcode$VJ_cdr3s_nt <- curr.contigs.LC$cdr3_nt
      curr.barcode$VDJ_chain_contig <- curr.contigs.HC$contig_id
      curr.barcode$VJ_chain_contig <- curr.contigs.LC$contig_id
      curr.barcode$VDJ_chain <- curr.contigs.HC$chain
      curr.barcode$VJ_chain <- curr.contigs.LC$chain
      curr.barcode$VDJ_vgene <- curr.contigs.HC$v_gene
      curr.barcode$VJ_vgene <- curr.contigs.LC$v_gene
      curr.barcode$VDJ_dgene <- curr.contigs.HC$d_gene
      curr.barcode$VDJ_jgene <- curr.contigs.HC$j_gene
      curr.barcode$VJ_jgene <- curr.contigs.LC$j_gene
      curr.barcode$VDJ_cgene <- curr.contigs.HC$c_gene
      curr.barcode$VJ_cgene <- curr.contigs.LC$c_gene
      curr.barcode$VDJ_raw_consensus_id <- curr.contigs.HC$raw_consensus_id
      curr.barcode$VJ_raw_consensus_id <- curr.contigs.LC$raw_consensus_id
      #adding raw sequences directly
      curr.barcode$VDJ_sequence_nt_raw <- curr.contigs.HC$raw_contig
      curr.barcode$VJ_sequence_nt_raw <- curr.contigs.LC$raw_contig

    } else { # this for cells with abberrant chain numbers

      contigs_pasted.HC <- stats::setNames(data.frame(matrix(ncol = ncol(curr.contigs.HC), nrow = 1)), names(curr.contigs.HC))
      contigs_pasted.LC <- stats::setNames(data.frame(matrix(ncol = ncol(curr.contigs.LC), nrow = 1)), names(curr.contigs.LC))

      #Heavy/b chain count
      if(HC_count == 0){
        contigs_pasted.HC[1,] <- ""
      } else if(HC_count == 1){
        contigs_pasted.HC[1,] <- curr.contigs.HC[1,]
      } else if(HC_count > 1){
        for(k in 1:ncol(curr.contigs.HC)){
          contigs_pasted.HC[1,k] <- paste0(curr.contigs.HC[,k], collapse = ";")
        }
      }

      #Light/a chain count
      if(LC_count == 0){
        contigs_pasted.LC[1,] <- ""
      } else if(LC_count == 1){
        contigs_pasted.LC[1,] <- curr.contigs.LC[1,]
      } else if(LC_count > 1){
        for(k in 1:ncol(curr.contigs.LC)){
          contigs_pasted.LC[1,k] <- paste0(curr.contigs.LC[,k], collapse = ";")
        }
      }

      #fill in the pasted info to curr.barcode
      if(HC_count > 0){
        curr.barcode$barcode <- curr.contigs.HC$barcode[1]
        curr.barcode$clonotype_id_10x <- curr.contigs.HC$raw_clonotype_id[1]
      } else {
        curr.barcode$barcode <- curr.contigs.LC$barcode[1]
        curr.barcode$clonotype_id_10x <- curr.contigs.LC$raw_clonotype_id[1]
      }
      cl_freq <- clonotypes[(clonotypes[,1] == curr.barcode$clonotype_id_10x),2]
      if(length(cl_freq) > 0){curr.barcode$clonotype_frequency <- clonotypes[(clonotypes[,1] == curr.barcode$clonotype_id_10x),2]} else{curr.barcode$clonotype_frequency <- 0}
      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count
      curr.barcode$VDJ_cdr3s_aa <- contigs_pasted.HC$cdr3
      curr.barcode$VJ_cdr3s_aa <- contigs_pasted.LC$cdr3
      curr.barcode$VDJ_cdr3s_nt <- contigs_pasted.HC$cdr3_nt
      curr.barcode$VJ_cdr3s_nt <- contigs_pasted.LC$cdr3_nt
      curr.barcode$VDJ_chain_contig <- contigs_pasted.HC$contig_id
      curr.barcode$VJ_chain_contig <- contigs_pasted.LC$contig_id
      curr.barcode$VDJ_chain <- contigs_pasted.HC$chain
      curr.barcode$VJ_chain <- contigs_pasted.LC$chain
      curr.barcode$VDJ_vgene <- contigs_pasted.HC$v_gene
      curr.barcode$VJ_vgene <- contigs_pasted.LC$v_gene
      curr.barcode$VDJ_dgene <- contigs_pasted.HC$d_gene
      curr.barcode$VDJ_jgene <- contigs_pasted.HC$j_gene
      curr.barcode$VJ_jgene <- contigs_pasted.LC$j_gene
      curr.barcode$VDJ_cgene <- contigs_pasted.HC$c_gene
      curr.barcode$VJ_cgene <- contigs_pasted.LC$c_gene
      curr.barcode$VDJ_raw_consensus_id <- stringr::str_split(contigs_pasted.HC$raw_consensus_id,";",simplify = T)[1]
      curr.barcode$VJ_raw_consensus_id <- stringr::str_split(contigs_pasted.LC$raw_consensus_id,";",simplify = T)[1] #Because we may have more than one consensus ID for light chains, we need to get one of them. Because the consensus ids are always the same for different light or heavy chains of the same cell, we can just take the first element of the str_split
      #adding raw sequences directly
      curr.barcode$VDJ_sequence_nt_raw <- contigs_pasted.HC$raw_contig
      curr.barcode$VJ_sequence_nt_raw <- contigs_pasted.LC$raw_contig

    } #end if HC | LC count > 1

    #Now on to the actual sequences
    reference_HC <- curr.references[names(curr.references) == curr.barcode$VDJ_raw_consensus_id]
    reference_LC <- curr.references[names(curr.references) == curr.barcode$VJ_raw_consensus_id]

    #HEAVY CHAIN / TRB
    if(HC_count == 1){
      if(append.raw.reference == T){
        if(length(as.character(reference_HC)) > 0){curr.barcode$VDJ_raw_ref <- as.character(reference_HC)}
      }

      if(trim.and.align == T){
        tryCatch({
          #find match in annotations
          HC_contig <- which(curr.annotations$contig_id == curr.barcode$VDJ_chain_contig) #This line is never reached if trim.and.align == F (either set by the user or set by the function if the all_contig_annotations.json was present in the input folders)
          #trim sequence
          curr.barcode$VDJ_sequence_nt_trimmed <- substr(curr.barcode$VDJ_sequence_nt_raw, as.numeric(curr.annotations$temp_start[HC_contig])+1, as.numeric(curr.annotations$temp_end[HC_contig])-1)
          #translate trimmed sequence
          if(nchar(curr.barcode$VDJ_sequence_nt_trimmed) > 1){
            curr.barcode$VDJ_sequence_aa <- as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VDJ_sequence_nt_trimmed))))
          } else {to_paste_aa <- ""}
          #align to reference and trim reference
          if(nchar(curr.barcode$VDJ_sequence_nt_trimmed) > 1){
            alignments <- Biostrings::pairwiseAlignment(curr.barcode$VDJ_sequence_nt_trimmed, as.character(reference_HC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
            curr.barcode$VDJ_trimmed_ref <- as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))]))
          } else {to_paste_ref_trimmed <-  ""}
        }, error=function(e){
          to_paste_ref_trimmed <- "ALIGNMENT ERROR"})
      }
    } else if(HC_count > 1){ #MORE THAN ONE HC
      #from the annotations extract sequence and paste

      if(append.raw.reference == T){
        if(length(as.character(reference_HC)) > 0){
          curr.barcode$VDJ_raw_ref <- paste0(lapply(reference_HC, function(x) return(as.character(unlist(x)))), collapse = ";")
        }
      }

      #Heavy/b
      to_paste <- curr.barcode$VDJ_sequence_nt_raw
      to_paste_ref_raw <- reference_HC
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      if(trim.and.align == T){
        tryCatch({
          #looping contigs in annotation
          for(l in 1:nrow(curr.annotations)){
            #looping over Hb contig ids (as there may be more than 1)
            for(c in 1:length(stringr::str_split(curr.barcode$VDJ_chain_contig, ";",simplify = T))){
              #find a match
              if(curr.annotations$contig_id[l] == stringr::str_split(curr.barcode$VDJ_chain_contig, ";",simplify = T)[c]){
                #trim sequence
                to_paste_trimmed <- append(to_paste_trimmed, substr(stringr::str_split(to_paste, ";",simplify = T)[c], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
                #translate trimmed sequence
                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  to_paste_aa <- append(to_paste_aa, as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)])))))
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
          curr.barcode$VDJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
          curr.barcode$VDJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
          curr.barcode$VDJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
        }, error=function(e){
          to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")})
      }
    }

    #Light/a
    if(LC_count == 1){

      if(append.raw.reference == T){
        if(length(as.character(reference_LC)) > 0){curr.barcode$VJ_raw_ref <- as.character(reference_LC)
        }
      }

      if(trim.and.align == T){
        tryCatch({
          #find match in annotations
          LC_contig <- which(curr.annotations$contig_id == curr.barcode$VJ_chain_contig)
          #trim sequence
          curr.barcode$VJ_sequence_nt_trimmed <- substr(curr.barcode$VJ_sequence_nt_raw, as.numeric(curr.annotations$temp_start[LC_contig])+1, as.numeric(curr.annotations$temp_end[LC_contig])-1)
          #translate trimmed sequence
          if(nchar(curr.barcode$VJ_sequence_nt_trimmed) > 1){
            curr.barcode$VJ_sequence_aa <- as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAStringSet(curr.barcode$VJ_sequence_nt_trimmed))))
          } else {to_paste_aa <- ""}
          #align to reference and trim reference
          if(nchar(curr.barcode$VJ_sequence_nt_trimmed) > 1){
            alignments <- Biostrings::pairwiseAlignment(curr.barcode$VJ_sequence_nt_trimmed, as.character(reference_LC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
            curr.barcode$VJ_trimmed_ref <- as.character(Biostrings::subject(alignments[which.max(Biostrings::score(alignments))]))
          } else {to_paste_ref_trimmed <-  ""}
        }, error=function(e){
          to_paste_ref_trimmed <- "ALIGNMENT ERROR"
        })
      }
    } else if(LC_count > 1){ #MORE THAN ONE LC

      if(append.raw.reference == T){
        if(length(as.character(reference_LC)) > 0){
          curr.barcode$VJ_raw_ref <- paste0(lapply(reference_LC, function(x) return(as.character(unlist(x)))), collapse = ";")
        }
      }

      to_paste <- curr.barcode$VJ_sequence_nt_raw
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      if(trim.and.align == T){
        tryCatch({
          #looping contigs in annotation
          for(l in 1:nrow(curr.annotations)){
            #looping over Hb contig ids (as there may be more than 1)
            for(c in 1:length(stringr::str_split(curr.barcode$VJ_chain_contig, ";",simplify = T))){
              #find a match
              if(curr.annotations$contig_id[l] == stringr::str_split(curr.barcode$VJ_chain_contig, ";",simplify = T)[c]){
                #trim sequence
                to_paste_trimmed <- append(to_paste_trimmed, substr(stringr::str_split(to_paste, ";",simplify = T)[c], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
                #translate trimmed sequence
                if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
                  to_paste_aa <- append(to_paste_aa, as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)])))))
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
          curr.barcode$VJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
          curr.barcode$VJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
          curr.barcode$VJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
        }, error=function(e){
          to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
        })
      }
    }
    return(curr.barcode)
  }

  #### Missing parameter defaults ####

  #Input
  if(missing(GEX.read.h5)) GEX.read.h5 <- FALSE
  #VDJ GEX integration
  if(missing(VDJ.combine)) VDJ.combine <- T
  if(missing(GEX.integrate)) GEX.integrate <- T
  if(GEX.integrate == F) VDJ.combine <- F #Default to avoid one GEX object and multiple VDJs
  if(missing(integrate.GEX.to.VDJ)) integrate.GEX.to.VDJ <- T
  if(missing(integrate.VDJ.to.GEX)) integrate.VDJ.to.GEX <- T
  if(missing(exclude.GEX.not.in.VDJ)) exclude.GEX.not.in.VDJ <- F
  if(missing(filter.overlapping.barcodes.GEX)) filter.overlapping.barcodes.GEX <- T
  if(missing(filter.overlapping.barcodes.VDJ)) filter.overlapping.barcodes.VDJ <- T
  #Stats
  if(missing(get.VDJ.stats)) get.VDJ.stats <- T
  #VDJ
  if(missing(append.raw.reference)) append.raw.reference <- T
  if(missing(select.excess.chains.by.umi.count)) select.excess.chains.by.umi.count <- F
  if(missing(excess.chain.confidence.count.threshold)) excess.chain.confidence.count.threshold <- 1000 #default is too high for any single chain to cross it. if set to 3 or lower, cells with double chains will be maintained if each double chain is over that umi threshold
  if(missing(trim.and.align)) trim.and.align <- F
  if(missing(parallel.processing)) parallel.processing <- "none"
  if(parallel.processing == "parlapply" | parallel.processing == "mclapply"){
    if(missing(numcores)) numcores <- parallel::detectCores()
    if(numcores > parallel::detectCores()){numCores <- parallel::detectCores()}
  } else{
    numcores <- 1
  }
  if(missing(gap.opening.cost)) gap.opening.cost <- 10
  if(missing(gap.extension.cost)) gap.extension.cost <- 4
  #GEX
  if(missing(exclude.on.cell.state.markers)) exclude.on.cell.state.markers <- "none"
  if(missing(exclude.on.barcodes)) exclude.on.barcodes <- "none"
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
  #FB
  if(missing(FB.ratio.threshold)) FB.ratio.threshold <- 2
  if(missing(FB.count.threshold)) FB.count.threshold <- 10
  if(missing(FB.exclude.pattern)) FB.exclude.pattern <- "not a pattern that will be filterd out"
  if(missing(subsample.barcodes)) subsample.barcodes = F
  if(missing(verbose)) verbose <- T

  #### Input organisation ####

  #Input presence checks will later be run by class checks. Not-provided inputs are set to "none" / class = "character"
  if(missing(Data.in)) Data.in <- "none"
  if(missing(Seurat.in)) Seurat.in <- "none"
  if(missing(FB.out.directory.list)) FB.out.directory.list <- "none"

  #Input processed Seurat object
  if(inherits(Seurat.in,"list")){
    stop("Please provide a single Seurat object as Seurat.in, and set VDJ.combine to TRUE if processing multiple VDJ samples")

  } else if (inherits(Seurat.in,"SeuratObject") | inherits(Seurat.in,"Seurat")){ #Seurat.in must be class Seurat.object

      if(!"sample_id" %in% names(Seurat.in@meta.data) | !"group_id" %in% names(Seurat.in@meta.data)) stop("Seurat.in objects need to contain sample_id and group_id columns")
      if(stringr::str_detect(Seurat.in@meta.data$sample_id, "s\\d")[1] == F) stop("Seurat.in objects sample_id column needs to follow sample naming scheme: s1 , s2, ... sn")

    #Setting related variables accoring to Seurat.input
    seurat.loaded <- T

    if(VDJ.combine == F){
      warning("Seurat input object provided, but VDJ.combine == FALSE. Setting VDJ.combine to TRUE")
      VDJ.combine <- T
    }
  }

  #No data.in input. Checking input cellranger paths
  if(inherits(Data.in,"character")){ #if Data.in was provided, this would be false as Data.in would be of class list
    samples.in <- "none"

    #Worst case: no input was provided at all
    if(missing(VDJ.out.directory.list) & missing(GEX.out.directory.list)){ stop("Please provide data input either as a as a list of local paths to VDJ.out.directory.list and/or GEX.out.directory.list or as list of R objects to Data.in (development only)")

      #VDJ but not GEX local paths provided
    } else if(!missing(VDJ.out.directory.list) & missing(GEX.out.directory.list)){
      GEX.out.directory.list <- "none"
      samples.paths.VDJ <- paste0(do.call("c",as.list(VDJ.out.directory.list)), collapse = " ; ") #For runtime parameters
      #Defining group_id and batches
      if(missing(group.id)) group.id <- 1:length(VDJ.out.directory.list)
      batches <- "none"

      #FB control
      if(!inherits(FB.out.directory.list[[1]], "character")){
        if(length(FB.out.directory.list) != length(VDJ.out.directory.list)){
          stop("Different number of input elements for VDJ and FB. If for some samples in VDJ have no adjunct FB information, please set the corresponding FB input list element to 'PLACEHOLDER'")
        }}

      #GEX but not VDJ local paths provided
    } else if(missing(VDJ.out.directory.list) & !missing(GEX.out.directory.list)){
      VDJ.out.directory.list <- "none"
      samples.paths.GEX <- paste0(do.call("c",as.list(GEX.out.directory.list)), collapse = " ; ") #For runtime parameters
      #Defining group_id and batches
      if(missing(group.id)) group.id <- 1:length(GEX.out.directory.list)
      batches <- "none"

      #FB control
      if(!inherits(FB.out.directory.list[[1]], "character")){
        if(length(FB.out.directory.list) != length(GEX.out.directory.list)){
          stop("Different number of input elements for GEX and FB. If for some samples in GEX have no adjunct FB information, please set the corresponding FB input list element to 'PLACEHOLDER'")
        }}

      #Both local paths provided
    } else if(!missing(VDJ.out.directory.list) & !missing(GEX.out.directory.list) & !seurat.loaded){
      if(length(VDJ.out.directory.list) != length(GEX.out.directory.list)){stop("Different number of input paths provided for GEX and VDJ. Please revise input")}

    #FB control
    if(!inherits(FB.out.directory.list[[1]], "character")){
      if(length(FB.out.directory.list) != length(GEX.out.directory.list)){
        stop("Different number of input elements for GEX/VDJ and FB. If for some samples in GEX/VDJ have no adjunct FB information, please set the corresponding FB input list element to 'PLACEHOLDER'")
      }
      }

    if(missing(group.id)) group.id <- 1:length(GEX.out.directory.list)
    batches <- "none"
    }

    #For runtime parameters
    samples.paths.VDJ <- paste0(do.call("c",as.list(VDJ.out.directory.list)), collapse = " ; ")
    samples.paths.GEX <- paste0(do.call("c",as.list(GEX.out.directory.list)), collapse = " ; ")


    #! Data.in input provided
  } else if(inherits(Data.in,"list")){ #This input is prioritized over local paths input

    GEX.out.directory.list <- "none"
    VDJ.out.directory.list <- "none"
    FB.out.directory.list <- "none"

    #figure out the input structure and reorder into a list were elements 1-n are all files for each sample
    samples.in <- list()
    samples.paths.VDJ <- c()
    samples.paths.GEX <- c()
    batches <- c()

    for(i in 1:length(Data.in)){ #First level
      if(names(Data.in[[i]])[1] == "VDJ"){ #check if we are already on a sample level
        samples.in[[length(samples.in)+1]] <- Data.in[[i]]
        samples.paths.VDJ <- c(samples.paths.VDJ, Data.in[[i]][[4]])
        samples.paths.GEX <- c(samples.paths.GEX, Data.in[[i]][[5]])
        batches <- c(batches, Data.in[[i]][[5]])
        Data.in[[i]] <- "None" #to limit ram usage
      } else {
        for(j in 1:length(Data.in[[i]])){ #Second level
          if(names(Data.in[[i]][[j]])[1] == "VDJ"){ #check if we are a sample level
            samples.in[[length(samples.in)+1]] <- Data.in[[i]][[j]]
            samples.paths.VDJ <- c(samples.paths.VDJ, Data.in[[i]][[j]][[4]])
            samples.paths.GEX <- c(samples.paths.GEX, Data.in[[i]][[j]][[5]])
            batches <- c(batches, Data.in[[i]][[j]][[5]])
            Data.in[[i]][[j]] <- "None" #to limit ram usage
          } else { #Data structure does not match expectations
            stop("Provided datastructure does not match required input format")
          }
        }
      }
    }

    samples.paths.GEX <- paste0(do.call("c",as.list(samples.paths.GEX)), collapse = " ; ")
    samples.paths.VDJ <- paste0(do.call("c",as.list(samples.paths.VDJ)), collapse = " ; ")

    if(missing(group.id)) group.id <- 1:length(samples.in)
    Data.in <- NULL
    data.in.loaded <- T
  }

  #### Load VDJ ####
    if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","\n Loading in data"))

  vdj_load_error <- tryCatch({

  if(data.in.loaded == F & VDJ.out.directory.list[[1]] != "none"){ #No Data.in => proceed with loading by VDJ.out.directory.list

    #Remove possible backslash at the end of the input path
    VDJ.out.directory.list <- lapply(VDJ.out.directory.list, function(x) return(gsub("/$", "", x)))

    #Check that all necessary files exist before loading
    file.check <- c()
    for(i in 1:length(VDJ.out.directory.list)){
      for(j in c("/clonotypes.csv","/concat_ref.fasta","/filtered_contig_annotations.csv","/filtered_contig.fasta")){
        if(file.exists(paste(VDJ.out.directory.list[[i]],j,sep="")) == F) file.check <- c(file.check, paste0("For VDJ input directory ",i, " file ", j, " was not found. \n"))
      }
    }
    if(trim.and.align == F) file.check <- file.check[!stringr::str_detect(file.check, "all_contig_annotations.json")] #If these files are not needed they are not reported as missing
    if(length(file.check) > 0){
      for(i in file.check){ message(i)}
      stop("VDJ loading of at least one file failed. Please revise input paths")}

      clonotype.list <- lapply(VDJ.out.directory.list, function(x) utils::read.table(paste(x,"/clonotypes.csv",sep=""), stringsAsFactors = FALSE,sep=",",header=T))
      reference.list <- lapply(VDJ.out.directory.list, function(x) seqinr::read.fasta(paste(x,"/concat_ref.fasta",sep=""), as.string = T,seqonly = F,forceDNAtolower = F))
      contig.table <- lapply(VDJ.out.directory.list, function(x) utils::read.csv(paste(x,"/filtered_contig_annotations.csv",sep=""),sep=",",header=T)) #better using the table format downstream

      #changes for compatibility with the cellranger multi pipeline
      if(all(file.exists(paste(VDJ.out.directory.list,"/metrics_summary.csv",sep="")))){
        metrics.table <- lapply(VDJ.out.directory.list, function(x) utils::read.csv(paste(x,"/metrics_summary.csv",sep=""),sep=",",header=T))
      } else {
        if(verbose) message("! metrics_summary.csv file not available in at least one of the VDJ input directories. Loading will be skipped \n")
      }

      #NEW in cellranger 6.1 => Loading the contigs independently from the annotations file. This allows to return full contig sequences despite not trimming and aligning
      #Adding the raw contigs to the contig table loaded from filtered_contig_annotation.csv above
      raw.contig.table <- lapply(VDJ.out.directory.list, function(x) seqinr::read.fasta(paste(x,"/filtered_contig.fasta",sep=""), as.string = T,seqonly = F,forceDNAtolower = F))
      for(ikj in 1:length(raw.contig.table)){
        #coercing this to a dataframe
        raw.contig.table[[ikj]] <- data.frame("contig_id" = names(raw.contig.table[[ikj]]), "raw_contig" = as.character(raw.contig.table[[ikj]]))
        #merge directly with contig dataframe to have annotations and sequence in one place for later
        contig.table[[ikj]] <- merge(contig.table[[ikj]], raw.contig.table[[ikj]], by = "contig_id", all.x = T, all.y = F) #making sure to only merge in raw contig sequences for contigs which are present in the table containing annotations
        if(sum(is.na(contig.table[[ikj]]$raw_contig)) > 0.5*nrow(contig.table[[ikj]])) warning("! Merging of raw contigs and filtered_contig_annotations showed unsusually low overlap \n")
      }

      #NEW in cellranger 6.1 =>
      #Reason: In Cellranger 6 the function Cellranger multi was introduced. The VDJ output of that does not contain the all_contig_annotations.json file.
      #Therefore we check if the file exists and skip it if not. If it does not exist, the function will not be able to perform trimming and aligning and a warning is issued

      if(trim.and.align == T){ #we only need this for trimming, so we skip the loading if that is not desired
        if(all(file.exists(paste(VDJ.out.directory.list,"/all_contig_annotations.json",sep="")))){
          annotations.list <- lapply(VDJ.out.directory.list, function(x) jsonlite::read_json(paste(x,"/all_contig_annotations.json",sep="")))

          # returns key features: featureRegions, and of featureRegions. Used for trimming where the V region starts and where the C region ends.
          annotations.table <- list()
          for(i in 1:length(annotations.list)){

            #get annotation table to make VDJ_barcode_iteration later on faster
            annotations.table[[i]] <- do.call(BiocGenerics::rbind, lapply(annotations.list[[i]], function(y){

                if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")]) == 1){
                  temp_start <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")][[1]]$contig_match_start
                } else {temp_start <- 10000 } #This is to cause substr() in trimming to return an empty string
                if(length(y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")]) == 1){
                  temp_end <- y$annotations[sapply(y$annotations, function(x) x$feature$region_type=="J-REGION")][[1]]$contig_match_end #!
                } else {temp_end <- 0 } #This is to cause substr() in trimming to return an empty string
                data.frame("contig_id" = y$contig_name,
                           "sequence" = y$sequence,
                           "temp_start" = temp_start,
                           "temp_end" = temp_end
                )}))### returns a dataframe with these four outputs.
          }
        } else { #in at least one directory the all_contig_annotations.json was not found
          warning("Warning: At least one VDJ input directory is missing the file all_contig_annotations.json. Without this, accurate trimming and aligning of sequences is not possible. Setting trim.and.align to FALSE and proceeding. For an alternate mode of aligment please refer to the function VDJ_call_MIXCR \n")
          trim.and.align <- F
          annotations.table <- as.list(rep("none", length(contig.table)))
        }
      } else if(trim.and.align == F){
        annotations.table <- as.list(rep("none", length(contig.table)))
      }

      #Reduce clonotype.list for an entry into VDJ_barcode_iteration function
      clonotype.list <- lapply(clonotype.list, function(x) return(x[,c(1,2,4)]))

      #change names so that the barcode function does not have to do that
      contig.table <- lapply(contig.table, function(x){
        x$raw_consensus_id <- gsub("_consensus_","_concat_ref_", x$raw_consensus_id)
        return(x)})

      #convert all columns of input objects to character
      contig.table <- lapply(contig.table, function(x) dplyr::mutate(.data = x, dplyr::across(dplyr::everything(), as.character)))
      if(!inherits(annotations.table[[1]],"character")){ #conditional as annotations may not have been loaded
        annotations.table <- lapply(annotations.table, function(x) dplyr::mutate(.data = x, dplyr::across(dplyr::everything(), as.character)))}

      vdj.loaded <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded VDJ data"))

  } else if(data.in.loaded == T){ #GET info from samples.in
    #get VDJs

    vdj.loaded <- F

      clonotype.list <- lapply(samples.in, function(x){return(x[[1]][[1]])})
      reference.list <- lapply(samples.in, function(x){return(x[[1]][[2]])})
      annotations.table <- lapply(samples.in, function(x){return(x[[1]][[3]])})
      contig.table <- lapply(samples.in, function(x){return(x[[1]][[4]])})
      metrics.table <- lapply(samples.in, function(x){return(x[[1]][[5]])})

      clonotype.list.class <- lapply(samples.in, function(x){return(class(x[[1]][[1]]))})

      if(!any(clonotype.list.class == "character")){ #All slots of data in contained VDJ info

      #clear VDJ part in object
      for(i in 1:length(samples.in)) samples.in[[i]][[1]] <- "loaded"

      #change names so that the barcode function does not have to do that
      contig.table <- lapply(contig.table, function(x){
        x$raw_consensus_id <- gsub("_consensus_","_concat_ref_", x$raw_consensus_id)
        return(x)})

      #Reduce clonotype.list for an entry into VDJ_barcode_iteration function
      clonotype.list <- lapply(clonotype.list, function(x) return(x[,c(1,2,4)]))

      #convert all columns of input objects to character
      contig.table <- lapply(contig.table, function(x) dplyr::mutate(.data = x, dplyr::across(dplyr::everything(), as.character)))
      if(!inherits(annotations.table[[1]],"character")){ #conditional as annotations may not have been loaded
        annotations.table <- lapply(annotations.table, function(x) dplyr::mutate(.data = x, dplyr::across(dplyr::everything(), as.character)))}

      vdj.loaded <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded VDJ data from Data.in"))

      } else {
        warning("At least on VDJ slot in data.in provided list did not contain VDJ information. Skipping VDJ processing")
        vdj.loaded <- F
      }
    }
  }, error = function(e){
    message("\n Loading VDJ failed")
    message(e)}) #end loading from Data.in

  #### Load GEX from Data.in or input paths ####

  gex_load_error <- tryCatch({

  #Choosing GEX paths
  if(data.in.loaded == F  & seurat.loaded == F & GEX.out.directory.list[[1]] != "none"){ #No Data.in or Seurat.in => proceed with loading by GEX.out.directory.list

    #Remove possible backslash at the end of the input path
    GEX.out.directory.list <- lapply(GEX.out.directory.list, function(x) return(gsub("/$", "", x)))

      if(GEX.read.h5 == F){
        if((stringr::str_detect(GEX.out.directory.list[[1]],"filtered_feature_bc_matrix") | #checking only for first path assuming that all sample inputs are processed with the same cellranger function
           stringr::str_detect(GEX.out.directory.list[[1]],"raw_feature_bc_matrix") |
           stringr::str_detect(GEX.out.directory.list[[1]],"sample_feature_bc_matrix"))) {
          #Nothing to append
        } else {
          if(dir.exists(paste(GEX.out.directory.list[[1]],"/filtered_feature_bc_matrix", sep = ""))) { #appending folder name of cellranger count output
                GEX.out.directory.list <- paste(GEX.out.directory.list, "/filtered_feature_bc_matrix", sep = "")
                if(verbose) message("\n Setting GEX directory to provided path/filtered_feature_bc_matrix")
          } else if (dir.exists(paste(GEX.out.directory.list[[1]], "/sample_feature_bc_matrix", sep = ""))) { #appending folder name of cellranger aggr or multi output
            GEX.out.directory.list <- paste(GEX.out.directory.list, "/sample_feature_bc_matrix", sep = "")
                if(verbose) message("\n Setting GEX directory to provided path/sample_feature_bc_matrix")
          } else {
            message("\n The GEX directory filtered_feature_bc_matrix or sample_feature_bc_matrix was not found at the given path. Trying to load GEX from raw input path")
            GEX.out.directory.list <- GEX.out.directory.list
          }
        }
      } else if (GEX.read.h5) {
            GEX.out.directory.list <- gsub("(/filtered_feature_bc_matrix$)|(/raw_feature_bc_matrix$)|(/sample_feature_bc_matrix$)", "", GEX.out.directory.list) #removing last path element if pointing to subfolders of sample output directory
            if(stringr::str_detect(GEX.out.directory.list[[1]],"\\.h5")){
              if(verbose) message("\n Setting GEX directory to provided .h5 file")
              #Nothing to do here
            } else {
              GEX.out.directory.list <- paste(GEX.out.directory.list, "/filtered_feature_bc_matrix.h5", sep = "")
              if(verbose) message("\n Setting GEX directory to provided path/filtered_feature_bc_matrix.h5")
            }
      }

      #Loading GEX
        if (!GEX.read.h5) {
          directory_read10x <- lapply(GEX.out.directory.list, function(x) tryCatch({Seurat::Read10X(data.dir = x)}, error = function(e){return(NULL)}))
        } else {
          directory_read10x <- lapply(GEX.out.directory.list,
                                      tryCatch({function(x) Seurat::Read10X_h5(filename = x, use.names = TRUE, unique.features = TRUE)}, error = function(e){return(NULL)}))
        }

      for(i in 1:length(directory_read10x)){
        if(is.null(directory_read10x[[i]])) stop("GEX loading failed for input path(s) ", i, "! Stopping run. Please revise input paths \n")
      }

      #NEW CELLRANGER 6.1. dealing with the possibility of Feature barcode information being included in the GEX directory.
      #=> Procedure: check if there is feature barcode info in the GEX that was just loaded => if not, proceed as normal
      #=> if yes, isolate the matrices for each input sample into two flat lists: GEX.list and FB.list
      fb.loaded <- F
      FB.list <- list()
      for(i in 1:length(directory_read10x)){ #iterating over main list
        if(inherits(directory_read10x[[i]],"list")){ #this returns true only if the GEX directory import contained more than one marices => i.e. there is a GEX and a FB matrix
          if(verbose) message(paste0("\n GEX input ", i, " contains multiple count matrices."))
          GEX_ind <- c() #open indices for GEX and FB list elements
          FB_ind <- c()
          for(j in 1:length(directory_read10x[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(directory_read10x[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              if(verbose) message(paste0("GEX input ", i, " element ", j, " contains > 100 features and will be loaded as GEX"))
              GEX_ind <- c(GEX_ind, j)
            } else if(nrow(directory_read10x[[i]][[j]]) < 100){
              if(verbose) message(paste0("GEX input ", i, " element ", j, " contains < 100 features and will be loaded as FB"))
              FB_ind <- c(FB_ind, j)
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("\n GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }

          if(length(FB_ind) > 0){
            FB.list[[length(FB.list)+1]] <- directory_read10x[[i]][[FB_ind]] #add this matrix to the FB list ! For GEX input where no FB was found this will become placeholder (see below)
            directory_read10x[[i]] <- directory_read10x[[i]][[-FB_ind]] #delete the matrix from the original list, so that only the GEX matrix remains
            #Moreover we want to prevent two user inputs with FB. Here we set fb.loaded to TRUE, so that the FB loading from disk module further down will be skipped.
            fb.loaded <- T
          }
        } else{
          #PLACEHOLDER MATRIX: can be converted to a seurat object and run through the whole function without needing extra IF conditions
          FB.list[[length(FB.list)+1]] <- Matrix::Matrix(data = c(rep(1001, 10), rep(1, 10)), nrow = 2, ncol = 10, sparse = TRUE, dimnames = list( c("No-FB-data", "column2"),LETTERS[11:20]))
        }
      }

      if(fb.loaded == T) {
        #got two lists: the FB.list and the directory_read10x which only contains GEX info, but is still nested. So we have to flatten it
        directory_read10x <- do.call(list, unlist(directory_read10x, recursive=FALSE))

        #Next we need to check column names / names of feature barcodes. Potentially more than one library of FBs was sequenced giving us more than one matrix but all with the same column names (because same FBs may have been used across samples). We therefore rename the feature barcodes individually to make sure they stay separated
        for(i in 1:length(FB.list)){
          if(ncol(FB.list[[i]]) > 0){
            #colnames(FB.list[[i]]) <- paste0("s", i, "_", colnames(FB.list[[i]]))
          } else {
            FB.list[[i]] <- Matrix::Matrix(data = c(rep(1001, 10), rep(1, 10)), nrow = 2, ncol = 10, sparse = TRUE, dimnames = list( c("No-FB-data", "column2"),LETTERS[11:20]))
            #colnames(FB.list[[i]]) <- paste0("s", i, "_", colnames(FB.list[[i]]))
            #PLACEHOLDER MATRIX: can be converted to a seurat object and run through the whole function without needing extra IF conditions
          }
        }
        #convert all to Seurat objects
        FB.list <- lapply(FB.list, function(x) Seurat::CreateSeuratObject(x))
      }

      #Continuing with GEX processing. Irregardless of FB data presence, the directory_read10x is a flat list of GEX matrices
      directory_read10x <- lapply(directory_read10x, function(x){
        rownames(x) <- toupper(rownames(x))
        return(x)})
      gex.list <- lapply(directory_read10x, function(x) Seurat::CreateSeuratObject(x))
      directory_read10x <- NULL

      #load the metrics file conditionally
      GEX.out.directory.list.metrics <- paste(GEX.out.directory.list,"/metrics_summary.csv",sep="")
      if(get.VDJ.stats == T & all(file.exists(GEX.out.directory.list.metrics))){
        gex.metrics.table <- lapply(GEX.out.directory.list.metrics, function(x) utils::read.csv(x,sep=",",header=T, ))
      } else {
        gex.metrics.table <- "none"
      }

      gex.loaded <- T
      if(fb.loaded == F){
        if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded GEX data"))
      } else {
        if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded GEX and FB data"))
      }


  } else if(data.in.loaded == T & inherits(samples.in, "list")){ #GET info from samples.in
      gex.list <- lapply(samples.in, function(x) return(x[[2]][[1]]))
      gex.list.class <- lapply(samples.in, function(x) return(class(x[[2]][[1]])))

      if(!any(gex.list.class == "character")){ #All slots of data in contained GEX info
      gex.list <- lapply(gex.list, function(x){
        rownames(x) <- toupper(rownames(x))
        return(x)})

      gex.list <- lapply(gex.list, function(x) Seurat::CreateSeuratObject(x))
      gex.metrics.table <- lapply(samples.in, function(x) return(x[[2]][[2]]))

      #clear GEX part in object
      for(i in 1:length(samples.in)){
        samples.in[[i]][[2]] <- "loaded"
      }

      gex.loaded <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded GEX data from Data.in"))

      } else {
        warning("At least on GEX slot in data.in provided list did not contain GEX information. Skipping GEX processing")
        gex.loaded <- F
      }
  } else {
    if(seurat.loaded == F) message("\n GEX not loaded")
  }

  }, error = function(e){
    message("\n Loading GEX failed")
    message(e)}) #end loading from Data.in

  #### Load GEX from Seurat object ####

  if(gex.loaded == F & seurat.loaded == T){
      gex.list <- list(Seurat.in) #changing to list type for compatibility
      Seurat.in <- NULL
      gex.loaded <- F #!
  }

  #### Load Feature barcodes ####

  #! Not checking if FB have already been loaded as part of GEX. If FB directory is supplied, any FB data loaded by GEX is therefore overwritten.
  # This should allow for more flexibility of the user without having to realign the whole GEX data
  if(inherits(samples.in,"character") & FB.out.directory.list[[1]] != "none" & fb.loaded == F){ #No Data.in or input => proceed with loading by BC.out.directory.list
    #Remove possible backslash at the end of the input path
    FB.out.directory.list <- lapply(FB.out.directory.list, function(x) return(gsub("/$", "", x)))

    gex_load_error <- tryCatch({suppressWarnings({
      #add the directory identifier

      if(stringr::str_detect(FB.out.directory.list[[1]], "filtered_feature_bc_matrix")){
        FB.out.directory.list <- FB.out.directory.list
        if(verbose) message("\n ! Feature barcode path was specified explicitely to filtered_feature_bc_matrix. For better rates we recommend using the raw_feature_bc_matrix folder content !")
      } else if (stringr::str_detect(FB.out.directory.list[[1]], "raw_feature_bc_matrix")){
        #Nothing to append
        FB.out.directory.list <- FB.out.directory.list
        if(verbose) message("\n Loading feature barcodes from raw_feature_bc_matrix folder")
      } else {
        if(dir.exists(paste(FB.out.directory.list[[1]],"/raw_feature_bc_matrix", sep = ""))){
          FB.out.directory.list <- paste0(FB.out.directory.list, "/raw_feature_bc_matrix")
          if(verbose) message("\n Loading feature barcodes from raw_feature_bc_matrix folder")
        } else if (dir.exists(paste(FB.out.directory.list[[1]],"/filtered_feature_bc_matrix", sep = ""))){
            FB.out.directory.list <- paste0(FB.out.directory.list, "/filtered_feature_bc_matrix")
            if(verbose) message("\n Loading feature barcodes from filtered_feature_bc_matrix folder")
        } else {
          FB.out.directory.list <- FB.out.directory.list
          if(verbose) message("\n Loading feature barcodes from raw input path")
        }
      }
      #Actually loading the data. Critical: if Cellranger 6.1.0 count was run with --libraries input containing both GEX and Feature Barcodes, reading it will result in a list of matrices instead of a single matrix. => the next section deals with this

      n_not_loaded <- 0
      directory_read10x <- lapply(FB.out.directory.list, function(x) tryCatch({Seurat::Read10X(data.dir=x)},
                                                                                error = function(e){
                                                                                  n_not_loaded <- n_not_loaded + 1
                                                                                  return(NULL)})) #Catching an error which occurs when trying to load the directory "PLACEHOLDER/filtered_feature_bc_matrix" if a placeholder was provided in input directories.
      #List elements of directory_read10x that failed loading (because a placeholder path was provided, will be NULL)
      #If all are null, we retry with a different input format
      if(all(is.null(directory_read10x))){
        warning("None of FB data could be loaded. Retrying with gene.column = 1")
        directory_read10x <- lapply(FB.out.directory.list, function(x) tryCatch({Seurat::Read10X(data.dir=x, gene.column = 1)},
                                                                        error = function(e){
                                                                          n_not_loaded <- n_not_loaded + 1
                                                                          return(NULL)}))
      }
      #If this failes to, we advise to revise input paths
      if(all(is.null(directory_read10x))) stop("FB data loading failed from provided paths. Please revise input paths")

      #Assuming that all FB from directories have been loaded correctly and that the NULL list items are "PLACEHOLDER" inputs
      #Replacing null values in the loaded list with placeholder matrices which can be processed without warnings and if conditions
      for(i in 1:length(directory_read10x)){
        if(is.null(directory_read10x[[i]])){
          directory_read10x[[i]] <- Matrix::Matrix(data = c(rep(1001, 10), rep(1, 10)), nrow = 2, ncol = 10, sparse = TRUE, dimnames = list( c("No-FB-data", "column2"),LETTERS[11:20]))
        }
      }

      print(length(directory_read10x))
      #dealing with possible mixed GEX FB inputs or multiple FB input matrices from the same directory
      for(i in 1:length(directory_read10x)){ #iterating over main list
        if(inherits(directory_read10x[[i]],"list")){ #this returns true only if the FB directory import contained more than one marices
          if(verbose) message(paste0("\n Feature barcode input ", i, " contains multiple count matrices."))
          to_del <- c()
          for(j in 1:length(directory_read10x[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(directory_read10x[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              if(verbose) message(paste0("Feature barcode input ", i, " element ", j, " contains > 100 features and likely corresponds to GEX data. This matrix will be removed from further FB processing"))
              to_del <- c(to_del, j) #Will be deleted later to not mess up the loop
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("\n GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }
          if(length(to_del) > 0){
            directory_read10x[[i]] <- directory_read10x[[i]][-to_del] #deleting
          }
        }
      }

      #now flatten the remaining list
      directory_read10x <- do.call(list, unlist(directory_read10x, recursive=FALSE))
      warning(paste0("At least one Feature barcode input contains multiple count matrices. Count matrices containing more than 100 features were filtered out from FB assignment as they most likely correspond to GEX data"))
      #Done => result should be a non-nested list of matrices only containing FB information.

      #We can now convert all to Seurat objects
      FB.list <- lapply(directory_read10x, function(x) Seurat::CreateSeuratObject(x))
      directory_read10x <- NULL

      fb.loaded <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded FB data"))

    })}, error = function(e){
      message("\n Loading FB failed")
      message(e)})

  } else if(inherits(samples.in,"list") & fb.loaded == F){ #GET info from samples.in

    FB.list <- lapply(samples.in, function(x) return(x[[3]][[1]]))
    FB.list.class <- lapply(samples.in, function(x) return(class(x[[3]][[1]])))

    if(!all(FB.list.class == "character")){ #in case that FB data was loaded for at least one sample
      #Replacing none values in the loaded list with placeholder matrices which can be processed without warnings and if conditions
      for(i in 1:length(FB.list)){
        if(inherits(FB.list[[i]],"character")){
          FB.list[[i]] <- Matrix::Matrix(data = c(rep(1001, 10), rep(1, 10)), nrow = 2, ncol = 10, sparse = TRUE, dimnames = list(c("No-FB-data", "column2"),LETTERS[11:20]))
        }
      }

      #Make Seurat
      FB.list <- lapply(FB.list, function(x) Seurat::CreateSeuratObject(x))
      fb.loaded <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Loaded FB data from data.in"))

    } else if(fb.loaded == F){
      FB.list <- "none"
    }
  } else if(fb.loaded == F){
    FB.list <- "none"
  }


  #### Call VDJ stats ####

  stats.done <- F
  if(get.VDJ.stats == T & vdj.loaded == T){
    if(verbose) message("\n Getting VDJ GEX stats")
    tryCatch({
      out.stats <- VDJ_GEX_stats_int(clonotype.list = clonotype.list,
                                     reference.list = reference.list,
                                     annotations.list = annotations.list,
                                     contig.table = contig.table,
                                     vdj.metrics = metrics.table,
                                     gex.metrics = gex.metrics.table,
                                     samples.paths.VDJ = stringr::str_split(samples.paths.VDJ, " ; ", simplify = T)[1,], #needing to split sample paths again as they will end up in separate rows in the stats dataframe
                                     verbose = verbose)
      stats.done <- T
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Got VDJ GEX stats"))

    }, error = function(e){e
      message(paste0("VDJ stats failed: ", e, "\n"))
      message()})
  }

  #### Select barcodes for processing ####

  #Proceed with selecting the barcodes which will be processed
  #Any barcode that is a cell in the VDJ or in the GEX it will be processed
  #because the filtered GEX is loaded, all barcodes from the GEX will be included in any case.
  #If a barcode is a cell in VDJ but is not present in the filtered GEX, it will be analyzed for VDJ irregardless

  barcodes_GEX <- list()
  barcodes_VDJ <- list()
  barcodes_GEX_raw <- list()
  barcodes_VDJ_raw <- list()

  if(gex.loaded == T & vdj.loaded == T & seurat.loaded == F){ #If both VDJ and GEX are available

    for(i in 1:length(gex.list)){
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])
      #raw barcodes GEX
      barcodes_GEX_raw[[i]] <- gsub("(^_)|(-\\d+.*$)","",barcodes_GEX[[i]])
      gex.list[[i]]@meta.data$orig_barcode <- barcodes_GEX_raw[[i]] #Add as column

      #VDJ filtering to make sure => with the new input strategy post Cellranger 6.1. this should not have any effect because of prefiltering by cellranger
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true"
                                                                  & tolower(contig.table[[i]]$high_confidence) == "true"
                                                                  & tolower(contig.table[[i]]$productive) == "true"
                                                                  & tolower(contig.table[[i]]$full_length) == "true")])
      #raw barcodes VDJ
      barcodes_VDJ_raw[[i]] <- gsub("(^_)|(-\\d+.*$)","",barcodes_VDJ[[i]])

      if(verbose) message(paste0("\n For sample ", i, ": ", length(barcodes_GEX_raw[[i]])," cell assigned barcodes in GEX, ", length(barcodes_VDJ_raw[[i]]), " cell assigned high confidence barcodes in VDJ. Overlap: ", sum(barcodes_GEX_raw[[i]] %in% barcodes_VDJ_raw[[i]])))

      vdj.gex.available <- barcodes_GEX_raw[[i]] %in% barcodes_VDJ_raw[[i]]
      gex.list[[i]] <- SeuratObject::AddMetaData(gex.list[[i]], vdj.gex.available, col.name = "VDJ_available")

      #remove all barcodes in GEX which are not present in VDJ (defaults to FALSE)
      if(exclude.GEX.not.in.VDJ == T){
        gex.list[[i]] <- subset(gex.list[[i]], cells = colnames(gex.list[[i]])[which(gex.list[[i]]$VDJ_available == T)])
        if(verbose) message(paste0("Removed ", length(vdj.gex.available)-sum(vdj.gex.available), " GEX cells that were not present in VDJ"))
      }
    }
  }

  if(gex.loaded == F & vdj.loaded == T & seurat.loaded == T){ #If VDJ is available and GEX is loaded from Seurat.in
    #Pulling single vectors from Seurat.in
    barcodes_GEX <- colnames(gex.list[[1]])
    barcodes_GEX_raw <- gsub("(^_)|(-\\d+.*$)","",barcodes_GEX)
    gex.list[[1]]@meta.data$orig_barcode <- barcodes_GEX_raw #Add as column

    barcodes_VDJ <- list()
    barcodes_VDJ_raw <- list()
    for(i in 1:length(contig.table)){##barcodes_VDJ holds the unique barcodes
      #VDJ filtering to make sure => with the new input strategy post Cellranger 6.1. this should not have any effect because of prefiltering by cellranger
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true"
                                                                  & tolower(contig.table[[i]]$high_confidence) == "true"
                                                                  & tolower(contig.table[[i]]$productive) == "true"
                                                                  & tolower(contig.table[[i]]$full_length) == "true")])
      #raw barcodes VDJ
      barcodes_VDJ_raw[[i]] <- sub("(^_)|(-\\d+.*$)","",barcodes_VDJ[[i]])
    }
    barcodes_VDJ_raw <- unlist(barcodes_VDJ_raw) #reducing to a single vector

    if(verbose) message(paste0("\n For input Seurat object: ", length(barcodes_GEX)," cells assigned barcodes in GEX"))

    vdj.gex.available <- barcodes_GEX_raw %in% barcodes_VDJ_raw
    if(verbose) message(paste0("For input Seurat object GEX and VDJ barcode overlap is: ", sum(vdj.gex.available)))
    gex.list[[1]] <- SeuratObject::AddMetaData(gex.list[[1]], vdj.gex.available, col.name = "VDJ_available")

    #remove all barcodes in GEX which are not present in VDJ (defaults to FALSE)
    if(exclude.GEX.not.in.VDJ == T){
      gex.list[[1]] <- subset(gex.list[[1]], cells = colnames(gex.list[[1]])[which(gex.list[[1]]$VDJ_available == T)])
      if(verbose) message(paste0("Removed ", length(vdj.gex.available)-sum(vdj.gex.available), " GEX cells that were not present in VDJ"))
      }
    }

  #If only GEX is processed
  if(gex.loaded == T & vdj.loaded == F){
    for(i in 1:length(gex.list)){
      barcodes_GEX[[i]] <- colnames(gex.list[[i]])
      barcodes_GEX_raw[[i]] <- gsub("(^_)|(-\\d+.*$)","",barcodes_GEX[[i]])
      gex.list[[i]]@meta.data$orig_barcode <- barcodes_GEX_raw[[i]] #Add as column
      if(verbose) message(paste0("\nFor sample ", i, ": ", length(barcodes_GEX[[i]])," cell assigned barcodes in GEX"))
    }
  }
  #If only VDJ is processed
  if(gex.loaded == F & vdj.loaded == T){
    for(i in 1:length(contig.table)){##barcodes_VDJ holds the unique barcodes
      barcodes_VDJ[[i]] <- unique(contig.table[[i]]$barcode[which(tolower(contig.table[[i]]$is_cell) == "true"
                                                                  & tolower(contig.table[[i]]$high_confidence) == "true"
                                                                  & tolower(contig.table[[i]]$productive) == "true"
                                                                  & tolower(contig.table[[i]]$full_length) == "true")])
      if(verbose) message(paste0("\nFor sample ", i, ": ", length(barcodes_VDJ[[i]]), " cells assigned with high confidence barcodes in VDJ"))
    }
  }

  #### Barcode overlap removal ####

  #GEX
  if(filter.overlapping.barcodes.GEX == T & gex.loaded == T){
    if(length(gex.list) > 1){
      barcodes_GEX_c <- do.call("c", lapply(gex.list, function(x) x$orig_barcode))
      unique_barcodes <- names(table(barcodes_GEX_c)[table(barcodes_GEX_c) == 1])
      for(i in 1:length(gex.list)){
        if(sum(gex.list[[i]]$orig_barcode %in% unique_barcodes) == 0){
          stop("\n Removal of non-unique barcodes in GEX has resulted in 0 cells in at least one GEX sample. Please check that no GEX input paths are duplicated")}
        gex.list[[i]] <- subset(gex.list[[i]], subset = orig_barcode %in% unique_barcodes)
      }
      if(verbose) message(paste0("\n Removed a total of ", length(unique(barcodes_GEX_c)) - length(unique_barcodes), " cells with non unique barcodes in GEX"))
    }
  }

  #VDJ
  if(filter.overlapping.barcodes.VDJ == T & vdj.loaded == T){
    if(length(barcodes_VDJ) > 1){
      barcodes_VDJ_c <- do.call("c", barcodes_VDJ)
      non_unique_barcodes <- names(table(barcodes_VDJ_c)[table(barcodes_VDJ_c) > 1])
      for(i in 1:length(barcodes_VDJ)){
        barcodes_VDJ[[i]] <- barcodes_VDJ[[i]][which(!barcodes_VDJ[[i]] %in% non_unique_barcodes)]
        if(length(barcodes_VDJ[[i]]) == 0){
          stop("Removal of non-unique barcodes in VDJ has resulted in 0 cells in at least one VDJ sample. Please check that no VDJ input paths are duplicated")}
      }
      if(verbose) message(paste0("\n Removed a total of ", length(non_unique_barcodes), " cells with non unique barcodes in VDJ"))
    }
  }

  if(seurat.loaded & filter.overlapping.barcodes.GEX){warning("Filtering of overlapping barcodes in GEX is not performed for Seurat.in input")}

  #### Filtering by input barcode list ####

  if(exclude.on.barcodes[1] != "none" & (gex.loaded == T | seurat.loaded == T)){
      barcodes_GEX_c <- do.call("c", lapply(gex.list, function(x) x$orig_barcode))
      to_keep <- barcodes_GEX_c[!barcodes_GEX_c %in% exclude.on.barcodes]
      for(i in 1:length(gex.list)){
        n_cells <- length(gex.list[[i]]$orig_barcode)
        gex.list[[i]] <- subset(gex.list[[i]], subset = orig_barcode %in% to_keep)
        if(verbose) message(paste0("\n In GEX input ", i ," removed a total of ", n_cells - length(gex.list[[i]]$orig_barcode), " cells with barcodes provided as exclude.on.barcodes"))
      }
    }

  if(exclude.on.barcodes[1] != "none" & vdj.loaded == T){
      for(i in 1:length(barcodes_VDJ)){
        n_cells <- length(barcodes_VDJ[[i]])
        barcodes_VDJ[[i]] <- barcodes_VDJ[[i]][which(!gsub("-\\d+.*$","",barcodes_VDJ[[i]]) %in% exclude.on.barcodes)]
        if(verbose) message(paste0("\n In VDJ sample ", i ," removed a total of ", n_cells - length(barcodes_VDJ[[i]]), " cells with barcodes provided as exclude.on.barcodes"))
      }
    }

  #### Exclude cells on cell state markers ####

  #handlers copied from GEX_phenotype, Thanks Alex!
  if(exclude.on.cell.state.markers[1] != "none" & gex.loaded == T){

    Cap<-function(x){
      temp<-c()
      for (i in 1:length(x)){
        s <- strsplit(x, ";")[[i]]
        temp[i]<-paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)), sep="", collapse=";")
      }
      return(temp)
    }

    #rename to match GEX_phenotype variables
    cell.state.markers <- exclude.on.cell.state.markers
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
            gex.list[[j]]@meta.data <- gex.list[[j]]@meta.data[-c(ncol(gex.list[[j]]@meta.data))] #remove that column again
            if(verbose) message(paste0("\n In GEX sample ", j ," excluded ", cells_unfiltered - ncol(gex.list[[j]])," cells based on ", cell.state.markers[i]))
            #If the Gene was not found in seurat object features
          } else{
            if(verbose) message(paste0("\n In GEX sample ", j ," failed to exclude cells based on: ", cell.state.markers[i], " Please check gene spelling"))
          }
        }
      }
      #larger error callback should be independent of user input.
    }, error = function(e){
      message("\n Exclusion based on cell markers failed")
      message(e)
    })
  }

  #### Subsample VDJ barcodes ####

  if(subsample.barcodes == T & vdj.loaded == T){
    #For development => shorten barcode list to shorten computational time during development
    if(verbose) message("Sampling 50 barcodes from all in VDJ per sample \n")
    for(i in 1:length(barcodes_VDJ)){
      barcodes_VDJ[[i]] <- sample(barcodes_VDJ[[i]],50)
    }
  }

  #### VDJ Processing per cell ####

  if(vdj.loaded == T){
    VDJ.proc.list <- list()
    for(i in 1:length(contig.table)){

      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " "," Starting VDJ barcode iteration ", i , " of ", length(contig.table), "..."))

      if(parallel.processing == "parlapply" | parallel.processing == "parLapply"){
        doParallel::stopImplicitCluster() #close any open clusters
        #open cluster for parallel computing
        cl <- parallel::makeCluster(numcores)
        if(verbose) message(paste0("Started parlapply cluster with ", numcores, " cores"))
        out.VDJ <- parallel::parLapply(cl, barcodes_VDJ[[i]],
                                       barcode_VDJ_iteration,
                                       contigs = contig.table[[i]],
                                       clonotypes = clonotype.list[[i]],
                                       references = reference.list[[i]],
                                       annotations = annotations.table[[i]],
                                       gap.extension.cost = gap.extension.cost,
                                       gap.opening.cost = gap.opening.cost,
                                       trim.and.align = trim.and.align,
                                       select.excess.chains.by.umi.count = select.excess.chains.by.umi.count,
                                       excess.chain.confidence.count.threshold = excess.chain.confidence.count.threshold)
        #close any open clusters
        doParallel::stopImplicitCluster()
      } else if(parallel.processing == "mclapply"){
        if(verbose) message(paste0("Started mcapply cluster with ", numcores, " cores"))
        out.VDJ <- parallel::mclapply(X = barcodes_VDJ[[i]],
                                      FUN = barcode_VDJ_iteration,
                                      contigs = contig.table[[i]],
                                      clonotypes = clonotype.list[[i]],
                                      references = reference.list[[i]],
                                      annotations = annotations.table[[i]],
                                      gap.extension.cost = gap.extension.cost,
                                      gap.opening.cost = gap.opening.cost,
                                      trim.and.align = trim.and.align,
                                      select.excess.chains.by.umi.count = select.excess.chains.by.umi.count,
                                      excess.chain.confidence.count.threshold = excess.chain.confidence.count.threshold)
      } else { #No parallel computing
        out.VDJ <- lapply(barcodes_VDJ[[i]], barcode_VDJ_iteration,
                          contigs = contig.table[[i]],
                          clonotypes = clonotype.list[[i]],
                          references = reference.list[[i]],
                          annotations = annotations.table[[i]],
                          gap.extension.cost = gap.extension.cost,
                          gap.opening.cost = gap.opening.cost,
                          trim.and.align = trim.and.align,
                          select.excess.chains.by.umi.count = select.excess.chains.by.umi.count,
                          excess.chain.confidence.count.threshold = excess.chain.confidence.count.threshold)
      }

      #bind list recieved from parLapply
      VDJ.proc <- dplyr::bind_rows(out.VDJ)
      VDJ.proc[VDJ.proc == ";"] <- "" #fix bug, where if two emtpy strings are concatenated, a ";" is left behind.
      VDJ.proc[is.na(VDJ.proc)] <- "" #Replace NA (empty values) with an empty string for format compatibility

      #update barcodes
      VDJ.proc$orig_barcode <- gsub("(^_)|(-\\d+.*$)","",VDJ.proc$barcode)
      VDJ.proc$barcode <- paste0("s",i,"_",VDJ.proc$orig_barcode)
      rownames(VDJ.proc) <- VDJ.proc$barcode
      VDJ.proc$sample_id <- paste0("s",i)
      VDJ.proc$group_id <- group.id[i]

      #update celltypes
      VDJ.proc$celltype[stringr::str_detect(paste0(VDJ.proc$VDJ_chain,VDJ.proc$VJ_chain), "TR")] <- "T cell"
      VDJ.proc$celltype[stringr::str_detect(paste0(VDJ.proc$VDJ_chain,VDJ.proc$VJ_chain), "IG")] <- "B cell"

      #Fill extra clonotype_id column with 10x default
      VDJ.proc$clonotype_id <- VDJ.proc$clonotype_id_10x

      #Add further columns to fill in in future updates
      VDJ.proc$specifity <- NA
      VDJ.proc$affinity <- NA

      VDJ.proc.list[[i]] <- VDJ.proc

      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Done with ", i , " of ", length(contig.table)))
    }
    VDJ.proc <- VDJ.proc.list
    VDJ.proc.list <- NULL
    gc()

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

  #### GEX Processing ####
  if(gex.loaded == T & seurat.loaded == F){ #Only execute this processing if the GEX object is not already processed
    if(verbose) message("\n Starting GEX pipeline")
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
                                    group.id = group.id,
                                    verbose = verbose)

    if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Done with GEX pipeline"))

  } else if (gex.loaded == F & seurat.loaded == F){
    GEX.proc <- "none"
  } else if (gex.loaded == F & seurat.loaded == T){

    if(verbose) message("\n Preparing Seurat.in object")

    GEX.proc <- gex.list
    #going through samples and renaming cells by that
    if(stringr::str_detect(colnames(GEX.proc[[1]])[1], "^s\\d_")){
      #input seurat object was already processed via platypus. not renaming
    } else {
    GEX.proc[[1]] <- SeuratObject::RenameCells(GEX.proc[[1]], new.names = paste0(GEX.proc[[1]]$sample_id, "_",gsub("(^_)|(-\\d+.*$)","",colnames(GEX.proc[[1]]))))
    }

    if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Done with Seurat.in"))
  }

  #### Feature Barcode preparation ####

  FB.processed <- F
  if(fb.loaded){
    if(verbose) message("\n Starting FB processing")
    tryCatch({
      #getting relevant info as a dataframe
      FB.list <- lapply(FB.list, function(x) return(SeuratObject::FetchData(x, rownames(x))))
      #filtering columns of the table based on pattern input
      for(i in 1:length(FB.list)){
        FB.list[[i]] <- FB.list[[i]][, stringr::str_detect(names(FB.list[[i]]), FB.exclude.pattern) == F]
      }
      #assigning FBs using also the user tunable threshold as an additional argument
      FB.list <- lapply(FB.list, pick_max_feature_barcode, FB.ratio.threshold, FB.count.threshold)
      #merge FB list but append sample id to barcode

      for(i in 1:length(FB.list)){
        if(stringr::str_detect(FB.list[[i]][1,1],"K")){ #This is LETTERS[11] which used as rownames for the place holder matrix
          if(verbose) message(paste0("\n For GEX/VDJ input ", i, " no FB data was loaded. In output FB_assignment column cells of this sampled are labelled 'Not assignable'"))
        }
      }

      FB.processed <- T #keeping track for VDJ
      if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Got FB assignment to cell barcodes"))
    }, error = function(e){
      message("Processing FB failed")
      message(e)
    })
  }

  #### Feature Barcode assignment to GEX ####
  #Here we merge feature barcodes into the GEX matrices.
  #Because users may choose to not combine individual samples, we need to deal with multiple seurat objects in the GEX.proc and VDJ.proc lists
  if(fb.loaded == T & inherits(GEX.proc,"list")){ #ensuring correct input

    if(verbose) message("\n Adding Feature barcode information to GEX")

    tryCatch({

      #make a copy of FB.list for VDJ below
      FB.list.copy.for.vdj <- FB.list

      #process FB further
      if(length(GEX.proc) == 1){ #Either only one sample or GEX.integrate == T
        if(length(unique(GEX.proc[[1]]$sample_id)) != length(FB.list)){ #no 1-1 corespondence between lengths...
          #This should not happen as every GEX library is associated with on FB library
          stop("\n FB assignment failed because number of FB input matrices does not match number of GEX input matrices")
        }

        #Format barcodes in FB and make dataframe
        for(i in 1:length(FB.list)){FB.list[[i]]$barcode <- paste0("s",i,"_", FB.list[[i]]$orig_barcode)}
        FB.list <- as.data.frame(dplyr::bind_rows(FB.list))

        #GET from GEX
        meta_to_merge <- SeuratObject::FetchData(GEX.proc[[1]], "orig_barcode") #getting a reference to merge into
        #now making the barcodes matching with the format in FB.list
        meta_to_merge$barcode <- paste0(stringr::str_split(rownames(meta_to_merge), "_", simplify = T)[,1], "_",meta_to_merge$orig_barcode)

        meta_to_merge <- merge(meta_to_merge, FB.list, by = "barcode", all.x = T, all.y = F, sort = F) #merging making sure to not add or remove any rows
        rownames(meta_to_merge) <- rownames(GEX.proc[[1]]@meta.data) #reconstitute the rownames.
        #make sure that there are no NAs
        meta_to_merge$FB_assignment[is.na(meta_to_merge$FB_assignment)] <- "Not assignable"
        GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], meta_to_merge[,"FB_assignment"], col.name = "FB_assignment") #add to object.
        #move the column to after sample_id as this column will probably be used a lot
        sample_id_index <- which(names(GEX.proc[[1]]@meta.data) == "sample_id")
        GEX.proc[[1]]@meta.data <- GEX.proc[[1]]@meta.data[,c(1:sample_id_index, ncol(GEX.proc[[1]]@meta.data), (sample_id_index+1):(ncol(GEX.proc[[1]]@meta.data)-1))]

      } else if(length(GEX.proc) > 1){ #With multiple GEX samples and GEX.integrate = FALSE
        if(length(GEX.proc) == length(FB.list)){ #if there is 1-1 correspondence

          for(i in 1:length(FB.list)){ #iterating over FB tables
            meta_to_merge <- SeuratObject::FetchData(GEX.proc[[i]], "orig_barcode") #getting a reference to merge into
            meta_to_merge <- merge(meta_to_merge, FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F, sort = F) #merging making sure to not add or remove any rows
            rownames(meta_to_merge) <- rownames(GEX.proc[[i]]@meta.data)#reconstitute the rownames.
            #make sure that there are no NAs
            meta_to_merge$FB_assignment[is.na(meta_to_merge$FB_assignment)] <- "Not assignable"
            GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], meta_to_merge[,"FB_assignment"], col.name = "FB_assignment") #add to object
            #move the column to after sample_id as this column will probably be used a lot
            sample_id_index <- which(names(GEX.proc[[i]]@meta.data) == "sample_id")
            GEX.proc[[i]]@meta.data <- GEX.proc[[i]]@meta.data[,c(1:sample_id_index, ncol(GEX.proc[[i]]@meta.data), (sample_id_index+1):(ncol(GEX.proc[[i]]@meta.data)-1))]
          }
        }else { #no 1-1 correspondence
          stop("\n FB assignment failed because number of FB input matrices does not match number of GEX input matrices")
        }
      }
      message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Added Feature barcode info to GEX"))
      meta_to_merge <- NULL
    }, error = function(e){
      message("Adding Feature barcode information to GEX failed")
      message(e)
    })
  } else {
    FB.processed <- F #keeping track for VDJ
  }

  #### Feature Barcode assignment to VDJ ####

  if(fb.loaded == T & (inherits(VDJ.proc,"list") | inherits(VDJ.proc,"data.frame"))){ #ensuring correct input

    if(verbose) message("\n Adding Feature barcode information to VDJ")
    tryCatch({
      if(inherits(VDJ.proc,"data.frame")){ #Either only one sample or VDJ.combine == T

        print(class(FB.list))
        if(inherits(FB.list, "data.frame")) FB.list <- FB.list.copy.for.vdj #reconstitute

        print(VDJ.proc$sample_id)
        if(length(unique(VDJ.proc$sample_id)) != length(FB.list)){ #no 1-1 corespondence between lengths...
          #This should not happen as every VDJ library is associated with on FB library
          stop("\n FB assignment failed because number of FB input matrices does not match number of VDJ input matrices")
        }

        #Format barcodes in FB and make dataframe
        for(i in 1:length(FB.list)){FB.list[[i]]$barcode <- paste0("s",i, "_", FB.list[[i]]$orig_barcode)}
        FB.list <- as.data.frame(dplyr::bind_rows(FB.list))

        #add the sample identifier to the orig_barcode column in VDJ similar to GEX
        VDJ.proc$barcode <- paste0(stringr::str_split(VDJ.proc$barcode, "_", simplify = T)[,1], "_",VDJ.proc$orig_barcode)
        VDJ.proc <- merge(VDJ.proc, FB.list, by = "barcode", all.x = T, all.y = F) #merging making sure to not add or remove any rows
        #make sure that there are no NAs
        VDJ.proc$FB_assignment[is.na(VDJ.proc$FB_assignment)] <- "Not assignable"
        #move the column to after sample_id as this column will probably be used a lot
        sample_id_index <- which(names(VDJ.proc) == "sample_id")
        VDJ.proc<- VDJ.proc[,c(1:sample_id_index, ncol(VDJ.proc), (sample_id_index+1):(ncol(VDJ.proc)-1))]

      } else if(inherits(VDJ.proc,"list")){
        if(length(VDJ.proc) == length(FB.list)){ #if there is 1-1 correspondence
          for(i in 1:length(FB.list)){ #iterating over FB tables
            VDJ.proc[[i]] <- merge(VDJ.proc[[i]], FB.list[[i]], by = "orig_barcode", all.x = T, all.y = F)
            #make sure that there are no NAs
            VDJ.proc[[i]]$FB_assignment[is.na(VDJ.proc[[i]]$FB_assignment)] <- "Not assignable"
            #move the column to after sample_id as this column will probably be used a lot
            sample_id_index <- which(names(VDJ.proc[[i]]) == "sample_id")
            VDJ.proc[[i]]<- VDJ.proc[[i]][,c(1:sample_id_index, ncol(VDJ.proc[[i]]), (sample_id_index+1):(ncol(VDJ.proc[[i]])-1))]
          }
        } else { #Should not happen
          stop("\n FB assignment failed because number of FB input matrices does not match number of VDJ input matrices")
        }
      }
      message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Added Feature barcode info to VDJ"))
      FB.list <- NULL
      FB.list.merge <- NULL
    }, error = function(e){
      message("Adding Feature barcode information to VDJ failed")
      message(e)
    })
  }

  #### VDJ GEX integration ####

  if(vdj.loaded == T & (gex.loaded == T | seurat.loaded == T)){
    if(verbose) message("\n Integrating VDJ and GEX ")
    #Combine GEX and VDJ
    #Multiple cases:
    if(inherits(VDJ.proc,"list") & length(GEX.proc) == 1){ #Multiple VDJ, one GEX

      for(i in 1:length(VDJ.proc)){ VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[1]])}

      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]],
                                                 colnames(GEX.proc[[1]]) %in% do.call(c ,lapply(VDJ.proc, function(x) x$barcode))
                                                 , col.name = "VDJ_available")

      #add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){

        #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
        seur_meta <- SeuratObject::FetchData(GEX.proc[[1]],
                                             vars = c("orig.ident","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
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

    } else if (inherits(VDJ.proc, "data.frame") & length(GEX.proc) == 1){ #one VDJ one GEX

      GEX.proc[[1]] <- SeuratObject::AddMetaData(GEX.proc[[1]], colnames(GEX.proc[[1]]) %in% VDJ.proc$barcode, col.name = "VDJ_available")
      VDJ.proc$GEX_available <- VDJ.proc$barcode %in% colnames(GEX.proc[[1]])
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]

      #add some GEX columns to VDJ table
      if(integrate.GEX.to.VDJ == T){
        #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
        seur_meta <- SeuratObject::FetchData(GEX.proc,
                                             vars = c("orig.ident","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
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

    } else if(inherits(VDJ.proc,"list") & length(GEX.proc) > 1){ #multiple VDJ multiple GEX
      for(i in 1:length(GEX.proc)){
        GEX.proc[[i]] <- SeuratObject::AddMetaData(GEX.proc[[i]], colnames(GEX.proc[[i]]) %in% VDJ.proc[[i]]$barcode, col.name = "VDJ_available")
        VDJ.proc[[i]]$GEX_available <- VDJ.proc[[i]]$barcode %in% colnames(GEX.proc[[i]])

        #add some GEX columns to VDJ table
        if(integrate.GEX.to.VDJ == T){
          if(verbose) message("Integrating multiple VDJ and GEX pairwise. ! VDJ and GEX input paths or data must be in the same order!")
          #get data from Seurat object to add to VDJ. In the future this could become an extra parameter
          seur_meta <- SeuratObject::FetchData(GEX.proc[[i]],
                                               vars = c("orig.ident","seurat_clusters","PC_1", "PC_2", "UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2"))
          seur_meta$barcode <- rownames(seur_meta)
          #merge to VDJ.proc
          VDJ.proc[[i]] <- merge(VDJ.proc[[i]], seur_meta, by = "barcode", all.x = T, all.y = F)
        }
        #add VDJ info to GEX as metadata columns
        if(integrate.VDJ.to.GEX == T){
          if(verbose) message("\n Integrating multiple VDJ and GEX pairwise. ! VDJ and GEX input paths or data must be in the same order!")
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

  } else if(vdj.loaded == T & gex.loaded == F & seurat.loaded == F){
    out.list <- list("VDJ" = VDJ.proc, "GEX" = "none")
  } else if(vdj.loaded == F & (gex.loaded == T | seurat.loaded == T)){
    if(length(GEX.proc) == 1){
      #Reduce GEX.proc to a seurat object
      GEX.proc <- GEX.proc[[1]]
    }
    out.list <- list("VDJ" = "none", "GEX" = GEX.proc)
  } else if(vdj.loaded == F & gex.loaded == F & Seurat.loaded == F){
    stop("Neither VDJ or GEX data loaded. Exiting")
  }

  #### Save run time parameters ####

  params <- c("sample.path.vdj" = paste0(samples.paths.VDJ, collapse = " / "),
              "samples.paths.GEX" = paste0(samples.paths.GEX, collapse = " / "),
              "FB.out.directory.list" = paste0(FB.out.directory.list, collapse = " / "),
              "GEX.read.h5" = GEX.read.h5,
              "VDJ.combine" = VDJ.combine,
              "GEX.integrate" = GEX.integrate,
              "integrate.GEX.to.VDJ" = integrate.GEX.to.VDJ,
              "integrate.VDJ.to.GEX" = integrate.VDJ.to.GEX,
              "exclude.GEX.not.in.VDJ" = exclude.GEX.not.in.VDJ,
              "filter.overlapping.barcodes.GEX" = filter.overlapping.barcodes.GEX,
              "filter.overlapping.barcodes.VDJ" = filter.overlapping.barcodes.VDJ,
              "exclude.on.cell.state.markers" = paste0(exclude.on.cell.state.markers,collapse = ";"),
              "exclude.on.barcodes (TRUE if barcodes provided)" = (exclude.on.barcodes[1] == "none"),
              "get.VDJ.stats" =  get.VDJ.stats,
              "numcores" = numcores,
              "trim.and.align" = trim.and.align,
              "append.raw.reference" = append.raw.reference,
              "select.excess.chains.by.umi.count" = select.excess.chains.by.umi.count,
              "excess.chain.confidence.count.threshold" = excess.chain.confidence.count.threshold,
              "gap.opening.cost," = gap.opening.cost,
              "gap.extension.cost" = gap.extension.cost,
              "parallel.processing" = parallel.processing,
              "integration.method" = integration.method,
              "VDJ.gene.filter" = VDJ.gene.filter,
              "mito.filter" = mito.filter,
              "norm.scale.factor" = norm.scale.factor,
              "n.feature.rna" = n.feature.rna,
              "n.count.rna.min" = n.count.rna.min,
              "n.count.rna.max" = n.count.rna.max,
              "n.variable.features" = n.variable.features,
              "cluster.resolution" = cluster.resolution,
              "neighbor.dim" = paste0(neighbor.dim, collapse = ";"),
              "mds.dim" = paste0(mds.dim, collapse = ";"),
              "subsample.barcodes" = subsample.barcodes,
              "group.id" = paste0(group.id, collapse = ";"),
              "FB.count.threshold" = FB.count.threshold,
              "FB.ratio.threshold" = FB.ratio.threshold)

  #### Compile output ####

  if(batches[[1]] != "none"){
    if(inherits(out.list[[1]],"data.frame")){ #check that vdj has been completed
      if(length(batches) == length(unique(out.list[[1]]$sample_id))){
        out.list[[1]]$batch_id <- as.character(out.list[[1]]$sample_id)
        for (i in 1:length(unique(out.list[[1]]$batch_id ))){
          out.list[[1]]$batch_id <- gsub(unique(out.list[[1]]$sample_id)[i], batches[i], out.list[[1]]$batch_id)}
      }
    }
    if(!inherits(out.list[[2]],"character") & GEX.integrate == T){ #check that gex has been completed and if there is only a single GEX object
      if(length(batches) == length(unique(out.list[[2]]$sample_id))){
        out.list[[2]]$batch_id <- as.character(out.list[[2]]$sample_id)
        for (i in 1:length(unique(out.list[[2]]$batch_id))){
          out.list[[2]]$batch_id <- gsub(unique(out.list[[2]]$sample_id)[i], batches[i], out.list[[2]]$batch_id)}
      }
    }
  } else if(batches[[1]] == "none"){#append the column anyways to keep formatting consistent
    if(inherits(out.list[[1]],"data.frame")){ #check that vdj has been completed
      out.list[[1]]$batches <- "Unspecified"
    }
    if(!inherits(out.list[[2]],"character") & GEX.integrate == T){ #check that gex has been completed
      out.list[[2]]$batches <- "Unspecified"
    }
  }

  #adding VDJ stats
  if(stats.done == T){
    if(verbose) message("\n Adding VDJ stats")
    out.list[[3]] <- out.stats
  } else {
    out.list[[3]] <- "VDJ stats failed"
  }

  #adding parameter info for reproducibility
  if(verbose) message(" Adding runtime params")
  out.list[[4]] <- params

  #finally add session info
  out.list[[5]] <- utils::sessionInfo()

  #rename for clarity
  names(out.list) <- c("VDJ", "GEX", "VDJ.GEX.stats", "Running params", "sessionInfo")

  if(verbose) message(paste0(stringr::str_split(Sys.time(), " ", simplify = T)[,2], " ","Done!"))
  if(verbose) message(Sys.time())

  return(out.list)
}
