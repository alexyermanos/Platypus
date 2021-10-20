#' Utility function for loading in local dataset as VDJ_GEX_matrix and PlatypusDB compatible R objects. Especially useful when wanting to integrate local and public datasets. This function only imports and does not make changes to format, row and column names. Exception: filtered_contig.fasta are appended to the filtered_contig_annotations.csv as a column for easy access
#' @param VDJ.out.directory.list List containing paths to VDJ output directories from cell ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#'@param GEX.out.directory.list List containing paths the outs/ directory of each sample or directly the raw or filtered_feature_bc_matrix folder. Order of list items must be the same as for VDJ. This outs directory may also contain Feature Barcode (FB) information. Do not specify FB.out.directory in this case.
#'@param FB.out.directory.list List of paths pointing at the outs/ directory of output of the Cellranger counts function which contain Feature barcode counts. Any input will overwrite potential FB data loaded from the GEX input directories. Length must match VDJ and GEX directory inputs. (in case of a single FB output directory for multiple samples, please specifiy this directory as many times as needed)
#'@param batches Integer vector. Defaults to all 1, yielding all samples with batch number "b1". Give a batch number to each sample (each entry in the VDJ/GEX input lists). This will be saved as element 5 in the sample list output.
#' @return Large nested list object containing all needed Cellranger outputs to run the VDJ_GEX_matrix function. Level 1 of the list are samples, level 2 are VDJ GEX and metadata information. (e.g. out[[1]][[1]] corresponds to VDJ data objects of sample 1)
#' @export
#' @examples
#' \dontrun{
#' VDJ.in <- list()
#' VDJ.in[[1]] <- c("~/VDJ/S1/")
#' VDJ.in[[2]] <- c("~/VDJ/S2/")
#' GEX.in <- list()
#' GEX.in[[1]] <- c("~/GEX/S1/")
#' GEX.in[[2]] <- c("~/GEX/S2/")
#' PlatypusDB_load_from_disk(VDJ.out.directory.list = VDJ.in, GEX.out.directory.list = GEX.in)
#' }
PlatypusDB_load_from_disk <- function(VDJ.out.directory.list,
                                      GEX.out.directory.list,
                                      FB.out.directory.list,
                                      batches){


  if(missing(GEX.out.directory.list)){GEX.out.directory.list <- "none"}
  if(missing(VDJ.out.directory.list)){VDJ.out.directory.list <- "none"}
  if(missing(FB.out.directory.list)){FB.out.directory.list <- "none"}

  if(GEX.out.directory.list[[1]] != "none" & VDJ.out.directory.list[[1]] != "none"){
    if(length(VDJ.out.directory.list) != length(GEX.out.directory.list)){stop("Different number of paths supplied for VDJ and GEX")}}

  if(GEX.out.directory.list[[1]] == "none" & missing(batches)){batches <- rep(1, length(VDJ.out.directory.list))}
  if(VDJ.out.directory.list[[1]] == "none" & missing(batches)){batches <- rep(1, length(GEX.out.directory.list))}
  if(GEX.out.directory.list[[1]] != "none" & VDJ.out.directory.list[[1]] != "none" & missing(batches)){
    batches <- rep(1, length(GEX.out.directory.list))
  }

  if(length(batches) != length(VDJ.out.directory.list) & length(batches) != length(GEX.out.directory.list)){
    stop("Batch input number does not match directory list input")
  }

  batches <- paste0("b", batches)

  #Remove possible backslash at the end of the input path
  for(k in 1:length(GEX.out.directory.list)){
    GEX.out.directory.list[[k]]<-  gsub("/$", "", GEX.out.directory.list[[k]])
  }
  for(k in 1:length(VDJ.out.directory.list)){
    VDJ.out.directory.list[[k]]<-  gsub("/$", "", VDJ.out.directory.list[[k]])
  }
  for(k in 1:length(FB.out.directory.list)){
    FB.out.directory.list[[k]]<-  gsub("/$", "", FB.out.directory.list[[k]])
  }

  params <- c(paste0(do.call("c",as.list(VDJ.out.directory.list)), collapse = " / "),
              paste0(do.call("c",as.list(GEX.out.directory.list)), collapse = " / "),
              paste0(do.call("c",as.list(FB.out.directory.list)), collapse = " / "))


  names(params) <- c("VDJ.out.directory.list",
                     "GEX.out.directory.list",
                     "FB.out.directory.list")

  gex.list <- list()
  out.list <- list() #open list containing matrices

  #### Load in VDJ ####
  vdj.loaded <- F
  if(VDJ.out.directory.list[[1]] != "none"){
    vdj_load_error <- tryCatch({

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
        warning("! metrics_summary.csv file not available in at least one of the VDJ input directories. Loading will be skipped \n")
        metrics.table <- as.list(rep("none",length(VDJ.out.directory.list)))
      }

      #NEW in cellranger 6.1 => We load the contigs independently from the annotations file. This allows us to return the full contig sequence (important for MIXCR and other tools later on) despite not trimming and aligning or potentially not having access to the all_contig_annotations.json file (see below)
      #We then add the raw contigs to the contig table loaded from filtered_contig_annotation.csv above
      VDJ.out.directory_raw_contigs <- paste(VDJ.out.directory.list,"/filtered_contig.fasta",sep="")
      raw.contig.table <- lapply(VDJ.out.directory_raw_contigs, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))

      for(ikj in 1:length(raw.contig.table)){
        #coercing this to a dataframe
        raw.contig.table[[ikj]] <- data.frame("contig_id" = names(raw.contig.table[[ikj]]), "raw_contig" = as.character(raw.contig.table[[ikj]]))
        #merge directly with contig dataframe to have annotations and sequence in one place for later
        contig.table[[ikj]] <- merge(contig.table[[ikj]], raw.contig.table[[ikj]], by = "contig_id", all.x = T, all.y = F) #making sure to only merge in raw contig sequences for contigs which are present in the table containing annotations
        if(sum(is.na(contig.table[[ikj]]$raw_contig)) > 0.5*nrow(contig.table[[ikj]])){
          warning("! Merging of raw contigs and filtered_contig_annotations showed unsusually low overlap \n")
        }
      }

      #New annotations read in
      if(all(file.exists(paste(VDJ.out.directory.list,"/all_contig_annotations.json",sep="")))){
        VDJ.out.directory_annotations <- paste(VDJ.out.directory.list,"/all_contig_annotations.json",sep="")
        annotations.list <- lapply(VDJ.out.directory_annotations, function(x) jsonlite::read_json(x))

        # pulls out the three important features: featureRegions, and of featureRegions. Used for trimming where the V region starts and where the C region ends.
        annotations.table <- list()
        for(i in 1:length(annotations.list)){
          #get annotation table to make VDJ_barcode_iteration later on faster
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
      } else { #in at least one directory the all_contig_annotations.json was not found
        warning("Warning: At least one VDJ input directory is missing the file all_contig_annotations.json. Without this, accurate trimming and aligning of sequences is not possible. Setting trim.and.align to FALSE and proceeding. For an alternate mode of aligment please refer to the function VDJ_call_MIXCR \n")
        trim.and.align <- F
        annotations.table <- as.list(rep("none", length(contig.table)))
      }



      #FOR AIRR COMPATIBILITY
      VDJ.out.directory_AIRR <- paste(VDJ.out.directory.list,"/airr_rearrangement.tsv",sep="")
      airr.table <- list()
      for(j in 1:length(VDJ.out.directory_AIRR)){
        if(file.exists(VDJ.out.directory_AIRR[[j]])){
          airr.table[[j]] <- utils::read.delim(VDJ.out.directory_AIRR[[j]], header = T)
        } else {
          warning(paste0("\n airr_rearrangements.tsv not found for sample ", j, ". AIRR compatibility for this sample will not be available"))
          airr.table[[j]] <- "none"
        }
      }

      vdj.loaded <- T
    }, error = function(e){
      message(paste0("Loading VDJ failed \n ",e))
})
  } else{
    clonotype.list <- as.list(rep("none",length(GEX.out.directory.list)))
    reference.list <- as.list(rep("none",length(GEX.out.directory.list)))
    annotations.table <- as.list(rep("none",length(GEX.out.directory.list)))
    contigs.table <- as.list(rep("none",length(GEX.out.directory.list)))
    VDJ.out.directory.list <- as.list(rep("none",length(GEX.out.directory.list)))
    metrics.table <- as.list(rep("none",length(GEX.out.directory.list)))
    airr.table <- as.list(rep("none",length(GEX.out.directory.list)))
  }

  #### Load in GEX ####
  gex.loaded <- F
  FB.loaded <- F
  if(GEX.out.directory.list[[1]] != "none"){
    gex_load_error <- tryCatch({
      #add the directory identifier
      if(stringr::str_detect(GEX.out.directory.list[[1]], "filtered_feature_bc_matrix") | stringr::str_detect(GEX.out.directory.list[[1]], "raw_feature_bc_matrix") | stringr::str_detect(GEX.out.directory.list[[1]], "sample_feature_bc_matrix")){
        #Nothing to append
        GEX.out.directory.list.p <- GEX.out.directory.list
      } else{

        if(dir.exists(paste(GEX.out.directory.list[[1]],"/filtered_feature_bc_matrix",sep=""))){ #checking only for first path assuming that all sample inputs are processed with the same cellranger function
          GEX.out.directory.list.p <- paste(GEX.out.directory.list,"/filtered_feature_bc_matrix",sep="")
          message("Setting GEX directory to provided path/filtered_feature_bc_matrix \n")
        } else if (dir.exists(paste(GEX.out.directory.list[[1]],"/sample_feature_bc_matrix",sep=""))){
          GEX.out.directory.list.p <- paste(GEX.out.directory.list,"/sample_feature_bc_matrix",sep="")
          message("Setting GEX directory to provided path/sample_feature_bc_matrix \n")
        } else {
          stop("The GEX directory filtered_feature_bc_matrix or sample_feature_bc_matrix was not found at the given path. Please revise GEX input paths")
        }
      }

      GEX.list <- lapply(GEX.out.directory.list.p, function(x) Seurat::Read10X(data.dir=x))

      #NEW CELLRANGER 6.1. dealing with the possibility of Feature barcode information being included in the GEX directory.
      #=> Procedure: check if there is feature barcode info in the GEX that was just loaded => if not, proceed as normal
      #=> if yes, isolate the matrices for each input sample into two flat lists: GEX.list and FB.list
      FB.list <- list()
      FB.loaded <- F
      #dealing with possible mixed GEX FB inputs or multiple FB input matrices from the same directory
      for(i in 1:length(GEX.list)){ #iterating over main list
        if(class(GEX.list[[i]]) == "list"){ #this returns true only if the GEX directory import contained more than one marices => i.e. there is a GEX and a FB matrix
          GEX_ind <- c() #open indices for GEX and FB list elements
          FB_ind <- c()
          for(j in 1:length(GEX.list[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(GEX.list[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              message(paste0("GEX input ", i, " element ", j, " contains > 100 features and will be loaded as GEX \n"))
              GEX_ind <- c(GEX_ind, j)
            } else if(nrow(GEX.list[[i]][[j]]) < 100){
              message(paste0("GEX input ", i, " element ", j, " contains < 100 features and will be loaded as FB \n"))
              FB_ind <- c(FB_ind, j)
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }
          if(length(FB_ind) > 0){
            FB.list <- c(FB.list, GEX.list[[i]][FB_ind]) #add this matrix to the FB list
            GEX.list[[i]] <- GEX.list[[i]][-FB_ind] #delete the matrix from the original list, so that only the GEX matrix remains
            #Moreover we want to prevent two user inputs with FB. Here we set FB.loaded to TRUE, so that the FB loading from disk module further down will be skipped.
            FB.loaded <- T
            #at this point we have two lists: the FB.list and the directory_read10x which only contains GEX info, but is still nested. So we have to flatten it
            GEX.list <- do.call(list, unlist(GEX.list, recursive=FALSE))
          }
        } else{
          FB.list[[i]] <- "none"
        }
      }

      #load the metrics file conditionally
      GEX.out.directory.list.metrics <- paste(GEX.out.directory.list,"/metrics_summary.csv",sep="")
      if(all(file.exists(GEX.out.directory.list.metrics))){
        GEX.metrics <- lapply(GEX.out.directory.list.metrics, function(x) utils::read.csv(x,sep=",",header=T, ))
      } else {
        GEX.metrics <- as.list(rep("none",length(GEX.out.directory.list)))
      }

      gex.loaded <- T
    }, error = function(e){
      message(paste0("Loading GEX failed \n", e))})
  } else{
    GEX.list <- as.list(rep("none",length(VDJ.out.directory.list)))
    GEX.metrics <- as.list(rep("none",length(VDJ.out.directory.list)))
    GEX.out.directory.list <- as.list(rep("none",length(VDJ.out.directory.list)))
    FB.list <- as.list(rep("none",length(VDJ.out.directory.list)))
  }

  #### Load in FB ####

  #! Not checking if FB have already been loaded as part of GEX. If FB directory is supplied, any FB data loaded by GEX is therefore overwritten. This should allow for more flexibility of the user without having to realign the whole GEX data
  if(FB.out.directory.list[[1]] != "none"){

    gex_load_error <- tryCatch({suppressWarnings({

      #add the directory identifier
      if(stringr::str_detect(FB.out.directory.list[[1]], "filtered_feature_bc_matrix")){
        FB.out.directory.list.p <- FB.out.directory.list
        warning("! Feature barcode path was specified explicitely to filtered_feature_bc_matrix. For better rates we recommend using the raw_feature_bc_matrix folder content ! \n ")
      } else if (stringr::str_detect(FB.out.directory.list[[1]], "raw_feature_bc_matrix")){
        #Nothing to append
        FB.out.directory.list.p <- FB.out.directory.list
        message("Loading feature barcodes from raw_feature_bc_matrix folder \n")
      } else{
        FB.out.directory.list.p <- paste(FB.out.directory.list,"/raw_feature_bc_matrix",sep="")
        message("Loading feature barcodes from raw_feature_bc_matrix folder \n")
      }
      #Actually loading the data. Critical: if Cellranger 6.1.0 count was run with --libraries input containing both GEX and Feature Barcodes, reading it will result in a list of matrices instead of a single matrix. => the next section deals with this
      n_not_loaded <- 0
      FB.list <- lapply(FB.out.directory.list.p, function(x) tryCatch({Seurat::Read10X(data.dir=x)},
                                                                                error = function(e){
                                                                                  n_not_loaded <- n_not_loaded + 1
                                                                                  return(NULL)})) #Catching an error which occurs when trying to load the directory "PLACEHOLDER/filtered_feature_bc_matrix" if a placeholder was provided in input directories.

      #next we check if all elements of the new loaded list are null => if so we stop and report to the user that apparently the paths he provided were not correct
      if(n_not_loaded == length(FB.list)){
        stop("FB data loading failed from provided paths. Please revise input paths")
      }

      #Replacing null values with "none" which will be picked up by the VGM function
      for(i in 1:length(FB.list)){
        if(is.null(FB.list[[i]])){
          FB.list[[i]] <- "none"
        }
      }



      #dealing with possible mixed GEX FB inputs or multiple FB input matrices from the same directory
      for(i in 1:length(FB.list)){ #iterating over main list
        if(class(FB.list[[i]]) == "list"){ #this returns true only if the FB directory import contained more than one marices
          to_del <- c()
          for(j in 1:length(FB.list[[i]])){ #Now iterating over the elements of this particular FB directory input
            if(nrow(FB.list[[i]][[j]]) > 100){ #Checking whether this matrix may contain cite seq or feature barcodes. If the number of features is over 100, the matrix almost certainly contains GEX information. We will discard this matrix
              message(paste0("Feature barcode input ", i, " element ", j, " contains > 100 features and likely corresponds to GEX data. This matrix will be removed from further FB processing \n"))
              to_del <- c(to_del, j) #Will be deleted later to not mess up the loop
            }
          }
          if(j > 2){ #If there are more two matrices for one sample. At the moment (18.8.21) we expect max 2 matrices per sample (1 GEX 1 FB). Potentially, this will change with future versions of cellrangers and the possibility to deal with CITE-seq data as well. For now we unfortunately stop the function if such input is provided
            stop(paste0("GEX loading error: for GEX directory input ", i, " the Read10x function returned more than 2 matrices. This is a likely a result of running Cellranger count or Cellranger multi with > 1 directory input for GEX or feature barcodes per sample or of having processed additional feature barcode data such as from Cite-seq. Currently this function is only capable of processing 2 output matrices from each GEX directory. Further compatibility may be added in the future."))
          }
          if(length(to_del) > 0){
            FB.list[[i]] <- FB.list[[i]][-j] #deleting
          }
        }
      }
      #now flatten the remaining list
      FB.list <- do.call(list, unlist(FB.list, recursive=FALSE))
      #Done => result should be a non-nested list of matrices only containing FB information. With this we can move forward

      FB.loaded <- T
    })}, error = function(e){
      message(paste0("Loading FB failed \n", e))})

  } else {
    if(FB.loaded == F){
    FB.loaded <- F
    if(gex.loaded){
      FB.list <- as.list(rep("none",length(GEX.out.directory.list)))
    } else {
      FB.list <- as.list(rep("none",length(VDJ.out.directory.list)))
    }
    } else {
      #Nothing => FB has been loaded from GEX, and no FB directory was supplied to overwrite that
    }
  }


  #combine into an output nested list on a by sample level
  out.list <- list()

      for(i in 1:length(VDJ.out.directory.list)){
    out.list[[i]] <- list(list(clonotype.list[[i]], reference.list[[i]], annotations.table[[i]], contig.table[[i]], metrics.table[[i]], airr.table[[i]]),list(GEX.list[[i]], GEX.metrics[[i]]),list(FB.list[[i]], FB.out.directory.list), VDJ.out.directory.list[[i]], GEX.out.directory.list[[i]], batches[i])


    names(out.list[[i]]) <- c("VDJ", "GEX","FB", "VDJ.path", "GEX.path", "Batch")
    names(out.list[[i]][[1]]) <- c("clonotypes.csv", "concat_ref.fasta", "all_contig_annotations.json", "filtered_contig_annotations.csv", "metrics_summary.csv", "airr_rearrangement.tsv")
    names(out.list[[i]][[2]]) <- c("GEX", "metrics_summary.csv")
    names(out.list[[i]][[3]]) <- c("FB", "FB.path") #Reason for adding this as a nested list: this gives some more room if we ever want to implement CITE seq support, so that the rest of the data structure can remain stable. Cite seq matrices could be saved as out.list[[i]][[3]][[2]] and the FB.path and CITE.path at out.list[[i]][[3]][[3]] and out.list[[i]][[3]][[4]], respectively
      }
  names(out.list) <- c(paste0("Sample", c(1:length(VDJ.out.directory.list))))

  return(out.list)
}




