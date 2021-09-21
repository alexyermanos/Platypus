#' Gives stats on number and quality of reads.
#' @param VDJ.out.directory List of paths with each element containing the path to the output of cellranger VDJ runs. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @param GEX.out.directory OPTIONAL list of paths with each element containing the path to the output of cellranger GEX runs. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @param sample.names OPTIONAL: an array of the same length as the input VDJ.out.directory list with custom names for each sample. If not provided samples will be numbered by processing order
#' @param metrics10x Whether to append metrics_summary.csv information provided by Cellranger for both VDJ and GEX. Defaults to T
#' @param save.csv Boolean. Defaults to TRUE. Whether to directly save the resuts as a comma delimited .csv file in the current working directory.
#' @param filename Character ending in .csv. Filename to save .csv as.
#' @return returns a single matrix where the rows are individual cells and the columns are repertoire features.
#' @export
#' @examples
#' \dontrun{
#' stats <- VDJ_GEX_stats(VDJ.out.directory = VDJ.out.directory.list
#' ,GEX.out.directory = GEX.out.directory.list,sample.names = c(1:4)
#' ,metrics10x = T,save.csv = T ,filename = "stats.csv")
#' }

VDJ_GEX_stats <- function(VDJ.out.directory,
                          GEX.out.directory,
                          sample.names,
                          metrics10x,
                          save.csv,
                          filename){


  if(missing(save.csv)) save.csv <- T
  if(missing(filename)) filename <- "VDJ_stats.csv"
  if(missing(metrics10x)) metrics10x <- F

  vdj.loaded <- F
  if(missing(VDJ.out.directory)==F){
    print("Reading in input files")
    VDJ.out.directory_clonotypes <- paste(VDJ.out.directory,"/clonotypes.csv",sep="")
    VDJ.out.directory_reference <- paste(VDJ.out.directory,"/concat_ref.fasta",sep="")
    VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/all_contig_annotations.csv",sep="")
    VDJ.out.directory_annotations <- paste(VDJ.out.directory,"/all_contig_annotations.json",sep="")
    VDJ.out.directory_metrics <- paste(VDJ.out.directory,"/metrics_summary.csv",sep="")

    #needed for VDJ_analyze module
    clonotype.list <- lapply(VDJ.out.directory_clonotypes, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))

    #needed for VDJ_per_cell_matrix module
    reference.list <- lapply(VDJ.out.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    annotations.list <- lapply(VDJ.out.directory_annotations, function(x) jsonlite::read_json(x))
    contig.list <- lapply(VDJ.out.directory_contigs, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
    VDJ.metrics.list <- lapply(VDJ.out.directory_metrics, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
    vdj.loaded <- T
  } else {
    stop("No VDJ.out.directory supplied")
  }

  gex.loaded <- F
  if(missing(GEX.out.directory) == F){
    if(length(GEX.out.directory) == length(VDJ.out.directory)){
      GEX.out.directory_metrics <- paste(GEX.out.directory,"/metrics_summary.csv",sep="")
      GEX.metrics.list <- lapply(GEX.out.directory_metrics, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
      gex.loaded <- T
    }else {
      print("Length of GEX out directory not the same as VDJ out directory. GEX will not be loaded")
    }
  } else {
    print("No GEX.out.directory supplied")
  }

  if(missing(sample.names)){sample.names <- c(1:length(VDJ.out.directory))}

  VDJ.stats.list <- list()
  for(k in 1:length(clonotype.list)){

    print(paste0("Starting with ", k, " of ", length(clonotype.list)))
    VDJ.stats <- c()

    #gsub to be able to process TCRs as well
    contig.list[[k]]$chain <- gsub(pattern = "TRB", replacement = "IGH", contig.list[[k]]$chain)
    contig.list[[k]]$chain <- gsub(pattern = "TRA", replacement = "IGL", contig.list[[k]]$chain)

    clonotype.list[[k]]$cdr3s_aa <- gsub(pattern = "TRB:", replacement = "IGH:", clonotype.list[[k]]$cdr3s_aa)
    clonotype.list[[k]]$cdr3s_aa <- gsub(pattern = "TRA:", replacement = "IGL:", clonotype.list[[k]]$cdr3s_aa)

    #info on sample
    VDJ.stats[length(VDJ.stats)+1] <- VDJ.out.directory[[k]]
    names(VDJ.stats)[length(VDJ.stats)] <- "Repertoir path"

    VDJ.stats[length(VDJ.stats)+1] <- sample.names[k]
    names(VDJ.stats)[length(VDJ.stats)] <- "Sample name"


    #Get number of unique barcodes
    VDJ.stats[length(VDJ.stats)+1] <- length(unique(contig.list[[k]]$barcode))
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr unique barcodes"

    #generate lookup table with HC and LC counts and stats per barcode
    print("Getting lookup tables")
    holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "%",
                                         width = 100, style = 3, file = "")
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
      utils::setTxtProgressBar(value = nr_bar/(length(unique(contig.list[[k]]$barcode)) + nrow(clonotype.list[[k]])),pb = holding_bar)
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
      utils::setTxtProgressBar(value = (nr_bar+l)/(length(unique(contig.list[[k]]$barcode)) + nrow(clonotype.list[[k]])),pb = holding_bar)
      clonotype_ids <- append(clonotype_ids, clonotype.list[[k]]$clonotype_id[l])

      nr_HC <- append(nr_HC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IGH:"))
      nr_LC <- append(nr_LC,stringr::str_count(clonotype.list[[k]]$cdr3s_aa[l], "IG(K|L):"))
    }
    lookup_stats_clono <- data.frame(clonotype_ids,nr_HC,nr_LC)
    names(lookup_stats_clono) <- c("clonotype_ids","nr_HC","nr_LC")

    close(holding_bar)

    #number of barcodes with
    #is cell == true
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr barcodes is_cell"

    #number of is.cell with 1 HC and 1 LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1HC 1LC"

    #number of cells with 1 HC and 0 LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 0 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1HC 0LC"

    #number of cells with 0 HC and 1 LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 0 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 0HC 1LC"

    #number of cells with 2 or more HC and 1 LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC == 1 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more HC 1LC"

    #number of cells with 1 HC and 2 or more LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC > 1 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1HC 2 or more LC"

    #number of cells with 2 or more HC and 2 or more LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC > 1 & lookup_stats$is_cell == "true",])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more HC 2 or more LC"

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
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells all true and 1HC 1LC"

    #number of clonotypes
    VDJ.stats[length(VDJ.stats)+1] <- nrow(clonotype.list[[k]])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes"

    #number of clonotypes with exactly 1HC 1LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC == 1 & lookup_stats_clono$nr_LC == 1,])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes 1HC 1LC"

    #number of clonotypes with  < 1HC 1LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC + lookup_stats_clono$nr_LC < 2,])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes < 1HC 1LC"

    #number of clonotypes with  > 1HC 1LC
    VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats_clono[lookup_stats_clono$nr_HC + lookup_stats_clono$nr_LC > 2,])
    names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes > 1HC 1LC"

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

  print("Compiling table")
  if(metrics10x == T){
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
      if(gex.loaded == T){

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
      print(paste0("Adding 10x metrix failed: ", e))})

    VDJ.stats.all <- cbind(VDJ.stats.all, VDJ.metrics.all)
  } else {}
  if(save.csv){utils::write.csv(VDJ.stats.all, file = filename)}

  return(VDJ.stats.all)
}
