#'AIRR to Platypus V3 VGM compatibility function
#'
#'@description Loads in and converts input AIRR-compatible tsv file(s) into the Platypus VGM object format.All compulsory AIRR data columns are needed. Additionally, the following columns are required: v_call, cell_id, clone_id. If trim.and.align is set to TRUE additionally the following columns are needed: v_sequence_start, j_sequence_end.
#' Note on TRUST4 input: TRUST4 (https://doi.org/10.1038/s41592-021-01142-n2) is a newly alignment tool for VDJ data by the Shirley lab. It is able to also extract VDJ sequences from 10x GEX data. We are actively testing TRUST4 as an alternative to Cellranger and can not give recommendations as of now. This function does support the conversion of TRUST4 airr output data into the Platypus VGM format. In that case, an extra column will be added describing whether the full length VDJ sequence was extracted for any given cell and chain.
#'@param AIRR.input Source of the AIRR table(s) as a list. There are 2 available input options: 1. 1. List with local paths to .tsv files / 3. List of AIRR tables loaded in as R objects within the current R environment.
#'@param get.VDJ.stats Boolean. Defaults to TRUE. Whether to generate summary statistics on repertoires and output those as output_VGM[[3]]
#'@param VDJ.combine Boolean. Defaults to TRUE. Whether to integrate repertoires. A sample identifier will be appended to each barcode both. Highy recommended for all later functions
#'@param trim.and.align Boolean. defaults to FALSE. Whether to trim VJ/VDJ seqs and add information from alignment in AIRR dataframe columns. ! No alignment is done here, instead, columns containing alignment information in the AIRR dataframes are reformatted.
#'@param filter.overlapping.barcodes.VDJ Boolean. defaults to TRUE. Whether to remove barcodes which are shared among samples in the GEX analysis. Shared barcodes normally appear at a very low rate.
#'@param group.id vector with integers specifying the group membership. c(1,1,2,2) would specify the first two elements of the input AIRR list are in group 1 and the third/fourth input elements will be in group 2.
#' @param verbose Writes runtime status to console. Defaults to FALSE
#'@return A VDJ_GEX_Matrix object used in Platypus V3 as an input to most analysis and plotting functions
#' @export
#' @examples
#' \donttest{
#' try({
#' VGM <- PlatypusDB_AIRR_to_VGM(AIRR.input =
#' list("~/pathto/s1/airr_rearrangement.tsv", "~pathto/s2/airr_rearrangement.tsv"),
#' VDJ.combine = TRUE, group.id = c(1,2), filter.overlapping.barcodes.VDJ = TRUE)
#' })
#' }

PlatypusDB_AIRR_to_VGM <- function(AIRR.input,
                                   get.VDJ.stats,
                                   VDJ.combine,
                                   trim.and.align,
                                   filter.overlapping.barcodes.VDJ,
                                   group.id,
                                   verbose){

  if(missing(verbose)) verbose <- FALSE

  clonotype_id_10x <- NULL
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL

  ############## Def of helper function
  #Helper function called in VDJ_GEX_matrix. Do not run as standalone!
  #FUN to call in parlapply mclapply or lapply
  AIRR_barcode_VDJ_iteration <- function(barcodes,
                                    airrs){


    curr.contigs <- airrs[airrs$cell_id == barcodes & tolower(as.character(airrs$productive)) == "true",]

    #set up data structure
    cols <- c("barcode","sample_id","group_id","clonotype_id","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa", "VJ_cdr3s_aa","VDJ_cdr3s_nt", "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain","VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw","VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_trimmed_ref","VJ_trimmed_ref")

    #trust 4 full lenght compatibility extra column
    if(all(c("consensus_count", "complete_vdj") %in% names(curr.contigs))){
      cols <- c(cols,"VDJ_complete_vdj","VJ_complete_vdj", "VDJ_consensus_count", "VJ_consensus_count")
      trust4 <- TRUE
    } else{
      trust4 <- FALSE
    }

    curr.barcode <- stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 1)), cols)

    #fill in information that do not need processing
    #Contig info on light/a and heavy/b chains is put into different columns (see cols)
    #If only one contig is available, the fields of the other are left blank
    #If more than two contigs of one chain (e.g. 2 TRB) are present, the elements will be pasted separated by a ";" into the relevant fields (in the case of TRB, into the Hb columns)

    #Get number of chains
    HC_count <- sum(stringr::str_count(curr.contigs$v_call, pattern = "(TRG|TRB|IGH)"))
    LC_count <- sum(stringr::str_count(curr.contigs$v_call, pattern = "(TRD|TRA|IGK|IGL)"))

    #In this case we need to make much less effort with pasting together, so we can save time
    if(HC_count == 1 & LC_count == 1){

      if(which(stringr::str_detect(curr.contigs$v_call, "(TRD|TRA|IGK|IGL)")) == 1){ #make row 1 the heavy chain in case it is not already
        curr.contigs <- curr.contigs[c(2,1),]}

      #fill in the pasted info to curr.barcode directly
      curr.barcode$barcode <- curr.contigs$cell_id[1]
      curr.barcode$clonotype_id_10x <- curr.contigs$clone_id[1]
      curr.barcode$sample_id <- ""
      curr.barcode$group_id <- ""
      if(stringr::str_detect(curr.contigs$v_call[1], "TR") | stringr::str_detect(curr.contigs$v_call[2], "TR")){curr.barcode$celltype <- "T cell"
      } else if(stringr::str_detect(curr.contigs$v_call[1], "IG") | stringr::str_detect(curr.contigs$v_call[2], "IG")) {curr.barcode$celltype <- "B cell"
      } else {curr.barcode$celltype <- "Unkown"}

      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count

      curr.barcode$VDJ_cdr3s_aa <- curr.contigs$junction_aa[1]
      curr.barcode$VJ_cdr3s_aa <- curr.contigs$junction_aa[2]
      curr.barcode$VDJ_cdr3s_nt <- curr.contigs$junction[1]
      curr.barcode$VJ_cdr3s_nt <- curr.contigs$junction[2]
      curr.barcode$VDJ_chain_contig <- curr.contigs$sequence_id[1]
      curr.barcode$VJ_chain_contig <- curr.contigs$sequence_id[2]
      curr.barcode$VDJ_chain <- substr(curr.contigs$c_call[1], 1,3)
      curr.barcode$VJ_chain <- substr(curr.contigs$c_call[2], 1,3)
      curr.barcode$VDJ_vgene <- curr.contigs$v_call[1]
      curr.barcode$VJ_vgene <- curr.contigs$v_call[2]
      curr.barcode$VDJ_dgene <- curr.contigs$d_call[1]
      curr.barcode$VDJ_jgene <- curr.contigs$j_call[1]
      curr.barcode$VJ_jgene <- curr.contigs$j_call[2]
      curr.barcode$VDJ_cgene <- curr.contigs$c_call[1]
      curr.barcode$VJ_cgene <- curr.contigs$c_call[2]
      curr.barcode$VDJ_raw_consensus_id <- ""
      curr.barcode$VJ_raw_consensus_id <- ""
      if(trust4){
        curr.barcode$VDJ_complete_vdj <- curr.contigs$complete_vdj[1]
        curr.barcode$VJ_complete_vdj <- curr.contigs$complete_vdj[2]
        curr.barcode$VDJ_consensus_count <- curr.contigs$consensus_count[1]
        curr.barcode$VJ_consensus_count <- curr.contigs$consensus_count[2]
      }

    } else { # this for cells with aberrant chain numbers

      contigs_pasted <- stats::setNames(data.frame(matrix(ncol = ncol(curr.contigs), nrow = length(unique(curr.contigs$v_call)))), names(curr.contigs)) #the dataframe may be one or two rows too long, this will not matter / ROW 1 = Heavy chain information / ROW 2 = Light chain information. This order is maintained even if one of the two chains is not present!
      if(nrow(contigs_pasted) == 3) contigs_pasted <- contigs_pasted[-3,]

      #Heavy/b chain count
      if(HC_count == 1){
        contigs_pasted[1,] <- curr.contigs[stringr::str_detect(curr.contigs$v_call, pattern = "(TRG|TRB|IGH)"),]
      } else if(HC_count == 0){
        contigs_pasted[1,] <- ""
      } else if(HC_count > 1){
        for(k in 1:ncol(curr.contigs)){
          contigs_pasted[1,k] <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$v_call, pattern = "(TRG|TRB|IGH)")), k], collapse = ";")
        }
      }
      ### Order of CDRs with multiple chains is determined here

      #Light/a chain count
      if(LC_count == 1){
        contigs_pasted[2,] <- curr.contigs[stringr::str_detect(curr.contigs$v_call, pattern = "(TRD|TRA|IGK|IGL)"),]
      } else if(LC_count == 0){
        contigs_pasted[2,] <- ""
      } else if(LC_count > 1){
        for(k in 1:ncol(curr.contigs)){
          contigs_pasted[2,k]  <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$v_call, pattern = "(TRD|TRA|IGK|IGL)")),k],collapse = ";")
        }
      }

      #fill in the pasted info to curr.barcode
      curr.barcode$barcode <- curr.contigs$cell_id[1]
      curr.barcode$clonotype_id_10x <- curr.contigs$clone_id[1]
      curr.barcode$sample_id <- ""
      curr.barcode$group_id <- ""
      if(stringr::str_detect(contigs_pasted$v_call[1], "TR") | stringr::str_detect(contigs_pasted$v_call[2], "TR")){curr.barcode$celltype <- "T cell"
      } else if(stringr::str_detect(contigs_pasted$v_call[1], "IG") | stringr::str_detect(contigs_pasted$v_call[2], "IG")) {curr.barcode$celltype <- "B cell"
      } else {curr.barcode$celltype <- "Unkown"}

      curr.barcode$Nr_of_VDJ_chains <- HC_count
      curr.barcode$Nr_of_VJ_chains <- LC_count

      curr.barcode$VDJ_cdr3s_aa <- contigs_pasted$junction_aa[1]
      curr.barcode$VJ_cdr3s_aa <- contigs_pasted$junction_aa[2]
      curr.barcode$VDJ_cdr3s_nt <- contigs_pasted$junction[1]
      curr.barcode$VJ_cdr3s_nt <- contigs_pasted$junction[2]
      curr.barcode$VDJ_chain_contig <- contigs_pasted$sequence_id[1]
      curr.barcode$VJ_chain_contig <- contigs_pasted$sequence_id[2]
      curr.barcode$VDJ_chain <- substr(contigs_pasted$v_call[1], 1,3)
      curr.barcode$VJ_chain <- substr(contigs_pasted$v_call[2], 1,3)
      curr.barcode$VDJ_vgene <- contigs_pasted$v_call[1]
      curr.barcode$VJ_vgene <- contigs_pasted$v_call[2]
      curr.barcode$VDJ_dgene <- contigs_pasted$d_call[1]
      curr.barcode$VDJ_jgene <- contigs_pasted$j_call[1]
      curr.barcode$VJ_jgene <- contigs_pasted$j_call[2]
      curr.barcode$VDJ_cgene <- contigs_pasted$v_call[1]
      curr.barcode$VJ_cgene <- contigs_pasted$v_call[2]
      curr.barcode$VDJ_raw_consensus_id <- ""
      curr.barcode$VJ_raw_consensus_id <- ""
      if(trust4){
        curr.barcode$VDJ_complete_vdj <- curr.contigs$complete_vdj[1]
        curr.barcode$VJ_complete_vdj <- curr.contigs$complete_vdj[2]
        curr.barcode$VDJ_consensus_count <- curr.contigs$consensus_count[1]
        curr.barcode$VJ_consensus_count <- curr.contigs$consensus_count[2]
      }

    } #end if HC | LC count > 1

    #HEAVY CHAIN / TRB
    #CHECK IF THERE IS 1 2 or 0 chains to process
    if(HC_count == 1){

      HC_index <- which(stringr::str_detect(curr.contigs$v_call, "(TRG|TRB|IGH)"))

      #extract match
      curr.barcode$VDJ_sequence_nt_raw <- curr.contigs$sequence[HC_index]
      if(trim.and.align){
      curr.barcode$VDJ_sequence_nt_trimmed <- substr(curr.barcode$VDJ_sequence_nt_raw, as.numeric(curr.contigs$v_sequence_start[HC_index]), as.numeric(curr.contigs$j_sequence_end[HC_index]))
      curr.barcode$VDJ_sequence_aa <- curr.contigs$sequence_aa[HC_index]
      curr.barcode$VDJ_trimmed_ref <- substr(curr.contigs$germline_alignment[HC_index], as.numeric(curr.contigs$v_sequence_start[HC_index]), as.numeric(curr.contigs$j_sequence_end[HC_index]))
      } else {
        curr.barcode$VDJ_sequence_nt_trimmed <- ""
        curr.barcode$VDJ_sequence_aa <- ""
        curr.barcode$VDJ_trimmed_ref <- ""
      }

    } else if(HC_count == 0){

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
      HC_rows <- curr.contigs[stringr::str_detect(curr.contigs$v_call, "(TRG|TRB|IGH)"),]
      #looping contigs in annotation
      for(l in 1:nrow(HC_rows)){
            #get sequence
            to_paste <- append(to_paste, HC_rows$sequence[l])
            if(trim.and.align == TRUE){

              to_paste_trimmed <- c(to_paste_trimmed, substr(HC_rows$sequence[l], as.numeric(HC_rows$v_sequence_start[l]), as.numeric(HC_rows$j_sequence_end[l])))

              to_paste_aa <- c(to_paste_aa, HC_rows$sequence_aa[l])

              to_paste_ref_trimmed <- c(to_paste_ref_trimmed, substr(HC_rows$germline_alignment[l], as.numeric(HC_rows$v_sequence_start[l]), as.numeric(HC_rows$j_sequence_end[l])))

            } else {
              to_paste_trimmed <- ""
              to_paste_aa <- ""
              to_paste_ref_trimmed <- ""
            }
          }
      curr.barcode$VDJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
      curr.barcode$VDJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VDJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VDJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }
    HC_rows <- NULL

    #Light/a
    if(LC_count == 1){

      LC_index <- which(stringr::str_detect(curr.contigs$v_call, "(TRD|TRA|IGK|IGL)"))

      #extract match
      curr.barcode$VJ_sequence_nt_raw <- curr.contigs$sequence[LC_index]
      if(trim.and.align){
        curr.barcode$VJ_sequence_nt_trimmed <- substr(curr.barcode$VJ_sequence_nt_raw, as.numeric(curr.contigs$v_sequence_start[LC_index]), as.numeric(curr.contigs$j_sequence_end[LC_index]))
        curr.barcode$VJ_sequence_aa <- curr.contigs$sequence_aa[LC_index]
        curr.barcode$VJ_trimmed_ref <- substr(curr.contigs$germline_alignment[LC_index], as.numeric(curr.contigs$v_sequence_start[2]), as.numeric(curr.contigs$j_sequence_end[LC_index]))
      } else {
        curr.barcode$VJ_sequence_nt_trimmed <- ""
        curr.barcode$VJ_sequence_aa <- ""
        curr.barcode$VJ_trimmed_ref <- ""
      }

    } else if(LC_count == 0){

      curr.barcode$VJ_sequence_nt_trimmed <- ""
      curr.barcode$VJ_sequence_aa <- ""
      curr.barcode$VJ_trimmed_ref <- ""

    } else if(LC_count > 1){ #MORE THAN ONE LC
      #from the annotations extract sequence and paste
      #Heavy/b
      to_paste <- c()
      to_paste_trimmed <- c()
      to_paste_aa <- c()
      to_paste_ref_trimmed <- c()
      LC_rows <- curr.contigs[stringr::str_detect(curr.contigs$v_call, "(TRD|TRA|IGK|IGL)"),]
      #looping contigs in annotation
      for(l in 1:nrow(LC_rows)){
        #get sequence
        to_paste <- append(to_paste, LC_rows$sequence[l])
        if(trim.and.align == TRUE){

          to_paste_trimmed <- c(to_paste_trimmed, substr(LC_rows$sequence[l], as.numeric(LC_rows$v_sequence_start[l]), as.numeric(LC_rows$j_sequence_end[l])))

          to_paste_aa <- c(to_paste_aa, LC_rows$sequence_aa[l])

          to_paste_ref_trimmed <- c(to_paste_ref_trimmed, substr(LC_rows$germline_alignment[l], as.numeric(LC_rows$v_sequence_start[l]), as.numeric(LC_rows$j_sequence_end[l])))

        } else {
          to_paste_trimmed <- ""
          to_paste_aa <- ""
          to_paste_ref_trimmed <- ""
        }
      }
      curr.barcode$VJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
      curr.barcode$VJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
      curr.barcode$VJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
      curr.barcode$VJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    }
    LC_rows <- NULL

    return(curr.barcode)
  }

  ##################################################################################### STOP barcode_VDJ_iteration

  #Gets statistics on VDJ and GEX
  AIRR_VDJ_stats_int <- function(airr.list){####START VDJ_GEX_stats

    ### VDJ stats - mainly info coming from the annotated contigs csv
    VDJ.stats.list <- list()
    for(k in 1:length(airr.list)){

      VDJ.stats <- c()

      #gsub to be able to process TCRs as well
      airr.list[[k]]$chain <- substr(airr.list[[k]]$v_call, 1,3)

      airr.list[[k]]$chain <- gsub(pattern = "TRB|TRG", replacement = "IGH", airr.list[[k]]$chain)
      airr.list[[k]]$chain <- gsub(pattern = "TRA|TRD", replacement = "IGL", airr.list[[k]]$chain)

      #info on sample
      VDJ.stats[length(VDJ.stats)+1] <- paste0("s", k)
      names(VDJ.stats)[length(VDJ.stats)] <- "Sample name"

      #Get number of unique barcodes
      VDJ.stats[length(VDJ.stats)+1] <- length(unique(airr.list[[k]]$cell_id))
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr unique barcodes"

      barcodes <- c()
      nr_HC <- c()
      nr_LC <- c()
      is_cell <-c()
      high_confidence <- c()
      productive <- c()
      full_length <- c()
      nr_bar <- 0
      for(j in unique(airr.list[[k]]$cell_id)){
        nr_bar <- nr_bar + 1
        #utils::setTxtProgressBar(value = nr_bar/(length(unique(contig.list[[k]]$barcode)) + nrow(clonotype.list[[k]])),pb = holding_bar)
        barcodes <- append(barcodes, j)

        productive <- append(productive, min(airr.list[[k]]$productive[which(airr.list[[k]]$cell_id == j)]))
        nr_HC <- append(nr_HC,stringr::str_count(paste0(airr.list[[k]]$chain[which(airr.list[[k]]$cell_id== j)],collapse = ""), "IGH"))
        nr_LC <- append(nr_LC,stringr::str_count(paste0(airr.list[[k]]$chain[which(airr.list[[k]]$cell_id == j)],collapse = ""), "IG(K|L)"))
      }

      lookup_stats <- data.frame(barcodes,nr_HC,nr_LC,productive)
      names(lookup_stats) <- c("barcodes","nr_HC","nr_LC","productive")

      #number of barcodes with
      #is cell ==TRUE
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$productive == 1,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr barcodes productive"

      #number of is.cell with 1 HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 1,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 1VJ"

      #number of cells with 1 HC and 0 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC == 0 ,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 0VJ"

      #number of cells with 0 HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 0 & lookup_stats$nr_LC == 1 ,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 0VDJ 1VJ"

      #number of cells with 2 or more HC and 1 LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC == 1 ,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more VDJ 1VJ"

      #number of cells with 1 HC and 2 or more LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC == 1 & lookup_stats$nr_LC > 1 ,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 1VDJ 2 or more VJ"

      #number of cells with 2 or more HC and 2 or more LC
      VDJ.stats[length(VDJ.stats)+1] <- nrow(lookup_stats[lookup_stats$nr_HC > 1 & lookup_stats$nr_LC > 1 ,])
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr cells 2 or more VDJ 2 or more VJ"

      #number of clonotypes
      VDJ.stats[length(VDJ.stats)+1] <- length(unique(airr.list[[k]]$clone_id))
      names(VDJ.stats)[length(VDJ.stats)] <- "Nr clonotypes"

      #percentages
      VDJ.stats.perc <- rep(NA, length(VDJ.stats)-2)
      VDJ.stats.perc[1] <- round(as.numeric(VDJ.stats[3]) / as.numeric(VDJ.stats[3]) *100, digits = 1) #for barcode items
      VDJ.stats.perc[c(2:13)] <- round(as.numeric(VDJ.stats[c(4:9)]) / as.numeric(VDJ.stats[4]) *100, digits = 1) #for barcode and is_cell items

      names(VDJ.stats.perc) <- paste0("% ", names(VDJ.stats[3:length(VDJ.stats)])) #set names to VDJ.stats.perc
      VDJ.stats <- c(VDJ.stats, VDJ.stats.perc) #combine vectors

      VDJ.stats.df <- as.data.frame(t(data.frame(VDJ.stats)))
      names(VDJ.stats.df) <- names(VDJ.stats)
      VDJ.stats.list[[k]] <- VDJ.stats.df

    }
    VDJ.stats.all <- do.call(rbind, VDJ.stats.list)

    return(VDJ.stats.all)
  }

  ########################################################################################### STOP VDJ_GEX_stats



  #combine all samples into one table
  if(missing(get.VDJ.stats)) get.VDJ.stats <- TRUE
  if(missing(VDJ.combine)) VDJ.combine <- TRUE
  if(missing(filter.overlapping.barcodes.VDJ)) filter.overlapping.barcodes.VDJ <- TRUE
  if(missing(trim.and.align)) trim.and.align <- TRUE

  vdj_loaded <- FALSE
  #get input as list
  if(inherits(AIRR.input,"list")){
    if(inherits(AIRR.input[[1]],"data.frame")){ #case 1.
      if(verbose) message("\n Dataframe input detected")
      airr.list <- AIRR.input
      airr.names <- "Input from R enviroment"
      vdj_loaded <- TRUE

    } else if(inherits(AIRR.input[[1]],"character")){ #case 3
      if(verbose) message("\n Local paths input detected. Loading in tables")
      airr.names <- paste0(unlist(AIRR.input),collapse = ";")
      vdj_load_error <- tryCatch({
        airr.list <- list()
        for(j in 1:length(AIRR.input)){
          if(file.exists(AIRR.input[[j]])){
            airr.list[[j]] <- utils::read.delim(AIRR.input[[j]], header = T)

            #trust 4 full lenght compatibility extra column
            if(all(c("consensus_count", "complete_vdj") %in% names(airr.list[[j]]))){
              if(verbose) message(paste0("\n TRUST4 input format detected for sample ", j, ". Appending extra columns for complete_vdj and consensus_count to output VGM"))
            }

            vdj_loaded <- TRUE
          } else {
            if(verbose) warning(paste0(" File not found for sample ", j, ". Skipping this sample..."))
            airr.list[[j]] <- "none"
          }
        }

      }, error = function(e){e
        message(e)})
      if(inherits(vdj_load_error,"error")){
        warning("Loading airr_rearrangement from disk failed")}
    }
  } else {
    stop("\n Please provide AIRR inputs as a list. e.g. for paths list(~/data/project/airr_rearrangement.tsv)")
  }

  if(vdj_loaded == FALSE){
    stop("\n None of the specified file paths were found")
  }


  if(missing(group.id)) group.id <- c(1:length(airr.list))

  #save runtime parameters for later
  #params <- do.call("rbind", as.list(environment()))
  params <- c(airr.names,
              VDJ.combine,
              filter.overlapping.barcodes.VDJ,
              paste0(group.id, collapse = ";"))
  names(params) <- c("AIRR.input",
                     "VDJ.combine",
                     "filter.overlapping.barcodes.VDJ",
                     "group.id")

  if(verbose) message("\n AIRR tables loaded     ")
  if(verbose) message(Sys.time())

  stats.done <- FALSE
  if(get.VDJ.stats == TRUE){
    tryCatch({

      out.stats <- AIRR_VDJ_stats_int(airr.list = airr.list)
      stats.done <- TRUE
      if(verbose) message("\n Got VDJ stats    ")
      if(verbose) message(Sys.time())

    }, error = function(e){e
      if(verbose) message("\n VDJ stats failed: ")
      out.stats <- "failed"
      if(verbose)message(e)})
  } else {
    out.stats <- "none"
  }


  #similarly to VDJ_GEX_matrix
  barcodes_VDJ <- list()
  for(i in 1:length(airr.list)){
  barcodes_VDJ[[i]] <- unique(airr.list[[i]]$cell_id)

  if(verbose) message(paste0("\n For sample ", i, ": ", length(barcodes_VDJ[[i]])," cell assigned barcodes in VDJ"))
  }

  #remove sample overlapping barcodes in VDJ
  if(filter.overlapping.barcodes.VDJ == TRUE){
    if(length(barcodes_VDJ) > 1){
      barcodes_VDJ_c <- do.call("c", barcodes_VDJ)
      non_unique_barcodes <- names(table(barcodes_VDJ_c)[table(barcodes_VDJ_c) > 1])
      for(i in 1:length(barcodes_VDJ)){
        barcodes_VDJ[[i]] <- barcodes_VDJ[[i]][which(!barcodes_VDJ[[i]] %in% non_unique_barcodes)]
      }
      if(verbose) message(paste0("\n Removed a total of ", length(non_unique_barcodes), " cells with non unique barcodes"))
    }
  }

  #VDJ Processing per cell
    VDJ.proc.list <- list()
    for(i in 1:length(airr.list)){

      if(verbose) message(paste0("\n Starting VDJ barcode iteration ", i , " of ", length(airr.list), "...    "))
      if(verbose) message(Sys.time())

      if(trim.and.align){#validate that start end end columns of V segments are available and numeric

        if(any(!c("v_cigar", "j_cigar") %in% names(airr.list[[i]]))){
          trim.and.align <- FALSE
          warning('Columns "v_cigar", "j_cigar" for trimming information where not found in current airr input table. Setting trim.and.align to FALSE')
        }
        if(all(c("v_cigar", "j_cigar") %in% names(airr.list[[i]]))) {
          #use CIGAR string for trimming
          #Assume standard CIGAR output: ...S...M...S (first number indicates soft clipped nucleotides before alignment, second number indicates matching nuclotides, third number indicate soft clipped nucleotides at the end of the sequence)
          airr.list[[i]]$v_sequence_start <- as.numeric(gsub("S","",stringr::str_extract(airr.list[[i]]$v_cigar, "^\\d+S")))+1 #first soft clipped value
          #This returns NA if there is no clipped sequence and the alignment starts from the first nucleotide. In that case we need no trimming
          airr.list[[i]]$v_sequence_start[is.na(airr.list[[i]]$v_sequence_start)] <- 0

          airr.list[[i]]$j_sequence_end <- nchar(airr.list[[i]]$sequence) - as.numeric(gsub("S","",stringr::str_extract(airr.list[[i]]$j_cigar, "\\d+S$"))) #last soft clipped value
          #This returns NA if there is no clipped sequence at the end and the cigar ends in ...M In that case we do not trim from this side and set the end to the length of the contig
          airr.list[[i]]$j_sequence_end[is.na(airr.list[[i]]$j_sequence_end)] <- nchar(airr.list[[i]]$sequence)[is.na(airr.list[[i]]$j_sequence_end)]
        }
      }

      VDJ.proc.list[[i]] <- lapply(barcodes_VDJ[[i]], AIRR_barcode_VDJ_iteration, airrs = airr.list[[i]])

      #bind list recieved from parLapply
      VDJ.proc.list[[i]] <- dplyr::bind_rows(VDJ.proc.list[[i]])
      VDJ.proc.list[[i]][VDJ.proc.list[[i]] == ";"] <- "" #fix bug, where if two emtpy strings are concatenated, a ";" is left behind.
      #filter empty rows
      nv <- nrow(VDJ.proc.list[[i]])
      VDJ.proc.list[[i]] <- subset(VDJ.proc.list[[i]], Nr_of_VDJ_chains > 0 | Nr_of_VJ_chains > 0)
      if(verbose) message(paste0("\n Filtered out  ", nv -  nrow(VDJ.proc.list[[i]]), " cells with 0 chains"))

      #update barcodes
      VDJ.proc.list[[i]]$orig_barcode <-  VDJ.proc.list[[i]]$barcode
      VDJ.proc.list[[i]]$barcode <- paste0("s",i,"_",VDJ.proc.list[[i]]$barcode)
      VDJ.proc.list[[i]]$sample_id <- paste0("s",i)
      VDJ.proc.list[[i]]$group_id <- group.id[i]

      #add frequency column (i.e. all cells in clonotype2 will have the same entry, that is the number of cells in clonotype2)

      if("clonotype_id" %in% names(VDJ.proc.list[[i]])){
        VDJ.proc.list[[i]]$clonotype_id <- "Undetermined"
        VDJ.proc.list[[i]]$clonotype_frequency <- "Undetermined"
      }

      #Add further columns to fill in in future updates
      VDJ.proc.list[[i]]$specifity <- NA
      VDJ.proc.list[[i]]$affinity <- NA


      #! Add other columns from the input AIRR frame
      #In essence, a user may want to bring in metadata for features from the AIRR table into the VGM.
      #We give this option by automatically merging VJ and VDJ info for each cell into new columns of the VGM

      #remove columns that were already used and will therefore not be appended as extra data
      airr.list[[i]]$chain_for_merge <- substr(airr.list[[i]]$v_call, 1,3)
      airr.list[[i]] <- airr.list[[i]][,-c(which(names(airr.list[[i]]) %in% c("sequence", "sequence_aa", "sequence_id", "clone_id", "v_call", "j_call", "d_call","c_call","productive", "junction", "junction_aa", "sequence_alignment", "germline_alignment", "is_cell", "consensus_count", "complete_vdj")))]

      #VJ
      VJ_airr <- subset(airr.list[[i]], stringr::str_detect(airr.list[[i]]$chain_for_merge, "(TRD|TRA|IG(K|L))"))
      VJ_airr <- VJ_airr[,-c(which(names(VJ_airr) %in% c("chain_for_merge")))]
      #rename to fit VGM
      names(VJ_airr) <- paste0("VJ_", names(VJ_airr))
      #rename for merge
      names(VJ_airr)[which(names(VJ_airr) == "VJ_cell_id")] <- "orig_barcode"
      #merge in
      VDJ.proc.list[[i]] <- merge(VDJ.proc.list[[i]], VJ_airr, by = "orig_barcode", all.x = T, all.y = F)
      VJ_airr <- NULL


      #VDJ
      VDJ_airr <- subset(airr.list[[i]], stringr::str_detect(airr.list[[i]]$chain_for_merge, "(TRG|TRB|IGH)"))
      VDJ_airr <- VDJ_airr[,-c(which(names(VDJ_airr) %in% c("chain_for_merge")))]
      #rename to fit VGM
      names(VDJ_airr) <- paste0("VDJ_", names(VDJ_airr))
      #rename for merge
      names(VDJ_airr)[which(names(VDJ_airr) == "VDJ_cell_id")] <- "orig_barcode"
      #merge in
      VDJ.proc.list[[i]] <- merge(VDJ.proc.list[[i]], VDJ_airr, by = "orig_barcode", all.x = T, all.y = F)
      VDJ_airr <- NULL

      if(verbose) message(paste0("\n \t Done with ", i , " of ", length(airr.list), "     "))
      if(verbose) message(Sys.time())
    }

    VDJ.proc <- VDJ.proc.list
    VDJ.proc.list <- NULL
    if (VDJ.combine == TRUE){ #processing all VDJ files together
      tryCatch({

        for(i in 1:length(VDJ.proc)){ #setting all AIRR columns to character as this has caused issues in the past
          for(j in 41:ncol(VDJ.proc[[i]])){
            VDJ.proc[[i]][,j] <- as.character(VDJ.proc[[i]][,j])
          }
        }

        VDJ.proc.comb <- dplyr::bind_rows(VDJ.proc)
      }, error = function(e){
        message(e)
        message("\n dplyr::bind_rows failed, returning list of VDJ dataframes")
        VDJ.proc.comb <- VDJ.proc
    })
    } else {
      VDJ.proc.comb <- VDJ.proc
    }

    if(verbose) message("\n Done     ")
    if(verbose) message(Sys.time())

    return(list("VDJ" = VDJ.proc.comb,
                "GEX" = "none",
                "VDJ.GEX.stats" = out.stats,
                "Running params" = params,
                "sessionInfo" = utils::sessionInfo()))
}
