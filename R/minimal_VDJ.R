#' Minimal version of the VDJ building part from VDJ_GEX_matrix. Adapted for Cellranger v6 and v7.


#'@description  Minimal version of the VDJ building part from VDJ_GEX_matrix. Adapted for Cellranger v6 and v7. Currently, Seurat objects need to be integrated by matching barcodes from the Suerat object's metadata with the barcodes of the VDJ dataframe.
#' @param VDJ.directory.list List containing paths to VDJ output directories from cell ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires. ! Neccessary files within this folder: all_contig_annotations.csv, filtered_contig_annotations.csv, all_contig_annotations.csv, filtered_contig.fasta, consensus.fasta, concat_ref.fasta,
#' @param trim.germlines bool - if TRUE, will trim the germlines to the consensus sequences (using BIostrings::pairwiseAlignment with the option "global-local" and a gap opening cost = gap.opening.cost).
#' @param gap.opening.cost float or Inf - the cost for opening a gap in Biostrings::pairwiseAlignment when aligning and trimming germline sequences. Default to gapless alignment (gap.opening.cost = Inf).
#' @param parallel bool - if TRUE, will execute the per-sample VDJ building and germline trimming routines in parallel (parallelized across samples).
#' @param num.cores integer or NULL - number of cores to be used of parallel = T. Defaults to all ncores = (available cores - 1) or the ncores = number of samples (depending which is smaller).
#' @param operating.system string - operating system for choosing the parallel processing method. Options include: "Windows", "Darwin", "Linux".
#' @return Returns the VDJ dataframe/ VGM[[1]] object. Integration with the Seurat object needs to be performed separately (matching by barcodes).
#' @export
#' @examples
#' \dontrun{
#' VDJ <- minimal_VDJ(VDJ_directories, trim.germlines = TRUE,
#' gap.opening.cost = Inf, parallel = TRUE)
#'}


minimal_VDJ <- function(VDJ.directory.list,
                        trim.germlines,
                        gap.opening.cost,
                        parallel,
                        num.cores,
                        operating.system
                        ){

  if(missing(VDJ.directory.list)) stop('Input list of VDJ sample paths from cellranger v7')
  if(missing(trim.germlines)) trim.germlines <- FALSE
  if(missing(gap.opening.cost)) gap.opening.cost <- Inf
  if(missing(parallel)) parallel <- FALSE
  if(missing(num.cores)) num.cores <- NULL
  if(missing(operating.system)){
        switch(Sys.info()[['sysname']],
               Windows= {message("Windows system detected")
                         operating.system <- "Windows"},
               Linux  = {message("Linux system detected")
                        operating.system <- "Linux"},
               Darwin = {message("MAC system detected")
                        operating.system <- "Darwin"})
  }
  barcode <- NULL
  chain <- NULL
  raw_consensus_id <- NULL

  trim_ref <- function(consensus, reference, gap.opening.cost = Inf){
    #Trims the reference based on the first and last codon
    #Arguments:
    #         consensus: Consensus sequence for the clonotype and chain
    #         reference: Full reference sequence for the clonotype and chain
    #         gap.opening.cost: Cost for opening a gap in the alignment. For gapless set to Inf ##As suggested by Anamay Samant.
    #Author: Evgenios Kladis

    if (is.na(reference) | is.na(consensus)){
      return ("")
    }

    globalAlign <- Biostrings::pairwiseAlignment(Biostrings::DNAString(consensus), Biostrings::DNAString(reference), type="global-local", gapOpening = gap.opening.cost)
    return(as.character(globalAlign@subject))
  }

  single_sample_parse <- function(VDJ.directory){
    #Obtains all necessary data from a VDJ directory: filtered_contig_annotation.csv (most receptor info and trimmed sequences), filtered_contig.fasta (raw receptors), consensus_annotations.csv (trimmed consensus sequences), consensus.fasta (raw consensus), concat_ref.fasta (raw germlines).
    #Arguments:
    #         VDJ.directory: Path to a single sample in the VDJ directory from Cellranger v6 and v7
    #Author: Tudor-Stefan Cotet

    out <- vector("list", 5)

    out[[1]] <- utils::read.csv(paste(VDJ.directory,"/filtered_contig_annotations.csv",sep=""),sep=",",header=T)
    out[[1]]$sequence_nt_trimmed <- paste0(out[[1]][["fwr1_nt"]],
                                           out[[1]][["cdr1_nt"]],
                                           out[[1]][["fwr2_nt"]],
                                           out[[1]][["cdr2_nt"]],
                                           out[[1]][["fwr3_nt"]],
                                           out[[1]][["cdr3_nt"]],
                                           out[[1]][["fwr4_nt"]])

    out[[2]] <- data.frame(Biostrings::readDNAStringSet(paste(VDJ.directory,"/filtered_contig.fasta",sep="")))
    colnames(out[[2]]) <- 'sequence_nt_raw'
    out[[2]]$contig_id <- rownames(out[[2]])

    out[[3]] <- utils::read.csv(paste(VDJ.directory,"/consensus_annotations.csv",sep=""),sep=",",header=T)
    out[[3]]$consensus_nt_trimmed <- paste0(out[[3]][["fwr1_nt"]],
                                            out[[3]][["cdr1_nt"]],
                                            out[[3]][["fwr2_nt"]],
                                            out[[3]][["cdr2_nt"]],
                                            out[[3]][["fwr3_nt"]],
                                            out[[3]][["cdr3_nt"]],
                                            out[[3]][["fwr4_nt"]])
    out[[3]]$consensus_id <- unlist(lapply(out[[3]]$consensus_id, function(x) gsub(x, pattern = "(consensus)(.*)", replacement = "\\1_\\2")))

    out[[4]] <- data.frame(Biostrings::readDNAStringSet(paste(VDJ.directory,"/consensus.fasta",sep="")))
    colnames(out[[4]]) <- 'consensus_nt_raw'
    out[[4]]$raw_consensus_id <- rownames(out[[4]])

    out[[5]] <- data.frame(Biostrings::readDNAStringSet(paste(VDJ.directory,"/concat_ref.fasta",sep="")))
    colnames(out[[5]]) <- 'germline_nt_raw'
    out[[5]]$raw_consensus_id <- rownames(out[[5]])
    out[[5]]$raw_consensus_id <- unlist(lapply(out[[5]]$raw_consensus_id, function(x) gsub(x, pattern = "_concat_ref_", replacement = "_consensus_")))


    names(out) <- c('filtered_contig_table', 'raw_contigs_table', 'trimmed_consensus_table', 'raw_consensus_table', 'raw_reference_table')
    return(out)
  }

  make_vdj_sample <- function(index,
                              VDJ.directory.list,
                              trim.germlines = TRUE,
                              gap.opening.cost = Inf){

    #Creates a VDJ dataframe for a single sample.
    #Arguments:
    #         index: Sample/iteration index.
    #         VDJ.directoty.list: List of directories to all samples from Cellranger v6/v7.
    #         trim.germlines: if TRUE, will trim the germlines to the consensus sequences (using BIostrings::pairwiseAlignment with the option "global-local" and a gap opening cost = gap.opening.cost).
    #         gap.opening.cost: the cost for opening a gap in Biostrings::pairwiseAlignment when aligning and trimming germline sequences. Default to gapless alignment (gap.opening.cost = Inf). Suggested by Anamay Samant.
    #Author: Tudor-Stefan Cotet, Aurora Desideri Perea, Anamay Samant


    #Obtain final column names
    colnames_unordered <- c('barcode','sample_id','VDJ_fwr1_aa','VDJ_fwr1_nt','VDJ_cdr1_aa','VDJ_cdr1_nt','VDJ_fwr2_aa','VDJ_fwr2_nt','VDJ_cdr2_aa','VDJ_cdr2_nt','VDJ_fwr3_aa','VDJ_fwr3_nt','VDJ_cdr3s_aa','VDJ_cdr3s_nt','VDJ_fwr4_aa','VDJ_fwr4_nt','VDJ_umis','VDJ_vgene','VDJ_dgene','VDJ_jgene','VDJ_cgene','VDJ_sequence_nt_trimmed','VDJ_consensus_nt_trimmed','VDJ_consensus_nt_raw','VDJ_germline_nt_raw','VDJ_germline_nt_trimmed','VDJ_sequence_aa_trimmed','VDJ_consensus_aa_trimmed','VDJ_consensus_aa_raw','VDJ_germline_aa_raw','VDJ_germline_aa_trimmed','VDJ_chain','VDJ_raw_consensus_id','VDJ_contig_id','VJ_fwr1_aa','VJ_fwr1_nt','VJ_cdr1_aa','VJ_cdr1_nt','VJ_fwr2_aa','VJ_fwr2_nt','VJ_cdr2_aa','VJ_cdr2_nt','VJ_fwr3_aa','VJ_fwr3_nt','VJ_cdr3s_aa','VJ_cdr3s_nt','VJ_fwr4_aa','VJ_fwr4_nt','VJ_umis','VJ_vgene','VJ_jgene','VJ_cgene','VJ_sequence_nt_trimmed','VJ_consensus_nt_trimmed','VJ_consensus_nt_raw','VJ_germline_nt_raw','VJ_germline_nt_trimmed','VJ_sequence_aa_trimmed','VJ_consensus_aa_trimmed','VJ_consensus_aa_raw','VJ_germline_aa_raw','VJ_germline_aa_trimmed','VJ_chain','VJ_raw_consensus_id','celltype','Nr_of_VDJ_chains','Nr_of_VJ_chains','clonotype_id','clonotype_frequency','clonotype_id_10x')

    #Obtain the final column names in a similar order to the previous VGM[[1]] objects.
    colnames_reordered <- c('barcode','sample_id','clonotype_id','clonotype_frequency','celltype','VDJ_chain','VJ_chain', 'Nr_of_VDJ_chains','Nr_of_VJ_chains', 'VDJ_fwr1_aa','VDJ_fwr1_nt','VDJ_cdr1_aa','VDJ_cdr1_nt','VDJ_fwr2_aa','VDJ_fwr2_nt','VDJ_cdr2_aa','VDJ_cdr2_nt','VDJ_fwr3_aa','VDJ_fwr3_nt','VDJ_cdr3s_aa','VDJ_cdr3s_nt','VDJ_fwr4_aa','VDJ_fwr4_nt','VDJ_umis','VDJ_vgene','VDJ_dgene','VDJ_jgene','VDJ_cgene','VDJ_sequence_nt_trimmed','VDJ_consensus_nt_trimmed','VDJ_consensus_nt_raw','VDJ_germline_nt_raw','VDJ_germline_nt_trimmed','VDJ_sequence_aa_trimmed','VDJ_consensus_aa_trimmed','VDJ_consensus_aa_raw','VDJ_germline_aa_raw','VDJ_germline_aa_trimmed','VDJ_raw_consensus_id','VDJ_contig_id','VJ_fwr1_aa','VJ_fwr1_nt','VJ_cdr1_aa','VJ_cdr1_nt','VJ_fwr2_aa','VJ_fwr2_nt','VJ_cdr2_aa','VJ_cdr2_nt','VJ_fwr3_aa','VJ_fwr3_nt','VJ_cdr3s_aa','VJ_cdr3s_nt','VJ_fwr4_aa','VJ_fwr4_nt','VJ_umis','VJ_vgene','VJ_jgene','VJ_cgene','VJ_sequence_nt_trimmed','VJ_consensus_nt_trimmed','VJ_consensus_nt_raw','VJ_germline_nt_raw','VJ_germline_nt_trimmed','VJ_sequence_aa_trimmed','VJ_consensus_aa_trimmed','VJ_consensus_aa_raw','VJ_germline_aa_raw','VJ_germline_aa_trimmed','VJ_raw_consensus_id','clonotype_id_10x')

    #Read and parse a single sample directory using single_sample_parse()
    VDJ.directory <- VDJ.directory.list[[index]]
    sample.dataframes <- single_sample_parse(VDJ.directory)

    #Creates and empty VDJ with all unique barcodes from the filtered_contigs dataframe
    unique_barcodes <- unique(sample.dataframes$filtered_contig_table$barcode)
    VDJ <- data.frame(matrix(nrow = length(unique_barcodes)))
    VDJ$barcode <- unique_barcodes
    VDJ$sample_id <- paste0("s", index)

    #Appednds useful information from the raw contig, raw consensus, trimmed consensus, and raw germlines dataframes
    sample.dataframes$filtered_contig_table <- sample.dataframes$filtered_contig_table |>
            dplyr::left_join(sample.dataframes$raw_contigs_table, by="contig_id") |>
            dplyr::left_join(sample.dataframes$trimmed_consensus_table[, c("consensus_id", "consensus_nt_trimmed")], by = c("raw_consensus_id" = "consensus_id")) |>
            dplyr::left_join(sample.dataframes$raw_consensus_table, by="raw_consensus_id") |>
            dplyr::left_join(sample.dataframes$raw_reference_table, by="raw_consensus_id")

    #Remove appended dataframes from memory  - will save time and RAM
    #sample.dataframes$raw_contigs_table <- NULL
    #sample.dataframes$trimmed_consensus_table <- NULL
    #sample.dataframes$raw_consensus_table <- NULL
    #sample.dataframes$raw_reference_table <- NULL

    #If trim.germlines is set to T, will trim and align the raw references to the consensus sequences to obtain trimmed germlines
    if(trim.germlines){
      germline_df <- sample.dataframes$filtered_contig_table[!is.na(sample.dataframes$filtered_contig_table$consensus_nt_trimmed),]
      germline_df <- germline_df |> dplyr::distinct(raw_consensus_id, .keep_all = TRUE)

      germline_df$germline_nt_trimmed <- mapply(function(x,y) trim_ref(x, y, gap.opening.cost), germline_df$consensus_nt_trimmed, germline_df$germline_nt_raw)
      sample.dataframes$filtered_contig_table <- sample.dataframes$filtered_contig_table |>
               dplyr::left_join(germline_df[, c('raw_consensus_id', 'germline_nt_trimmed')], by="raw_consensus_id")
    }else{
      sample.dataframes$filtered_contig_table$germline_nt_trimmed <- NA
    }


    #Adds sequences as amino acids as well
    sequence_columns <- c('sequence_nt_trimmed','sequence_nt_raw','consensus_nt_trimmed','consensus_nt_raw','germline_nt_trimmed','germline_nt_raw')
    for(col in sequence_columns){
      new_col <- stringr::str_replace(col, '_nt_', '_aa_')
      sample.dataframes$filtered_contig_table[[new_col]] <- ifelse(is.na(sample.dataframes$filtered_contig_table[[col]]), NA, as.character(suppressWarnings(Biostrings::translate(Biostrings::DNAStringSet(sample.dataframes$filtered_contig_table[[col]])))))
    }

    #Separately match VDJ/heavy chains
    contigs_vdj <- subset(sample.dataframes$filtered_contig_table, chain %in% c("IGH","TRG","TRGB")) |>
         dplyr::select("fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
                "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
                "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
                "fwr4", "fwr4_nt", "umis",
                "v_gene", "d_gene", "j_gene", "c_gene",
                "sequence_nt_trimmed", "sequence_nt_trimmed",
                "consensus_nt_trimmed", "consensus_nt_raw",
                "germline_nt_raw", "germline_nt_trimmed",
                "sequence_aa_trimmed", "sequence_aa_trimmed",
                "consensus_aa_trimmed", "consensus_aa_raw",
                "germline_aa_raw", "germline_aa_trimmed",
                "chain",
                "raw_clonotype_id", "raw_consensus_id", "contig_id", "barcode")
    contigs_vdj <- contigs_vdj[order(contigs_vdj$umis, decreasing = TRUE),]
    colnames(contigs_vdj)[1:33] <- paste0('VDJ_', colnames(contigs_vdj)[1:33])

    #Separately match VJ/light chains
    contigs_vj <- subset(sample.dataframes$filtered_contig_table, chain %in% c('TRD','TRA','IGK','IGL')) |>
          dplyr::select("fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
                 "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
                 "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
                 "fwr4", "fwr4_nt", "umis",
                 "v_gene", "j_gene", "c_gene",
                 "sequence_nt_trimmed", "sequence_nt_trimmed",
                 "consensus_nt_trimmed", "consensus_nt_raw",
                 "germline_nt_raw", "germline_nt_trimmed",
                 "sequence_aa_trimmed", "sequence_aa_trimmed",
                 "consensus_aa_trimmed", "consensus_aa_raw",
                 "germline_aa_raw", "germline_aa_trimmed",
                 "chain",
                 "raw_clonotype_id", "raw_consensus_id", "contig_id", "barcode")
    contigs_vj <- contigs_vj[order(contigs_vj$umis, decreasing = TRUE),]
    colnames(contigs_vj)[1:32] <- paste0('VJ_', colnames(contigs_vj)[1:32])

    #Filter out aberrant cells with 2 or more chains by selecting the ones with highest UMIs.
    VDJ <- VDJ |>
           dplyr::left_join(contigs_vdj, by="barcode") |>
           dplyr::distinct(barcode, .keep_all = TRUE) |>
           dplyr::left_join(contigs_vj, by="barcode") |>
           dplyr::distinct(barcode, .keep_all = TRUE)

    VDJ <- VDJ[, 2:(length(VDJ)-1)]

    #Appends cell types
    VDJ$celltype[stringr::str_detect(paste0(VDJ$VDJ_chain,VDJ$VJ_chain), "TR")] <- "T cell"
    VDJ$celltype[stringr::str_detect(paste0(VDJ$VDJ_chain,VDJ$VJ_chain), "IG")] <- "B cell"

    #Add number of VDJ and VJ chains
    VDJ$Nr_of_VDJ_chains <- ifelse(is.na(VDJ$VDJ_chain), 0, 1)
    VDJ$Nr_of_VJ_chains <- ifelse(is.na(VDJ$VJ_chain), 0, 1)

    #Appends sample IDs to the barcode, as done in the previous versions of VGM[[1]]/VDJ_GEX_matrix
    VDJ$barcode <- paste0(VDJ$sample_id, rep('_', nrow(VDJ)), VDJ$barcode)

    #Adds clonotype IDs from the VDJ and VJ raw clonotype IDs (can be found in filtered_contig_annotations.csv)
    VDJ$clonotype_id <- VDJ$VDJ_raw_clonotype_id
    VDJ[is.na(VDJ$clonotype_id),]$clonotype_id <- VDJ[is.na(VDJ$clonotype_id),]$VJ_raw_clonotype_id

    #Clonotype counts/frequencies
    VDJ$clonotype_frequency <- unlist(lapply(VDJ$clonotype_id, function(x) length(which(VDJ$clonotype_id == x))))

    #Added to match the previous versions, but clonotype_id_10x = clonotype_id.
    VDJ$clonotype_id_10x <- VDJ$clonotype_id

    #Removes unwanted columns
    VDJ$VDJ_raw_clonotype_id <- NULL
    VDJ$VJ_raw_clonotype_id <- NULL

    #Predefined column names and order matching the old versions of VGM[[1]]. Also reorders by clonotype IDs
    colnames(VDJ) <- colnames_unordered
    VDJ <- VDJ[order(nchar(VDJ$clonotype_id), VDJ$clonotype_id), colnames_reordered]

    #Removes row index
    rownames(VDJ) <- NULL

    return(VDJ)
  }

  #Partial function definition for lapply
  partial_function <- function(x) {make_vdj_sample(x,
                                                   VDJ.directory.list = VDJ.directory.list,
                                                   trim.germlines = trim.germlines,
                                                   gap.opening.cost = gap.opening.cost
                                                   )}

  #Parallel execution via mclapply (Linux or Mac), parLapply (Windows).
  if(parallel){
    if(is.null(num.cores)){
      cores <- parallel::detectCores() - 1
      cores <- min(c(cores, length(VDJ.directory.list)))
    }else{
      cores <- num.cores
    }

    if(operating.system %in% c('Linux', 'Darwin')){
      VDJ_list <- parallel::mclapply(1:length(VDJ.directory.list), partial_function, mc.cores = cores)

    }else{
      cl <- parallel::makeCluster(cores)
      pdb_list <- parallel::parLapply(cl, 1:length(VDJ.directory.list), partial_function)
      parallel::stopCluster(cl)
    }
  }else{
      VDJ_list <- lapply(1:length(VDJ.directory.list), function(x) partial_function(x))
  }

  #Rbinds all per-sample dataframes into a single VDJ
  VDJ <- do.call("rbind", VDJ_list)

  return(VDJ)
}
