#' Analyzes and processes the repertoire sequencing data from cellranger vdj. This provides information on the single-cell level for each clone, as opposed to the output from VDJ_analyze.
#' @title VDJ_per_cell
#' @param VDJ.out.directory Character vector with each element containing the path to the output of cellranger vdj runs. This corresponds to the same object used for the VDJ_analyze function. Multiple repertoires to be integrated in a single transcriptome should be supplied as multiple elements of the character vector. This can be left blank if supplying the clonotypes and contig files directly as input. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @param clonotype.list Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, the i'th entry in the directory supplied to GEX.list should correspond to the i'th element in the clonotype.list object.
#' @param contig.list List of dataframe based on the all_contigs.csv file from cellranger vdj output. If 10x sequencing was not used then this object should be formatted with the same columns as the 10x object.
#' @param fasta.list Contains the full-length sequence information in the same format as filtered_contig.fasta file from the output of cellranger.
#' @param reference.list Contains the reference sequence information in the same format as concat_ref.fasta file from the output of cellranger.
#' @param filtered.contigs Logical indicating if the filtered contigs file should be used. TRUE will read VDJ information from only the filtered output of cellranger. FALSE will read the all contigs file from cellranger. Default set to TRUE (filtered output)
#' @return Returns a list of dataframes containing
#' @export
#' @examples
#' \dontrun{
#' check_VDJ_per_cell <- VDJ_per_cell(clonotype.list = output.from.VDJ_analyze, VDJ.out.directory = "path/to/cellranger/outs/")
#' }
VDJ_per_cell <- function(clonotype.list,
         VDJ.out.directory,
         contig.list,
         fasta.list,
         reference.list,
         filtered.contigs,
         annotations.json) {
  require(seqinr)
  require(jsonlite)
  require(dplyr)
  require(msa)
  if(missing(VDJ.out.directory)) print("No output directory supplied. Assuming clonotype and contig are provided as list objects")
  if(missing(contig.list)) print("No contig.list supplied. Assuming contigs should be extracted from working directory")
  if(missing(filtered.contigs)) filtered.contigs <- TRUE

  print("Reading in output files")
  ### need to also read in the fastas
  if(missing(VDJ.out.directory)==F) {
    if(filtered.contigs==T) VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/filtered_contig_annotations.csv",sep="")
    if(filtered.contigs==F) VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/all_contig_annotations.csv",sep="")

    contig.list <- lapply(VDJ.out.directory_contigs, function(x) utils::read.csv(x, stringsAsFactors = FALSE,sep=",",header=T))
    VDJ.out.directory_fasta <- paste(VDJ.out.directory,"/filtered_contig.fasta",sep="")
    fasta.list <- lapply(VDJ.out.directory_fasta, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    VDJ.out.directory_reference <- paste(VDJ.out.directory,"/concat_ref.fasta",sep="")
    reference.list <- lapply(VDJ.out.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    annotations.json <- read_json(paste0(VDJ.out.directory, "/all_contig_annotations.json"), simplifyVector = T)
  }
  VDJ.per.cell <- list()
  ### each element in VDJ.per.cell is a dataframe containing per cell information
  for(i in 1:length(clonotype.list)) {
    VDJ.per.cell[[i]] <- list()
    temp_trb <- which(contig.list[[i]]$chain == "TRB")
    temp_tra <- which(contig.list[[i]]$chain == "TRA")
    contig.list[[i]]$chain[temp_trb] <- "IGH"
    contig.list[[i]]$chain[temp_tra] <- "IGK"
    holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "=", width = NA, style = 1, file = "")
    print(paste(i,"from",length(clonotype.list),"repertoires"))
    for(j in 1:nrow(clonotype.list[[i]])) {
      utils::setTxtProgressBar(value = j/nrow(clonotype.list[[i]]),pb = holding_bar)
      temp.cell.number <- clonotype.list[[i]]$frequency[j]
      temp_concat_ref_HC <- gsub(pattern = "_consensus_",replacement = "_concat_ref_",x = contig.list[[i]]$raw_consensus_id[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGH")][1])
      temp_concat_ref_LC <- gsub(pattern = "_consensus_",replacement = "_concat_ref_",x = contig.list[[i]]$raw_consensus_id[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGK")][1])
      ref_HC_seq <- as.character(reference.list[[i]][which(names(reference.list[[i]])==temp_concat_ref_HC)])
      ref_LC_seq <- as.character(reference.list[[i]][which(names(reference.list[[i]])==temp_concat_ref_LC)])

      # Adding VDJ trimmed sequence
      selected_contigs_HC <- annotations.json[which(annotations.json$info$raw_clonotype_id == clonotype.list[[i]]$clonotype_id[j]
                              & sapply(annotations.json$annotations, function(x) any(x$feature$chain == "IGH")) & annotations.json$productive & annotations.json$is_cell),]
      seq_trimmmed_HC <- substr(selected_contigs_HC$sequence,
                              sapply(selected_contigs_HC$annotations, function(x) x$contig_match_start[x$feature$region_type == "L-REGION+V-REGION"]),
                              sapply(selected_contigs_HC$annotations, function(x) x$contig_match_end[x$feature$region_type == "J-REGION"]))
      # Same for LC
      selected_contigs_LC <- annotations.json[which(annotations.json$info$raw_clonotype_id == clonotype.list[[i]]$clonotype_id[j]
                              & sapply(annotations.json$annotations, function(x) any(x$feature$chain == "IGK")) & annotations.json$productive & annotations.json$is_cell),]
      seq_trimmmed_LC <- substr(selected_contigs_LC$sequence,
                              sapply(selected_contigs_LC$annotations, function(x) x$contig_match_start[x$feature$region_type == "L-REGION+V-REGION"]),
                              sapply(selected_contigs_LC$annotations, function(x) x$contig_match_end[x$feature$region_type == "J-REGION"]))
      if (length(ref_HC_seq) == 0) {
        ref_trimmed_HC <- rep("", temp.cell.number)
        ref_HC_seq <- rep("", temp.cell.number)
        seq_trimmmed_HC <- rep("", temp.cell.number)
      }
      if(length(ref_LC_seq) == 0) {
        ref_trimmed_LC <- rep("", temp.cell.number)
        ref_LC_seq <- rep("", temp.cell.number)
        seq_trimmmed_LC <- rep("", temp.cell.number)
      } else {
        # HC germline trimmed sequence
        alignments <- pairwiseAlignment(unique(seq_trimmmed_HC), ref_HC_seq, type = "global-local")
        ref_trimmed_HC <- as.character(alignments@subject[which.max(alignments@score)])
        # LC germline trimmed sequence
        alignments <- pairwiseAlignment(unique(seq_trimmmed_LC), ref_LC_seq, type = "global-local")
        ref_trimmed_LC <- as.character(alignments@subject[which.max(alignments@score)])
      }
      empty <- rep("", temp.cell.number)
      VDJ.per.cell[[i]][[j]] <- data.frame(barcode=empty, clonotype_id=rep(clonotype.list[[i]]$clonotype_id[j],temp.cell.number),
                                           isotype_hc=empty,isotype_lc=empty, HC_vgene=empty,HC_jgene=empty,LC_vgene=empty,
                                           LC_jgene=empty,full_HC_sequence=empty, full_LC_sequence=empty,
                                           full_HC_germline= ref_HC_seq,
                                           full_LC_germline= ref_LC_seq,
                                           trimmed_HC_sequence=seq_trimmmed_HC,
                                           trimmed_LC_sequence=seq_trimmmed_LC,
                                           trimmed_HC_germline=ref_trimmed_HC,
                                           trimmed_LC_germline=ref_trimmed_LC,
                                           umi_count=empty, contig_id_hc=empty,
                                           contig_id_lc=empty,
                                           stringsAsFactors = F)
      temp_str_split <- stringr::str_split(string = clonotype.list[[i]]$nt_clone_ids[j],pattern = ";")[[1]]
      temp_barcode_list <- list()
      for(k in 1:length(temp_str_split)){
        temp_barcode_list[[k]] <- as.character(unique(contig.list[[i]]$barcode[which(contig.list[[i]]$raw_clonotype_id==temp_str_split[k]
                                                                                    & contig.list[[i]]$productive=="True")]))
      }
      VDJ.per.cell[[i]][[j]]$barcode <- unlist(temp_barcode_list)

      #### now need to be flexible with the barcodes
      for(k in 1:length(VDJ.per.cell[[i]][[j]]$barcode)){
        VDJ.per.cell[[i]][[j]]$isotype_hc[k] <- names(which.max(table(contig.list[[i]]$c_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$isotype_lc[k] <- names(which.max(table(contig.list[[i]]$c_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$umi_count[k] <- names(which.max(table(contig.list[[i]]$umis[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$contig_id_lc[k] <- names(which.max(table(contig.list[[i]]$contig_id[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$contig_id_hc[k] <- names(which.max(table(contig.list[[i]]$contig_id[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$HC_vgene[k] <- names(which.max(table(contig.list[[i]]$v_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$HC_jgene[k] <- names(which.max(table(contig.list[[i]]$j_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$LC_vgene[k] <- names(which.max(table(contig.list[[i]]$v_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$LC_jgene[k] <- names(which.max(table(contig.list[[i]]$j_gene[which(contig.list[[i]]$is_cell=="True" & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
        VDJ.per.cell[[i]][[j]]$full_HC_sequence[k] <- as.character(fasta.list[[i]][which(names(fasta.list[[i]])==VDJ.per.cell[[i]][[j]]$contig_id_hc[k])])
        VDJ.per.cell[[i]][[j]]$full_LC_sequence[k] <- as.character(fasta.list[[i]][which(names(fasta.list[[i]])==VDJ.per.cell[[i]][[j]]$contig_id_lc[k])])
      }
    }
  }
  return(VDJ.per.cell)
}
