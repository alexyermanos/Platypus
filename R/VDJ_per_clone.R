#' Analyzes and processes the repertoire sequencing data from cellranger vdj. This provides information on the single-cell level for each clone, as opposed to the output from VDJ_analyze.
#' @title VDJ_per_clone
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
#' check_VDJ_per_clone <- VDJ_per_clone(clonotype.list = output.from.VDJ_analyze, VDJ.out.directory = "path/to/cellranger/outs/")
#' }
VDJ_per_clone <- function(clonotype.list,
                               VDJ.out.directory,
                               contig.list,
                               fasta.list,
                               reference.list,
                               filtered.contigs,
                               annotations.json,
                               JSON) {
  require(seqinr)
  require(jsonlite)
  require(dplyr)
  require(msa)
  if(missing(VDJ.out.directory)) print("No output directory supplied. Assuming clonotype and contig are provided as list objects")
  if(missing(contig.list)) print("No contig.list supplied. Assuming contigs should be extracted from working directory")
  if(missing(filtered.contigs)) filtered.contigs <- TRUE
  if(missing(JSON)) JSON <- FALSE

  print("Reading in output files")
  ### need to also read in the fastas
  if(missing(VDJ.out.directory)==F){
    if(filtered.contigs==T)VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/filtered_contig_annotations.csv",sep="")
    if(filtered.contigs==F)VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/all_contig_annotations.csv",sep="")

    contig.list <- lapply(VDJ.out.directory_contigs, function(x) utils::read.csv(x, stringsAsFactors = FALSE,sep=",",header=T))
    VDJ.out.directory_fasta <- paste(VDJ.out.directory,"/filtered_contig.fasta",sep="")
    fasta.list <- lapply(VDJ.out.directory_fasta, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    VDJ.out.directory_reference <- paste(VDJ.out.directory,"/concat_ref.fasta",sep="")
    reference.list <- lapply(VDJ.out.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    if (JSON==TRUE) annotations.json <- read_json(paste0(VDJ.out.directory, "/all_contig_annotations.json"), simplifyVector = T)
  }

  VDJ.per.cell <- list()
  ### each element in VDJ.per.cell is a dataframe containing per cell information
  tryCatch({

    for(i in 1:length(clonotype.list)){
      VDJ.per.cell[[i]] <- list()
      temp_trb <- which(contig.list[[i]]$chain=="TRB")
      temp_tra <- which(contig.list[[i]]$chain=="TRA")
      contig.list[[i]]$chain[temp_trb] <- "IGH"
      contig.list[[i]]$chain[temp_tra] <- "IGK"
      holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "=",
                                           width = NA, style = 1, file = "")
      print(paste(i,"from",length(clonotype.list),"repertoires"))
      for(j in 1:nrow(clonotype.list[[i]])){
        utils::setTxtProgressBar(value = j/nrow(clonotype.list[[i]]),pb = holding_bar)

        temp.cell.number <- clonotype.list[[i]]$frequency[j]
        print(temp.cell.number)
        temp_concat_ref_HC <- gsub(pattern = "_consensus_",replacement = "_concat_ref_",x = contig.list[[i]]$raw_consensus_id[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH")][1])
        temp_concat_ref_LC <- gsub(pattern = "_consensus_",replacement = "_concat_ref_",x = contig.list[[i]]$raw_consensus_id[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGK")][1])
        temp_concat_ref_HC_seq <- as.character(reference.list[[i]][which(names(reference.list[[i]])==temp_concat_ref_HC)])
        temp_concat_ref_LC_seq <- as.character(reference.list[[i]][which(names(reference.list[[i]])==temp_concat_ref_LC)])

        VDJ.per.cell[[i]][[j]] <- data.frame(barcode=rep("",temp.cell.number),clonotype_id=rep(clonotype.list[[i]]$clonotype_id[j],temp.cell.number),isotype_hc=rep("",temp.cell.number),isotype_lc=rep("",temp.cell.number),
                                             HC_vgene=rep("",temp.cell.number),HC_jgene=rep("",temp.cell.number),LC_vgene=rep("",temp.cell.number),LC_jgene=rep("",temp.cell.number),full_HC_sequence=rep("",temp.cell.number),
                                             full_LC_sequence=rep("",temp.cell.number),full_HC_germline=rep(temp_concat_ref_HC_seq,temp.cell.number),full_LC_germline=rep(temp_concat_ref_LC_seq,temp.cell.number),
                                             trimmed_LC_sequence=rep("",temp.cell.number), trimmed_HC_sequence=rep("",temp.cell.number),trimmed_HC_germline=rep("",temp.cell.number), trimmed_LC_germline=rep("",temp.cell.number),
                                             umi_count=rep("",temp.cell.number),contig_id_hc=rep("",temp.cell.number),contig_id_lc=rep("",temp.cell.number),stringsAsFactors = F)
        temp_str_split <- stringr::str_split(string = clonotype.list[[i]]$nt_clone_ids[j],pattern = ";")[[1]]
        temp_barcode_list <- list()
        for(k in 1:length(temp_str_split)){
          temp_barcode_list[[k]] <- as.character(unique(contig.list[[i]]$barcode[which(contig.list[[i]]$raw_clonotype_id==temp_str_split[k] & str_detect(contig.list[[i]]$productive, "(?i)true"))]))
        }
        VDJ.per.cell[[i]][[j]]$barcode <- unlist(temp_barcode_list)

        #### now need to be flexible with the barcodes
        for(k in 1:length(VDJ.per.cell[[i]][[j]]$barcode)){

          VDJ.per.cell[[i]][[j]]$isotype_hc[k] <- names(which.max(table(contig.list[[i]]$c_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$isotype_lc[k] <- names(which.max(table(contig.list[[i]]$c_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$umi_count[k] <- names(which.max(table(contig.list[[i]]$umis[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$contig_id_lc[k] <- names(which.max(table(contig.list[[i]]$contig_id[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$contig_id_hc[k] <- names(which.max(table(contig.list[[i]]$contig_id[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))

          VDJ.per.cell[[i]][[j]]$HC_vgene[k] <- names(which.max(table(contig.list[[i]]$v_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$HC_jgene[k] <- names(which.max(table(contig.list[[i]]$j_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$LC_vgene[k] <- names(which.max(table(contig.list[[i]]$v_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))
          VDJ.per.cell[[i]][[j]]$LC_jgene[k] <- names(which.max(table(contig.list[[i]]$j_gene[which(str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH" & contig.list[[i]]$productive=="True" & contig.list[[i]]$barcode==VDJ.per.cell[[i]][[j]]$barcode[k])])))

          VDJ.per.cell[[i]][[j]]$full_HC_sequence[k] <- as.character(fasta.list[[i]][which(names(fasta.list[[i]])==VDJ.per.cell[[i]][[j]]$contig_id_hc[k])])
          VDJ.per.cell[[i]][[j]]$full_LC_sequence[k] <- as.character(fasta.list[[i]][which(names(fasta.list[[i]])==VDJ.per.cell[[i]][[j]]$contig_id_lc[k])])

          if(JSON==TRUE){

            # Adding VDJ trimmed sequence
            selected_contig <- annotations.json[annotations.json$contig_name == VDJ.per.cell[[i]][[j]]$contig_id_hc[k],]
            info <- selected_contig$annotations[[1]]
            vdj_seq_start <- info$contig_match_start[info$feature$region_type == "L-REGION+V-REGION"]
            vdj_seq_end <- info$contig_match_end[info$feature$region_type == "J-REGION"]
            seq_trimmed <- substr(selected_contig$sequence, vdj_seq_start, vdj_seq_end)
            ref_trimmed <- as.character(pairwiseAlignment(seq_trimmed, VDJ.per.cell[[i]][[j]]$full_HC_sequence[1], type = "local")@subject)
            VDJ.per.cell[[i]][[j]]$trimmed_HC_sequence[k] <- seq_trimmed
            VDJ.per.cell[[i]][[j]]$trimmed_HC_germline[k] <- ref_trimmed

            # Same for light chain
            selected_contig <- annotations.json[annotations.json$contig_name == VDJ.per.cell[[i]][[j]]$contig_id_lc[k],]
            info <- selected_contig$annotations[[1]]
            vdj_seq_start <- info$contig_match_start[info$feature$region_type == "L-REGION+V-REGION"]
            vdj_seq_end <- info$contig_match_end[info$feature$region_type == "J-REGION"]
            seq_trimmed <- substr(selected_contig$sequence, vdj_seq_start, vdj_seq_end)
            ref_trimmmed <- as.character(pairwiseAlignment(seq_trimmed, VDJ.per.cell[[i]][[j]]$full_LC_sequence[1], type = "local")@subject)
            VDJ.per.cell[[i]][[j]]$trimmed_LC_sequence[k] <- seq_trimmed
            VDJ.per.cell[[i]][[j]]$trimmed_LC_germline[k] <- ref_trimmmed
          }
        }
      }
    }
  }, error=function(e){})
  return(VDJ.per.cell)
}
