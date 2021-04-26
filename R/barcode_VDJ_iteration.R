
#' Helper function called in VDJ_GEX_matrix. Do not run as standalone!
#' @param All_PARAMS Inputs provided by VDJ_GEX_matrix
#' @return List of tables, with each table containing information one one barcode
#' @export
#' @examples
#' \dontrun{
#' Do not run as standalone!
#' }

#FUN to call in parlapply mclapply or lapply
barcode_VDJ_iteration <- function(barcodes, contigs, references, annotations, gap.opening.cost, gap.extension.cost,trim.and.align){
  
  require(stringr)
  require(Biostrings)
  
  #Get all the info needed to shrink data usage and search times later in the function
  #Filtering out non productive or non full length contigs from cell. This is neccessary, as a cell labled as productive and full length may still have associated contigs not fullfilling these criteria.
  curr.contigs <- contigs[which(contigs$barcode == barcodes & tolower(contigs$is_cell) == "true" & tolower(contigs$high_confidence) == "true" & tolower(contigs$productive) == "true" & tolower(contigs$full_length) == "true"),]
  if(curr.contigs$raw_clonotype_id[1] != ''){
    curr.references <- references[which(str_detect(names(references), curr.contigs$raw_clonotype_id[1]))]} else {curr.references <- ""}
  
  #getting the relevant annotations
  curr.annotations <- annotations[str_detect(annotations$contig_id, barcodes),]
  
  #set up data structure
  cols <- c("barcode","sample_id","group_id","clonotype_id_10x","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa", "VJ_cdr3s_aa","VDJ_cdr3s_nt", "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain","VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw","VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_trimmed_ref","VJ_trimmed_ref")
  curr.barcode <- stats::setNames(data.frame(matrix(ncol = length(cols), nrow = 1)), cols)
  
  #fill in information that do not need processing
  #Contig info on light/a and heavy/b chains is put into different columns (see cols)
  #If only one contig is available, the fields of the other are left blank
  #If more than two contigs of one chain (e.g. 2 TRB) are present, the elements will be pasted separated by a ";" into the relevant fields (in the case of TRB, into the Hb columns)
  
  #open temporary data format. This is used as a bridge between raw contigs and the curr.barcode table. It provides the flexibility needed to process both cells with 0, 1 or more Heavy/b or Light/a chains
  contigs_pasted <- setNames(data.frame(matrix(ncol = ncol(curr.contigs), nrow = length(unique(curr.contigs$chain)))), names(curr.contigs)) #the dataframe may be one or two rows too long, this will not matter / ROW 1 = Heavy chain information / ROW 2 = Light chain information. This order is maintained even if one of the two chains is not present!
  
  #Heavy/b chain count
  if(stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRB|IGH)") == 1){
    contigs_pasted[1,] <- curr.contigs[which(str_detect(curr.contigs$chain, pattern = "(TRB|IGH)") == T),]
  } else if(stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRB|IGH)") > 1){
    for(k in 1:ncol(curr.contigs)){
      contigs_pasted[1,k] <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$chain, pattern = "(TRB|IGH)") == T), k], collapse = ";")
    }
  } else {
    contigs_pasted[1,] <- ""
  }
  
  #Light/a chain count
  if(stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRA|IG(K|L))") == 1){
    contigs_pasted[2,] <- curr.contigs[which(stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))") == T),]
  } else if(stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRA|IG(K|L))") > 1){
    for(k in 1:ncol(curr.contigs)){
      contigs_pasted[2,k]  <- paste0(curr.contigs[which(stringr::str_detect(curr.contigs$chain, pattern = "(TRA|IG(K|L))") == T),k],collapse = ";")
    }
  } else {
    contigs_pasted[2,] <- ""
  }
  
  #fill in the pasted info to curr.barcode
  curr.barcode$barcode <- curr.contigs$barcode[1]
  curr.barcode$clonotype_id_10x <- curr.contigs$raw_clonotype_id[1]
  curr.barcode$sample_id <- curr.contigs$raw_clonotype_id[1]
  curr.barcode$group_id <- curr.contigs$raw_clonotype_id[1]
  if(stringr::str_detect(paste0(contigs_pasted$chain, collapse = " "), "TR") == T){
    curr.barcode$celltype <- "T cell"
  } else if(stringr::str_detect(paste0(contigs_pasted$chain, collapse = " "), "IG") == T) {
    curr.barcode$celltype <- "B cell"
  } else {
    curr.barcode$celltype <- "Unkown"
  }
  curr.barcode$Nr_of_VDJ_chains <- stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRB|IGH)")
  curr.barcode$Nr_of_VJ_chains <- stringr::str_count(paste0(curr.contigs$chain, collapse= " "), pattern = "(TRA|IG(K|L))")
  
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
  
  #Now on to the actual sequences
  reference_HC <- curr.references[which(names(curr.references) == gsub("_consensus_","_concat_ref_",curr.barcode$VDJ_raw_consensus_id))]
  reference_LC <- curr.references[which(names(curr.references) == gsub("_consensus_","_concat_ref_",curr.barcode$VJ_raw_consensus_id))]
  
  if(trim.and.align == T){
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
          #trim sequence
          to_paste_trimmed <- append(to_paste_trimmed, substr(curr.annotations$sequence[l], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
          #translate trimmed sequence
          if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
            to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
          } else {to_paste_aa <- ""}
          #align to reference and trim reference
          tryCatch({
            if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
              alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_HC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
              to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(subject(alignments[which.max(score(alignments))])))
            } else {
              to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
            }
          }, error=function(e){
            to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
          })
        }
      }
    }
    curr.barcode$VDJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
    curr.barcode$VDJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
    curr.barcode$VDJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
    curr.barcode$VDJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    
    #Light/a
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
          #trim sequence
          to_paste_trimmed <- append(to_paste_trimmed, substr(curr.annotations$sequence[l], as.numeric(curr.annotations$temp_start[l])+1, as.numeric(curr.annotations$temp_end[l])-1))
          #translate trimmed sequence
          if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
            to_paste_aa <- append(to_paste_aa, as.character(Biostrings::translate(DNAStringSet(to_paste_trimmed[length(to_paste_trimmed)]))))
          } else {to_paste_aa <- ""}
          #align to reference and trim reference
          tryCatch({
            if(nchar(to_paste_trimmed[length(to_paste_trimmed)]) > 1){
              alignments <- Biostrings::pairwiseAlignment(to_paste_trimmed[length(to_paste_trimmed)], as.character(reference_LC), type = "local", gapOpening = gap.opening.cost, gapExtension = gap.extension.cost)
              to_paste_ref_trimmed <- append(to_paste_ref_trimmed, as.character(subject(alignments[which.max(score(alignments))])))
            } else {
              to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "")
            }
          }, error=function(e){
            to_paste_ref_trimmed <- append(to_paste_ref_trimmed, "ALIGNMENT ERROR")
          })
        }
      }
    }
    curr.barcode$VJ_sequence_nt_raw <- paste0(to_paste, collapse = ";")
    curr.barcode$VJ_sequence_nt_trimmed <- paste0(to_paste_trimmed, collapse = ";")
    curr.barcode$VJ_sequence_aa <- paste0(to_paste_aa, collapse = ";")
    curr.barcode$VJ_trimmed_ref <- paste0(to_paste_ref_trimmed, collapse = ";")
    
  }else {
    curr.barcode$VDJ_sequence_nt_raw <- NA
    curr.barcode$VDJ_sequence_nt_trimmed <- NA
    curr.barcode$VDJ_sequence_aa <- NA
    curr.barcode$VDJ_trimmed_ref <- NA
    
    curr.barcode$VJ_sequence_nt_raw <- NA
    curr.barcode$VJ_sequence_nt_trimmed <- NA
    curr.barcode$VJ_sequence_aa <- NA
    curr.barcode$VJ_trimmed_ref <- NA
  }
  return(curr.barcode)
}
