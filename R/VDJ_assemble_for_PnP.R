#' Assembles sequences from MIXCR output into inserts for expression in PnP cells. ! ALWAYS VALIDATE INDIVIDUAL SEQUENCE IN GENEIOUS OR OTHER SOFTWARE BEFORE ORDERING ! Check notes on column content below ! Only cells with 1 VDJ and 1 VJ sequence are considered. Warnings are issued if sequences do not pass necessary checks
#' @param VDJ.mixcr.matrix Output dataframe from the VDJ_call_MIXCR function or a dataframe generated using the VDJ_GEX_matrix function and supplemented with MIXCR information (Needed columns: All Framework and CDR sequences)
#' @param id.column Character. Column name of VDJ.mixcr.matrix to use as ID for the assembled sequences. Defaults to "barcode"
#' @param species Character. Which IgKC sequence to use. Can be "human" or "mouse". Defaults to "mouse"
#' @param manual_IgKC Character. Manual overwrite for sequence used as IgKC.
#' @param manual_2A Character. Manual overwrite for sequence used as Furine 2A site.
#' @param manual_VDJLeader Character. Manual overwrite for sequence used as VDJ Leader and signal peptide.
#' @param write.to.disk Boolean. Defaults to TRUE. Whether to save assembled sequences to working directory
#' @param filename Character. Output file name for .fasta and .csv files if write.to.disk == T. Defaults to PnP_assembled_seqs.fasta/.csv
#' @return Returns the input VGM matrix with one additional column containing the assembles sequences. If write.to.disk == T writes a CSV containing key columns of the VGM as well as a .FASTA file to the current working director (getwd())
#' ! Important notes on column content:
#' 1. The column "seq_length_check" contains either "passed" or "FAILED". If FAILED, this means that at least one of the sequences (e.g. FRL1) was shorter than 9NTs and therefore considered invalid. Please check for missing sequences if you find any warnings
#' 2. The column "seq_codon_check" is deemed "passed" if all CDR and FR input sequences of a cell contain only full codons (i.e. are divisible by 3)
#' 3. The column "PnP_assembled_seqs" contains the assembled sequences / inserts for PnP expression. These should be validated manually in Geneious or other software and can then be ordered to be synthesized.
#' 4. The column "PnP_assembled_annotations" contains a string of annotations for the respective assembled sequence. The structure is | [Sequence element] -> [index (starting from 1) of last nucleotide of the sequence element] ...
#' 5. The column "PnP_assembled_translations" contains the amino acid translation of the full contig that will result from the assembled insert in the backbone PnP vector. Please note: the sequences in the PnP_assembled_translation resulted from pasting the VJ leader sequence (contained in the PnP vector backbone), the PnP_assembled_seqs (The insert itself) and a surrogate stop codon ATAA. If correct, the translation should only contain one * (stop codon) at the very end. For reference: VJLeader sequence: ATGGATTTTCAGGTGCAGATTTTCAGCTTCCTGCTAATCAGCGCTTCAGTTATAATGTCCCGGGGG
#' 6. The column "seq_VJCDR3_check" is deemed "passed" if the translated sequence of the input VJ CDR3 is found in the translated assembled sequence. If this test fails, there is likely an issue with the VJ segment
#' 7. The column "seq_Fur2A_check" is deemed "passed" if correct AA sequence of the 2A site is found in the translated assembled sequence. If this test fails, and the seq_VJCDR3_test was passed, there is likely an issue at the border between VJ and IgKC/2A sequences
#' 8. The column "seq_VDJCDR3_check" is deemed "passed" if the translated sequence of the input VDJ CDR3 is found in the translated assembled sequence.
#' 9. The column "seq_splicesite_check" is deemed passed if the last 6 nucleotides of the assembled sequence are one of the following: "TCCTCA", "TCTTCA","TCGTCA","TCATCA".
#' @export
#' @examples
#' \dontrun{
#'
#' VGM_with_PnP_seq <- VDJ_assemble_for_PnP(VDJ.mixcr.matrix = VDJ_call_MIXCR.output
#' , id.column = "barcode",species = "mouse", manual_IgKC = "none", manual_2A = "none"
#' , manual_VDJLeader = "none", write.to.disk = T, filename = "PnP_seq_example")
#'
#'}

VDJ_assemble_for_PnP <- function(VDJ.mixcr.matrix,
                                 id.column,
                                 species,
                                 manual_IgKC,
                                 manual_2A,
                                 manual_VDJLeader,
                                 write.to.disk,
                                 filename){


  platypus.version <- "v3"
  if(missing(species)) species <- "mouse"
  if(missing(manual_IgKC)) manual_IgKC <- "none"
  if(missing(manual_2A)) manual_2A <- "none"
  if(missing(manual_VDJLeader)) manual_VDJLeader <- "none"
  if(missing(write.to.disk)) write.to.disk <- T
  if(missing(id.column)) id.column <- "barcode"
  if(missing(filename)) filename <- "PnP_assembled_seqs"

  VDJ.matrix <- VDJ.mixcr.matrix

  if(all(c("Nr_of_VJ_chains", "Nr_of_VDJ_chains") %in% names(VDJ.matrix))){
    cat("\n Excluded cells with more or less than 1 VDJ 1 VJ chain")
    VDJ.matrix <- subset(VDJ.matrix, Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
  }

  #check data
  if(any(!c("VDJ_nSeqFR1", "VDJ_nSeqFR2","VDJ_nSeqFR3","VDJ_nSeqFR4","VDJ_nSeqCDR1","VDJ_nSeqCDR2","VDJ_nSeqCDR3","VJ_nSeqFR1", "VJ_nSeqFR2","VJ_nSeqFR3","VJ_nSeqFR4","VJ_nSeqCDR1","VJ_nSeqCDR2","VJ_nSeqCDR3") %in% names(VDJ.matrix))){
    stop(paste0("At least one of the neccessary columns in the input matrix are missing. Neccessary columns are: VDJ_nSeqFR1, VDJ_nSeqFR2,VDJ_nSeqFR3,VDJ_nSeqFR4,VDJ_nSeqCDR1,VDJ_nSeqCDR2,VDJ_nSeqCDR3,VJ_nSeqFR1, VJ_nSeqFR2,VJ_nSeqFR3,VJ_nSeqFR4,VJ_nSeqCDR1,VJ_nSeqCDR2,VJ_nSeqCDR3"))
  }
  if(!id.column %in% names(VDJ.matrix)) stop("id.column not found in input matrix")

  #choose static sequences
  if(manual_IgKC == "none"){
    if(species == "human"){ IgKC <- "CGAACTGTGGCTGCACCATCTGTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACTGCCTCTGTTGTGTGCCTGCTGAATAACTTCTATCCCAGAGAGGCCAAAGTACAGTGGAAGGTGGATAACGCCCTCCAATCGGGTAACTCCCAGGAGAGTGTCACAGAGCAGGACAGCAAGGACAGCACCTACAGCCTCAGCAGCACCCTGACGCTGAGCAAAGCAGACTACGAGAAACACAAAGTCTACGCCTGCGAAGTCACCCATCAGGGCCTGAGCTCGCCCGTCACAAAGAGCTTCAACAGGGGAGAGTGT"
    } else if(species == "mouse") { IgKC <- "CGGGCCGACGCGGCCCCAACTGTATCCATCTTCCCACCATCCAGTGAGCAGTTAACATCTGGAGGTGCCTCAGTCGTGTGCTTCTTGAACAACTTCTACCCCAAAGACATCAATGTCAAGTGGAAGATTGATGGCAGTGAACGACAAAATGGCGTCCTGAACAGTTGGACTGATCAGGACAGCAAAGACAGCACCTACAGCATGAGCAGCACCCTCACGTTGACCAAGGACGAGTATGAACGACATAACAGCTATACCTGTGAGGCCACTCACAAGACATCAACTTCACCCATTGTCAAGAGCTTCAACAGGAATGAGTGT"
    } else {
      stop("Please enter either mouse or human as the species parameter")
    }
  } else{
    IgKC <- manual_IgKC
  }

  if(manual_2A == "none"){

    Fur_2A <- "AGGAAAAGACGACACAAACAGAAAATTGTGGCACCGGTGAAACAGACTTTGAATTTTGACCTTCTCAAGTTGGCGGGAGACGTCGAGTCCAACCCTGGGCCC"
  } else{
    Fur_2A <- manual_2A
  }

  if(manual_VDJLeader == "none"){
    VDJLeader <- "ATGATGGTGTTAAGTCTTCTGTACCTGTTGACAGCACTTCCGGGTGAGTGTTTCCATTTCATACATGTGCCATGAGGATTTTTCAAAATGTGTGATTGACAGATTTGATTCTTTTTGTCTAAAGGTATCCTGTCA"
  } else{
    VDJLeader <- manual_VDJLeader
  }

  #Needed for contig check later
  VJLeader <- "ATGGATTTTCAGGTGCAGATTTTCAGCTTCCTGCTAATCAGCGCTTCAGTTATAATGTCCCGGGGG"

  cat("\n Got sequences; starting assembly...")

  seq_frame <- VDJ.matrix[,which(names(VDJ.matrix) %in% c(id.column, "VDJ_nSeqFR1", "VDJ_nSeqFR2","VDJ_nSeqFR3","VDJ_nSeqFR4","VDJ_nSeqCDR1","VDJ_nSeqCDR2","VDJ_nSeqCDR3","VJ_nSeqFR1", "VJ_nSeqFR2","VJ_nSeqFR3","VJ_nSeqFR4","VJ_nSeqCDR1","VJ_nSeqCDR2","VJ_nSeqCDR3"))]

  #check all input sequences for missing ones or too short ones (< 9NTs)

  seq_length_check <- NULL
  seq_frame$seq_length_check <- "passed"

  for(i in 1:nrow(seq_frame)){
    if(any(nchar(seq_frame[i,2:ncol(seq_frame)]) < 6)){
      warning(paste0("At least one FR or CDR3 seq of cell id ",seq_frame[i,1]), " is less than 9 nt long. Please check for possible issues or missing sequences")
      seq_frame$seq_length_check[i] <- "FAILED"
    }
  }

  seq_codon_check <- NULL
  seq_frame$seq_codon_check <- "passed"

  for(i in 1:nrow(seq_frame)){
    if(sum(nchar(seq_frame[i,2:ncol(seq_frame)]) %% 3 != 0) > 2){ #This returns true only if more than two sequences are not divisible by three. Reason: The VDJ_FR4 and VJ_FR4 by MIXCR always contain one extra nucleotide at the end. This is trimmed of in the next section
      warning(paste0("At least one FR or CDR3 seq of cell id ",seq_frame[i,1]), " does contain partial codons (i.e. sequence length not divisible by 3")
      seq_frame$seq_codon_check[i] <- "FAILED"
    }
  }

  #trim last nucleotide of relevant sequences
  VJ_nSeqFR4.lastnttrimmed <- NULL
  VDJ_nSeqFR4.lastnttrimmed <- NULL
  seq_frame$VJ_nSeqFR4.lastnttrimmed <- substr(seq_frame$VJ_nSeqFR4, start = 0, stop = nchar(seq_frame$VJ_nSeqFR4)-1)
  seq_frame$VDJ_nSeqFR4.lastnttrimmed <- substr(seq_frame$VDJ_nSeqFR4, start = 0, stop = nchar(seq_frame$VDJ_nSeqFR4)-1)

  #paste sequences together
  #ORDER: FRL1 CDRL1 FRL2 CDRL2 FRL3 CDRL3 FRL4_lastnttrimmed IgKC Fur_2A VDJLeader FRH1 CDRH1 FRH2 CDRH2 FRH3 CDRH3 FRH4_lastnttrimmed

  pasted_seqs <- paste0(seq_frame$VJ_nSeqFR1,seq_frame$VJ_nSeqCDR1,seq_frame$VJ_nSeqFR2,seq_frame$VJ_nSeqCDR2,seq_frame$VJ_nSeqFR3,seq_frame$VJ_nSeqCDR3,seq_frame$VJ_nSeqFR4.lastnttrimmed,IgKC, Fur_2A,VDJLeader,seq_frame$VDJ_nSeqFR1,seq_frame$VDJ_nSeqCDR1,seq_frame$VDJ_nSeqFR2,seq_frame$VDJ_nSeqCDR2,seq_frame$VDJ_nSeqFR3,seq_frame$VDJ_nSeqCDR3,seq_frame$VDJ_nSeqFR4.lastnttrimmed)
#
  #get annotations to make checking easier later
  nchar_frame <- seq_frame
  nchar_frame$IgKC <- IgKC
  nchar_frame$Fur_2A <- Fur_2A
  nchar_frame$VDJLeader <- VDJLeader

  #order
  nchar_frame <- nchar_frame[,c("VJ_nSeqFR1","VJ_nSeqCDR1", "VJ_nSeqFR2","VJ_nSeqCDR2","VJ_nSeqFR3","VJ_nSeqCDR3","VJ_nSeqFR4.lastnttrimmed","IgKC","Fur_2A","VDJLeader","VDJ_nSeqFR1","VDJ_nSeqCDR1", "VDJ_nSeqFR2","VDJ_nSeqCDR2","VDJ_nSeqFR3","VDJ_nSeqCDR3","VDJ_nSeqFR4.lastnttrimmed")]

  for(i in 1:nrow(nchar_frame)){ #getting length of each element and summing these cumulatively
    nchar_frame[i,] <- cumsum(nchar(nchar_frame[i,]))
  }

  pasted_annotations <- paste0(
  "VJ_FR1 -> ",nchar_frame$VJ_nSeqFR1," | VJ_CDRL -> ",nchar_frame$VJ_nSeqCDR1,
  " | VJ_FR2 -> ",nchar_frame$VJ_nSeqFR2," | VJ_CDR2 -> ",nchar_frame$VJ_nSeqCDR2,
  " | VJ_FR3 -> ",nchar_frame$VJ_nSeqFR3," | VJ_CDR3 -> ",nchar_frame$VJ_nSeqCDR3,
  " | VJ_FRL4.lastnttrimmed -> ",nchar_frame$VJ_nSeqFR4.lastnttrimmed,
  " | IgKC -> ",nchar_frame$IgKC, " | Fur 2A -> ",nchar_frame$Fur_2A,
  " | VDJLeader -> ",nchar_frame$VDJLeader,
  " | VDJ_FRH1 -> ",nchar_frame$VDJ_nSeqFR1," | VDJ_CDR1 -> ",nchar_frame$VDJ_nSeqCDR1,
  " | VDJ_FR2 -> ",nchar_frame$VDJ_nSeqFR2," | VDJ_CDR2 -> ",nchar_frame$VDJ_nSeqCDR2,
  " | VDJ_FR3 -> ",nchar_frame$VDJ_nSeqFR3," | VDJ_CDR3 -> ",nchar_frame$VDJ_nSeqCDR3,
  " | VDJ_FR4.lastnttrimmed -> ",nchar_frame$VDJ_nSeqFR4.lastnttrimmed)

  PnP_assembled_seqs <- NULL
  PnP_assembled_annotations <- NULL
  seq_frame$PnP_assembled_seqs <- pasted_seqs
  seq_frame$PnP_assembled_annotations <- pasted_annotations

  #Check sequences for contig integrity
  #This is a bit tricky, because these sequences do not start with a start codon (this is located in the vector backbone VJ leader sequence)
  #we therefore append the pasted sequences to that leader peptide and than check for contig integrity
  leader_seqs_pasted <- paste0(VJLeader, seq_frame$PnP_assembled_seqs)

  #now translate these sequences
  trans_test <- as.character(Biostrings::translate(Biostrings::DNAStringSet(leader_seqs_pasted)))

  #for verification
  trans_VJ_CDR3s <- as.character(Biostrings::translate(Biostrings::DNAStringSet(seq_frame$VJ_nSeqCDR3)))
  trans_VDJ_CDR3s <- as.character(Biostrings::translate(Biostrings::DNAStringSet(seq_frame$VDJ_nSeqCDR3)))

  PnP_assembled_translations <- NULL
  seq_frame$PnP_assembled_translations <- trans_test

  cat("\n Please note: the sequences in the PnP_assembled_translations column resulted from pasting the VJ leader sequence (contained in the PnP vector backbone) and the PnP_assembled_seqs (The insert itself)")

  seq_VJCDR3_check <- NULL
  seq_frame$seq_VJCDR3_check <- "passed"
  seq_Fur2A_check <- NULL
  seq_frame$seq_Fur2A_check <- "passed"
  seq_VDJCDR3_check <- NULL
  seq_frame$seq_VDJCDR3_check <- "passed"
  seq_splicesite_check <- NULL
  seq_frame$seq_splicesite_check <- "passed"

  for(i in 1:nrow(seq_frame)){
    #Check for VJ CDR3 correct translation
    if(!stringr::str_detect(seq_frame$PnP_assembled_translations[i], pattern = trans_VJ_CDR3s[i])){
      warning(paste0("Correct VJ CDR3 AA seq not found in the assembled sequence for cell id ", seq_frame[i,1]), "!")
      seq_frame$seq_VJCDR3_check[i] <- "FAILED"
    }
    #Check for 2A correct translation
    if(!stringr::str_detect(seq_frame$PnP_assembled_translations[i], pattern = as.character(Biostrings::translate(Biostrings::DNAStringSet(Fur_2A))))){
      warning(paste0("Correct Furine 2A AA seq not found in the assembled sequence for cell id ", seq_frame[i,1]), "!")
      seq_frame$seq_Fur2A_check[i] <- "FAILED"
    }
    #Check for VDJ CDR3 correct translation
    if(!stringr::str_detect(seq_frame$PnP_assembled_translations[i], pattern = trans_VDJ_CDR3s[i])){
      warning(paste0("Correct VDJ CDR3 AA seq not found in the assembled sequence for cell id ", seq_frame[i,1]), "!")
      seq_frame$seq_VDJCDR3_check[i] <- "FAILED"
    }
    #Check for the presence of the correct two codons at the end of the VDJ FR4 region ! Using the untranslated sequence
    if(!substr(seq_frame$PnP_assembled_seqs[i], start = nchar(seq_frame$PnP_assembled_seqs[i])-5, stop = 10000) %in% c("TCCTCA", "TCTTCA","TCGTCA","TCATCA")){
      warning(paste0("Splicing site not found at end of VDJ FR4 in the assembled sequence for cell id ", seq_frame[i,1]), "! Valid splicing site requires a TCA at the sequence end")
      seq_frame$seq_splicesite_check[i] <- "FAILED"
    }
  }

  cat("\n Assembly and checks done")

  #Assemble VGM output
  cat("\n Adding additional columns to VDJ.mixcr.matrix input")
  VDJ.matrix <- cbind(VDJ.matrix, seq_frame)

  if(write.to.disk == F){
    cat("\n Done")
    return(VDJ.matrix)
  } else {
    cat("\n Building .FASTA and .csv file")
    fasta_names <- paste0("Seq ID: ", seq_frame[,1], " | seq_length_check: ", seq_frame$seq_length_check," | seq_codon_check: ", seq_frame$seq_codon_check, " | seq_VJCDR3_check: ", seq_frame$seq_VJCDR3_check," | seq_Fur2A_check: ", seq_frame$seq_Fur2A_check," | seq_VDJCDR3_check: ", seq_frame$seq_VDJCDR3_check," | seq_splicesite_check: ", seq_frame$seq_splicesite_check, " Annotations: ", seq_frame$PnP_assembled_annotations)

    seqinr::write.fasta(as.list(unlist(seq_frame$PnP_assembled_seqs)), names = fasta_names, file.out = paste0(filename,".fasta"))

    write.csv(seq_frame, file = paste0(filename,".csv"))

    cat("\n Done \n")
    return(VDJ.matrix)
  }
}
