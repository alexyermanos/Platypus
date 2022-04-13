#' Platypus V2 utility for full germline sequence via MiXCR
#'
#'@description Only Platypus v2. Extracts the full-length germline sequence as determined by cellranger. This function returns an object that now contains the reference germline for each of the clones. If multiple clones (as determined by cellranger) have been merged using the VDJ_clonotype function then these sequences may have distinct germline sequences despite being in the same clonal family (nested list). This is particularly possible when homology thresholds were used to determine the clonotypes.
#' @param VDJ.per.clone The output from the VDJ_per_clone function. This object should have information regarding the contigs and clonotype_ids for each cell.
#' @param mixcr.directory  The directory containing an executable version of MiXCR. This must be downloaded separately and is under a separate license.
#' @param extract.VDJRegion Default is TRUE. Future iterations will allow for distinct gene regions to be extracted.
#' @param species Either "mus" or "hsa" for mouse and human respectively. Default is set to mouse.
#' @return Returns a dataframe containing repertoire information, such as isotype, CDR sequences, mean number of UMIs. This output can be supplied to furhter packages VDJ_extract_sequences and VDJ_GEX_integrate
#' @export
#' @examples
#' \dontrun{
#' VDJ_extract_germline(VDJ.per.clone=VDJ.per.clone.output
#' ,mixcr.directory="~/Downloads/mixcr-3.0.12/mixcr"
#' ,extract.VDJRegion=T,species = "mmu")
#'}
#'
VDJ_extract_germline <- function(VDJ.per.clone,
                                 mixcr.directory,
                                 extract.VDJRegion,species){
  output.list <- list()
  sequence_output_list <- list()
  if(missing(extract.VDJRegion)) extract.VDJRegion <- T
  if(missing(species)) species <- "mus"

  for(i in 1:length(VDJ.per.clone)){
    sequence_output_list[[i]] <- list()

    clone_number <- length(VDJ.per.clone[[i]])
    output.list[[i]] <- data.frame(clonotype_id=rep("",clone_number),
                                   full_HC_germline=rep("",clone_number),
                                   full_HC_germline=rep("",clone_number),stringsAsFactors = F)
    for(j in 1:length(VDJ.per.clone[[i]])){
      output.list[[i]]$full_HC_germline[j] <- names(which.max(table(VDJ.per.clone[[i]][[j]]$full_HC_germline)))
      output.list[[i]]$full_LC_germline[j] <- names(which.max(table(VDJ.per.clone[[i]][[j]]$full_LC_germline)))
      output.list[[i]]$clonotype_id[j] <- names(which.max(table(VDJ.per.clone[[i]][[j]]$clonotype_id)))
    }

    if(extract.VDJRegion==T){
      seqinr::write.fasta(sequences = as.list(output.list[[i]]$full_HC_germline),names = output.list[[i]]$clonotype_id,file.out = "temphc.fasta")
      system(paste(mixcr.directory," align -OsaveOriginalReads=true -s ", species," temphc.fasta tempmixcrhc.out.vdjca",sep=""))
      system(paste(mixcr.directory," exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrhc.out.vdjca tempmixcrhc.out.txt",sep=""))
      seqinr::write.fasta(sequences = as.list(output.list[[i]]$full_HC_germline),names = output.list[[i]]$clonotype_id,file.out = "templc.fasta")
      system(paste(mixcr.directory," align -OsaveOriginalReads=true -s ", species," templc.fasta tempmixcrlc.out.vdjca",sep=""))
      system(paste(mixcr.directory," exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrlc.out.vdjca tempmixcrlc.out.txt",sep=""))
      temp.mixcr.hc <- utils::read.table(file ="tempmixcrhc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
      temp.mixcr.lc <- utils::read.table(file ="tempmixcrlc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
      sequence_output_list[[i]][[1]] <- temp.mixcr.hc
      sequence_output_list[[i]][[2]] <- temp.mixcr.lc

      sequence_output_list[[i]][[1]]$VDJ.NT.HC <- paste(sequence_output_list[[i]][[1]]$nSeqFR1,
                                                        sequence_output_list[[i]][[1]]$nSeqCDR1,
                                                        sequence_output_list[[i]][[1]]$nSeqFR2,
                                                        sequence_output_list[[i]][[1]]$nSeqCDR2,
                                                        sequence_output_list[[i]][[1]]$nSeqFR3,
                                                        sequence_output_list[[i]][[1]]$nSeqCDR3,
                                                        sequence_output_list[[i]][[1]]$nSeqFR4,sep="")
      sequence_output_list[[i]][[1]]$VDJ.AA.HC <- paste(sequence_output_list[[i]][[1]]$aaSeqFR1,
                                                        sequence_output_list[[i]][[1]]$aaSeqCDR1,
                                                        sequence_output_list[[i]][[1]]$aaSeqFR2,
                                                        sequence_output_list[[i]][[1]]$aaSeqCDR2,
                                                        sequence_output_list[[i]][[1]]$aaSeqFR3,
                                                        sequence_output_list[[i]][[1]]$aaSeqCDR3,
                                                        sequence_output_list[[i]][[1]]$aaSeqFR4,sep="")
      sequence_output_list[[i]][[1]]$VDJ.NT.LC <- paste(sequence_output_list[[i]][[2]]$nSeqFR1,
                                                        sequence_output_list[[i]][[2]]$nSeqCDR1,
                                                        sequence_output_list[[i]][[2]]$nSeqFR2,
                                                        sequence_output_list[[i]][[2]]$nSeqCDR2,
                                                        sequence_output_list[[i]][[2]]$nSeqFR3,
                                                        sequence_output_list[[i]][[2]]$nSeqCDR3,
                                                        sequence_output_list[[i]][[2]]$nSeqFR4,sep="")
      sequence_output_list[[i]][[1]]$VDJ.AA.LC <- paste(sequence_output_list[[i]][[2]]$aaSeqFR1,
                                                        sequence_output_list[[i]][[2]]$aaSeqCDR1,
                                                        sequence_output_list[[i]][[2]]$aaSeqFR2,
                                                        sequence_output_list[[i]][[2]]$aaSeqCDR2,
                                                        sequence_output_list[[i]][[2]]$aaSeqFR3,
                                                        sequence_output_list[[i]][[2]]$aaSeqCDR3,
                                                        sequence_output_list[[i]][[2]]$aaSeqFR4,sep="")
      sequence_output_list[[i]][[1]]$VDJ.AA.HC.LC <- paste(sequence_output_list[[i]][[1]]$VDJ.AA.HC,sequence_output_list[[i]][[1]]$VDJ.AA.LC,sep="")
      sequence_output_list[[i]][[1]]$VDJ.NT.HC.LC <- paste(sequence_output_list[[i]][[1]]$VDJ.NT.HC,sequence_output_list[[i]][[1]]$VDJ.NT.LC,sep="")

      system("rm temphc.fasta")
      system("rm templc.fasta")
      system("rm tempmixcrhc.out.vdjca")
      system("rm tempmixcrlc.out.vdjca")
      system("rm tempmixcrhc.out.txt")
      system("rm tempmixcrlc.out.txt")
    }
  }

  return(sequence_output_list)
}
