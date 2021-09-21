#' Extracts information on the VDJRegion level using MiXCR. This function assumes the user can run an executable instance of MiXCR and is elgible to use MiXCR as determined by license agreements. The VDJRegion corresponds to the recombined heavy and light chain loci starting from framework region 1 (FR1) and extending to frame work region 4 (FR4). This can be useful for extracting full-length sequences ready to clone and further calculating somatic hypermutation occurances.
#' @param VDJ.per.clone The output from the VDJ_per_clone function. This object should have information regarding the contigs and clonotype_ids for each cell.
#' @param mixcr.directory The directory containing an executable version of MiXCR. This must be downloaded separately and is under a separate license.
#' @param species Either "mmu" for mouse or "hsa" for human. These use the default germline genes for both species contained in MIXCR.
#' @return Returns a nested list containing VDJRegion information as determined by MIXCR. The outer list corresponds to the individual repertoires in the same structure as the input  VDJ.per.clone. The inner list corresponds to each clonal family, as determined by either the VDJ_clonotype function or the defaul nucleotide clonotyping produced by cellranger.Each element in the inner list corresponds to a dataframe containing repertoire information such as isotype, CDR sequences, mean number of UMIs. This output can be supplied to further package functions such as VDJ_extract_sequences and VDJ_GEX_integrate.
#' @seealso VDJ_extract_sequences
#' @export
#' @examples
#' \dontrun{
#'call_MIXCR(VDJ.per.clone = VDJ.per.clone.output
#',mixcr.directory = "~/Downloads/mixcr-3.0.12/mixcr",species = "mmu")
#'}
call_MIXCR <- function(VDJ.per.clone,mixcr.directory,species){

  if(missing(mixcr.directory)) print("No mixcr directory supplied. Assuming working directory holds mixcr executable")
  ### need to also read in the fastas
  for(i in 1:length(VDJ.per.clone)){
    temp.seq.hc <- list()
    temp.seq.lc <- list()
    temp.name.hc <- list()
    temp.name.lc <- list()
    for(j in 1:length(VDJ.per.clone[[i]])){
      temp.seq.hc[[j]] <- VDJ.per.clone[[i]][[j]]$full_HC_sequence ## now need to check this sequence is correct
      temp.seq.lc[[j]] <- VDJ.per.clone[[i]][[j]]$full_LC_sequence
      temp.name.hc[[j]] <- VDJ.per.clone[[i]][[j]]$contig_id_hc
      temp.name.lc[[j]] <- VDJ.per.clone[[i]][[j]]$contig_id_lc
    }
    temp.seq.hc_unlist <- unlist(temp.seq.hc)
    temp.seq.lc_unlist <- unlist(temp.seq.lc)
    temp.name.hc_unlist <- unlist(temp.name.hc)
    temp.name.lc_unlist <- unlist(temp.name.lc)

    seqinr::write.fasta(sequences = as.list(temp.seq.hc_unlist),names = temp.name.hc_unlist,file.out = "temphc.fasta")
    system(paste(mixcr.directory," align -OsaveOriginalReads=true -s ", species," temphc.fasta tempmixcrhc.out.vdjca",sep=""))
    system(paste(mixcr.directory," exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrhc.out.vdjca tempmixcrhc.out.txt",sep=""))
    seqinr::write.fasta(sequences = as.list(temp.seq.lc_unlist),names = temp.name.lc_unlist,file.out = "templc.fasta")
    system(paste(mixcr.directory," align -OsaveOriginalReads=true -s ", species," templc.fasta tempmixcrlc.out.vdjca",sep=""))
    system(paste(mixcr.directory," exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrlc.out.vdjca tempmixcrlc.out.txt",sep=""))
    temp.mixcr.hc <- utils::read.table(file ="tempmixcrhc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
    temp.mixcr.lc <- utils::read.table(file ="tempmixcrlc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)

    system("rm temphc.fasta")
    system("rm templc.fasta")
    system("rm tempmixcrhc.out.vdjca")
    system("rm tempmixcrlc.out.vdjca")
    system("rm tempmixcrhc.out.txt")
    system("rm tempmixcrlc.out.txt")

    ## now need to fill VDJ_per_clone

      for(j in 1:length(VDJ.per.clone[[i]])){
        VDJ.per.clone[[i]][[j]]$FRH1.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH2.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH3.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH4.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH1.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH2.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH3.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL1.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL2.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL3.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL4.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL1.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL2.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL3.NT <- rep("",nrow(VDJ.per.clone[[i]][[j]]))

        VDJ.per.clone[[i]][[j]]$FRH1.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH2.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH3.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRH4.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH1.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH2.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRH3.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL1.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL2.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL3.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$FRL4.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL1.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL2.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$CDRL3.AA <- rep("",nrow(VDJ.per.clone[[i]][[j]]))

        VDJ.per.clone[[i]][[j]]$bestVHlignment <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$bestVLAlignment <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$bestJHAlignment <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$bestJLAlignment <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$isotype <- rep("",nrow(VDJ.per.clone[[i]][[j]]))

        VDJ.per.clone[[i]][[j]]$AA.Jmutations.hc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$AA.Jmutations.lc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$NT.Jmutations.hc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$NT.Jmutations.lc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))

        VDJ.per.clone[[i]][[j]]$AA.Vmutations.hc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$AA.Vmutations.lc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$NT.Vmutations.hc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))
        VDJ.per.clone[[i]][[j]]$NT.Vmutations.lc <- rep("",nrow(VDJ.per.clone[[i]][[j]]))


        for(k in 1:nrow(VDJ.per.clone[[i]][[j]])){
          tryCatch({
          index.hc <- which(temp.mixcr.hc$descrsR1==VDJ.per.clone[[i]][[j]]$contig_id_hc[k])
          index.lc <- which(temp.mixcr.lc$descrsR1==VDJ.per.clone[[i]][[j]]$contig_id_lc[k])
          VDJ.per.clone[[i]][[j]]$FRH1.NT[k] <- temp.mixcr.hc$nSeqFR1[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH2.NT[k] <- temp.mixcr.hc$nSeqFR2[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH3.NT[k] <- temp.mixcr.hc$nSeqFR3[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH4.NT[k] <- temp.mixcr.hc$nSeqFR4[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH1.NT[k] <- temp.mixcr.hc$nSeqCDR1[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH2.NT[k] <- temp.mixcr.hc$nSeqCDR2[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH3.NT[k] <- temp.mixcr.hc$nSeqCDR3[index.hc]
          VDJ.per.clone[[i]][[j]]$FRL1.NT[k] <- temp.mixcr.lc$nSeqFR1[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL2.NT[k] <- temp.mixcr.lc$nSeqFR2[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL3.NT[k] <- temp.mixcr.lc$nSeqFR3[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL4.NT[k] <- temp.mixcr.lc$nSeqFR4[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL1.NT[k] <- temp.mixcr.lc$nSeqCDR1[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL2.NT[k] <- temp.mixcr.lc$nSeqCDR2[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL3.NT[k] <- temp.mixcr.lc$nSeqCDR3[index.lc]

          VDJ.per.clone[[i]][[j]]$FRH1.AA[k] <- temp.mixcr.hc$aaSeqFR1[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH2.AA[k] <- temp.mixcr.hc$aaSeqFR2[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH3.AA[k] <- temp.mixcr.hc$aaSeqFR3[index.hc]
          VDJ.per.clone[[i]][[j]]$FRH4.AA[k] <- temp.mixcr.hc$aaSeqFR4[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH1.AA[k] <- temp.mixcr.hc$aaSeqCDR1[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH2.AA[k] <- temp.mixcr.hc$aaSeqCDR2[index.hc]
          VDJ.per.clone[[i]][[j]]$CDRH3.AA[k] <- temp.mixcr.hc$aaSeqCDR3[index.hc]
          VDJ.per.clone[[i]][[j]]$FRL1.AA[k] <- temp.mixcr.lc$aaSeqFR1[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL2.AA[k] <- temp.mixcr.lc$aaSeqFR2[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL3.AA[k] <- temp.mixcr.lc$aaSeqFR3[index.lc]
          VDJ.per.clone[[i]][[j]]$FRL4.AA[k] <- temp.mixcr.lc$aaSeqFR4[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL1.AA[k] <- temp.mixcr.lc$aaSeqCDR1[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL2.AA[k] <- temp.mixcr.lc$aaSeqCDR2[index.lc]
          VDJ.per.clone[[i]][[j]]$CDRL3.AA[k] <- temp.mixcr.lc$aaSeqCDR3[index.lc]

          VDJ.per.clone[[i]][[j]]$bestVHAlignment[k] <- temp.mixcr.hc$bestVAlignment[index.hc]
          VDJ.per.clone[[i]][[j]]$bestVLAlignment[k] <- temp.mixcr.lc$bestVAlignment[index.lc]
          VDJ.per.clone[[i]][[j]]$bestJHAlignment[k] <- temp.mixcr.hc$bestJAlignment[index.hc]
          VDJ.per.clone[[i]][[j]]$bestJLAlignment[k] <- temp.mixcr.lc$bestJAlignment[index.lc]
          VDJ.per.clone[[i]][[j]]$isotype[k] <- temp.mixcr.hc$bestCAlignment[index.hc]

          VDJ.per.clone[[i]][[j]]$AA.Vmutations.hc[k] <- temp.mixcr.hc$aaMutationsVRegion[index.hc]
          VDJ.per.clone[[i]][[j]]$AA.Vmutations.lc[k] <- temp.mixcr.lc$aaMutationsVRegion[index.lc]
          VDJ.per.clone[[i]][[j]]$NT.Vmutations.hc[k] <- temp.mixcr.hc$nMutationsVRegion[index.hc]
          VDJ.per.clone[[i]][[j]]$NT.Vmutations.lc[k] <- temp.mixcr.lc$nMutationsVRegion[index.lc]

          VDJ.per.clone[[i]][[j]]$AA.Jmutations.hc[k] <- temp.mixcr.hc$aaMutationsJRegion[index.hc]
          VDJ.per.clone[[i]][[j]]$AA.Jmutations.lc[k] <- temp.mixcr.lc$aaMutationsJRegion[index.lc]
          VDJ.per.clone[[i]][[j]]$NT.Jmutations.hc[k] <- temp.mixcr.hc$nMutationsJRegion[index.hc]
          VDJ.per.clone[[i]][[j]]$NT.Jmutations.lc[k] <- temp.mixcr.lc$nMutationsJRegion[index.lc]


          }, error=function(e){})
        }#k loop
        VDJ.per.clone[[i]][[j]]$VDJ.AA.HC <- paste(VDJ.per.clone[[i]][[j]]$FRH1.AA,VDJ.per.clone[[i]][[j]]$CDRH1.AA,VDJ.per.clone[[i]][[j]]$FRH2.AA,VDJ.per.clone[[i]][[j]]$CDRH2.AA,VDJ.per.clone[[i]][[j]]$FRH3.AA,VDJ.per.clone[[i]][[j]]$CDRH3.AA,VDJ.per.clone[[i]][[j]]$FRH4.AA,sep="")
        VDJ.per.clone[[i]][[j]]$VDJ.AA.LC <- paste(VDJ.per.clone[[i]][[j]]$FRL1.AA,VDJ.per.clone[[i]][[j]]$CDRL1.AA,VDJ.per.clone[[i]][[j]]$FRL2.AA,VDJ.per.clone[[i]][[j]]$CDRL2.AA,VDJ.per.clone[[i]][[j]]$FRL3.AA,VDJ.per.clone[[i]][[j]]$CDRL3.AA,VDJ.per.clone[[i]][[j]]$FRL4.AA,sep="")
        VDJ.per.clone[[i]][[j]]$VDJ.NT.HC <- paste(VDJ.per.clone[[i]][[j]]$FRH1.NT,VDJ.per.clone[[i]][[j]]$CDRH1.NT,VDJ.per.clone[[i]][[j]]$FRH2.NT,VDJ.per.clone[[i]][[j]]$CDRH2.NT,VDJ.per.clone[[i]][[j]]$FRH3.NT,VDJ.per.clone[[i]][[j]]$CDRH3.NT,VDJ.per.clone[[i]][[j]]$FRH4.NT,sep="")
        VDJ.per.clone[[i]][[j]]$VDJ.NT.LC <- paste(VDJ.per.clone[[i]][[j]]$FRL1.NT,VDJ.per.clone[[i]][[j]]$CDRL1.NT,VDJ.per.clone[[i]][[j]]$FRL2.NT,VDJ.per.clone[[i]][[j]]$CDRL2.NT,VDJ.per.clone[[i]][[j]]$FRL3.NT,VDJ.per.clone[[i]][[j]]$CDRL3.NT,VDJ.per.clone[[i]][[j]]$FRL4.NT,sep="")
        VDJ.per.clone[[i]][[j]]$VDJ.NT.HC.LC <- paste(VDJ.per.clone[[i]][[j]]$VDJ.NT.HC,VDJ.per.clone[[i]][[j]]$VDJ.NT.LC,sep="_")
        VDJ.per.clone[[i]][[j]]$VDJ.AA.HC.LC <- paste(VDJ.per.clone[[i]][[j]]$VDJ.AA.HC,VDJ.per.clone[[i]][[j]]$VDJ.AA.LC,sep="_")

      }
  }
  return(VDJ.per.clone)
}


