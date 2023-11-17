#' Making the trimmed reference and concatenating fr1-fr4
#'
#' @description Function that takes the VDJ and the fr1-fr4 sequence per antibody
#' Based on the ref argument, if TRUE it also returns the returns in the VDJ/VJ_ref.nt/aa the trimmed reference based
#' on the alignement with the consensus.
#' @param  VDJ VDJ or vgm[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param n_clones integer, denoting the top n clones to get the reference. If NA it is performed in all clones
#' @param sample list of sample names, with the same order as they were accessed to make the VGM
#' @param ref bool, denoting whether or not we trim the reference of the antibodies.
#' @param path_tOData str, denoting the folder containing the VDJ folder with VDJ information per sample
#' @return $vdj: VDJ containing the VDJ/VJ_ref.nt/aa columns if ref = TRUE and the full_VDJ, full_VJ columns with the fr1-fr4. $clones: clone_ids for which a reference was made.
#' @examples
#' \dontrun{
#' samples = c('LCMV', 'TNFR')
#' vgm = read("VGM.RData")
#' n_clones = 20
#' result = VDJ_extract_germline_consensus_ref(vgm$VDJ, n_clones,
#' samples, ref = TRUE,
#' path_toData="../Data/")
#' VDJ = result[1]$vdj
#' clone_counts = result[2]$clones
#' }
#'
VDJ_extract_germline_consensus_ref<- function(VDJ, n_clones = NA, samples = NA, ref = TRUE, path_toData = "../Data/"){
  #SUPLEMENTARY FUNCTIONS
  #1
  trim_ref <- function(consensus, reference){
    #Trims the reference based on the first and last codon
    #Arguments:
    #         consensus: Consensus sequence for the clonotype and chain
    #         reference: Full reference sequence for the clonotype and chain

    if (is.na(reference) | is.na(consensus)){
      return ("")
    }
    #Since I make consensus with the concatenation of fr1-fr4 consensus, it shouldnt need to be trimmed since it is already at apropriate shape
    s1 <- DNAString(consensus)
    s2 <- DNAString(reference)
    globalAlign <- pairwiseAlignment(s1, s2, type="global-local",gapOpening = Inf)
    new_reference <- as.character(globalAlign@subject)
    return (new_reference)
  }
  #2
  ref_per_cell = function(row, consensus_data, examined_refs){
    #Functions that does the reference trimming for one cell on both chains
    #Arguemens:
    #         row: Vgm row for the cell to be examines
    #         consensus_data: The consensus information for the specific cell
    #         examined_refs: A names list with existed computed references per clone, so that we don't trim the same reference more than once.
    #Output
    #         output: A list with [[1]] containing the VDJ.nt VDJ.aa VJ.aa VJ.nt
    #                             [[2]] containing the examined_refs updated with new trimmed references



    #Subset the vgm on one specific cell and extract the barcode, to map the results to the original vgm

    vgm_barcode = row
    barcode = vgm_barcode["barcode"]
    #if(nrow(VDJ[VDJ$barcode == barcode,])!=0){
    output = list(VDJ_ref.aa = "", VDJ_ref.nt = "", VJ_ref.aa = "", VJ_ref.nt = "")
    for (chain in c("VDJ","VJ")){
      #Due to the issue with the consensus id of th VGM, take it from the reference id
      con_id = str_sub(vgm_barcode[[paste(chain,"_raw_consensus_id",sep="")]],-1,-1)
      clonotype_cons = vgm_barcode[[paste(chain,"_raw_consensus_id",sep="")]]
      #In case the cell has only light or only heavy chain, skip the process
      if(clonotype_cons == ""){
        reference_n = ""
      }
      else{
        if(!clonotype_cons %in% names(examined_refs)){
          #print(paste(clonotype_cons,names(examined_refs)))
          #Get the first sequence in the clonotype and start trimming the ref
          reference_n = vgm_barcode[[paste(chain, "_raw_ref", sep="")]]
          consensus_n = consensus_data[consensus_data["consensus_id"] == paste(clone, "_consensus",con_id,sep=""),]
          consensus_n = paste(consensus_n["fwr1_nt"], consensus_n["cdr1_nt"], consensus_n["fwr2_nt"], consensus_n["cdr2_nt"],
                              consensus_n["fwr3_nt"], consensus_n["cdr3_nt"], consensus_n["fwr4_nt"],sep = "")
          #Trim reference
          reference_n = trim_ref(consensus_n, reference_n)
          examined_refs[clonotype_cons] = reference_n
        }
        else{
          reference_n = examined_refs[[clonotype_cons]]
        }
        output[paste(chain,"_ref.nt",sep = "")] = reference_n
        output[paste(chain,"_ref.aa",sep = "")] = translate_DNA(reference_n)
        #VDJ[VDJ$barcode == barcode,paste(chain,"_ref_trimmed_nt",sep = "")] = reference_n
        #VDJ[VDJ$barcode == barcode,paste(chain,"_ref_trimmed_aa",sep = "")] = translate_DNA(reference_n)
      }
    }

    #}
    return(list(output, examined_refs))
  }


  #3
  translate_DNA<- function(sequence){
    #Gets a nucleotide sequence with also "-" and turns it into an aminonacid sequence
    if (sequence == ""){
      return("")
    }
    genetic_code <- list(
      "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
      "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
      "TAT"="Y", "TAC"="Y", "TAA"="*", "TAG"="*",
      "TGT"="C", "TGC"="C", "TGA"="*", "TGG"="W",
      "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
      "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
      "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
      "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
      "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
      "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
      "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
      "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
      "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
      "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
      "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
      "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
    )

    codons <- strsplit(sequence, "(?<=.{3})", perl=TRUE)[[1]] #Break into codons

    for (codon_id in 1:length(codons)){
      if(nchar(codons[codon_id]) < 3){ #If codon is less than 2 chars, ignore it
        codons[codon_id] = ""
      }else if (grepl("-", codons[codon_id], fixed = TRUE)){ #If it contains missing nucleotides ignore it
        codons[codon_id] = "-" #Maybe i should remove it altogether and not leave a dash
      } else{
        codons[codon_id] = genetic_code[[codons[codon_id]]]
      }
    }
    sequence <- paste(codons, collapse="")
    return(sequence)
  }

  #4
  VDJ_merge_chain <- function(VDJ, samples, path_toData) {
    #Merges cell entries per chain, into one entry per cell.
    #Author: Aurora
    #Binding all contig annotations files in a unique file
    all_contig_annotations = data.frame()
    for (x in samples) {
      x = paste(path_toData,"VDJ/",x,"/",sep="")
      contig_annotations <- read_csv(paste(x,"filtered_contig_annotations.csv",sep="")) %>%
        #Removing final -1 from barcodes
        mutate(VDJ_barcode = substr(barcode,1,nchar(barcode)-2))
      all_contig_annotations = rbind(all_contig_annotations, contig_annotations)#CHANGE BECAUSE OF DIFFERENT collumns LCMV TNFR2 files
    }
    #Dividing into heavy and light chain
    contigs_HC <- subset(all_contig_annotations, chain == "IGH") %>%
      #Selecting columns of interest
      select("fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
             "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
             "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
             "fwr4", "fwr4_nt", "umis", "VDJ_barcode")
    #Adding prefix to differentiate
    colnames(contigs_HC)[0:15] <- paste('HC', colnames(contigs_HC)[0:15], sep = '_')

    #Repeating for light chain
    contigs_LC <- subset(all_contig_annotations, chain %in% c("IGK","IGL")) %>%
      select("fwr1", "fwr1_nt", "cdr1", "cdr1_nt",
             "fwr2", "fwr2_nt", "cdr2", "cdr2_nt",
             "fwr3", "fwr3_nt", "cdr3", "cdr3_nt",
             "fwr4", "fwr4_nt", "umis", "VDJ_barcode")
    colnames(contigs_LC)[0:15] <- paste('LC', colnames(contigs_LC)[0:15], sep = '_')

    #Removing sample nr from barcode from VDJ
    VDJ %<>% mutate(VDJ_barcode = sub(".*_","",barcode))

    #Joining columns of interest to the initial VDJ
    VDJ_HC_contigs <- left_join(VDJ, contigs_HC, by="VDJ_barcode")
    VDJ_contigs <- left_join(VDJ_HC_contigs, contigs_LC, by="VDJ_barcode") %>%
      select(-VDJ_barcode)

    VDJ_contigs$full_VDJ <- paste0(VDJ_contigs$HC_fwr1_nt,
                                   VDJ_contigs$HC_cdr1_nt,
                                   VDJ_contigs$HC_fwr2_nt,
                                   VDJ_contigs$HC_cdr2_nt,
                                   VDJ_contigs$HC_fwr3_nt,
                                   VDJ_contigs$HC_cdr3_nt,
                                   VDJ_contigs$HC_fwr4_nt)

    VDJ_contigs$full_VJ <- paste0(VDJ_contigs$LC_fwr1_nt,
                                  VDJ_contigs$LC_cdr1_nt,
                                  VDJ_contigs$LC_fwr2_nt,
                                  VDJ_contigs$LC_cdr2_nt,
                                  VDJ_contigs$LC_fwr3_nt,
                                  VDJ_contigs$LC_cdr3_nt,
                                  VDJ_contigs$LC_fwr4_nt)

    return(VDJ_contigs)
  }


  #Master function code


  if(ref == FALSE){
    #Just keep the 1 light 1 heavy chain and find the FULL VDJ
    VDJ <- VDJ[grepl(";",VDJ$VDJ_chain_contig) == FALSE,]
    VDJ <- VDJ[grepl(";",VDJ$VJ_chain_contig) == FALSE,]
    VDJ = VDJ_merge_chain(VDJ, samples, path_toData)
    return(VDJ)
  }
  if("VDJ_ref.nt" %in% names(VDJ) | "VDJ_ref.aa" %in% names(VDJ) | "VJ_ref.nt" %in% names(VDJ) | "VJ_ref.aa" %in% names(VDJ)){
    message("Reference is already trimmed")
    return(VDJ)
  }
  sample_id = 0

  VDJ["VDJ_ref.nt"] = "None"
  VDJ["VJ_ref.nt"] = "None"
  VDJ["VDJ_ref.aa"] = "None"
  VDJ["VJ_ref.aa"] = "None"
  VDJ["rank_post_filter"] = -1
  clone_counts_all = c()
  for (i in samples){
    #Path to the sample data
    path = paste(path_toData, "VDJ/",i,"/",sep = "")
    #Examined sample_id
    sample_id = sample_id + 1
    #read the consensus
    consensus = read.csv(paste(path,"consensus_annotations.csv",sep=""), sep=",")
    vgm_subset = VDJ[VDJ$sample_id == paste("s",sample_id,sep = ""),]
    #Filtering for one light one heavy chain
    VDJ <- VDJ[grepl(";",VDJ$VDJ_chain_contig) == FALSE,]
    VDJ <- VDJ[grepl(";",VDJ$VJ_chain_contig) == FALSE,]

    clone_counts = table(vgm_subset$clonotype_id)

    total_clones = length(clone_counts)
    used_clones = n_clones
    if(is.na(n_clones)){
      message("Computing ref for all ", total_clones," clones")
      used_clones = total_clones
    }else if(n_clones>total_clones){
      message("Wanted clones exceed the clones in the VDJ. Computing ref for the top ", total_clones," clones instead")
      used_clones = total_clones
    }
    topn_clones = names(clone_counts[order(-clone_counts)][1:used_clones])
    clone_counts_all[[i]] = topn_clones
    print(clone_counts_all)
    #for each clonotype find the reference
    clone_rank = 1
    for (clone in topn_clones){
      VDJ[(VDJ$sample_id == paste("s",sample_id,sep = "") & VDJ$clonotype_id == clone), "rank_post_filter"] = clone_rank
      #Fetching the light and heavy reference IS IT REALLY 2-LIGHT AND 1-HEAVY?
      vgm_clone = vgm_subset[(vgm_subset$clonotype_id == clone),]#CHANGE NOT TO BE DONE ONCE PER ITEM
      consensus_curr = consensus[consensus$clonotype_id == clone,]
      #In case clone does not exist
      if (is.null(vgm_clone)){
        next
      }
      #Trim consensus and map to reference
      examined_refs = c()
      for (row_n in rownames(vgm_clone)){
        vgm_cell = vgm_clone[row_n,]
        out = ref_per_cell(vgm_cell,consensus_curr,examined_refs)
        examined_refs = out[[2]]
        output = out[[1]]
        for(name in names(output)){
          VDJ[VDJ$barcode == vgm_cell$barcode,name] = output[name]
        }
      }
      clone_rank = clone_rank + 1
    }
  }
  VDJ = VDJ_merge_chain(VDJ, samples, path_toData)
  return(list(vdj = VDJ,clones = clone_counts_all))
}
