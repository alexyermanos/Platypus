#' MiXCR wrapper for Platypus V3 VDJ object
#'
#'@description Extracts information on the VDJRegion level using MiXCR on WINDOWS, MAC and UNIX systems for input from both Platypus v2 (VDJ.per.clone) or v3 (Output of VDJ_GEX_matrix) This function assumes the user can run an executable instance of MiXCR and is elgible to use MiXCR as determined by license agreements. ! FOR WINDOWS USERS THE EXECUTABLE MIXCR.JAR HAS TO PRESENT IN THE CURRENT WORKING DIRECTORY ! In case the function is not able to call mixcr in "Windows" mode, please try setting the operating.system to "Linux".  The VDJRegion corresponds to the recombined heavy and light chain loci starting from framework region 1 (FR1) and extending to frame work region 4 (FR4). This can be useful for extracting full-length sequences ready to clone and further calculating somatic hypermutation occurrences.
#' @param VDJ For platypus.version = "v2" the output from the VDJ_per_clone function. This object should have information regarding the contigs and clonotype_ids for each cell. For platypus.version = "v3" the VDJ dataframe output of the VDJ_GEX_matrix function (VDJ.GEX.matri.output[[1]])
#' @param operating.system Can be either "Windows", "Darwin" (for MAC) or "Linux". If left empty this is detected automatically
#' @param custom.cmd.call This function calls Mixcr via a terminal. The commands are build as follows: FOR DARWIN and LINUX machines: paste0(custom.cmd.call," ", mixcr.directory,"/mixcr.jar ... {mixcr parameters assembled by function}"). FOR WINDOWS machines: paste0(custom.cmd.call," mixcr.jar ... {mixcr parameters assembled by function}") custom.cmd.call defaults to "java -jar". In case different MIXCR installations, changing this to meet your systems requirements may be necessary.
#' @param mixcr.directory The directory containing an executable version of MiXCR. FOR WINDOWS USERS THIS IS SET TO THE CURRENT WORKING DIRECTORY (please paste the content of the MIXCR folder after unzipping into your working directory. Make sure, that mixcr.jar is not within any subfolders.)
#' @param species Either "mmu" for mouse or "hsa" for human. These use the default germline genes for both species contained in MIXCR. Defaults to "hsa"
#' @param simplify Only relevant when platypus.version = "v3". Boolean. Defaults to TRUE. If FALSE the full MIXCR output and computed SHM column is appended to the VDJ If TRUE only the framework and CDR3 region columns and computed SHM column is appended. To discriminate between VDJ and VJ chains, prefixes are added to all MIXCR output columns
#' @param platypus.version Character. Defaults to "v3". Can be "v2" or "v3" dependent on the input format
#' @return For platypus.version = "v3" returns input VDJ dataframe supplemented with MIXCR output information. For platypus.version = "v2" returns a nested list containing VDJRegion information as determined by MIXCR. The outer list corresponds to the individual repertoires in the same structure as the input  VDJ.per.clone. The inner list corresponds to each clonal family, as determined by either the VDJ_clonotype function or the defaul nucleotide clonotyping produced by cellranger.Each element in the inner list corresponds to a dataframe containing repertoire information such as isotype, CDR sequences, mean number of UMIs. This output can be supplied to further package functions such as VDJ_extract_sequences and VDJ_GEX_integrate.
#' @seealso VDJ_extract_sequences
#' @export
#' @examples
#' \dontrun{
#'#For platypus version 2
#'VDJ_call_MIXCR(VDJ = VDJ.per.clone.output,
#'mixcr.directory = "~/Downloads/mixcr-3.0.12/mixcr",species = "mmu")
#'
#'#For platypus version 3 on a Windows system
#'VDJ_call_MIXCR(VDJ = VDJ_GEX_matrix.output[[1]],
#'mixcr.directory = "WILL BE SET TO CURRENT WORKING DIRECTORY",
#'species = "mmu", platypus.version = "v3", simplify = TRUE)
#'}

VDJ_call_MIXCR <- function(VDJ,
                           operating.system,
                           custom.cmd.call,
                           mixcr.directory,
                           species,
                           simplify,
                           platypus.version){
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL
  descrsR1 <- NULL

    if(missing(simplify)) simplify <- T
    if(missing(species)) species <- "hsa"
    if(missing(platypus.version)) platypus.version <- "v3"
    if(missing(custom.cmd.call)) custom.cmd.call <- "java -jar"

    if(missing(operating.system)){
      switch(Sys.info()[['sysname']],
             Windows= {message("Windows system detected")
                       operating.system <- "Windows"},
             Linux  = {message("Linux system detected")
                      operating.system <- "Linux"},
             Darwin = {message("MAC system detected")
                      operating.system <- "Darwin"})

    }

    switch(platypus.version,
           v3 = {if(!inherits(VDJ,"data.frame")){stop("When selecting platypus.version = 'v3', please input the VDJ matrix of the output of the VDJ_GEX_matrix function (usually VDJ.GEX.matrix.output[[1]]")}},
           v2 = {if(!inherits(VDJ,"list")){stop("When selecting platypus.version = 'v2', please input the list output of VDJ.per.clone")}},
           {stop("Please input either 'v2' or 'v3' as platypus.version")})

    #Everything set up, now going through all four possibilities
    if(platypus.version == "v3" & operating.system == "Windows"){

    mixcr.directory <- getwd()

    system("cmd.exe", input = paste("cd ", getwd(),sep=""))
    system("cmd.exe", input = paste("del temphc.fasta"))
    system("cmd.exe", input = paste("del templc.fasta"))
    system("cmd.exe", input = paste("del tempmixcrhc.out.vdjca"))
    system("cmd.exe", input = paste("del tempmixcrlc.out.vdjca"))
    system("cmd.exe", input = paste("del tempmixcrhc.out.txt"))
    system("cmd.exe", input = paste("del tempmixcrlc.out.txt"))

    #Deal with aberrant cells first
    aberrants <- subset(VDJ, Nr_of_VDJ_chains == 2 | Nr_of_VJ_chains == 2)
    bcs_VDJ <- c()
    bcs_VJ <- c()
    VDJ_raw <- c()
    VJ_raw <- c()
    if(nrow(aberrants)>0){
    for(i in 1:nrow(aberrants)){
        bcs_VDJ <- c(bcs_VDJ, rep(aberrants$barcode[i],aberrants$Nr_of_VDJ_chains[i]))
        VDJ_raw <- c(VDJ_raw, unlist(stringr::str_split(aberrants$VDJ_sequence_nt_raw[i], ";", simplify = T)[1,]))
        bcs_VJ <- c(bcs_VJ, rep(aberrants$barcode[i],aberrants$Nr_of_VJ_chains[i]))
        VJ_raw <- c(VJ_raw, unlist(stringr::str_split(aberrants$VJ_sequence_nt_raw[i], ";", simplify = T)[1,]))
    }
    }

    normals <- subset(VDJ, Nr_of_VDJ_chains < 2 & Nr_of_VJ_chains < 2)

  ### need to also read in the fastas
    temp.seq.hc_unlist <- c(unlist(normals$VDJ_sequence_nt_raw),VDJ_raw)
    temp.seq.lc_unlist <- c(unlist(normals$VJ_sequence_nt_raw),VJ_raw)
    temp.name.hc_unlist <- c(unlist(normals$barcode),bcs_VDJ)
    temp.name.lc_unlist <- c(unlist(normals$barcode),bcs_VJ)

    seqinr::write.fasta(sequences = as.list(temp.seq.hc_unlist),names = temp.name.hc_unlist,file.out = "temphc.fasta")

    system("cmd.exe", input = paste0("echo Hello world"))

    print(paste0(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"')))

    system("cmd.exe", input = paste0(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"')))

           system("cmd.exe", input = paste0("java -jar mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"')))

    system("cmd.exe", input = paste0(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"')))

    system("cmd.exe",input =  paste0(custom.cmd.call, " ", "mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrhc.out.vdjca ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.txt",'"')))

    seqinr::write.fasta(sequences = as.list(temp.seq.lc_unlist),names = temp.name.lc_unlist,file.out = "templc.fasta")

    system("cmd.exe", input = paste0(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," templc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrlc.out.vdjca",'"')))

    system("cmd.exe", input = paste0(custom.cmd.call, " ", "mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrlc.out.vdjca ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrlc.out.txt",'"')))

    temp.mixcr.hc <- utils::read.table(file ="tempmixcrhc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
    temp.mixcr.lc <- utils::read.table(file ="tempmixcrlc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)

    system("cmd.exe", input = paste("del temphc.fasta"))
    system("cmd.exe", input = paste("del templc.fasta"))
    system("cmd.exe", input = paste("del tempmixcrhc.out.vdjca"))
    system("cmd.exe", input = paste("del tempmixcrlc.out.vdjca"))
    system("cmd.exe", input = paste("del tempmixcrhc.out.txt"))
    system("cmd.exe", input = paste("del tempmixcrlc.out.txt"))

    ## now need to fill VDJ
    if(simplify == F){

      aberrants_mixcr <- subset(temp.mixcr.hc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

      to_merge_hc <- rbind(subset(temp.mixcr.hc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)
      #add SHM measures
      to_merge_hc$SHM <- stringr::str_count(to_merge_hc$bestVAlignment,"S") + stringr::str_count(to_merge_hc$bestJAlignment,"S") + stringr::str_count(to_merge_hc$bestDAlignment,"S")
      #append prefix to names
      names(to_merge_hc) <- paste0("VDJ_",names(to_merge_hc))
      #add barcode for merge
      to_merge_hc$barcode <- to_merge_hc$VDJ_descrsR1


      aberrants_mixcr <- subset(temp.mixcr.lc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

      to_merge_lc <- rbind(subset(temp.mixcr.lc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)
      #add SHM measures
      to_merge_lc$SHM <- stringr::str_count(to_merge_lc$bestVAlignment,"S") + stringr::str_count(to_merge_lc$bestJAlignment,"S")
      #append prefix to names
      names(to_merge_lc) <- paste0("VJ_",names(to_merge_lc))
      #add barcode for merge
      to_merge_lc$barcode <- to_merge_lc$VJ_descrsR1

    } else {
      aberrants_mixcr <- subset(temp.mixcr.hc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

      to_merge_hc <- rbind(subset(temp.mixcr.hc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)[,c("nSeqFR1","nSeqCDR1", "nSeqFR2","nSeqCDR2", "nSeqFR3","nSeqCDR3", "nSeqFR4", "aaSeqFR1","aaSeqCDR1", "aaSeqFR2","aaSeqCDR2", "aaSeqFR3","aaSeqCDR3", "aaSeqFR4","bestVAlignment","bestJAlignment","bestDAlignment","descrsR1")]
      #add SHM measures
      to_merge_hc$SHM <- stringr::str_count(to_merge_hc$bestVAlignment,"S") + stringr::str_count(to_merge_hc$bestJAlignment,"S") + stringr::str_count(to_merge_hc$bestDAlignment,"S")
      #append prefix to names
      names(to_merge_hc) <- paste0("VDJ_",names(to_merge_hc))
      #add barcode for merge
      to_merge_hc$barcode <- to_merge_hc$VDJ_descrsR1

      aberrants_mixcr <- subset(temp.mixcr.lc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
      aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

      to_merge_lc <- rbind(subset(temp.mixcr.lc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)[,c("nSeqFR1","nSeqCDR1", "nSeqFR2","nSeqCDR2", "nSeqFR3","nSeqCDR3", "nSeqFR4", "aaSeqFR1","aaSeqCDR1", "aaSeqFR2","aaSeqCDR2", "aaSeqFR3","aaSeqCDR3", "aaSeqFR4","bestVAlignment","bestJAlignment","descrsR1")]
      #add SHM measures
      to_merge_lc$SHM <- stringr::str_count(to_merge_lc$bestVAlignment,"S") + stringr::str_count(to_merge_lc$bestJAlignment,"S")

      #append prefix to names
      names(to_merge_lc) <- paste0("VJ_",names(to_merge_lc))
      #add barcode for merge
      to_merge_lc$barcode <- to_merge_lc$VJ_descrsR1
    }

    ncol_raw <- ncol(VDJ)
    VDJ <- merge(VDJ, to_merge_hc, by = "barcode", all.x = T)
    VDJ <- merge(VDJ, to_merge_lc, by = "barcode", all.x = T)

    #cleanup
    for(i in (ncol_raw+1):ncol(VDJ)){
      if(names(VDJ)[i] %in% c("VDJ_SHM", "VJ_SHM")){
        #leave these numeric columns as they are
      } else {
        #replace NA from merging with empty strings similarly to other parts of the VDJ
        VDJ[is.na(VDJ[,i]),i] <- ""
      }
    }
    return(VDJ)


    } else if (platypus.version == "v2" & operating.system == "Windows") {

      mixcr.directory <- getwd()

      system("cmd.exe", input = paste("cd ", getwd(),sep=""))
      system("cmd.exe", input = paste("del temphc.fasta"))
      system("cmd.exe", input = paste("del templc.fasta"))
      system("cmd.exe", input = paste("del tempmixcrhc.out.vdjca"))
      system("cmd.exe", input = paste("del tempmixcrlc.out.vdjca"))
      system("cmd.exe", input = paste("del tempmixcrhc.out.txt"))
      system("cmd.exe", input = paste("del tempmixcrlc.out.txt"))


      ### need to also read in the fastas
      for(i in 1:length(VDJ.per.clone)){
        temp.seq.hc <- list()
        temp.seq.lc <- list()
        temp.name.hc <- list()
        temp.name.lc <- list()
        for(j in 1:length(VDJ.per.clone[[i]])){
          temp.seq.hc[[j]] <- VDJ.per.clone[[i]][[j]]$full_HC_sequence
          temp.seq.lc[[j]] <- VDJ.per.clone[[i]][[j]]$full_LC_sequence
          temp.name.hc[[j]] <- VDJ.per.clone[[i]][[j]]$contig_id_hc
          temp.name.lc[[j]] <- VDJ.per.clone[[i]][[j]]$contig_id_lc
        }
        temp.seq.hc_unlist <- unlist(temp.seq.hc)
        temp.seq.lc_unlist <- unlist(temp.seq.lc)
        temp.name.hc_unlist <- unlist(temp.name.hc)
        temp.name.lc_unlist <- unlist(temp.name.lc)

        seqinr::write.fasta(sequences = as.list(temp.seq.hc_unlist),names = temp.name.hc_unlist,file.out = "temphc.fasta")

        system("cmd.exe", input = paste(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"'),sep=""))


        system("cmd.exe", input = paste(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.vdjca",'"'),sep=""))


        system("cmd.exe",input =  paste(custom.cmd.call, " ", "mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrhc.out.vdjca ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrhc.out.txt",'"'),sep=""))

        seqinr::write.fasta(sequences = as.list(temp.seq.lc_unlist),names = temp.name.lc_unlist,file.out = "templc.fasta")

        system("cmd.exe", input = paste(custom.cmd.call, " ", "mixcr.jar align -OsaveOriginalReads=true -s ", species," templc.fasta ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrlc.out.vdjca",'"'),sep=""))

        system("cmd.exe", input = paste(custom.cmd.call, " ", "mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrlc.out.vdjca ",paste0('"',gsub("\\\\","/",mixcr.directory),"/tempmixcrlc.out.txt",'"'),sep=""))

        temp.mixcr.hc <- utils::read.table(file ="tempmixcrhc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
        temp.mixcr.lc <- utils::read.table(file ="tempmixcrlc.out.txt", sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)

        system("cmd.exe", input = paste("del temphc.fasta"))
        system("cmd.exe", input = paste("del templc.fasta"))
        system("cmd.exe", input = paste("del tempmixcrhc.out.vdjca"))
        system("cmd.exe", input = paste("del tempmixcrlc.out.vdjca"))
        system("cmd.exe", input = paste("del tempmixcrhc.out.txt"))
        system("cmd.exe", input = paste("del tempmixcrlc.out.txt"))


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




      } else if (platypus.version == "v3" & operating.system %in% c("Darwin", "Linux")) {

      #Deal with aberrant cells first
      aberrants <- subset(VDJ, Nr_of_VDJ_chains == 2 | Nr_of_VJ_chains == 2)
      bcs_VDJ <- c()
      bcs_VJ <- c()
      VDJ_raw <- c()
      VJ_raw <- c()
      if(nrow(aberrants)>0){
        for(i in 1:nrow(aberrants)){
          bcs_VDJ <- c(bcs_VDJ, rep(aberrants$barcode[i],aberrants$Nr_of_VDJ_chains[i]))
          VDJ_raw <- c(VDJ_raw, unlist(stringr::str_split(aberrants$VDJ_sequence_nt_raw[i], ";", simplify = T)[1,]))
          bcs_VJ <- c(bcs_VJ, rep(aberrants$barcode[i],aberrants$Nr_of_VJ_chains[i]))
          VJ_raw <- c(VJ_raw, unlist(stringr::str_split(aberrants$VJ_sequence_nt_raw[i], ";", simplify = T)[1,]))
        }
      }

      normals <- subset(VDJ, Nr_of_VDJ_chains < 2 & Nr_of_VJ_chains < 2)

      ### need to also read in the fastas
      temp.seq.hc_unlist <- c(unlist(normals$VDJ_sequence_nt_raw),VDJ_raw)
      temp.seq.lc_unlist <- c(unlist(normals$VJ_sequence_nt_raw),VJ_raw)
      temp.name.hc_unlist <- c(unlist(normals$barcode),bcs_VDJ)
      temp.name.lc_unlist <- c(unlist(normals$barcode),bcs_VJ)

      seqinr::write.fasta(sequences = as.list(temp.seq.hc_unlist),names = temp.name.hc_unlist,file.out = paste0("temphc.fasta"))
      system(paste0(custom.cmd.call, " ", mixcr.directory,"/mixcr.jar align -OsaveOriginalReads=true -s ", species," temphc.fasta tempmixcrhc.out.vdjca",sep=""))
      system(paste0(custom.cmd.call, " ", mixcr.directory,"/mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrhc.out.vdjca tempmixcrhc.out.txt",sep=""))
      seqinr::write.fasta(sequences = as.list(temp.seq.lc_unlist),names = temp.name.lc_unlist,file.out = paste0("templc.fasta"))
      system(paste0(custom.cmd.call, " ", mixcr.directory,"/mixcr.jar align -OsaveOriginalReads=true -s ", species," templc.fasta tempmixcrlc.out.vdjca",sep=""))
      system(paste0(custom.cmd.call, " ", mixcr.directory,"/mixcr.jar exportAlignments --preset full -descrsR1 -vAlignment -dAlignment -jAlignment -aaMutations VRegion -aaMutations JRegion -nMutations VRegion -nMutations JRegion tempmixcrlc.out.vdjca tempmixcrlc.out.txt",sep=""))
      temp.mixcr.hc <- utils::read.table(file =paste0("tempmixcrhc.out.txt"), sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)
      temp.mixcr.lc <- utils::read.table(file =paste0("tempmixcrlc.out.txt"), sep="\t", header = T, stringsAsFactors=FALSE, fill=TRUE)

      system("rm temphc.fasta")
      system("rm templc.fasta")
      system("rm tempmixcrhc.out.vdjca")
      system("rm tempmixcrlc.out.vdjca")
      system("rm tempmixcrhc.out.txt")
      system("rm tempmixcrlc.out.txt")

      ## now need to fill VDJ
      ## now need to fill VDJ
      if(simplify == F){

        aberrants_mixcr <- subset(temp.mixcr.hc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

        to_merge_hc <- rbind(subset(temp.mixcr.hc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)
        #add SHM measures
        to_merge_hc$SHM <- stringr::str_count(to_merge_hc$bestVAlignment,"S") + stringr::str_count(to_merge_hc$bestJAlignment,"S") + stringr::str_count(to_merge_hc$bestDAlignment,"S")
        #append prefix to names
        names(to_merge_hc) <- paste0("VDJ_",names(to_merge_hc))
        #add barcode for merge
        to_merge_hc$barcode <- to_merge_hc$VDJ_descrsR1


        aberrants_mixcr <- subset(temp.mixcr.lc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

        to_merge_lc <- rbind(subset(temp.mixcr.lc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)
        #add SHM measures
        to_merge_lc$SHM <- stringr::str_count(to_merge_lc$bestVAlignment,"S") + stringr::str_count(to_merge_lc$bestJAlignment,"S")
        #append prefix to names
        names(to_merge_lc) <- paste0("VJ_",names(to_merge_lc))
        #add barcode for merge
        to_merge_lc$barcode <- to_merge_lc$VJ_descrsR1

      } else {
        aberrants_mixcr <- subset(temp.mixcr.hc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

        to_merge_hc <- rbind(subset(temp.mixcr.hc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)[,c("nSeqFR1","nSeqCDR1", "nSeqFR2","nSeqCDR2", "nSeqFR3","nSeqCDR3", "nSeqFR4", "aaSeqFR1","aaSeqCDR1", "aaSeqFR2","aaSeqCDR2", "aaSeqFR3","aaSeqCDR3", "aaSeqFR4","bestVAlignment","bestJAlignment","bestDAlignment","descrsR1")]
        #add SHM measures
        to_merge_hc$SHM <- stringr::str_count(to_merge_hc$bestVAlignment,"S") + stringr::str_count(to_merge_hc$bestJAlignment,"S") + stringr::str_count(to_merge_hc$bestDAlignment,"S")
        #append prefix to names
        names(to_merge_hc) <- paste0("VDJ_",names(to_merge_hc))
        #add barcode for merge
        to_merge_hc$barcode <- to_merge_hc$VDJ_descrsR1

        aberrants_mixcr <- subset(temp.mixcr.lc, descrsR1 %in% descrsR1[duplicated(descrsR1)])
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::group_by(descrsR1) %>% dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~paste0(.x, collapse = ";")))
        aberrants_mixcr <- aberrants_mixcr %>% dplyr::mutate(dplyr::across(.cols = dplyr::everything(), .fns = ~gsub("(^;)|(;$)","",.x)))

        to_merge_lc <- rbind(subset(temp.mixcr.lc, !descrsR1 %in% descrsR1[duplicated(descrsR1)]), aberrants_mixcr)[,c("nSeqFR1","nSeqCDR1", "nSeqFR2","nSeqCDR2", "nSeqFR3","nSeqCDR3", "nSeqFR4", "aaSeqFR1","aaSeqCDR1", "aaSeqFR2","aaSeqCDR2", "aaSeqFR3","aaSeqCDR3", "aaSeqFR4","bestVAlignment","bestJAlignment","descrsR1")]
        #add SHM measures
        to_merge_lc$SHM <- stringr::str_count(to_merge_lc$bestVAlignment,"S") + stringr::str_count(to_merge_lc$bestJAlignment,"S")

        #append prefix to names
        names(to_merge_lc) <- paste0("VJ_",names(to_merge_lc))
        #add barcode for merge
        to_merge_lc$barcode <- to_merge_lc$VJ_descrsR1
      }

      ncol_raw <- ncol(VDJ)
      VDJ <- merge(VDJ, to_merge_hc, by = "barcode", all.x = T)
      VDJ <- merge(VDJ, to_merge_lc, by = "barcode", all.x = T)

      #cleanup
      for(i in (ncol_raw+1):ncol(VDJ)){
        if(names(VDJ)[i] %in% c("VDJ_SHM", "VJ_SHM")){
          #leave these numeric columns as they are
        } else {
          #replace NA from merging with empty strings similarly to other parts of the VDJ
          VDJ[is.na(VDJ[,i]),i] <- ""
        }
      }
      return(VDJ)

    } else if (platypus.version == "v2" & operating.system %in% c("Darwin", "Linux")) {


      if(missing(mixcr.directory)) message("No mixcr directory supplied. Assuming working directory holds mixcr executable")
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
}


