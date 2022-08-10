#' MiXCR wrapper for Platypus V3 VDJ object. In addition to the VDJ_call_MIXCR function, the output also contains the concatenated sequences
#' from FR1 all the way to FR2 for both the VDJ and VJ. 
#'
#'@description Extracts information on the VDJRegion level using MiXCR on WINDOWS, MAC and UNIX systems for input from both Platypus v2 (VDJ.per.clone) or v3 (Output of VDJ_GEX_matrix) 
#'This function assumes the user can run an executable instance of MiXCR and is elgible to use MiXCR as determined by license agreements. 
#'! FOR WINDOWS USERS THE EXECUTABLE MIXCR.JAR HAS TO PRESENT IN THE CURRENT WORKING DIRECTORY ! 
#'The VDJRegion corresponds to the recombined heavy and light chain loci starting from framework region 1 (FR1) and extending to frame work region 4 (FR4). 
#'This can be useful for extracting full-length sequences ready to clone and further calculating somatic hypermutation occurrences.
#'In addition to the VDJ_call_MIXCR function, the output also contains the concatenated sequences from FR1 all the way to FR2 for both the VDJ and VJ. 
#' @param VDJ For platypus.version = "v2" the output from the VDJ_per_clone function. This object should have information regarding the contigs and clonotype_ids for each cell. For platypus.version = "v3" the VDJ dataframe output of the VDJ_GEX_matrix function (VDJ.GEX.matri.output[[1]])
#' @param operating.system Can be either "Windows", "Darwin" (for MAC) or "Linux". If left empty this is detected automatically
#' @param mixcr.directory The directory containing an executable version of MiXCR. FOR WINDOWS USERS THIS IS SET TO THE CURRENT WORKING DIRECTORY (please paste the content of the MIXCR folder after unzipping into your working directory. Make sure, that mixcr.jar is not within any subfolders.)
#' @param species Either "mmu" for mouse or "hsa" for human. These use the default germline genes for both species contained in MIXCR. Defaults to "hsa"
#' @param simplify Only relevant when platypus.version = "v3". Boolean. Defaults to TRUE. If FALSE the full MIXCR output and computed SHM column is appended to the VDJ If TRUE only the framework and CDR3 region columns and computed SHM column is appended. To discriminate between VDJ and VJ chains, prefixes are added to all MIXCR output columns
#' @param platypus.version Character. Defaults to "v3". Can be "v2" or "v3" dependent on the input format


VDJ_call_MIXCR_full <- function (VDJ , 
                                 mixcr.directory,
                                 species, 
                                 platypus.version, 
                                 operating.system, 
                                 simplify){
  
  if(missing(simplify)) simplify <- T
  if(missing(species)) species <- "hsa"
  if(missing(platypus.version)) platypus.version <- "v3"
  
  if(missing(operating.system)){
    switch(Sys.info()[['sysname']],
           Windows= {message("Windows system detected")
             operating.system <- "Windows"},
           Linux  = {message("Linux system detected")
             operating.system <- "Linux"},
           Darwin = {message("MAC system detected")
             operating.system <- "Darwin"})
    
  }
  
  VDJ_mixcr_out <- Platypus::VDJ_call_MIXCR(VDJ = VDJ , 
                                  mixcr.directory = mixcr.directory, 
                                  species = species, 
                                  platypus.version = platypus.version, 
                                  operating.system = operating.system, 
                                  simplify = simplify) 
  
  VDJ_mixcr_out$VDJ_nt_mixcr <- paste0(VDJ_mixcr_out$VDJ_nSeqFR1,VDJ_mixcr_out$VDJ_nSeqCDR1, VDJ_mixcr_out$VDJ_nSeqFR2,VDJ_mixcr_out$VDJ_nSeqCDR2, VDJ_mixcr_out$VDJ_nSeqFR3,VDJ_mixcr_out$VDJ_nSeqCDR3, VDJ_mixcr_out$VDJ_nSeqFR4 )
  
  VDJ_mixcr_out$VJ_nt_mixcr <- paste0(VDJ_mixcr_out$VJ_nSeqFR1,VDJ_mixcr_out$VJ_nSeqCDR1, VDJ_mixcr_out$VJ_nSeqFR2,VDJ_mixcr_out$VJ_nSeqCDR2, VDJ_mixcr_out$VJ_nSeqFR3,VDJ_mixcr_out$VJ_nSeqCDR3, VDJ_mixcr_out$VJ_nSeqFR4 )
  
  VDJ_mixcr_out$VDJ_aa_mixcr <- paste0(VDJ_mixcr_out$VDJ_aaSeqFR1,VDJ_mixcr_out$VDJ_aaSeqCDR1, VDJ_mixcr_out$VDJ_aaSeqFR2,VDJ_mixcr_out$VDJ_aaSeqCDR2, VDJ_mixcr_out$VDJ_aaSeqFR3,VDJ_mixcr_out$VDJ_aaSeqCDR3, VDJ_mixcr_out$VDJ_aaSeqFR4 )
  
  VDJ_mixcr_out$VJ_aa_mixcr <- paste0(VDJ_mixcr_out$VJ_aaSeqFR1,VDJ_mixcr_out$VJ_aaSeqCDR1, VDJ_mixcr_out$VJ_aaSeqFR2,VDJ_mixcr_out$VJ_aaSeqCDR2, VDJ_mixcr_out$VJ_aaSeqFR3,VDJ_mixcr_out$VJ_aaSeqCDR3, VDJ_mixcr_out$VJ_aaSeqFR4 )
  
  
  
  return(VDJ_mixcr_out)
  
}