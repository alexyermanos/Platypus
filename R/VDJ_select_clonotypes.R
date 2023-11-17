#'Select clonotypes
#'
#'@description For prediction of antibody structures from a big data set it might be of interest to select the most expanded clonotypes for prediction.
#'This function can select the top most expanded clonotypes based on the desired clone strategy.
#'Among the most expanded clonotypes the cells are ranked according to the UMI count and then the top unique sequences are selected to use for prediction.
#'The function's input is the Platypus VGM object. In order to integrate UMI counts to the data, the raw data which is the output of the PlatypusDB_fetch() function is needed in addition.
#'From the selected clonotypes the germline reference sequences are obtained by calling MIXCR. This requires a local installation of MIXCR on your computer.
#'!FOR WINDOWS USERS THE EXECUTABLE MIXCR.JAR HAS TO PRESENT IN THE CURRENT WORKING DIRECTORY !
#'
#'The output of the VDJ_select_clonotypes function can directly be used for structure prediction by the AlphaFold_prediction() function.
#'
#'
#'@param VGM The platypus vgm object his used as an input for the function.
#'@param raw.data In order to integrate the UMI counts per cell, the raw data has to be specified as a second input to the function which is the output of the PlatypusDB_fetch() function.
#'@param clone.strategy The desired clone strategy can be specified as a string. Possible options are 10x.default, cdr3.nt, cdr3.aa, VDJJ.VJJ, VDJJ.VJJ.cdr3length, VDJJ.VJJ.cdr3length.cdr3.homology, VDJJ.VJJ.cdr3length.VDJcdr3.homology, cdr3.homology, VDJcdr3.homology. 10x.default is used as default. cdr3.aa will convert the default cell ranger clonotyping to amino acid bases. 'VDJJ.VJJ' groups B cells with identical germline genes (V and J segments for both heavy chain and light chain. Those arguments including 'cdr3length' will group all sequences with identical VDJ and VJ CDR3 sequence lengths. Those arguments including 'cdr3.homology' will additionally impose a homology requirement for CDRH3 and CDRL3 sequences.'CDR3.homology',or 'CDRH3.homology' will group sequences based on homology only (either of the whole CDR3 sequence or of the VDJ CDR3 sequence respectively). All homology calculations are performed on the amino acid level.
#'@param VDJ.VJ.1chain If the VDJ.VJ.1chain argument is set to TRUE only cells with  one VDJ and one VJ sequences are included in the selection.
#'@param donut.plot If set to TRUE a donut plot for visualization of the clonotypes is returned.
#'@param clonotypes.per.sample By default the top clonotypes are selected per sample. If the top clonotypes over all samples are desired the clonotypes.per.sample argument can be set to FALSE.
#'@param top.clonotypes Specify the number of top clonotypes that will be selected either per sample if clonotypes.per.sample = T or over all if clonotypes.per.sample = F.
#'@param seq.per.clonotype Specify the number of unique sequences per clonotype that are selected. The clonotypes are ordered according to UMI expression.
#'@param mixcr.directory The path to the directory containing an executable version of MIXCR.
#'@param species Either "mmu" for mouse or "hsa" for human. These use the default germline genes for both species contained in MIXCR. Default is set to "hsa".
#'@param platypus.version Character. Defaults to "v3". Can be "v2" or "v3" dependent on the input format
#'@param operating.system Can be either "Windows", "Darwin" (for MAC) or "Linux". If left empty this is detected automatically
#'@param simplify Only relevant when platypus.version = "v3". Boolean. Defaults to TRUE. If FALSE the full MIXCR output and computed SHM column is appended to the VDJ If TRUE only the framework and CDR3 region columns and computed SHM column is appended. To discriminate between VDJ and VJ chains, prefixes are added to all MIXCR output columns
#' @return ADD DESCRIPTION OF RETURN VALUE HERE
#' @export
#' @examples
#' \dontrun{
#'
#' ADD EXAMPLE CODE HERE
#'
#' }


VDJ_select_clonotypes <- function(VGM,
                                  raw.data,
                                  clone.strategy,
                                  VDJ.VJ.1chain,
                                  donut.plot,
                                  clonotypes.per.sample,
                                  top.clonotypes,
                                  seq.per.clonotype,
                                  mixcr.directory,
                                  species,
                                  platypus.version,
                                  operating.system,
                                  simplify
                                  ) {



  if(missing(VGM)) {stop("Please specify a VGM object to run this function")}
  if(missing(raw.data)) {stop("Please specify the raw input data to integrate the UMI counts to the VDJ dataframe")}
  if(missing(clone.strategy)) {clone.strategy <- "10x.default"}
  if(missing(VDJ.VJ.1chain)) {VDJ.VJ.1chain <- F}
  if(missing(donut.plot)) {donut.plot <- F}
  if(missing(clonotypes.per.sample)) {clonotypes.per.sample <- T}
  if(missing(top.clonotypes)) {top.clonotypes <- 25}
  if(missing(seq.per.clonotype)) {seq.per.clonotype <- 3}
  if(missing(mixcr.directory)) {stop("Please specify a path to the local mixcr installation ('mixcr.jar' file)")}
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


  names.top.clonotypes <- NULL
  VGM.top.clonotypes <- NULL
  unique.clonotypes <- NULL
  VGM.top.clonotypes.mixcr.out <- NULL
  chain <- NULL
  barcode <- NULL
  umis <- NULL
  reads <- NULL
  sample_id_clonotype <- NULL
  sample_id <- NULL
  Nr_of_VDJ_chains <- NULL
  VDJ_aa_mixcr <- NULL
  Nr_of_VJ_chains <- NULL
  VJ_aa_mixcr <- NULL

  #Integrate the Heavy Chain UMI and read count to the data
  df_full <- data.frame()
  ii <- 0
  for(sam in names(raw.data)){
    ii <- ii + 1
    df <- raw.data[[sam]]$VDJ$filtered_contig_annotations.csv
    df <- dplyr::filter(df,chain == "IGH")
    df <- dplyr::distinct(df,barcode, .keep_all = TRUE)
    df <- dplyr::select(df,barcode,umis,reads)
    df$orig_barcode <- df$barcode
    df$barcode <- paste0(rep(paste0("s",ii,"_"), nrow(df)),df$orig_barcode)
    df$orig_barcode <- NULL
    df_full <- dplyr::bind_rows(df_full,df)

  }
  VGM[[1]] <- dplyr::left_join(VGM[[1]],df_full, by = "barcode")




  ## Determine clonotypes based on a selected clonotype
  VGM[[1]] <- Platypus::VDJ_clonotype(VDJ = VGM[[1]],
                        VDJ.VJ.1chain = VDJ.VJ.1chain,
                        clone.strategy = clone.strategy)

  ## Show Donut Plot
  if(donut.plot) {
  plots <- Platypus::VDJ_clonal_donut(VDJ = VGM[[1]],
                            counts.to.use = paste0("clonotype_id_",clone.strategy))
  print(plots)
  }


  ##Select the top clonotpyes by user default
  if(clonotypes.per.sample){
    VGM.top.clonotypes <- tidyr::unite(VGM[[1]],sample_id_clonotype,sample_id,!!as.symbol(paste0("clonotype_id_",clone.strategy)),remove = FALSE)
    names.top.clonotypes <- paste0(rep("clonotype",top.clonotypes),1:top.clonotypes)
    VGM.top.clonotypes <- dplyr::filter(VGM.top.clonotypes,!!as.symbol(paste0("clonotype_id_",clone.strategy)) %in% names.top.clonotypes)
  }

 if(!clonotypes.per.sample){
   VGM.top.clonotypes <- tidyr::unite(VGM[[1]],sample_id_clonotype,sample_id,!!as.symbol(paste0("clonotype_id_",clone.strategy)),remove = FALSE)
   VGM.top.clonotypes <- dplyr::arrange(VGM.top.clonotypes,plyr::desc(!!as.symbol(paste0("clonotype_frequency_",clone.strategy))))
   unique.clonotypes <- unique(VGM.top.clonotypes$sample_id_clonotype)[1:top.clonotypes]
   VGM.top.clonotypes <- dplyr::filter(VGM.top.clonotypes,sample_id_clonotype %in% unique.clonotypes)
 }



  ## Run MIXCR alignment
  VDJ_mixcr_out <- Platypus::VDJ_call_MIXCR(VDJ = VGM.top.clonotypes ,
                                            mixcr.directory = mixcr.directory,
                                            species = species,
                                            platypus.version = platypus.version,
                                            operating.system = operating.system,
                                            simplify = simplify)

  VDJ_mixcr_out$VDJ_nt_mixcr <- paste0(VDJ_mixcr_out$VDJ_nSeqFR1,VDJ_mixcr_out$VDJ_nSeqCDR1, VDJ_mixcr_out$VDJ_nSeqFR2,VDJ_mixcr_out$VDJ_nSeqCDR2, VDJ_mixcr_out$VDJ_nSeqFR3,VDJ_mixcr_out$VDJ_nSeqCDR3, VDJ_mixcr_out$VDJ_nSeqFR4 )

  VDJ_mixcr_out$VJ_nt_mixcr <- paste0(VDJ_mixcr_out$VJ_nSeqFR1,VDJ_mixcr_out$VJ_nSeqCDR1, VDJ_mixcr_out$VJ_nSeqFR2,VDJ_mixcr_out$VJ_nSeqCDR2, VDJ_mixcr_out$VJ_nSeqFR3,VDJ_mixcr_out$VJ_nSeqCDR3, VDJ_mixcr_out$VJ_nSeqFR4 )

  VDJ_mixcr_out$VDJ_aa_mixcr <- paste0(VDJ_mixcr_out$VDJ_aaSeqFR1,VDJ_mixcr_out$VDJ_aaSeqCDR1, VDJ_mixcr_out$VDJ_aaSeqFR2,VDJ_mixcr_out$VDJ_aaSeqCDR2, VDJ_mixcr_out$VDJ_aaSeqFR3,VDJ_mixcr_out$VDJ_aaSeqCDR3, VDJ_mixcr_out$VDJ_aaSeqFR4 )

  VDJ_mixcr_out$VJ_aa_mixcr <- paste0(VDJ_mixcr_out$VJ_aaSeqFR1,VDJ_mixcr_out$VJ_aaSeqCDR1, VDJ_mixcr_out$VJ_aaSeqFR2,VDJ_mixcr_out$VJ_aaSeqCDR2, VDJ_mixcr_out$VJ_aaSeqFR3,VDJ_mixcr_out$VJ_aaSeqCDR3, VDJ_mixcr_out$VJ_aaSeqFR4 )



  #Select the most expanded antibodies per clonotype


  VGM.top.clonotypes.mixcr.out <- data.frame()
  for(ClonoTyp in unique(VDJ_mixcr_out$sample_id_clonotype)) {
    df <- dplyr::filter(VDJ_mixcr_out, sample_id_clonotype == ClonoTyp & Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)
    df <- df[stringr::str_count(df$VDJ_aa_mixcr,"[*_]") <= 1,]
    df <- df[stringr::str_count(df$VJ_aa_mixcr,"[*_]") <= 1,]
    df <- dplyr::arrange(df,plyr::desc("umis"))
    df <- dplyr::distinct(df,VDJ_aa_mixcr,VJ_aa_mixcr, .keep_all = TRUE)

    if(nrow(df) == 0) {}
    else if(nrow(df) < seq.per.clonotype) {VGM.top.clonotypes.mixcr.out <- dplyr::bind_rows(VGM.top.clonotypes.mixcr.out,df[1:nrow(df),])}
    else {VGM.top.clonotypes.mixcr.out <- dplyr::bind_rows(VGM.top.clonotypes.mixcr.out,df[1:seq.per.clonotype,])}
  }


  #VGM.top.clonotypes.mixcr.out$sample_id_clonotype <- NULL
  return(VGM.top.clonotypes.mixcr.out)


  } # End of function




