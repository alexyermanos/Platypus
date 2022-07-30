#'Utility function for bulk data to standard Platypus format conversion
#'
#' @description The VDJ_bulk_to_vgm function converts bulk output files from MIXCR or MAF into
#' a vgm-format compatible with most downstream Platypus functions
#' used for VDJ repertoire analysis.
#' @param VDJ.bulk.out.directory.list List containing paths to bulk VDJ output files from MIXCR or MAF. TRUST4 (and TRUST4.FULL) require an RDS file as input
#' @param input.type Character vector. Defaults to "MIXCR". "MIXCR", "MAF", "TRUST4", and "TRUST4.FULL" are supported. "TRUST4.FULL" contains TRUST additional columns, which were not originally supported by vgm: "cdr1", "cdr2", "v_cigar", "d_cigar", "j_cigar", "v_identity", "j_identity", "complete_vdj".
#' @param integrate.MIXCR.output Boolean. Defaults to TRUE. Whether to include in the VGM output additional MiXCR (49-78) columns.
#' @param vgm.expanded Boolean. Defaults to TRUE. Whether to include vgm[[9]] in the output list, where vgm[[9]] is the expanded version of vgm[[1]] having 1 line per read. For some Platypus functions, only vgm[[9]] (and not vgm[[1]]) may be compatible.
#' @param clone.strategy Character vector to specify the clonotyping strategy. Defaults to "cdr3.aa". Note that MIXCR input comes with clonotypes already assigned, and therefore clone.strategy should be specified only when the user wants to change the clonoyping strategy, and if no clone.strategy is provided, re-clonotyping will not be performed. Meanwhile, MAF inputs do not come with the clonotypes pre-assigned. Hence, if no clone.strategy is specified, "cdr3.aa" will be used as the default clonotyping strategy. The clonotyping strategies available in this function are: "cdr3.aa", "VDJJ.VJJ", "VDJJ.VJJ.cdr3length".
#' @param group.id Numeric vector. Defaults to NA. The user can specify to which group does each file belong to (e.g. a group could correspond to some specific treatment). The length of this numeric vector should match the number of samples in the VDJ.bulk.out.directory.list input.
#' @param cell.type Character vector. Defaults to NA. Cell type (e.g., "Bcell") of the MIXCR or MAF file that is provided as input.
#' @param batches Numeric vector. Defaults to NA. An additional grouping parameter that can be specified by the user. The length of this numeric vector should match the number of samples in the VDJ.bulk.out.directory.list input.
#' @param best.match.only Boolean. Whether only the highest scoring gene (V,J,D,C gene should) should be included in the output, or all matching genes in MIXCR should be included (MAF outputs: for the same read we can only have one possible V,J,D or C gene). Defaults to TRUE.
#' @return a VGM object (vgm.bulk.list). vgm.bulk.list[[1]]: each line correspond to a clonotype. vgm.bulk.list[[9]] (if vgm.expanded==TRUE): each line correspond to a read. The other (2-8) entries of the list are left empty for compatibility with Platypus functions.
#' @export
#' @examples
#' \dontrun{
#' Run from local directory using MIXCR/MAF bulk VDJ-repertoire files as inputs:
#' VDJ.bulk.out.directory.list <- list()
#' VDJ.bulk.out.directory.list[[1]] <- c("~/MIXCR_vdj_cdr3_clonotyping/C4.txt")
#' VDJ.bulk.out.directory.list[[2]] <- c("~/MIXCR_vdj_cdr3_clonotyping/C6.txt")
#' bulk.vgm.MIXCR <- VDJ_bulk_to_vgm(VDJ.bulk.out.directory.list = VDJ.bulk.out.directory.list,
#' input.type = 'MIXCR',
#' integrate.MIXCR.output = TRUE,
#' group.id = c(1,2),
#' cell.type = "Bcells",
#' batches = c(1,1),
#' vgm.expanded = TRUE,
#' best.match.only = FALSE)
#'
#' To re-clonotype MIXCR samples based on e.g., the CDR3 a.a. sequence:
#' bulk.vgm.MIXCR <- VDJ_bulk_to_vgm(VDJ.bulk.out.directory.list = VDJ.bulk.out.directory.list,
#' input.type = 'MIXCR',
#' integrate.MIXCR.output = TRUE,
#' group.id = c(1,2),
#' cell.type = "Bcells",
#' batches = c(1,1),
#' vgm.expanded = TRUE,
#' best.match.only = FALSE,
#' clone.strategy = "cdr3.aa")
#' }


VDJ_bulk_to_vgm <- function(VDJ.bulk.out.directory.list,
                        input.type,
                        integrate.MIXCR.output,
                        vgm.expanded,
                        clone.strategy,
                        group.id,
                        cell.type,
                        batches,
                        best.match.only) {

  #For CRAN checks
  cloneId <- NULL
  cloneID <- NULL
  nMutationsVRegion <- NULL
  nMutationsJRegion <- NULL
  nMutationsDRegion <- NULL
  V_count <- NULL
  J_count <- NULL
  D_count <- NULL
  SHM <- NULL
  SMH <- NULL


  # options(dplyr.summarise.inform = FALSE)

  if(missing(input.type)) {
    input.type <- "MIXCR"
    print("No input.type provided, standard: MIXCR")}

  if(missing(integrate.MIXCR.output)) {
    integrate.MIXCR.output <- TRUE}

  if(missing(vgm.expanded)) {
    vgm.expanded <- TRUE}

  if(missing(clone.strategy)) {
    want.clonotyping = FALSE
    clone.strategy = "cdr3.aa" #MAF does not give clonotype IDs, so the clonotyping function is always called in MAF standard is cdr3.aa  want.clonotyping = FALSE
    if(input.type=="MAF" | input.type=="") message("No clone.strategy provided. Data will be clonotyped by CDR3")
    if(input.type=="MIXCR") message("No clone.strategy provided. 10x clonotyping will be kept")
    if(input.type=="TRUST4" | input.type=="TRUST4.FULL") message("No clone.strategy provided. Sequence_id values will be used as clonotypes")
    }
  else {want.clonotyping = TRUE}

  if(missing(group.id)){
    group.id <- NA}

  if(missing(cell.type)){
    cell.type <- NA}

  if(missing(batches)){
    batches <- NA}

  if(length(group.id) != length(VDJ.bulk.out.directory.list) && !is.na(group.id)) {
    stop("numeric group.id vector should be of the same length as
          the number of samples in VDJ.bulk.out.directory.list.")}

  if(length(batches) != length(VDJ.bulk.out.directory.list) && !is.na(batches)) {
    stop("numeric batches vector should be of the same length as
          the number of samples in VDJ.bulk.out.directory.list.")}

  if(missing(best.match.only)){
    best.match.only <- TRUE
  }

  if(input.type == 'MIXCR'){
    VDJ.dataframe.list <- list()
    data.df <- data.frame()
    # Store data frames of MIXCR-outputs in the VDJ.dataframe.list;
    # combine the data frames in data_df; sample ids ("s1", "s2" ...) are assigned based
    # on the user-defined order in the VDJ.bulk.out.directory.list.
    # group ids & batches are assigned based on user-defined vector input to the function (defaults to NA).
    for(i in 1:length(VDJ.bulk.out.directory.list)) {
      VDJ.dataframe.list[[i]] <- utils::read.table(VDJ.bulk.out.directory.list[[i]], sep="\t", header = TRUE)
      VDJ.dataframe.list[[i]]$sample_id <- paste0("s", i)
      VDJ.dataframe.list[[i]]$group_id <- group.id[i]
      VDJ.dataframe.list[[i]]$batches <- batches[i]
      data.df <- rbind(data.df, VDJ.dataframe.list[[i]])
    }

    # Then instead of using a for-loop to fill in the vgm[[1]], vectorized operations are used to
    # speed up the computation time.

    # The first 48 columns exactly match the formatting of vgm[[1]] generated after running
    # VDJ_GEX_matrix for single-cell data. In case of unavailable data, NA is used.

    i <- substr(data.df$allCHitsWithScore, 1, 3) %in% c("IGH", "") # subset the heavy chains
    j <- substr(data.df$allCHitsWithScore, 1, 3) %in% c("IGL", "IGK") # subset the light chains


    #1:
    barcode <- rep(NA, nrow(data.df))

    #2:
    sample_id <- data.df$sample_id

    #3:
    group_id <- data.df$group_id

    #4:
    clonotype_id_10x <- data.df$cloneId
    clonotype_id_10x <- paste0("clonotype", clonotype_id_10x)

    #5:
    celltype <- rep(cell.type, nrow(data.df))

    #6:
    Nr_of_VDJ_chains <- rep(NA, nrow(data.df))
    Nr_of_VDJ_chains[i] <- 1

    #7:
    Nr_of_VJ_chains <- rep(NA, nrow(data.df))
    Nr_of_VJ_chains[j]<- 1

    #8:
    VDJ_cdr3s_aa <- rep(NA, nrow(data.df))
    VDJ_cdr3s_aa[i] <- data.df$aaSeqCDR3[i]

    #9:
    VJ_cdr3s_aa <- rep(NA, nrow(data.df))
    VJ_cdr3s_aa[j] <- data.df$aaSeqCDR3[j]

    #10:
    VDJ_cdr3s_nt <- rep(NA, nrow(data.df))
    VDJ_cdr3s_nt[i] <- data.df$nSeqCDR3[i]

    #11:
    VJ_cdr3s_nt <- rep(NA, nrow(data.df))
    VJ_cdr3s_nt[j] <- data.df$nSeqCDR3[j]

    #12:
    VDJ_chain_contig = rep(NA, nrow(data.df))

    #13:
    VJ_chain_contig = rep(NA, nrow(data.df))

    #14:
    VDJ_chain <- rep(NA, nrow(data.df))
    VDJ_chain[i] <- substr(data.df$allCHitsWithScore, 1, 3)[i]

    #15:
    VJ_chain <- rep(NA, nrow(data.df))
    VJ_chain[j] <- substr(data.df$allCHitsWithScore, 1, 3)[j]

    #16:
    VDJ_vgene <- rep(NA, nrow(data.df))
    VDJ_vgene[i] <- gsub("\\*[^,]*", "", data.df$allVHitsWithScore[i])
    VDJ_vgene[i] <- gsub(",", ";", VDJ_vgene[i])

    #17:
    VJ_vgene <- rep(NA, nrow(data.df))
    VJ_vgene[j] <- gsub("\\*[^,]*", "", data.df$allVHitsWithScore[j])
    VJ_vgene[j] <- gsub(",", ";", VJ_vgene[j])

    #18:
    VDJ_dgene <- rep(NA, nrow(data.df))
    if("allDHitsWithScore" %in% colnames(data.df)) {
      VDJ_dgene[i] <- gsub("\\*[^,]*", "", data.df$allDHitsWithScore[i])
      VDJ_dgene[i] <- gsub(",", ";", VDJ_dgene[i])}

    #19:
    VDJ_jgene <- rep(NA, nrow(data.df))
    if("allJHitsWithScore" %in% colnames(data.df)) {
      VDJ_jgene[i] <- gsub("\\*[^,]*", "", data.df$allJHitsWithScore[i])
      VDJ_jgene[i] <- gsub(",", ";", VDJ_jgene[i])}

    #20:
    VJ_jgene <- rep(NA, nrow(data.df))
    if("allJHitsWithScore" %in% colnames(data.df)) {
      VJ_jgene[j] <- gsub("\\*[^,]*", "", data.df$allJHitsWithScore[j])
      VJ_jgene[j] <- gsub(",", ";", VJ_jgene[j])}

    #21:
    VDJ_cgene <- rep(NA, nrow(data.df))
    VDJ_cgene[i] <- gsub("\\*[^,]*", "", data.df$allCHitsWithScore[i])
    VDJ_cgene[i] <- gsub(",", ";", VDJ_cgene[i])

    #22:
    VJ_cgene <- rep(NA, nrow(data.df))
    VJ_cgene[j] <- gsub("\\*[^,]*", "", data.df$allCHitsWithScore[j])
    VJ_cgene[j] <- gsub(",", ";", VJ_cgene[j])

    #23:
    VDJ_sequence_nt_raw <- rep(NA, nrow(data.df))

    #24:
    VJ_sequence_nt_raw <- rep(NA, nrow(data.df))

    #25:
    VDJ_sequence_nt_trimmed = rep(NA, nrow(data.df))
    names_n <- c("nSeqFR1", "nSeqCDR1", "nSeqFR2", "nSeqCDR2",
                 "nSeqFR3", "nSeqCDR3", "nSeqFR4")
    if(all(names_n %in% colnames(data.df))){
      VDJ_sequence_nt_trimmed[i] <- paste0(data.df$nSeqFR1[i],
                                           data.df$nSeqCDR1[i],
                                           data.df$nSeqFR2[i],
                                           data.df$nSeqCDR2[i],
                                           data.df$nSeqFR3[i],
                                           data.df$nSeqCDR3[i],
                                           data.df$nSeqFR4[i])}

    #26:
    VJ_sequence_nt_trimmed = rep(NA, nrow(data.df))
    if(all(names_n %in% colnames(data.df))){
      VJ_sequence_nt_trimmed[j] <- paste0(data.df$nSeqFR1[j],
                                          data.df$nSeqCDR1[j],
                                          data.df$nSeqFR2[j],
                                          data.df$nSeqCDR2[j],
                                          data.df$nSeqFR3[j],
                                          data.df$nSeqCDR3[j],
                                          data.df$nSeqFR4[j])}

    #27:
    VDJ_sequence_aa <- rep(NA, nrow(data.df))
    names_a <- c("aaSeqFR1", "aaSeqCDR1", "aaSeqFR2", "aaSeqCDR2", "aaSeqFR3",
                 "aaSeqCDR3", "aaSeqFR4")
    if(all(names_a %in% colnames(data.df))){
      VDJ_sequence_aa[i] <- paste0(data.df$aaSeqFR1[i],
                                   data.df$aaSeqCDR1[i], data.df$aaSeqFR2[i],
                                   data.df$aaSeqCDR2[i], data.df$aaSeqFR3[i],
                                   data.df$aaSeqCDR3[i], data.df$aaSeqFR4[i])}

    #28:
    VJ_sequence_aa <- rep(NA, nrow(data.df))
    if(all(names_a %in% colnames(data.df))){
      VJ_sequence_aa[j] <- paste0(data.df$aaSeqFR1[j],
                                  data.df$aaSeqCDR1[j], data.df$aaSeqFR2[j],
                                  data.df$aaSeqCDR2[j], data.df$aaSeqFR3[j],
                                  data.df$aaSeqCDR3[j], data.df$aaSeqFR4[j])}

    #29:
    VDJ_trimmed_ref <- rep(NA, nrow(data.df))
    if(all(names_n %in% colnames(data.df))){
      VDJ_trimmed_ref[i] <- paste0(data.df$nSeqFR1[i],
                                   data.df$nSeqCDR1[i], data.df$nSeqFR2[i],
                                   data.df$nSeqCDR2[i], data.df$nSeqFR3[i],
                                   data.df$nSeqCDR3[i], data.df$nSeqFR4[i])}

    #30:
    VJ_trimmed_ref <- rep(NA, nrow(data.df))
    if(all(names_n %in% colnames(data.df))){
      VJ_trimmed_ref[j] <- paste0(data.df$nSeqFR1[j],
                                  data.df$nSeqCDR1[j], data.df$nSeqFR2[j],
                                  data.df$nSeqCDR2[j], data.df$nSeqFR3[j],
                                  data.df$nSeqCDR3[j], data.df$nSeqFR4[j])}

    #31:
    VDJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #32:
    VJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #33:
    orig_barcode <- rep(NA, nrow(data.df))

    #34:
    clonotype_frequency <- data.df$cloneCount

    #35:
    specifity <- rep(NA, nrow(data.df))

    #36:
    affinity <- rep(NA, nrow(data.df))

    #37:
    GEX_available <- rep("FALSE", nrow(data.df))

    #38:
    orig.ident <- rep(NA, nrow(data.df))

    #39:
    orig_barcode_GEX  <- rep(NA, nrow(data.df))

    #40:
    seurat_clusters <- rep(NA, nrow(data.df))

    #41:
    PC_1  <-  rep(NA, nrow(data.df))

    #42:
    PC_2 <-  rep(NA, nrow(data.df))

    #43:
    UMAP_1 = rep(NA, nrow(data.df))

    #44:
    UMAP_2 <- rep(NA, nrow(data.df))

    #45:
    tSNE_1 <- rep(NA, nrow(data.df))

    #46:
    tSNE_2 <- rep(NA, nrow(data.df))

    #47:
    batches <- data.df$batches

    #48:
    clonotype_id <- data.df$cloneId
    clonotype_id <- paste0("clonotype", clonotype_id)

    # Columns 49-78 below match the 30 additional columns that are added to
    # vgm[[1]] after running VDJ_call_MIXCR:

    #49:
    VDJ_nSeqFR1 <- rep(NA, nrow(data.df))
    if("nSeqFR1" %in% colnames(data.df)) {
      VDJ_nSeqFR1[i] <- data.df$nSeqFR1[i]}

    #50:
    VDJ_nSeqCDR1 <- rep(NA, nrow(data.df))
    if("nSeqCDR1" %in% colnames(data.df)) {
      VDJ_nSeqCDR1[i] <- data.df$nSeqCDR1[i]}

    #51:
    VDJ_nSeqFR2 <- rep(NA, nrow(data.df))
    if("nSeqFR2" %in% colnames(data.df)) {
      VDJ_nSeqFR2[i] <- data.df$nSeqFR2[i]}

    #52:
    VDJ_nSeqCDR2 <- rep(NA, nrow(data.df))
    if("nSeqCDR2" %in% colnames(data.df)) {
      VDJ_nSeqCDR2[i] <- data.df$nSeqCDR2[i]}

    #53:
    VDJ_nSeqFR3 <- rep(NA, nrow(data.df))
    if("nSeqFR3" %in% colnames(data.df)) {
      VDJ_nSeqFR3[i] <- data.df$nSeqFR3[i]}

    #54:
    VDJ_nSeqCDR3 <- rep(NA, nrow(data.df))
    VDJ_nSeqCDR3[i] <- data.df$nSeqCDR3[i] # CDR3 is the minimum information we expect to always have

    #55:
    VDJ_nSeqFR4 <- rep(NA, nrow(data.df))
    if("nSeqFR4" %in% colnames(data.df)) {
      VDJ_nSeqFR4[i] <- data.df$nSeqFR4[i]}

    #56:
    VDJ_aaSeqFR1 <- rep(NA, nrow(data.df))
    if("aaSeqFR1" %in% colnames(data.df)) {
      VDJ_aaSeqFR1[i] <- data.df$aaSeqFR1[i]}

    #57:
    VDJ_aaSeqCDR1 <- rep(NA, nrow(data.df))
    if("aaSeqCDR1" %in% colnames(data.df)) {
      VDJ_aaSeqCDR1[i] <- data.df$aaSeqCDR1[i]}

    #58:
    VDJ_aaSeqFR2 <- rep(NA, nrow(data.df))
    if("aaSeqFR2" %in% colnames(data.df)) {
      VDJ_aaSeqFR2[i] <- data.df$aaSeqFR2[i]}

    #59:
    VDJ_aaSeqCDR2 <- rep(NA, nrow(data.df))
    if("aaSeqCDR2" %in% colnames(data.df)) {
      VDJ_aaSeqCDR2[i] <- data.df$aaSeqCDR2[i]}

    #60:
    VDJ_aaSeqFR3 <- rep(NA, nrow(data.df))
    if("aaSeqFR3" %in% colnames(data.df)) {
      VDJ_aaSeqFR3[i] <- data.df$aaSeqFR3[i]}

    #61:
    VDJ_aaSeqCDR3 <- rep(NA, nrow(data.df))
    VDJ_aaSeqCDR3[i] <- data.df$aaSeqCDR3[i]

    #62:
    VDJ_aaSeqFR4 <- rep(NA, nrow(data.df))
    if("aaSeqFR4" %in% colnames(data.df)) {
      VDJ_aaSeqFR4[i] <- data.df$aaSeqFR4[i]}

    #63:
    VDJ_SHM <- rep(NA, nrow(data.df))
    names_SMH <- c("nMutationsVRegion", "nMutationsJRegion",
                   "nMutationsDRegion")
    if(all(names_SMH  %in% colnames(data.df))) {
      VDJ_SHM[i] <- data.df[i,] %>%
        dplyr::select(cloneId, nMutationsVRegion,
                      nMutationsJRegion, nMutationsDRegion) %>%
        dplyr::mutate(V_count = stringr::str_count(nMutationsVRegion, "\\w+"),
                      J_count = stringr::str_count(nMutationsJRegion, "\\w+"),
                      D_count = stringr::str_count(nMutationsDRegion, "\\w+"),
                      SMH = V_count + J_count + D_count) %>% dplyr::pull(SMH)
    }

    #64:
    VJ_nSeqFR1 <- rep(NA, nrow(data.df))
    if("nSeqFR1" %in% colnames(data.df)){
      VJ_nSeqFR1[j] <- data.df$nSeqFR1[j]}

    #65:
    VJ_nSeqCDR1 <- rep(NA, nrow(data.df))
    if("nSeqCDR1" %in% colnames(data.df)){
      VJ_nSeqCDR1[j] <- data.df$nSeqCDR1[j]}

    #66:
    VJ_nSeqFR2 <- rep(NA, nrow(data.df))
    if("nSeqFR2" %in% colnames(data.df)){
      VJ_nSeqFR2[j] <- data.df$nSeqFR2[j]}

    #67:
    VJ_nSeqCDR2 <- rep(NA, nrow(data.df))
    if("nSeqCDR2" %in% colnames(data.df)){
      VJ_nSeqCDR2[j] <- data.df$nSeqCDR2[j]}

    #68:
    VJ_nSeqFR3 <- rep(NA, nrow(data.df))
    if("nSeqFR3" %in% colnames(data.df)){
      VJ_nSeqFR3[j] <- data.df$nSeqFR3[j]}

    #69:
    VJ_nSeqCDR3 <- rep(NA, nrow(data.df))
    VJ_nSeqCDR3[j] <- data.df$nSeqCDR3[j]

    #70:
    VJ_nSeqFR4 <- rep(NA, nrow(data.df))
    if("nSeqFR4" %in% colnames(data.df)){
      VJ_nSeqFR4[j] <- data.df$nSeqFR4[j]}

    #71:
    VJ_aaSeqFR1 <- rep(NA, nrow(data.df))
    if("aaSeqFR1" %in% colnames(data.df)){
      VJ_aaSeqFR1[j] <- data.df$aaSeqFR1[j]}

    #72:
    VJ_aaSeqCDR1 <- rep(NA, nrow(data.df))
    if("aaSeqCDR1" %in% colnames(data.df)){
      VJ_aaSeqCDR1[j] <- data.df$aaSeqCDR1[j]}

    #73:
    VJ_aaSeqFR2 <- rep(NA, nrow(data.df))
    if("aaSeqFR2" %in% colnames(data.df)){
      VJ_aaSeqFR2[j] <- data.df$aaSeqFR2[j]}

    #74:
    VJ_aaSeqCDR2 <- rep(NA, nrow(data.df))
    if("aaSeqCDR2" %in% colnames(data.df)){
      VJ_aaSeqCDR2[j] <- data.df$aaSeqCDR2[j]}

    #75:
    VJ_aaSeqFR3 <- rep(NA, nrow(data.df))
    if("aaSeqFR3" %in% colnames(data.df)){
      VJ_aaSeqFR3[j] <- data.df$aaSeqFR3[j]}

    #76:
    VJ_aaSeqCDR3 <- rep(NA, nrow(data.df))
    VJ_aaSeqCDR3[j] <- data.df$aaSeqCDR3[j]

    #77:
    VJ_aaSeqFR4 <- rep(NA, nrow(data.df))
    if("aaSeqFR4" %in% colnames(data.df)){
      VJ_aaSeqFR4[j] <- data.df$aaSeqFR4[j]}

    #78:
    VJ_SHM <- rep(NA, nrow(data.df))
    if(all(names_SMH %in% colnames(data.df))) {
      VJ_SHM[j] <- data.df[j,] %>%
        dplyr::select(cloneId, nMutationsVRegion,
                      nMutationsJRegion, nMutationsDRegion) %>%
        dplyr::mutate(V_count = stringr::str_count(nMutationsVRegion, "\\w+"),
                      J_count = stringr::str_count(nMutationsJRegion, "\\w+"),
                      D_count = stringr::str_count(nMutationsDRegion, "\\w+"),
                      SMH = V_count + J_count + D_count) %>% dplyr::pull(SMH)
    }

    if (best.match.only == TRUE){
      VDJ_vgene <- gsub("\\;.*", "", VDJ_vgene)
      VJ_vgene <- gsub("\\;.*", "", VJ_vgene)
      VDJ_dgene <- gsub("\\;.*", "", VDJ_dgene)
      VDJ_jgene <- gsub("\\;.*", "", VDJ_jgene)
      VJ_jgene <- gsub("\\;.*", "", VJ_jgene)
      VDJ_cgene <- gsub("\\;.*", "", VDJ_cgene)
      VJ_cgene <- gsub("\\;.*", "", VJ_cgene)
    }

    # Combining all 78 lists into one data frame:
    vgm.bulk <- data.frame(barcode, sample_id, group_id, clonotype_id_10x,
                           celltype, Nr_of_VDJ_chains, Nr_of_VJ_chains,
                           VDJ_cdr3s_aa, VJ_cdr3s_aa, VDJ_cdr3s_nt,
                           VJ_cdr3s_nt, VDJ_chain_contig, VJ_chain_contig,
                           VDJ_chain, VJ_chain, VDJ_vgene, VJ_vgene,VDJ_dgene,
                           VDJ_jgene, VJ_jgene, VDJ_cgene, VJ_cgene,
                           VDJ_sequence_nt_raw, VJ_sequence_nt_raw,
                           VDJ_sequence_nt_trimmed, VJ_sequence_nt_trimmed,
                           VDJ_sequence_aa, VJ_sequence_aa, VDJ_trimmed_ref,
                           VJ_trimmed_ref, VDJ_raw_consensus_id,
                           VJ_raw_consensus_id, orig_barcode,
                           clonotype_frequency, specifity, affinity,
                           GEX_available, orig.ident, orig_barcode_GEX,
                           seurat_clusters, PC_1, PC_2, UMAP_1, UMAP_2,
                           tSNE_1, tSNE_2, batches, clonotype_id,
                           VDJ_nSeqFR1, VDJ_nSeqCDR1, VDJ_nSeqFR2,
                           VDJ_nSeqCDR2, VDJ_nSeqFR3, VJ_nSeqCDR3,
                           VDJ_nSeqFR4, VDJ_aaSeqFR1, VDJ_aaSeqCDR1,
                           VDJ_aaSeqFR2, VDJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VDJ_aaSeqCDR3, VDJ_aaSeqFR4, VDJ_SHM,
                           VJ_nSeqFR1, VJ_nSeqCDR1, VJ_nSeqFR2,
                           VJ_nSeqCDR2, VJ_nSeqFR3, VJ_nSeqCDR3,
                           VJ_nSeqFR4, VJ_aaSeqFR1, VJ_aaSeqCDR1,
                           VJ_aaSeqFR2, VJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VJ_aaSeqCDR3, VJ_aaSeqFR4, VJ_SHM)


    # Change all white spaces to NA:
    # vgm.bulk[vgm.bulk == ""] <- NA

    # Return columns based on whether the user specifies integrate.MIXCR.output to be T(rue) or F(false).
    if (integrate.MIXCR.output == TRUE) {
      vgm.bulk <- vgm.bulk[, 1:78]
    }  else if (integrate.MIXCR.output == FALSE) {
      vgm.bulk <- vgm.bulk[, 1:48]
    }

    # create an empty list vgm.bulk.list in we will store our two data frames:
    # vgm.bulk.list[[1]] contains aggregated clonotype-frequencies stored in a single clonotype_frequency row.

    # vgm.bulk.list[[9]] will contain "per-cell" information: rows will be added for each bulk clonotype based on
    # their aggregated clonotype_frequency number. I.e.: if we have clonotype0 with 8000 cells, vgm.bulk.list[[9]]
    # will contain 8000 rows for each cell of clonotype0. This is done to achieve better compatibility with some
    # pre-existing Platypus functions used for analyzing sc-VDJ sequence data.

    vgm.bulk.list <- list()
    vgm.bulk.list[[1]] <- vgm.bulk

    if (vgm.expanded == TRUE) {
      # creating vgm.bulk[[9]]:
      vgm.per.cell <- vgm.bulk %>%
        utils::type.convert(as.is = TRUE) %>%
        tidyr::uncount(clonotype_frequency, .remove = FALSE)
      vgm.bulk.list[[9]] <- vgm.per.cell
    }

    # if we want to re-clonotype:
    if(want.clonotyping == TRUE & vgm.expanded == TRUE) {
      vgm.per.cell <- Platypus::VDJ_clonotype(VDJ = vgm.bulk.list[[9]],
                                                 clone.strategy = clone.strategy,
                                                 global.clonotype = F,
                                                 VDJ.VJ.1chain = F,
                                                 hierarchical = 'none')

      # overwrite pre-existing columns created on the previous clonotyping strategy:
      clonotype_id_ <- paste0("clonotype_id_", clone.strategy)
      clonotype_frequency_ <- paste0("clonotype_frequency_", clone.strategy)

      vgm.per.cell$clonotype_id <- vgm.per.cell[[clonotype_id_]]
      vgm.per.cell$clonotype_id_10x <- vgm.per.cell[[clonotype_id_]]
      vgm.per.cell$clonotype_frequency <-  vgm.per.cell[[clonotype_frequency_]]

      vgm.bulk.list[[9]] <- vgm.per.cell[, 1:78]
      vgm.bulk.list[[1]] <- vgm.bulk.list[[9]] %>%
        dplyr::group_by(clonotype_id, sample_id) %>%
        dplyr::distinct()
    } else if(want.clonotyping == TRUE & vgm.expanded == FALSE){
      print("Re-clonotyping can only be carried out if both want.clonotyping and vgm.expanded arguments are set to TRUE")
    }

  } else if (input.type == "MAF"){

    # Store data frames of MAF-outputs in the VDJ.dataframe.list:
    VDJ.dataframe.list <- list()
    data.df <- data.frame()
    for(i in 1:length(VDJ.bulk.out.directory.list)) {
      VDJ.dataframe.list[[i]] <-  utils::read.delim(VDJ.bulk.out.directory.list[[i]])
      VDJ.dataframe.list[[i]]$sample_id <- paste0("s", i)
      VDJ.dataframe.list[[i]]$group_id <- group.id[i]
      VDJ.dataframe.list[[i]]$batches <- batches[i]
      data.df <- rbind(data.df, VDJ.dataframe.list[[i]])
    }


    # initialize an empty data frame to store vgm[[1]]:
    vgm.bulk <- data.frame()

    # retrieving and appending information to the vgm.bulk data frame:

    #1
    barcode = data.df$Rev.BC_After_Correction

    #2 MAF outputs are missing the sample_id column; therefore, we assign sample ids ("s1", "s2" ...) based on the user-defined order in the VDJ.bulk.out.directory.list:
    sample_id <- data.df$sample_id


    #3 group_id column; group.id vector with user-defined integers specifying the group membership (defaults to NA):
    group_id <- data.df$group_id


    #4 clonotype_id_10x column:
    clonotype_id_10x <- rep(NA, nrow(data.df))

    #5
    celltype <- rep(cell.type, nrow(data.df))

    #6
    Nr_of_VDJ_chains <- rep(1, nrow(data.df))

    #7
    Nr_of_VJ_chains <- rep(0, nrow(data.df))

    #8 VDJ_cdr3s_aa column:
    VDJ_cdr3s_aa = data.df$CDR3.AA

    #9
    VJ_cdr3s_aa = rep(NA, nrow(data.df))

    #10 VDJ_cdr3s_nt column:
    VDJ_cdr3s_nt <- rep(NA, nrow(data.df))

    #11
    VJ_cdr3s_nt <- rep(NA, nrow(data.df))

    #12
    VDJ_chain_contig <- rep(NA, nrow(data.df))

    #13
    VJ_chain_contig <- rep(NA, nrow(data.df))

    #14
    VDJ_chain <- substr(data.df$VGene, 1,3)

    #15
    VJ_chain <- rep(NA, nrow(data.df))

    #16 VDJ_vgene column:
    VDJ_vgene <- data.df$VGene

    #17
    VJ_vgene <- rep(NA, nrow(data.df))

    #18 VDJ_dgene column:
    VDJ_dgene <- rep(NA, nrow(data.df))

    #19 VDJ_jgene column:
    VDJ_jgene <- data.df$JGene

    #20
    VJ_jgene <- rep(NA, nrow(data.df))

    #21
    VDJ_cgene <- toupper(paste0(substr(data.df$Isotype,2,3),"H",substr(data.df$Isotype,4,4)))

    #22
    VJ_cgene <- rep(NA, nrow(data.df))

    #23 VDJ_sequence_nt_raw column:
    VDJ_sequence_nt_raw <- rep(NA, nrow(data.df))

    #24
    VJ_sequence_nt_raw <- rep(NA, nrow(data.df))

    #25
    VDJ_sequence_nt_trimmed = data.df$Corrected_NucleotideSeq

    #26
    VJ_sequence_nt_trimmed <- rep(NA, nrow(data.df))

    #27 VDJ_sequence_aa column:
    VDJ_sequence_aa <- data.df$Corrected_VDJ.AA

    #28
    VJ_sequence_aa <- rep(NA, nrow(data.df))

    #29
    VDJ_trimmed_ref = rep(NA, nrow(data.df))

    #30
    VJ_trimmed_ref <- rep(NA, nrow(data.df))

    #31
    VDJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #32
    VJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #33
    orig_barcode <- rep(NA, nrow(data.df))

    #34
    clonotype_frequency <- rep(NA, nrow(data.df))

    #35
    specifity <- rep(NA, nrow(data.df))

    #36
    affinity <- rep(NA, nrow(data.df))

    #37
    GEX_available <- rep("FALSE", nrow(data.df))

    #38
    orig.ident <- rep(NA, nrow(data.df))

    #39
    orig_barcode_GEX <- rep(NA, nrow(data.df))
    #40
    seurat_clusters <- rep(NA, nrow(data.df))

    #41
    PC_1 <- rep(NA, nrow(data.df))

    #42
    PC_2 <- rep(NA, nrow(data.df))

    #43
    UMAP_1 <- rep(NA, nrow(data.df))

    #44
    UMAP_2 <- rep(NA, nrow(data.df))

    #45
    tSNE_1 <- rep(NA, nrow(data.df))

    #46
    tSNE_2 <- rep(NA, nrow(data.df))

    #47
    batches <- data.df$batches

    #48
    clonotype_id <- rep(NA, nrow(data.df))

    #49:
    VDJ_nSeqFR1 <- rep(NA, nrow(data.df))

    #50:
    VDJ_nSeqCDR1 <- rep(NA, nrow(data.df))

    #51:
    VDJ_nSeqFR2 <- rep(NA, nrow(data.df))

    #52:
    VDJ_nSeqCDR2 <- rep(NA, nrow(data.df))

    #53:
    VDJ_nSeqFR3 <- rep(NA, nrow(data.df))

    #54:
    VDJ_nSeqCDR3 <- rep(NA, nrow(data.df))

    #55:
    VDJ_nSeqFR4 <- rep(NA, nrow(data.df))

    #56:
    VDJ_aaSeqFR1 <- rep(NA, nrow(data.df))

    #57:
    VDJ_aaSeqCDR1 <- rep(NA, nrow(data.df))

    #58:
    VDJ_aaSeqFR2 <- rep(NA, nrow(data.df))

    #59:
    VDJ_aaSeqCDR2 <-rep(NA, nrow(data.df))

    #60:
    VDJ_aaSeqFR3 <- rep(NA, nrow(data.df))

    #61:
    VDJ_aaSeqCDR3 <- rep(NA, nrow(data.df))

    #62:
    VDJ_aaSeqFR4 <- rep(NA, nrow(data.df))

    #63:
    VDJ_SHM <- data.df$CB_VRegion_SHM_Tot

    #64:
    VJ_nSeqFR1 <- rep(NA, nrow(data.df))

    #65:
    VJ_nSeqCDR1 <- rep(NA, nrow(data.df))

    #66:
    VJ_nSeqFR2 <- rep(NA, nrow(data.df))

    #67:
    VJ_nSeqCDR2 <- rep(NA, nrow(data.df))

    #68:
    VJ_nSeqFR3 <- rep(NA, nrow(data.df))

    #69:
    VJ_nSeqCDR3 <- rep(NA, nrow(data.df))

    #70:
    VJ_nSeqFR4 <- rep(NA, nrow(data.df))

    #71:
    VJ_aaSeqFR1 <- rep(NA, nrow(data.df))

    #72:
    VJ_aaSeqCDR1 <- rep(NA, nrow(data.df))

    #73:
    VJ_aaSeqFR2 <- rep(NA, nrow(data.df))

    #74:
    VJ_aaSeqCDR2 <- rep(NA, nrow(data.df))

    #75:
    VJ_aaSeqFR3 <- rep(NA, nrow(data.df))

    #76:
    VJ_aaSeqCDR3 <- rep(NA, nrow(data.df))

    #77:
    VJ_aaSeqFR4 <- rep(NA, nrow(data.df))

    #78:
    VJ_SHM <- rep(NA, nrow(data.df))


    # start of appending the required information to a local data frame. The first 48 columns that
    # are added do this data frame exactly match the formatting of vgm[[1]] generated after running
    # VDJ_GEX_matrix for single-cell data. In case of unavailable data, NA is used.

    vgm.bulk <- data.frame(barcode, sample_id, group_id, clonotype_id_10x,
                           celltype, Nr_of_VDJ_chains, Nr_of_VJ_chains,
                           VDJ_cdr3s_aa, VJ_cdr3s_aa, VDJ_cdr3s_nt,
                           VJ_cdr3s_nt, VDJ_chain_contig, VJ_chain_contig,
                           VDJ_chain, VJ_chain, VDJ_vgene, VJ_vgene,VDJ_dgene,
                           VDJ_jgene, VJ_jgene, VDJ_cgene, VJ_cgene,
                           VDJ_sequence_nt_raw, VJ_sequence_nt_raw,
                           VDJ_sequence_nt_trimmed, VJ_sequence_nt_trimmed,
                           VDJ_sequence_aa, VJ_sequence_aa, VDJ_trimmed_ref,
                           VJ_trimmed_ref, VDJ_raw_consensus_id,
                           VJ_raw_consensus_id, orig_barcode,
                           clonotype_frequency, specifity, affinity,
                           GEX_available, orig.ident, orig_barcode_GEX,
                           seurat_clusters, PC_1, PC_2, UMAP_1, UMAP_2,
                           tSNE_1, tSNE_2, batches, clonotype_id,
                           VDJ_nSeqFR1, VDJ_nSeqCDR1, VDJ_nSeqFR2,
                           VDJ_nSeqCDR2, VDJ_nSeqFR3, VJ_nSeqCDR3,
                           VDJ_nSeqFR4, VDJ_aaSeqFR1, VDJ_aaSeqCDR1,
                           VDJ_aaSeqFR2, VDJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VDJ_aaSeqCDR3, VDJ_aaSeqFR4, VDJ_SHM,
                           VJ_nSeqFR1, VJ_nSeqCDR1, VJ_nSeqFR2,
                           VJ_nSeqCDR2, VJ_nSeqFR3, VJ_nSeqCDR3,
                           VJ_nSeqFR4, VJ_aaSeqFR1, VJ_aaSeqCDR1,
                           VJ_aaSeqFR2, VJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VJ_aaSeqCDR3, VJ_aaSeqFR4, VJ_SHM)


    if (integrate.MIXCR.output == TRUE) {
      vgm.bulk <- vgm.bulk[, 1:78]
    }
    else if (integrate.MIXCR.output == FALSE) {
      vgm.bulk <- vgm.bulk[, 1:48]}

    # create an empty list vgm.bulk.list in we will store our two data frames:
    # vgm.bulk.list[[1]] contains aggregated clonotype-frequencies stored in a single clonotype_frequency row.
    # vgm.bulk.list[[9]] will contain "per-cell" information: rows will be added for each bulk clonotype based on
    # their aggregated clonotype_frequency number. I.e.: if we have clonotype0 with 8000 cells, vgm.bulk.list[[9]] will
    # contain 8000 rows for each cell of clonotype0. This is done to achieve better compatibility with some pre-existing Platypus
    # functions used for analyzing sc-VDJ sequence data.

    vgm.bulk.list <- list()

    vgm.clonotyped <- VDJ_clonotype(VDJ = vgm.bulk, clone.strategy = clone.strategy,
                                       global.clonotype = F, VDJ.VJ.1chain = F, hierarchical = 'none') # TRY ALSO hierarchical = 'single.chains':
    # Not filtering cells with counts other than 1VDJ 1VJ chain and integrating these cells hierarchically into clonotypes

    ### LINES BELOW NO LONGER NECESSARY THANKD TO VDJ_clonotype UPDATE
    # clonotype_id_ <- paste0("clonotype_id_", clone.strategy)
    # clonotype_frequency_ <- paste0("clonotype_frequency_", clone.strategy)
    # vgm.clonotyped$clonotype_id <- vgm.clonotyped[[clonotype_id_]]
    # vgm.clonotyped$clonotype_id_10x <- vgm.clonotyped[[clonotype_id_]]
    # vgm.clonotyped$clonotype_frequency <-  vgm.clonotyped[[clonotype_frequency_]]

    ### HERE THE OLD VERSION OF vgm[[9]] -> vgm[[1]] compression
    #
    # vgm.bulk.list[[1]] <- vgm.clonotyped %>%
    #   add_count(VDJ_vgene, VDJ_jgene, VDJ_sequence_nt_trimmed, VDJ_sequence_aa, VDJ_cgene, VDJ_SHM, VDJ_cdr3s_aa) %>%
    #   dplyr::group_by(clonotype_id, sample_id) %>%
    #   dplyr::mutate(VDJ_vgene = VDJ_vgene[n == max(n)][1],
    #                 VDJ_jgene = VDJ_jgene[n == max(n)][1],
    #                 VDJ_sequence_nt_trimmed = VDJ_sequence_nt_trimmed[n == max(n)][1],
    #                 VDJ_sequence_aa = VDJ_sequence_aa[n == max(n)][1],
    #                 VDJ_cgene = VDJ_cgene[n == max(n)][1],
    #                 VDJ_SHM = VDJ_SHM[n == max(n)][1],
    #                 VDJ_cdr3s_aa = VDJ_cdr3s_aa[n == max(n)][1]
    #   ) %>%
    #   dplyr::select(-n, -barcode)
    # vgm.bulk.list[[1]] <- vgm.bulk.list[[1]] %>% dplyr::distinct()
    #
    # vgm.bulk.list[[1]] <- cbind(barcode = NA, vgm.bulk.list[[1]])
    # vgm.bulk.list[[1]] <- as.data.frame(vgm.bulk.list[[1]])


    ##### Creating vgm.bulk[[1]] (one line = one clonotype):

    most_frequent <- function(x) {
      if(is.numeric(x))
        return(as.numeric(names(which.max(table(x)))))
      else
        return(names(which.max(table(x))))
      }

    vgm.clonotyped[is.na(vgm.clonotyped)] <- ""

    columns_to_summarise <- c()
    for (n in names(vgm.clonotyped)) {
      if (length(unique(vgm.clonotyped[[n]]))>1 & !(n %in% c("clonotype_id","sample_id","barcode")))
        columns_to_summarise <- c(columns_to_summarise, n)}

    vgm.bulk.list[[1]] <- vgm.clonotyped %>%
      dplyr::select(-barcode) %>%
      dplyr::group_by(sample_id, clonotype_id) %>%
      # summarise(across(everything(), most_frequent))  ### TOO SLOW
      # summarise(across(columns_to_summarise, most_frequent), across(), .groups = "keep") ### Throws a warning
      # summarise(across(where( ~ length(unique(.x))>1 && ! .x %in% c("clonotype_id","sample_id","barcode")), most_frequent), across(), .groups = "keep") # warnings
      dplyr::summarise(dplyr::across(.cols = dplyr::all_of(columns_to_summarise) , most_frequent), dplyr::across(), .groups = "keep") %>%
      dplyr::distinct()

    vgm.bulk.list[[1]] <- cbind(barcode = NA, vgm.bulk.list[[1]])
    vgm.bulk.list[[1]] <- as.data.frame(vgm.bulk.list[[1]])

    vgm.bulk.list[[1]] <- vgm.bulk.list[[1]][order(as.numeric(gsub(".*?([0-9]+).*", "\\1", vgm.bulk.list[[1]]$sample_id)), as.numeric(gsub(".*?([0-9]+).*", "\\1", vgm.bulk.list[[1]]$clonotype_id))),]

    vgm.bulk.list[[1]][vgm.bulk.list[[1]] == ""] <- NA

    if (vgm.expanded == TRUE) {
      vgm.bulk.list[[9]] <- vgm.clonotyped
      vgm.bulk.list[[9]][vgm.bulk.list[[9]] == ""] <- NA
    }


  } else if (input.type == "TRUST4" | input.type == "TRUST4.FULL") {

    VDJ.dataframe.list <- list()
    data.df <- data.frame()
    for(i in 1:length(VDJ.bulk.out.directory.list)) {

      if (is.data.frame(VDJ.bulk.out.directory.list[[i]]))
        VDJ.dataframe.list[[i]] <- VDJ.bulk.out.directory.list[[i]]
      else
        message("ERROR, TRUST4 only supports .RDS as input.type!")
        # VDJ.dataframe.list[[i]] <-  utils::read.delim(VDJ.bulk.out.directory.list[[i]])
      VDJ.dataframe.list[[i]]$sample_id <- paste0("s", i)
      VDJ.dataframe.list[[i]]$group_id <- group.id[i]
      VDJ.dataframe.list[[i]]$batches <- batches[i]
      data.df <- rbind(data.df, VDJ.dataframe.list[[i]])
    }


    # initialize an empty data frame to store vgm[[1]]:
    vgm.bulk <- data.frame()

    i <- substr(data.df$c_call, 1, 3) %in% c("IGH", "TRB", "TRD","") # subset the heavy chains
    j <- substr(data.df$c_call, 1, 3) %in% c("IGL", "IGK", "TRA", "TRG") # subset the light chains
    # t <- substr(data.df$c_call, 1, 1) %in% c("T") # T cells


    # retrieving and appending information to the vgm.bulk data frame:

    #1
    barcode = rep(NA, nrow(data.df))

    #2
    sample_id <- data.df$sample_id

    #3 group_id column; group.id vector with user-defined integers specifying the group membership (defaults to NA):
    group_id <- data.df$group_id

    #4 clonotype_id_10x column:
    clonotype_id_10x <- rep(NA, nrow(data.df))

    #5
    celltype <- rep(cell.type, nrow(data.df))

    #6:
    Nr_of_VDJ_chains <- rep(NA, nrow(data.df))
    Nr_of_VDJ_chains[i] <- 1

    #7:
    Nr_of_VJ_chains <- rep(NA, nrow(data.df))
    Nr_of_VJ_chains[j]<- 1

    #8 VDJ_cdr3s_aa column:
    VDJ_cdr3s_aa <- rep(NA, nrow(data.df))
    VDJ_cdr3s_aa[i] <- data.df$junction_aa[i]

    #9:
    VJ_cdr3s_aa <- rep(NA, nrow(data.df))
    VJ_cdr3s_aa[j] <- data.df$junction_aa[j]

    #10:
    VDJ_cdr3s_nt <- rep(NA, nrow(data.df))
    VDJ_cdr3s_nt[i] <- data.df$junction[i]

    #11:
    VJ_cdr3s_nt <- rep(NA, nrow(data.df))
    VJ_cdr3s_nt[j] <- data.df$junction[j]

    #12:
    VDJ_chain_contig = rep(NA, nrow(data.df))

    #13:
    VJ_chain_contig = rep(NA, nrow(data.df))

    #14:
    VDJ_chain <- rep(NA, nrow(data.df))
    VDJ_chain[i] <- substr(data.df$c_call, 1, 3)[i]

    #15:
    VJ_chain <- rep(NA, nrow(data.df))
    VJ_chain[j] <- substr(data.df$c_call, 1, 3)[j]

    #16:
    VDJ_vgene <- rep(NA, nrow(data.df))
    VDJ_vgene[i] <- data.df$v_call[i]

    #17:
    VJ_vgene <- rep(NA, nrow(data.df))
    VJ_vgene[j] <- data.df$v_call[j]

    #18:
    VDJ_dgene <- rep(NA, nrow(data.df))
    VDJ_dgene[i] <- data.df$d_call[i]

    #19:
    VDJ_jgene <- rep(NA, nrow(data.df))
    VDJ_jgene[i] <- data.df$j_call[i]

    #20:
    VJ_jgene <- rep(NA, nrow(data.df))
    VDJ_jgene[j] <- data.df$j_call[j]

    #21:
    VDJ_cgene <- rep(NA, nrow(data.df))
    VDJ_cgene[i] <- data.df$c_call[i]

    #22:
    VJ_cgene <- rep(NA, nrow(data.df))
    VJ_cgene[j] <- data.df$c_call[j]

    #23:
    VDJ_sequence_nt_raw <- rep(NA, nrow(data.df))
    VDJ_sequence_nt_raw[i] <- data.df$sequence[i]

    #24:
    VJ_sequence_nt_raw <- rep(NA, nrow(data.df))
    VJ_sequence_nt_raw[j] <- data.df$sequence[j]

    #25:
    VDJ_sequence_nt_trimmed <- rep(NA, nrow(data.df))
    VDJ_sequence_nt_trimmed[i] <- data.df$sequence[i]

    #26:
    VJ_sequence_nt_trimmed <- rep(NA, nrow(data.df))
    VJ_sequence_nt_trimmed[j] <- data.df$sequence[j]

    #27:
    VDJ_sequence_aa <- rep(NA, nrow(data.df))

    #28:
    VJ_sequence_aa <- rep(NA, nrow(data.df))

    #29:
    VDJ_trimmed_ref <- rep(NA, nrow(data.df))

    #30:
    VJ_trimmed_ref <- rep(NA, nrow(data.df))

    #31
    VDJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #32
    VJ_raw_consensus_id <- rep(NA, nrow(data.df))

    #33
    orig_barcode <- rep(NA, nrow(data.df))

    #34
    clonotype_frequency <- rep(NA, nrow(data.df))

    #35
    specifity <- rep(NA, nrow(data.df))

    #36
    affinity <- rep(NA, nrow(data.df))

    #37
    GEX_available <- rep("FALSE", nrow(data.df))

    #38
    orig.ident <- rep(NA, nrow(data.df))

    #39
    orig_barcode_GEX <- rep(NA, nrow(data.df))

    #40
    seurat_clusters <- rep(NA, nrow(data.df))

    #41
    PC_1 <- rep(NA, nrow(data.df))

    #42
    PC_2 <- rep(NA, nrow(data.df))

    #43
    UMAP_1 <- rep(NA, nrow(data.df))

    #44
    UMAP_2 <- rep(NA, nrow(data.df))

    #45
    tSNE_1 <- rep(NA, nrow(data.df))

    #46
    tSNE_2 <- rep(NA, nrow(data.df))

    #47
    batches <- data.df$batches

    #48
    clonotype_id <- data.df$sequence_id   ### In TRUST4, sequence_id is not unique for each read. We use it as clonotype when no clonotyping strategy is specified as argument

    #49:
    VDJ_nSeqFR1 <- rep(NA, nrow(data.df))

    #50:
    VDJ_nSeqCDR1 <- rep(NA, nrow(data.df))

    #51:
    VDJ_nSeqFR2 <- rep(NA, nrow(data.df))

    #52:
    VDJ_nSeqCDR2 <- rep(NA, nrow(data.df))

    #53:
    VDJ_nSeqFR3 <- rep(NA, nrow(data.df))

    #54:
    VDJ_nSeqCDR3 <- rep(NA, nrow(data.df))

    #55:
    VDJ_nSeqFR4 <- rep(NA, nrow(data.df))

    #56:
    VDJ_aaSeqFR1 <- rep(NA, nrow(data.df))

    #57:
    VDJ_aaSeqCDR1 <- rep(NA, nrow(data.df))

    #58:
    VDJ_aaSeqFR2 <- rep(NA, nrow(data.df))

    #59:
    VDJ_aaSeqCDR2 <-rep(NA, nrow(data.df))

    #60:
    VDJ_aaSeqFR3 <- rep(NA, nrow(data.df))

    #61:
    VDJ_aaSeqCDR3 <- rep(NA, nrow(data.df))

    #62:
    VDJ_aaSeqFR4 <- rep(NA, nrow(data.df))

    #63:
    VDJ_SHM <- rep(NA, nrow(data.df))

    #64:
    VJ_nSeqFR1 <- rep(NA, nrow(data.df))

    #65:
    VJ_nSeqCDR1 <- rep(NA, nrow(data.df))

    #66:
    VJ_nSeqFR2 <- rep(NA, nrow(data.df))

    #67:
    VJ_nSeqCDR2 <- rep(NA, nrow(data.df))

    #68:
    VJ_nSeqFR3 <- rep(NA, nrow(data.df))

    #69:
    VJ_nSeqCDR3 <- rep(NA, nrow(data.df))

    #70:
    VJ_nSeqFR4 <- rep(NA, nrow(data.df))

    #71:
    VJ_aaSeqFR1 <- rep(NA, nrow(data.df))

    #72:
    VJ_aaSeqCDR1 <- rep(NA, nrow(data.df))

    #73:
    VJ_aaSeqFR2 <- rep(NA, nrow(data.df))

    #74:
    VJ_aaSeqCDR2 <- rep(NA, nrow(data.df))

    #75:
    VJ_aaSeqFR3 <- rep(NA, nrow(data.df))

    #76:
    VJ_aaSeqCDR3 <- rep(NA, nrow(data.df))

    #77:
    VJ_aaSeqFR4 <- rep(NA, nrow(data.df))

    #78:
    VJ_SHM <- rep(NA, nrow(data.df))

    ######## TRUST4 Additional columns

    cdr1 <- data.df$cdr1

    cdr2 <- data.df$cdr2

    v_cigar <- data.df$v_cigar

    d_cigar <- data.df$d_cigar

    j_cigar <- data.df$j_cigar

    v_identity <- data.df$v_identity

    j_identity <- data.df$j_identity

    complete_vdj <- data.df$complete_vdj

    consensus_count <- data.df$consensus_count

    # start of appending the required information to a local data frame. The first 48 columns that
    # are added do this data frame exactly match the formatting of vgm[[1]] generated after running
    # VDJ_GEX_matrix for single-cell data. In case of unavailable data, NA is used.

    vgm.bulk <- data.frame(barcode, sample_id, group_id, clonotype_id_10x,
                           celltype, Nr_of_VDJ_chains, Nr_of_VJ_chains,
                           VDJ_cdr3s_aa, VJ_cdr3s_aa, VDJ_cdr3s_nt,
                           VJ_cdr3s_nt, VDJ_chain_contig, VJ_chain_contig,
                           VDJ_chain, VJ_chain, VDJ_vgene, VJ_vgene,VDJ_dgene,
                           VDJ_jgene, VJ_jgene, VDJ_cgene, VJ_cgene,
                           VDJ_sequence_nt_raw, VJ_sequence_nt_raw,
                           VDJ_sequence_nt_trimmed, VJ_sequence_nt_trimmed, #
                           VDJ_sequence_aa, VJ_sequence_aa, VDJ_trimmed_ref,
                           VJ_trimmed_ref, VDJ_raw_consensus_id,
                           VJ_raw_consensus_id, orig_barcode,
                           clonotype_frequency, specifity, affinity,
                           GEX_available, orig.ident, orig_barcode_GEX,
                           seurat_clusters, PC_1, PC_2, UMAP_1, UMAP_2,
                           tSNE_1, tSNE_2, batches, clonotype_id,   ### 48
                           VDJ_nSeqFR1, VDJ_nSeqCDR1, VDJ_nSeqFR2,
                           VDJ_nSeqCDR2, VDJ_nSeqFR3, VJ_nSeqCDR3,
                           VDJ_nSeqFR4, VDJ_aaSeqFR1, VDJ_aaSeqCDR1,
                           VDJ_aaSeqFR2, VDJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VDJ_aaSeqCDR3, VDJ_aaSeqFR4, VDJ_SHM,
                           VJ_nSeqFR1, VJ_nSeqCDR1, VJ_nSeqFR2,
                           VJ_nSeqCDR2, VJ_nSeqFR3, VJ_nSeqCDR3,
                           VJ_nSeqFR4, VJ_aaSeqFR1, VJ_aaSeqCDR1,
                           VJ_aaSeqFR2, VJ_aaSeqCDR2, VJ_aaSeqFR3,
                           VJ_aaSeqCDR3, VJ_aaSeqFR4, VJ_SHM,   ### 78
                           ### Starting TRUST4 additional columns
                           cdr1, cdr2, v_cigar,
                           d_cigar, j_cigar, v_identity,
                           j_identity, complete_vdj, consensus_count) ### 87


    if (integrate.MIXCR.output == TRUE & input.type=="TRUST4") {
      vgm.bulk <- vgm.bulk[, 1:78]
    }
    else if (integrate.MIXCR.output == FALSE & input.type=="TRUST4") {
      vgm.bulk <- vgm.bulk[, 1:48]
    }
    else if (integrate.MIXCR.output == TRUE & input.type=="TRUST4.FULL") {
      vgm.bulk <- vgm.bulk[, 1:87]
    }
    else if (integrate.MIXCR.output == FALSE & input.type=="TRUST4.FULL") {
      vgm.bulk <- vgm.bulk[, c(1:48,79:87)]
    }

    # create an empty list vgm.bulk.list in we will store our two data frames:
    # vgm.bulk.list[[1]] contains aggregated clonotype-frequencies stored in a single clonotype_frequency row.
    # vgm.bulk.list[[9]] will contain "per-cell" information: rows will be added for each bulk clonotype based on
    # their aggregated clonotype_frequency number. I.e.: if we have clonotype0 with 8000 cells, vgm.bulk.list[[9]] will
    # contain 8000 rows for each cell of clonotype0. This is done to achieve better compatibility with some pre-existing Platypus
    # functions used for analyzing sc-VDJ sequence data.

    vgm.bulk.list <- list()

    if (want.clonotyping) {
      vgm.clonotyped <- VDJ_clonotype(VDJ = vgm.bulk, clone.strategy = clone.strategy,
                                         global.clonotype = F, VDJ.VJ.1chain = F, hierarchical = 'none') # TRY ALSO hierarchical = 'single.chains':
      # Not filtering cells with counts other than 1VDJ 1VJ chain and integrating these cells hierarchically into clonotypes
      }
    else {
      vgm.clonotyped <- vgm.bulk
    }


    most_frequent <- function(x) {
      if(is.numeric(x))
        return(as.numeric(names(which.max(table(x)))))
      else
        return(names(which.max(table(x))))
    }

    vgm.clonotyped[is.na(vgm.clonotyped)] <- ""

    columns_to_summarise <- c()
    for (n in names(vgm.clonotyped)) {
      if (length(unique(vgm.clonotyped[[n]]))>1 & !(n %in% c("clonotype_id","sample_id","barcode")))
        columns_to_summarise <- c(columns_to_summarise, n)}

    vgm.bulk.list[[1]] <- vgm.clonotyped %>%
      dplyr::select(-barcode) %>%
      dplyr::group_by(sample_id, clonotype_id) %>%
      dplyr::summarise(dplyr::across(.cols = dplyr::all_of(columns_to_summarise) , most_frequent), dplyr::across(), .groups = "keep") %>%
      dplyr::distinct()

    vgm.bulk.list[[1]] <- cbind(barcode = NA, vgm.bulk.list[[1]])
    vgm.bulk.list[[1]] <- as.data.frame(vgm.bulk.list[[1]])

    vgm.bulk.list[[1]] <- vgm.bulk.list[[1]][order(as.numeric(gsub(".*?([0-9]+).*", "\\1", vgm.bulk.list[[1]]$sample_id)), as.numeric(gsub(".*?([0-9]+).*", "\\1", vgm.bulk.list[[1]]$clonotype_id))),]

    vgm.bulk.list[[1]][vgm.bulk.list[[1]] == ""] <- NA

    if (vgm.expanded == TRUE) {
      vgm.bulk.list[[9]] <- vgm.clonotyped
      vgm.bulk.list[[9]][vgm.bulk.list[[9]] == ""] <- NA
      }

      ### TESTING #################################

      # vgm_bulk_cl12 <- vgm_bulk[[9]][vgm_bulk[[9]]$clonotype_id %in% c("clonotype1","clonotype2","clonotype3","clonotype4","clonotype5"),]
      # vgm_bulk_cl12[is.na(vgm_bulk_cl12)] <- ""
      #
      # done <- vgm_bulk_cl12 %>%
      #   dplyr::group_by(sample_id, clonotype_id) %>%
      #   summarise(across(everything(), mf))
      # done <- as.data.frame(done)
      #
      # for (si in unique(vgm_bulk_cl12$sample_id)){
      #   for (ci in unique(vgm_bulk_cl12[vgm_bulk_cl12$sample_id==si,]$clonotype_id)) {
      #
      #     # this_clone = vgm_bulk_cl12[vgm_bulk_cl12$sample_id==si & vgm_bulk_cl12$clonotype_id==ci,]
      #     print(names(which.max(table(vgm_bulk_cl12[vgm_bulk_cl12$sample_id==si & vgm_bulk_cl12$clonotype_id==ci, "VDJ_sequence_nt_trimmed"]))))
      #     print(names(which.max(table(vgm_bulk_cl12[vgm_bulk_cl12$sample_id==si & vgm_bulk_cl12$clonotype_id==ci, "VDJ_sequence_nt_trimmed"]))) == done[done$sample_id==si & done$clonotype_id==ci, "VDJ_sequence_nt_trimmed"])
      #   }
      # }

      #############################################

      # SAME AS ABOVE BUT WITH A FOR LOOP
      # for (si in unique(vgm_bulk[[1]]$sample_id)){
      #   for (ci in unique(vgm_bulk[[1]][vgm_bulk[[1]]$sample_id==si,]$clonotype_id)) {
      #
      #     this_clone = vgm_bulk[[1]][vgm_bulk[[1]]$sample_id==si & vgm_bulk[[1]]$clonotype_id==ci,]
      #
      #     for (co in names(vgm_bulk[[1]])) {
      #       a <- this_clone[[co]]
      #       if (!all(is.na(a)))
      #         vgm_bulk[[1]][vgm_bulk[[1]]$sample_id==si & vgm_bulk[[1]]$clonotype_id==ci, co] <- rep(names(which.max(na.omit(table(a)))), length(a))
      #       # else we keep all the NA
      #
      #       }
      #     }
      #   }

    }

  return(vgm.bulk.list)
}
