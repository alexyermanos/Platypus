#' Returns a list of clonotype dataframes following additional clonotyping. This function works best following filtering to ensure that each clone only has one heavy chain and one light chain.
#' @param VDJ For platypus v2 output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire. For platypus v3 VDJ output from the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]])
#' @param clone.strategy (Updated keywords, previous format is also functional) String describing the clonotyping strategy. Possible options include 'cdr3.nt', 'cdr3.aa','VDJJ.VJJ','VDJJ.VJJ.cdr3lengths','VDJJ.VJJ.cdr3length.cdr3homology', 'VDJJ.VJJ.cdr3length.VDJcdr3homology', 'cdr3.homology',or 'VDJcdr3.homology'. 'cdr3.aa' will convert the default cell ranger clonotyping to amino acid based. 'Hvj.Lvj' groups B cells with identical germline genes (V and J segments for both heavy chain and light chain. Those arguments including 'CDR3length' will group all sequences with identical CDRH3 and CDRL3 sequence lengths. Those arguments including 'CDR3homology' will additionally impose a homology requirement for CDRH3 and CDRL3 sequences.'CDR3homology',or 'CDRH3homology' will group sequences based on homology only (either of the whole CDR3 sequence or of the CDRH3 sequence respictevely).
#' All homology calculations are performed on the amino acid level.
#' @param homology.threshold Numeric value between 0 and 1 corresponding to the homology threshold forn the clone.strategy arguments that require a homology threshold. Default value is set to 70 percent sequence homology. For 70 percent homology, 0.3 should be supplied as input.
#' @param hierarchical Boolean. Defaults to FALSE. This is an extention specifically for cells with aberrant numbers of chains (i.e. 0VDJ 1VJ, 1VDJ 0VJ, 0VDJ 2VJ, 2VDJ 0VJ). Cells with 2VDJ 2VJ are filtered out as these are most likely doublets. Aberrant cells are clonotyped hierarchically in post, following this procedure: 1. define clonotypes classically with all cells containing exactly 1VDJ 1VJ chains. 2. For cells with only a single chain (either VDJ or VJ), check if any clone exists, which matches the clonotyping criteria for this chain. If true, add this cell to that clone. If false, create a new clone containing that cell. In case that more than 1 existing clone matches the aberrant cell, the cell is assigned to the most frequent existing clone. Two reasons are behind this decision: 2.1. The aberrant cells is numerically more likely to be a part of the more frequent existing clone. 2.2 In case of a wrong assignment, the effect of the error is lower, if an already expanded clone is increase by one count, rather than a existing non-expanded clone being assigned a second entry and thereby resulting as expanded. 3. For cells with 3 chains, verify the clonotyping criteria on both combinations of chains (i.e. VDJ1 - VJ1, VDJ2-VJ1 in case of a cell with 2VDJ 1VJ).
#' @param global.clonotype Logical specifying whether clonotyping should occur across samples or only within a single sample.
#' @param VDJ.VJ.1chain Logical specifying whether cells with multiple VDJ and VJ chains should be removed from the clonotyping. Can be either T or F for those definitions not requiring germline genes or homology thresholds, as calculating the later is difficult when multiple chains are present.
#' @param output.format String specifies function output format. Options are "vgm" (default), "dataframe.per.sample", "clone.level.dataframes", or "phylo.dataframe". "vgm" will update the existing $clonotype_id column of the input vgm, which is the output from VDJ_GEX_matrix. "dataframe.per.sample" will return a list of VDJ dataframes, where each dataframe contains the cell-level information for a given sample. "clone.level.dataframes" will convert the per.cell matrix to a clonal dataframe, in which cells of the same clone will be merged into a single row. "dataframe.per.clone" will generate nested lists of dataframes, where each dataframe contains cell-level information of a given clone.
#' @param platypus.version Default is "v3". To use the output of VDJ_GEX_matrix function, one should change this argument to "v3".
#' @return Returns a list of clonotype dataframes where each list element matches the  repertoire index in the input clonotype.list object. The dataframes will be updated with clonal frequencies based on the new clonotyping definition.
#' @export
#' @examples
#' reclonotyped_vgm <- VDJ_clonotype(VDJ=Platypus::small_vgm[[1]],
#'  clone.strategy="VDJJ.VJJ",
#'  homology.threshold=".3", platypus.version = "v3")
#'
VDJ_clonotype <- function(VDJ,
                          clone.strategy,
                          homology.threshold,
                          hierarchical,
                          VDJ.VJ.1chain,
                          global.clonotype,
                          output.format,
                          platypus.version){
  #for CRAN checks
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL
  ccombs1 <- NULL


  if(missing(platypus.version)) platypus.version <- "v3"
  if(missing(output.format)) output.format <- "vgm"
  if(missing(global.clonotype)) global.clonotype <- F
  if(missing(clone.strategy)) clone.strategy <- "cdr3.aa"
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
  if(missing(VDJ)) stop("Please provide input data as VDJ")
  if(missing(hierarchical)) hierarchical <- F

  #Making cloning stategy fitting with VDJ / VJ naming scheme
  #This way the old keyworks will still work and this update should not break any old code
  #remember for renaming later
  clone.strategy.as.input <- clone.strategy
  switch(clone.strategy,
         VDJJ.VJJ.cdr3length.cdr3homology = {clone.strategy <- 'hvj.lvj.CDR3length.CDR3homology'},
         VDJJ.VJJ.cdr3length.VDJcdr3homology  = {clone.strategy <- 'hvj.lvj.CDR3length.CDRH3homology'},
         cdr3.homology = {clone.strategy <- 'CDR3.homology'},
         VDJcdr3.homology = {clone.strategy <- 'CDRH3.homology'},
         VDJJ.VJJ = {clone.strategy <- "hvj.lvj"},
         VDJJ.VJJ.cdr3 = {clone.strategy <- "hvj.lvj.cdr3"},
         VDJJ.VJJ.cdr3lengths = {clone.strategy <- "hvj.lvj.cdr3lengths"})

  if(platypus.version=="v2"){####START v2

    #compatibility with input naming scheme
    clonotype.list <- VDJ
    VDJ <- NULL

    output.clonotype <- list()
    if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) message("No homology threshold supplied. Clonotyping based on 70% amino acid homology.")
    if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) homology.threshold<-0.3  # Setting default homology threshold

    #Possible strategy options:'cdr3.aa','hvj.lvj','hvj.lvj.cdr3lengths','hvj.lvj.cdr3length.cdr3homology', 'hvj.lvj.CDR3length.CDRH3homology', 'CDR3homology',or 'CDRH3homology'.
    for(i in 1:length(clonotype.list)){
      if(clone.strategy=="cdr3.nt"){
        unique_clones <- unique(clonotype.list[[i]]$CDR3_nt_pasted)
        clonotype.list[[i]]$new_clone_unique <- clonotype.list[[i]]$CDR3_nt_pasted
      }
      else if(clone.strategy=="cdr3.aa"){
        unique_clones <- unique(clonotype.list[[i]]$CDR3_aa_pasted)
        clonotype.list[[i]]$new_clone_unique <- clonotype.list[[i]]$CDR3_aa_pasted
      }
      else if(clone.strategy=="hvj.lvj"){
        unique_clones <- unique(paste(clonotype.list[[i]]$HC_vgene,
                                      clonotype.list[[i]]$HC_jgene,
                                      clonotype.list[[i]]$LC_vgene,
                                      clonotype.list[[i]]$LC_jgene,sep="_"))
        clonotype.list[[i]]$new_clone_unique <- paste(clonotype.list[[i]]$HC_vgene,
                                                      clonotype.list[[i]]$HC_jgene,
                                                      clonotype.list[[i]]$LC_vgene,
                                                      clonotype.list[[i]]$LC_jgene,sep="_")
      }
      else if(clone.strategy=="hvj.lvj.cdr3lengths"){
        unique_clones <- unique(paste(clonotype.list[[i]]$HC_vgene,
                                      clonotype.list[[i]]$HC_jgene,
                                      clonotype.list[[i]]$LC_vgene,
                                      clonotype.list[[i]]$LC_jgene,
                                      nchar(clonotype.list[[i]]$CDRH3_aa),
                                      nchar(clonotype.list[[i]]$CDRL3_aa),sep="_"))
        clonotype.list[[i]]$new_clone_unique <- paste(clonotype.list[[i]]$HC_vgene,
                                                      clonotype.list[[i]]$HC_jgene,
                                                      clonotype.list[[i]]$LC_vgene,
                                                      clonotype.list[[i]]$LC_jgene,
                                                      nchar(clonotype.list[[i]]$CDRH3_aa),
                                                      nchar(clonotype.list[[i]]$CDRL3_aa),sep="_")

      }
      else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology" | clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){  #taking into account both cases
        clones_temp <- (paste(clonotype.list[[i]]$HC_vgene,
                              clonotype.list[[i]]$HC_jgene,
                              clonotype.list[[i]]$LC_vgene,
                              clonotype.list[[i]]$LC_jgene,
                              nchar(clonotype.list[[i]]$CDRH3_aa),
                              nchar(clonotype.list[[i]]$CDRL3_aa),sep="_"))
        clonotype.list[[i]]$new_clone_unique <- clones_temp
        unique_clones <- unique(clones_temp)
        for(j in 1:length(unique_clones)){
          original_clone_indices <- which(clones_temp==unique_clones[j])
          ### calculate distance for all within each
          if (length(original_clone_indices) >= 2){
            #different vl_distance depending on the strategy
            vh_distance <- stringdist::stringdistmatrix(clonotype.list[[i]]$CDRH3_aa[original_clone_indices],clonotype.list[[i]]$CDRH3_aa[original_clone_indices],method = "lv")/nchar(clonotype.list[[i]]$CDRH3_aa[original_clone_indices])
            if (clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){
              vl_distance <- stringdist::stringdistmatrix(clonotype.list[[i]]$CDRL3_aa[original_clone_indices],clonotype.list[[i]]$CDRL3_aa[original_clone_indices],method = "lv")/nchar(clonotype.list[[i]]$CDRL3_aa[original_clone_indices])
            }else{
              vl_distance <- 0
            }
            combined_distance <- vh_distance + vl_distance
            diag(combined_distance) <- NA
            hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
            hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
            # paste j and cluster
            clonotype.list[[i]]$new_clone_unique[original_clone_indices] <- paste(clonotype.list[[i]]$new_clone_unique[original_clone_indices],j,hclust_combined_cut)
            # need to account for the fact that hclust will not work if we have only 1 object to cluster. So assign value manually to the groups with one object.
          }else{
            clonotype.list[[i]]$new_clone_unique[original_clone_indices] <- paste(clonotype.list[[i]]$new_clone_unique[original_clone_indices],j,"1")
          }
        }
        unique_clones <- unique(clonotype.list[[i]]$new_clone_unique)
      }
      else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){
        vh_distance <- stringdist::stringdistmatrix(clonotype.list[[i]]$CDRH3_aa, clonotype.list[[i]]$CDRH3_aa, method = "lv")/nchar(clonotype.list[[i]]$CDRH3_aa)
        if(clone.strategy=="CDR3.homology"){
          vl_distance <- stringdist::stringdistmatrix(clonotype.list[[i]]$CDRL3_aa, clonotype.list[[i]]$CDRL3_aa, method = "lv")/nchar(clonotype.list[[i]]$CDRL3_aa)
        }else{
          vl_distance <- 0
        }
        combined_distance <- vh_distance + vl_distance
        diag(combined_distance) <- NA
        hclust_combined <- stats::hclust(stats::as.dist(combined_distance))
        hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
        # paste j and cluster
        clonotype.list[[i]]$new_clone_unique <- paste(hclust_combined_cut)
        unique_clones <- unique(clonotype.list[[i]]$new_clone_unique)
      }
      clone_number <-length(unique_clones)
      output.clonotype[[i]] <- data.frame(clonotype_id=paste("clonotype",1:clone_number,sep=""),frequency=rep(NA,clone_number),proportion=rep("",clone_number),cdr3s_aa=rep("",clone_number),cdr3s_nt=rep("",clone_number),HC_count=rep("",clone_number),IGK_count=rep("",clone_number),IGL_count=rep("",clone_number),LC_count=rep("",clone_number),CDRH3_aa=rep("",clone_number),CDRL3_aa=rep("",clone_number),CDRH3_nt=rep("",clone_number),CDRL3_nt=rep("",clone_number),CDR3_aa_pasted=rep("",clone_number),CDR3_nt_pasted=rep("",clone_number),HC_cgene=rep("",clone_number),HC_vgene=rep("",clone_number),HC_dgene=rep("",clone_number),HC_jgene=rep("",clone_number),LC_cgene=rep("",clone_number),LC_vgene=rep("",clone_number),LC_jgene=rep("",clone_number),barcodes=rep("",clone_number),nt_clone_ids=rep("",clone_number),new_unique_clone=unique_clones,nt_clone_cdrh3s=rep("",clone_number),nt_clone_cdrl3s=rep("",clone_number),stringsAsFactors = F)

      for(j in 1:length(unique_clones)){
        output.clonotype[[i]]$frequency[j] <- sum(clonotype.list[[i]]$frequency[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])
        output.clonotype[[i]]$proportion[j] <- output.clonotype[[i]]$frequency[j]/sum(clonotype.list[[i]]$frequency)
        output.clonotype[[i]]$cdr3s_aa[j] <- names(which.max(table(clonotype.list[[i]]$cdr3s_aa[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$cdr3s_nt[j] <- names(which.max(table(clonotype.list[[i]]$cdr3s_nt[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$HC_count[j] <- names(which.max(table(clonotype.list[[i]]$HC_count[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$IGK_count[j] <- names(which.max(table(clonotype.list[[i]]$IGK_count[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$IGL_count[j] <- names(which.max(table(clonotype.list[[i]]$IGL_count[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$LC_count[j] <- names(which.max(table(clonotype.list[[i]]$LC_count[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))

        output.clonotype[[i]]$CDRH3_aa[j] <- names(which.max(table(clonotype.list[[i]]$CDRH3_aa[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$CDRL3_aa[j] <- names(which.max(table(clonotype.list[[i]]$CDRL3_aa[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$CDRH3_nt[j] <- names(which.max(table(clonotype.list[[i]]$CDRH3_nt[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$CDRL3_nt[j] <- names(which.max(table(clonotype.list[[i]]$CDRL3_nt[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))

        output.clonotype[[i]]$CDR3_aa_pasted[j] <- names(which.max(table(clonotype.list[[i]]$CDR3_aa_pasted[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$CDR3_nt_pasted[j] <- names(which.max(table(clonotype.list[[i]]$CDR3_nt_pasted[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$HC_cgene[j] <- names(which.max(table(clonotype.list[[i]]$HC_cgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$HC_vgene[j] <- names(which.max(table(clonotype.list[[i]]$HC_vgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$HC_dgene[j] <- names(which.max(table(clonotype.list[[i]]$HC_dgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$HC_jgene[j] <- names(which.max(table(clonotype.list[[i]]$HC_jgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$LC_cgene[j] <- names(which.max(table(clonotype.list[[i]]$LC_cgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$LC_vgene[j] <- names(which.max(table(clonotype.list[[i]]$LC_vgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$LC_jgene[j] <- names(which.max(table(clonotype.list[[i]]$LC_jgene[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])])))
        output.clonotype[[i]]$barcodes[j] <- gsub(toString(clonotype.list[[i]]$barcodes[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])]),pattern=", ", replacement = ";")
        output.clonotype[[i]]$nt_clone_ids[j] <- gsub(toString(clonotype.list[[i]]$clonotype_id[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])]),pattern = ", ", replacement = ";")
        output.clonotype[[i]]$nt_clone_cdrh3s[j] <- gsub(toString(clonotype.list[[i]]$CDRH3_nt[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])]),pattern = ", ",replacement = ";")
        output.clonotype[[i]]$nt_clone_cdrl3s[j] <- gsub(toString(clonotype.list[[i]]$CDRL3_aa[which(clonotype.list[[i]]$new_clone_unique==unique_clones[j])],sep=";"),pattern = ", ",replacement = ";")
      }
    }
    return(output.clonotype)
  }####STOP v2
  if(platypus.version=="v3"){####START v3

    #compatibility with input naming scheme
    VDJ.GEX.matrix <- list()
    VDJ.GEX.matrix[[1]] <- VDJ
    VDJ <- NULL

    if(hierarchical == F){ #Standard clonotyping for all cells including those with abberant chain numbers


      if(global.clonotype==F){ # loop through each repertoire individually
        repertoire.number <- unique(VDJ.GEX.matrix[[1]]$sample_id)
        sample_dfs <- list()
        for(i in 1:length(repertoire.number)){ ####START sample loop
          sample_dfs[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==repertoire.number[i]),]
          ####only include clones with one heavy and one light chain
          if(VDJ.VJ.1chain== T){
            sample_dfs[[i]]<- sample_dfs[[i]][which(sample_dfs[[i]]$Nr_of_VDJ_chains==1 & sample_dfs[[i]]$Nr_of_VJ_chains==1),]
          }####STOP strict

          ####Clonotyping strategies
          if(clone.strategy=="10x.default"){ ####START cdr3.nt

            sample_dfs[[i]]$new_clonal_feature <- sample_dfs[[i]]$clonotype_id_10x
          } ####STOP cdr3.nt
          if(clone.strategy=="cdr3.nt"){ ####START cdr3.nt

            sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_nt,
                                                         sample_dfs[[i]]$VJ_cdr3s_nt)
          } ####STOP cdr3.nt
          else if(clone.strategy=="cdr3.aa"){ ####START cdr3.aa
            sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                         sample_dfs[[i]]$VJ_cdr3s_aa)
          } ####STOP cdr3.aa
          else if(clone.strategy=="hvj.lvj"){ ####START hvj.lvj
            sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        sample_dfs[[i]]$VJ_jgene,
                                                        sample_dfs[[i]]$VJ_jgene,sep="_")
          }####STOP hvj.lvj
          else if(clone.strategy=="hvj.lvj.cdr3"){ ####START hvj.lvj.cdr3
            sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                        sample_dfs[[i]]$VJ_cdr3s_aa,
                                                        sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        sample_dfs[[i]]$VJ_jgene,
                                                        sample_dfs[[i]]$VJ_jgene,sep="_")
          }####STOP hvj.lvj
          else if(clone.strategy=="hvj.lvj.cdr3lengths"){ ####START hvj.lvj.cdr3lengths
            sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        sample_dfs[[i]]$VJ_vgene,
                                                        sample_dfs[[i]]$VJ_jgene,
                                                        nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                                        nchar(sample_dfs[[i]]$VJ_cdr3s_aa),sep="_")

          } ####STOP hvj.lvj.cdr3lengths   / START Homology based clonotyping
          else if(clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology" | clone.strategy=="Hvj.Lvj.CDR3length.CDRH3homology"){  #taking into account both cases
            clones_temp <- (paste(sample_dfs[[i]]$VDJ_vgene,
                                  sample_dfs[[i]]$VDJ_jgene,
                                  sample_dfs[[i]]$VJ_vgene,
                                  sample_dfs[[i]]$VJ_jgene,
                                  nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                  nchar(sample_dfs[[i]]$VJ_cdr3s_aa),sep="_"))
            sample_dfs[[i]]$new_clonal_feature <- clones_temp
            unique_clones <- unique(clones_temp)
            for(j in 1:length(unique_clones)){
              original_clone_indices <- which(clones_temp==unique_clones[j])
              ### calculate distance for all within each
              if(length(original_clone_indices) >= 2){
                #different vl_distance depending on the strategy

                #Deal with the possibility of a missing chain and verify that a nchar length is present
                if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                  #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                  nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
                  nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
                } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                  nchars_vh <- rep(1,length(original_clone_indices))
                } else {
                  nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
                }

                vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vh
                if (clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology"){

                  if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                    #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                    nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                    nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
                  } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                    nchars_vl <- rep(1,length(original_clone_indices))
                  } else {
                    nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                  }

                  vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vl
                }else{
                  vl_distance <- 0
                }
                combined_distance <- vh_distance + vl_distance
                diag(combined_distance) <- NA
                hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
                hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
                # paste j and cluster
                sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,hclust_combined_cut)
                # need to account for the fact that hclust will not work if we have only 1 object to cluster. So assign value manually to the groups with one object.
              }else{
                sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,"1")
              }
            }
            unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)
          }
          else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){

            #Deal with the possibility of a missing chain and verify that a nchar length is present
            if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
              #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
              nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
              nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
            } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
              nchars_vh <- rep(1,length(sample_dfs[[i]]$VDJ_cdr3s_aa))
            } else {
              nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
            }

            vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VDJ_cdr3s_aa, method = "lv")/nchars_vh
            if(clone.strategy=="CDR3.homology"){

              #Deal with the possibility of a missing chain and verify that a nchar length is present
              if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
                #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
                nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
                nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
              } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
                nchars_vl <- rep(1,length(sample_dfs[[i]]$VJ_cdr3s_aa))
              } else {
                nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
              }

              vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa, method = "lv")/nchars_vl
            }else{
              vl_distance <- 0
            }
            combined_distance <- vh_distance + vl_distance
            diag(combined_distance) <- NA
            hclust_combined <- stats::hclust(stats::as.dist(combined_distance))
            hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
            # paste j and cluster
            sample_dfs[[i]]$new_clonal_feature <- paste(hclust_combined_cut)
            unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)
          }

          ####START recalculating clonotype_id and clonal_frequency
          #place holders
          sample_dfs[[i]]$new_clonotype_id <- rep(NA,nrow(sample_dfs[[i]]))
          sample_dfs[[i]]$new_clonal_frequency <- rep(NA,nrow(sample_dfs[[i]]))
          sample_dfs[[i]]$new_clonal_rank <- rep(NA,nrow(sample_dfs[[i]]))


          unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
          unique.clonal.frequencies <- rep(NA,length(unique.clonal.features))
          for(j in 1:length(unique.clonal.features)){####START assigning new frequency
            unique.clonal.frequencies[j] <- length(which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j]))
            sample_dfs[[i]]$new_clonal_frequency[which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j])] <- unique.clonal.frequencies[j]
          }####STOP assigning new frequency

          #assigning new new_clonal_rank
          sample_dfs[[i]] <-sample_dfs[[i]][with(sample_dfs[[i]], order(-new_clonal_frequency)), ]
          unique.clone.frequencies <- unique(sample_dfs[[i]]$new_clonal_frequency)
          for(j in 1:length(unique.clone.frequencies)){
            sample_dfs[[i]]$new_clonal_rank[which(sample_dfs[[i]]$new_clonal_frequency==unique.clone.frequencies[j])] <- j
          }####STOP assigning new new_clonal_rank

          #new clonotype id
          unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
          for(j in 1:length(unique.clonal.features)){
            sample_dfs[[i]]$new_clonotype_id[which(sample_dfs[[i]]$new_clonal_feature == unique.clonal.features[j])] <- paste0("clonotype",j)
          }####STOP assigning new clonotype id
        }####STOP sample loop
      }####STOP global.clonotype==F
      else if(global.clonotype==T){####START global.clonotype == T
        sample_dfs <- VDJ.GEX.matrix[[1]]
        sample_dfs$clonotype_id_10x <- paste0(sample_dfs$clonotype_id_10x,"_",sample_dfs$sample_id)

        if(VDJ.VJ.1chain==T){
          sample_dfs <- sample_dfs[which(sample_dfs$Nr_of_VDJ_chains==1 & sample_dfs$Nr_of_VJ_chains==1), ]}

        if(clone.strategy=="10x.default"){ ####START cdr3.nt
          sample_dfs$new_clonal_feature <- sample_dfs$clonotype_id_10x
        } ####STOP cdr3.nt
        if(clone.strategy=="cdr3.nt"){ ####START cdr3.nt
          sample_dfs$new_clonal_feature <- paste0(sample_dfs$VDJ_cdr3s_nt,
                                                  sample_dfs$VJ_cdr3s_nt)
        } ####STOP cdr3.nt
        else if(clone.strategy=="cdr3.aa"){ ####START cdr3.aa
          sample_dfs$new_clonal_feature <- paste0(sample_dfs$VDJ_cdr3s_aa,
                                                  sample_dfs$VJ_cdr3s_aa)
        } ####STOP cdr3.aa
        else if(clone.strategy=="hvj.lvj"){ ####START hvj.lvj
          sample_dfs$new_clonal_feature <- paste(sample_dfs$VDJ_vgene,
                                                 sample_dfs$VDJ_jgene,
                                                 sample_dfs$VJ_vgene,
                                                 sample_dfs$VJ_jgene,sep="_")
        }####STOP hvj.lvj
        else if(clone.strategy=="hvj.lvj.cdr3"){ ####START hvj.lvj.cdr3
          sample_dfs$new_clonal_feature <- paste(sample_dfs$VDJ_cdr3s_aa,
                                                 sample_dfs$VJ_cdr3s_aa,
                                                 sample_dfs$VDJ_vgene,
                                                 sample_dfs$VDJ_jgene,
                                                 sample_dfs$VJ_jgene,
                                                 sample_dfs$VJ_jgene,sep="_")
        }####STOP hvj.lvj.cdr3
        else if(clone.strategy=="hvj.lvj.cdr3lengths"){ ####START hvj.lvj.cdr3lengths
          sample_dfs$new_clonal_feature <- paste(sample_dfs$VDJ_vgene,
                                                 sample_dfs$VDJ_jgene,
                                                 sample_dfs$VJ_vgene,
                                                 sample_dfs$VJ_jgene,
                                                 nchar(sample_dfs$VDJ_cdr3s_aa),
                                                 nchar(sample_dfs$VJ_cdr3s_aa),sep="_")
        }####STOP hvj.lvj.cdr3lengths
        else if(clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology" | clone.strategy=="Hvj.Lvj.CDR3length.CDRH3homology"){  #taking into account both cases
          clones_temp <- (paste(sample_dfs$VDJ_vgene,
                                sample_dfs$VDJ_jgene,
                                sample_dfs$VJ_vgene,
                                sample_dfs$VJ_jgene,
                                nchar(sample_dfs$VDJ_cdr3s_aa),
                                nchar(sample_dfs$VJ_cdr3s_aa),sep="_"))
          sample_dfs$new_clonal_feature <- clones_temp
          unique_clones <- unique(clones_temp)
          for(j in 1:length(unique_clones)){
            original_clone_indices <- which(clones_temp==unique_clones[j])
            ### calculate distance for all within each
            if(length(original_clone_indices) >= 2){
              #different vl_distance depending on the strategy

              #Deal with the possibility of a missing chain and verify that a nchar length is present
              if(any(nchar(sample_dfs$VDJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                nchars_vh <- nchar(sample_dfs$VDJ_cdr3s_aa[original_clone_indices])
                nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
              } else if(all(nchar(sample_dfs$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                nchars_vh <- rep(1,length(original_clone_indices))
              } else {
                nchars_vh <- nchar(sample_dfs$VDJ_cdr3s_aa[original_clone_indices])
              }

              vh_distance <- stringdist::stringdistmatrix(sample_dfs$VDJ_cdr3s_aa[original_clone_indices],sample_dfs$VDJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vh
              if (clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology"){

                if(any(nchar(sample_dfs$VJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                  #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                  nchars_vl <- nchar(sample_dfs$VJ_cdr3s_aa[original_clone_indices])
                  nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
                } else if(all(nchar(sample_dfs$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                  nchars_vl <- rep(1,length(original_clone_indices))
                } else {
                  nchars_vl <- nchar(sample_dfs$VJ_cdr3s_aa[original_clone_indices])
                }

                vl_distance <- stringdist::stringdistmatrix(sample_dfs$VJ_cdr3s_aa[original_clone_indices],sample_dfs$VJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vl
              }else{
                vl_distance <- 0
              }
              combined_distance <- vh_distance + vl_distance
              diag(combined_distance) <- NA
              hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
              hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
              # paste j and cluster
              sample_dfs$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,hclust_combined_cut)
              # need to account for the fact that hclust will not work if we have only 1 object to cluster. So assign value manually to the groups with one object.
            }else{
              sample_dfs$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,"1")
            }
          }
          unique_clones <- unique(sample_dfs$new_clonal_feature)
        }
        else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){

          #Deal with the possibility of a missing chain and verify that a nchar length is present
          if(any(nchar(sample_dfs$VDJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs$VDJ_cdr3s_aa) == 0)){
            #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
            nchars_vh <- nchar(sample_dfs$VDJ_cdr3s_aa)
            nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
          } else if(all(nchar(sample_dfs$VDJ_cdr3s_aa) == 0)){
            nchars_vh <- rep(1,length(sample_dfs$VDJ_cdr3s_aa))
          } else {
            nchars_vh <- nchar(sample_dfs$VDJ_cdr3s_aa)
          }

          vh_distance <- stringdist::stringdistmatrix(sample_dfs$VDJ_cdr3s_aa, sample_dfs$VDJ_cdr3s_aa, method = "lv")/nchars_vh
          if(clone.strategy=="CDR3.homology"){

            #Deal with the possibility of a missing chain and verify that a nchar length is present
            if(any(nchar(sample_dfs$VJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs$VJ_cdr3s_aa) == 0)){
              #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
              nchars_vl <- nchar(sample_dfs$VJ_cdr3s_aa)
              nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
            } else if(all(nchar(sample_dfs$VJ_cdr3s_aa) == 0)){
              nchars_vl <- rep(1,length(sample_dfs$VJ_cdr3s_aa))
            } else {
              nchars_vl <- nchar(sample_dfs$VJ_cdr3s_aa)
            }

            vl_distance <- stringdist::stringdistmatrix(sample_dfs$VJ_cdr3s_aa, sample_dfs$VJ_cdr3s_aa, method = "lv")/nchars_vl
          }else{
            vl_distance <- 0
          }
          combined_distance <- vh_distance + vl_distance
          diag(combined_distance) <- NA
          hclust_combined <- stats::hclust(stats::as.dist(combined_distance))
          hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
          # paste j and cluster
          sample_dfs$new_clonal_feature <- paste(hclust_combined_cut)
          unique_clones <- unique(sample_dfs$new_clonal_feature)
        }

        ####START recalculating clonotype_id and clonal_frequency
        #definde placeholder columns
        sample_dfs$new_clonotype_id <- rep(NA,nrow(sample_dfs))
        sample_dfs$new_clonal_frequency <- rep(NA,nrow(sample_dfs))
        sample_dfs$new_clonal_rank <- rep(NA,nrow(sample_dfs))

        unique.clonal.features <- unique(sample_dfs$new_clonal_feature)

        unique.clonal.frequencies <- rep(NA,length(unique.clonal.features))
        for(j in 1:length(unique.clonal.features)){####START assigning new frequency
          unique.clonal.frequencies[j] <- length(which(sample_dfs$new_clonal_feature==unique.clonal.features[j]))
          sample_dfs$new_clonal_frequency[which(sample_dfs$new_clonal_feature==unique.clonal.features[j])] <- unique.clonal.frequencies[j]
        }####STOP assigning new frequency

        #assigning new new_clonal_rank
        sample_dfs <-sample_dfs[with(sample_dfs, order(-new_clonal_frequency)), ]
        unique.clone.frequencies <- unique(sample_dfs$new_clonal_frequency)
        for(j in 1:length(unique.clone.frequencies)){
          sample_dfs$new_clonal_rank[which(sample_dfs$new_clonal_frequency==unique.clone.frequencies[j])] <- j
        }###

        #new clonotype id
        unique.clonal.features <- unique(sample_dfs$new_clonal_feature)
        for(j in 1:length(unique.clonal.features)){
          sample_dfs$new_clonotype_id[which(sample_dfs$new_clonal_feature == unique.clonal.features[j])] <- paste0("clonotype",j)
        }
      }####STOP global.clonotype==T

    } ####stop hierarchical == F

    if(hierarchical){ #Clonotyping were cells with abberant chain numbers are added sequentially to existing clones to improve their integration
      #Preamble to this non elegant solution: All ifs and elses below are coded out explicitely for each case. A helper function with variable inputs would be nicer. For now this should ensure that the code is easily understandable.
      if(global.clonotype==F){ # loop through each repertoire individually

        repertoire.number <- unique(VDJ.GEX.matrix[[1]]$sample_id)
        sample_dfs <- list()
        for(i in 1:length(repertoire.number)){ ####START sample loop
          sample_dfs[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==repertoire.number[i]),]

          prior_filtering <- nrow(sample_dfs[[i]])
          sample_dfs[[i]] <- subset(sample_dfs[[i]],  (Nr_of_VDJ_chains > 0 | Nr_of_VJ_chains > 0) & sample_dfs[[i]]$Nr_of_VDJ_chains + sample_dfs[[i]]$Nr_of_VJ_chains < 4)
          if(nrow(sample_dfs[[i]]) > 0){
          message(paste0("Filtered out ", prior_filtering - nrow(sample_dfs[[i]]), " cells containing more than one VDJ AND VJ chain, as these likely correspond to doublets"))}

          ####only include clones with one heavy and one light chain
          if(VDJ.VJ.1chain== T){message("Hierarchical clonotyping is specifically designed to incorporate cells with abberand numbers of chains. Filtering for 1VDJ 1VJ chain thereby defeats its purpose. Function will continue without filtering.")}####STOP strict

          #Prepwork to increase function speed
          #find cells with 1VJ chain only
          aberant_cells <- subset(sample_dfs[[i]], Nr_of_VDJ_chains != 1 | Nr_of_VJ_chains != 1)
          onlyVJ_ind <- which(aberant_cells$Nr_of_VJ_chains > 0 & aberant_cells$Nr_of_VDJ_chains == 0)
          onlyVDJ_ind <- which(aberant_cells$Nr_of_VDJ_chains > 0 & aberant_cells$Nr_of_VJ_chains == 0)
          #cells with more than one VJ chain and one VDJ chain
          multVJ_ind <- which(aberant_cells$Nr_of_VJ_chains > 1 & aberant_cells$Nr_of_VDJ_chains == 1)
          #cells with more than one VDJ chain and one VJ chain
          multVDJ_ind <- which(aberant_cells$Nr_of_VDJ_chains > 1 & aberant_cells$Nr_of_VJ_chains == 1)
          #add the new column already
          aberant_cells$new_clonal_feature <- NA

          #now filter out the rest of abberant cells from the sample_dfs
          sample_dfs[[i]] <- subset(sample_dfs[[i]],Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)

          ####Clonotyping strategies
          if(clone.strategy=="10x.default"){ ####START cdr3.nt
            sample_dfs[[i]]$new_clonal_feature <- sample_dfs[[i]]$clonotype_id_10x
          } ####STOP cdr3.nt
          if(clone.strategy=="cdr3.nt"){ ####START cdr3.nt

            sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_nt,
                                                         sample_dfs[[i]]$VJ_cdr3s_nt)

            n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_cdr3s_nt[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,2])))
              } else {
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_nt[cel]))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_nt[cel]
              }
            }
            }
            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_cdr3s_nt[cel], ";")){ #catches a cell with two vDJ chains but no VJ chain (very very rare)
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,2])))
                } else {
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VDJ_cdr3s_nt[cel]))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- aberant_cells$VDJ_cdr3s_nt[cel]
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- aberant_cells$VDJ_cdr3s_nt[cel]
                VJs <- stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,]

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2])

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_nt[cel],aberant_cells$VJ_cdr3s_nt[cel])
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,]
                VJs <- aberant_cells$VJ_cdr3s_nt[cel]

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2])

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_nt[cel],aberant_cells$VJ_cdr3s_nt[cel])
                }
              }
            }


          } ####STOP cdr3.nt
          else if(clone.strategy=="cdr3.aa"){ ####START cdr3.aa
            sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                         sample_dfs[[i]]$VJ_cdr3s_aa)


            n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VJ_cdr3s_aa[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
                } else {
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_aa[cel]))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_aa[cel]
                }
              }
            }


            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_cdr3s_aa[cel], ";")){ #catches a cell with two vDJ chains but no VJ chain (very very rare)
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
                } else {
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VDJ_cdr3s_aa[cel]))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- aberant_cells$VDJ_cdr3s_aa[cel]
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- aberant_cells$VDJ_cdr3s_aa[cel]
                VJs <- stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2])

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel],aberant_cells$VJ_cdr3s_aa[cel])
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]
                VJs <- aberant_cells$VJ_cdr3s_aa[cel]

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2])

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel],aberant_cells$VJ_cdr3s_aa[cel])
                }
              }
            }

          } ####STOP cdr3.aa
          else if(clone.strategy=="hvj.lvj"){ ####START hvj.lvj
            sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        sample_dfs[[i]]$VJ_vgene,
                                                        sample_dfs[[i]]$VJ_jgene,sep="_")


            n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }


            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- paste0(aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                VJs <- c(paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- c(paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
                VJs <- paste0(aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }


          }####STOP hvj.lvj
          else if(clone.strategy=="hvj.lvj.cdr3"){ ####START hvj.lvj.cdr3
            sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                        sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        sample_dfs[[i]]$VJ_cdr3s_aa,
                                                        sample_dfs[[i]]$VJ_vgene,
                                                        sample_dfs[[i]]$VJ_jgene,sep="_")

            n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(aberant_cells$VJ_cdr3s_aa[cel], "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }


            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                VJs <- c(paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- c(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
                VJs <- paste0(aberant_cells$VJ_cdr3s_aa[cel], "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }


          }####STOP hvj.lvj
          else if(clone.strategy=="hvj.lvj.cdr3lengths"){ ####START hvj.lvj.cdr3lengths
            sample_dfs[[i]]$new_clonal_feature <- paste(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                                        sample_dfs[[i]]$VDJ_vgene,
                                                        sample_dfs[[i]]$VDJ_jgene,
                                                        nchar(sample_dfs[[i]]$VJ_cdr3s_aa),
                                                        sample_dfs[[i]]$VJ_vgene,
                                                        sample_dfs[[i]]$VJ_jgene, sep="_")

            n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }


            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_", aberant_cells$VDJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- comb1
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                VJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
                VJs <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- ccombs[1]
                }
              }
            }

          } ####STOP hvj.lvj.cdr3lengths   / START Homology based clonotyping
          else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology" | clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){  #taking into account both cases

            clones_temp <- (paste(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                  sample_dfs[[i]]$VDJ_vgene,
                                  sample_dfs[[i]]$VDJ_jgene,
                                  nchar(sample_dfs[[i]]$VJ_cdr3s_aa),
                                  sample_dfs[[i]]$VJ_vgene,
                                  sample_dfs[[i]]$VJ_jgene,sep="_"))
            sample_dfs[[i]]$new_clonal_feature <- clones_temp
            unique_clones <- unique(clones_temp)
            for(j in 1:length(unique_clones)){
              original_clone_indices <- which(clones_temp==unique_clones[j])
              ### calculate distance for all within each
              if(length(original_clone_indices) >= 2){
                #different vl_distance depending on the strategy

                #Deal with the possibility of a missing chain and verify that a nchar length is present
                if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                  #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                  nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
                  nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
                } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                  nchars_vh <- rep(1,length(original_clone_indices))
                } else {
                  nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
                }

                vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vh
                if (clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){

                  if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                    #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                    nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                    nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
                  } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                    nchars_vl <- rep(1,length(original_clone_indices))
                  } else {
                    nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                  }

                  vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vl
                }else{
                  vl_distance <- 0
                }
                combined_distance <- vh_distance + vl_distance
                diag(combined_distance) <- NA
                hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
                hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
                # paste j and cluster
                sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,hclust_combined_cut)
                # need to account for the fact that hclust will not work if we have only 1 object to cluster. So assign value manually to the groups with one object.
              }else{
                sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,"1")
              }
            }
            unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)


            #Now deal with cells of abberant clone numbers
            #Generate a surrogate clonal feature including classical features + the new homology based one
            sur_clonal_feature <- paste(clones_temp, "_", sample_dfs[[i]]$new_clonal_feature)

            #Proceed very similarly to above during length and V J gene clonotyping with the added step of a stringdist
            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])

                  clone_matches <- c(which(stringr::str_detect(sur_clonal_feature, comb1)), which(stringr::str_detect(sur_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sur_clonal_feature, comb1))
                }
                #We now have one or a set of matching clones / Now we check their homology via stringdist

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                  if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #in the case that homology is only calculated on the VDJ chain, here we cannot check homology, because the aberrant query cell does not contain a VDJ chain
                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone
                  } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we check for homology with all matching clones for the VJ chain
                    if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){
                      dists1 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                      dists2 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                      dists <- c(dists1,dists2) #concatenate
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                      }
                    } else {
                      dists <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], aberant_cells$VJ_cdr3s_aa[cel]) / nchar(aberant_cells$VJ_cdr3s_aa[cel]) #get distances
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                      }
                    }
                  }
                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                  #We still check if the CDR3 homology threshold is fullfilled here
                  if(any(stringdist::stringdist(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                    aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                  } else { #open a new clone
                    aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                  }
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                }
              }
            }


            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                #check if the light chain matches any already existing clone
                #check if the light chain matches any already existing clone
                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                  comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
                } else {
                  comb1 <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_", aberant_cells$VDJ_jgene[cel])
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
                }

                #We now have one or a set of matching clones / Now we check their homology via stringdist

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                  if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology

                    if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                      dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                      dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                      dists <- c(dists1,dists2) #concatenate
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                      }
                    } else {
                      dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                      }
                    }

                  } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we now also only check for VDJ homology, because there is no VJ chain

                    if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                      dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                      dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                      dists <- c(dists1,dists2) #concatenate
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                      }
                    } else {
                      dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VJD_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                      }
                    }
                  }
                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                  #We still check if the CDR3 homology threshold is fullfilled here
                  if(any(stringdist::stringdist(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                    aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                  } else { #open a new clone
                    aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                  }
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){
                #get combinations
                VDJs <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                VJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology / given that we know that this cell has 1 VDJ chain, we can skip a few ifs

                      dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                      }

                  } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we have to make sure that homology is fullfilled for both the VDJ chain and one of the VJ chains

                      dists1 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                       paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])) / nchar(paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]))
                      dists2 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                       paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])) / nchar(paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]))
                      dists <- c(dists1,dists2) #concatenate
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                      }
                    }
                  } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the aberrant query clone
                  if(stringdist::stringdist(aberant_cells$VDJ_cdr3s_aa[cel],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(aberant_cells$VDJ_cdr3s_aa[cel]) <= homology.threshold
                     & any(stringdist::stringdist(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                    aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                  } else { #open a new clone
                    aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                  }

                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){
                #get combinations
                VDJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
                VJs <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                #check if any combination of VDJ and VJ chains matches any already existing clone
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                  if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology / given that we know that this cell has 2 VDJ chain, we str_split those

                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    dists <- c(dists1,dists2) #concatenate
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel],"_", aberant_cells$VJ_cdr3s_aa[cel])
                    }

                  } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we have to make sure that homology is fullfilled for both the VDJ chain and one of the VJ chains

                    dists1 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                     paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_", aberant_cells$VJ_cdr3s_aa[cel])) / nchar(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_", aberant_cells$VJ_cdr3s_aa[cel]))
                    dists2 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                     paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_", aberant_cells$VJ_cdr3s_aa[cel])) / nchar(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_", aberant_cells$VJ_cdr3s_aa[cel]))
                    dists <- c(dists1,dists2) #concatenate
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                    }
                  }
                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the aberrant query clone
                  if(stringdist::stringdist(aberant_cells$VJ_cdr3s_aa[cel],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(aberant_cells$VJ_cdr3s_aa[cel]) <= homology.threshold
                    & any(stringdist::stringdist(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                      aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                    } else { #open a new clone
                      aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                    }

                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                }
              }
            }

          }
          else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){

            #Deal with the possibility of a missing chain and verify that a nchar length is present
            if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
              #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
              nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
              nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
            } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
              nchars_vh <- rep(1,length(sample_dfs[[i]]$VDJ_cdr3s_aa))
            } else {
              nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
            }

            vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VDJ_cdr3s_aa, method = "lv")/nchars_vh
            if(clone.strategy=="CDR3.homology"){

              #Deal with the possibility of a missing chain and verify that a nchar length is present
              if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
                #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
                nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
                nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
              } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
                nchars_vl <- rep(1,length(sample_dfs[[i]]$VJ_cdr3s_aa))
              } else {
                nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
              }

              vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa, method = "lv")/nchars_vl
            }else{
              vl_distance <- 0
            }
            combined_distance <- vh_distance + vl_distance
            diag(combined_distance) <- NA
            hclust_combined <- stats::hclust(stats::as.dist(combined_distance))
            hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
            # paste j and cluster
            sample_dfs[[i]]$new_clonal_feature <- paste(hclust_combined_cut)
            unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)


            #Just check for stringdist
            #check cells with only one VJ chain and nothing else
            if(length(onlyVJ_ind) > 0){
              for(cel in onlyVJ_ind){

                  if(clone.strategy=="CDRH3.homology"){ #in the case that homology is only calculated on the VDJ chain, here we cannot check homology, so we check if the light chain matches to an existing clone exaclty

                    if(stringr::str_detect(aberant_cells$VJ_cdr3s_aa[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                      clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
                    } else {
                      clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_aa[cel]))
                    }

                    if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                      aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone
                    } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                      aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                    } else { #no clone found with the light chain of this cell => open a new clone
                      aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_aa[cel]
                    }

                  } else if(clone.strategy=="CDR3homology"){ #here we check for homology with all matching clones for the VJ chain
                    if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){
                      dists1 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                      dists2 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                      #get the minimum of both for each element
                      dists <- c()
                      for(k in 1:length(dists1)){
                        if(dists1[k] < dists2[k]){
                          dists <- c(dists, dists1[k])
                        } else{
                          dists <- c(dists, dists2[k])
                        }
                      }

                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VJ_cdr3s_aa[cel])
                      }
                    } else {
                      dists <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, aberant_cells$VJ_cdr3s_aa[cel]) / nchar(aberant_cells$VJ_cdr3s_aa[cel]) #get distances
                      if(any(dists <= homology.threshold)){
                        aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                      } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                        aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VJ_cdr3s_aa[cel])
                      }
                    }
                  }
                }
              }

            #check cells with only one VDJ chain and nothing else
            if(length(onlyVDJ_ind) > 0){
              for(cel in onlyVDJ_ind){
                if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                  if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    #get the minimum of both for each element
                    dists <- c()
                    for(k in 1:length(dists1)){
                      if(dists1[k] < dists2[k]){
                        dists <- c(dists, dists1[k])
                      } else{
                        dists <- c(dists, dists2[k])
                      }
                    }

                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  } else {
                    dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  }


                } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                  if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    #get the minimum of both for each element
                    dists <- c()
                    for(k in 1:length(dists1)){
                      if(dists1[k] < dists2[k]){
                        dists <- c(dists, dists1[k])
                      } else{
                        dists <- c(dists, dists2[k])
                      }
                    }

                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  } else {
                    dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  }
                }
              }
            }

            #check cells with 2 VJ chains and 1 VDJ chain
            if(length(multVJ_ind) > 0){
              for(cel in multVJ_ind){

                if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                    dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }

                } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                  #get combinations
                  VDJs <- aberant_cells$VDJ_cdr3s_aa[cel]
                  VJs <- c(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])

                  ccombs <- expand.grid(VDJs, VJs)
                  ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                  pasted_sample_dfs <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa)

                  dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[1]) / nchar(ccombs[1])
                  dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[2]) / nchar(ccombs[2])
                    #get the minimum of both for each element
                    dists <- c()
                    for(k in 1:length(dists1)){
                      if(dists1[k] < dists2[k]){
                        dists <- c(dists, dists1[k])
                      } else{
                        dists <- c(dists, dists2[k])
                      }
                    }

                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- ccombs1
                    }
                }
              }
            }

            #check cells with 2 VDJ chains and 1 VJ chain
            if(length(multVDJ_ind) > 0){
              for(cel in multVDJ_ind){

                if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    #get the minimum of both for each element
                    dists <- c()
                    for(k in 1:length(dists1)){
                      if(dists1[k] < dists2[k]){
                        dists <- c(dists, dists1[k])
                      } else{
                        dists <- c(dists, dists2[k])
                      }
                    }

                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                    }

                } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                  #get combinations
                  VDJs <- c(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                  VJs <- aberant_cells$VJ_cdr3s_aa[cel]

                  ccombs <- expand.grid(VDJs, VJs)
                  ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                  pasted_sample_dfs <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa)

                  dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[1]) / nchar(ccombs[1])
                  dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[2]) / nchar(ccombs[2])
                  #get the minimum of both for each element
                  dists <- c()
                  for(k in 1:length(dists1)){
                    if(dists1[k] < dists2[k]){
                      dists <- c(dists, dists1[k])
                    } else{
                      dists <- c(dists, dists2[k])
                    }
                  }

                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- ccombs1
                  }
                }
              }

            }
          }

          #### Combine the aberant cells dataframe with the rest
          sample_dfs[[i]] <- rbind(sample_dfs[[i]], aberant_cells)
          unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)

          ####START recalculating clonotype_id and clonal_frequency
          #place holders
          sample_dfs[[i]]$new_clonotype_id <- rep(NA,nrow(sample_dfs[[i]]))
          sample_dfs[[i]]$new_clonal_frequency <- rep(NA,nrow(sample_dfs[[i]]))
          sample_dfs[[i]]$new_clonal_rank <- rep(NA,nrow(sample_dfs[[i]]))


          unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
          unique.clonal.frequencies <- rep(NA,length(unique.clonal.features))

          for(j in 1:length(unique.clonal.features)){####START assigning new frequency
            unique.clonal.frequencies[j] <- length(which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j]))
            sample_dfs[[i]]$new_clonal_frequency[which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j])] <- unique.clonal.frequencies[j]
          }####STOP assigning new frequency

          #assigning new new_clonal_rank
          sample_dfs[[i]] <-sample_dfs[[i]][with(sample_dfs[[i]], order(-new_clonal_frequency)), ]
          unique.clone.frequencies <- unique(sample_dfs[[i]]$new_clonal_frequency)
          for(j in 1:length(unique.clone.frequencies)){
            sample_dfs[[i]]$new_clonal_rank[which(sample_dfs[[i]]$new_clonal_frequency==unique.clone.frequencies[j])] <- j
          }####STOP assigning new new_clonal_rank

          #new clonotype id
          unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
          for(j in 1:length(unique.clonal.features)){
            sample_dfs[[i]]$new_clonotype_id[which(sample_dfs[[i]]$new_clonal_feature == unique.clonal.features[j])] <- paste0("clonotype",j)
          }####STOP assigning new clonotype id
        }####STOP sample loop
      }####STOP global.clonotype==F

      else if(global.clonotype==T){####START global.clonotype == T#
        sample_dfs <- list() #to avoid recoding. This is a bodge
        sample_dfs[[1]] <- VDJ.GEX.matrix[[1]]
        i <- 1
        sample_dfs$clonotype_id_10x <- paste0(sample_dfs$clonotype_id_10x,"_",sample_dfs$sample_id)

        prior_filtering <- nrow(sample_dfs[[i]])
        sample_dfs[[i]] <- subset(sample_dfs[[i]],  (Nr_of_VDJ_chains > 0 | Nr_of_VJ_chains > 0) & sample_dfs[[i]]$Nr_of_VDJ_chains + sample_dfs[[i]]$Nr_of_VJ_chains < 4)
        if(nrow(sample_dfs[[i]]) > 0){
          message(paste0("Filtered out ", prior_filtering - nrow(sample_dfs[[i]]), " cells containing more than one VDJ AND VJ chain, as these likely correspond to doublets"))}

        ####only include clones with one heavy and one light chain
        if(VDJ.VJ.1chain== T){message("Hierarchical clonotyping is specifically designed to better incorporate cells with abberand numbers of chains. Filtering for 1VDJ 1VJ chain thereby defeats its purpose. Function will continue with out filtering. For standard clonotyping with filtering set hierarchical = FALSE. ")}####STOP strict

        #Prepwork to increase function speed
        #find cells with 1VJ chain only
        aberant_cells <- subset(sample_dfs[[i]], Nr_of_VDJ_chains != 1 | Nr_of_VJ_chains != 1)
        onlyVJ_ind <- which(aberant_cells$Nr_of_VJ_chains > 0 & aberant_cells$Nr_of_VDJ_chains == 0)
        onlyVDJ_ind <- which(aberant_cells$Nr_of_VDJ_chains > 0 & aberant_cells$Nr_of_VJ_chains == 0)
        #cells with more than one VJ chain and one VDJ chain
        multVJ_ind <- which(aberant_cells$Nr_of_VJ_chains > 1 & aberant_cells$Nr_of_VDJ_chains == 1)
        #cells with more than one VDJ chain and one VJ chain
        multVDJ_ind <- which(aberant_cells$Nr_of_VDJ_chains > 1 & aberant_cells$Nr_of_VJ_chains == 1)
        #add the new column already
        aberant_cells$new_clonal_feature <- NA

        #now filter out the rest of abberant cells from the sample_dfs
        sample_dfs[[i]] <- subset(sample_dfs[[i]],Nr_of_VDJ_chains == 1 & Nr_of_VJ_chains == 1)

        if(clone.strategy=="10x.default"){ ####START cdr3.nt
          sample_dfs[[i]]$new_clonal_feature <- sample_dfs[[i]]$clonotype_id_10x
        } ####STOP cdr3.nt
        if(clone.strategy=="cdr3.nt"){ ####START cdr3.nt

          sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_nt,
                                                       sample_dfs[[i]]$VJ_cdr3s_nt)

          n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_cdr3s_nt[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,2])))
              } else {
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_nt[cel]))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_nt[cel]
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_cdr3s_nt[cel], ";")){ #catches a cell with two vDJ chains but no VJ chain (very very rare)
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,2])))
              } else {
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VDJ_cdr3s_nt[cel]))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- aberant_cells$VDJ_cdr3s_nt[cel]
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- aberant_cells$VDJ_cdr3s_nt[cel]
              VJs <- stringr::str_split(aberant_cells$VJ_cdr3s_nt[cel], ";", simplify = T)[1,]

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2])

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_nt[cel],aberant_cells$VJ_cdr3s_nt[cel])
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- stringr::str_split(aberant_cells$VDJ_cdr3s_nt[cel], ";", simplify = T)[1,]
              VJs <- aberant_cells$VJ_cdr3s_nt[cel]

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2])

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_nt[cel],aberant_cells$VJ_cdr3s_nt[cel])
              }
            }
          }


        } ####STOP cdr3.nt
        else if(clone.strategy=="cdr3.aa"){ ####START cdr3.aa
          sample_dfs[[i]]$new_clonal_feature <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                       sample_dfs[[i]]$VJ_cdr3s_aa)


          n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_cdr3s_aa[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
              } else {
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_aa[cel]))

              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone


              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone

                aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_aa[cel]
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_cdr3s_aa[cel], ";")){ #catches a cell with two vDJ chains but no VJ chain (very very rare)
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
              } else {
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VDJ_cdr3s_aa[cel]))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone


                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone


                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- aberant_cells$VDJ_cdr3s_aa[cel]

              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- aberant_cells$VDJ_cdr3s_aa[cel]
              VJs <- stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2])

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel],aberant_cells$VJ_cdr3s_aa[cel])
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]
              VJs <- aberant_cells$VJ_cdr3s_aa[cel]

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2])

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel],aberant_cells$VJ_cdr3s_aa[cel])
              }
            }
          }

        } ####STOP cdr3.aa
        else if(clone.strategy=="hvj.lvj"){ ####START hvj.lvj
          sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                      sample_dfs[[i]]$VDJ_jgene,
                                                      sample_dfs[[i]]$VJ_vgene,
                                                      sample_dfs[[i]]$VJ_jgene,sep="_")


          n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- paste0(aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
              VJs <- c(paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- c(paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
              VJs <- paste0(aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }


        }####STOP hvj.lvj
        else if(clone.strategy=="hvj.lvj.cdr3"){ ####START hvj.lvj.cdr3
          sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                      sample_dfs[[i]]$VDJ_vgene,
                                                      sample_dfs[[i]]$VDJ_jgene,
                                                      sample_dfs[[i]]$VJ_cdr3s_aa,
                                                      sample_dfs[[i]]$VJ_vgene,
                                                      sample_dfs[[i]]$VJ_jgene,sep="_")

          n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(aberant_cells$VJ_cdr3s_aa[cel], "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
              VJs <- c(paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- c(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
              VJs <- paste0(aberant_cells$VJ_cdr3s_aa[cel], "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }


        }####STOP hvj.lvj
        else if(clone.strategy=="hvj.lvj.cdr3lengths"){ ####START hvj.lvj.cdr3lengths
          sample_dfs[[i]]$new_clonal_feature <- paste(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                                      sample_dfs[[i]]$VDJ_vgene,
                                                      sample_dfs[[i]]$VDJ_jgene,
                                                      nchar(sample_dfs[[i]]$VJ_cdr3s_aa),
                                                      sample_dfs[[i]]$VJ_vgene,
                                                      sample_dfs[[i]]$VJ_jgene, sep="_")

          n_new_clones <- length(unique(sample_dfs[[i]]$new_clonal_feature))

          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_", aberant_cells$VDJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- comb1
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
              VJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
              VJs <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone

              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone

                aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- ccombs[1]
              }
            }
          }

        } ####STOP hvj.lvj.cdr3lengths   / START Homology based clonotyping
        else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology" | clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){  #taking into account both cases

          clones_temp <- (paste(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa),
                                sample_dfs[[i]]$VDJ_vgene,
                                sample_dfs[[i]]$VDJ_jgene,
                                nchar(sample_dfs[[i]]$VJ_cdr3s_aa),
                                sample_dfs[[i]]$VJ_vgene,
                                sample_dfs[[i]]$VJ_jgene,sep="_"))
          sample_dfs[[i]]$new_clonal_feature <- clones_temp
          unique_clones <- unique(clones_temp)
          for(j in 1:length(unique_clones)){
            original_clone_indices <- which(clones_temp==unique_clones[j])
            ### calculate distance for all within each
            if(length(original_clone_indices) >= 2){
              #different vl_distance depending on the strategy

              #Deal with the possibility of a missing chain and verify that a nchar length is present
              if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
                nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
              } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices]) == 0)){
                nchars_vh <- rep(1,length(original_clone_indices))
              } else {
                nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices])
              }

              vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VDJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vh
              if (clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){

                if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                  #get nchar of heavy chains and surrogate missing once with the mean of existing ones
                  nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                  nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
                } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices]) == 0)){
                  nchars_vl <- rep(1,length(original_clone_indices))
                } else {
                  nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices])
                }

                vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],sample_dfs[[i]]$VJ_cdr3s_aa[original_clone_indices],method = "lv")/nchars_vl
              }else{
                vl_distance <- 0
              }
              combined_distance <- vh_distance + vl_distance
              diag(combined_distance) <- NA
              hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
              hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
              # paste j and cluster
              sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,hclust_combined_cut)
              # need to account for the fact that hclust will not work if we have only 1 object to cluster. So assign value manually to the groups with one object.
            }else{
              sample_dfs[[i]]$new_clonal_feature[original_clone_indices] <- paste(sample_dfs[[i]]$new_clonal_feature[original_clone_indices],j,"1")
            }
          }
          unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)


          #Now deal with cells of abberant clone numbers
          #Generate a surrogate clonal feature including classical features + the new homology based one
          sur_clonal_feature <- paste(clones_temp, "_", sample_dfs[[i]]$new_clonal_feature)

          #Proceed very similarly to above during length and V J gene clonotyping with the added step of a stringdist
          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2])

                clone_matches <- c(which(stringr::str_detect(sur_clonal_feature, comb1)), which(stringr::str_detect(sur_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_",aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sur_clonal_feature, comb1))
              }
              #We now have one or a set of matching clones / Now we check their homology via stringdist

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #in the case that homology is only calculated on the VDJ chain, here we cannot check homology, because the aberrant query cell does not contain a VDJ chain
                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone
                } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we check for homology with all matching clones for the VJ chain
                  if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){
                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    dists <- c(dists1,dists2) #concatenate
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                    }
                  } else {
                    dists <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches], aberant_cells$VJ_cdr3s_aa[cel]) / nchar(aberant_cells$VJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                    }
                  }
                }
              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                #We still check if the CDR3 homology threshold is fullfilled here
                if(any(stringdist::stringdist(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                } else { #open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
                }
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VJ_cdr3s_aa[cel])
              }
            }
          }


          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              #check if the light chain matches any already existing clone
              #check if the light chain matches any already existing clone
              if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                comb1 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1])
                comb2 <- paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_" ,stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2])
                clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1)), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb2)))
              } else {
                comb1 <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_", aberant_cells$VDJ_jgene[cel])
                clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, comb1))
              }

              #We now have one or a set of matching clones / Now we check their homology via stringdist

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology

                  if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    dists <- c(dists1,dists2) #concatenate
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  } else {
                    dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  }

                } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we now also only check for VDJ homology, because there is no VJ chain

                  if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                    dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                    dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                    dists <- c(dists1,dists2) #concatenate
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  } else {
                    dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VJD_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                    if(any(dists <= homology.threshold)){
                      aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                    } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                      aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                    }
                  }
                }
              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                #We still check if the CDR3 homology threshold is fullfilled here
                if(any(stringdist::stringdist(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                } else { #open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                }
              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){
              #get combinations
              VDJs <- paste0(nchar(aberant_cells$VDJ_cdr3s_aa[cel]), "_",aberant_cells$VDJ_vgene[cel], "_",aberant_cells$VDJ_jgene[cel])
              VJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VJ_jgene[cel], ";", simplify = T)[1,2]))

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology / given that we know that this cell has 1 VDJ chain, we can skip a few ifs

                  dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[clone_matches[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel])
                  }

                } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we have to make sure that homology is fullfilled for both the VDJ chain and one of the VJ chains

                  dists1 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                   paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])) / nchar(paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]))
                  dists2 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                   paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])) / nchar(paste0(aberant_cells$VDJ_cdr3s_aa[cel], "_", stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]))
                  dists <- c(dists1,dists2) #concatenate
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                  }
                }
              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the aberrant query clone
                if(stringdist::stringdist(aberant_cells$VDJ_cdr3s_aa[cel],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(aberant_cells$VDJ_cdr3s_aa[cel]) <= homology.threshold
                   & any(stringdist::stringdist(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                } else { #open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                }

              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){
              #get combinations
              VDJs <- c(paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,1], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,1]),paste0(nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]), "_",stringr::str_split(aberant_cells$VDJ_vgene[cel], ";", simplify = T)[1,2], "_", stringr::str_split(aberant_cells$VDJ_jgene[cel], ";", simplify = T)[1,2]))
              VJs <- paste0(nchar(aberant_cells$VJ_cdr3s_aa[cel]), "_", aberant_cells$VJ_vgene[cel], "_",aberant_cells$VJ_jgene[cel])

              ccombs <- expand.grid(VDJs, VJs)
              ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

              #check if any combination of VDJ and VJ chains matches any already existing clone
              clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[1])),which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature,ccombs[2])))

              if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone
                if(clone.strategy=="hvj.lvj.CDR3length.CDRH3homology"){ #checking homology / given that we know that this cell has 2 VDJ chain, we str_split those

                  dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                  dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                  dists <- c(dists1,dists2) #concatenate
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(comb1,"_", aberant_cells$VDJ_cdr3s_aa[cel],"_", aberant_cells$VJ_cdr3s_aa[cel])
                  }

                } else if(clone.strategy=="hvj.lvj.cdr3length.cdr3homology"){ #here we have to make sure that homology is fullfilled for both the VDJ chain and one of the VJ chains

                  dists1 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                   paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_", aberant_cells$VJ_cdr3s_aa[cel])) / nchar(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], "_", aberant_cells$VJ_cdr3s_aa[cel]))
                  dists2 <- stringdist::stringdist(paste0(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches],"_",sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]),
                                                   paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_", aberant_cells$VJ_cdr3s_aa[cel])) / nchar(paste0(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2], "_", aberant_cells$VJ_cdr3s_aa[cel]))
                  dists <- c(dists1,dists2) #concatenate
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[rep(clone_matches,2)[which.min(dists)]] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                  }
                }
              } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the aberrant query clone
                if(stringdist::stringdist(aberant_cells$VJ_cdr3s_aa[cel],unique(sample_dfs[[i]]$VJ_cdr3s_aa[clone_matches]))/nchar(aberant_cells$VJ_cdr3s_aa[cel]) <= homology.threshold
                   & any(stringdist::stringdist(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,],unique(sample_dfs[[i]]$VDJ_cdr3s_aa[clone_matches]))/nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,]) <= homology.threshold)){ #make sure to account for the rare case of 2 VJ chains and 0 VDJ chains
                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #assign the clone
                } else { #open a new clone
                  aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
                }

              } else { #no clone found with the light chain of this cell => open a new clone
                aberant_cells$new_clonal_feature[cel] <- paste0(ccombs[1],"_", aberant_cells$VDJ_cdr3s_aa[cel], "_", aberant_cells$VJ_cdr3s_aa[cel])
              }
            }
          }

        }
        else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){

          #Deal with the possibility of a missing chain and verify that a nchar length is present
          if(any(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
            #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
            nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
            nchars_vh[which(nchars_vh == 0)] <- mean(nchars_vh[nchars_vh > 0])
          } else if(all(nchar(sample_dfs[[i]]$VDJ_cdr3s_aa) == 0)){
            nchars_vh <- rep(1,length(sample_dfs[[i]]$VDJ_cdr3s_aa))
          } else {
            nchars_vh <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
          }

          vh_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VDJ_cdr3s_aa, method = "lv")/nchars_vh
          if(clone.strategy=="CDR3.homology"){

            #Deal with the possibility of a missing chain and verify that a nchar length is present
            if(any(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0) & !all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
              #get nchar number of heavy chains and surrogate missing once with the mean of existing ones
              nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
              nchars_vl[which(nchars_vl == 0)] <- mean(nchars_vl[nchars_vl > 0])
            } else if(all(nchar(sample_dfs[[i]]$VJ_cdr3s_aa) == 0)){
              nchars_vl <- rep(1,length(sample_dfs[[i]]$VJ_cdr3s_aa))
            } else {
              nchars_vl <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
            }

            vl_distance <- stringdist::stringdistmatrix(sample_dfs[[i]]$VJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa, method = "lv")/nchars_vl
          }else{
            vl_distance <- 0
          }
          combined_distance <- vh_distance + vl_distance
          diag(combined_distance) <- NA
          hclust_combined <- stats::hclust(stats::as.dist(combined_distance))
          hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
          # paste j and cluster
          sample_dfs[[i]]$new_clonal_feature <- paste(hclust_combined_cut)
          unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)


          #Just check for stringdist
          #check cells with only one VJ chain and nothing else
          if(length(onlyVJ_ind) > 0){
            for(cel in onlyVJ_ind){

              if(clone.strategy=="CDRH3.homology"){ #in the case that homology is only calculated on the VDJ chain, here we cannot check homology, so we check if the light chain matches to an existing clone exaclty

                if(stringr::str_detect(aberant_cells$VJ_cdr3s_aa[cel], ";")){ #catches a cell with two vJ chains but no VDJ chain (very rare)
                  clone_matches <- c(which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])), which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])))
                } else {
                  clone_matches <- which(stringr::str_detect(sample_dfs[[i]]$new_clonal_feature, aberant_cells$VJ_cdr3s_aa[cel]))
                }

                if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

                  aberant_cells$new_clonal_feature[cel] <- names(which.max(table(sample_dfs[[i]]$new_clonal_feature[clone_matches]))) #Assigning the aberrant query clone to the most frequent matching clone
                } else if(length(unique(sample_dfs[[i]]$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
                  aberant_cells$new_clonal_feature[cel] <- unique(sample_dfs[[i]]$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
                } else { #no clone found with the light chain of this cell => open a new clone
                  aberant_cells$new_clonal_feature[cel] <- aberant_cells$VJ_cdr3s_aa[cel]
                }

              } else if(clone.strategy=="CDR3homology"){ #here we check for homology with all matching clones for the VJ chain
                if(stringr::str_detect(aberant_cells$VJ_vgene[cel], ";")){
                  dists1 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                  dists2 <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                  #get the minimum of both for each element
                  dists <- c()
                  for(k in 1:length(dists1)){
                    if(dists1[k] < dists2[k]){
                      dists <- c(dists, dists1[k])
                    } else{
                      dists <- c(dists, dists2[k])
                    }
                  }

                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VJ_cdr3s_aa[cel])
                  }
                } else {
                  dists <- stringdist::stringdist(sample_dfs[[i]]$VJ_cdr3s_aa, aberant_cells$VJ_cdr3s_aa[cel]) / nchar(aberant_cells$VJ_cdr3s_aa[cel]) #get distances
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VJ_cdr3s_aa[cel])
                  }
                }
              }
            }
          }

          #check cells with only one VDJ chain and nothing else
          if(length(onlyVDJ_ind) > 0){
            for(cel in onlyVDJ_ind){
              if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                  dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                  dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                  #get the minimum of both for each element
                  dists <- c()
                  for(k in 1:length(dists1)){
                    if(dists1[k] < dists2[k]){
                      dists <- c(dists, dists1[k])
                    } else{
                      dists <- c(dists, dists2[k])
                    }
                  }

                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                  }
                } else {
                  dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                  }
                }


              } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                if(stringr::str_detect(aberant_cells$VDJ_vgene[cel], ";")){
                  dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                  dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                  #get the minimum of both for each element
                  dists <- c()
                  for(k in 1:length(dists1)){
                    if(dists1[k] < dists2[k]){
                      dists <- c(dists, dists1[k])
                    } else{
                      dists <- c(dists, dists2[k])
                    }
                  }

                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                  }
                } else {
                  dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                  if(any(dists <= homology.threshold)){
                    aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                  } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                    aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                  }
                }
              }
            }
          }

          #check cells with 2 VJ chains and 1 VDJ chain
          if(length(multVJ_ind) > 0){
            for(cel in multVJ_ind){

              if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                dists <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, aberant_cells$VDJ_cdr3s_aa[cel]) / nchar(aberant_cells$VDJ_cdr3s_aa[cel]) #get distances
                if(any(dists <= homology.threshold)){
                  aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                }

              } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                #get combinations
                VDJs <- aberant_cells$VDJ_cdr3s_aa[cel]
                VJs <- c(stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,1], stringr::str_split(aberant_cells$VJ_cdr3s_aa[cel], ";", simplify = T)[1,2])

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                pasted_sample_dfs <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa)

                dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[1]) / nchar(ccombs[1])
                dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[2]) / nchar(ccombs[2])
                #get the minimum of both for each element
                dists <- c()
                for(k in 1:length(dists1)){
                  if(dists1[k] < dists2[k]){
                    dists <- c(dists, dists1[k])
                  } else{
                    dists <- c(dists, dists2[k])
                  }
                }

                if(any(dists <= homology.threshold)){
                  aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                  aberant_cells$new_clonal_feature[cel] <- ccombs1
                }
              }
            }
          }

          #check cells with 2 VDJ chains and 1 VJ chain
          if(length(multVDJ_ind) > 0){
            for(cel in multVDJ_ind){

              if(clone.strategy=="CDRH3.homology"){ #check homology for the VDJ chain

                dists1 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1])
                dists2 <- stringdist::stringdist(sample_dfs[[i]]$VDJ_cdr3s_aa, stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2]) / nchar(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                #get the minimum of both for each element
                dists <- c()
                for(k in 1:length(dists1)){
                  if(dists1[k] < dists2[k]){
                    dists <- c(dists, dists1[k])
                  } else{
                    dists <- c(dists, dists2[k])
                  }
                }

                if(any(dists <= homology.threshold)){
                  aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                  aberant_cells$new_clonal_feature[cel] <- paste0(aberant_cells$VDJ_cdr3s_aa[cel])
                }

              } else if(clone.strategy=="CDR3homology"){ #check homology for VDJ again, because there is no VJ chain

                #get combinations
                VDJs <- c(stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,1], stringr::str_split(aberant_cells$VDJ_cdr3s_aa[cel], ";", simplify = T)[1,2])
                VJs <- aberant_cells$VJ_cdr3s_aa[cel]

                ccombs <- expand.grid(VDJs, VJs)
                ccombs <- paste0(ccombs[,1], ccombs[,2]) #make it an array

                pasted_sample_dfs <- paste0(sample_dfs[[i]]$VDJ_cdr3s_aa, sample_dfs[[i]]$VJ_cdr3s_aa)

                dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[1]) / nchar(ccombs[1])
                dists1 <- stringdist::stringdist(pasted_sample_dfs, ccombs[2]) / nchar(ccombs[2])
                #get the minimum of both for each element
                dists <- c()
                for(k in 1:length(dists1)){
                  if(dists1[k] < dists2[k]){
                    dists <- c(dists, dists1[k])
                  } else{
                    dists <- c(dists, dists2[k])
                  }
                }

                if(any(dists <= homology.threshold)){
                  aberant_cells$new_clonal_feature[cel] <- sample_dfs[[i]]$new_clonal_feature[which.min(dists)] #Assigning the aberrant the closest clone (which.min only returns the index of the first minimal value if this value appears more than once)
                } else { #No pattern match found which also is within the homology threshold -> open a new clonotype
                  aberant_cells$new_clonal_feature[cel] <- ccombs1
                }
              }
            }

          }
        }

        #### Combine the aberant cells dataframe with the rest
        sample_dfs[[i]] <- rbind(sample_dfs[[i]], aberant_cells)
        unique_clones <- unique(sample_dfs[[i]]$new_clonal_feature)

        ####START recalculating clonotype_id and clonal_frequency
        #place holders
        sample_dfs[[i]]$new_clonotype_id <- rep(NA,nrow(sample_dfs[[i]]))
        sample_dfs[[i]]$new_clonal_frequency <- rep(NA,nrow(sample_dfs[[i]]))
        sample_dfs[[i]]$new_clonal_rank <- rep(NA,nrow(sample_dfs[[i]]))


        unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
        unique.clonal.frequencies <- rep(NA,length(unique.clonal.features))
        for(j in 1:length(unique.clonal.features)){####START assigning new frequency
          unique.clonal.frequencies[j] <- length(which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j]))
          sample_dfs[[i]]$new_clonal_frequency[which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j])] <- unique.clonal.frequencies[j]
        }####STOP assigning new frequency

        #assigning new new_clonal_rank
        sample_dfs[[i]] <-sample_dfs[[i]][with(sample_dfs[[i]], order(-new_clonal_frequency)), ]
        unique.clone.frequencies <- unique(sample_dfs[[i]]$new_clonal_frequency)
        for(j in 1:length(unique.clone.frequencies)){
          sample_dfs[[i]]$new_clonal_rank[which(sample_dfs[[i]]$new_clonal_frequency==unique.clone.frequencies[j])] <- j
        }####STOP assigning new new_clonal_rank

        #new clonotype id
        unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
        for(j in 1:length(unique.clonal.features)){
          sample_dfs[[i]]$new_clonotype_id[which(sample_dfs[[i]]$new_clonal_feature == unique.clonal.features[j])] <- paste0("clonotype",j)
        }####STOP assigning new clonotype id
      }

    } ####stop hierarchical == T

    if(output.format=="dataframe.per.sample"){
      return(sample_dfs)
    }
    else if(output.format=="vgm"){
      if(!global.clonotype) VDJ.GEX.matrix <- do.call("rbind",sample_dfs)
      if(global.clonotype) VDJ.GEX.matrix <- sample_dfs[[1]]

      #shift new columns to the front and rename
      clono_10x_index <- which(names(VDJ.GEX.matrix) == "clonotype_id_10x")
      VDJ.GEX.matrix<- VDJ.GEX.matrix[,c(1:clono_10x_index, ((ncol(VDJ.GEX.matrix)-3):ncol(VDJ.GEX.matrix)), (clono_10x_index+1):(ncol(VDJ.GEX.matrix)-4))]

      names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonotype_id")] <- paste0("clonotype_id_",clone.strategy.as.input)
      names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_feature")] <- paste0("clonal_feature_",clone.strategy.as.input)
      names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_frequency")] <- paste0("clonotype_frequency_",clone.strategy.as.input)
      names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_rank")] <- paste0("clonal_rank_",clone.strategy.as.input)

      return(VDJ.GEX.matrix)
    }

    #convert into a data frame where each clone has a single row
    else if(output.format=="clone.level.dataframes" & global.clonotype == F){####START clone.level.dataframes
      clone.dataframe.list <- list()
      for(i in 1:length(sample_dfs)){
        sample_dfs[[i]]$VDJ_VJ_trimmed <- paste0(sample_dfs[[i]]$VDJ_sequence_nt_trimmed,sample_dfs[[i]]$VJ_sequence_nt_trimmed)
        clones_unique <- (sample_dfs[[i]][!duplicated(sample_dfs[[i]]$clonotype_id),])
        clones_unique$count.VDJ_VJ_trimmed_majority <- rep(NA,nrow(clones_unique))
        clones_unique$VDJ_VJ_trimmed_majority <- rep(NA,nrow(clones_unique))
        clones_unique$VDJ_trimmed_majority <- rep(NA,nrow(clones_unique))
        clones_unique$count.unique.trimVH.trimVL <- rep(NA,nrow(clones_unique))

        for (k in 1:nrow(clones_unique)){
          cells.per.clone <- sample_dfs[[i]][sample_dfs[[i]]$clonotype_id %in% clones_unique$clonotype_id[k], ]
          cells.per.clone.stats.VDJ_VJ <- sort(table(cells.per.clone$VDJ_VJ_trimmed),decreasing = T)
          cells.per.clone.stats.VDJ <- sort(table(cells.per.clone$VDJ_sequence_nt_trimmed),decreasing = T)
          cells.per.clone.stats.isotype <- sort(table(cells.per.clone$VDJ_cgene),decreasing = T)

          clones_unique$count.VDJ_VJ_trimmed_majority[k] <- cells.per.clone.stats.VDJ_VJ[1]
          clones_unique$VDJ_VJ_trimmed_majority[k] <- names(cells.per.clone.stats.VDJ_VJ)[1]
          clones_unique$count.unique.trimVH.trimVL[k] <- length(unique(cells.per.clone$VDJ_VJ_trimmed))
          clones_unique$count.VDJ_trimmed_majority[k] <- cells.per.clone.stats.VDJ[1]
          clones_unique$VDJ_trimmed_majority[k] <- names(cells.per.clone.stats.VDJ)[1]
          clones_unique$count.unique.trimVH[k] <- length(unique(cells.per.clone$VDJ_sequence_nt_trimmed))
          clones_unique$VDJ_cgene[k] <- names(cells.per.clone.stats.isotype)[1]
        }
        clone.dataframe.list[[i]] <- clones_unique

      }
      return(clone.dataframe.list)

    }####STOP clone.level.dataframes per Sample

    else if(output.format=="clone.level.dataframes" & global.clonotype == T){####START clone.level.dataframes
      sample_dfs$VDJ_VJ_trimmed <- paste0(sample_dfs$VDJ_sequence_nt_trimmed,sample_dfs$VJ_sequence_nt_trimmed)
      clones_unique <- (sample_dfs[!duplicated(sample_dfs$clonotype_id),])
      clones_unique$count.VDJ_VJ_trimmed_majority <- rep(NA,nrow(clones_unique))
      clones_unique$VDJ_VJ_trimmed_majority <- rep(NA,nrow(clones_unique))
      clones_unique$VDJ_trimmed_majority <- rep(NA,nrow(clones_unique))
      clones_unique$count.unique.trimVH.trimVL <- rep(NA,nrow(clones_unique))

      for (k in 1:nrow(clones_unique)){

        cells.per.clone <- sample_dfs[sample_dfs$clonotype_id %in% clones_unique$clonotype_id[k], ]
        cells.per.clone.stats.VDJ_VJ <- sort(table(cells.per.clone$VDJ_VJ_trimmed),decreasing = T)
        cells.per.clone.stats.VDJ <- sort(table(cells.per.clone$VDJ_sequence_nt_trimmed),decreasing = T)
        cells.per.clone.stats.isotype <- sort(table(cells.per.clone$VDJ_cgene),decreasing = T)

        clones_unique$count.VDJ_VJ_trimmed_majority[k] <- cells.per.clone.stats.VDJ_VJ[1]
        clones_unique$VDJ_VJ_trimmed_majority[k] <- names(cells.per.clone.stats.VDJ_VJ)[1]
        clones_unique$count.unique.trimVH.trimVL[k] <- length(unique(cells.per.clone$VDJ_VJ_trimmed))
        clones_unique$count.VDJ_trimmed_majority[k] <- cells.per.clone.stats.VDJ[1]
        clones_unique$VDJ_trimmed_majority[k] <- names(cells.per.clone.stats.VDJ)[1]
        clones_unique$count.unique.trimVH[k] <- length(unique(cells.per.clone$VDJ_sequence_nt_trimmed))
        clones_unique$VDJ_cgene[k] <- names(cells.per.clone.stats.isotype)[1]
      }
      return(clones_unique)
    }####STOP clone.level.dataframes (global)

    else if(output.format=="phylo.dataframes"){####START phylo.dataframes
      phylo.dataframe.list <- list()
      for(i in 1:length(sample_dfs)){
        if(VDJ.VJ.1chain==T){
          sample_dfs[[i]] <- sample_dfs[[i]][which(sample_dfs[[i]]$clonotype_id!="clonotypeNA"),]
        }
        phylo.dataframe.list[[i]] <- split(sample_dfs[[i]],sample_dfs[[i]]$clonotype_id)
        ## need to split the dataframe by clonotype_id column into lists
      }
      return(phylo.dataframe.list)
    }####STOP phylo.dataframes
  }####STOP v3
}

