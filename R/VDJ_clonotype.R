#' Returns a list of clonotype dataframes following additional clonotyping. This function works best following filtering to ensure that each clone only has one heavy chain and one light chain.
#' @param clonotype.list Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @param clone.strategy String describing the clonotyping strategy. Possible options include 'cdr3.nt', 'cdr3.aa','hvj.lvj','hvj.lvj.cdr3lengths','Hvj.Lvj.CDR3length.CDR3homology', 'Hvj.Lvj.CDR3length.CDRH3homology', 'CDR3homology',or 'CDRH3homology'. 'CDR3aa' will convert the default cell ranger clonotyping to amino acid based. 'Hvj.Lvj' groups B cells with identical germline genes (V and J segments for both heavy chain and light chain. Those arguments including 'CDR3length' will group all sequences with identical CDRH3 and CDRL3 sequence lengths. Those arguments including 'CDR3homology' will additionally impose a homology requirement for CDRH3 and CDRL3 sequences.'CDR3homology',or 'CDRH3homology' will group sequences based on homology only (either of the whole CDR3 sequence or of the CDRH3 sequence respictevely).
#' All homology calculations are performed on the amino acid level.
#' @param homology.threshold Numeric value between 0 and 1 corresponding to the homology threshold forn the clone.strategy arguments that require a homology threshold. Default value is set to 70 percent sequence homology. For 70 percent homology, 0.3 should be supplied as input.
#' @param platypus.version Default is "v2" for compatibility. To use the output of VDJ_GEX_matrix function, one should change this argument to "v3".
#' @param VDJ.GEX.matrix Output from the VDJ.GEX.matrix function. The output object should have the VDJ information (e.g., the original VDJ_GEX_matrix call should have had cellranger's VDJ output supplied as input).
#' @param output.dataframe Logical that specifies whether the output should return a list element of dataframes, where each dataframe contains the cell-level information for a given clone. By default this is set to false, which will instead update the "clonotype_id" column in the VDJ.GEX.matrix[[1]] element.
#' @param global.clonotype Logical specifying whether clonotyping should occur across samples or only within a single sample.
#' @param VDJ.VJ.1chain Logical specifying whether cells with multiple VDJ and VJ chains should be removed from the clonotyping. Can be either T or F for those definitions not requiring germline genes or homology thresholds, as calculating the later is difficult when multiple chains are present.
#' @return Returns a list of clonotype dataframes where each list element matches the  repertoire index in the input clonotype.list object. The dataframes will be updated with clonal frequencies based on the new clonotyping definition.
#' @export
#' @examples
#' \dontrun{
#' VDJ_clonotype(clonotype.list=VDJ.analyze.output,
#'  clone.strategy="Hvj.Lvj.CDR3length.CDR3homology",
#'  homology.threshold=".3")
#'  }
#'
VDJ_clonotype <- function(clonotype.list,
                          clone.strategy,
                          homology.threshold,
                          platypus.version,
                          VDJ.GEX.matrix,
                          output.dataframe,
                          global.clonotype,
                          VDJ.VJ.1chain){
  require(stringdist)
  if(missing(platypus.version)) platypus.version <- "v2"
  if(missing(VDJ.GEX.matrix)) VDJ.GEX.matrix <- list()
  if(missing(output.dataframe)) output.dataframe <- FALSE
  if(missing(global.clonotype)) global.clonotype <- FALSE
  if(missing(clone.strategy)) clone.strategy <- "cdr3.nt"
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T

  if(platypus.version=="v2"){####START v2
    output.clonotype <- list()
    if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) print("No homology threshold supplied. Clonotyping based on 70% amino acid homology.")
    if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) homology.threshold<-0.3   # Setting default homology threshold

    #Possible strategy options:'cdr3.aa','hvj.lvj','hvj.lvj.cdr3lengths','Hvj.Lvj.CDR3length.CDR3homology', 'Hvj.Lvj.CDR3length.CDRH3homology', 'CDR3homology',or 'CDRH3homology'.
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
      else if(clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology" | clone.strategy=="Hvj.Lvj.CDR3length.CDRH3homology"){  #taking into account both cases
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
            if (clone.strategy=="Hvj.Lvj.CDR3length.CDR3homology"){
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
    require(parallel)
    if(global.clonotype==F){ # loop through each repertoire individualually
      repertoire.number <- unique(VDJ.GEX.matrix[[1]]$sample_id)
      sample_dfs <- list()
      for(i in 1:length(repertoire.number)){ ####START sample loop
        sample_dfs[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==repertoire.number[i]),]
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

        } ####STOP hvj.lvj.cdr3lengths
        ####START homology based clonotyping - to add still for v3
        ####STOP homology based clonotyping - to add still for v3

        ####START recalculating clonotype_id and clonal_frequency
        if(VDJ.VJ.1chain==T) sample_dfs[[i]]$new_clonal_feature[which(sample_dfs[[i]]$Nr_of_VDJ_chains!=1 & sample_dfs[[i]]$Nr_of_VJ_chains!=1)] <- NA
        sample_dfs[[i]]$new_clonal_frequency <- rep(NA,nrow(sample_dfs[[i]]))
        sample_dfs[[i]]$new_clonal_rank <- rep(NA,nrow(sample_dfs[[i]]))
        sample_dfs[[i]]$new_tied_clonal_rank <- rep(NA,nrow(sample_dfs[[i]]))

        unique.clonal.features <- unique(sample_dfs[[i]]$new_clonal_feature)
        unique.clonal.frequencies <- rep(NA,length(unique.clonal.features))
        for(j in 1:length(unique.clonal.features)){####START assigning new frequency
          unique.clonal.frequencies[j] <- length(which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j]))
          sample_dfs[[i]]$new_clonal_frequency[which(sample_dfs[[i]]$new_clonal_feature==unique.clonal.features[j])] <- unique.clonal.frequencies[j]
        }####STOP assigning new frequency

        ####START new clonotype_id
        new.clonal.order <- unique.clonal.features[order(unique.clonal.frequencies,decreasing = T)]
        for(j in 1:length(new.clonal.order)){
          sample_dfs[[i]]$new_clonal_rank[which(sample_dfs[[i]]$new_clonal_feature==new.clonal.order[j])] <- j

        }
        sample_dfs[[i]]$clonotype_id <- paste0("clonotype",sample_dfs[[i]]$new_clonal_rank)
        ####STOP new clonotype_id

        ####START clonal rank including ties
        #e.g. if two different clones have same number of cells, they have the same rank
        unique.clone.frequencies <- sort(unique(sample_dfs[[i]]$new_clonal_frequency),decreasing = T)
        for(j in 1:length(unique.clone.frequencies)){
          sample_dfs[[i]]$new_tied_clonal_rank[which(sample_dfs[[i]]$new_clonal_frequency==unique.clone.frequencies[j])] <- j
        }
        ####STOP clonal rank including ties

      }####STOP sample loop
    }####STOP global.clonotype==F
    else if(global.clonotype==T){####START global.clonotype==T

    }####STOP global.clonotype==T

    if(output.dataframe==T){
      return(sample_dfs)
    }
    else if(output.dataframe==F){

      VDJ.GEX.matrix[[1]] <- do.call("rbind",sample_dfs)
      return(VDJ.GEX.matrix)
    }


  }####STOP v3
}
