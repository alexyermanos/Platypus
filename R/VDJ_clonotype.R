#' Returns a list of clonotype dataframes following additional clonotyping. This function works best following filtering to ensure that each clone only has one heavy chain and one light chain.
#' @param clonotype.list Output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire.
#' @param clone.strategy String describing the clonotyping strategy. Possible options include 'cdr3.aa','hvj.lvj','hvj.lvj.cdr3lengths','Hvj.Lvj.CDR3length.CDR3homology', 'Hvj.Lvj.CDR3length.CDRH3homology', 'CDR3homology',or 'CDRH3homology'. 'CDR3aa' will convert the default cell ranger clonotyping to amino acid based. 'Hvj.Lvj' groups B cells with identical germline genes (V and J segments for both heavy chain and light chain. Those arguments including 'CDR3length' will group all sequences with identical CDRH3 and CDRL3 sequence lengths. Those arguments including 'CDR3homology' will additionally impose a homology requirement for CDRH3 and CDRL3 sequences.'CDR3homology',or 'CDRH3homology' will group sequences based on homology only (either of the whole CDR3 sequence or of the CDRH3 sequence respictevely).
#' All homology calculations are performed on the amino acid level.
#' @param homology.threshold Numeric value between 0 and 1 corresponding to the homology threshold forn the clone.strategy arguments that require a homology threshold. Default value is set to 70 percent sequence homology. For 70 percent homology, 0.3 should be supplied as input.
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
                          homology.threshold){
  require(stringdist)

  output.clonotype <- list()
  if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) print("No homology threshold supplied. Clonotyping based on 70% amino acid homology.")
  if(missing(homology.threshold) & grepl(clone.strategy,pattern = "homology")) homology.threshold<-0.3   # Setting default homology threshold

  #Possible strategy options:'cdr3.aa','hvj.lvj','hvj.lvj.cdr3lengths','Hvj.Lvj.CDR3length.CDR3homology', 'Hvj.Lvj.CDR3length.CDRH3homology', 'CDR3homology',or 'CDRH3homology'.
  for(i in 1:length(clonotype.list)){
    if(clone.strategy=="cdr3.aa"){
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
}
