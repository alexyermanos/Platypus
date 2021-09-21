#' Processes and organizes the repertoire sequening data from cellranger vdj and returns a list of dataframes, where each dataframe corresponds to an individual repertoire. The function will return split CDR3 sequences, germline gene information, filter out those clones with either incomplete information or doublets (multiple CDR3 sequences for a given chain). This function should be called once for desired integrated repertoire and transcriptome. For example, if there are 3 VDJ libraries and 3 GEX libraries and the goal is to analyze all three GEX libraries together (e.g. one UMAP/tSNE reduction) this then function should be called one time and the three VDJ directories should be provided as input to the single function call.
#' @param VDJ.out.directory Character vector with each element containing the path to the output of cellranger vdj runs. Multiple repertoires to be integrated in a single transcriptome should be supplied as multiple elements of the character vector. This can be left blank if supplying the clonotypes and contig files directly as input. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @param filter.1HC.1LC Logical indicating whether only those clones containing 1 VH/TRB and VL/TRA should be maintined for furhter analysis. Default is set to TRUE, which restricts the analysis to only clones with exactly 1 heavy chain and 1 light chain (or 1 beta + 1 alpha in the case of T cells).
#' @param clonotype.list List of dataframes containing clonotyping information for each repertoire. The column names should correspond to the clonotypes.csv file from cellranger vdj output.
#' @param contig.list List of dataframes containing the contig information for each repertoire.  The column names should correspond to the all_contigs.csv file from cellranger vdj output.
#' @param filtered.contigs Logical indicating if the filtered contigs file should be used. TRUE will read VDJ information from only the filtered output of cellranger. FALSE will read the all contigs file from cellranger. Default set to TRUE (filtered output)
#' @return Returns a list of dataframes where each dataframe corresponds to one input directory. If only one file is supplied, the output list will only contain one element. This output can be supplied as input to other functions including VDJ_per_clone, VDJ_network, VDJ_germline_genes, VDJ_expansion, visualize_clones_GEX, VDJ_phylo, VDJ_clonotype. Germline gene information is based on the majority of cells within each clonotype. For example, if the majority of cells in clonotype1 have the IGHG1 isotype then then entire clonal family will be determined as IGHG1. For a cell-specific investigation, the output of this function can be supplied to the function VDJ_per_clone, which will provide isotype, sequence, germline gene, etc information for each cell within the each clone.
#' @export
#' @examples
#' \dontrun{
#' example.vdj.analyze <- VDJ_analyze(
#' VDJ.out.directory = "~/path/to/cellranger/vdj/outs/", filter.1HC.1LC = T)
#' }
#'
VDJ_analyze <- function(VDJ.out.directory,
                        filter.1HC.1LC,
                        clonotype.list,
                        contig.list,
                        filtered.contigs){
  if(missing(VDJ.out.directory)) print("No output directory supplied. Assuming clonotype and contig are provided as list objects")
  if(missing(filtered.contigs)) filtered.contigs <- TRUE
  if(missing(filter.1HC.1LC)) filter.1HC.1LC <- TRUE

  if(missing(VDJ.out.directory)==F){
    VDJ.out.directory_clonotypes <- paste(VDJ.out.directory,"/clonotypes.csv",sep="")
    if(filtered.contigs==T)VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/filtered_contig_annotations.csv",sep="")
    if(filtered.contigs==F)VDJ.out.directory_contigs <- paste(VDJ.out.directory,"/all_contig_annotations.csv",sep="")


    clonotype.list <- lapply(VDJ.out.directory_clonotypes, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
    contig.list <- lapply(VDJ.out.directory_contigs, function(x) utils::read.table(x, stringsAsFactors = FALSE,sep=",",header=T))
  }
  for(i in 1:length(clonotype.list)){

    clonotype.list[[i]]$HC_count <- stringr::str_count(pattern = "IGH:",string = clonotype.list[[i]]$cdr3s_aa)
    clonotype.list[[i]]$HC_count <- clonotype.list[[i]]$HC_count + stringr::str_count(pattern = "TRB:",string = clonotype.list[[i]]$cdr3s_aa)


    clonotype.list[[i]]$IGK_count <- stringr::str_count(pattern = "IGK:",string = clonotype.list[[i]]$cdr3s_aa)
    clonotype.list[[i]]$IGL_count <- stringr::str_count(pattern = "IGL:",string = clonotype.list[[i]]$cdr3s_aa)
    clonotype.list[[i]]$LC_count <- clonotype.list[[i]]$IGK_count + clonotype.list[[i]]$IGL_count
    clonotype.list[[i]]$LC_count <- clonotype.list[[i]]$LC_count + stringr::str_count(pattern = "TRA:",string = clonotype.list[[i]]$cdr3s_aa)

    if(filter.1HC.1LC==T){ ## Filter out those clones that do not contain exactly one heavy/beta and one light/alpha chain.
      clonotype.list[[i]] <- clonotype.list[[i]][clonotype.list[[i]]$HC_count==1 & clonotype.list[[i]]$LC_count==1,]
    }

    clonotype.list[[i]]$CDRH3_aa <- rep("",nrow(clonotype.list[[i]]))
    clonotype.list[[i]]$CDRL3_aa <- rep("",nrow(clonotype.list[[i]]))
    clonotype.list[[i]]$CDRH3_nt <- rep("",nrow(clonotype.list[[i]]))
    clonotype.list[[i]]$CDRH3_nt <- rep("",nrow(clonotype.list[[i]]))
    for(j in 1:nrow(clonotype.list[[i]])){
      clonotype.list[[i]]$CDRH3_aa[j] <- gsub(pattern = "IGH:",x=stringr::str_split(clonotype.list[[i]]$cdr3s_aa[j],pattern = ";")[[1]][1],replacement = "")
      clonotype.list[[i]]$CDRL3_aa[j] <- gsub(pattern="IGL:",x=gsub(pattern = "IGK:",x=stringr::str_split(clonotype.list[[i]]$cdr3s_aa[j],pattern = ";")[[1]][2],replacement = ""),replacement = "")
      clonotype.list[[i]]$CDRH3_nt[j] <- gsub(pattern = "IGH:",x=stringr::str_split(clonotype.list[[i]]$cdr3s_nt[j],pattern = ";")[[1]][1],replacement = "")
      clonotype.list[[i]]$CDRL3_nt[j] <- gsub(pattern="IGL:",x=gsub(pattern = "IGK:",x=stringr::str_split(clonotype.list[[i]]$cdr3s_nt[j],pattern = ";")[[1]][2],replacement = ""),replacement = "")
    }
    clonotype.list[[i]]$CDR3_aa_pasted <- paste(clonotype.list[[i]]$CDRH3_aa,clonotype.list[[i]]$CDRL3_aa,sep="")
    clonotype.list[[i]]$CDR3_nt_pasted <- paste(clonotype.list[[i]]$CDRH3_nt,clonotype.list[[i]]$CDRL3_nt,sep="")

  clonotype.list[[i]]$HC_cgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$HC_vgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$HC_dgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$HC_jgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$LC_cgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$LC_vgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$LC_jgene <- rep("",nrow(clonotype.list[[i]]))
  clonotype.list[[i]]$barcodes <- rep("",nrow(clonotype.list[[i]]))

  temp_trb <- which(contig.list[[i]]$chain=="TRB")
  temp_tra <- which(contig.list[[i]]$chain=="TRA")
  contig.list[[i]]$chain[temp_trb] <- "IGH"
  contig.list[[i]]$chain[temp_tra] <- "IGK"
  holding_bar <- utils::txtProgressBar(min = 0, max = 1, initial = 0, char = "=",
                                width = NA, style = 1, file = "")
  print(paste(i,"from",length(clonotype.list),"repertoires"))
  for(j in 1:nrow(clonotype.list[[i]])){
    utils::setTxtProgressBar(value = j/nrow(clonotype.list[[i]]),pb = holding_bar)
    tryCatch({
      clonotype.list[[i]]$HC_cgene[j] <- names(which.max(table(contig.list[[i]]$c_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH")])))
      clonotype.list[[i]]$HC_vgene[j] <- names(which.max(table(contig.list[[i]]$v_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH" & contig.list[[i]]$v_gene!="None")])))
      clonotype.list[[i]]$HC_dgene[j] <- names(which.max(table(contig.list[[i]]$d_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH")])))
      clonotype.list[[i]]$HC_jgene[j] <- names(which.max(table(contig.list[[i]]$j_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain=="IGH")])))
      clonotype.list[[i]]$LC_cgene[j] <- names(which.max(table(contig.list[[i]]$c_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH")])))
      clonotype.list[[i]]$LC_vgene[j] <- names(which.max(table(contig.list[[i]]$v_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH")])))
      clonotype.list[[i]]$LC_jgene[j] <- names(which.max(table(contig.list[[i]]$j_gene[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true") & contig.list[[i]]$chain!="IGH")])))
      clonotype.list[[i]]$barcodes[j] <- gsub(toString(unique(contig.list[[i]]$barcode[which(contig.list[[i]]$raw_clonotype_id==clonotype.list[[i]]$clonotype_id[j] & stringr::str_detect(contig.list[[i]]$is_cell, "(?i)true"))])),pattern = ", ",replacement = ";")
    }, error=function(e){})

    }
    clonotype.list[[i]]$nt_clone_ids <- clonotype.list[[i]]$clonotype_id
  }

  return(clonotype.list)
}
