#' Expand the aberrant cells in a VDJ dataframe by converting them into additional rows


#'@description  Expand the aberrant cells in a VDJ dataframe by converting them into additional rows. Aberrant cells consist of cells with more than 1 VDJ or VJ chain.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param chain.to.expand string, 'VDJ' to expand VDJ aberrants, 'VJ' to expand VJ aberrants, 'VDJ.VJ' for both.
#' @param add.barcode.prefix boolean - if T, a new barcode will be added for each expanded aberrant.
#' @param additional.VDJ.features vector of strings - VDJ_expand_aberrants will only expand across the sequence columns of VDJ. If you have additional columns with aberrant cell features (e.g., both 'yes' and 'no' binders for a single sequence), where the aberrants are VDJ-specific, include them here.
#' @param additional.VJ.features vector of strings - VDJ_expand_aberrants will only expand across the sequence columns of VDJ. If you have additional columns with aberrant cell features (e.g., both 'yes' and 'no' binders for a single sequence), where the aberrants are VJ-specific, include them here.
#' @param add.CDR3aa boolean - if T, will create a new column 'CDR3aa' with pasted VDJ_cdr3s_aa and VJ_cdr3s_aa.
#' @param add.expanded.number boolean - if T, will add the number of new cells resulting from an aberrant one.
#' @param recalculate.clonotype.frequency boolean - if T, will recalculate the clonotype frequencies for the resulting, expanded VDJ.
#' @return Returns a VDJ format dataframe in which cells with more than one VDJ or VJ chain are split into multiple rows each containing only one VDJ VJ chain combination.
#' @export
#' @examples
#' VDJ_expand_aberrants(VDJ = Platypus::small_vgm[[1]],
#' chain.to.expand='VDJ.VJ',
#' add.barcode.prefix=TRUE, recalculate.clonotype.frequency=FALSE)
#'

VDJ_expand_aberrants  <- function(VDJ,
                                  chain.to.expand,
                                  add.barcode.prefix,
                                  additional.VDJ.features,
                                  additional.VJ.features,
                                  add.CDR3aa,
                                  add.expanded.number,
                                  recalculate.clonotype.frequency){

  if(missing(VDJ)) stop('Please input your data as a VDJ matrix')
  if(missing(chain.to.expand)) chain.to.expand <- 'VDJ.VJ'
  if(missing(add.barcode.prefix)) add.barcode.prefix <- TRUE
  if(missing(additional.VDJ.features)) additional.VDJ.features <- NULL
  if(missing(additional.VJ.features)) additional.VJ.features <- NULL
  if(missing(add.CDR3aa)) add.cdr3aa <- TRUE
  if(missing(add.expanded.number)) add.expanded.number <- TRUE
  if(missing(recalculate.clonotype.frequency)) recalculate.clonotype.frequency <- TRUE

  get_sequence_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';'){
    if(split.x==TRUE) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
    if(split.y==TRUE) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

    ccombs <- expand.grid(x,y)
    ccombs<-paste0(ccombs[,1], ccombs[,2])
    ccombs <- paste0(ccombs, collapse=collapse.by)

    return(ccombs)
  }

  VDJ_columns <- list('Nr_of_VDJ_chains', 'VDJ_cdr3s_aa', 'VDJ_cdr3s_nt', 'VDJ_chain_contig', 'VDJ_chain', 'VDJ_vgene', 'VDJ_dgene', 'VDJ_jgene', 'VDJ_cgene', 'VDJ_sequence_nt_raw', 'VDJ_sequence_nt_trimmed', 'VDJ_sequence_aa', 'VDJ_trimmed_ref')
  VJ_columns <- list('Nr_of_VJ_chains', 'VJ_cdr3s_aa', 'VJ_cdr3s_nt', 'VJ_chain_contig', 'VJ_chain', 'VJ_vgene', 'VJ_jgene', 'VJ_cgene', 'VJ_sequence_nt_raw', 'VJ_sequence_nt_trimmed', 'VJ_sequence_aa', 'VJ_trimmed_ref')

  VDJ_columns <- unlist(c(VDJ_columns, additional.VDJ.features))
  VJ_columns <- unlist(c(VJ_columns, additional.VJ.features))


  VDJ.matrix <- VDJ
  VDJ <- NULL

  VDJ_columns <- base::intersect(VDJ_columns, colnames(VDJ.matrix))
  VJ_columns <- base::intersect(VDJ_columns, colnames(VDJ.matrix))


  if(chain.to.expand=='VDJ'){
    VDJ.matrix <- tidyr::separate_rows(VDJ.matrix, tidyselect::all_of(VDJ_columns), sep=';', convert=TRUE)


  }else if(chain.to.expand=='VJ'){
    VDJ.matrix <- tidyr::separate_rows(VDJ.matrix, tidyselect::all_of(VJ_columns), sep=';', convert=TRUE)

  }else if(chain.to.expand=='VDJ.VJ'){
    VDJ.matrix <- tidyr::separate_rows(VDJ.matrix, tidyselect::all_of(VDJ_columns), sep=';', convert=TRUE)
    VDJ.matrix <- tidyr::separate_rows(VDJ.matrix, tidyselect::all_of(VJ_columns), sep=';', convert=TRUE)

  }else stop('Chain not found.')

  VDJ.matrix <- as.data.frame(VDJ.matrix)

  if(add.cdr3aa){
    VDJ.matrix$CDR3aa <- mapply(function(x,y) get_sequence_combinations(x,y, split.x=TRUE, split.y=TRUE), VDJ.matrix$VDJ_cdr3s_aa, VDJ.matrix$VJ_cdr3s_aa)
  }

  if(add.expanded.number){
    VDJ.matrix$expanded_number <-  rep(NA, nrow(VDJ.matrix))
  }

  if(add.barcode.prefix){
    ind <- which(duplicated(VDJ.matrix$barcode))
    ind <- unname(split(ind, cumsum(c(1, diff(ind) != 1))))

    for(i in 1:length(ind)){
      barcode_ind <- 1
      split_barcode <- stringr::str_split(VDJ.matrix$barcode[ind[[i]][1]-1], '_', simplify=TRUE)
      VDJ.matrix$barcode[ind[[i]][1]-1] <- paste0(split_barcode[1], '_', barcode_ind, '_', split_barcode[2])

      if(add.expanded.number){
        VDJ.matrix$expanded_number[ind[[i]][1]-1] <- length(ind[[i]]) + 1
        VDJ.matrix$expanded_number[ind[[i]]] <- length(ind[[i]]) + 1
      }

      for(j in 1:length(ind[[i]])){
        barcode_ind <- barcode_ind + 1
        split_barcode <- stringr::str_split(VDJ.matrix$barcode[ind[[i]][j]], '_', simplify=TRUE)
        VDJ.matrix$barcode[ind[[i]][j]] <- paste0(split_barcode[1], '_', barcode_ind, '_', split_barcode[2])
      }
    }
  }

  if(recalculate.clonotype.frequency){
    repertoire_numbers <- unique(VDJ.matrix$sample_id)
    for(rep in repertoire_numbers){
      clonotypes <- VDJ.matrix$clonotype_id[which(VDJ.matrix$sample_id==rep)]
      frequencies <- unlist(lapply(clonotypes, function(x) length(which(clonotypes==x))))
      VDJ.matrix$clonotype_frequency[which(VDJ.matrix$sample_id==rep)] <- frequencies
      VDJ.matrix[which(VDJ.matrix$sample_id==rep),] <- VDJ.matrix[which(VDJ.matrix$sample_id==rep),][order(VDJ.matrix$clonotype_frequency[which(VDJ.matrix$sample_id==rep)], decreasing=TRUE),]
    }
  }

  VDJ.matrix$X <- 1:nrow(VDJ.matrix)
  return(VDJ.matrix)

}
