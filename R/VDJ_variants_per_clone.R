#' Returns statistics and plots to examine diversity of any sequence or metadata item within clones on a by sample level or global level
#' @param VDJ.matrix VDJ dataframe generated using the VDJ_GEX_matrix. VDJ.matrix supplemented with with MIXCR information is also valid
#' @param variants.of Character vector. Defaults to c("VDJ_cdr3s_aa", "VJ_cdr3s_aa"). Column name(s) of VDJ.matrix to examine variants of. If more than one name is given, these columns will be pasted together. The default will therefore return statistics on the number of variants of VDJ and VJ cdr3s in every clone 
#' @param clonotypes.col Column name of the VDJ.matrix column containing clonotype information. Defaults to "clonotype_id_10x". This is useful if alternative clonotyping strategies have been used and are stored in other columns
#' @param stringDist.method Character. Passed to Biostrings::strinDist. Method to calculate distance between variants of a clone. Defaults to "levenshtein". Other options are "hamming", "quality". If "hamming" variants of a clone will be shortened from the end to the shortest variant to make all input sequences the same length.
#' @param split.by Character. Defaults to "sample_id". Column name of VDJ.matrix to split the analysis by. This is necessary, if clonotyping was done on a per sample level (e.g. "clonotype1" in sample 1 is not the same sequence as "clonotype1" in sample 2). If clonotyping was done across samples and no splitting is necessary input "none" 
#' @param platypus.version Character. Only "v3" available.  
#' @return Returns a list of dataframes. Each dataframe contains the statistics of one split.by element (in most cases, one sample)
#' @export
#' @examples
#' \dontrun{
#'
#'variants_per_clone <- VDJ_variants_per_clone(VDJ.matrix = VDJ_comb[[1]],variants.of = c("VDJ_cdr3s_aa", "VJ_cdr3s_aa"),stringDist.method = "levenshtein", split.by = "sample_id")
#'take a sneakpeak
#'for(i in 1:length(variants_per_clone)){print(variants_per_clone[[i]][1:10,])}
#'}

VDJ_variants_per_clone <- function(VDJ.matrix, 
                                   variants.of,
                                   clonotypes.col,
                                   stringDist.method,
                                   split.by,
                                   platypus.version){
  
  #function to use in lapply iterating over clonotypes
  find_vars <- function(clonotypes, dat){
    
    #for clones 
    cur_clone <- subset(dat, clonotype == clonotypes) 
    n_var <- length(unique(cur_clone$pasted)) #number of variants
    n_tot <- nrow(cur_clone) #nr of total sequences for that clone
    most_frequent <- names(which.max(table(cur_clone$pasted)))
    if(n_var > 1){
      #trim sequences in case for hamming method
      if(stringDist.method == "hamming" & length(unique(nchar(unique(cur_clone$pasted)))) > 1){
        cur_clone$pasted <- substr(cur_clone$pasted, 0, min(nchar(unique(cur_clone$pasted))))
      }
      dist_mat <- Biostrings::stringDist(unique(cur_clone$pasted), method = stringDist.method)
      m_stringdist <- mean(dist_mat)
      max_stringdist <- max(dist_mat)
      min_stringdist <- min(dist_mat)
    } else {
      m_stringdist <- NA
      max_stringdist <- NA
      min_stringdist <- NA
    }
    
    #for isotypes
    n_iso_var <- length(unique(cur_clone$VDJ_cgene)) #number of variants
    most_frequent_iso <- names(which.max(table(cur_clone$VDJ_cgene)))
    
    out <- c(as.character(cur_clone$split[1]),clonotypes, n_tot, n_var, m_stringdist, max_stringdist, min_stringdist,most_frequent, n_iso_var, most_frequent_iso)
    names(out) <- c("split_group","clonotype_id","n sequences", "n variants", "mean variants distance", "max variants distance", "min variants distance","most frequent variant", "n VDJ cgene variants", "most frequent VDJ cgene")
    return(out)
  }
  
platypus.version <- "v3"
if(missing(split.by)) split.by <- "sample_id"
if(missing(variants.of)) variants.of <- c("VDJ_cdr3s_aa", "VJ_cdr3s_aa")
if(missing(clonotypes.col)) clonotypes.col <- "clonotype_id_10x"
if(missing(stringDist.method)) stringDist.method <- "levenshtein"
  
  #get data
if(any(!variants.of %in% names(VDJ.matrix))){stop("At least one of the input columns does not exist in the VDJ.matrix. Please input valid column names to variants.of")}
if(!clonotypes.col %in% names(VDJ.matrix)){stop("clonotypes.col not found in VDJ.matrix. Please input valid column name")}
if(!split.by %in% names(VDJ.matrix) & split.by != "none"){stop("split.by not found in VDJ.matrix. Please input valid column name")}
  
#remove any rows that do not contain an entry for a given feature
to_remove <- c()
for(n in 1:nrow(VDJ.matrix)){
  if("" %in% VDJ.matrix[n,c(variants.of)]){
    to_remove <- c(to_remove, n)}
}
if(length(to_remove) > 0){
  VDJ.matrix <- VDJ.matrix[-to_remove,]
}

#get dataframe of clonotype ids and pasted variants.of columns
grouping <- data.frame("clonotype" = VDJ.matrix[, clonotypes.col])
if(length(variants.of) > 1){
  grouping$pasted <- do.call(paste, c(VDJ.matrix[,c(variants.of)], sep="/"))
} else {
  grouping$pasted <- VDJ.matrix[, c(variants.of)]
}
  
#add isotype information
if("VDJ_cgene" %in% names(VDJ.matrix)){
  grouping$VDJ_cgene <- VDJ.matrix$VDJ_cgene
} else {
  print("VDJ_cgene column not found; Isotype info will not be attached")
}

#loop over split.by
if(split.by == "none"){grouping$split <- "1"
} else {grouping$split <- VDJ.matrix[,c(split.by)]}

out_list <- list()
for(j in 1:length(unique(grouping$split))){
  
  grouping_sub <- subset(grouping, split == unique(grouping$split)[j])
  
  out_vars <- lapply(X = unique(grouping_sub$clonotype), FUN = find_vars, grouping_sub)
  out_vars <- as.data.frame(bind_rows(out_vars))
  
  for(i in 3:(ncol(out_vars)-1)){
    if(names(out_vars)[i] != "most frequent variant"){
    out_vars[,i] <- as.numeric(out_vars[,i])}
  }
  out_vars[,"mean variants distance"] <- round(out_vars[,"mean variants distance"],2)
  
  out_vars <- out_vars[order(out_vars[,"n variants"], decreasing = T),]
  
  out_list[[j]] <- out_vars
  names(out_list)[[j]] <- unique(grouping$split)[j]
  }
  return(out_list)
}



