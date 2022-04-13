#' Integrates antigen-specific information into the VDJ/VDJ.GEX.matrix[[1]] object

#'@description Integrate antigen-specific information from a list of antigen dataframes or antigen csv file paths. The antigen data should contain either the clonotypes, cell barcodes, or sequences with the specific column names of the VDJ/VDJ.GEX.matrix[[1]] object. These columns will be used to rematch the binder information at the cell, sequence, or clonotype level into the main VDJ.GEX.matrix[[1]].
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param antigen.data.list list of antigen csv file paths or antigen dataframes for the specific antigen datasets. To ease matching, the column names by which we will match should be the same as the column names in the original VDJ/VDJ.GEX.matrix[[1]] object.
#' @param antigen.features vector of columns of antigen features to be integrated from the antigen csv files into the VDJ/VDJ.GEX.matrix[[1]] object. The vector can also use unique, short-hand names of the columns to add (e.g., 'affinity' for 'octet.affinity.[nM]').
#' @param binder.threshold list or nested list of threshold values and specific features by which to define binders in the VDJ.
#' For example, if binder.threshold=list(list('affinity', 0.2), list('elisa', 0.8)), we will have two new binder columns: binders_affinity if the values are greater than 0.2, binders_elisa if they are greater than 0.8.
#' @param VDJ.VJ.1chain boolean, if T will remove aberrant cells (more than 1 VDJ of VJ chain), if F it will keep them in the VDJ when matching antigen data.
#' @param match.by string, represents the method by which to match the antigen data and integrate it into the VDJ/VDJ.GEX.matrix[[1]] object. 'clonotype' will match by 'clonotype_id' (needs to be present in the antigen data), 'clonotype.v3' will match by v3 cellranger clonotypes (you need a v3_clonotypes column in the VDJ/VDJ.GEX.matrix[[1]], 'cdr3.aa' by VDJ and VJ cdr3s amino acid sequences, 'cdrh3.aa' by VDJ cdr3s amino acid sequences, 'VDJ.VJ.aa' by full VDJ and VJ aa sequences, 'VDJ.VJ.nt' by trimmed nt VDJ and VJ sequences (must run VDJ_call_MIXCR first on the VDJ),'cdr3.nt' by VDJ and VJ cdr3s as nucleotides, 'cdrh3.nt.' by VDJ cdr3s as nucleotides, 'absolut' will match the VDJ_cdr3s_aa with the CDR3 column in Absolut! datasets.
#' @param matching.type string, either 'exact' for exact sequence matching if the match.by parameter is a sequence type, or 'homology' for homology matching (matches if the Levehnstein distance is less than the distance.threshold parameter).
#' @param distance.threshold integer, maximum string distance value by which to match sequences in the antigen data and sequences in the VDJ object (to further integrate the antigen data).
#' @param sample.id boolean, if T then will also match by the 'sample_id' column in the antigen dataframes.
#' @param aberrant.chosen.sequences boolean, if T will add a column of the chosen aberrant sequences (which matched a sequence in the antigen data) if matching by sequence (and VDJ.VJ.1chain=F).
#' @param output.format string, 'vgm' - returns the full VDJ object, 'dataframe.per.sample' - list of VDJ dataframes for each sample.
#' @return Either the original VDJ dataframe with additional columns of the antigen features integrated, a list of VDJ dataframes per sample.
#' @export

#' @examples
#' \dontrun{
#' VDJ_antigen_integrate_v2(VDJ,antigen.directory.list=antigen.directory.list,
#' antigen.feature=c('elisa', 'affinity'),VDJ.VJ.1chain=T,
#' match.by='clonotype',sample.id=T, output.format='vgm')
#'}


VDJ_antigen_integrate <- function(VDJ,
                                 antigen.data.list,
                                 antigen.features,
                                 binder.threshold,
                                 VDJ.VJ.1chain,
                                 match.by,
                                 matching.type,
                                 distance.threshold,
                                 sample.id,
                                 aberrant.chosen.sequences,
                                 output.format){


  if(missing(VDJ)) stop('Please provide input data as VDJ')
  if(missing(antigen.data.list)) stop('Please provide a list of antigen data file(s) or a list of antigen dataframes')
  if(missing(antigen.features)) stop('Please provide the features to be extracted form the antigen data')
  if(missing(binder.threshold)) binder.threshold <- list()
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- F
  if(missing(match.by)) match.by <- 'clonotype'
  if(missing(matching.type)) matching.type <- 'exact'
  if(missing(distance.threshold)) distance.threshold <- 3
  if(missing(sample.id)) sample.id <- T
  if(missing(aberrant.chosen.sequences)) aberrant.chosen.sequences <- T
  if(missing(output.format)) output.format <- 'vgm'

  get_sequence_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';'){
   if(split.x==T) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
   if(split.y==T) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

   ccombs <- expand.grid(x,y)
   ccombs<-paste0(ccombs[,1], ccombs[,2])
   ccombs <- paste0(ccombs, collapse=collapse.by)

   return(ccombs)
  }


  extract_features_from_clonotype_matchings <- function(sample_df, antigen_df, clonotype_column,
                                                        features=antigen_df_columns, new_names=new_column_names,
                                                        binder_threshold=binder.threshold, sample_id=sample.id){

    if(sample_id){
      if(!(unique(sample_df$sample_id) %in% unique(antigen_df$sample_id))){
        for(i in 1:length(new_names)){
          sample_df[new_names[[i]]] <- rep(NA, nrow(sample_df))

           if(length(binder_threshold)!=0){
            for(threshold in binder_threshold){

              if(is.list(threshold) & length(threshold)==2 & stringr::str_detect(new_names[[i]], threshold[[1]])){
                new_binder_feature <- paste0(new_names[i], '_', 'binders')
                sample_df$new_binder_feature <- rep(NA, nrow(sample_df))
                names(sample_df)[names(sample_df) == "new_binder_feature"] <- new_binder_feature
              }
            }
           }
         }
      return(sample_df)
     }
      antigen_df <- antigen_df[which(antigen_df$sample_id==unique(sample_df$sample_id)),]
    }

    all_clonotypes <- sample_df[,clonotype_column]

    for(i in 1:length(features)){

      matched_feature_values <- sapply(all_clonotypes, function(x) if(x!='clonotypeNA') antigen_df[features[[i]]][which(antigen_df[clonotype_column]==x),]  else NA)

      matched_max_feature_values <- sapply(matched_feature_values, function(x) if(inherits(x[1],'numeric')) {if(!is.na(x[1])) x[which.max(x)] else NA} else if(!inherits(x[1],'numeric')){if(!is.na(x[1])) paste(unlist(x), collapse=';') else NA})
      sample_df$new_feature <- matched_max_feature_values

      if(length(binder_threshold)!=0){
        for(threshold in binder_threshold){
          if(is.list(threshold) & length(threshold)==2 & stringr::str_detect(new_names[[i]], threshold[[1]])){
            new_binder_feature <- paste0(new_names[i], '_', 'binders')
            sample_df$new_binder_feature <- rep(NA, nrow(sample_df))
            sample_df$new_binder_feature[which(sample_df$new_feature >= threshold[[2]])] <- 'yes'
            sample_df$new_binder_feature[which(sample_df$new_feature < threshold[[2]])] <- 'no'
            names(sample_df)[names(sample_df) == "new_binder_feature"] <- new_binder_feature
          }
        }
      }
      names(sample_df)[names(sample_df)=='new_feature'] <- new_names[i]
    }

    return(sample_df)
  }


  extract_features_from_sequence_matchings <- function(sample_df, antigen_df, sequence_column,
                            matching_type=matching.type, distance_threshold=distance.threshold, features=antigen_df_columns,
                            new_names=new_column_names, aberrant_chosen_sequences = aberrant.chosen.sequences, binder_threshold=binder.threshold, sample_id=sample.id){

    if(sample_id){
      if(!(unique(sample_df$sample_id) %in% unique(antigen_df$sample_id))){
        for(i in 1:length(new_names)){
          sample_df[new_names[[i]]] <- rep(NA, nrow(sample_df))

          if(aberrant_chosen_sequences){
            aberrant_name <- paste0(new_names[i], '_aberrant_sequences')
            sample_df[aberrant_name] <- rep(NA, nrow(sample_df))
          }

           if(length(binder_threshold)!=0){
            for(threshold in binder_threshold){
              if(is.list(threshold) & length(threshold)==2 & stringr::str_detect(new_names[[i]], threshold[[1]])){
                new_binder_feature <- paste0(new_names[i], '_', 'binders')
                sample_df$new_binder_feature <- rep(NA, nrow(sample_df))
                names(sample_df)[names(sample_df) == "new_binder_feature"] <- new_binder_feature
              }
            }
           }
         }
      return(sample_df)
    }
      antigen_df <- antigen_df[which(antigen_df$sample_id==unique(sample_df$sample_id)),]
   }


    sample_sequences <- sample_df[,sequence_column]
    if(any(stringr::str_detect(sample_sequences, ';'))){ sample_sequences <-  unlist(lapply(sample_sequences, function(x) stringr::str_split(x,';'))) }
    sample_sequences <- unique(sample_sequences)
    antigen_sequences <- antigen_df[,sequence_column]

    distance_matrix <- stringdist::stringdistmatrix(sample_sequences, antigen_sequences, method='lv')

    if(matching_type=='exact'){
      matched_indices_from_distance_matrix <- lapply(1:nrow(distance_matrix), function(x) if(min(distance_matrix[x, ]) == 0) which(distance_matrix[x, ]==min(distance_matrix[x, ])) else NA)
      }else if(matching_type=='homology'){
      matched_indices_from_distance_matrix <- lapply(1:nrow(distance_matrix), function(x) if(min(distance_matrix[x, ]) <= distance_threshold) which(distance_matrix[x, ]==min(distance_matrix[x, ])) else NA)
    }


    for(i in 1:length(features)){
      if(aberrant_chosen_sequences){
        sample_df$aberrant_chosen_sequences <- rep(NA, nrow(sample_df))
      }
      aberrant_indices <- which(sample_df$Nr_of_VDJ_chains==2 | sample_df$Nr_of_VJ_chains==2)
      matched_feature_values <- lapply(matched_indices_from_distance_matrix, function(x) if(!is.na(x[1])) unlist(antigen_df[,features[[i]]][x]) else NA)

      #matched_max_indices <- mapply(function(x,y) if(inherits(y[1])=='numeric') {if(!is.na(x[1])) x[which.max(y)] else NA} else if(inherits(y)!='numeric{'){ x[1]}, matched_indices_from_distance_matrix, matched_feature_values)
      matched_values <- sapply(matched_feature_values, function(x) if(inherits(x[1],'numeric')) {if(!is.na(x[1])) x[which.max(x)] else NA} else if(!inherits(x,'numeric')){if(!is.na(x[1])) paste0(x, collapse=';') else NA})
      class <- unique(lapply(matched_values, function(x) class(x)))

      matched_indices <- which(!is.na(matched_values))
      matched_sequences <- sample_sequences[matched_indices]
      matched_values <- matched_values[matched_indices]

      sample_df$new_feature <- rep(NA, nrow(sample_df))

      for(j in 1:length(matched_sequences)){
        sequence <- matched_sequences[j]

        matched_na_indices <- which(stringr::str_detect(sample_df[,sequence_column], sequence) & is.na(sample_df$new_feature))
        pre_matched_indices <- which(stringr::str_detect(sample_df[,sequence_column], sequence) & !is.na(sample_df$new_feature))

        sample_df$new_feature[matched_na_indices] <- matched_values[j]

        if(aberrant_chosen_sequences & length(intersect(matched_na_indices, aberrant_indices))!=0){
          sample_df$aberrant_chosen_sequences[intersect(matched_na_indices, aberrant_indices)] <- sequence
        }

        if(class=='numeric' & length(intersect(pre_matched_indices, aberrant_indices))!=0){
          if(aberrant_chosen_sequences){
            sample_df$aberrant_chosen_sequences[intersect(pre_matched_indices, aberrant_indices)][which(sample_df$new_feature[intersect(pre_matched_indices, aberrant_indices)] < matched_values[j])] <- sequence
          }
          sample_df$new_feature[intersect(pre_matched_indices, aberrant_indices)][which(sample_df$new_feature[intersect(pre_matched_indices, aberrant_indices)] < matched_values[j])] <- matched_values[j]

        }else if(class!='numeric' & length(intersect(pre_matched_indices, aberrant_indices))!=0){
          if(aberrant_chosen_sequences){
            sample_df$aberrant_chosen_sequences[intersect(pre_matched_indices, aberrant_indices)] <- paste(sample_df$aberrant_chosen_sequences[intersect(pre_matched_indices, aberrant_indices)],'/', sequence)
          }
          sample_df$new_feature[intersect(pre_matched_indices, aberrant_indices)] <- paste(sample_df$new_feature[intersect(pre_matched_indices, aberrant_indices)], '/', matched_values[j])

        }
      }

      if(length(binder_threshold)!=0){
        for(threshold in binder_threshold){
          if(is.list(threshold) & length(threshold)==2 & stringr::str_detect(new_names[[i]], threshold[[1]])){
            new_binder_feature <- paste0(new_names[i], '_', 'binders')
            sample_df$new_binder_feature <- rep(NA, nrow(sample_df))
            sample_df$new_binder_feature[which(sample_df$new_feature >= threshold[[2]])] <- 'yes'
            sample_df$new_binder_feature[which(sample_df$new_feature < threshold[[2]])] <- 'no'
            names(sample_df)[names(sample_df) == "new_binder_feature"] <- new_binder_feature
          }
        }
      }

      names(sample_df)[names(sample_df)=='new_feature'] <- new_names[i]
      if(aberrant_chosen_sequences){
        names(sample_df)[names(sample_df)=='aberrant_chosen_sequences'] <- paste0(new_names[i], '_aberrant_sequences')
      }
    }
     return(sample_df)
  }


  VDJ.matrix <- VDJ
  VDJ <- NULL

  antigen_dfs <- list()
  if(inherits(antigen.data.list[[1]],'data.frame')){
    antigen_dfs <- antigen.data.list
    antigen_names <- names(antigen.data.list)
    if(is.null(antigen_names)) {antigen_names <- paste0(rep('antigen_', length(antigen_dfs)), 1:length(antigen_dfs))}
  }else{
    antigen_dfs <- lapply(antigen.data.list, function(x) utils::read.csv(x,sep=',',header=T))
    antigen_names <- lapply(antigen.data.list, function(x) stringr::str_match(x, pattern='([A-Z0-9]+)')[1])
  }

  if(VDJ.VJ.1chain==T){
    VDJ.matrix <- VDJ.matrix[which(VDJ.matrix$Nr_of_VDJ_chains==1 & VDJ.matrix$Nr_of_VJ_chains==1),]
  }

  sample_dfs <- list()

  if(sample.id==T){
    repertoire.number <- unique(VDJ.matrix$sample_id)

    for(i in 1:length(repertoire.number)){
      sample_dfs[[i]] <- VDJ.matrix[which(VDJ.matrix$sample_id==repertoire.number[i]),]
    }

  }else{
    sample_dfs[[1]] <- VDJ.matrix
  }


  for(i in 1:length(antigen_dfs)){
    col_names <- tolower(colnames(antigen_dfs[[i]]))
    antigen_df_columns <- list()
    new_column_names <- list()

    for(j in 1:length(antigen.features)){
      antigen_df_columns[j] <- colnames(antigen_dfs[[i]])[stringr::str_which(col_names, tolower(antigen.features[[j]]))]
      new_column_names[j] <-  paste0(antigen_names[i], '_', antigen.features[j])
    }

    if(match.by=='clonotype'){
       if(is.null(antigen_dfs[[i]]$clonotype_id)) stop('Please make sure the antigen data shares the same sequence/clonotype column names as VDJ')
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_clonotype_matchings(x, antigen_dfs[[i]], clonotype_column='clonotype_id'))

    }else if(match.by=='clonotype.v3'){
       #if(is.null(antigen_dfs[[i]]$clonotype_id)) stop('Please make sure the antigen data shares the same sequence/clonotype column names as VDJ')

       names(antigen_dfs[[i]])[names(antigen_dfs[[i]])=='clonotype_id'] <- 'v3_clonotypes'
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_clonotype_matchings(x, antigen_dfs[[i]], clonotype_column='v3_clonotypes'))

    }else if(match.by=='barcode'){
       if(is.null(antigen_dfs[[i]]$barcode)) stop('Please make sure the antigen data shares the same sequence/clonotype column names as VDJ')
       if(any(is.list(antigen_dfs[[i]]$barcode))) stop('Please ensure your antigen data is in long format - avoid lists of cell barcodes per sequence rows')

       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_clonotype_matchings(x, antigen_dfs[[i]], clonotype_column='barcode'))

    }else if(match.by=='cdr3.aa'){
       for(j in 1:length(sample_dfs)){
        sample_dfs[[j]]$CDR3aa <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T,split.y=T), sample_dfs[[j]]$VDJ_cdr3s_aa, sample_dfs[[j]]$VJ_cdr3s_aa)
       }

       if(!('cdr3aa' %in% tolower(colnames(antigen_dfs[[i]])))){
        antigen_dfs[[i]]$CDR3aa <- paste0(antigen_dfs[[i]]$VDJ_cdr3s_aa, antigen_dfs[[i]]$VJ_cdr3s_aa)
        if(is.null(antigen_dfs[[i]]$CDR3aa)) stop('Please make sure the antigen data shares the same sequence/clonotype column names as VDJ')
       }

       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='CDR3aa'))

    }else if(match.by=='cdr3.nt'){
       for(j in 1:length(sample_dfs)){
         sample_dfs[[j]]$CDR3nt <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T,split.y=T), sample_dfs[[j]]$VDJ_cdr3s_nt, sample_dfs[[j]]$VJ_cdr3s_nt)
       }
       if(!('cdr3nt' %in% tolower(colnames(antigen_dfs[[i]])))){
         antigen_dfs[[i]]$CDR3nt <- paste0(antigen_dfs[[i]]$VDJ_cdr3s_nt, antigen_dfs[[i]]$VJ_cdr3s_nt)
         if(is.null(antigen_dfs[[i]]$CDR3nt)) stop('Please make sure the antigen data shares the same sequence/clonotype column names as VDJ')
       }

       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='CDR3aa'))

   }else if(match.by=='cdrh3.aa'){
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='VDJ_cdr3s_aa'))

   }else if(match.by=='cdrh3.nt'){
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='VDJ_cdr3s_nt'))

   }else if(match.by=='VDJ.VJ.aa'){
       for(j in 1:length(sample_dfs)){
         sample_dfs[[j]]$VDJ.VJ.aa <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T,split.y=T), sample_dfs[[j]]$VDJ_sequence_aa, sample_dfs[[j]]$VJ_sequence_aa)
       }

       antigen_dfs[[i]]$VDJ.VJ.aa <- paste0(antigen_dfs[[i]]$VDJ_sequence_aa, antigen_dfs[[i]]$VJ_sequence_aa)
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='VDJ.VJ.aa'))


   }else if(match.by=='VDJ.VJ.nt'){
       for(j in 1:length(sample_dfs)){
        sample_dfs[[j]]$VDJ.VJ.nt <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T,split.y=T), sample_dfs[[j]]$VDJ_sequence_nt_trimmed, sample_dfs[[j]]$VJ_sequence_nt_trimmed)
       }

       antigen_dfs[[i]]$VDJ.VJ.nt <- paste0(antigen_dfs[[i]]$VDJ_sequence_nt_trimmed, antigen_dfs[[i]]$VJ_sequence_nt_trimmed)
       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='VDJ.VJ.nt'))

   }else if(match.by=='absolut'){

       names(antigen_dfs[[i]])[names(antigen_dfs[[i]])=='CDR3'] <- 'VDJ_cdr3s_aa'
       antigen_dfs[[i]]$Energy <- abs(antigen_dfs[[i]]$Energy)

       sample_dfs <- lapply(sample_dfs, function(x) extract_features_from_sequence_matchings(x, antigen_dfs[[i]], sequence_column='VDJ_cdr3s_aa'))
   }else{
     stop('Not implemented yet')
   }
  }

  if(output.format=='dataframe.per.sample'){
    return(sample_dfs)

  }else if(output.format=='vgm'){
    if(sample.id) VDJ.matrix <- do.call('rbind',sample_dfs)
    if(!sample.id) VDJ.matrix <- sample_dfs[[1]]

    return(VDJ.matrix)
  }
}
