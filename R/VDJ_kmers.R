#' Calculates and plots kmers distributions and frequencies.

#' @param VDJ VDJ dataframe output from the VDJ_GEX_matrix function.
#' @param sequence.column Character vector. One or more sequence column names from the VDJ for kmer counting. if more than one column is provided (e.g. c("VDJ_cdr3s_aa","VJ_cdr3s_aa")) these columns will be pasted together before counting the kmers.
#' @param grouping.column Character. Column name of a column to group kmer counting by. This could be "sample_id" to group each kmer by the sample.
#' @param kmer.k Integer. Length k of each kmer.
#' @param max.kmers Integer. Maximum number of kmers to be plotted in the output barplots.
#' @param specific.kmers Character vector. Specific kmers to be plotted in the output barplots.
#' @param plot.format Character. The output plot format: 'barplot' for barplots of kmer frequency per group, 'pca' for group-level PCA reduction across the kmer vectors, 'density' for kmer count density plots.
#' @param as.proportions Boolean. If TRUE, will return the kmer barplot as proportions instead of absolute counts.
#' @return Returns a ggplot with the kmer analysis depedning on the plot.format parameter
#' @export
#' @examples
#' \donttest{
#' try({
#'  VDJ_kmers(VDJ = Platypus::small_vgm[[1]],
#'  sequence.column = c("VDJ_cdr3s_aa"), grouping.column = "sample_id", kmer.k = 2, max.kmers = 5)
#'  })
#' }
#'

VDJ_kmers <- function(VDJ,
                      sequence.column,
                      grouping.column,
                      kmer.k,
                      max.kmers,
                      specific.kmers,
                      plot.format, #barplot, histogram
                      as.proportions){



  if(missing(VDJ)) stop('Please input your VDJ matrix for the kmer analysis')
  if(missing(sequence.column)) sequence.column <- 'VDJ_cdr3s_aa'
  if(missing(grouping.column)) grouping.column <- 'sample_id'
  if(missing(kmer.k)) kmer.k <- 5
  if(missing(max.kmers)) max.kmers <- 30
  if(missing(specific.kmers)) specific.kmers <- NULL
  if(missing(plot.format)) plot.format <- 'barplot' #or spectrum or pca
  if(missing(as.proportions)) as.proportions <- FALSE

  get_feature_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';', combine.sequences=FALSE){
   if(split.x==TRUE) x <- stringr::str_split(x, split.by ,simplify=TRUE)[1,]
   if(split.y==TRUE) y <- stringr::str_split(y, split.by ,simplify=TRUE)[1,]

   ccombs <- expand.grid(x,y)
   if(!combine.sequences){
     ccombs<-paste0(ccombs[,1], ' ', ccombs[,2])
   }else{
     ccombs<-paste0(ccombs[,1], ccombs[,2])

   }
   ccombs <- paste0(ccombs, collapse=collapse.by)

   return(ccombs)
  }

  if(length(sequence.column) == 2){
    VDJ[[paste0(sequence.column, collapse = ';')]] <- mapply(function(x,y) get_feature_combinations(x, y, split.x=TRUE, split.y=TRUE, combine.sequences=TRUE), VDJ[sequence.column[1]], VDJ[sequence.column[2]])
    sequence.column <- paste0(sequence.column, collapse = ';')
  }

  if(stringr::str_detect(sequence.column, '_aa')){
    aa <- TRUE
  }else if(stringr::str_detect(sequence.column, '_nt')){
    aa <- FALSE
  }

  if(length(grouping.column) > 1){
    new_name <- paste0(grouping.column, collapse = '; ')

    for(col in grouping.column){
      VDJ <- VDJ[which(!is.na(VDJ[col]) & !is.null(VDJ[col]) & VDJ[col] != ''), ]
    }

    VDJ[[new_name]] <- do.call(paste, c(VDJ[, c(grouping.column)], sep=" / "))
    grouping.column <- new_name
  }


  get_kmer_counts_per_group <- function(VDJ){
    unique_groups <- unique(VDJ[[grouping.column]])
    out_dfs <- vector(mode = 'list', length(unique_groups))

    for(i in 1:length(unique_groups)){
      sequences <- unname(unlist(VDJ[sequence.column][which(VDJ[[grouping.column]] == unique_groups[i]),]))
      #if(any(stringr::str_detect(sequences, ';'))) {
      #  sequences <- sapply(sequences, function(x) unlist(stringr::str_split(x, ';')))
      #}

      if(!aa){
        out_dfs[[i]] <- as.data.frame(colSums(kmer::kcount(ape::as.DNAbin(Biostrings::DNAStringSet(sequences)), k = kmer.k)))
      }else{
        out_dfs[[i]] <- as.data.frame(colSums(kmer::kcount(ape::as.AAbin(Biostrings::AAStringSet(sequences)), k = kmer.k)))
      }

      colnames(out_dfs[[i]]) <- 'counts'
      out_dfs[[i]]$kmers <- rownames(out_dfs[[i]])
      out_dfs[[i]]$group <- unlist(unique_groups[i])
    }

    out_dfs <- do.call('rbind', out_dfs)
    return(out_dfs)
  }

  plot_kmers <- function(kmer_df){

    #For CRAN checks
    kmers <- NULL
    counts <- NULL
    total_counts <- NULL
    group <- NULL
    PC1 <- NULL
    PC2 <- NULL

    if(plot.format == 'barplot'){
      if(is.null(specific.kmers)){
        temp_df <- kmer_df %>%
                   dplyr::group_by(kmers) %>%
                   dplyr::mutate(total_counts = sum(counts)) %>%
                   dplyr::distinct(kmers, .keep_all = TRUE) %>%
                   dplyr::arrange(dplyr::desc(total_counts))

        specific.kmers <- temp_df[1:max.kmers,]$kmers
        kmer_df <- kmer_df[kmer_df$kmers %in% specific.kmers,]

      }else{
        kmer_df <- kmer_df[kmer_df$kmers %in% specific.kmers,]
      }

      if(as.proportions){
        kmer_df <- kmer_df %>%
                   dplyr::group_by(kmers) %>%
                   dplyr::mutate(counts = (counts / sum(counts)))
      }

      plot <- ggplot2::ggplot() +
              ggplot2::geom_bar(data = kmer_df, ggplot2::aes(y = counts, x = kmers, fill = group), stat = "identity", width=0.6, color="black") +
              ggplot2::theme_bw() +
              ggplot2::theme_classic() +
              ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1))

      if(as.proportions){
        plot <- plot + ggplot2::labs(title = paste0(kmer.k, '-kmer distribution across ' , grouping.column), x = paste0(kmer.k, '-kmers'), y = paste0(kmer.k, '-kmer relative proportions'), fill = paste0(grouping.column))

      }else{
        plot <- plot + ggplot2::labs(title = paste0(kmer.k, '-kmer distribution across ' , grouping.column), x = paste0(kmer.k, '-kmers'), y = paste0(kmer.k, '-kmer frequency'), fill = paste0(grouping.column))
      }


    }else if(plot.format == 'density'){
      #requireNamespace('ggridges')
      plot <- ggplot2::ggplot() +
              ggridges::geom_density_ridges(data = kmer_df, ggplot2::aes(x = counts, y = group, fill = group), alpha = 1) +
              ggplot2::theme_bw() +
              ggplot2::theme_classic() +
              ggplot2::labs(title = paste0(kmer.k, '-kmer density for ', sequence.column, ' grouped by ',  grouping.column), x = paste0(kmer.k, '-kmer frequency'), y = paste0(grouping.column), fill = paste0(grouping.column))

    }else if(plot.format == 'pca'){
      unique_groups <- unique(kmer_df$group)
      counts_list <- vector(mode = 'list', length = length(unique_groups))

      for(i in 1:length(unique_groups)){
        counts_list[[i]] <- kmer_df$counts[kmer_df$group == unique_groups[i]]
      }
      counts_df <- do.call('rbind', counts_list)
      rownames(counts_df) <- unique_groups
      pca_out <- stats::prcomp(t(counts_df))$rotation[,1:2] %>% as.data.frame()
      pca_out$group <- rownames(pca_out)
      plot <- ggplot2::ggplot(data = pca_out) +
                  ggplot2::geom_vline(xintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_hline(yintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, fill = group), size = 6, alpha = 1, shape = 21, colour = 'black') +
                  #ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2), shape = 1, size = 6, alpha = 1, colour = 'black') +
                  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                  ggplot2::labs(title = paste0('PCA reduction of ', sequence.column, ' ', kmer.k, '-kmers'), color = paste0(grouping.column))

    }else{
      stop('Unrecognized plotting format!')
    }

    return(plot)
  }


  output_plot <- VDJ %>% get_kmer_counts_per_group() %>% plot_kmers()

  return(output_plot)
}
