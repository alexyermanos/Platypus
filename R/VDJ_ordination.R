#' Performs ordination/dimensionality reduction for a species incidence matrix, depending on the species selected in the feature.columns parameter.

#' @param VDJ VDJ dataframe output from the VDJ_build function.
#' @param feature.columns Character vector. One or more column names from the VDJ to indicate the unique species for the incidence/count matrix. if more than one column is provided (e.g. c("VDJ_cdr3_aa","VJ_cdr3_aa")) these columns will be pasted together before metric calculation.
#' @param grouping.column Character. Column name of a column to group the ordination by. This could be "sample_id" to reduce across each sample. Indicative of 'sites' in a typical community data matrix/incidence matrix used in community ecology analyses (species by sites).
#' @param method Character. The ordination method; choose from either: PCA - 'pca', t-SNE - 'tsne', UMAP - 'umap', PCOA/MDS - 'mds', DCA - 'dca'.
#' @param reduction.level Character. Whether to reduce across groups ('groups'), features/sequences ('features'), or both ('both').
#' @param VDJ.VJ.1chain Boolean defaults to TRUE. Whether to filter out aberrant cells (more than 1 VDJ or VJ chain).
#' @param umap.n.neighbours Integer. Control the number of UMAP KNN neighbours when method = 'umap'.
#' @param umap.n.neighbours Integer. Control the t-SNE perplexity when method = 'tsne'.
#' @param tsne.perplexity Integrer. Defaults to 1
#' @return Returns a ggplot with the ordination analysis performer across features, groups, or both
#' @export
#' @examples
#' plot <- VDJ_ordination(VDJ = Platypus::small_vdj
#' ,feature.columns = c("VDJ_cdr3_aa"), grouping.column = "sample_id"
#' ,method = "pca", reduction.level = 'groups')
#'

VDJ_ordination <- function(VDJ,
                           feature.columns,
                           grouping.column,
                           method,
                           reduction.level,
                           VDJ.VJ.1chain,
                           umap.n.neighbours,
                           tsne.perplexity){

  if(missing(VDJ)) stop('VDJ matrix not found. Please input the VDJ/VGM[[1]] matrix!')
  if(missing(feature.columns)) feature.columns <- 'VDJ_cdr3_aa'
  if(missing(grouping.column)) grouping.column <- 'sample_id'
  if(missing(method)) method <- 'mds'
  if(missing(reduction.level)) reduction.level <- 'both'
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- TRUE
  if(missing(umap.n.neighbours)) umap.n.neighbours <- 3
  if(missing(tsne.perplexity)) tsne.perplexity <- 1


  get_abundances <- function(VDJ, feature.columns, grouping.column, VDJ.VJ.1chain){

    if(length(feature.columns) > 1){
      combine.features <- TRUE
    }else{
      combine.features <- FALSE
    }

    abundance_df <- VDJ_abundances(VDJ,
                                   feature.columns = feature.columns,
                                   proportions = 'absolute',
                                   grouping.column = grouping.column,
                                   max.groups = NULL,
                                   specific.groups = 'none',
                                   sample.column = 'none',
                                   VDJ.VJ.1chain = VDJ.VJ.1chain,
                                   treat.incomplete.groups = 'exclude',
                                   treat.incomplete.features = 'exclude',
                                   combine.features = combine.features,
                                   treat.combined.features = 'exclude',
                                   treat.combined.groups = 'exclude',
                                   specific.feature.colors = NULL,
                                   output.format = 'abundance.df')

    return(abundance_df)
  }


  abundance_to_incidence_df <- function(abundance_df){
    groups <- unique(abundance_df$group)
    groups <- groups[order(nchar(groups), groups)]

    species <- unique(abundance_df$unique_feature_values)

    incidence_matrix <- matrix(0, length(groups), length(species))
    colnames(incidence_matrix) <- species
    rownames(incidence_matrix) <- groups

    #TO DO: method to avoid loop
    for(i in 1:nrow(abundance_df)){
      incidence_matrix[abundance_df$group[i], abundance_df$unique_feature_value[i]] <- abundance_df$feature_value_counts[i]
    }

    incidence_df <- as.data.frame(incidence_matrix)
    #incidence_df$group <- rownames(incidence_df)


    return(incidence_df)
  }


  compute_ordination <- function(vdj, feature.columns, grouping.column, VDJ.VJ.1chain, method, reduction.level, umap.n.neighbours, tsne.perplexity){

    #For CRAN checks
    PC1 <- NULL
    PC2 <- NULL
    NMDS1 <- NULL
    NMDS2 <- NULL
    DCA1 <- NULL
    DCA2 <- NULL
    DIM1 <- NULL
    DIM2 <- NULL
    group <- NULL

     incidence_df <- vdj %>%
                     get_abundances(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                     abundance_to_incidence_df()

     if(method == 'pca'){

       if(reduction.level == 'groups' | reduction.level == 'both'){
         groups_ordination <- as.data.frame(stats::prcomp(t(incidence_df))$rotation[,1:2])
         groups_ordination <- groups_ordination %>% dplyr::rename(DIM1 = PC1, DIM2 = PC2)
         groups_ordination$method <- 'PCA reduction'
       }


       if(reduction.level == 'features' | reduction.level == 'both'){
         features_ordination <- as.data.frame(stats::prcomp(incidence_df)$rotation[,1:2])
         features_ordination <- features_ordination %>% dplyr::rename(DIM1 = PC1, DIM2 = PC2)
         features_ordination$method <- 'PCA reduction'

       }

     }else if(method == 'tsne'){

       if(reduction.level == 'groups' | reduction.level == 'both'){
         #requireNamespace('Rtsne')
         groups_ordination <- Rtsne::Rtsne(incidence_df, perplexity = tsne.perplexity, check_duplicates = FALSE)$Y
         groups_ordination <- as.data.frame(groups_ordination)

         colnames(groups_ordination) <- c('DIM1', 'DIM2')
         rownames(groups_ordination) <- rownames(incidence_df)
         groups_ordination$method <- 't-SNE reduction'
       }


       if(reduction.level == 'features' | reduction.level == 'both'){
         #requireNamespace('Rtsne')
         features_ordination <- Rtsne::Rtsne(t(incidence_df), perplexity = tsne.perplexity, check_duplicates = FALSE)$Y
         features_ordination <- as.data.frame(features_ordination)

         colnames(features_ordination) <- c('DIM1', 'DIM2')
         rownames(features_ordination) <- colnames(incidence_df)
         features_ordination$method <- 't-SNE reduction'

       }
     }else if(method == 'umap'){

       if(reduction.level == 'groups' | reduction.level == 'both'){
         #requireNamespace('umap')
         groups_ordination <- umap::umap(incidence_df, n_neighbors = umap.n.neighbours)$layout
         groups_ordination <- as.data.frame(groups_ordination)

         colnames(groups_ordination) <- c('DIM1', 'DIM2')
         rownames(groups_ordination) <- rownames(incidence_df)
         groups_ordination$method <- 'UMAP reduction'
       }

       if(reduction.level == 'features' | reduction.level == 'both'){
         features_ordination <- umap::umap(t(incidence_df), n_neighbors = umap.n.neighbours)$layout
         features_ordination <- as.data.frame(features_ordination)

         colnames(features_ordination) <- c('DIM1', 'DIM2')
         rownames(features_ordination) <- colnames(incidence_df)
         features_ordination$method <- 'UMAP reduction'
       }
     }else if(method == 'pcoa' | method == 'mds' | method == 'metamds'){
       invisible(utils::capture.output(mds <- vegan::metaMDS(incidence_df) %>% suppressWarnings()))

       if(reduction.level == 'groups' | reduction.level == 'both'){
         groups_ordination <- as.data.frame(vegan::scores(mds, display = 'sites'))
         groups_ordination <- groups_ordination %>% dplyr::rename(DIM1 = NMDS1, DIM2 = NMDS2)
         groups_ordination$method <- 'MDS/PCoA ordination'
       }

       if(reduction.level == 'features' | reduction.level == 'both'){
         features_ordination <- as.data.frame(vegan::scores(mds, display = 'species'))
         features_ordination <- features_ordination %>% dplyr::rename(DIM1 = NMDS1, DIM2 = NMDS2)
         features_ordination$method <- 'MDS/PCoA ordination'
       }

     }else if(method == 'dca'){
       dca <- vegan::decorana(incidence_df)

       if(reduction.level == 'groups' | reduction.level == 'both'){
         groups_ordination <- as.data.frame(vegan::scores(dca, display = 'sites'))
         groups_ordination$DCA3 <- NULL
         groups_ordination$DCA4 <- NULL
         groups_ordination <- groups_ordination %>% dplyr::rename(DIM1 = DCA1, DIM2 = DCA2)
         groups_ordination$method <- 'Detrended Correspondence Analysis (DCA)'
       }

       if(reduction.level == 'features' | reduction.level == 'both'){
         features_ordination <- as.data.frame(vegan::scores(dca, display = 'species'))
         features_ordination$DCA3 <- NULL
         features_ordination$DCA4 <- NULL
         features_ordination <- features_ordination %>% dplyr::rename(DIM1 = DCA1, DIM2 = DCA2)
         features_ordination$method <- 'Detrended Correspondence Analysis (DCA)'
       }
     }else{
       stop('Dimensionality reduction/ordination method not recognized/implemented!')
     }

     if(reduction.level == 'groups'){
       return(groups_ordination)

     }else if(reduction.level == 'features'){
       return(features_ordination)

     }else if(reduction.level == 'both'){
       out <- list(groups_ordination, features_ordination)
       names(out) <- c('groups', 'features')
       return(out)
     }
  }

  plot_ordination <- function(ordination_df, reduction.level){

    DIM1 <- NULL
    DIM2 <- NULL
    group <- NULL

    if(reduction.level == 'groups'){
      groups_ordination <- ordination_df
      groups_ordination$group <- rownames(groups_ordination)

      out_plot <- ggplot2::ggplot(data = groups_ordination) +
                  ggplot2::geom_vline(xintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_hline(yintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_point(ggplot2::aes(x = DIM1, y = DIM2, color = group), size = 6, alpha = 0.8) +
                  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"), plot.title = ggplot2::element_text(size=10)) +
                  ggplot2::labs(title = unique(groups_ordination$method), color = grouping.column)



    }else if(reduction.level == 'features'){
      features_ordination <- ordination_df

      out_plot <- ggplot2::ggplot(data = features_ordination) +
                  ggplot2::geom_vline(xintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_hline(yintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_point(ggplot2::aes(x = DIM1, y = DIM2), size = 1) +
                  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                  ggplot2::labs(title = unique(features_ordination$method))


    }else if(reduction.level == 'both'){
      features_ordination <- ordination_df$features
      groups_ordination <- ordination_df$groups
      groups_ordination$group <- rownames(groups_ordination)


      out_plot <- ggplot2::ggplot() +
                  ggplot2::geom_vline(xintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_hline(yintercept = c(0), color = "grey70", linetype = 2, size = 0.75) +
                  ggplot2::geom_point(data = features_ordination, ggplot2::aes(x = DIM1, y = DIM2), size = 1) +
                  ggplot2::geom_point(data = groups_ordination, ggplot2::aes(x = DIM1, y = DIM2, color = group), size = 6, alpha = 0.8) +
                  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                  ggplot2::labs(title = unique(features_ordination$method), color = grouping.column)
    }

    return(out_plot)
  }



  output_plot <- VDJ %>%
                 compute_ordination(feature.columns = feature.columns,
                                    grouping.column = grouping.column,
                                    VDJ.VJ.1chain = VDJ.VJ.1chain,
                                    method = method,
                                    reduction.level = reduction.level,
                                    umap.n.neighbours = umap.n.neighbours,
                                    tsne.perplexity = tsne.perplexity) %>%
                 plot_ordination(reduction.level = reduction.level)


  return(output_plot)
}
