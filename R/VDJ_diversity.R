#' Calculates and plots common diversity and overlap measures for repertoires and alike. Requires the vegan package

#' @param VDJ VDJ dataframe output from the VDJ_GEX_matrix function.
#' @param feature.columns Character vector. One or more column names from the VDJ of which diversity or overlap metrics are calculated. if more than one column is provided (e.g. c("VDJ_cdr3s_aa","VJ_cdr3s_aa")) these columns will be pasted together before metric calculation.
#' @param grouping.column Character. Column name of a column to group metrics by. This could be "sample_id" to calculate the metric for each sample. This column is required if metric = "simpson". If so, the simpson overlap index will be calculated pairwise for all combinations of elements in the grouping.column. Defaults to "none".
#' @param metric Character. Diversity or overlap metric to calculate. Can be c("richness", "bergerparker", "simpson", "ginisimpson", "shannon", "shannonevenness", "jaccard"). Defaults to "shannon". If jaccard is selected, a heatmap with the pairwise comparisons between all groups is returned. If any of the others is selected, a dotplot is returned
#' @param subsample.to.same.n Boolean defaults to TRUE. Whether to subsample larger groups down to the size of the smallest group
#' @param VDJ.VJ.1chain Boolean defaults to TRUE. Whether to filter out aberrant cells (more than 1 VDJ or VJ chain).
#' @return Returns a ggplot with the calculated metric for each group (if provided).
#' @export
#' @examples
#'
#' #Calculate shannon index for VDJ CDR3s by sample
#' plot <- VDJ_diversity(VDJ = Platypus::small_vgm[[1]],
#' ,feature.columns = c("VDJ_cdr3s_aa"), grouping.column = "sample_id"
#' ,metric = "shannon")
#'
#' #Calculate Gini-simpson and Simpson index for VDJ and VJ CDR3s by sample
#' VDJ_diversity(VDJ = Platypus::small_vgm[[1]],
#' ,feature.columns = c("VDJ_cdr3s_aa","VJ_cdr3s_aa"), grouping.column = "sample_id"
#' ,metric = "ginisimpson")
#'
#' #Calculate Jaccard index of J gene usage between two samples
#' VDJ_diversity(VDJ = Platypus::small_vgm[[1]],
#',feature.columns = c("VDJ_jgene"), grouping.column = "sample_id"
#',metric = "jaccard")
#'

#TO ADD: alpha-diversity: bootstrap, jackknife, qstat
#        beta: pearson, yuleq (other correlation metrics)
#        more plotting options for alpha and beta (other than barplots/heatmaps)
#        option to rarefy some alpha diversity metrics (e.g., chao1)


VDJ_diversity <- function(VDJ,
                          feature.columns,
                          grouping.column,
                          metric,
                          VDJ.VJ.1chain,
                          subsample.to.same.n){


  if(missing(VDJ)) stop('Please input your VDJ/VGM[[1]] matrix for the diversity analysis')
  if(missing(feature.columns)) feature.columns <- 'VDJ_cdr3s_aa'
  if(missing(grouping.column)) grouping.column <- 'sample_id'
  if(missing(metric)) metric <- 'richness'
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
  if(missing(subsample.to.same.n)) subsample.to.same.n <- T

  ############## UTILITY 1: internal wrapper for VDJ_abundances (get cell counts of unique features per grouping column) ##############
  get_abundances <- function(VDJ, feature.columns, grouping.column, VDJ.VJ.1chain){

    if(length(feature.columns) > 1){
      combine.features <- T
    }else{
      combine.features <- F
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


  ############## UTILITY 2: subsample groups to get to sample size = minimum size across all groups (for consistency in the diversity analyses) ##############
  subsample_abundance_df <- function(abundance_df){

    if(subsample.to.same.n){

      min_sample_n <- min(abundance_df$group_frequency)
      unique_groups <- unique(abundance_df$group)
      final_abundance_df <- vector(mode = 'list', length = length(unique_groups))

      for(i in 1:length(unique_groups)){
        subset_df <- abundance_df[abundance_df$group == unique_groups[i],]

        if(unique(subset_df$group_frequency) == min_sample_n){
          final_abundance_df[[i]] <- subset_df
          next

        }else{
          values <- subset_df$unique_feature_values
          counts <- subset_df$feature_value_counts
          sampled_values <- c(table(sample(rep(values, counts), min_sample_n)))

          final_abundance_df[[i]] <- data.frame(group = unique_groups[i], sample = unique(subset_df$sample), sample_frequency = unique(subset_df$sample_frequency),
                                                group_frequency = min_sample_n, unique_feature_values = names(sampled_values), feature_value_counts = unname(sampled_values),
                                                feature_name = unique(subset_df$feature_name))
        }
      }

      return(do.call('rbind', final_abundance_df))

    }else{
      return(abundance_df)
    }
  }


  ############## UTILITY 3: convert an abundance dataframe to an incidence one (groups x all unique features across all groups - zero if a group misses a feature, otherwise the number of cells for a specific feature of that specific group) ##############
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


  ############## ALPHA DIVERSITY - computes the following alpha diversity metrics: richness, bergerparker, simpson, ginisimpson, inversesimpson, shannon, chao1, ace ##############
  compute_alpha_diversity <- function(vdj, feature.columns, grouping.column, VDJ.VJ.1chain){

    #missing variable definitions
    group <- NULL
    unique_feature_values <- NULL
    feature_value_counts <- NULL

    abundance_df <- vdj %>%
                    get_abundances(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                    subsample_abundance_df()


    if(metric == 'richness'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = length(unique_feature_values))

      metric_df$metric_name <- 'Species richness'

    }else if(metric == 'bergerparker'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = max(feature_value_counts) / sum(feature_value_counts))

      metric_df$metric_name <- 'Berger-Parker index'

    }else if(metric == 'simpson'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = vegan::diversity(feature_value_counts, 'simpson'))

      metric_df$metric_name <- 'Simpson index'

    }else if(metric == 'ginisimpson'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = 1 - vegan::diversity(feature_value_counts, 'simpson'))

      metric_df$metric_name <- 'Gini-Simpson index'

    }else if(metric == 'inversesimpson'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = vegan::diversity(feature_value_counts, 'invsimpson'))

      metric_df$metric_name <- 'Inverse Simpson index'

    }else if(metric == 'shannon'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = vegan::diversity(feature_value_counts, 'shannon'))

      metric_df$metric_name <- 'Shannon index'

    }else if(metric == 'chao1'){

      metric_df <- abundance_df %>%
                   abundance_to_incidence_df() %>%
                   vegan::estimateR()

      metric_df <- data.frame(group = colnames(metric_df), metric = round(metric_df['S.chao1',],2), se = round(metric_df['se.chao1',],2))
      metric_df$metric_name <- 'Chao1 index'

    }else if(metric == 'ace'){

      metric_df <- abundance_df %>%
                   abundance_to_incidence_df() %>%
                   vegan::estimateR()

      metric_df <- data.frame(group = colnames(metric_df), metric = round(metric_df['S.ACE',],2), se = round(metric_df['se.ACE',],2))
      metric_df$metric_name <- 'ACE index'

    }else if(metric == 'coverage'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = length(which(feature_value_counts == 1)) / sum(feature_value_counts))

      metric_df$metric_name <- "Good's coverage"

    }else{
      stop('Diversity metric was not found/implemented!')
    }

    metric_df$colors <- grDevices::rainbow(length(metric_df$group))

    return(metric_df)
  }


  ############## SPECIES EVENNESS - computes the following evenness metrics: shannonevenness, simpsonevenness, bulla, camargo, smithwilson, heip  ##############
  compute_evenness <- function(vdj, feature.columns, grouping.column, VDJ.VJ.1chain){

    #missing variable definitions
    group <- NULL
    unique_feature_values <- NULL
    feature_value_counts <- NULL

    abundance_df <- vdj %>%
                    get_abundances(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                    subsample_abundance_df()


    if(metric == 'shannonevenness'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = vegan::diversity(feature_value_counts, 'shannon') / log(length(unique_feature_values)))

      metric_df$metric_name <- 'Shannon evenness'

    }else if(metric == 'simpsonevenness' | metric == 'pielou'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = 1/(vegan::diversity(feature_value_counts, 'simpson') * length(unique_feature_values)))

      metric_df$metric_name <- 'Simpson evenness'


    }else if(metric == 'bulla'){
      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = bulla_evenness(feature_value_counts))

      metric_df$metric_name <- 'Bulla evenness'


    }else if(metric == 'camargo'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = camargo_evenness(feature_value_counts))

      metric_df$metric_name <- 'Camargo evenness'


    }else if(metric == 'smithwilson'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = smithwilson_evenness(feature_value_counts))

      metric_df$metric_name <- 'Smith-Wilson evenness'

    }else if(metric == 'heip'){

      metric_df <- abundance_df %>%
                   dplyr::group_by(group) %>%
                   dplyr::summarise(metric = (exp(vegan::diversity(feature_value_counts, 'shannon')) - 1) / (length(unique_feature_values) - 1))

      metric_df$metric_name <- 'Heip evenness'


    }else{
      stop('Evenness metric was not found/implemented!')
    }

    metric_df$colors <- grDevices::rainbow(length(metric_df$group))

    return(metric_df)
  }



  ############## EVENESS METRICS ##############
  bulla_evenness <- function(abundance_vector, ignore.zeros = T){

    #Modified from https://rdrr.io/github/microbiome/microbiome/src/R/evenness.R
    if(ignore.zeros){
      abundance_vector <- abundance_vector[abundance_vector > 0]
    }

    S <- sum(abundance_vector > 0)
    rel_abundance <- abundance_vector / sum(abundance_vector)
    O <- sum(pmin(rel_abundance, 1/S))

    evenness <- (O - 1/S) / (1 - 1/S)

    return(evenness)

  }

  camargo_evenness <- function(abundance_vector, ignore.zeros = T){

    #Modified from https://rdrr.io/github/microbiome/microbiome/src/R/evenness.R
    if(ignore.zeros){
      abundance_vector <- abundance_vector[abundance_vector > 0]
    }

    S <- sum(abundance_vector > 0)

    out <- 0
    for(i in 1:(S - 1)){
      out <- out + sum(abs(abundance_vector[(i+1):S] - abundance_vector[i]))
    }

    evenness <- 1 - out/(sum(abundance_vector) * S)

    return(evenness)
  }

  smithwilson_evenness <- function(abundance_vector, ignore.zeros = T){

    #missing variable definitions
    x <- NULL

    #Modified from https://rdrr.io/github/microbiome/microbiome/src/R/evenness.R
    if(ignore.zeros){
      abundance_vector <- abundance_vector[abundance_vector > 0]
    }

    total_abundance <- sum(abundance_vector)
    log_abundance <- log(abundance_vector)
    log_abundance[abundance_vector == 0] <- 0

    S <- sum(abundance_vector > 0)

    a <- log_abundance / S
    b <- sum(x) #!!!
    c <- (log_abundance - b) ^ 2 / S
    c[abundance_vector == 0 ] <- 0
    d <- sum(c)

    evenness <- (1 - 2/pi * atan(d))

    return(evenness)
  }


  ############## PLOTTING FUNCTION FOR ALPHA DIVERSITY AND EVENNESS METRICS - barplots for the specified alpha diversity/evenness metric ##############
  plot_alpha_diversity <- function(metric_df){

    #missing variable definitions
    group <- NULL
    colors <- NULL
    se <- NULL

    title <- unique(metric_df$metric_name)
    plot_out <- ggplot2::ggplot(metric_df, ggplot2::aes(x = group, y = metric, fill = colors)) +
                ggplot2::geom_bar(show.legend = F, stat = "identity", width=0.6, color="black") +
                ggplot2::labs(title = title, x = "", y = title) +
                ggplot2::theme(panel.background = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "none")

    if('se' %in% colnames(metric_df)){
      plot_out <- plot_out +
                  ggplot2::geom_errorbar(ggplot2::aes(ymin = metric - se, ymax = metric + se), width= 0.1, size = 0.6)

    }

    return(plot_out)
  }


  ############## BETA DIVERSITY - computes the following beta diversity metrics: anderberg, braycurtis, cosine, euclidean, hamming, jaccard, manhattan, morisitahorn, sorensen ##############
  compute_beta_diversity <- function(vdj, feature.columns, grouping.column, VDJ.VJ.1chain){

    abundance_df <- vdj %>%
                    get_abundances(feature.columns, grouping.column, VDJ.VJ.1chain)

    if(metric == 'anderberg'){
      metric_df <- anderberg_similarity(abundance_df)

    }else if(metric == 'braycurtis'){
      metric_df <- braycurtis_index(abundance_df)


    }else if(metric == 'cosine'){
      metric_df <- cosine_similarity(abundance_df)


    }else if(metric == 'euclidean'){
      metric_df <- euclidean_distance(abundance_df)


    }else if(metric == 'hamming'){
      metric_df <- hamming_distance(abundance_df)


    }else if(metric == 'jaccard'){
      metric_df <- jaccard_similarity(abundance_df)


    }else if(metric == 'manhattan'){
      metric_df <- manhattan_distance(abundance_df)


    }else if(metric == 'morisitahorn'){
      metric_df <- morisitahorn_index(abundance_df)


    }else if(metric == 'sorensen'){
      metric_df <- sorensen_similarity(abundance_df)

    }else{
      stop('Beta diversity metric was not found/implemented!')
    }

    return(metric_df)
  }


  ############## BETA DIVERSITY METRICS ##############
  anderberg_similarity <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Anderberg similarity'

    for(i in 1:nrow(combs)){
      a <- df$unique_feature_values[df$group == combs[i,1]]
      b <- df$unique_feature_values[df$group == combs[i,2]]
      a_b <- intersect(a, b)

      combs$metric[i] <- 1 - ( length(a_b) / (2 * length(a) + 2 * length(b) - 3 * length(a_b)) )
    }

    return(combs)
  }


  braycurtis_index <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Bray-Curtis index'

    df <- abundance_to_incidence_df(df)

    for(i in 1:nrow(combs)){
      a <- as.numeric(df[combs[i,1],])
      b <- as.numeric(df[combs[i,2],])

      m <- sum(a)
      n <- sum(b)

      c <- sapply(1:length(a), function(i) min(c(a[i], b[i])))
      c <- sum(c)

      combs$metric[i] <- 1 - 2 * (c / (m + n))
    }

    return(combs)
  }

  cosine_similarity <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Cosine similarity'

    df <- abundance_to_incidence_df(df)

    for(i in 1:nrow(combs)){
      a <- as.numeric(df[combs[i,1],])
      b <- as.numeric(df[combs[i,2],])

      combs$metric[i] <- (a %*% b) / (norm(a, type = '2') * norm(b, type = '2'))
    }

    return(combs)

  }

  euclidean_distance <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Euclidean distance'

    df <- abundance_to_incidence_df(df)

    for(i in 1:nrow(combs)){
      a <- as.numeric(df[combs[i,1],])
      b <- as.numeric(df[combs[i,2],])

      dist <- mapply(function(x,y) (x - y) ^ 2, a, b)

      combs$metric[i] <- sqrt(sum(dist))
    }

    return(combs)
  }


  hamming_distance <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Hamming distance'

    for(i in 1:nrow(combs)){
      a <- df$unique_feature_values[df$group == combs[i,1]]
      b <- df$unique_feature_values[df$group == combs[i,2]]
      a_b <- intersect(a, b)

      combs$metric[i] <- length(a) + length(b) - 2 * length(a_b)
    }

    return(combs)
  }


  jaccard_similarity <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Jaccard similarity'

    for(i in 1:nrow(combs)){
      a <- df$unique_feature_values[df$group == combs[i,1]]
      b <- df$unique_feature_values[df$group == combs[i,2]]

      intersection <- length(intersect(a, b))
      union <- length(a) + length(b) - intersection
      combs$metric[i] <- intersection / union
    }

    return(combs)

  }

  manhattan_distance <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Manhattan distance'

    df <- abundance_to_incidence_df(df)

    for(i in 1:nrow(combs)){
      a <- as.numeric(df[combs[i,1],])
      b <- as.numeric(df[combs[i,2],])

      combs$metric[i] <- sum(abs(a-b))
    }

    return(combs)

  }

  morisitahorn_index <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Morisita-Horn index'

    df <- abundance_to_incidence_df(df)

    for(i in 1:nrow(combs)){
      a <- as.numeric(df[combs[i,1],])
      b <- as.numeric(df[combs[i,2],])

      m <- sum(a)
      n <- sum(b)

      c <- sapply(1:length(a), function(i) (a[i] / m) * (b[i] / n))
      c <- sum(c)
      d <- sapply(1:length(a), function(i) (a[i] / m) ^ 2)
      d <- sum(d)
      e <- sapply(1:length(a), function(i) (b[i] / n) ^ 2)
      e <- sum(e)

      combs$metric[i] <- 1 - (2 * (c / (d + e)))
    }

    return(combs)
  }

  sorensen_similarity <- function(df){

    groups <- unique(df$group)
    combs <- data.frame(t(utils::combn(groups, m = 2, simplify = TRUE)))

    combs$metric <- NA
    combs$metric_name <- 'Sørensen-Dice similarity'

    for(i in 1:nrow(combs)){
      a <- df$unique_feature_values[df$group == combs[i,1]]
      b <- df$unique_feature_values[df$group == combs[i,2]]
      a_b <- intersect(a, b)

      combs$metric[i] <- (2 * length(a_b)) / (length(a) + length(b))
    }

    return(combs)
  }


  ############## PLOTTING FUNCITON FOR BETA DIVERSITY - plots beta diversity/overlap/dissimilariy heatmaps ##############
  plot_beta_diversity <- function(metric_df){

    title_out <- unique(metric_df$metric_name)
    metric_df$metric <- as.numeric(round(metric_df$metric, 3))

    plot_out <- ggplot2::ggplot(metric_df, ggplot2::aes(x = metric_df[,2], y = metric_df[,1], fill = as.numeric(metric))) +
                ggplot2::geom_tile() +
                ggplot2::geom_text(ggplot2::aes(label = metric), size = 4)+
                ggplot2::theme(panel.background = ggplot2::element_blank(),
                            axis.text = ggplot2::element_text(size = 30), axis.line.x = ggplot2::element_blank(),axis.line.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                            text = ggplot2::element_text(size=30), legend.key = ggplot2::element_rect(colour = "white"), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, size = 25),
                            plot.subtitle = ggplot2::element_text(size = 15),axis.text.x = ggplot2::element_text(angle = 60,vjust = 1, hjust=1, size = 12), axis.text.y = ggplot2::element_text(size = 12)) +
                ggplot2::labs(title = title_out, x = "", y = "", fill = "") +
                ggplot2::scale_fill_viridis_c()

    return(plot_out)
  }



  ############## MAIN FUNCTION CALL - categorize the metrics, call the respective subroutines and their plotting function, output plots
  if(metric %in% c('richness', 'bergerparker', 'simpson', 'ginisimpson', 'inversesimpson', 'shannon', 'chao1', 'ace')){
    output_plot <- VDJ %>%
                   compute_alpha_diversity(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                   plot_alpha_diversity()


  }else if(metric %in% c('shannonevenness', 'simpsonevenness', 'bulla', 'camargo', 'smithwilson', 'heip')){
    output_plot <- VDJ %>%
                   compute_evenness(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                   plot_alpha_diversity()


  }else if(metric %in% c('anderberg', 'braycurtis', 'cosine', 'euclidean', 'hamming', 'jaccard', 'manhattan', 'morisitahorn', 'sorensen')){
    output_plot <- VDJ %>%
                   compute_beta_diversity(feature.columns, grouping.column, VDJ.VJ.1chain) %>%
                   plot_beta_diversity()


  }else{
    stop('Diversity metric not recognized/implemented!')
  }


  return(output_plot)
}
