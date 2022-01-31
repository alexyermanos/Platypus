#' Creates phylogenetic trees from a VDJ dataframe


#'@description Creates phylogenetic trees as tidytree dataframes from an input VDJ dataframe. The resulting phylogenetic trees can be plotted using VDJ_phylogenetic_trees_plot. Both of these functions require the tidytree and ggtree packages.
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param sequence.type string - sequences which will be used when creating the phylogenetic trees. 'cdr3' for CDR3s of both VDJs and VJs, 'cdrh3' for VDJ CDR3s, 'VDJ.VJ' for pasted full sequences of both VDJ and VJ, 'VDJ' for full VDJ sequences, 'VJ' for full VJ.
#' @param as.nucleotide boolean - if T, will only consider the DNA sequences specified by sequence.type, else it will consider the amino acid ones.
#' @param trimmed boolean - in the case of full VDJ or VJ nt sequences, if the trimmed sequences should be consider (trimmed=T), or raw ones. You need to call MIXCR first on the VDJ dataframe using VDJ_call_MIXCR().
#' @param include.germline boolean - if T, a germline sequence will be included in the trees (root), obtained by pasting the VDJ_trimmed_ref and VJ_trimmed_ref sequences. You need to call MIXCR first on the VDJ dataframe using VDJ_call_MIXCR().
#' @param global.clonotype boolean - if T, will ignore samples from the sample_id column, creating global clonotypes.
#' @param VDJ.VJ.1chain boolean - if T, will remove aberrant cells from the VDJ matrix.
#' @param additional.feature.columns list of strings or NULL - VDJ column names which will comprise the per-sequence features to be included in the tidytree dataframe, which will be used to label nodes/ determines their color/ size etc. See also the VDJ_phylogenetic_trees_plot function.
#' @param filter.na.columns list of strings - VDJ columns names: if a phylogenetic tree/tidytree dataframe has all elements = NA in that feature, that tree will be completely removed.
#' @param maximum.lineages integer or 'all' - maximum number of clonotypes to create trees for. If 'all', will create trees for all clonotypes.
#' @param minimum.sequences integer - lower bound of sequences for a tree. Defaults to 3. Trees with a lower number will be automatically removed.
#' @param maximum.sequences integer - upper bound of sequences for a tree. Additional sequences will be removed, after being ordered by their total frequency.
#' @param tree.algorithm string - the algorithm used when constructing the phylogenetic trees. 'nj' for Neighbour-Joining, 'bionj', 'fastme.bal', and 'fastme.ols'
#' @param tree.level string - level at which to build phylogenetic trees. 'intraclonal' - tree per clonotype, per sample, 'global.clonotype' - global clonotype trees (include.germline must be F), irrespective of sample, 'combine.first.trees' will combine the trees for the most expanded clonotypes, per sample (include.germline must be F).
#' @param no.trees.combined integer - number of trees to combine if tree.level='combine.first.trees'.
#' @param germline.scale.factor numeric - as germlines are incredibly distant from their closest neighbours (in the tree), this controls the scale factor for the germline tree branch length for more intelligible downstream plotting.
#' @param output.format string - 'tree.df.list' returns a nested list of tidytree dataframes, per clonotype and per sample; 'lineage.df.list' returns a list of lineage dataframes - unique sequences per clonotype,
#' @param parallel string - parallelization method to be used to accelerate computations, 'none', 'mclapply', or 'parlapply'.

#' @return Nested list of tidytree dataframes or lineage dataframes.
#' @export
#' @examples
#' \dontrun{
#' VDJ_phylogenetic_trees(VDJ=VDJ, sequence.type='VDJ.VJ',
#' trimmed=TRUE, as.nucleotide=TRUE, include.germline=TRUE,
#' additional.feature.columns=NULL, tree.level='intraclonal',
#' output.format='tree.df.list')
#'}


VDJ_phylogenetic_trees <- function(VDJ,
                                   sequence.type,
                                   as.nucleotide,
                                   trimmed,
                                   include.germline,
                                   global.clonotype,
                                   VDJ.VJ.1chain,
                                   additional.feature.columns,
                                   filter.na.columns,
                                   maximum.lineages,
                                   minimum.sequences,
                                   maximum.sequences,
                                   tree.algorithm,
                                   tree.level,
                                   no.trees.combined,
                                   germline.scale.factor,
                                   output.format,
                                   parallel){

 if(missing(output.format)) output.format <- 'tree.df.list'
 if(missing(VDJ)) stop('Please input your data as VDJ')
 if(missing(sequence.type)) sequence.type <- 'VDJ.VJ'
 if(missing(as.nucleotide)) as.nucleotide <- T
 if(missing(trimmed)) trimmed <- T
 if(missing(include.germline)) include.germline <- T
 if(missing(global.clonotype)) global.clonotype <- F
 if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
 if(missing(additional.feature.columns)) additional.feature.columns <- NULL
 if(missing(filter.na.columns)) filter.na.columns <- NULL
 if(missing(maximum.lineages)) maximum.lineages <- 50
 if(missing(minimum.sequences)) minimum.sequences <- 4
 if(missing(maximum.sequences)) maximum.sequences <- 20
 if(missing(tree.algorithm) & output.format=='tree.df.list') tree.algorithm <- 'nj'
 if(missing(tree.level)) tree.level <- 'intraclonal'
 if(missing(germline.scale.factor)) germline.scale.factor <- 0.0001
 if(missing(no.trees.combined)) no.trees.combined <- 5
 if(missing(parallel)) parallel <- 'none'

 options(warn=-1)

  get_sequence_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';'){
    if(split.x==T) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
    if(split.y==T) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

    ccombs <- expand.grid(x,y)
    ccombs<-paste0(ccombs[,1], ccombs[,2])
    ccombs <- paste0(ccombs, collapse=collapse.by)

    return(ccombs)
  }


  transform_clonotype_to_lineage_df <- function(clonotype_df){


    if(sequence.type=='cdr3'){
      if(as.nucleotide==T){
        combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_cdr3s_nt, clonotype_df$VJ_cdr3s_nt)
        clonotype_df$lineage_sequences <- combined_sequences
      }else{
        combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_cdr3s_aa, clonotype_df$VJ_cdr3s_aa)
        clonotype_df$lineage_sequences <- combined_sequences
      }
    }else if(sequence.type=='cdrh3'){
      if(as.nucleotide==T){
        clonotype_df$lineage_sequences <- clonotype_df$VDJ_cdr3s_nt
      }else if(as.nucleotide==F){
        clonotype_df$lineage_sequences <- clonotype_df$VDJ_cdr3s_aa
      }
    }else if(sequence.type=='VDJ.VJ'){
      if(as.nucleotide==T){
        if(trimmed==T){
          combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_sequence_nt_trimmed, clonotype_df$VJ_sequence_nt_trimmed)
          if(length(combined_sequences)==0) stop('Make sure you used trim.and.align=T to obtain your clonotype_df.')
        }else{
          combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_sequence_nt_raw, clonotype_df$VJ_sequence_nt_raw)
        }
        clonotype_df$lineage_sequences <- combined_sequences

      }else if(as.nucleotide==F){
        if(trimmed==T){
          combined_sequences <- mapply(function(x,y) get_sequence_combinations(x,y,split.x=T, split.y=T), clonotype_df$VDJ_sequence_aa, clonotype_df$VJ_sequence_aa)
          if(length(combined_sequences)==0) stop('Make sure you used trim.and.align=T to obtain your clonotype_df.')
        }else{
          stop('Not available for raw aa sequences.')
        }
        clonotype_df$lineage_sequences <- combined_sequences
      }
    }else if(sequence.type!='VDJ.VJ' & sequence.type!='cdr3' & sequence.type!='cdrh3'){
      if(as.nucleotide==T){
        if(trimmed==T){
          column_name <- paste0(sequence.type, '_sequence_nt_trimmed')
        }else{
          column_name <- paste0(sequence.type, '_sequence_nt_raw')
        }
        if(!(column_name %in% colnames(clonotype_df))) stop('Make sure you used trim.and.align=T to obtain your clonotype_df.')
        clonotype_df$lineage_sequences <-  clonotype_df$column_name
     }else{
        if(trimmed==T){
          column_name <- paste0(sequence.type, '_sequence_aa')
        }else{
          stop('Not available for raw aa sequences.')
        }
        if(!(column_name %in% colnames(clonotype_df))) stop('Make sure you used trim.and.align=T to obtain your clonotype_df.')
        clonotype_df$lineage_sequences <-  clonotype_df$column_name
      }
    }


    all_sequences <- unlist(lapply(clonotype_df$lineage_sequences, function(x) unlist(stringr::str_split(x, ';'))))
    unique_sequences <- unlist(unique(all_sequences))
    sequence_frequency <- unlist(lapply(unique_sequences, function(x) length(all_sequences[all_sequences==x])))
    unique_sequences <- unique_sequences[order(sequence_frequency, decreasing=T)]
    sequence_frequency <- sequence_frequency[order(sequence_frequency, decreasing=T)]
    cell_barcodes <- lapply(unique_sequences, function(x) clonotype_df$barcode[which(stringr::str_detect(clonotype_df$lineage_sequences, x))])

    lineage_df <- data.frame(lineage_sequences=unique_sequences, sequence_frequency=sequence_frequency, barcodes=matrix(cell_barcodes))


    features_to_select <- c('clonotype_id', 'clonotype_frequency', 'sample_id', unlist(additional.feature.columns))
    for(feature in features_to_select){
      feature_list <- list()
      for(i in 1:length(cell_barcodes)){
        feature_list[[i]] <- unique(unlist(lapply(cell_barcodes[[i]], function(x)clonotype_df[feature][which(clonotype_df$barcode==x),])))[1]
      }
      lineage_df[feature] <- list(feature_list)
    }

    if(nrow(lineage_df)<minimum.sequences) {
      return(NULL)
      #break
    }

    for(col in filter.na.columns){
      if(all(is.na(unlist(lineage_df[,col])))  | all(unlist(lineage_df[,col])=='')) {
        return(NULL)
      }
    }

    if(!is.null(maximum.sequences)){
      if(nrow(lineage_df)>maximum.sequences) lineage_df <- lineage_df[1:maximum.sequences,]
    }

    if(include.germline==T){
      if(sequence.type=='VJ' & as.nucleotide==T & trimmed==T){
        if(!('VJ_trimmed_ref' %in% colnames(clonotype_df))) stop('Please run VDJ.GEX.matrix with trim.and.align=T first.')
        germline_seq <- unique(clonotype_df$VJ_trimmed_ref)[1]
      }else if(sequence.type=='VDJ' & as.nucleotide==T & trimmed==T){
        if(!('VDJ_trimmed_ref' %in% colnames(clonotype_df))) stop('Please run VDJ.GEX.matrix with trim.and.align=T first.')
        germline_seq <- unique(clonotype_df$VDJ_trimmed_ref)[1]
      }else if(sequence.type=='VDJ.VJ' & as.nucleotide==T & trimmed==T){
        if(!('VDJ_trimmed_ref' %in% colnames(clonotype_df))) stop('Please run VDJ.GEX.matrix with trim.and.align=T first.')
        VDJ <- unique(clonotype_df$VDJ_trimmed_ref)
        VJ <- unique(clonotype_df$VDJ_trimmed_ref)
        germline_seq <- paste0(VDJ,VJ)[1]
      }else{
        germline_seq <- 'Not found'
      }

      lineage_df$germline <- rep('no', nrow(lineage_df))
      lineage_df <- rbind(lineage_df, NA)
      lineage_df$lineage_sequences[nrow(lineage_df)] <- germline_seq
      lineage_df$germline[which(lineage_df$lineage_sequences==germline_seq)] <- 'yes'
    }

    lineage_df$sequence_id <- as.character(c(1:nrow(lineage_df)))

    return(lineage_df)
  }


  transform_lineage_to_tree_df <- function(lineage_df){

    if(tree.algorithm == 'nj'){
      output_tree <- tidytree::as_tibble(ape::nj(stringdist::stringdistmatrix(lineage_df$lineage_sequences, lineage_df$lineage_sequences)))
      if(include.germline==T){
        output_tree <- phytools::reroot(ape::as.phylo(output_tree), node.number = output_tree$node[which(output_tree$label == nrow(lineage_df))])
        output_tree <- tidytree::as_tibble(output_tree)
      }
    }

    if(tree.algorithm == 'bionj'){
      output_tree <- tidytree::as_tibble(ape::bionj(stringdist::stringdistmatrix(lineage_df$lineage_sequences, lineage_df$lineage_sequences)))
      if(include.germline==T){
        output_tree <- phytools::reroot(ape::as.phylo(output_tree), node.number = output_tree$node[which(output_tree$label == nrow(lineage_df))])
        output_tree <- tidytree::as_tibble(output_tree)
      }
    }

    if(tree.algorithm == 'fastme.bal'){
      output_tree <- tidytree::as_tibble(ape::fastme.bal(stringdist::stringdistmatrix(lineage_df$lineage_sequences, lineage_df$lineage_sequences)))
      if(include.germline==T){
        output_tree <- phytools::reroot(ape::as.phylo(output_tree), node.number = output_tree$node[which(output_tree$label == nrow(lineage_df))])
        output_tree <- tidytree::as_tibble(output_tree)

      }
    }

    if(tree.algorithm == 'fastme.ols'){
      output_tree <- tidytree::as_tibble(ape::fastme.ols(stringdist::stringdistmatrix(lineage_df$lineage_sequences, lineage_df$lineage_sequences)))
      if(include.germline==T){
        output_tree <- phytools::reroot(ape::as.phylo(output_tree), node.number = output_tree$node[which(output_tree$label == nrow(lineage_df))])
        output_tree <- tidytree::as_tibble(output_tree)
      }
    }

    output_tree <- dplyr::full_join(output_tree, lineage_df, by=c('label'='sequence_id'))
    if(include.germline & !(is.null(germline.scale.factor))){
      output_tree$branch.length[which(output_tree$germline=='yes')] <- germline.scale.factor*output_tree$branch.length[which(output_tree$germline=='yes')]
      output_tree$branch.length[which(output_tree$node==nrow(output_tree))] <- germline.scale.factor*output_tree$branch.length[which(output_tree$node==nrow(output_tree))]

    }

    return(output_tree)
  }


  if(class(VDJ)=='data.frame'){
    VDJ.GEX.matrix <- list()
    VDJ.GEX.matrix[[1]] <- VDJ
    VDJ <- NULL

    sample_dfs <- list()
    lineage_dfs <- list()

    if(VDJ.VJ.1chain){
      VDJ <- VDJ[which(VDJ$Nr_of_VDJ_chains==1 & VDJ$Nr_of_VJ_chains==1),]
    }

    if(global.clonotype==F){
      repertoire.number <- unique(VDJ.GEX.matrix[[1]]$sample_id)

      for(i in 1:length(repertoire.number)){
        sample_dfs[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id==repertoire.number[i]),]
        sample_dfs[[i]] <- sample_dfs[[i]]
      }
    }else{
      sample_dfs[[1]] <-VDJ.GEX.matrix[[1]]
    }

    for(i in 1:length(sample_dfs)){
      lineage_dfs[[i]] <- list()
      sample_dfs[[i]] <- sample_dfs[[i]][order(sample_dfs[[i]]$clonotype_frequency, decreasing=T), ]
      clonotype_dfs <-  split(sample_dfs[[i]], factor(sample_dfs[[i]]$clonotype_id, levels=unique(sample_dfs[[i]]$clonotype_id)))

      if(maximum.lineages!='all') clonotype_dfs <- clonotype_dfs[1:maximum.lineages]

      if(parallel=='mclapply'){

        #require(parallel)
        cores <- parallel::detectCores() - 1
        lineage_dfs[[i]] <- parallel::mclapply(clonotype_dfs, transform_clonotype_to_lineage_df, mc.cores=cores)

      }else if(parallel=='parlapply'){

        #require(doParallel)
        no_cores <- parallel::detectCores() - 1
        doParallel::registerDoParallel(cores=no_cores)
        cl <- parallel::makeCluster(no_cores, type="FORK")
        lineage_dfs[[i]] <- parallel::parLapply(cl, clonotype_dfs, transform_clonotype_to_lineage_df)
        parallel::stopCluster(cl)

      }else{

        lineage_dfs[[i]] <- lapply(clonotype_dfs, function(x) transform_clonotype_to_lineage_df(x))
      }
    }

  }else if(class(VDJ)=='list'){
    clonotype_dfs <- VDJ
    VDJ <- NULL
    for(i in 1:length(clonotype_dfs)){
      lineage_dfs <- list()

      if(maximum.lineages!='all') clonotype_dfs[[i]] <- clonotype_dfs[[i]][1:maximum.lineages]

      if(parallel=='mclapply'){

        #require(parallel)
        cores <- parallel::detectCores() - 1
        lineage_dfs[[i]] <- parallel::mclapply(clonotype_dfs[[i]], transform_clonotype_to_lineage_df, mc.cores=cores)

      }else if(parallel=='parlapply'){

        #require(doParallel)
        no_cores <- parallel::detectCores() - 1
        doParallel::registerDoParallel(cores=no_cores)
        cl <- parallel::makeCluster(no_cores, type="FORK")
        lineage_dfs[[i]] <- parallel::parLapply(cl, clonotype_dfs[[i]], transform_clonotype_to_lineage_df)
        parallel::stopCluster(cl)

      }else{

        lineage_dfs[[i]] <- lapply(clonotype_dfs[[i]], function(x) transform_clonotype_to_lineage_df(x))
      }
    }

  }else stop('Incompatible input data')


  lineage_dfs <- lapply(lineage_dfs, function(x) x[!sapply(x,is.null)])
  lineage_dfs <- lapply(lineage_dfs, function(x) unname(x))


  if(output.format=='lineage.df.list'){
    return(lineage_dfs)

  }else if(output.format=='tree.df.list'){
    #require('tidytree')
    #require('ggtree')
    if(tree.level=='intraclonal'){
      output_trees <- list()
      for(i in 1:length(lineage_dfs)){
        output_trees[[i]] <- list()

        if(parallel=='mclapply'){

          #require(parallel)
          cores <- parallel::detectCores() - 1
          output_trees[[i]] <- parallel::mclapply(lineage_dfs[[i]], transform_lineage_to_tree_df, mc.cores=cores)

        }else if(parallel=='parlapply'){

          #require(doParallel)
          no_cores <- parallel::detectCores() - 1
          doParallel::registerDoParallel(cores=no_cores)
          cl <- parallel::makeCluster(no_cores, type="FORK")
          output_trees[[i]] <- parallel::parLapply(cl, lineage_dfs[[i]], transform_lineage_to_tree_df)
          parallel::stopCluster(cl)

        }else{

          output_trees[[i]] <- lapply(lineage_dfs[[i]], function(x) transform_lineage_to_tree_df(x))
        }
      }

      return(output_trees)
    }else if(tree.level=='global.clonotype' & include.germline==F){

      temp_dfs <- list()

      for(i in length(lineage_dfs)){
        temp_dfs[[i]] <- do.call('rbind', lineage_dfs[[i]])
      }

      combined_df <- do.call('rbind', temp_dfs)

      lineage_dfs <-  split(combined_df, combined_df$clonotype_id)

      output_trees <- list()

      for(i in 1:length(lineage_dfs)){
        output_trees[[i]] <- list()

        if(parallel=='mclapply'){

          #require(parallel)
          cores <- parallel::detectCores() - 1
          output_trees[[i]] <- parallel::mclapply(lineage_dfs[[i]], transform_lineage_to_tree_df, mc.cores=cores)

        }else if(parallel=='parlapply'){

          #require(doParallel)
          no_cores <- parallel::detectCores() - 1
          doParallel::registerDoParallel(cores=no_cores)
          cl <- parallel::makeCluster(no_cores, type="FORK")
          output_trees[[i]] <- parallel::parLapply(cl, lineage_dfs[[i]], transform_lineage_to_tree_df)
          parallel::stopCluster(cl)

        }else{

          output_trees[[i]] <- lapply(lineage_dfs[[i]], function(x) transform_lineage_to_tree_df(x))
        }
      }

    }else if(tree.level=='combine.first.trees' & include.germline==F){
      output_trees <- list()

      for (i in 1:length(lineage_dfs)){
        output_trees[[i]] <- list()
        if(!is.null(no.trees.combined)){
          combined_lineage <- do.call('rbind', lineage_dfs[[i]][1:no.trees.combined])
        }else{
          combined_lineage <- do.call('rbind', lineage_dfs[[i]])
        }
        combined_lineage$sequence_id <- as.character(c(1:nrow(combined_lineage)))
        output_trees[[i]][[1]] <- transform_lineage_to_tree_df(combined_lineage)
      }

    }else stop('Unavailable tree configuration')
  }
  return(output_trees)
}
