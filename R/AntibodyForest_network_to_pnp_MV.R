#' From Networks to PnP Sequences
#' Extracts from Abforest network function. Initial Version; extracts n most expanded Antibody Sequences per Network from Network/Tree
#' @param AntibodyForest.networks List of samples (sample_ids need to be unique), which is a list of clonotypes/igraph objects (AbForest output containing mixcr alignemnt results in network/igraph objects)
#' @param metric.function String. Name of vertex attribute to be selected by. It will be selected for nodes having either "yes", TRUE or the highest numeric. (function that can describe metric better (i.e. should be applyable to network and return T or F if it should be included; later on))
#' @param n.per.sample Integer or Vector of Integers with length of nbr of samples. Defaults to 1.
#' @param n.per.network Integer. Defaults to 1. Enter this value only if metric.function (e.g column name in iGraph function: cell_number == expansion) is a number; will select top n highest values. How many of top expanded clones from each network should be at most picked. E.g., if set to 2, for the top 2 expanded clones PnP sequence constructs will be extracted.
#' @param filter.isotype not implemented.
#' @param filter.VJcgene not implemented.
#' @param species Character. Which IgKC sequence to use. Can be "human" or "mouse". Defaults to "mouse".
#' @param manual_IgKC Character. Manual overwrite for sequence used as IgKC.
#' @param manual_2A Character. Manual overwrite for sequence used as Furine 2A site.
#' @param manual_VDJLeader Character. Manual overwrite for sequence used as VDJ Leader and signal peptide.
#' @param write.to.disk Boolean. Defaults to TRUE. Whether to save assembled sequences to working directory.
#' @param filename Character. Output file name for .fasta and .csv files if write.to.disk == T. Defaults to PnP_assembled_seqs.fasta/.csv.
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter.
#' @return A data.frame containing PnP sequences of selected sequences of networks.
#' @export
#' @examples
#' \dontrun{
#' To return results for a non default column
#' physiologically_relevant <- to_something(VDJ = VDJ_GEX_matrix.output[[1]]
#' , column.to.plot = "VJ_jgene", normalization factor = 20)
#'}

#-----------------------------------------------------------------------------
# To do:
# add input on which Ig isotype to select
# add constant sequences for IgKappa and Iglambda
# and adapt the sequence input for pnp assemble depending on that
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Example of how the code was run for the OAS Flu Project, where the VDJ signal and the 2A sequence were exchanged:
# AntibodyForest_network_to_pnp_MV(output_ABforest,
#                                  metric.function = "most_expanded",
#                                  n.per.sample = 30, # Vector of Integers. Select top n networks of each sample to select antibody sequences from.
#                                  n.per.network = 1, # How many sequences per network
#                                  order.networks.per.sample.by.expansion = T, # Boolean. For Top N networks, networks will be ordered per sample (among networks of one sample) based on cell number (cell count)
#                                  species = "mouse",
#                                  manual_IgKC = "none",
#                                  manual_2A = "CGCAAACGACGGGGAAGCGGAGCTACTAACTTCAGCCTGCTGAAGCAGGCTGGAGACGTGGAGGAGAACCCTGGACCT",
#                                  manual_VDJLeader= "ATGATGGTGTTAAGTCTTCTGTACCTGTTGACAGCCCTTCCGGGTATCCTGTCA",
#                                  write.to.disk = T,
#                                  filename = desiredname,
#                                  whichgroup = groupname
# )


#-----------------------------------------------------------------------------



# Name to AntibodyForest_..
AntibodyForest_network_to_pnp_MV <- function(AntibodyForest.networks, #! Include all parameters you mentioned in the documentation and mention all parameters you include here in the documentation !
                                    metric.function, # String. Needs to be a vertex attribute in iGraph Objects which is either defined as a boolean or as a string of either "yes" or "no" ("yes" will be selected). Later on function input will be implemented to be applied internally to iGraph object to allow for more freedom in metrics
                                    n.per.sample, # Vector of Integers. Select top n networks of each sample to select antibody sequences from.
                                    n.per.network, #
                                    order.networks.per.sample.by.expansion, # Boolean. For Top N networks, networks will be ordered per sample (among networks of one sample) based on cell number (cell count)
                                    filter.isotype, # vector of string(s); if wanted to subselect for specific isotype; if set to FALSE no filtering will occur. if set to c('IgK','IgL) kappa and lambda Ig; default is no filtering
                                    filter.VJcgene, # vector of string(s); if wanted to subselect for specific light chain types; if set to FALSE no filtering will occur. if set to c('lambda','kappa) kappa and lambda Ig # default is no filtering
                                    species,
                                    manual_IgKC,
                                    manual_2A,
                                    manual_VDJLeader,
                                    write.to.disk,
                                    filename,
                                    whichgroup,
                                    platypus.version,
                                    removeintron,
                                    intron){

  #-----------------------------------------------------------------------------
  # Basic set-ups
  #-----------------------------------------------------------------------------
  cell_number <- NULL
  sample_id <- NULL
  clonotype_id <- NULL
  cellsInClonotype <- NULL
  .data <- NULL
  . <- NULL
  # Missing inputs, default values and small checks
  if(missing(AntibodyForest.networks)) stop("Please provide AntibodyForest.networks input for this function")
  if(missing(metric.function)){
    metric.function <- "most_expanded" #Default
    print("metric.function parameter was set to most_expanded")} #let the user know in case they did not read the documentation (who does that anyways ;) )
  if(missing(n.per.sample)){
    n.per.sample <- c(10,5,1)
    print("n.per.sample parameter was set to c(10,5,1)")}
  if(missing(n.per.network)){
    n.per.network <- 1
    print("n.per.network parameter was set to 1")}
  if (missing(order.networks.per.sample.by.expansion)) {
    order.networks.per.sample.by.expansion <- T}
  if(missing(species)) species <- "mouse"
  if(missing(manual_IgKC)) manual_IgKC <- "none"
  if(missing(manual_2A)) manual_2A <- "none"
  if(missing(manual_VDJLeader)) manual_VDJLeader <- "none"
  if(missing(write.to.disk)) write.to.disk <- T
  if(missing(filename)) filename <- "PnP_assembled_seqs"
  if(missing(removeintron)) removeintron <- F
  if(missing(intron)) intron <- ""
  platypus.version <- "v3" #in case your function is backwards compatible, you will have to include an if clause for "v2" and "v3". Default should be "v2" at the time of writing.


  metric_function_params <- c("most_expanded", "hub") # simple metrics that are callable by names
  # from me: added @tree
  metric_function_params <- unique(c(metric_function_params, names(igraph::vertex.attributes(AntibodyForest.networks[[1]][[1]]@tree)))) # add to those any column name that is contained in network

  # if n.per.sample is a single integer expand to nbr of samples

  # if not a custom can skip custom function calculation --> need to check actually in tree!! ideally to ensure one can add also a column by itself
  custom_function <- T
  if(metric.function %in% metric_function_params){
    custom_function <- F
  }
  # here also check that AntibodyForest.networks contains mixcr alignment output. e.g. nSeqFR1.... # just check for one
  if (!"VDJ_nSeqFR1" %in% metric_function_params) {
    stop("mixcr alignment output seems to be missing in iGraph objects of networks. Please provide mixcr output to allow processing to PnP sequences.")
  }



  #  Variable Initialisation

  subset_col <- NULL # subset_column to be added to networks based on which variants are selected from network
  node_dfs <- list() # container for selected nodes of networks
  node_dfs_single <- list() # container for selected nodes of networks
  node_dfs_total <- list()
    #table_with_cellnbrOrdered <- NULL # nodes table per sample which will be ordered by nbr of cells in corresponding clonotpye
  node_df_subset <- NULL # data.frame that will contain node subsets based on metrics and to select Abs
  node_df_all <- NULL # data.frame that will collect all selected/subsetted nodes to be input into pnp assemble function

  #-----------------------------------------------------------------------------
  # Calculate Metrics needed for each network (ideally by using Abforest_metrics)
  # add unique meaningful naming scheme (i.e. s001_c0001_n0001 for sample, clonotype and node label; node label is ordered by expansion usually)
  #-----------------------------------------------------------------------------
  # loop over each sample and each network and calculate based on metric subset_col
  if (custom_function) {
    for (sample_name in names(AntibodyForest.networks)) {
      for (sample_name in names(AntibodyForest.networks)) {
        # pass; add custom func
        # create 'subset_col' in each iGraph object as vertex attribute which is based on metrics function used
      }
    }
  }



  #-----------------------------------------------------------------------------
  # Subset networks for wanted sequences
  # Create plot (if wanted) that shows which sequences will be selected
  #-----------------------------------------------------------------------------
  # Perform following actions for each sample
  for(i in 1:length(AntibodyForest.networks)){
    # First get for each network a data frame

    #this part is from me: loop through each clonotype to get a data frame from each network. old code was with lapply
    for (num_clone in 1:n.per.sample[1]){
      node_dfs_single[[i]] <- igraph::as_data_frame(AntibodyForest.networks[[i]][[num_clone]]@tree, what='vertices')
      node_dfs_total[[num_clone]] <- node_dfs_single[[i]]
      #print(node_dfs_single[[i]])
      # Collapse the list of data.frames per sample to one big data.frame per sample from which we will select

    }
    node_dfs[[i]] <- do.call(rbind, node_dfs_total)

    #old code
    # from me : added [[1]]@tree
    #node_dfs[[i]] <- lapply(AntibodyForest.networks[[i]][[1]]@tree, function(x) igraph::as_data_frame(x, what='vertices'))

    if (order.networks.per.sample.by.expansion) { #THIS PART (depends on metrics though) AND THE SELECTION OF NETWORKS SHOULD BE MOVED ABOVE METRIC CALCULATION TO REDUCE COMPUTATION (of metrics)
      # Get ordering vector to select then topN networks
      node_dfs[[i]] <- node_dfs[[i]] %>% dplyr::group_by(sample_id, clonotype_id) %>% dplyr::mutate(cellsInClonotype = sum(cell_number, na.rm = T))
      node_dfs[[i]] <- node_dfs[[i]] %>% dplyr::arrange(dplyr::desc(cellsInClonotype))
    }
    # Select nodes/antibody/cell variants based on metric and n.per.sample
    node_df_subset <- NULL
    node_df_subset <- node_dfs[[i]] %>% dplyr::filter(clonotype_id %in% unique(node_dfs[[i]]$clonotype_id)[1:n.per.sample[i]])
    if (custom_function) {
      # use selection_col
      stop("not impemented yet")
    } else {
      # use column name from 'metric.function' input
      if (is.character(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::filter(.data[[metric.function]] == "yes")
      } else if (is.logical(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::filter(.data[[metric.function]] == TRUE)
      } else if (is.integer(node_df_subset[[metric.function]]) | is.double(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id, clonotype_id) %>%
          dplyr::slice_max(., order_by = .data[[metric.function]], n = n.per.network, with_ties = F) %>%
          dplyr::filter(!is.na(.data[[metric.function]])) %>% dplyr::ungroup() # to ensure no NA values are kept need to remove them (NAs could arise from germlines!)
      }
    }

    # collect subset dfs
    node_df_all <- rbind(node_df_all, node_df_subset)
  }


  print(node_df_all)



  #-----------------------------------------------------------------------------
  # If Mixcr output is not provided yet; run alignment via call_mixcr
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # If several isotypes: need to filter and seperate dataframes to only contain one because of PnP_assemble function
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # Create PnP Sequences using PnP_assemble function
  #-----------------------------------------------------------------------------
  # after filtering data frame to wanted cells, we can assembly PnP sequences
  # create a unique identifier consisting of sampleID, clonotypeID and network label
  #old code
  #new_node_df_all <- node_df_all %>% dplyr::mutate(pnp_id = paste0(sample_id,"_", clonotype_id, "_l",label), .before =1) %>%
    #as.data.frame() %>% dplyr::mutate_if(is.list,as.character)
  #new code
  new_node_df_all <- node_df_all %>% dplyr::mutate(pnp_id = paste0(whichgroup,"_", clonotype_id), .before =1) %>%
    as.data.frame() %>% dplyr::mutate_if(is.list,as.character)
  print(new_node_df_all)

  print(new_node_df_all$VJ_nSeqCDR3)


  PnP_seq_output <-# network contains mixcr outputs as lists of 1 char but PnP assemble functions needs them as characters
    VDJ_assemble_for_PnP_github_MV(VDJ.mixcr.matrix = new_node_df_all,
                         id.column = "pnp_id",
                         species = "mouse",
                         manual_IgKC = manual_IgKC,
                         manual_2A = manual_2A,
                         manual_VDJLeader = manual_VDJLeader,
                         write.to.disk = write.to.disk,
                         filename = filename,
                         verbose = T,
                         removeintron = removeintron,
                         intron = intron)
  #-----------------------------------------------------------------------------
  #
  #-----------------------------------------------------------------------------





  #  11. Output
  #      In case you return a list make sure to include the output structure in your documentation under @return
  return(PnP_seq_output)
}

#Once done with writing the function there are a couple more things to do:

#  12. Dependencies and imports
#     If you are using functions from other package e.g. ggplot2::ggplot these packages need to be references in the
#     DESCRIPTION file of the R package. If you use a new package that no other function of Platypus uses, please let me
#     or Alex know, so we can update the DESCRIPTION file.
#     => Again make sure that you always write package::function() in case its not from the base package.

#  13. Documenting
#      Please make sure that your code is understandable to others. Add comments to important functions or lines

#  14. Debugging
#      Please debug your code. Test all possible scenarios and the limit of your functions. If needed, provide warnings
#      and recommendations to the user so that the function stays stable.

#  15. Operating systems
#      Some packages and external software (e.g. doParallel, MIXCR, velocity) behave differently depending on the
#      operating system. In case you are using one of these, include the following snippet to make sure you know what the
#      computer is running...

#switch(Sys.info()[['sysname']],
#       Windows= {print("Windows system detected")
#         operating.system <- "Windows"},
#       Linux  = {print("Linux system detected")
#         operating.system <- "Linux"},
#       Darwin = {print("MAC system detected")
#         operating.system <- "Darwin"})

#...and provide the user with a warning:
#if(operating.system == "Darwin") print("Question: How can you tell which one of your friends has the new gold iMac?
#Answer: Don???t worry, he will let you know")
