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
AntibodyForest_network_to_pnp_DL <- function(AntibodyForest.networks, #! Include all parameters you mentioned in the documentation and mention all parameters you include here in the documentation !
                                             metric.function, # String. Needs to be a vertex attribute in iGraph Objects which is either defined as a boolean or as a string of either "yes" or "no" ("yes" will be selected). Later on function input will be implemented to be applied internally to iGraph object to allow for more freedom in metrics
                                             n.per.sample, # Vector of Integers. Select top n networks of each sample to select antibody sequences from.
                                             n.per.network, #
                                             order.networks.per.sample.by.expansion, # Boolean. For Top N networks, networks will be ordered per sample (among networks of one sample) based on cell number (cell count)
                                             igg.default, #Boolean. If true, networks are filtered so that only those remain that have less than 50% IGHM/IGHD
                                             custom.network.filter, # (List of) Expression. Needs the keyword
                                             custom.node.filter,
                                             filter.isotype, # vector of string(s); if wanted to subselect for specific isotype; if set to FALSE no filtering will occur. if set to c('IgK','IgL) kappa and lambda Ig; default is no filtering
                                             filter.VJcgene, # vector of string(s); if wanted to subselect for specific light chain types; if set to FALSE no filtering will occur. if set to c('lambda','kappa) kappa and lambda Ig # default is no filtering
                                             custom.function,
                                             custom.function.args,
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

  sampleVec <- NULL
  sample_id <- NULL
  clonotype_id <- NULL
  cellsInClonotype <- NULL
  .data <- NULL
  . <- NULL
  # Missing inputs, default values and small checks
  if(missing(AntibodyForest.networks)){
    stop("Please provide AntibodyForest.networks input for this function")
  }
  if(missing(metric.function)){
    metric.function <- FALSE
    print("No filtering via metric.function applied")}
  if(!is.character(metric.function)&!is.logical(metric.function)){
    print("Metric.function filters via column names of the resulting data object. Metric.function was set to false.")
    metric.function <- FALSE
  }
    # metric.function <- "most_expanded" #Default
    # print("metric.function parameter was set to most_expanded")} #let the user know in case they did not read the documentation (who does that anyways ;) )

  if(missing(n.per.sample)){
    sampleVec <- rep(10, times<-length(AntibodyForest.networks))
    print(paste0("n.per.sample parameter was set to c(",paste(sampleVec, collapse=", "),")"))
    n.per.sample <- sampleVec
  }

  if(length(n.per.sample)!=length(AntibodyForest.networks)){ #if the user has not supplied the exact amount of values in the vector (only c(1,1), even though there are three sample_id's)
    if(length(n.per.sample)<length(AntibodyForest.networks)){ #if the user has not supplied enough values
      indices<- ((length(n.per.sample)+1):length(AntibodyForest.networks)) #check how many values are necessary
      n.per.sample[indices]<-n.per.sample[1] #fill these indices with the first value
      warning(paste0("Not enough values in n.per.sample \nn.per.sample was expanded to c(",paste0(n.per.sample,collapse=", "),")"))
    }
    else{
      n.per.sample<-n.per.sample[1:length(AntibodyForest.networks)]
      warning(paste0("Too much values in n.per.sample \nn.per.sample was set to c(",paste(n.per.sample,collapse=", "),")"))
    }
  }
  #n.per.sample <- c(10,5,1)
  #print("n.per.sample parameter was set to c(10,5,1)")}
  if(missing(n.per.network)){
    n.per.network <- 1
    print("n.per.network parameter was set to 1")}

  if (missing(order.networks.per.sample.by.expansion)) {
      order.networks.per.sample.by.expansion <- T}

  if(missing(igg.default)){
    igg.default<-TRUE
  }

  if(missing(custom.network.filter)){
    custom.network.filter<-FALSE
    print("custom.network.filter was set to FALSE - No custom network filtering applied")
  }
  #if((all(unlist(lapply(custom.network.filter,is.expression)))==FALSE)&(!custom.network.filter[1]==FALSE)){
  if(!is.logical(custom.node.filter)){
    if((any(unlist(lapply(custom.network.filter,is.expression)))==FALSE)){
      warning("Some elements of custom.network.filter were not of type expression - No custom network filtering applied")
      custom.network.filter<-NULL
      custom.network.filter<-FALSE
    }
  }
  if(missing(custom.node.filter)){
    custom.node.filter<-FALSE
    print("custom.node.filter was set to FALSE - No custom node filtering applied")
  }
  #if((all(unlist(lapply(custom.node.filter,is.expression)))==FALSE)&(!custom.node.filter[1]==FALSE)){
  if(!is.logical(custom.network.filter)){
    if((any(unlist(lapply(custom.network.filter,is.expression)))==FALSE)){
      warning("Some elements of custom.node.filter were not of type expression - No custom node filtering applied")
      custom.node.filter<-NULL
      custom.node.filter<-FALSE
    }
  }

   if(missing(filter.isotype)){
     filter.isotype <- FALSE
     print("filter.isotype was set to FALSE - No isotype filtering applied")  }

  if(missing(filter.VJcgene)){
    filter.VJcgene <- FALSE
    print("filter.VJcgene was set to FALSE - No light chain filtering applied")}


  if(missing(custom.function)){
    custom.function <- FALSE
    print("No custom.function supplied - custom.function was set to FALSE")}


  if(missing(custom.function.args)){
    if(custom.function==FALSE){
    print("No additional arguments supplied to custom.function")}
      }

  if(missing(species)) species <- "mouse"
  if(missing(manual_IgKC)) manual_IgKC <- "none"
  if(missing(manual_2A)) manual_2A <- "none"
  if(missing(manual_VDJLeader)) manual_VDJLeader <- "none"
  if(missing(write.to.disk)) write.to.disk <- T
  if(missing(filename)) filename <- "PnP_assembled_seqs"
  if(missing(removeintron)) removeintron <- F
  if(missing(intron)) intron <- ""
  platypus.version <- "v3" #in case your function is backwards compatible, you will have to include an if clause for "v2" and "v3". Default should be "v2" at the time of writing.


  #Checks if mixcr alignement is there
  if (!"VDJ_nSeqFR1" %in% names(AntibodyForest.networks[[1]][[1]]@node_features)) { #This check is very ARBitratry and may not be the best choice
    stop("mixcr alignment output seems to be missing in iGraph objects of networks. Please provide mixcr output to allow processing to PnP sequences.")
  }


  ############################################################
  ### Variable Initialisation
  ############################################################

  node_dfs_single <- list() # container for selected nodes of 1 clonotype/networks
  node_dfs_total <- list() #container for selected nodes of all clonotype/networks ()
  node_dfs <- list() # container for selected nodes of networks (REDUNDANT) (succeeds node_dfs_total and passes immediately to node_df_subset)

  node_df_subset <- NULL # data.frame that will contain node subsets based on metrics and to select Abs
  node_df_all <- NULL # data.frame that will collect all selected/subsetted nodes to be input into pnp assemble function

  established_isotypes <- c("IGHG1","IGHG2a","IGHG2b","IGHG2c","IGHG3","IGHA","IGHM","IGHD","IGHE","Unknown") #For filter.isotype
  established_lc_designation<-(c("IGK","IGL","IGKC","IGLC1","IGLC2","IGLC3","kappa","lambda")) #For filter.VJcgene

  network.filter.report <- data.frame(matrix(nrow = 0, ncol=5)) #Report so you can check which networks got filtered
  names(network.filter.report) <- c("Filter","Expression","sample_id","clonotype_id","Value")

  node.filter.report <- data.frame(matrix(nrow = 0, ncol=6)) #Report so you can check which nodes got filtered
  names(node.filter.report) <- c("Filter","Expression","sample_id","clonotype_id","label","Value")

  for(i in 1:length(AntibodyForest.networks)){ #For every sample that was put in i <- sample

    ###########################################################
    ### Which networks/clonotypes are we going to get?
    ###########################################################
    network.ensemble <- NULL #Empty list to input the chosen clonotype in

    if(igg.default==TRUE){
      igg.expression <- expression((length(grep(paste(c("IGHM","IGHD"),collapse="|"),network$VDJ_cgene))/length(network$VDJ_cgene))<0.5)
      igg.filter.result <- AntibodyForest_network_filter(AntibodyForest.networks[[i]], list(igg.expression))
      igg.filter.selection <- igg.filter.result[[1]]
    }

    if((is.expression(custom.network.filter)|is.expression(custom.network.filter[[1]]))){ #If something was passed into the custom network filter, do the filtering | is custom.network.filter OR the first element in the list custom.network.filter an expression

      network.filter.result <- AntibodyForest_network_filter(AntibodyForest.networks[[i]],custom.network.filter)
      network.selection <- network.filter.result[[1]]
      network.filter.report <- rbind(network.filter.report,network.filter.result[[2]])

      if(igg.default==TRUE){
        network.selection<-network.selection&igg.filter.selection
      }

      network.ensemble <- names(AntibodyForest.networks[[i]])[network.selection]
    }else{ #No filter on
      if(igg.default==TRUE){# but igg default
        network.ensemble <- names(AntibodyForest.networks[[i]])[igg.filter.selection]
      }else{ #No filter on
        network.ensemble <- names(AntibodyForest.networks[[i]])
        #n.per.sample comes as a last step, as it is possible that some clonotypes fulfill the network criteria but no node criteria

      }
    }


    ###########################################################
    ### Which nodes from the before chosen clonotypes are we going to get?
    ###########################################################
    node_dfs_total <- list() #While it was initialized before, it needs to be reset again for every sample, as the way the loop was built, it's possible that the past entries do not get overwritten

    for(ii in network.ensemble){ #for every clonotype/network chosen from AntibodyForest.networks[[i]]

      node.selection <-    logical(length(igraph::V(AntibodyForest.networks[[i]][[ii]]@tree))) #This is the "end vector" that gets compared to.

      custom.node.filter.Select <- logical(length(igraph::V(AntibodyForest.networks[[i]][[ii]]@tree))) #This one gets overwritten during the filtering process

      #Because of this next line, we do a loop instead of fancy list operations. I did not find a way to do it.
      node_dfs_single[[i]] <- igraph::as_data_frame(AntibodyForest.networks[[i]][[ii]]@tree, what='vertices') #gets continuously overwritten for every sample
      node_dfs_single[[i]]$cellsInClonotype <- sum(node_dfs_single[[i]]$cell_number) #new column

      ### Are we filtering with the custom.node.filter?
      if((is.expression(custom.node.filter)|is.expression(custom.node.filter[[1]]))){ #If something was passed into the custom node filter, do the filtering

        node.filter.result <- AntibodyForest_node_filter(node_dfs_single[[i]],custom.node.filter)
        custom.node.filter.Select <- node.filter.result[[1]]
        node.filter.report <- rbind(node.filter.report,node.filter.result[[2]])
        #maybe sort the report?

      }else{
        custom.node.filter.Select <- !logical(length(igraph::V(AntibodyForest.networks[[i]][[ii]]@tree)))

      }

      ###########################################################
      ### Isotype filtering (filter.isotype)
      ###########################################################

      if(!filter.isotype[1]==FALSE){ #If FALSE, then filter.isotype was missing

        filter.isotype.flag <- TRUE # Raise Flag that there should be a filtering process
        if(!(is.character(filter.isotype))){ #If the filter.isotype is not of type character
          warning("Provided filter.isotype argument is not a string and hence no filtering is applied - \nInitialised designations: IGHG1, IGHG2a, IGHG2b, IGHG2c, IGHG3, IGHM, IGHD, IGHE, Unknown")
          filter.isotype.flag <- FALSE #flag down - no filtering process
          filter.isotype.Select <- !logical(length(node_dfs_single[[i]]$VDJ_cgene))

        }

        if(filter.isotype.flag){ #If filtering process is OK'd, do it
          #Filter function
          filter.isotype.Select <- AntibodyForest_isotype_filter(node_dfs_single,i, filter.isotype, established_isotypes)
        }


      } else{
        filter.isotype.Select <- !logical(length(node_dfs_single[[i]]$VDJ_cgene)) #Set all to TRUE so nothing gets filtered due to filter.isotype
      }

      ###########################################################
      ### VJcgene filtering (filter.VJcgene)
      ###########################################################

      if(!filter.VJcgene[1]==FALSE){ #Filtering? Y/N?
        filter.VJcgene.flag <- TRUE #Raise flag so there is a filtering process
        if(!is.character(filter.VJcgene)){ #Is it even a valid type?
          warning("Provided filter.VJcgene argument is not a string and hence no filtering is applied - \nInitialised designations: IGHG1, IGHG2a, IGHG2b, IGHG2c, IGHG3, IGHM, IGHD, IGHE, Unknown")
          filter.VJcgene.flag <- FALSE #Take flag down
          filter.VJcgene.Select <- !logical(length((node_dfs_single[[i]]$label)))
        }

        #put here function filter.VJcgene
        if(!filter.VJcgene.flag==FALSE){ #If filtering process is okay
          #Filter function
          filter.VJcgene.Select <- AntibodyForest_VJcgene_filter(node_dfs_single,i, filter.VJcgene, established_lc_designation,filter.VJcgene.flag)

         }
      }else{
        filter.VJcgene.Select <- !logical(length((node_dfs_single[[i]]$label))) #sets all to TRUE, so nothing gets thrown out due to filter.VJcgene
        }


      ###########################################################
      ### Selection
      ###########################################################


      node.selection <- custom.node.filter.Select&filter.isotype.Select&filter.VJcgene.Select

      node_dfs_total[[length(node_dfs_total)+1]] <- node_dfs_single[[i]][node.selection,] #get the selected nodes


    } # for(ii in network.ensemble)

    node_dfs[[i]] <- do.call(rbind, node_dfs_total)

    node_df_subset <- NULL
    node_df_subset <- node_dfs[[i]]

    # basically create new column(s) and already filter on those or
    if(!is.logical(custom.function)){
      tryCatch(
      node_df_subset <- do.call(custom.function,append(list(nodes=node_df_subset,ABForest=AntibodyForest.networks),custom.function.args)),
      #node_df_subset <- custom.function(node_df_subset,AntibodyForest.networks[[i]],custom.funtion.args) #HERE I SHOULD PASS ABFOREST TOO
      error = function(e)
        print("Custom.function did not work"))
    }

    # if one used the custom.node.filter, one can filter now for "NodeFilterMatches", which is the amount of TRUE values one got from the applied filters
    if(nrow(node.filter.report)>0){
      node_df_subset$NodeFilterMatches <- 0

      for(t in 1:nrow(node_df_subset)){
        node_df_subset[t,"NodeFilterMatches"]<-length(node.filter.report[(node.filter.report$Expression!="ALL"
                                                                               &node.filter.report$sample_id==node_df_subset[t,"sample_id"]
                                                                               &node.filter.report$clonotype_id==node_df_subset[t,"clonotype_id"]
                                                                               &node.filter.report$label==as.numeric(node_df_subset[t,"label"])),"Value"])
      }
    }

    if (order.networks.per.sample.by.expansion) { #THIS PART (depends on metrics though) AND THE SELECTION OF NETWORKS SHOULD BE MOVED ABOVE METRIC CALCULATION TO REDUCE COMPUTATION (of metrics)
      # Get ordering vector to select then topN networks
      #node_dfs[[i]] <- node_dfs[[i]] %>% dplyr::group_by(sample_id, clonotype_id) %>% dplyr::mutate(cellsInClonotype = sum(cell_number, na.rm = T))
      node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id, clonotype_id) %>% dplyr::arrange(dplyr::desc(cellsInClonotype))
    }



    # use column name from 'metric.function' input
    if(!metric.function==FALSE){
      if(!(metric.function%in%names(node_df_subset))){
        stop("The supplied column name via argument column.function does not exist in the resulting data object.")
      }
      if (is.character(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::filter(.data[[metric.function]] == "yes")
      } else if (is.logical(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::filter(.data[[metric.function]] == TRUE)
      } else if (is.integer(node_df_subset[[metric.function]]) | is.double(node_df_subset[[metric.function]])) {
        node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id, clonotype_id) %>%
          dplyr::slice_max(., order_by = .data[[metric.function]], n = n.per.network, with_ties = F) %>%
          dplyr::filter(!is.na(.data[[metric.function]])) %>% dplyr::ungroup() %>% dplyr::arrange(nchar(.data$clonotype_id),.data$clonotype_id) # to ensure no NA values are kept need to remove them (NAs could arise from germlines!)
      }
    }

    # Only 10 clonotypes/networks, or specified amount of the clonotypes/networks
    if(n.per.sample[i]<length(unique(node_df_subset$clonotype_id))){
    node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id) %>% dplyr::filter(clonotype_id %in% unique(.$clonotype_id)[1:n.per.sample[i]])
    }else{
      node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id) %>% dplyr::filter(clonotype_id %in% unique(.$clonotype_id))
    }
    # Only 1 node per clonotype/network
    ## if the certain column name is there - order by it
    node_df_subset <- node_df_subset %>% dplyr::group_by(sample_id, clonotype_id) %>% dplyr::slice_head(.,n = n.per.network)

    # collect subset dfs
    node_df_all <- rbind(node_df_all, node_df_subset)

  } #for(i in 1:length(AntibodyForest.networks))

  #print(node_df_all)

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
    VDJ_assemble_for_PnP_github_DL(VDJ.mixcr.matrix = new_node_df_all,
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

#-----------------------------------------------------------------------------
# Additional functions
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# AntibodyForest_network_filter
#-----------------------------------------------------------------------------

#' Filter function for the AntibodyForest_network_to_pnp function
#' Input: AntibodyForest for 1 sample and a (list of) expression(s), which are conditions one wants to filter the clonotypes/networks which are in the ABForest.
#' @param AntibodyForest.networks.sample One sample (sample_ids need to be unique) from AntibodyForest.networks, which is a list of clonotypes/igraph objects (AbForest output containing mixcr alignemnt results in network/igraph objects)
#' @param custom.network.filter List of expressions, which are to evaluated using eval(). Multiple expressions are "OR'd". If one wishes to "AND" multiple expressions, one has to write that explicitely in one condition
#' @param network.selection Boolean vector. Is iterated upon over the conditions and is returned at the end.
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter.
#' @return A boolean vector with TRUE at the indexes of the clonotypes that passed all the conditions set in custom.network.filter and FALSE on all the other ones


AntibodyForest_network_filter <- function(AntibodyForest.networks.sample,custom.network.filter){ #THIS INPUT IS ALREADY AntibodyForest.networks[[i]]
  #network.preselection <- logical(length=length(AntibodyForest.networks.sample))
  for(z in 1:length(custom.network.filter)){
    tempVar <- deparse(custom.network.filter[[z]], width.cutoff = 500)
    tempVar <- gsub("network","igraph::V(network)",tempVar)
    tempVar <- sub("expression","", tempVar)
    #custom.network.filter[[z]]<-NULL
    custom.network.filter[[z]]<- parse(text=tempVar, keep.source = FALSE)
  }

  network.preselection <- FALSE
  network.selection <-    logical(length=length(AntibodyForest.networks.sample))

  network_filter_report <- data.frame(matrix(nrow = 0, ncol=5))
  names(network_filter_report) <- c("Filter","Expression","sample_id","clonotype_id","Value")

  iii <- 1 #loop count

  for(i in 1:length(custom.network.filter)){ # go through all of the conditions
    for(ii in 1:length(AntibodyForest.networks.sample)){ #for every network/clonotype
      network <- AntibodyForest.networks.sample[[ii]]@tree #making the network ready for the expression
      nodes.selection <- eval(custom.network.filter[[i]]) #should give a boolean out for every node if it fits the expression or not
      if(!is.logical(nodes.selection)){
        stop("Your expression did not result in a logical value ('TRUE'\\\'FALSE')")
      }
      if(length(nodes.selection)>1 & ii == 1){
        warning("Your expression resulted in a longer logical vector than expected.\n\"Any()\" was applied on the logical vector.")
      }
      network.preselection[ii] <- any(nodes.selection)#ticks for the current condition i if there's a true node or not
      network_filter_report[iii,] <- c("Network_Filter",deparse(custom.network.filter[[i]], width.cutoff = 500),igraph::vertex_attr(network,"sample_id")[1],igraph::vertex_attr(network,"clonotype_id")[1],network.preselection[ii])
      iii <- iii+1
    }#For loop end - through all the networks
    network.selection <- network.selection|network.preselection
  }#For loop end - through all the conditions

  network_filter_endreport <- data.frame(Filter = "Network_Filter",
                                         Expression = "ALL",
                                         sample_id = unique(network_filter_report$sample_id),
                                         clonotype_id = unique(network_filter_report$clonotype_id),
                                         Value = network.selection)

  network_filter_report <- rbind(network_filter_report,network_filter_endreport)

  To_Return <- list(network.selection,network_filter_report)
  return(To_Return)
} #function end

#-----------------------------------------------------------------------------
# AntibodyForest_node_filter
#-----------------------------------------------------------------------------

#' Filter function for the AntibodyForest_network_to_pnp function
#' Input: Dataframe with node features for 1 specific clonotype of 1 specific sample and a (list of) expression(s), which are conditions one wants to filter the clonotypes/networks which are in the ABForest.
#' @param nodes Dataframe of node_features extracted from the igraph objects from the AntibodyForest object
#' @param custom.node.filter List of expressions, which are to evaluated using eval()
#' @param node.selection Boolean vector. Is iterated upon over the conditions specified in custom.node.filter and is returned at the end.
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter.
#' @return A boolean vector equaling the number of nodes - which nodes pass the filter


AntibodyForest_node_filter <- function(nodes,custom.node.filter){
  node.preselection <- logical(length=nrow(nodes))
  node.selection <-    logical(length=nrow(nodes))

  node_filter_report <- data.frame(matrix(nrow = 0, ncol=6))
  names(node_filter_report) <- c("Filter","Expression","sample_id","clonotype_id","label","Value")


  for(i in 1:length(custom.node.filter)){ # go through all of the conditions
    node_filter_prereport <- node_filter_report
    node.preselection <- eval(custom.node.filter[[i]])
    if(!is.logical(node.preselection)){
      stop("Your expression did not result in a logical value ('TRUE'\\\'FALSE')")
    }
    node.selection <- node.selection|node.preselection

    node_filter_prereport<- data.frame(Filter = "Node_Filter",
                                       Expression = deparse(custom.node.filter[[i]], width.cutoff = 500),
                                       sample_id = nodes[,c("sample_id")],
                                       clonotype_id = nodes[,c("clonotype_id")],
                                       label = nodes[,c("label")],
                                       Value = node.preselection)
    node_filter_report <- rbind(node_filter_report,node_filter_prereport)

  }#For loop end - through all the conditions

  node_filter_endreport <- data.frame(Filter = "Node_Filter",
                                      Expression = "ALL",
                                      sample_id = nodes[,c("sample_id")],
                                      clonotype_id = nodes[,c("clonotype_id")],
                                      label = nodes[,c("label")],
                                      Value = node.selection)
  node_filter_report <- rbind(node_filter_report,node_filter_endreport)

  To_Return <- list(node.selection,node_filter_report)

  return(To_Return)
} #function end

#-----------------------------------------------------------------------------
# AntibodyForest_isotype_filter
#-----------------------------------------------------------------------------

#' Filter function for the AntibodyForest_network_to_pnp function
#' Input: Dataframe with node features for 1 specific clonotype of 1 specific sample and a vector of character strings, which represent the $VDJ_cgene features one wants to filter
#' @param node_df_single Dataframe of node_features extracted from the igraph objects from the AntibodyForest object - already eventually processed by custom.node.filter
#' @param filter.isotype Vector of character strings. Should be values found in the $VDJ_cgene column
#' @param established.isotypes Character string of established VDJ_cgene designations. If the strings in filter.isotype deviate from the established_isotypes, it will give a warning.
#' @return A boolean vector equaling the number of nodes - which nodes pass the filter

AntibodyForest_isotype_filter <- function(node_dfs_single,i, filter.isotype, established_isotypes){
  #Check if the provided filter is a part of the established_isotypes or something else
  #There will be a warning, but it nevertheless tries to filter

  #Warnings

  if(!any(grepl(paste(filter.isotype,collapse="|"),established_isotypes, ignore.case = TRUE))){
    #Do any strings have subsets of the established designations in them but are not the designations themselves?
    #This distinction is important because it then does not make it into the warning due to the way the code was written
    subsetStringSelect <- logical(length(established_isotypes)) #Establish an empty logical vector
    for(k in 1:length(established_isotypes)){ #find the designations which are subsets for the specific filter.VJcgene
      subsetStringSelect[k]<-subsetStringSelect[k]|any(grepl(established_isotypes[k],filter.isotype))
    }
    #Put them into this vector
    subsetStrings <- paste0(filter.isotype[grepl(paste(established_isotypes[subsetStringSelect],collapse="|"),filter.isotype)&established_isotypes[subsetStringSelect]!=filter.isotype],collapse = " ,")
    if(subsetStrings>1){
      warning(paste0("filter.isotype '",paste(grep(paste(established_isotypes, collapse="|"),filter.isotype, ignore.case = TRUE, invert=TRUE, value=TRUE), collapse = ", "),
                     "' and '",
                     subsetStrings,
                     "' does not match any prespecified isotype designation \nData structure might be empty after filtering \nInitialised designations: IGHG1, IGHG2a, IGHG2b, IGHG2c, IGHG3, IGHA, IGHM, IGHD, IGHE, Unknown"))
    } else{
      warning(paste0("filter.isotype '",paste(grep(paste(established_isotypes, collapse="|"),filter.isotype, ignore.case = TRUE, invert=TRUE, value=TRUE), collapse = ", "),
                     "' does not match any prespecified isotype designation \nData structure might be empty after filtering \nInitialised designations: IGHG1, IGHG2a, IGHG2b, IGHG2c, IGHG3, IGHA, IGHM, IGHD, IGHE, Unknown"))
    }
  }

  #Actual filtering
  filter.isotype.PreSelect <- NULL #running variable vector
  filter.isotype.Select <- logical(length(node_dfs_single[[i]]$VDJ_cgene)) #Vector you end up

  #for loop to check over the provided criteria
  for(isoI in filter.isotype){
    filter.isotype.PreSelect <- grepl(isoI,node_dfs_single[[i]]$VDJ_cgene)
    filter.isotype.Select<- (filter.isotype.Select|filter.isotype.PreSelect)
  }

return(filter.isotype.Select)
}#function end

#-----------------------------------------------------------------------------
# AntibodyForest_VJcgene_filter
#-----------------------------------------------------------------------------

#' Filter function for the AntibodyForest_network_to_pnp function
#' Input: Dataframe with node features for 1 specific clonotype of 1 specific sample and a vector of character strings, which represent the $VDJ_cgene features one wants to filter
#' @param node_df_single Dataframe of node_features extracted from the igraph objects from the AntibodyForest object - already eventually processed by custom.node.filter
#' @param filter.VJcgene Vector of character strings. Should be values found in the $VJ_cgene or $VJchain column
#' @param established_lc_designation Character string of established light chain designations. If the strings in filter.VJcgene deviate from the established_lc_designations, it will give a warning.
#' @return A boolean vector equaling the number of nodes - which nodes pass the filter


AntibodyForest_VJcgene_filter <- function(node_dfs_single,i, filter.VJcgene, established_lc_designation,filter.VJcgene.flag){
  #Check if the provided filter is a part of the established_lc_designations or something else
  #There will be a warning, but it nevertheless tries to filter

  #Warnings
  if(!any(grepl(paste(filter.VJcgene, collapse="|"),established_lc_designation, ignore.case = TRUE))){

    #Do any strings have subsets of the established designations in them but are not the designations themselves?
    subsetStringSelect <- logical(length(established_lc_designation)) #Establish an empty logical vector
    for(k in 1:length(established_lc_designation)){ #find the designations which are subsets for the specific filter.VJcgene
      subsetStringSelect[k]<-subsetStringSelect[k]|any(grepl(established_lc_designation[k],filter.VJcgene))
    }
    #Put them into this vector
    subsetStrings <- paste0(filter.VJcgene[grepl(paste(established_lc_designation[subsetStringSelect],collapse="|"),filter.VJcgene)&established_lc_designation[subsetStringSelect]!=filter.VJcgene],collapse = " ,")
    if(subsetStrings>1){
      warning(paste0("filter.VJcgene '",
                     paste(grep(paste(established_lc_designation, collapse="|"),filter.VJcgene, ignore.case = TRUE, invert=TRUE, value=TRUE),collapse = ", "),
                     "' and '",
                     subsetStrings,
                     "' does not match any prespecified light chain designation \nData structure might be empty after filtering \nInitialised designations: IGK, IGL, IGKC, IGLC1, IGLC2,IGLC3, kappa, lambda"))
    } else{
      warning(paste0("filter.VJcgene '",
                     paste(grep(paste(established_lc_designation, collapse="|"),filter.VJcgene, ignore.case = TRUE, invert=TRUE, value=TRUE),collapse = ", "),
                     "' does not match any prespecified light chain designation \nData structure might be empty after filtering \nInitialised designations: IGK, IGL, IGKC, IGLC1, IGLC2,IGLC3, kappa, lambda"))
    }
  }
  #Substituting kappa with IGK
  if(any(grepl("kappa",filter.VJcgene, ignore.case = TRUE))){
    print(" 'kappa' was substitued with 'IGK' ")
    filter.VJcgene <- unique(gsub("kappa","IGK",filter.VJcgene))
  }
  #Substituting lambda with IGL
  if(any(grepl("lambda",filter.VJcgene, ignore.case = TRUE))){
    print(" 'lambda' was substitued with 'IGL' ")
    filter.VJcgene <- unique(gsub("lambda","IGL",filter.VJcgene))
  }

  #Initialising flags
  #Two flags due to two possible columns that one could filter with
  filter.VJcgene.vjchainflag <- TRUE
  filter.VJcgene.vjcgeneflag <- TRUE

  #Checking which columns are available
  if(!all(c("VJ_chain","VJ_cgene")%in%names(node_dfs_single[[i]]))){
    #one or both columns are missing
    if(!any(c("VJ_chain","VJ_cgene")%in%names(node_dfs_single[[i]]))){ #both columns are missing
      warning("Both columns (VJ_chain, VJ_cgene) necessary for filtering are missing \nPlease adjust your AntibodyForest call and add the column names to node.feature \nProcessing continues without light chain filtering")
      filter.VJcgene.flag <- FALSE
      #Here set a flag so it exits this filtering process
    }
    else if(!"VJ_chain"%in%names(node_dfs_single[[i]])){
      warning("Column VJ_chain is missing \nFor optimal filtering results please include columns VJ_chain and VJ_cgene in the AntibodyForest call")
      #Here only filtering with VJ_chain
      filter.VJcgene.vjchainflag <- FALSE

    }
    else if(!"VJ_cgene"%in%names(node_dfs_single[[i]])){
      warning("Column VJ_cgene is missing \nFor optimal filtering results please include columns VJ_chain and VJ_cgene in the AntibodyForest call")
      #Here only filtering with VJ_cgene
      filter.VJcgene.vjcgeneflag <- FALSE
    }
  }



  #Filtering process
  if(filter.VJcgene.flag){ #If at least one column is around

    #Initialising selection vectors
    filter.VJcgene.PreSelect <- NULL
    filter.VJcgene.Select <- logical(length((node_dfs_single[[i]]$label)))

    if(filter.VJcgene.vjchainflag==TRUE){ #vjchainflag is raised
      for(isoF1 in filter.VJcgene[which(nchar(filter.VJcgene)<4)]){ #VJchain is either IGL/IGK, so nchar <4
        # if(!any(grepl(isoF1,established_lc_designation))){
        #   #print(paste0("filter.VJcgene ",isoF1," does not match any prespecified light chain designation - Data structure might be empty \nInitialised designations: IGK, IGL, IGKC, IGLC1, IGLC2, IGLC3, kappa, lambda"))
        #   }
        filter.VJcgene.PreSelect <- grepl(isoF1,node_dfs_single[[i]]$VJ_chain)
        filter.VJcgene.Select<- (filter.VJcgene.Select|filter.VJcgene.PreSelect)
      }
    }
    if(filter.VJcgene.vjcgeneflag==TRUE){ #VJcgeneflag is raised
      for(isoF2 in filter.VJcgene[which(nchar(filter.VJcgene)>=4)]){ #VJcgene is normally >4
        # if(!any(grepl(isoF2,established_lc_designation))){
        #   #print(paste0("filter.VJcgene ",isoF2," does not match any prespecified light chain designation - Data structure might be empty \nInitialised designations: IGK, IGL, IGKC, IGLC1, IGLC2, IGLC3, kappa, lambda"))
        #   }
        filter.VJcgene.PreSelect <- grepl(isoF2,node_dfs_single[[i]]$VJ_cgene)
        filter.VJcgene.Select<- (filter.VJcgene.Select|filter.VJcgene.PreSelect)
      }
    }
  }else{
    filter.VJcgene.Select <- !logical(length((node_dfs_single[[i]]$label)))
  }

  return(filter.VJcgene.Select)

}
#-----------------------------------------------------------------------------
# Next function?
#-----------------------------------------------------------------------------


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
