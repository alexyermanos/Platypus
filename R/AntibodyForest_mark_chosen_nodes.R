#' From PnP Sequences to ABForest plots (targets the plot.ready slot)
#' Takes identifiers from extracted PnP sequences and marks then the respective nodes in the ABForest plot
#' @param AntibodyForest.networks List of samples (sample_ids need to be unique), which is a list of clonotypes/igraph objects (AbForest output containing mixcr alignemnt results in network/igraph objects)
#' @param Seq_output Dataframe with extracted sequences. Output of VDJ_assemble_for_PnP_github_MV or correspondingly the AntibodyForest_network_to_pnp function
#' @param Frame.color Character. Specifying with which color one wants the nodes to be marked.
#' @param Plot.Suppress Boolean. Whether the plots should be plotted on function call.
#' @param specific.node.colors Vector of strings. Colors for the legend for the specific isotypes. Normally already defined.
#' @return A data.frame containing PnP sequences of selected sequences of networks.
#' @export
#' @examples
#' \dontrun{
#' what
#' To return results for a non default column
#' physiologically_relevant <- to_something(VDJ = VDJ_GEX_matrix.output[[1]]
#' , column.to.plot = "VJ_jgene", normalization factor = 20)
#'}





AntibodyForest_mark_chosen_nodes <- function(AntibodyForest.networks, Seq_output,Frame.color, Plot.Suppress, specific.node.colors){
  #is "plot_ready" there?
  if(missing(AntibodyForest.networks)){
    stop("Argument \"AntibodyForest.networks\" is missing. Please supply the output of the \"AntibodyForests\" call")
  }
  if(missing(Seq_output)){
    stop("Argument \"Seq_output\" is missing. Please supply the output of the \"AntibodyForest_network_to_pnp\" call")
  }

  if(missing(Frame.color)){
    print("Argument \"Frame.color\" is missing. Nodes present in the \"AntibodyForest_network_to_pnp\" will be marked with a\"red\" outline")
    Frame.color<-"red"
  }
  if(missing(Plot.Suppress)){
    print("Argument \"Plot.Suppress\" is missing. Plot.Suppress is set to FALSE and the resulting objects will be plotted. ")
    Plot.Suppress <- FALSE
  }

  if(methods::.hasSlot(AntibodyForest.networks[[1]][[1]],"plot_ready")==FALSE){
    stop("The supplied AntibodyForest object does not have a slot \"plot_ready\"\n
         Please use the function \"AntibodyForests_plot\" on the object that resulted from the \"AntibodyForests\" call. ")
  }

  PlotList <- NULL
  sample_id <- NULL

  for(i in 1:length(unique(Seq_output$sample_id))){
    for(ii in 1:length(unique(subset(Seq_output, sample_id==unique(Seq_output$sample_id)[i])$clonotype_id))){
      clonotype_ensemble <- unique(subset(Seq_output, sample_id==unique(Seq_output$sample_id)[i])$clonotype_id)
      k <- unique(Seq_output$sample_id)[i] #this is to make sure, that the sample matches between the outmost loop and treecodes

      treecodes <-igraph::vertex_attr(AntibodyForest.networks[[k]][[clonotype_ensemble[ii]]]@plot_ready,"cell_barcodes") #all the codes from the current tree
      cellbarc <- Seq_output$cell_barcodes[which((Seq_output$sample_id==k)&(Seq_output$clonotype_id==clonotype_ensemble[ii]))] #all the codes from the sequences

      treee <- (AntibodyForest.networks[[k]][[clonotype_ensemble[ii]]]@plot_ready) #the unmarked tree

      # igraph::vertex_attr(treee,"frame.width")
      # igraph::V(treee)$frame.width<-1
      # igraph::V(treee)$frame.width[which(treecodes%in%cellbarc==TRUE)]<- 3 #Framewidth of chosen node
      # igraph::V(treee)$frame.color[which(treecodes%in%cellbarc==TRUE)]<- "black" #Framecolor of chosen node

      igraph::vertex_attr(treee,"frame.width") #initializing the frame.width attribute in that treee
      igraph::vertex_attr(treee,"frame.color") #initializing the frame.color attribute in that treee

      igraph::V(treee)$frame.width<-1 #setting all the frame.widths to 1
      igraph::V(treee)$frame.color<-"black"
      igraph::V(treee)$frame.width[which(treecodes%in%cellbarc)]<- 3 #Framewidths of chosen nodes
      igraph::V(treee)$frame.color[which(treecodes%in%cellbarc)]<- Frame.color #Framecolor of chosen nodes

      if(!Plot.Suppress){
        plot(treee)
        graphics::title(paste0(k," - ",clonotype_ensemble[ii]))
        graphics::legend(
          "left",
          legend = unlist(names(specific.node.colors)),
          pt.bg  = unlist(specific.node.colors),
          pch    = 21,
          cex    = 1,
          bty    = "n",
          title  = "VDJ_cgene"
        )
      }

      AntibodyForest.networks[[k]][[clonotype_ensemble[ii]]]@plot_ready<-treee



    }
  }
  #Markdown complains

  #return(AntibodyForest.networks)
}
