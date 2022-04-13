#' Phylogenetic tree plotting
#'
#' @description Function to plot phylogenetic trees obtained from VDJ_phylogenetic_trees

#'@description !Requires the ggtree package to be loaded! Plots trees from function VDJ_phylogenetic_trees
#' @param tree.dfs nested list of tidytree dataframes obtained from VDJ_phylogenetic_trees with output.format='tree.df.list'. tree.dfs[[1]][[2]] represent a tree dataframe for the first sample, second clonotype.
#' @param color.by string - VDJ or tree df column name which will be used to color the tree nodes.
#' @param size.by string or NULL - VDJ or tree df column name which determines the node size. If NULL, node sizes will be equal.
#' @param shape.by string or NULL - VDJ or tree df column name which determines the node shape. If NULL, node sizes will be equal.
#' @param specific.leaf.colors named list or NULL - if NULL, colors will be automatically selected for each node according to its color.by value.
#' @return nested list of ggtree plot objects for each sample and each clonotype.
#' @export
#' @examples
#' \dontrun{
#' VDJ_phylogenetic_trees_plot(tree.dfs,color.by='clonotype_id', size.by='sequence_frequency')
#'}

VDJ_phylogenetic_trees_plot <- function(tree.dfs,
                           color.by,
                           size.by,
                           shape.by,
                           specific.leaf.colors){

 if(missing(tree.dfs)) stop('Please input a nested list of tidytree dataframes from VDJ_clonal_lineages(ouput.format=tree.df.list)')
 if(missing(color.by)) color.by <- 'clonotype_id'
 if(missing(size.by)) size.by <- 'sequence_frequency'
 if(missing(shape.by)) shape.by <- NULL
 if(missing(specific.leaf.colors)) specific.leaf.colors <- NULL

 #! No require call inside functions. Please update package DESCRIPTION instead !
 #require(ggtree)
 #require(tidytree)

 output_plots <- list()

 for(i in 1:length(tree.dfs)){
   output_plots[[i]] <- list()
   for(j in 1:length(tree.dfs[[i]])){
     label_df <- tree.dfs[[i]][[j]][which(!is.na(tree.dfs[[i]][[j]]$label) & tree.dfs[[i]][[j]]$germline=='no'),]
     label_df[color.by] <- as.character(unlist(label_df[color.by]))
    # label_df[shape.by] <- as.character(label_df[shape.by])
     label_df[size.by] <- as.integer(unlist(label_df[size.by]))

      #Trycatch to make sure that ggtree is loaded
     tryCatch({
        #R CRAN checks throw an note concerning the %<+% operator.
        #One way is to formally import this via Roxygen2 as
        #@importFrom ggree %<+%. This is problematic as ggtree has therefore to be installed
        #upon loading of Platypus. Given that this package is still young and that is showed some clashes with rvcheck,
        #we want to avoid ggtree being necessary to load Platypus.
        #To therefore mask the function from CRAN checks, we provide it as a string and eval that as a command.
        #Not elegant but it works...
        eval(parse(text = "output_plots[[i]][[j]] <- ggtree::ggtree(ape::as.phylo(tree.dfs[[i]][[j]])) %<+% label_df + ggtree::geom_tippoint(ggplot2::aes_string(color=color.by, size=size.by, shape=shape.by)) + ggplot2::ggtitle(paste0(unique(as.character(label_df$sample_id)), '_', unique(as.character(label_df$clonotype_id)))) + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.text = ggplot2::element_blank()) + ggplot2::labs(color=paste0(color.by), size=paste0(size.by), shape=paste0(shape.by))"))

     }, error=function(e){
        print(e)
        stop("ggtree-based plot could not be generated. Please ensure that ggtree is loaded (library(ggtree)) and the %<+% operator is available in the current R session")
     })

     if(!is.null(specific.leaf.colors)){
       output_plots[[i]][[j]] <- output_plots[[i]][[j]] + ggplot2::scale_color_manual(values = specific.leaf.colors)
     }
   }
 }

 return(output_plots)
}
