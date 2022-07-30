
#' Function to cutomise the Dot Plot of CellPhoneDB analysis results.
#' @description Function to cutomise the Dot Plot of CellPhoneDB analysis results.
#' @param vgm.input Output of the VDJ_GEX_matrix function. Mandatory. Object where to save the dotplot.
#' @param selected.rows Strings vector. Defaults to NULL. Vector of rows to plot (interacting genes pair), one per line.
#' @param selected.columns Strings vector. Defaults to  NULL. Vector of columns (interacting groups) to plot, one per line
#' @param threshold.type Character vector. Defaults to NULL. Possible arguments: "pvalue", "log2means", "pvalue_topn", "log2means_topn". Which thresholding system the user wants to use.
#' @param threshold.value Numerical. Defaults to NULL. Value below/above (depending on whether it's pvalue or log2means) which genes to plot are selected.
#' @param project.name Character vector. Defaults to NULL. Subfolder where to find the output of the CellPhoneDB analysis and where to save the dot_plot output plot.
#' @param filename Character vector. Defaults to "selection_plot.pdf". Name of the file where the dot_plot output plot will be saved.
#' @param width Numerical. Defaults to 8. Width of the plot.
#' @param height Numerical. Defaults to 10. Height of the plot.
#' @param text.size Numerical. Defaults to 12. Font size of the plot axes
#' @param return.vector Logical. Defaults to FALSE. If set to TRUE, it includes the vector of genes_pairs present in the dot_plot in the VDJ_GEX_matrix[[10]] list.
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter
#' @return VDJ_GEX_matrix object with output dot plot added to VDJ_GEX_matrix[[10]] list.
#' @export
#' @examples
#' \dontrun{
#' vgm_cellphonedb<-dot_plot(vgm.input=vgm_cellphonedb,
#' selected.columns = c("group_1_3|group_2_6", "group_1_3|group_3_9", "group_2_6|group_1_3",
#'"group_2_6|group_3_9", "group_3_9|group_1_3", "group_3_9|group_2_6"),
#'threshold.type="pvalue_topn", threshold.value=50,
#'project.name = "test", height = 12, width=6, text.size=10, return.vector=TRUE)
#'}

dot_plot = function(vgm.input, #vgm where to save the dotplot

                    selected.rows, #Defaults to NULL #vector of rows to plot (interacting genes pair), one per line
                    selected.columns , #Defaults to  NULL #vector of columns (interacting groups) to plot, one per line
                    threshold.type, #Defaults to NULL. Possible arguments: "pvalue", "log2means", "pvalue_topn", "log2means_topn". Which thresholding system the user wants to use.
                    threshold.value, #Defaults to NULL #value below/above (depending on whether it's pvalue or log2means) which select the genes to plot
                    project.name , #Defaults to  NULL
                    filename, #Defaults to "selection_plot.pdf"
                    width, #Defaults to 8
                    height, #Defaults to 10
                    text.size, #Defaults to 12
                    return.vector, #Defaults to FALSE. If set to TRUE, it includes the vector of genes_pairs present in the dot_plot in the vgm CellPhoneDB list item.
                    platypus.version #Defaults to v3"

){

  #For CRAN checks
  count <- NULL
  clusters <- NULL
  pair <- NULL
  pvalue <- NULL

  platypus.version <- "v3"

  #set default value for selected.rows
  if(missing(selected.rows)){
    selected.rows=NULL
  }

  #set default value for selected.columns
  if(missing(selected.columns)){
    selected.columns = NULL
  }

  #set default value for pvalue.threshold
  if(missing(threshold.type)){
    threshold.type = NULL
  }

  #set default value for pvalue.threshold
  if(missing(threshold.value)){
    threshold.value = NULL
  }

  #set default value for filename
  if(missing(filename)){
    filename = "selection_plot.pdf"
  }

  #set default value for project.name
  if(missing(project.name)){
    project.name = NULL
  }

  #set default value for width
  if(missing(width)){
    width = 8
  }

  #set default value for height
  if(missing(height)){
    height = 10
  }

  #set default value for text.size
  if(missing(text.size)){
    text.size = 12
  }

  #set default value for return.vector
  if(missing(return.vector)){
    return.vector = FALSE
  }

  means_separator = "\t"
  pvalues_separator = "\t"

  #get files directories
  my_directory<-getwd()

  if(!is.null(project.name)){
    means_path<-paste(my_directory, "/out/",project.name,"/means.txt", sep ="")
    pvalues_path<-paste(my_directory, "/out/",project.name,"/pvalues.txt", sep ="")
    plot_path<-paste(my_directory, "/out/",project.name, sep ="")
  } else {
    means_path<-paste(my_directory, "/out/means.txt", sep ="")
    pvalues_path<-paste(my_directory, "/out/pvalues.txt", sep ="")
    plot_path<-paste(my_directory, "/out", sep ="")
  }

  all_pval = utils::read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)
  all_means = utils::read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected.rows)){
    selected.rows = intr_pairs
  }

  if(is.null(selected.columns)){
    selected.columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected.rows, intr_pairs), selected.columns]
  sel_means = all_means[match(selected.rows, intr_pairs), selected.columns]
  intr_pairs<-selected.rows
  #select only rows where at least one interacting pair has p-value velow the threshold
  if(!(is.null(threshold.value))&!(is.null(threshold.type))){

    #according to threshold.type values are selected in different ways
    if(threshold.type == "pvalue"){
      sel_col<-NULL
      sel_pval$interacting_pairs<-selected.rows
      #create a column counting how many times there's a pvalue below the threshold in the row
      sel_pval$count<-rep(0, length(selected.rows))
      for(g in 1:length(selected.rows)){
        for(c in 1:length(selected.columns)){
          if(!(is.na(sel_pval[g,c])) & sel_pval[g,c]<=threshold.value){
            sel_pval$count[g]<-sel_pval$count[g]+1
          }
        }
      }
      #only rows with at least one pvalue below the threshold are selected
      sel_pval<-dplyr::filter(sel_pval, count>0)
      sel_pval$count<-NULL
      selected.rows<-NULL
      selected.rows<-sel_pval$interacting_pairs
      sel_pval$interacting_pairs<-NULL
      sel_means = sel_means[match(selected.rows, intr_pairs), ]
    }

    #take the same approach to select the means
    if(threshold.type == "log2means"){
      sel_col<-NULL
      sel_means$interacting_pairs<-selected.rows
      sel_means$count<-rep(0, length(selected.rows))
      for(g in 1:length(selected.rows)){
        for(c in 1:length(selected.columns)){
          if(!(is.na(sel_means[g,c])) & sel_means[g,c]>=threshold.value){
            sel_means$count[g]<-sel_means$count[g]+1
          }
        }
      }
      sel_means<-dplyr::filter(sel_means, count>0)
      sel_means$count<-NULL
      selected.rows<-NULL
      selected.rows<-sel_means$interacting_pairs
      sel_means$interacting_pairs<-NULL
      sel_pval = sel_pval[match(selected.rows, intr_pairs), ]
    }

    if(threshold.type == "pvalue_topn"){
      sel_col<-NULL
      sel_pval$interacting_pairs<-selected.rows

      #pivot the data.frame so that we have all the pvalues in a single column with the respective interacting groups pair
      sel_piv<-tidyr::pivot_longer(sel_pval, dplyr::all_of(selected.columns), names_to = "group_pair", values_to = "temp")
      sel_piv<-cbind(pvalue=sel_piv$temp, sel_piv)
      sel_piv$temp<-NULL

      #order data.frame accoridng to the pvalue values
      sel_piv<-sel_piv[do.call(order,sel_piv),]

      # take only the first n rows
      sel_piv<-sel_piv[1:threshold.value,]
      sel_pval<-sel_pval[match(unique(sel_piv$interacting_pairs), sel_pval$interacting_pairs),]

      selected.rows<-NULL
      selected.rows<-sel_pval$interacting_pairs
      sel_pval$interacting_pairs<-NULL
      sel_means = sel_means[match(selected.rows, intr_pairs), ]
    }

    if(threshold.type == "log2means_topn"){
      sel_col<-NULL
      sel_means$interacting_pairs<-selected.rows

      #pivot the data.frame so that we have all the pvalues in a single column with the respective interacting groups pair
      sel_piv<-tidyr::pivot_longer(sel_means, dplyr::all_of(selected.columns), names_to = "group_pair", values_to = "temp")
      sel_piv<-cbind(mean=sel_piv$temp, sel_piv)
      sel_piv$temp<-NULL

      #order data.frame accoridng to the pvalue values
      sel_piv<-sel_piv[do.call("order",c(sel_piv, list(decreasing=TRUE))),]

      # take only the first n rows
      sel_piv<-sel_piv[1:threshold.value,]
      sel_means<-sel_means[match(unique(sel_piv$interacting_pairs), sel_means$interacting_pairs),]

      selected.rows<-NULL
      selected.rows<-sel_means$interacting_pairs
      sel_means$interacting_pairs<-NULL
      sel_pval = sel_pval[match(selected.rows, intr_pairs), ]
    }
  }

  df_names = expand.grid(selected.rows, selected.columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- grDevices::colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  #make the dot plot, including only the selcted rows and columns
  selection_dot_plot<-ggplot2::ggplot(plot.data,ggplot2::aes(x=clusters,y=pair)) +
    ggplot2::geom_point(ggplot2::aes(size=-log10(pvalue),color=mean)) +
    ggplot2::scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   axis.text=ggplot2::element_text(size=text.size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   axis.text.y = ggplot2::element_text(size=text.size, colour = "black"),
                   axis.title=ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(size = 0.7, linetype = "solid", colour = "black"))

  #save the plot in the output directory
  ggplot2::ggsave(filename, path=plot_path, width = width, height = height, limitsize=FALSE)

  #include the plot in the vgm
  vgm.input[[10]]$selection_dot_plot<-selection_dot_plot

  #if argument return.vector==TRUE, include the vector of interacting gene pairs in the vgm
  if(return.vector==TRUE){
    vgm.input[[10]]$interacting_pairs<-selected.rows
  }
  return(vgm.input)
}
