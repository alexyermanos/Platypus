#' Yields overlap heatmap and datatable of features or combined features for different samples or groups
#' @param VDJ VDJ output of the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]])
#' @param feature.columns A character array of column names of which the overlap should be displayed. The content of these columns is pasted together (separated by "/"). E.g. if the overlap in cells germline gene usage is desired, the input could be c("VDJ_jgene","VDJ_dgene","VDJ_vgene"). These columns would be pasted and compared across the grouping variable.
#' @param grouping.column A column which acts as a grouping variable. If repertoires are to be compared use the sample_id column.
#' @param plot.type Character. Either "ggplot" or "pheatmap". Defaults to Pheatmap
#' @param pvalues.label.size Numeric. Defaults to 4. Is passed on to ggplot theme
#' @param axis.label.size Numeric. Defaults to 4. Is passed on to ggplot theme
#' @param add.barcode.table Boolean. Defaults to T. Whether to generate a dataframe with frequencies and barcodes of cells with overlapping features. This is useful to e.g. analyze deferentially expressed genes between cells of two samples or groups expressing the same VDJ or VJ chain
#' @return A list of a ggplot (out[[1]]), the source table or matrix for the plot out[[2]] and a table containing additional information in case that add.barcode.table was set to TRUE (out[[3]])
#' @export
#' @examples
#' #To test the overlap of CDR3s between multiple samples
#' overlap <- VDJ_overlap_heatmap(VDJ = Platypus::small_vgm[[1]]
#' ,feature.columns = c("VDJ_cdr3s_aa"),
#' grouping.column = "sample_id", axis.label.size = 15
#' , plot.type = "ggplot")
#'

VDJ_overlap_heatmap <- function(VDJ,
                                feature.columns,
                                grouping.column,
                                plot.type,
                                pvalues.label.size,
                                axis.label.size,
                                add.barcode.table){
  group <- NULL
  overlap <- NULL
  overlap_lab <- NULL

  #VERSION is set for now:
  platypus.version <- "v3"

  if(missing(pvalues.label.size)) pvalues.label.size <- 4
  if(missing(axis.label.size)) axis.label.size <- 4
  if(missing(add.barcode.table)) add.barcode.table <- T
  if(missing(plot.type)) plot.type <- "pheatmap"

  #remove any rows that do not contain an entry for a given feature
  to_remove <- c()
  for(n in 1:nrow(VDJ)){
    if("" %in% VDJ[n,c(feature.columns)]){
      to_remove <- c(to_remove, n)}
  }
  if(length(to_remove) > 0){
  VDJ <- VDJ[-to_remove,]
  }
  grouping <- data.frame("group" = VDJ[, grouping.column])
  if(NA %in% grouping$group) stop("NA values in grouping columns. Please choose another column or replace NA values")
  if(length(unique(grouping$group)) < 3 & plot.type == "pheatmap"){
    message("\n Pheatmap plot not possible with less than 3 groups. Returning ggplot")
    plot.type <- "ggplot"
  }

  if(length(feature.columns) > 1){
    grouping$pasted <- do.call(paste, c(VDJ[,c(feature.columns)], sep="/"))
  } else {
    grouping$pasted <- VDJ[, c(feature.columns)]
  }

  sample.names <- unique(grouping[,1])

  df.list <- list()
  for(i in 1:length(unique(grouping[,1]))){
    df.list[[i]] <- unique(subset(grouping, grouping[,1] == unique(grouping[,1])[i])[,2])#get unique values of pasted feature columns per grouping / per repertoire
  }
  names(df.list) <- sample.names

  if(length(sample.names) > 2){
  combs <- as.data.frame(t(utils::combn(as.character(sample.names), m = 2,simplify = TRUE)))#get combinations to test

  combs[,1] <- ordered(as.factor(combs[,1]), levels = (sample.names))
  combs[,2] <- ordered(as.factor(combs[,2]), levels = (sample.names))

  } else {
    combs <- data.frame(sample.names[1], sample.names[2])
  }
  combs$overlap <- NA
  combs$items.overlapping <- NA
  ov_temp_list <- list()
  for(i in 1:nrow(combs)){

    if(all(is.na(df.list[[which(names(df.list) == combs[i,1])]])) == F & all(is.na(df.list[[which(names(df.list) == combs[i,2])]])) == F){
    combs$overlap[i] <- sum(df.list[[which(names(df.list) == combs[i,1])]] %in% df.list[[which(names(df.list) == combs[i,2])]])
    ov_temp <- df.list[[which(names(df.list) == combs[i,1])]][which(df.list[[which(names(df.list) == combs[i,1])]] %in% df.list[[which(names(df.list) == combs[i,2])]])]
    combs$items.overlapping[i] <- paste0(ov_temp, collapse = ";")

    } else { #if any of the two vectors to compare is empty (returning a 0 would be misleading)
      combs$overlap[i] <- NA
      ov_temp <- ""
    }

    if(add.barcode.table == T){
      ov_temp_list[[i]] <- ov_temp
    }
  }

  #now add a third dataframe with frequencies and barcodes of the overlapping elements
  if(add.barcode.table == T){
    if(!"barcode" %in% names(VDJ)) stop("'barcode' column must be present in input dataframe to add barcode table")

    ov_all <- do.call("c", ov_temp_list)
    if(length(ov_all) > 1){
    ov_df <- data.frame("overlapping_items" = ov_all)
    if(length(feature.columns) > 1){
      for(i in 1:length(feature.columns)){
        ov_df[,ncol(ov_df) + 1] <- stringr::str_split(ov_all, "/", simplify = T)[,i]
      }
      names(ov_df)[2:ncol(ov_df)] <- feature.columns
    }

    nc_start <- ncol(ov_df)
    #open some more columns
    for(i in 1:length(sample.names)){
      ov_df[,ncol(ov_df) + 1] <- NA
      ov_df[,ncol(ov_df) + 1] <- NA
      names(ov_df)[c(ncol(ov_df)-1,ncol(ov_df))] <- c(paste0("Freq in ", sample.names[i]),paste0("Barcodes ", sample.names[i]))
    }

    lookup <- cbind(grouping, VDJ[,"barcode"])
    names(lookup)[3] <- "barcode"

    sample_count <- 1
    for(j in seq((1+nc_start),(nc_start + 2*length(sample.names)), 2)){ #this way I get the correct column index directly

      curr_group <- subset(lookup, group == sample.names[sample_count])
      
      for(i in 1:nrow(ov_df)){
        
          ov_df[i,j] <- sum(curr_group$pasted == ov_df$overlapping_items[i]) #frequency per group
          ov_df[i,j+1] <- paste0(curr_group$barcode[which(curr_group$pasted == ov_df$overlapping_items[i])],collapse = ";")
      }
      sample_count <- sample_count + 1
    }
    } else {
      ov_df <- "none"
  }
  } else {
    ov_df <- "none"
  }

  #update: new labeling strategy: to show NA labels for non possible combinations, we reformat the overlap column and print that instead

  combs$overlap_lab <- as.character(combs$overlap)
  combs$overlap_lab[is.na(combs$overlap)] <- "NA"
  if(plot.type == "ggplot"){

  plot_out <- ggplot2::ggplot(combs, ggplot2::aes(x = combs[,1], y = combs[,2],fill=overlap)) + ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label=overlap_lab), size = pvalues.label.size)+ ggplot2::scale_fill_gradient2(low="navy", mid="white", high="red", limits=range(combs$overlap)) + ggplot2::theme(panel.background = ggplot2::element_blank(),axis.text = ggplot2::element_text(size = 30), axis.line.x = ggplot2::element_blank(),axis.line.y = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), text = ggplot2::element_text(size=30), legend.key = ggplot2::element_rect(colour = "white"), legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5, size = 25), plot.subtitle = ggplot2::element_text(size = 15),axis.text.x = ggplot2::element_text(angle = 60,vjust = 1, hjust=1, size = axis.label.size),axis.text.y = ggplot2::element_text(size = axis.label.size)) + ggplot2::labs(title = "", x = "", y = "", subtitle = paste0("Overlap features: " ,paste0(feature.columns, collapse = " ; ")), fill = "") + ggplot2::scale_y_discrete(limits=rev)

  } else {
    #getting symmetric matrix for pheatmap
    allcombs <- c(unlist(as.character(combs[,1])), unlist(as.character(combs[,2]))) #get the unique values from the combinations table
    allcombs <- ordered(as.factor(allcombs), levels = (sample.names)) #order those
    pheat_map <- matrix(data = NA, nrow = length(unique(allcombs)), ncol = length(unique(allcombs))) #make a symmetric matrix template
    colnames(pheat_map) <- unique(allcombs[order(allcombs)]) #rename
    rownames(pheat_map) <- unique(allcombs[order(allcombs)])
    for(i in 1:nrow(combs)){ #assign the corresponding values twice (once for upper triangle and once for lower)
      #upper triangle
      pheat_map[which(colnames(pheat_map) == combs[i,1]), which(rownames(pheat_map) == combs[i,2])] <- as.numeric(combs[i,5])
      #lower triangle
      pheat_map[which(colnames(pheat_map) == combs[i,2]), which(rownames(pheat_map) == combs[i,1])] <- as.numeric(combs[i,5])
    }
    plot_out <- pheatmap::pheatmap(pheat_map,main = paste0("Overlap features: " ,paste0(feature.columns, collapse = " ; ")), border_color = "white", scale = "none", cluster_rows = F, cluster_cols = F,display_numbers = T, number_format = "%.0f", angle_col = 315)

    combs <- pheat_map

  }
  return(list(plot_out,combs,ov_df))
}
