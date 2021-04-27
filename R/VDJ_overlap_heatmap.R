#' Yields overlap heatmap and datatable of features or combined features for different samples or groups
#' @param VDJ.matrix VDJ matrix from the VDJ_GEX_matrix function 
#' @param feature.columns A character array of column names of which the overlap should be displayed. The content of these columns is pasted together (separated by "/"). E.g. if the overlap in cells germline gene usage is desired, the input could be c("VDJ_jgene","VDJ_dgene","VDJ_vgene"). These columns would be pasted and compared across the grouping variable.
#' @param grouping.column A column which acts as a grouping variable. If repertoires are to be compared use the sample_id column. 
#' @param pvalues.label.size Numeric. Defaults to 4. Is passed on to ggplot theme
#' @param axis.label.size Numeric. Defaults to 4. Is passed on to ggplot theme
#' @param add.barcode.table Boolean. Defaults to T. Whether to generate a dataframe with frequencies and barcodes of cells with overlapping features. This is useful to e.g. analyze deferentially expressed genes between cells of two samples or groups expressing the same VDJ or VJ chain
#' #' @param platypus.version Character. At the moment this function runs only on the output of the VDJ_GEX_matrix function meaning that it is exclusively part of Platypus "v3". With further updates the functionality will be extended.  
#' @return A list of a ggplot (out[[1]]) and the table from which the ggplot was generated (out[[2]])
#' @export
#' @examples
#' \dontrun{
#' #To test the overlap of barcodes between multiple samples
#' barcode_overlap <- VDJ_overlap_heatmap(VDJ.matrix = VDJ_GEX_matrix_output[[1]], feature.columns = c("orig_barcode") ,grouping.column = "repertoire", axis.label.size = 15)
#' }

VDJ_overlap_heatmap <- function(VDJ.matrix.output,
                                feature.columns,
                                grouping.column,
                                pvalues.label.size,
                                axis.label.size,
                                add.barcode.table,
                                platypus.version){
  
  #VERSION is set for now:
  platypus.version <- "v3"
  
  if(missing(pvalues.label.size)) pvalues.label.size <- 4
  if(missing(axis.label.size)) axis.label.size <- 4
  if(missing(add.barcode.table)) add.barcode.table <- T
  
  grouping <- data.frame("group" = VDJ.matrix.output[, grouping.column])
  if(length(feature.columns) > 1){
    grouping$pasted <- do.call(paste, c(VDJ.matrix.output[, c(feature.columns)], sep="/"))
  } else {
    grouping$pasted <- VDJ.matrix.output[, c(feature.columns)]
  }
  
  sample.names <- unique(grouping[,1])
  df.list <- list()
  for(i in 1:length(unique(grouping[,1]))){
    df.list[[i]] <- unique(subset(grouping, grouping[,1] == unique(grouping[,1])[i])[,2])#get unique values of pasted feature columns per grouping / per repertoire 
  }
  names(df.list) <- sample.names
  
  if(length(sample.names) > 2){
  combs <- as.data.frame(t(combn(sample.names, m = 2,simplify = TRUE)))#get combinations to test

  combs[,1] <- ordered(as.factor(combs[,1]), levels = rev(sample.names))
  combs[,2] <- ordered(as.factor(combs[,2]), levels = sample.names)
  
  } else {
    combs <- data.frame(sample.names[1], sample.names[2])
  }
  combs$overlap <- NA
  combs$items.overlapping <- NA
  ov_temp_list <- list()
  for(i in 1:nrow(combs)){
    print(i)
    combs$overlap[i] <- sum(df.list[[which(names(df.list) == combs[i,1])]] %in% df.list[[which(names(df.list) == combs[i,2])]])
    ov_temp <- df.list[[which(names(df.list) == combs[i,1])]][which(df.list[[which(names(df.list) == combs[i,1])]] %in% df.list[[which(names(df.list) == combs[i,2])]])]
    combs$items.overlapping[i] <- paste0(ov_temp, collapse = ";")
    
    if(add.barcode.table == T){
      ov_temp_list[[i]] <- ov_temp
    }
  }
  
  #now add a third dataframe with frequencies and barcodes of the overlapping elements
  if(add.barcode.table == T){
    if(!"barcode" %in% names(VDJ.matrix.output)) stop("'barcode' column must be present in input dataframe to add barcode table")
  
    ov_all <- do.call(c, ov_temp_list)
    ov_df <- data.frame("overlapping_items" = ov_all)
    if(length(feature.columns) > 1){
      for(i in 1:length(feature.columns)){
        ov_df[,ncol(ov_df) + 1] <- str_split(ov_all, "/", simplify = T)[,i]
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
    
    lookup <- cbind(grouping, VDJ.matrix.output[,"barcode"])
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
  
  plot_out <- ggplot(combs, aes(x = combs[,2], y = combs[,1],fill=overlap)) + geom_tile() +geom_text(aes(label=overlap), size = pvalues.label.size)+ scale_fill_gradient2(low="navy", mid="white", high="red", limits=range(combs$overlap)) + theme(panel.background = element_blank(),axis.text = element_text(size = 30), axis.line.x = element_blank(),axis.line.y = element_blank(), axis.ticks = element_blank(), text = element_text(size=30), legend.key = element_rect(colour = "white"), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 25), plot.subtitle = element_text(size = 15),axis.text.x = element_text(angle = 60,vjust = 1, hjust=1, size = axis.label.size),axis.text.y = element_text(size = axis.label.size)) + labs(title = "", x = "", y = "", subtitle = paste0("Overlap features: " ,paste0(feature.columns, collapse = " ; ")), fill = "")
  
  print(plot_out)
  return(list(plot_out,combs,ov_df))  
}


