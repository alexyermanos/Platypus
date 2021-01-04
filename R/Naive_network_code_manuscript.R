VDJ_mutational_network <- function(lineage.list,distance.calculation,network.algorithm,weighted.edges,weighted.germline,random.seed,resolve.ties){
  library(data.table)
  library(scales)
  library(reshape2)
  library(igraph)
  library(stringdist)
  if(missing(distance.calculation)) distance.calculation <- "lv"
  if(missing(network.algorithm)) network.algorithm <- "naive"
  if(missing(weighted.edges)) weighted.edges <- F
  if(missing(weighted.germline)) weighted.germline <- F
  if(missing(random.seed)) random.seed <- 1
  if(missing(resolve.ties)) resolve.ties <- "random"
  
  
  # global parameters
  set.seed(random.seed)
  num.networks <- length(lineage.list)
  cells.per.network <- rep(0,num.networks)
  variants.per.network <- rep(0,num.networks)
  variant.sequences <- list()
  network.objects <- list()
  distance.matrices <- list()
  adj.matrices <- list()
  return.list <- list()
  germline.index.list <- list()
  
  # nested list parameters
  cells.per.variant <- list()
  cell.indicies.per.variant <- list()
  pasted.variant.names <- list()
  
  new.variant.names <- list()
  
  ## total isotype states and total transcriptome states 
  
  temp_rbind <- (do.call("rbind",lineage.list))
  #unique.clusters <- 
  #unique.isotypes <- 
  isotype.per.cell <- list()
  transcriptome.cluster.per.cell <- list()
  isotype.per.variant <- list()
  transcriptome.cluster.per.variant <- list()
  
  ## for each lineage loop
  for(i in 1:length(lineage.list)){
    ## record global parameters
    cells.per.network[i] <- length(lineage.list[[i]]$Seq)
    variant.sequences[[i]] <- unique(lineage.list[[i]]$Seq)
    variants.per.network[i] <- length(variant.sequences[[i]])
    return.list[[i]] <- list()
    
    ## per cell parameters
    ## need to extract cluster - remove all things before the first two _
    clusters <- gsub(".*_","",lineage.list[[i]]$Name)
    clusters_integer <- as.integer(gsub("[^0-9.]", "",  clusters))
    clusters_integer <- clusters_integer - 1
    clusters_integer[is.na(clusters_integer)==T] <- "Unknown"
    clusters_integer[which(grepl(lineage.list[[i]]$Name,pattern = "germline")==T)] <- "Germline"
    transcriptome.cluster.per.cell[[i]] <- clusters_integer
    print(transcriptome.cluster.per.cell[[i]])
    isotypes <- sub(".*?_", "", lineage.list[[i]]$Name)
    isotypes1 <- sub(".*?_", "", isotypes)
    isotypes2 <- gsub("_.*","",isotypes1)
    isotypes2[which(grepl(lineage.list[[i]]$Name,pattern = "germline")==T)] <- "Germline"
    isotype.per.cell[[i]] <- isotypes2
    #print(isotype.per.cell[[i]])
    isotype.per.variant[[i]] <- rep("",variants.per.network[i])
    transcriptome.cluster.per.variant[[i]] <- rep("",variants.per.network[i])
    
    ## generate nested list for per variant object
    cells.per.variant[[i]] <- rep(0,variants.per.network[i])
    cell.indicies.per.variant[[i]] <- list()
    pasted.variant.names[[i]] <- rep("",variants.per.network[i])
    
    new.variant.names[[i]] <- rep("",variants.per.network[i])
    
    ## extracting per variant information
    for(j in 1:length(variant.sequences[[i]])){
      cell.indicies.per.variant[[i]][[j]] <- which(lineage.list[[i]]$Seq==variant.sequences[[i]][j])
      cells.per.variant[[i]][j] <- length(cell.indicies.per.variant[[i]][[j]])
      pasted.variant.names[[i]][j] <- paste(lineage.list[[i]]$Name[cell.indicies.per.variant[[i]][[j]]],sep="",collapse = ";")
      isotype.per.variant[[i]][j] <- paste(isotype.per.cell[[i]][cell.indicies.per.variant[[i]][[j]]],sep="",collapse = ";")
      transcriptome.cluster.per.variant[[i]][j] <- paste(transcriptome.cluster.per.cell[[i]][cell.indicies.per.variant[[i]][[j]]],sep="",collapse = ";")
      
      
    }
    ## now including lineage number, 
    new.variant.names[[i]] <- paste("L",i,"_","v",1:length(new.variant.names[[i]]),"_","c",cells.per.variant[[i]],sep="")
    
    ## now need to change variant to germline 
    new.variant.names[[i]][which(grepl(pasted.variant.names[[i]],pattern = "germline")==T)] <- paste("L",i,"_Germline",sep="")
    
    
    # Simple case, no pruning
    if(network.algorithm=="naive"){
      distance.matrices[[i]] <- stringdistmatrix(variant.sequences[[i]],variant.sequences[[i]],method=distance.calculation)
      rownames(distance.matrices[[i]]) <- new.variant.names[[i]]
      colnames(distance.matrices[[i]]) <- new.variant.names[[i]]
      
      adj.matrices[[i]] <- matrix(0,nrow=nrow(distance.matrices[[i]]),ncol=ncol(distance.matrices[[i]]))
      rownames(adj.matrices[[i]]) <- new.variant.names[[i]]
      colnames(adj.matrices[[i]]) <- new.variant.names[[i]]
      germline.index <- which(grepl(new.variant.names[[i]],pattern = "Germline")==T)
      germline.index.list[[i]] <- germline.index
      #distance.matrices[[i]][upper.tri(distance.matrices[[i]])] <- Inf
      diag(distance.matrices[[i]]) <- Inf
      all.nodes <- 1:nrow(distance.matrices[[i]])
      current.nodes.in.network <- germline.index
      nodes.not.in.network <- all.nodes[-current.nodes.in.network]
      j <- 1
      while(length(nodes.not.in.network)>0){
        print(j)
        if(j==1){         
          x <- distance.matrices[[i]][germline.index,]
          closest.nodes <- which(x==min(x))
          
          if(weighted.germline==F) gl.dist <- 1
          else if(weighted.germline==T) gl.dist <- min(distance.matrices[[i]][germline.index,])
          for(k in 1:length(closest.nodes)){
            
            adj.matrices[[i]][germline.index,closest.nodes[k]] <- gl.dist
            adj.matrices[[i]][closest.nodes[k],germline.index] <- gl.dist
            
          }
          current.nodes.in.network <- unique(c(current.nodes.in.network,closest.nodes))
          nodes.not.in.network <- all.nodes[-current.nodes.in.network]
        }## germline set
        if(j>1){ ############################################################################
          print(paste('current nodes in network are',paste(current.nodes.in.network,collapse = ';')))
          print(paste('current nodes not in network are',paste(nodes.not.in.network,collapse = ';')))          
          x <- distance.matrices[[i]]
          x[nodes.not.in.network,] <- Inf
          x[,current.nodes.in.network] <- Inf
          
          next.min.distance <- min(x)
          
          next.min.distance.indices <- which(x == min(x), arr.ind = TRUE)
          print((next.min.distance.indices))
          print(class(next.min.distance.indices))
          new.unique.newnodes <- unique(next.min.distance.indices[,2])
          unique.column.nodes <- (which(duplicated(next.min.distance.indices[,2])==F))
          
          if(resolve.ties=='first'){
            duplicated.index <- which(duplicated(next.min.distance.indices[,2])==T)
            closest.nodes.to.current <- matrix(next.min.distance.indices[unique.column.nodes,],nrow=length(unique.column.nodes))
          }
          else if(resolve.ties=='random'){
            new.unique.newnodes <- unique(next.min.distance.indices[,2])
            if(length(next.min.distance.indices[,2])==length(unique(next.min.distance.indices[,2]))){
              closest.nodes.to.current <- (next.min.distance.indices)
            }
            else{
              keep.row <- rep(NA,length(new.unique.newnodes))
              for(m in 1:length(new.unique.newnodes)){
                temp.which <- which(next.min.distance.indices[,2]==new.unique.newnodes[m]) # 1 and 3
                print(temp.which)
                if(length(temp.which)>1){
                  sample.which <- sample(temp.which,size = 1)
                  keep.row[m] <- sample.which
                }
                else{
                  keep.row[m] <- temp.which
                }
                
                
              }
              closest.nodes.to.current <- matrix(next.min.distance.indices[keep.row,],nrow=length(new.unique.newnodes))
            }
            print(closest.nodes.to.current)
            print(paste('nrow of closest.nodes.to.current:',nrow(closest.nodes.to.current)))
            
            ## need to get row and column indices of each one, then update the adj matrix
            if(weighted.edges==F) next.min.distance <- 1
            else if(weighted.edges==T) next.min.distance <- min(x)
            
            adj.row.index <- rep(NA,nrow(closest.nodes.to.current))
            adj.col.index <- rep(NA,nrow(closest.nodes.to.current))
            
            for(k in 1:nrow(closest.nodes.to.current)){
              print(paste('which(rownames(distance.matrices[[i]])==rownames(x)[closest.nodes.to.current[k,1]])',rownames(x)[closest.nodes.to.current[k,1]]))
              
              adj.row.index[k] <- closest.nodes.to.current[k,1]
              adj.col.index[k] <- closest.nodes.to.current[k,2]
              adj.matrices[[i]][adj.row.index[k],adj.col.index[k]] <- next.min.distance
              adj.matrices[[i]][adj.col.index[k],adj.row.index[k]] <- next.min.distance
              
            }
          }
          ## need to move the rows to in network
          current.nodes.in.network <- unique(c(current.nodes.in.network,adj.col.index))
          print(paste('current nodes in network are',paste(current.nodes.in.network,collapse = ';')))
          nodes.not.in.network <- all.nodes[-current.nodes.in.network]
          print(paste('current nodes not in network are',paste(nodes.not.in.network,collapse = ';')))          
          #if(j==15) nodes.not.in.network <- NULL
        }## end of j not 1
        j <- j + 1
        
      }
    }
    ## Your algorithm
    else if(network.algorithm=="pruning"){
      distance.matrices[[i]] <- stringdistmatrix(variant.sequences[[i]],variant.sequences[[i]],method=distance.calculation)
      rownames(distance.matrices[[i]]) <- new.variant.names[[i]]
      colnames(distance.matrices[[i]]) <- new.variant.names[[i]]
      germline.index <- which(grepl(new.variant.names[[i]],pattern = "Germline")==T)
      
      ## Your algorithm
      
      
      
    }
    diag(adj.matrices[[i]]) <- NA
    network.objects[[i]] <- igraph::graph_from_adjacency_matrix(adj.matrices[[i]], mode = c("undirected"), weighted = NULL, diag = FALSE,add.colnames = NULL, add.rownames = NA)
  }
  return.list[[1]] <- network.objects
  return.list[[2]] <- distance.matrices
  return.list[[3]] <- adj.matrices
  return.list[[4]] <- cells.per.network
  return.list[[5]] <- variants.per.network
  return.list[[6]] <- variant.sequences
  return.list[[7]] <- cells.per.variant
  return.list[[8]] <- cell.indicies.per.variant
  return.list[[9]] <- pasted.variant.names
  return.list[[10]] <- new.variant.names
  return.list[[11]] <- germline.index.list
  return.list[[12]] <- isotype.per.cell
  return.list[[13]] <- transcriptome.cluster.per.cell
  return.list[[14]] <- isotype.per.variant
  return.list[[15]] <- transcriptome.cluster.per.variant
  
  return(return.list)
}

construct_structures<-function(list1,list2,list3,list4,opt){
  df <- rbindlist(Map(data.frame, list1, list2,
                      MoreArgs=list(stringsAsFactors=FALSE)),fill=TRUE)
  list_dfs<-tapply(as.list(df), gl(ncol(df)/2, 2), as.data.frame)
  list_dfs<-lapply(list_dfs, na.omit)
  list_dfs<-lapply(list_dfs, setNames, nm = c('sequences','count_uniq'))
  list_dfs<-lapply(list_dfs, FUN = function(x) {rownames(x)  <- NULL;x })
  list3 <- lapply(list3,function(x) { strsplit(x,";") })
  names_cols<-names(table(unlist(lapply(list3, function(x) lapply(x,unique)))))
  count_el<-lapply(list3,function(x) lapply(x, function(z) {
    z<-t(as.data.frame(table(z)))
    rownames(z)<-NULL
    colnames(z) <- z[1,]
    z<-z[-1,]
  }))
  df_iso<-dplyr::bind_rows(count_el)
  df_iso <- as.data.frame(sapply(df_iso, as.numeric)) 
  df_iso[is.na(df_iso)] <- 0
  lengths_lists<-as.vector(unique(unname(unlist(sapply(list_dfs,lengths)))))
  list_iso<-split(df_iso,rep(1:length(lengths_lists),lengths_lists )) 
  list_iso<-lapply(list_iso, FUN = function(x) {rownames(x)  <- NULL;x })
  merged.list.of.dfs<-mapply(data.frame, list_dfs, list_iso, SIMPLIFY = FALSE)
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {names(x) <- gsub("X", "", names(x), fixed = TRUE);x})
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {
    x[,paste(colnames(x[,3:ncol(x)]), "Ratio", sep="")]<-x[,3:ncol(x)]/x[,2]
    x
  })
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x){
    x$Names<- 1:nrow(x)
    x$Names <- sub("^", "B", x$Names)
    col_germ<-unname(which(sapply(names(x), function(x) any(grepl("Germline", x)))))
    row_germ<-which(apply(x[,col_germ][1], 1, function(r) r==1))
    x$Names[row_germ]<-"Germline"
    x$SeqLen<-nchar(as.vector(x$sequences))
    x
  })
  merged.list.of.dfs<-mapply(data.frame, merged.list.of.dfs, list4, SIMPLIFY = FALSE)
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {names(x) <- c(names(x[-length(x)]), "Names_L_v_c");x})
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {rownames(x) <- x$Names_L_v_c;x})
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {names(x) <- gsub("X", "", names(x), fixed = TRUE);x})
  merged.list.of.dfs<-lapply(merged.list.of.dfs, function(x) {ind<-sapply(x, function(z) sum(z==0)) != nrow(x);x[,ind]})
  return(merged.list.of.dfs)
}

add_edge_weights<-function(edges,weight,data,mat){
  
  data<-as.matrix(data)
  weights<-c()
  for (j in 1:nrow(edges)){
    row<-edges[j,]$From
    col<-edges[j,]$To
    row_f<-names(which(data[,"Names"]=="Germline", arr.ind=TRUE))
    if (row==row_f){
      if (weight==TRUE){
        mat[row,col]<-mat[row,col]-abs(as.integer(data[row,'SeqLen'])-as.integer(data[col,'SeqLen']))
      }
      else{
        mat[row,col]<-1
      }
    }
    weights<-c(weights,mat[unlist(row),unlist(col)])
    
  }
  return(weights)
}

vertex_scale<-function(data, uniq_seq,scaleByClonalFreq){
  
  if(scaleByClonalFreq==TRUE){
    node_size<-mapply(as.numeric(uniq_seq), rep(2.5,length(uniq_seq)), FUN = function(x, y) {10+y*x/8})
  }else{
    print('ko')
    node_size<-6+0.1*nrow(data)
  }
  
  label_size <- 0.08*node_size
  return(list(node_size = node_size, label_size = label_size))
}

create_graph<-function(graph,dfs,distMat,germ_index,clonal_frequency,scaleByClonalFreq,weight,opt){
  graphs<-c()
  
  graphs<-lapply(dfs, function(k) {
    #print(paste(parent.frame()$i))
    progr_flag<-FALSE
    columns <- grep('Ratio', colnames(k), value=TRUE)
    columns_x<-gsub('.{5}$', '', columns)
    pos_germ<-which(grepl("Germline",columns_x))
    colrs_sp<-c()
    if(opt=="isotype"){
      #print('enter')
      colors<-c("gray85","red","springgreen2","brown","white","yellow","purple","darkorange","pink","darkgreen")
      names(colors)<-c("Germline","IGHA","IGHG1","IGHG2B","None","IGHG2C","IGHM","IGHG3","IGHD","IGHE")
      
      for (co in names(colors)){
        if (any(grepl(co,columns))){
          colrs_sp<-c(colrs_sp,colors[co])
        }
      }
      colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
    }else{
      if(any(grepl("Unknown",columns_x)) | nchar(trimws(columns_x[pos_germ]))<=8){
        progr_flag<-TRUE
        first_cols<-c("gray85", "#A9A9A9")
        colors<-c(first_cols,scales::hue_pal()(30))
        names_temp<-c("0",as.character(c(1:29)))
        names(colors)<-c(c("Germline","Unknown"),names_temp)
        
      }else{
        colors<-columns_x
        pos_germ<-which(grepl("Germline",columns_x))
        colors[pos_germ]<-gsub("^.*?Germline"," ",colors[pos_germ])
        colors<-gsub(" ", "", colors, fixed = TRUE)
        names(colors)<-c("0",1:length(colors))
        names(colors)[pos_germ]<-"Germline"
        
      }
      
      colrs_sp<-colors[order(match(names(colors),columns_x))]
    }
    difer<-length(colrs_sp)-length(columns_x)
    if(length(colrs_sp)>length(columns_x)){
      colrs_sp<-head(colrs_sp, -difer)
    }
    values <- lapply(1:nrow(k), function(x) unlist(k[x,columns]))
    asd<-lapply(values, function(x) names(which((unlist(x) != 0))))
    max_length <- max(unlist(lapply (asd, FUN = length)))
    asd <- lapply (asd, function (x) {length (x) <- max_length; return (x)})
    matter <- do.call("rbind",asd)
    mapping<-c()
    for (co in names(colors)){
      mapping<-c(mapping,columns[grepl(co,columns_x)])
    }
    mapping<-unique(mapping)
    mapping<-mapping[order(match(mapping,columns))]
    col_cor<-colrs_sp
    names(col_cor)<-mapping
    matter[] <- col_cor[unlist(matter)]
    col_mat<-matter[]
    col_mat<-split(col_mat,1:nrow(col_mat))
    col_mat<-lapply(col_mat, function(x) x[!is.na(x)])
    shapes<-c()
    shapes<-lapply(col_mat, function(x) ifelse(length(x)>1,x<-"pie","circle"))
    colrs <-lapply(values, function(x) colrs_sp)
    col_mat<-lapply(col_mat, function(x) ifelse(length(x)>1,x<-"#01FF01",x))
    edges_df<-igraph::get.edgelist(graph[[as.integer(parent.frame()$i)]])
    edges_df<-as.data.frame(edges_df)
    names(edges_df)<-c('From','To')
    igraph::E(graph[[as.integer(parent.frame()$i)]])$weight <-add_edge_weights(edges_df,weight,k,distMat[[parent.frame()$i]])
    igraph::E(graph[[as.integer(parent.frame()$i)]])$color <- "black"
    if(clonal_frequency==TRUE & scaleByClonalFreq==TRUE){
      uniq_seq <- k$count_uniq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$name<-uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$scale<- uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$label<-uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$order<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$clone_freq<-uniq_seq
    }else if(scaleByClonalFreq==TRUE){
      uniq_seq <- k$count_uniq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$name<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$scale<- uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$label<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$order<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$clone_freq<-uniq_seq
    }else if(clonal_frequency==TRUE){
      uniq_seq <- k$count_uniq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$name<-uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$scale<- 1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$clone_freq<-uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$label<-uniq_seq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$order<-1:nrow(k)
    }else{
      uniq_seq <- k$count_uniq
      igraph::V(graph[[as.integer(parent.frame()$i)]])$name<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$scale<- NA
      igraph::V(graph[[as.integer(parent.frame()$i)]])$label<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$order<-1:nrow(k)
      igraph::V(graph[[as.integer(parent.frame()$i)]])$clone_freq<-uniq_seq
    }
    igraph::V(graph[[as.integer(parent.frame()$i)]])$size<-vertex_scale(k,igraph::V(graph[[as.integer(parent.frame()$i)]])$scale,scaleByClonalFreq)$node_size
    igraph::V(graph[[as.integer(parent.frame()$i)]])$shape<-unlist(shapes)
    igraph::V(graph[[as.integer(parent.frame()$i)]])$color<-unlist(col_mat)
    igraph::V(graph[[as.integer(parent.frame()$i)]])$pie<-values
    igraph::V(graph[[as.integer(parent.frame()$i)]])$pie.color<-colrs
    igraph::V(graph[[as.integer(parent.frame()$i)]])$arrow.size<-0.03
    igraph::V(graph[[as.integer(parent.frame()$i)]])$edge.label.font <-0.5
    igraph::V(graph[[as.integer(parent.frame()$i)]])$edge.label.cex <- 1
    igraph::V(graph[[as.integer(parent.frame()$i)]])$edge.arrow.size <- 1
    igraph::V(graph[[as.integer(parent.frame()$i)]])$label.cex <- vertex_scale(k,igraph::V(graph[[as.integer(parent.frame()$i)]])$scale,scaleByClonalFreq)$label_size
    igraph::V(graph[[as.integer(parent.frame()$i)]])$label.color<-"black"
    coords<-igraph::layout_as_tree(graph[[as.integer(parent.frame()$i)]],root = germ_index[[parent.frame()$i]], circular = FALSE, rootlevel = numeric(), mode = "all", flip.y = TRUE)
    graph[[as.integer(parent.frame()$i)]]$layout<-coords
    leg_pos<-"topright"
    leg_bty<-"n"
    if(!(opt=="isotype") & progr_flag){
      stings_only<-c("Unknown","Germline")
      n_cols <- names(colrs_sp)[! names(colrs_sp) %in% stings_only]
      ord<-as.character(sort(as.numeric(n_cols), decreasing=F))
      if(any(grepl("Unknown", columns_x)) & any(grepl("Germline", columns_x))){
        new_ord<-c(ord,stings_only)
      }else{
        new_ord<-c(ord,"Germline")
      }
      leg_names<-list(new_ord)
      leg_fill<-list(unname(colrs_sp[order(match(names(colrs_sp),new_ord))]))
    }else{
      leg_names<-list(names(colrs_sp))
      leg_fill<-list(unname(colrs_sp))
    }
    leg_border<-"black"
    leg_xpd<-TRUE
    inset<-c(-0.0005,0)
    if(opt=="isotype"){
      args.legend<-list(title="Isotypes")
    }else{
      args.legend<-list(title="Clusters")
    }
    title.cex<-1
    legend_param<-c(toString(leg_pos),toString(leg_bty),leg_names,leg_fill,toString(leg_border),inset,toString(leg_xpd),args.legend,title.cex)
    graph[[as.integer(parent.frame()$i)]]$progr_flag<-progr_flag
    graphs[[parent.frame()$i]]<-list("network" = graph[[as.integer(parent.frame()$i)]],"legend"=legend_param)
    #k<-graphs[[parent.frame()$i]]
    graphs[[parent.frame()$i]]
  })
  return(graphs)
}

PlotGraphs<-function(graphs,no_arg,topdf,filename){
  
  if(is.null((no_arg))){
    temp<-graphs[sapply(graphs,function(x) "network" %in% names(x))]
    if(length(temp)==0){
      temp<-lapply(graphs,function(x) x[sapply(x,function(p) "network" %in% names(p))])
    }
  }else{
    el_list<-sapply(graphs,function(x) x[no_arg])
    temp<-el_list[sapply(el_list,function(x) "network" %in% names(x))]
    if(length(el_list)==0){
      temp<-lapply(el_list,function(x) x[sapply(x,function(p) "network" %in% names(p))])
    }
  }
  
  if ("igraph" %in% unlist(lapply(temp,function(x) lapply(x,class))) || "igraph" %in% unlist(lapply(temp,function(x) lapply(x,function(a) lapply(a,class))))){
    
    if(topdf==TRUE){
      write_to_pdf(filename)
      lapply(temp, function(p){
        if(is.null(p$network)){
          p<-p[sapply(p,function(x) "network" %in% names(x))][[1]]
        }
        plot(p$network)
        if (!is.null(p$legend)){
          legend(p$legend[[1]],bty = p$legend[[2]],legend=p$legend[[3]],fill=p$legend[[4]],border=p$legend[[5]],inset=p$legend[[6]],xpd=p$legend[[7]],title=p$legend$title,cex=p$legend[[10]])
        }
        
        
      })
      
      dev.off()
    }else{
      lapply(temp, function(p){
        if(is.null(p$network)){
          p<-p[sapply(p,function(x) "network" %in% names(x))][[1]]
        }
        plot(p$network)
        if (!is.null(p$legend)){
          legend(p$legend[[1]],bty = p$legend[[2]],legend=p$legend[[3]],fill=p$legend[[4]],border=p$legend[[5]],inset=p$legend[[6]],xpd=p$legend[[7]],title=p$legend$title,cex=p$legend[[10]])
        }
        
        
      })
    }
  }else{
    if(topdf==TRUE){
      tryCatch(
        expr = {
          print_res<-lapply(graphs, function(x){
            x[[no_arg]]})
          
          gg_class<-c("gg","ggplot")
          #if(all(gg_class==unlist(lapply(print_res, class)))){
          if(all(match( gg_class, unlist(lapply(print_res, class))))){
            write_to_pdf(filename)
            invisible(lapply(print_res, print))
            dev.off()
          }else{
            stop('Error: Wrong arguments! Try again!')
          }
          
        },
        error = function(e){
          message('Error: Wrong arguments! Try again!')
        }
      )
      
      
    }else{
      tryCatch(
        expr = {
          print_res<-lapply(graphs, function(x){
            x[[no_arg]]
          })
        },
        error = function(e){
          message('Error: Wrong arguments! Try again!')
        }
      )
    }
  }
  
}

write_to_pdf<-function(number){
  pdf(width=15, height=10,file=paste0(number,".pdf"))
  
}
#save(Manuscript_LCMV_clonal_lineages,file='~/PHD/lcmv_10/Manuscript_LCMV_clonal_lineages.RData')

load('Manuscript_LCMV_clonal_lineages.RData')
testing_network <-Manuscript_LCMV_clonal_lineages[[1]][2:2]
### Now need to pick out the virus specific clones 
## could getthe barcodes of the specific ones 
####which Manuscript_LCMV_clonal_lineages has clonotypes 18,15,23,2,46,28
#16,14,23,2,33,27 ## 33 and 27 rows in original, but maybe earlier ones are merged?
## just to be sure, will re look up barcodes
# new_clonotype_index
## getting two barcodes from rows 27 and 33
repertoire_list[[1]]$barcodes[27] ## AACTGGTGTATGCTTG-1
repertoire_list[[1]]$barcodes[33] ## AACACGTGTGTGTGCC-1
for(i in 1:length(Manuscript_LCMV_clonal_lineages[[1]])){
  if(grepl(pattern = "AACACGTGTGTGTGCC-1",x = Manuscript_LCMV_clonal_lineages[[1]][[i]]$Name)==T)print(i)
}## 26 and 31

manuscript_lineage_ID <- c(16,14,23,2,31,26)
manuscript_lineages <-Manuscript_LCMV_clonal_lineages[[1]][c(16,14,23,2,31,26)]
network.function.seed1 <- VDJ_mutational_network(manuscript_lineages,random.seed=1,resolve.ties = "random")

clonal_frequency<-T # FALSE
scaleByClonalFreq<-TRUE
weight<-FALSE
dfs<-construct_structures(network.function.seed1[[6]],network.function.seed1[[7]],network.function.seed1[[14]],network.function.seed1[[10]],"isotype")
#dfs[1]
graphs<-create_graph(network.function.seed1[[1]],dfs,network.function.seed1[[2]],network.function.seed1[[11]],clonal_frequency,scaleByClonalFreq,weight,"isotype")

Manuscript_LCMV_clonal_lineages[[1]][[23]]$Name
## double check NP34, L28

network.function.seed1 <- VDJ_mutational_network(testing_network,random.seed=1,resolve.ties = "random")
network.function.seed3 <- VDJ_mutational_network(testing_network,random.seed=3,resolve.ties = "random")
plot(network.function.seed1[[1]][[1]])
clonal_frequency<-T # FALSE
scaleByClonalFreq<-TRUE
weight<-TRUE
topdf<-F # produce pdf, FALSE: see plot on console plots

#isotypes
dfs<-construct_structures(network.function.seed1[[6]],network.function.seed1[[7]],network.function.seed1[[14]],network.function.seed1[[10]],"isotype")
#dfs[1]
graphs<-create_graph(network.function.seed1[[1]],dfs,network.function.seed1[[2]],network.function.seed1[[11]],clonal_frequency,scaleByClonalFreq,weight,"isotype")


#clusters
dfs<-construct_structures(network.function.seed1[[6]],network.function.seed1[[7]],network.function.seed1[[15]],network.function.seed1[[10]],"cluster")
graphs<-create_graph(network.function.seed1[[1]],dfs,network.function.seed1[[2]],network.function.seed1[[11]],clonal_frequency,scaleByClonalFreq,weight,"cluster")

#For ploting a specific network e.g in position 2 of list
#PlotGraphs(graphs[2],NULL,topdf,"random_graphs_clon_sc")

#plots all networks in list
pdf("~/PHD/lcmv_10/submission/submission_figures/Network_Xale2_cluster.pdf",height=10,width=10)
PlotGraphs(graphs,NULL,topdf,"random_graphs_clon_sc_cluster")
dev.off()
#plot(network.function.seed1[[1]][[1]],vertex.label=1:33)
#plot(network.function.seed3[[1]][[1]],vertex.label=1:33)

network.function.seed1
clonal_frequency<-T # FALSE
scaleByClonalFreq<-TRUE
weight<-TRUE
topdf<-F # produce pdf, FALSE: see plot on console plots

#isotypes
dfs<-construct_structures(network.function.seed1[[6]],network.function.seed1[[7]],network.function.seed1[[14]],network.function.seed1[[10]],"isotype")
#dfs[1]
graphs<-create_graph(network.function.seed1[[1]],dfs,network.function.seed1[[2]],network.function.seed1[[11]],clonal_frequency,scaleByClonalFreq,weight,"isotype")
  
  
#clusters
dfs<-construct_structures(network.function.seed1[[6]],network.function.seed1[[7]],network.function.seed1[[15]],network.function.seed1[[10]],"cluster")
graphs<-create_graph(network.function.seed1[[1]],dfs,network.function.seed1[[2]],network.function.seed1[[11]],clonal_frequency,scaleByClonalFreq,weight,"cluster")
  
#For ploting a specific network e.g in position 2 of list
#PlotGraphs(graphs[2],NULL,topdf,"random_graphs_clon_sc")
  
#plots all networks in list
pdf("~/PHD/lcmv_10/submission/submission_figures/Networks_cluster_merged_PDF_2020_12_31.pdf",height=10,width=10)
PlotGraphs(graphs,NULL,topdf,"random_graphs_clon_sc_cluster")
dev.off()


## these should be correct now just need to change the colors of them. 