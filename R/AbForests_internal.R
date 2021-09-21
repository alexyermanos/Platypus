#Internal functions

.MAIN<-function(sequences,isotypes,count,distance_mat,clonal_frequency,scaleByClonalFreq,weight,tie_flag,scaleByclocloseness_metr,scaleBybetweenness,opt,random.seed,alg_opt,number.of.clusters){
  cells.per.network<-length(sequences)
  isotype.per.variant<-c()
  transcriptome.cluster.per.variant<-c()
  pre_data<-.PREPROCESSING(sequences,isotypes,distance_mat,count)
  new.variant.names<-rownames(pre_data$dist_mat)
  adj_mat<-.ADJ_MAT(pre_data$dist_mat)
  germline.index<-which(apply(pre_data$countData, 1, function(r) any(grepl("Germline",r))))
  if(alg_opt=="two-step"){
    edgeList<-.FIND_EDGES_v1(pre_data$dist_mat,pre_data$countData)
    edgeList<-.DELETE_MULTIPLE_PARENTS(edgeList$node_list,edgeList$edges_final,pre_data$dist_mat,pre_data$countData,tie_flag,weight,random.seed)
  }else{
    edgeList<-.EDGES(pre_data$dist_mat,adj_mat,germline.index,weight,random.seed)
  }
  graph<-.CREATE_GRAPH(edgeList$edges_final,edgeList$node_list,pre_data$countData,pre_data$dist_mat,clonal_frequency,scaleByClonalFreq,weight,scaleBybetweenness,scaleByclocloseness_metr,opt,number.of.clusters)
  if(opt=="cluster" | opt=="isotype"){
    pre_data$countData<-pre_data$countData[order(match(rownames(pre_data$countData),igraph::V(graph$graph)$order)), ]
    if(opt=="isotype"){
      isotype.per.variant<-pre_data$countData$item_per_category
      transcriptome.cluster.per.variant<-pre_data$countData$empty_per_category
    }else{
      transcriptome.cluster.per.variant<-pre_data$countData$item_per_category
      isotype.per.variant<-pre_data$countData$empty_per_category
    }
    variants.per.network<-nrow(pre_data$countData)
    variant.sequences<-pre_data$countData$sequences
    cells.per.variant<-pre_data$countData$count_uniq
    isotype.per.cell<-isotypes
    cell.indicies.per.variant<-list()
    list_temp<-list()
    for(j in 1:length(variant.sequences)){
      for(i in 1:length(sequences)){
        if (sequences[i]==variant.sequences[j]){
          list_temp <- c(list_temp,i)
        }
      }
      cell.indicies.per.variant[[j]]<-list_temp
      list_temp<-c()
    }
    if(alg_opt=="two-step"){
      adj.matrix<-igraph::as_adjacency_matrix(graph$graph,sparse=FALSE)
      idx <- which(adj.matrix != 0, arr.ind=TRUE)
      adj.matrix_int<-cbind(idx, adj.matrix[idx])
      max_row<-max(adj.matrix_int[,1])
      if(is.infinite(max_row)==TRUE & is.infinite(max_row)==TRUE){
        adj.matrix<-matrix(NA,nrow=1,ncol=1)
        colnames(adj.matrix) <- rownames(adj.matrix)<- paste("L",count,"_Germline",sep="")
      }else{
        max_column<-max(adj.matrix_int[,2])
        adj.matrix<-matrix(rep(0,max_row*max_column),nrow=max(max_row,max_column),ncol=max(max_row,max_column))
        for( i in 1:length(adj.matrix_int[,1])){
          adj.matrix[adj.matrix_int[,1][i],adj.matrix_int[,2][i]]<-adj.matrix_int[,3][i]
        }
        rownames(adj.matrix)<-colnames(adj.matrix)<-edgeList$node_list
        adj.matrix<-adj.matrix[order(match(rownames(adj.matrix),rownames(pre_data$countData))), ]
        adj.matrix<-adj.matrix[,order(match(colnames(adj.matrix),rownames(adj.matrix))) ]
        times_arg<-nrow(adj.matrix)
        rownames(adj.matrix) <- paste("L",rep(count,times_arg),"_","v",rownames(adj.matrix),"_","c",pre_data$countData$count_uniq,sep="")
        rownames(adj.matrix)[germline.index] <- paste("L",count,"_Germline",sep="")
        colnames(adj.matrix) <- rownames(adj.matrix)
        diag(adj.matrix)<-NA
      }
    }else{
      adj.matrix<-edgeList$adj_mat
      if(any(!(is.na(adj.matrix)))){
        rownames(adj.matrix)[germline.index] <- paste("L",count,"_Germline",sep="")
        colnames(adj.matrix) <- rownames(adj.matrix)
        diag(adj.matrix)<-NA
      }else{
        adj.matrix<-matrix(NA,nrow=1,ncol=1)
        colnames(adj.matrix) <- rownames(adj.matrix)<- paste("L",count,"_Germline",sep="")
      }

    }

  }else{
    adj.matrix<-NULL
    distance.matrix<-pre_data$dist_mat
    cells.per.network<-NULL
    variants.per.network<-NULL
    variant.sequences<-NULL
    cells.per.variant<-NULL
    cell.indicies.per.variant<-NULL
    new.variant.names<-NULL
    germline.index<-NULL
    isotype.per.variant<-NULL
    transcriptome.cluster.per.variant<-NULL
  }

  return(list("network" = graph$graph,"legend"=graph$leg,count.rand= edgeList$count_rand,adj.matrix=adj.matrix,distance.matrix=pre_data$dist_mat,cells.per.network=cells.per.network,variants.per.network=variants.per.network,variant.sequences=variant.sequences,cells.per.variant=cells.per.variant,cell.indicies.per.variant=cell.indicies.per.variant,new.variant.names=new.variant.names,germline.index=germline.index,isotype.per.variant=isotype.per.variant,transcriptome.cluster.per.variant=transcriptome.cluster.per.variant))
}

.PREPROCESSING<-function(sequences,isotypes,distance_mat,count){
  data <- cbind.data.frame(sequences, isotypes)
  data <- cbind(data,.ISOTYPE_ENCODING(data))
  countData<-.COUNT_DATA_INFO(data)
  if(is.null(distance_mat)){
    dist_mat<-.CREATE_DIST_MAT(countData,count)
  }else{
    dist_mat<-distance_mat
  }
  return(list(countData=countData,dist_mat=dist_mat))
}


.ISOTYPE_ENCODING<-function(data){
  id <- rownames(data)
  data <- cbind(id=id, data)
  one_hot_en <- stats::reshape(data,idvar="id",timevar="isotypes", direction="wide", v.names="isotypes", sep="")
  one_hot_en[which(is.na(one_hot_en), arr.ind=TRUE)] <- 20
  sel_col<-one_hot_en[grepl("isotypes", colnames(one_hot_en))]
  sel_col[which(sel_col != 20, arr.ind=TRUE)] <- 1
  sel_col[which(sel_col != 1, arr.ind=TRUE)] <- 0
  return(sel_col)
}


.COUNT_DATA_INFO<-function(data){

  data<-data[ , !(names(data) %in% c("isotypes"))]
  count_data_uni<-data %>%
    dplyr::group_by(data$sequences) %>%
    dplyr::summarise(count_uniq = dplyr::n(),.groups = 'drop')

  data[grepl("isotypes", colnames(data))][] <- lapply(data[grepl("isotypes", colnames(data))], as.numeric)
  count_data_sum<-stats::aggregate(. ~ sequences, data, sum)
  count_data_sum$list_item <- lapply(1:nrow(count_data_sum), function(x) unlist(unname(mapply(rep, names(count_data_sum[2:ncol(count_data_sum)]), unname(count_data_sum[x,2:ncol(count_data_sum)])))))
  count_data_sum$list_item<-count_data_sum$list_item[unlist(unname(lapply(count_data_sum$list_item,length))>0)]
  count_data_sum$list_item<-lapply(count_data_sum$list_item, lapply, function(z) gsub('^.{8}', '',z))
  item_per_category<-count_data_sum$list_item
  count_data_sum$list_item<-NULL
  count_data_sum$empty_item<-lapply(1:nrow(count_data_sum), function(x) unlist(unname(mapply(rep, rep("Unknown",ncol(count_data_sum)-1), unname(count_data_sum[x,2:ncol(count_data_sum)])))))
  empty_per_category<-count_data_sum$empty_item
  count_data_sum$empty_item<-NULL
  count_data_sum$sequences<-NULL
  count_data<-cbind(count_data_uni, count_data_sum)
  count_data<-as.data.frame(count_data)
  colnames(count_data)[1]<-gsub('^.{5}', '', colnames(count_data)[1])
  list1 <- list(count_data[,colnames(count_data[,3:ncol(count_data)])])
  list2 <- list(count_data[,colnames(count_data[2])])
  count_data[,paste(colnames(count_data[,3:ncol(count_data)]), "Ratio", sep="")] <- unlist(list1)/unlist(list2)
  count_data$Names<- 1:nrow(count_data)
  col_germ<-unname(which(sapply(names(count_data), function(x) any(grepl("germline", x)))))
  row_germ<-which(apply(count_data[,col_germ][1], 1, function(r) r==1))
  count_data$Names[row_germ]<-"Germline"
  count_data$SeqLen<-nchar(as.vector(count_data$sequences))
  count_data$item_per_category<-item_per_category
  count_data$empty_per_category<-empty_per_category
  return(count_data)
}


.CREATE_DIST_MAT <-function(data,count){
  dist_mat <- stringdist::stringdistmatrix(data$sequences,method="lv")
  dist_mat_df <- reshape2::melt(as.matrix(dist_mat), varnames = c("row", "col"))
  mat<-matrix(dist_mat_df$value, nrow=length(data$sequences), ncol=length(data$sequences))
  rownames(mat) <- paste("L",rep(count,nrow(data)),"_","v",data$Names,"_","c",data$count_uniq,sep="")
  rownames(mat)[which(apply(data, 1, function(r) any(grepl("Germline",r))))] <- paste("L",count,"_Germline",sep="")
  colnames(mat) <- rownames(mat)
  diag(mat) <- Inf
  return(mat)
}

.ADJ_MAT<-function(mat){
  adj.mat<- matrix(0,nrow=nrow(mat),ncol=ncol(mat))
  rownames(adj.mat) <- rownames(mat)
  colnames(adj.mat) <- rownames(adj.mat)
  return(adj.mat)
}

.REMOVE_DUPLICATES<-function(node_list){
  temp_list<-duplicated(node_list)
  node_list <- node_list[temp_list == FALSE]
  return(node_list)
}


.FIND_MIN<-function(mat){
  min_el<-min(mat[mat>0])
  pos<-which(mat == min(mat[mat>0]))
  return(list(pos=pos,min_el=min_el))
}


.IS_CONTAINED<-function(vec1,vec2){
  x=vector(length = length(vec1))
  for (i in 1:length(vec1)) {
    x[i] = vec1[i] %in% vec2
    if(length(which(vec1[i] %in% vec2)) == 0) vec2 else
      vec2=vec2[-match(vec1[i], vec2)]
  }
  y=(x==T)
  return(y)
}


.LIST_TO_DF<-function(edges){
  df <- data.frame(matrix(unlist(edges), nrow=length(edges),byrow=T),stringsAsFactors=FALSE)
  colnames(df) <- "Edges"
  edges_df <- data.frame(df,do.call(rbind,stringr::str_split(df$Edges,",")))
  edges_df<-edges_df[ , !(names(edges_df) %in% c("Edges"))]
  colnames(edges_df) <- c("From","To")
  return(edges_df)
}

.FIND_EDGES_v1<-function(mat,data){
  germ_row<-which(apply(data, 1, function(r) any(grepl("Germline",r))))
  node<-c(germ_row)
  level_next<-c()
  edges_final <- c()
  change_min<-c()
  node_list<-c()
  change_min_final<-c()
  l_next<-c()
  min_pos<-tryCatch(expr={
    .FIND_MIN(mat[node,])$pos
  },
  warning = function(w){
    return("NA")
  }
  )

  edges_final<-c(edges_final,apply(expand.grid(node, min_pos), 1, paste, collapse=","))
  if (all(ifelse(unique(min_pos) == unique(rep("NA",length(min_pos))), TRUE, FALSE))){
    fer<-unlist(strsplit(edges_final, ","))
    node_list<-c(node_list,as.integer(fer[1]))
    return(list(edges_final = edges_final,node_list = node_list))
  }
  mat[node,min_pos]<-0
  mat[min_pos,node]<-0
  mat[min_pos,min_pos]<-0
  node_list<-c(node,min_pos)
  i<-0
  while(length(node_list)<nrow(mat)){
    mat[node_list,node_list]<-0
    numbers_mat<-1:nrow(mat)
    numbers_mat<-numbers_mat[!(numbers_mat %in% node_list)]
    temp_io<-c()
    i<-i+1
    min_values<-c()
    node_sel<-c()
    for(j in seq_along(node_list)){
      combs<-expand.grid(node_list[j],numbers_mat,stringsAsFactors = FALSE)
      io<-.FIND_MIN(sapply(seq(nrow(combs)), function(i) mat[combs$Var1[i], combs$Var2[i]]))$pos
      io_min_el<-.FIND_MIN(sapply(seq(nrow(combs)), function(i) mat[combs$Var1[i], combs$Var2[i]]))$min_el
      temp_io<-c(temp_io,combs$Var2[io])
      min_values<-c(min_values,mat[node_list[j],combs$Var2[io]])
      node_sel<-c(node_sel,rep(j,length(io)))
    }
    io<-.FIND_MIN(min_values)$pos
    comma_el<-rep(",",length(node_list[node_sel[io]]))
    edges_final<-c(edges_final,paste0(node_list[node_sel[io]], comma_el,temp_io[io]))
    edges_final<-.REMOVE_DUPLICATES(edges_final)
    node_list<-c(node_list,temp_io[io])
    node_list<-.REMOVE_DUPLICATES(node_list)
  }

  return(list(edges_final = edges_final, node_list = unname(node_list),dist_mat=mat))
}


.EDGES<-function(dist_mat,adj_mat,germline.index,weight,random.seed){
  all.nodes <- 1:nrow(dist_mat)
  edges_final_final<-c()
  set.seed(random.seed,kind = "Mersenne-Twister")
  if(length(all.nodes)>1){
    current.nodes.in.network <- germline.index
    nodes.not.in.network <- all.nodes[-current.nodes.in.network]
    j <- 1
    count_rand<-0
    while(length(nodes.not.in.network)>0){
      if(j==1){
        x <- dist_mat[germline.index,]
        closest.nodes <- which(x==min(x))
        if(weight==FALSE) gl.dist <- 1
        else if(weight==TRUE) gl.dist <- min(dist_mat[germline.index,])
        for(k in 1:length(closest.nodes)){

          adj_mat[germline.index,closest.nodes[k]] <- gl.dist
          adj_mat[closest.nodes[k],germline.index] <- gl.dist

        }
        current.nodes.in.network <- unique(c(current.nodes.in.network,closest.nodes))
        nodes.not.in.network <- all.nodes[-current.nodes.in.network]
      }
      if(j>1){
        x <- dist_mat
        x[nodes.not.in.network,] <- Inf
        x[,current.nodes.in.network] <- Inf

        next.min.distance <- min(x)
        next.min.distance.indices <- which(x == min(x), arr.ind = TRUE)
        new.unique.newnodes <- unique(next.min.distance.indices[,2])
        unique.column.nodes <- (which(duplicated(next.min.distance.indices[,2])==F))

        new.unique.newnodes <- unique(next.min.distance.indices[,2])
        if(length(next.min.distance.indices[,2])==length(unique(next.min.distance.indices[,2]))){
          closest.nodes.to.current <- (next.min.distance.indices)
        }else{
          keep.row <- rep(NA,length(new.unique.newnodes))
          for(m in 1:length(new.unique.newnodes)){
            temp.which <- which(next.min.distance.indices[,2]==new.unique.newnodes[m])
            if(length(temp.which)>1){
              sample.which <- sample(temp.which,size = 1)
              count_rand<-count_rand+1
              keep.row[m] <- sample.which
            }
            else{
              keep.row[m] <- temp.which
            }
          }
          closest.nodes.to.current <- matrix(next.min.distance.indices[keep.row,],nrow=length(new.unique.newnodes))
        }
        if(weight==FALSE) next.min.distance <- 1
        else if(weight==TRUE) next.min.distance <- min(x)
        adj.row.index <- rep(NA,nrow(closest.nodes.to.current))
        adj.col.index <- rep(NA,nrow(closest.nodes.to.current))


        for(k in 1:nrow(closest.nodes.to.current)){
          adj.row.index[k] <- closest.nodes.to.current[k,1]
          adj.col.index[k] <- closest.nodes.to.current[k,2]
          adj_mat[adj.row.index[k],adj.col.index[k]] <- next.min.distance
          adj_mat[adj.col.index[k],adj.row.index[k]] <- next.min.distance
        }
        current.nodes.in.network <- unique(c(current.nodes.in.network,adj.col.index))
        nodes.not.in.network <- all.nodes[-current.nodes.in.network]
      }
      j <- j + 1
    }
    diag(adj_mat) <- NA
    net_obj<- igraph::graph_from_adjacency_matrix(adj_mat, mode = c("undirected"), weighted = NULL, diag = FALSE,add.colnames = NULL, add.rownames = NA)
    edges_final<-igraph::as_ids(igraph::E(net_obj))
    edges_final<-gsub("\\|", ",", edges_final)
    edges_final_1<-sapply(strsplit(edges_final, ","), "[", 1)
    edges_final_2<-sapply(strsplit(edges_final, ","), "[", 2)
    edges_final_1<-.NET_TO_COLON_FORMAT(edges_final_1,germline.index)
    edges_final_2<-.NET_TO_COLON_FORMAT(edges_final_2,germline.index)
    edges_final<-paste(edges_final_1,edges_final_2,sep=",")
    edges_final<-.REMOVE_DUPLICATES(edges_final)
    node_list<-current.nodes.in.network
    for(edge in edges_final){
      fer<-unlist(strsplit(edge, ","))
      pos1<-match(as.numeric(fer[1]),node_list)
      pos2<-match(as.numeric(fer[2]),node_list)
      if(pos1>pos2){
        edges_final_final<-c(edges_final_final,paste(fer[2],fer[1],sep=","))
      }else{
        edges_final_final<-c(edges_final_final,edge)
      }
    }
    adj_mat[is.na(adj_mat)] <- 0
    node_list<-as.numeric(node_list)
    node_list<-.REMOVE_DUPLICATES(node_list)
  }else{
    node_list <-all.nodes
    edges_final_final<-"1,1"
    adj_mat<-NA
    count_rand<-0
  }

  return(list(edges_final = edges_final_final, node_list =node_list ,adj_mat=adj_mat,count_rand=count_rand))
}

.NET_TO_COLON_FORMAT<-function(edges_final,germline.index){
  edges_final<-sub("_(.*)_.*", "\\1", edges_final)
  edges_final<-gsub(".*v","",edges_final)
  germ_pos<-which(sapply(edges_final, function(r) any(grepl("Germline",r))))
  edges_final[germ_pos]<-germline.index
  return(edges_final)
}

.DELETE_MULTIPLE_PARENTS<-function(nodes,edges,mat,data,flag,weight,random.seed){
  set.seed(random.seed,kind = "Mersenne-Twister")
  count_rand<-0
  for (node in nodes){
    if(is.null(node)) break
    match_level<-paste0(",",as.character(node),"$")
    found_in_edges<-data.table::like(edges, match_level,ignore.case = FALSE, fixed = FALSE)
    if (all(sum(found_in_edges, na.rm = TRUE))){
      pos_t<-c(which(found_in_edges %in% c(T)))
      values_l<-c()
      inter_edges<-c()
      match_min_el_mult<-c()
      for (min_el in pos_t){
        match_min_el<-gsub(",.*$", "", edges[min_el])
        match_min_el_mult<-c(match_min_el_mult,match_min_el)
        inter_edges<-c(inter_edges,edges[min_el])
      }
      if (length(inter_edges)>1){
        if(flag=='rand'){
          idx <- sample(match_min_el_mult,1)
          match_level<-paste0("^",idx,",")
          found_in_edges<-data.table::like(inter_edges,match_level,ignore.case = FALSE, fixed = FALSE)
          if (all(sum(found_in_edges, na.rm = TRUE))){
            idx<-c(which(found_in_edges %in% c(T)))
          }
          count_rand<-count_rand+1
        }else if (flag=='full'){
          next
        }else if(flag=="most_expanded" |flag=="least_expanded" ){
          w_t<-1:length(inter_edges)
          uni_seq<-unlist(lapply(match_min_el_mult, function(x) unlist(data[x,"count_uniq"])))
          if (flag=="most_expanded"){
            idx<-which(uni_seq == max(uni_seq[uni_seq>0]))
          }else{
            idx<-which(uni_seq == min(uni_seq[uni_seq>0]))
          }
          if( length(idx) > 1 ){
            idx <- sample(idx,1)
            match_level<-paste0("^",idx,",")
            found_in_edges<-data.table::like(inter_edges,match_level,ignore.case = FALSE, fixed = FALSE)
            if (all(sum(found_in_edges, na.rm = TRUE))){
              idx<-c(which(found_in_edges %in% c(T)))
            }
            count_rand<-count_rand+1
          }

        }else{
          first_half<-0
          paths<-list()
          j<-0
          for (el in match_min_el_mult){
            j<-j+1
            path_to_germ<-c()
            first_half<-el
            while (is.null(first_half)==FALSE) {
              match_el<-paste0(",",as.character(first_half),"$")
              if((length(match_el))>1){
                match_el<-paste(match_el,collapse="|")
              }
              found_in_list_edges<-grepl(match_el,edges)
              recur_list<-c(which(found_in_list_edges %in% c(T)))
              path_to_germ<-c(path_to_germ,edges[recur_list])
              first_half<- gsub(",.*$", "", as.character(edges[recur_list]))
              if (identical(first_half, character(0))){
                break
              }

            }
            paths[[j]] <- path_to_germ
          }
          path_size<-c()
          if (flag=='close_to_germ' | flag=='far_from_germ'){
            path_size<-rapply(paths, length, how="unlist")
          }else{
            weights<-.ADD_EDGE_WEIGHTS(.LIST_TO_DF(edges),weight,data,mat)
            df<-.LIST_TO_DF(edges)
            df$label <- paste(df$From, df$To,sep=",")
            for (k in 1:length(paths)){
              idx_row<-Reduce(union, lapply(paths[[k]], function(a) which(rowSums(df == a) > 0)))
              path_size<-c(path_size,sum(sapply(idx_row, function(z) weights[z])))
            }
          }

          if (flag=='close_to_germ' |flag=='close_path_to_germ'){
            idx<-which(path_size == min(path_size[path_size>0]))
          }else{
            idx<-which(path_size == max(path_size[path_size>0]))
          }

          if( length(idx) > 1 ){
            idx <- sample(idx,1)
            match_level<-paste0("^",idx,",")
            found_in_edges<-data.table::like(inter_edges,match_level,ignore.case = FALSE, fixed = FALSE)
            if (all(sum(found_in_edges, na.rm = TRUE))){
              idx<-c(which(found_in_edges %in% c(T)))
            }
            count_rand<-count_rand+1
          }
        }
        inter_edges<-inter_edges[inter_edges!=inter_edges[idx]]
        edges<-edges[!(edges %in% inter_edges)]
      }
    }

  }
  nodes<-.REMOVE_DUPLICATES(nodes)
  return(list(edges_final = edges, node_list = nodes,count_rand=count_rand))
}

.CREATE_GRAPH<-function(edges,nodes,data,distMat,clonal_frequency,scaleByClonalFreq,weight,scaleBybetweenness,scaleByclocloseness_metr,opt,number.of.clusters){
  progr_flag<-FALSE
  edges_df<-.LIST_TO_DF(edges)
  if(nrow(edges_df)>1){
    edges_df<-edges_df[!grepl("NA", edges_df$To),]
    rownames(edges_df)<-NULL
    nodes<-nodes[!is.na(nodes)]
  }
  if ("NA" %in% edges_df[,"To"]) {
    nodes<-append(nodes, "NA")
  }
  columns <- grep('Ratio', colnames(data), value=TRUE)
  columns_x<-gsub('.{5}$', '', columns)
  columns_x<-gsub('^.{8}', '',   columns_x)
  pos_germ<-which(grepl("germline",columns_x))
  colrs_sp<-c()
  if(opt=="isotype"){
    colors<-c("gray85","red","springgreen2","brown","white","yellow","purple","darkorange","pink","darkgreen")
    names(colors)<-c("germline","IGHA","IGHG1","IGHG2B","Unknown","IGHG2C","IGHM","IGHG3","IGHD","IGHE") #None

    for (co in names(colors)){
      if (any(grepl(co,columns))){
        colrs_sp<-c(colrs_sp,colors[co])
      }
    }
    colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
  }else{
    if(any(grepl(c("germline"),columns_x)) & any(grepl(c("Unknown"),columns_x))) {
      progr_flag<-TRUE
      first_cols<-c("gray85", "#A9A9A9")
      if(any(grepl("[A-Za-z]", number.of.clusters))==TRUE | length(number.of.clusters)==0){
        if(length(number.of.clusters)==0){
          colors<-first_cols
        }else{
          colors<-c(first_cols,scales::hue_pal()(length(number.of.clusters)))
        }
        names_temp<-number.of.clusters
        names(colors)<-c(c("germline","Unknown"),names_temp)
        colrs_sp<-base::subset(colors, names(colors) %in% columns_x)
        colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
      }else{
        colors<-c(first_cols,scales::hue_pal()(max(as.numeric(number.of.clusters))+1))
        names_temp<-as.character(c(0:(max(as.numeric(number.of.clusters)))))#as.character(c(1:30))
        names(colors)<-c(c("germline","Unknown"),names_temp)
        colrs_sp<-base::subset(colors, names(colors) %in% columns_x)
        colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
      }

    }else if (any(grepl("germline",columns_x))){
      first_cols<-c("gray85")
      if(any(grepl("#[A-Za-z]", number.of.clusters))==TRUE){
        colors<-columns_x
        pos_germ<-which(grepl("germline",columns_x))
        colors[pos_germ]<-gsub("^.*?germline"," ",colors[pos_germ])

        colors<-gsub(" ", "", colors, fixed = TRUE)
        names(colors)<-c(0:(length(colors)-1))
        names(colors)[pos_germ]<-"germline"
        colrs_sp<-colors
        colors<-append(colors, "#A9A9A9")
        names(colors)[length(colors)]<-"Unknown"


      }else if(any(grepl("[A-Za-z]", number.of.clusters))==TRUE){
        colors<-c(first_cols,scales::hue_pal()(length(number.of.clusters)))
        names_temp<-number.of.clusters
        names(colors)<-c("germline",names_temp)
        colrs_sp<-base::subset(colors, names(colors) %in% columns_x)
        colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
        colors<-append(colors, "#A9A9A9", after=1)
        names(colors)[2]<-"Unknown"

      }else{
        colors<-c(first_cols,scales::hue_pal()(max(as.numeric(number.of.clusters))+1))
        names_temp<-as.character(c(0:(max(as.numeric(number.of.clusters)))))
        names(colors)<-c("germline",names_temp)
        colrs_sp<-base::subset(colors, names(colors) %in% columns_x)
        colrs_sp<-colrs_sp[order(match(names(colrs_sp),columns_x))]
        colors<-append(colors, "#A9A9A9", after=1)
        names(colors)[2]<-"Unknown"
      }

    }


  }
  difer<-length(colrs_sp)-length(columns_x)
  if(length(colrs_sp)>length(columns_x)){
    colrs_sp<-utils::head(colrs_sp, -difer)
  }
  sequences <- lapply(nodes, function(x) unlist(data[x,1]))
  values <- lapply(nodes, function(x) unlist(data[x,columns]))

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
  if (!"NA" %in% edges_df[,"To"]) {
    g <- igraph::graph_from_data_frame(edges_df, vertices = nodes)
    g<-igraph::simplify(g)
    g$edges_list<-edges_df
    igraph::E(g)$weight <-.ADD_EDGE_WEIGHTS(edges_df,weight,data,distMat)
    igraph::E(g)$color <- "black"
  }else{
    pos_na<-which(nodes %in% "NA")
    g <- igraph::graph_from_data_frame(edges_df, vertices = nodes)
    edges<-gsub(",", "|",edges)
    g<-igraph::delete_edges(g,edges)
    g<- igraph::delete_vertices(g,pos_na )
    nodes<-nodes[-pos_na]
    shapes<-shapes[-pos_na]
    col_mat<-col_mat[-pos_na]
    values<-values[-pos_na]
    colrs<-colrs[-pos_na]

  }
  if(clonal_frequency==TRUE & (scaleBybetweenness==TRUE | scaleByclocloseness_metr==TRUE )){
    uniq_seq <- lapply(nodes, function(x) unlist(data[x,"count_uniq"]))
    igraph::V(g)$name<- uniq_seq
    igraph::V(g)$scale<- nodes
    igraph::V(g)$label<-uniq_seq
    igraph::V(g)$order<-nodes
    igraph::V(g)$clone_freq<-uniq_seq
  }else if(clonal_frequency==TRUE & scaleByClonalFreq==TRUE){
    uniq_seq <- lapply(nodes, function(x) unlist(data[x,"count_uniq"]))
    igraph::V(g)$name<-uniq_seq
    igraph::V(g)$scale<- uniq_seq
    igraph::V(g)$label<-uniq_seq
    igraph::V(g)$order<-nodes
    igraph::V(g)$clone_freq<-uniq_seq
  }else if(scaleByClonalFreq==TRUE){
    uniq_seq <- lapply(nodes, function(x) unlist(data[x,"count_uniq"]))
    igraph::V(g)$name<-nodes
    igraph::V(g)$scale<- uniq_seq
    igraph::V(g)$label<-nodes
    igraph::V(g)$order<-nodes
    igraph::V(g)$clone_freq<-uniq_seq
  }else if(clonal_frequency==TRUE){
    uniq_seq <- lapply(nodes, function(x) unlist(data[x,"count_uniq"]))
    igraph::V(g)$name<-uniq_seq
    igraph::V(g)$scale<- nodes
    igraph::V(g)$clone_freq<-uniq_seq
    igraph::V(g)$label<-uniq_seq
    igraph::V(g)$order<-nodes
  }else{
    uniq_seq <- lapply(nodes, function(x) unlist(data[x,"count_uniq"]))
    igraph::V(g)$name<-nodes
    igraph::V(g)$scale<-nodes
    igraph::V(g)$label<-nodes
    igraph::V(g)$order<-nodes
    igraph::V(g)$clone_freq<-uniq_seq
  }
  igraph::V(g)$size<-.VERTEX_SCALE(g,igraph::E(g),data,igraph::V(g)$scale,igraph::E(g)$weight,scaleByClonalFreq,scaleBybetweenness,scaleByclocloseness_metr)$node_size
  igraph::V(g)$shape<-unlist(shapes)
  igraph::V(g)$color<-unlist(col_mat)
  igraph::V(g)$pie<-values
  igraph::V(g)$pie.color<-colrs
  igraph::V(g)$seq<-sequences
  igraph::V(g)$arrow.size<-0.03
  g$object_colors<-colrs_sp
  g$all_colors<-colors
  igraph::V(g)$edge.label.font <-0.5
  igraph::V(g)$edge.label.cex <- 1
  igraph::V(g)$edge.arrow.size <- 1
  igraph::V(g)$label.cex <- .VERTEX_SCALE(g,igraph::E(g),data,igraph::V(g)$scale,igraph::E(g)$weight,scaleByClonalFreq,scaleBybetweenness,scaleByclocloseness_metr)$label_size
  igraph::V(g)$label.color<-"black"
  coords<-igraph::layout_as_tree(g,root = igraph::V(g)[1], circular = FALSE, rootlevel = numeric(), mode = "all", flip.y = TRUE)
  g$layout<-coords #
  leg_pos<-"topright"
  leg_bty<-"n"
  if(!(opt=="isotype") & progr_flag){
    stings_only<-c("Unknown","germline")
    n_cols <- names(colrs_sp)[! names(colrs_sp) %in% stings_only]
    ord<-as.character(sort(as.numeric(n_cols), decreasing=F))
    if(any(grepl("Unknown", columns_x)) & any(grepl("germline", columns_x))){
      new_ord<-c(ord,stings_only)
    }else{
      new_ord<-c(ord,"germline")
    }
    leg_names<- list(names(colors))
    leg_fill<-list(colors)
    }else{
    leg_names<- list(names(colors))
    leg_fill<-list(colors)
  }
  leg_border<-"black"
  leg_xpd<-TRUE
  inset<-c(-0.0005,0)
  if(opt=="isotype"){
    args.legend<-list(title="Isotypes")
  }else if (opt=="cluster"){
    args.legend<-list(title="Transcriptional cluster")
  }else{
    args.legend<-list(title=opt)
  }
  title.cex<-1.5
  legend_param<-c(toString(leg_pos),toString(leg_bty),leg_names,leg_fill,toString(leg_border),inset,toString(leg_xpd),args.legend,title.cex)
  g$progr_flag<-progr_flag

  return(list(graph=g,legend=legend_param))

}


.VERTEX_SCALE<-function(g,edges,data, uniq_seq,weights,scaleByClonalFreq,scaleBybetweenness,scaleByclocloseness_metr){

  if(scaleByClonalFreq==TRUE){
    node_size<-mapply(as.numeric(uniq_seq), rep(2,length(uniq_seq)), FUN = function(x, y) {8+y*x/10})
  }else if(scaleBybetweenness==TRUE && !(any(is.na(.CENTRALITY_METRICS(g,uniq_seq,edges)$betweenessCentrality))==TRUE)){
    node_size<-2.5+.CENTRALITY_METRICS(g,uniq_seq,edges)$betweenessCentrality/30
  }else if (scaleByclocloseness_metr==TRUE && !(any(is.na(.CENTRALITY_METRICS(g,uniq_seq,edges)$clocloseness_metr))==TRUE)){
    node_size<-20+.CENTRALITY_METRICS(g,uniq_seq,edges)$clocloseness_metr/30
  }else{
    node_size<-6+0.1*nrow(data)
  }

  label_size <- 0.15*node_size
  return(list(node_size = node_size, label_size = label_size))
}

.ADD_EDGE_WEIGHTS<-function(edges,weight,data,mat){
  edges<-edges[!grepl("NA",edges$To),]
  edges$From <- as.numeric(as.character(edges$From))
  edges$To <- as.numeric(as.character(edges$To))
  weights<-c()
  for (j in 1:nrow(edges)){
    row<-edges[j,]$From
    col<-edges[j,]$To
    if (length(row==which(data$Names=="Germline", arr.ind=TRUE))>0| length(row==which(data$Names=="germline", arr.ind=TRUE))>0){
      if (weight==TRUE){
        mat[row,col]<-mat[row,col]-abs(data$SeqLen[row]-data$SeqLen[col])
      }
      else{
        mat[row,col]<-1
      }
    }
    weights<-c(weights,mat[unlist(row),unlist(col)])

  }

  return(weights)
}


.LONGEST_PATH_FROM_GERMLINE<- function(graph,type,weights){
  shortest.paths<-igraph::all_shortest_paths(graph, from=igraph::V(graph), to = igraph::V(graph), mode ="out", weights = weights)
  if(type=="weighted"){
    score<-igraph::distances(graph, v = igraph::V(graph), to = igraph::V(graph), mode = "out", weights = weights)
  }else{
    score<-igraph::distances(graph, v = igraph::V(graph), to = igraph::V(graph), mode = "out", weights = NA)
  }
  score_idx<- unname(score[1,])
  max_pos<-which(score_idx == max(score_idx), arr.ind=TRUE)
  if (length(max_pos)>1){
    paths<-c()
    max_scores<-c()
    i<-0
    for (j in max_pos){
      i<-i+1
      paths[[i]]<-shortest.paths$res[[j]]
      max_scores[[i]] <-score_idx[j]
    }

    return(list(shortest_paths = paths, score =max_scores))
  }
  return(list(shortest_paths = shortest.paths$res[[max_pos]], score =score_idx[max_pos]))
}


.COLOR_PATH<-function(graph,path,colrs_sp,opt,flag){

  if(is.vector(path)){
    colors_path<- scales::hue_pal()(30)
    for (i in 1:length(path)){
      igraph::E(graph, path=path[[i]])$color<-colors_path[i]
    }
  }else{
    igraph::E(graph, path=path)$color <- "magenta"
  }

  leg_pos<-"topright"
  leg_bty<-"n"
  if(!(opt=="isotype") & flag ){
    stings_only<-c("Unknown","germline")
    n_cols <- names(colrs_sp)[! names(colrs_sp) %in% stings_only]
    ord<-as.character(sort(as.numeric(n_cols), decreasing=F))
    new_ord<-c(ord,stings_only)
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
  return(list("network"=graph,legend=legend_param))
}



.AVERAGE_DAUGHTERS<-function(g,nodes){

  degree<-igraph::degree(g, v = as.integer(nodes), mode = "out",loops = FALSE, normalized = FALSE)
  avg_daughters<-sum(degree)/length(nodes)
  std_daughters<-stats::sd(degree)
  min_daughters<-min(degree)
  max_daughters<-max(degree)

  df<-as.data.frame(table(degree))
  p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'degree', y = 'Freq', fill= 'degree')) + ggplot2::xlab("Vertex degree") + ggplot2::ylab("Frequency")+ ggplot2::geom_bar(width=0.4,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+ #ggtitle("Degree Distribution of daughter cells")+
    ggplot2::coord_cartesian()+
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add=c(0.6,0.5)))+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(2), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
          panel.background = ggplot2::element_blank(),legend.title=ggplot2::element_text(size=12,face="bold")) + ggplot2::scale_fill_discrete(name = "Vertex Degrees")#+ scale_y_continuous(expand = c(0, 0), limits = c(0, 20))



  return(list(avg_daughters=avg_daughters,std_daughters=std_daughters,min_daughters=min_daughters,max_daughters=max_daughters,p=p))
}

.CLUSTER_COEFFICIENTS<-function(g){
  global_cluster_coef<-igraph::transitivity(g, type = "local",vids=NULL)
  avg_cluster_coef<-igraph::transitivity(g, type = "average")
  return(list(global_cluster_coef=global_cluster_coef,avg_cluster_coef=avg_cluster_coef))
}


.GRAPH_STRENGTH_DISTRIBUTION<- function(graph,weights){
  if (!igraph::is.igraph(graph)) {
    stop("Not a graph object")
  }
  edge_strengths <- igraph::graph.strength(graph,mode = "out",loops=FALSE,weights=weights)
  df<-as.data.frame(table(edge_strengths))
  p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'edge_strengths', y = 'Freq',fill= edge_strengths)) + ggplot2::xlab("Weighted Vertex degree") + ggplot2::ylab("Frequency")+ggplot2::geom_bar(width=0.4,stat="identity", position="dodge",colour = "black",show.legend = FALSE)+  #+ ggtitle("Weighted Degree Distribution of daughter cells")
    ggplot2::coord_cartesian()+
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add=c(0.6,0.5)))+
    ggplot2::scale_y_continuous(expand = c(0,0)) +
   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(2), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5,colour = "black",face="bold"),
          panel.background = ggplot2::element_blank(),legend.title=ggplot2::element_text(size=12,face="bold")) +ggplot2::scale_fill_discrete(name = "Edge Strengths")#+ scale_y_continuous(expand = c(0, 0), limits = c(0, 20))


  return(list(edge_strengths=edge_strengths,p=p))
}


.CENTRALITY_METRICS<-function(g,nodes,edges){
  if(length(edges)==0){
    betweenessCentrality <- igraph::betweenness(g,as.integer(nodes), directed = TRUE, weights = NULL, nobigint = FALSE)
    edgeCentrality<-NA
    clocloseness_metr<-igraph::closeness(g,as.integer(nodes), mode = c("all"), weights = NULL, normalized = TRUE)
  }else{
    if (sum(edges$weight>0)==length(edges$weight)){
      betweenessCentrality <- igraph::betweenness(g,as.integer(nodes), directed = TRUE, weights = edges$weight, nobigint = FALSE)
      edgeCentrality<-igraph::edge_betweenness(g, edges, directed = TRUE, weights = edges$weight)
      clocloseness_metr<-igraph::closeness(g,as.integer(nodes), mode = c("all"), weights = edges$weight, normalized = TRUE)
    }else{
      betweenessCentrality<-NA
      edgeCentrality<-NA
      clocloseness_metr<-NA
    }
  }

  return(list(betweenessCentrality=betweenessCentrality,edgeCentrality=edgeCentrality,clocloseness_metr=clocloseness_metr))
}


.AVERAGE_ISOTYPES<-function(g,nodes,colrs_sp,opt){
  isotype_table<-do.call("rbind", nodes$pie)
  count_isotypes<-rep(0,ncol(isotype_table))

  for(i in 1:nrow(isotype_table)){
    count_isotypes[which(isotype_table[i,]!=0,arr.ind = T)]<-count_isotypes[which(isotype_table[i,]!=0,arr.ind = T)]+1
  }
  names(count_isotypes)<-colnames(isotype_table)
  names(count_isotypes)<-gsub('.{5}$', '', names(count_isotypes))
  names(count_isotypes)<-gsub('^.{8}', '',   names(count_isotypes))
  if(!(opt=="isotype")){
    names(count_isotypes)<-names(colrs_sp)
  }
  avg_isotypes<-count_isotypes/length(nodes)
  df <- data.frame(Combinations=names(count_isotypes), Freq=count_isotypes)
  df<-cbind(df,colrs_sp)
  df<-df[match(names(colrs_sp), df$Combinations),]
  df[] <- lapply( df, factor)
  df$colrs_sp <- as.factor(df$colrs_sp)

  if(opt=="isotype"){
    p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Combinations', y = 'Freq', fill='colrs_sp')) + ggplot2::xlab("")+ ggplot2::ylab("Frequency")+ggplot2::geom_bar(stat="identity", position="dodge",colour="black",width=0.4,show.legend = FALSE)+ #+ ggplot2::xlab(" Isotypes")+ ggtitle("Distribution of isotypes")
      ggplot2::coord_cartesian()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 25, vjust = 0.5,colour = "black",face="bold"),
         axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
            axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 60),
            axis.text = ggplot2::element_text(size = 14),
            panel.background = ggplot2::element_blank(), legend.title=ggplot2::element_text(size=12,face="bold"))+ ggplot2::scale_fill_manual(values = levels(df$colrs_sp))
  }else{
    p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Combinations', y = 'Freq', fill='colrs_sp')) + ggplot2::xlab("")+ ggplot2::ylab("Frequency")+ ggplot2::geom_bar(stat="identity", position="dodge",colour="black",width=0.4,show.legend = FALSE)+ #+ ggplot2::xlab(" Clusters")ggtitle("Distribution of Transcriptional Clusters")+
      ggplot2::coord_cartesian()+
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 25, vjust = 0.5,colour = "black",face="bold"),
          axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
            axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
            axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 60),
            axis.text = ggplot2::element_text(size = 14),
            panel.background = ggplot2::element_blank(), legend.title=ggplot2::element_text(size=12,face="bold"))+ ggplot2::scale_fill_manual(values = levels(df$colrs_sp))
  }

  return(list(avg_isotypes=avg_isotypes,p=p))

}


.GENERATE_COMBINATIONS<-function(x,g,isotype_info){
  combinations<- expand.grid(x,x,stringsAsFactors = FALSE)
  replacement<-names(isotype_info)
  colnames(combinations)<-c('Parent','Child')

  return((list(combinations = combinations, replacement = replacement)))
}


.ISOTYPE_POSITION<-function(g,nodes,list_edges,isotype_info,opt){

  isotype_table<-do.call("rbind", igraph::V(g)$pie)
  isotype_table <- cbind(isotype_table, Names=nodes)
  isotype_table<-as.data.frame(isotype_table)
  isotype_table <- isotype_table %>%
    dplyr::arrange(isotype_table$Names)

  isotype_table$Names<-NULL
  x<-names(isotype_table)[]
  x<-gsub('.{5}$', '', x)
  x<-gsub('^.{8}', '',   x)
  x<-unlist(lapply(x, function(z){ z[!is.na(z) & z != ""]}))
  first_char <- substr(x, start=1, stop=1)
  split <- strsplit(first_char, split="")
  if(any(split=="#")){
    isotype_info<- gsub("^.{0,1}", "",isotype_info )
  }
  if(!(opt=="isotype")){
    x_temp<-c()
    for (co in names(isotype_info)){
      if (any(grepl(co,x))){
        x_temp<-c(x_temp,names(isotype_info)[which(grepl(co,x))])
      }
    }
    x<-x_temp
  }
  x<-gsub('.{5}$', '', colnames(isotype_table))
  x<-gsub('^.{8}', '',  x)
  combs_x<-.GENERATE_COMBINATIONS(x,g,isotype_info)
  combs_x$replacement<-combs_x$replacement[order(match(combs_x$replacement,x))]
  combs<-combs_x$combinations
  par<-c()
  if(opt=="cluster" & any(grepl(paste0(x, collapse="|"),combs)) ){
    i<-0
    for (parent in x){
      i<-i+1
      if(length(which(combs$Parent == parent))>0){
        pos<-which(combs$Parent == parent)
        par[pos]<-c(combs_x$replacement[i])[1]
        }
    }
    combs$Parent<-par
    i<-0
    par<-c()
    for (parent in x){
      i<-i+1
      if(length(which(combs$Child == parent))>0){
        pos<-which(combs$Child == parent)
        par[pos]<-c(combs_x$replacement[i])[1]
      }
    }
    combs$Child<-par
  }
  combs<-dplyr::filter(combs, combs$Child != "germline")
  colnames(isotype_table)<-x
  count_list<-c()
  if(nrow(list_edges)>1){
    for (m in 1:nrow(combs)){
      tmp_P<-combs$Parent[m]
      tmp_C<-combs$Child[m]
      count<-0

      for(i in 1:nrow(list_edges)){
        for(j in 1:ncol(list_edges)){
          k<-list_edges[i,j]
        }
        temp<-t(sapply(list_edges[i,], function(k) isotype_table[as.numeric(k),]))
        combs_x$replacement<-combs_x$replacement[order(match(combs_x$replacement,x))]
        colnames(temp)<-combs_x$replacement
        temp<-data.table::data.table(temp)
        temp<-as.data.frame(temp)
        if(temp[[tmp_P]][1]>0){
          if(temp[[tmp_C]][2]>0){
            count<-count+1

          }
        }
      }
      count_list<-c(count_list,count)
    }
    combs$label <- paste(combs$Parent, "-", combs$Child)
    combs$Freq<-count_list
    combs<-combs[combs$Freq != 0, ]
    combs<- unique(combs)
    df <- data.frame(Combinations=combs$label, freq=combs$Freq)
    if(opt=="isotype"){
      p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Combinations', y = 'freq',fill='Combinations' )) + ggplot2::xlab("")+ ggplot2::ylab("Frequency")+ggplot2::geom_bar(stat="identity", position="dodge",colour = "black",show.legend = FALSE)+ #+ + ggplot2::xlab("Isotype Combinations") ggtitle("Isotype Directionality")
        ggplot2::coord_cartesian()+
        ggplot2::scale_y_continuous(expand = c(0,0)) +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,colour = "black",face="bold"),
              axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
              axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
              axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 60),
              axis.text = ggplot2::element_text(size = 12),
              panel.background = ggplot2::element_blank(), legend.title=ggplot2::element_text(size=12,face="bold")
        )

    }else{
      p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Combinations', y = 'freq',fill='Combinations' )) + ggplot2::xlab("")+ ggplot2::ylab("Frequency")+ggplot2::geom_bar(stat="identity", position="dodge",colour = "black",show.legend = FALSE)+  # xlab("Clusters Combinations") ++ ggtitle("Transcriptional Clusters Directionality")
        ggplot2::coord_cartesian()+
        ggplot2::scale_y_continuous(expand = c(0,0)) +
       ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.2,colour = "black",face="bold"),
             axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
              axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
              axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 60),
              axis.text = ggplot2::element_text(size = 12),
              panel.background = ggplot2::element_blank(), legend.title=ggplot2::element_text(size=12,face="bold")
        )

    }
  }else{
    combs <- combs[combs$Parent %in% combs_x$replacement, ]
    combs$Child<-"NA"
    combs<-unique(combs)
    p<-NULL
  }

  return((list(combs = combs,p=p)))

}


.ISOTYPE_TRANSITIONS<-function(combs){

  if('Freq' %in% names(combs)){
    colnames(combs)[colnames(combs) == 'Freq'] <- 'thickness'
    g <- igraph::graph_from_data_frame(combs)
    igraph::V(g)$label.cex <- 0.4#0.4
    igraph::V(g)$color<-"white"
    igraph::E(g)$width <- 0.1*combs$thickness + 0.5
    igraph::E(g)$arrow.width <- 0.05*combs$thickness +0.5
    igraph::E(g)$color<-"magenta"
    igraph::E(g)$arrow.size<-0.4
    igraph::E(g)$label<-combs$thickness
    el <- t(apply(igraph::get.edgelist(g),1,sort))
    curve<-0.3
    igraph::E(g)$curved <- 0
    igraph::E(g)[duplicated(el) | duplicated(el,fromLast =TRUE)]$curved <- curve
    coords<-igraph::layout.circle(g)
    g$layout<-coords*0.5
  }else{
    g <- igraph::make_empty_graph(directed=TRUE)
    g <- igraph::add.vertices(g, 1)
    igraph::V(g)$label.cex <- 0.4
    igraph::V(g)$color<-"white"
    igraph::V(g)$frame.color<-"blue"
    order_labels<-list(paste(combs$Parent,sep="\\n",collapse = "\n"))
    igraph::V(g)$name <- order_labels
    coords<-igraph::layout.circle(g)
    g$layout<-coords
  }
  return(list("network"=g))
}



.CLONAL_FREQUENCY_EDGES_FROM_GERM<-function(g,weights){
  shortest.paths<-igraph::all_shortest_paths(g, from=igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  count_table<- t(sapply(shortest.paths, lengths))
  count_table<-count_table[1,]

  no_edges_from_germ<-Map("-", count_table, rep(1,length(count_table) ))
  no_edges_from_germ<-unlist(no_edges_from_germ, use.names=FALSE)

  mean_no_edges_from_germ<-mean(no_edges_from_germ)


  ratio<-no_edges_from_germ/unlist(igraph::V(g)$clone_freq)

  node<-utils::head(names(rapply(shortest.paths, function(x) utils::tail(x, 1))),-1)
  node<-as.numeric(gsub("[^.]*.(.*)", "\\1", node))
  df <- data.frame(Edges_from_germ=no_edges_from_germ, Clonal_freq=unlist(igraph::V(g)$clone_freq, use.names=FALSE))
  p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Edges_from_germ', y = 'Clonal_freq'))+ ggplot2::geom_point(colour = "blue") + ggplot2::xlab("No of Edges from Germline") + ggplot2::ylab("Clonal Frequency")+ #+ ggtitle("Number of edges from Germline vs Clonal frequency")
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          panel.background = ggplot2::element_blank())


  return(list(ratio=ratio,mean_no_edges_from_germ=mean_no_edges_from_germ,p=p))
}


.CLONAL_FREQUENCY_PATH_FROM_GERM<-function(g,weights){
  shortest.paths<-igraph::all_shortest_paths(g, from=igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  score<-igraph::distances(g, v = igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  score_idx<- unname(score[1,])
  mean_score_idx<-mean(score_idx)

  ratio<-score_idx/unlist(igraph::V(g)$clone_freq)

  df <- data.frame(Path_length=score_idx, Clonal_freq=unlist(igraph::V(g)$clone_freq, use.names=FALSE))
  p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Path_length', y = 'Clonal_freq'))+ ggplot2::geom_point(colour = "blue") + ggplot2::xlab("Path length") + ggplot2::ylab("Clonal Frequency") + #+ ggtitle("Total Path length from Germline vs Clonal frequency")
   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(3), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          panel.background = ggplot2::element_blank())


  return(list(ratio=ratio,mean_score_idx=mean_score_idx,p=p))

}


.NODE_DEGREE_EDGES_FROM_GERM<-function(g,weights){
  node_degree<-igraph::degree(g, v = igraph::V(g), mode = c("out"),loops = FALSE, normalized = FALSE)
  shortest.paths<-igraph::all_shortest_paths(g, from=igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  count_table<- t(sapply(shortest.paths, lengths))
  count_table<-count_table[1,]
  no_edges_from_gem<-Map("-", count_table, rep(1,length(count_table) ))
  no_edges_from_gem<-unlist(no_edges_from_gem, use.names=FALSE)
  node<-utils::head(names(rapply(shortest.paths, function(x) utils::tail(x, 1))),-1)
  node<-as.numeric(gsub("[^.]*.(.*)", "\\1", node))
  df <- data.frame(Edges_from_germ=no_edges_from_gem,node_degree=unname(node_degree))
  p<-ggplot2::ggplot(df, ggplot2::aes_string(x = 'Edges_from_germ', y = 'node_degree', label='node'))+ ggplot2::geom_point(colour = "blue") + ggplot2::xlab("No of Edges from Germline") + ggplot2::ylab("Node Degree") +
   ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(2), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          panel.background = ggplot2::element_blank())


  return(p)
}


.NODE_DEGREE_PATH_FROM_GERM<-function(g,weights){
  node_degree<-igraph::degree(g, v = igraph::V(g), mode = c("out"),loops = FALSE, normalized = FALSE)
  shortest.paths<-igraph::all_shortest_paths(g, from=igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  score<-igraph::distances(g, v = igraph::V(g), to = igraph::V(g), mode = c("out"), weights = weights)
  score_idx<- unname(score[1,])

  node<-utils::head(names(rapply(shortest.paths, function(x) utils::tail(x, 1))),-1)
  node<-as.numeric(gsub("[^.]*.(.*)", "\\1", node))

  df <- data.frame(Path_length=score_idx,node_degree=unname(node_degree))
  p<- ggplot2::ggplot(df, ggplot2::aes_string(x = 'Path_length', y = 'node_degree', label='node'))+ ggplot2::geom_point(colour = "blue") + ggplot2::xlab("Path length") + ggplot2::ylab("Node Degree") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = ggplot2::element_line(colour = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.8,colour = "black",face="bold"),
          axis.title.y = ggplot2::element_text(size = ggplot2::rel(2), angle = 90),
          axis.title.x = ggplot2::element_text(size = ggplot2::rel(2), angle = 0),
          axis.text = ggplot2::element_text(size = 15),
          panel.background = ggplot2::element_blank())

  return(p)
}


.WRITE_TO_PDF<-function(number){
  grDevices::pdf(width=15, height=12,file=paste0(number,".pdf"))
}


.CREATE_DF<-function(list,metrics,number){
  nets<-sapply(list,function(x) x[1])
  allData<-data.frame(do.call(rbind, metrics))
  allData<-allData[, -grep("Plot", colnames(allData))]
  allData<-allData[, -grep("network", colnames(allData))]
  allData$Repertoire.id<-as.list(rep(paste0("",number),length(list)))
  allData$Number.of.sequences<-lapply(nets, function(x) length(igraph::V(x)))
  return(allData)
}

.FIND_MATCHES<-function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))

  ans = rep_len(as.character(fill), length(x))
  empty = seq_along(x)

  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]
    empty = empty[!greps]
  }

  return(ans)
}

.REFERENCE_INFO<-function() return(list(type_iso=c("Germline","IGHA","IGHG1","IGHG2B","IGHG2C","IGHM","IGHG3","IGHD","IGHE","Unknown")))

.MODIFY_INPUT<-function(k,type_iso,cdr3,opt){

  cont<-grepl(paste(c("barcode","sample_id","group_id","clonotype_id_10x","celltype","Nr_of_VDJ_chains","Nr_of_VJ_chains","VDJ_cdr3s_aa","VJ_cdr3s_aa","VDJ_cdr3s_nt",
                      "VJ_cdr3s_nt","VDJ_chain_contig","VJ_chain_contig","VDJ_chain",
                      "VJ_chain","VDJ_vgene","VJ_vgene","VDJ_dgene","VDJ_jgene","VJ_jgene","VDJ_cgene","VJ_cgene","VDJ_sequence_nt_raw","VJ_sequence_nt_raw",
                      "VDJ_sequence_nt_trimmed","VJ_sequence_nt_trimmed","VDJ_sequence_aa","VJ_sequence_aa","VDJ_trimmed_ref","VJ_trimmed_ref","VDJ_raw_consensus_id",
                      "VJ_raw_consensus_id","orig_barcode.x","clonotype_frequency","GEX_available","orig.ident","orig_barcode.y","seurat_clusters","PC_1","PC_2","UMAP_1",
                      "UMAP_2","tSNE_1","tSNE_2","clonotype_id","new_clonal_feature","new_clonal_frequency","new_clonal_rank","new_tied_clonal_rank"), collapse = "|"), names(k))

  if(any(cont)==TRUE){
    k<-subset(k, k$Nr_of_VDJ_chains ==1 & k$Nr_of_VJ_chains==1)
    if(length(unique(unlist(unname(lengths(k)))))==0){
      return(NULL)
    }
    k$ref<-paste(k$VDJ_trimmed_ref,k$VJ_trimmed_ref,sep='')
    tt <- table(k$ref)
    k[k==""]<-"Unknown"
    if(opt=="isotype" | opt=="cluster"){
      k$seurat_clusters<-as.character(k$seurat_clusters)
      k$seurat_clusters[is.na(k$seurat_clusters)] <- "Unknown"
      final_df<-data.frame(Seq=paste(k$VDJ_sequence_nt_trimmed,k$VJ_sequence_nt_trimmed,sep=''),
                           Name=paste(k$clonotype_id_10x,k$VDJ_cgene,k$orig_barcode.x,"cluster",k$seurat_clusters,sep='_'))

      final_df$Name<-gsub("(.*)\\_(.*)", "\\1\\2", final_df$Name)
      final_df<-rbind(final_df,c(names(which.max(tt)),"germline_clusterUnknown"))
    }else{
      final_df<-data.frame(Seq=paste(k$VDJ_sequence_nt_trimmed,k$VJ_sequence_nt_trimmed,sep=''),
                           Name=k[opt])
      final_df<-rbind(final_df,c(names(which.max(tt)),"germline"))
    }
    k<-final_df
  }
  if(opt=="isotype"){
    res_iso<- unname(which(sapply(k, function(x) any(grepl(paste(type_iso,collapse="|"), x)))))
    names(k)[res_iso]<-"match"
    k$isotype<-.FIND_MATCHES(k$match, type_iso, ignore.case = TRUE)
    if(!(k$match[which(k$isotype=="Germline")]==k$isotype[which(k$isotype=="Germline")])){
      k$isotype[which(k$isotype=="Germline")]<-k$match[which(k$isotype=="Germline")]
    }
    k$match<-NULL
    res_clus<- unname(which(sapply(k, function(x) any(grepl("cluster", x)))))

    if(length(res_clus)>0){
      if(!res_iso==res_clus){
        if (is.integer(res_clus) && length(res_clus) == 0){
          res_clus_1<- unname(which(sapply(k, function(x) any(grepl("(#)", x)))))
          if (!(is.integer(res_clus_1) && length(res_clus_1) == 0)){
            k[,res_clus_1]<-NULL
          }
        }else{
          k[,res_clus]<-NULL
        }
      }
    }
    names_clus<- unname(which(sapply(names(k), function(x) any(grepl("Cluster", x)))))
    if(length(names_clus)>0){
      k[names_clus]<-NULL
    }

  }else if (opt=="cluster"){
    res_clus<- unname(which(sapply(k, function(x) any(grepl("cluster", x)))))
    if (is.integer(res_clus) && length(res_clus) == 0){
      res_clus_1<- unname(which(sapply(k, function(x) any(grepl("(#)", x)))))
      row_germ<-unname(which(sapply(unname(unlist(k[names(k)[res_clus_1][1]])), function(r) which(any(grepl("germline",r))))==1))
      k[res_clus_1][[1]][-row_germ]<- gsub("^.*?#","#",k[res_clus_1][[1]][-row_germ])
      names(k)[res_clus_1]<-"cluster"
    }else{
      row_germ<-unname(which(sapply(unname(unlist(k[names(k)[res_clus][[1]]])), function(r) which(any(grepl("germline",r))))==1))
      k[res_clus][[1]][row_germ]<-"germline"
      k[res_clus][[1]][-row_germ]<-gsub("^.*?cluster","",k[res_clus][[1]][-row_germ])
      names(k)[res_clus]<-"cluster"
      k$cluster[!nzchar(k$cluster)]<-"Unknown"
    }
    res_iso<- unname(which(sapply(k, function(x) any(grepl(paste(c("Germline","IGHA","IGHG1","IGHG2B","IGHG2C","IGHM","IGHG3","IGHD","IGHE"),collapse="|"), x)))))
    if (!(is.integer(res_iso) && length(res_iso) == 0)){
      k[res_iso]<-NULL
    }
  }else{
    pos_col<- unname(which(sapply(names(k), function(x) any(grepl(opt, x)))))
    if(!any(grepl("germline", k[pos_col]))){
      one_string<-apply(k,1,paste,collapse=" ")
      row_germ<-grep("germline", one_string)
      k[pos_col][[1]][row_germ]<-"germline"
    }

  }
  if(opt=="isotype"){
    pos_col<-which(names(k)=="isotype" )
    names(pos_col)<-"isotype"
    row_germ<-unname(which(sapply(unname(unlist(k[names(k)[pos_col][[1]]])), function(r) which(any(grepl("germline",r))))==1))
    k[pos_col][[1]][row_germ]<-"germline"
    if(length(is.na(k[pos_col][[1]]))>0){
      pos_NA<-which(is.na(k[pos_col][[1]]))
    }
    k[pos_col][[1]][pos_NA]<-"Unknown"
  }else if (opt=="cluster"){
    pos_col<-which(names(k)=="cluster" )
    names(pos_col)<-"cluster"
  }
  k<-lapply(k, function(z) gsub("[[:punct:]]","",z))
  id<-which(sapply(k, function(z) any(grepl("^[A-Za-z]+$", z))))
  if(!(opt=="isotype") & !(opt=="cluster") & (length(names(k))>2)){
    id<-id[-which(grepl("isotype",names(id)))]
    k<-k[id]

  }else if ( (opt=="isotype") |(opt=="cluster")) {
    k<-k[unname(unique(c(id,pos_col)))]
  }

  if(any(grep("cdr3", names(k)))){
    if(cdr3>0){
      col_pos<-which(sapply(names(k), function(z) any(grepl("[H].*[L]", z))))
      cols<-c("HC_Seq","LC_Seq")
      names(k)[grep("cdr3", names(k))] <- cols
      k[col_pos]<-NULL
      k$Seq<-paste(k$HC_Seq, k$LC_Seq,sep="")
      k[cols]<-NULL
    }else{
      k<-k[-grep("cdr3", names(k))]
      if(opt=="isotype"){
        res_seq<-which(sapply(k, function(z) any(grepl("[0-9][A-Za-z]", z))))
        if(!(names(-res_seq)=="Seq")){
          names(k)[-res_seq] <- "Seq"
        }
      }else if(opt=="cluster"){
        res_seq<-which(sapply(k, function(z) any(grepl("^[A-Za-z]+$", z))))
        if(!(names(res_seq)=="Seq")){
          names(k)[res_seq] <- "Seq"
        }
      }else{
        res_seq<-which(sapply(k, function(z) any(grepl("[A-Za-z][0-9]", z))))
        if(!(names(-res_seq)=="Seq")){
          names(k)[-res_seq] <- "Seq"
        }
      }
    }
    if(opt=="cluster"){
      names_clus<- unname(which(sapply(names(k), function(x) any(grepl("cluster", x)))))
      k[names_clus]<- lapply(k[names_clus], function(x) paste0("#",x))
      row_germ<-unname(which(sapply(unname(unlist(k[names(k)[names_clus][[1]]])), function(r) which(any(grepl("germline",r))))==1))
      k[names_clus][[1]][row_germ]<-substring(k[names_clus][[1]][row_germ], 2)
      k[names_clus][[1]][row_germ]<-sub("(germline)", "\\1#\\2",k[names_clus][[1]][row_germ])
    }
    temp<-c()
    if(cdr3>0){
      if(opt=="isotype"){
        clus_num<- unname(which(sapply(names(k), function(x) any(grepl("isotype", x)))))
        if(clus_num==1){
          temp<-append(temp,k[clus_num+1])
          temp<-append(temp,k[clus_num])
          names(temp)<-c("Seq","isotype")
        }
      }else if (opt=="cluster"){
        clus_num<- unname(which(sapply(names(k), function(x) any(grepl("cluster", x)))))
        if(clus_num==1){
          temp<-append(temp,k[clus_num+1])
          temp<-append(temp,k[clus_num])
          names(temp)<-c("Seq","cluster")
        }

      }else{
        clus_num<- unname(which(sapply(names(k), function(x) any(grepl(opt, x)))))
        if(clus_num==1){
          temp<-append(temp,k[clus_num+1])
          temp<-append(temp,k[clus_num])
          names(temp)<-c("Seq",opt)
        }
      }
      k<-temp
    }else{
      if(!(opt=="cluster") & !(opt=="isotype")){
        clus_num<- unname(which(sapply(names(k), function(x) any(grepl(opt, x)))))
        if(clus_num==1){
          temp<-append(temp,k[clus_num+1])
          temp<-append(temp,k[clus_num])
          names(temp)<-c("Seq",opt)
        }
        k<-temp

      }

    }
    pos_none<-which(is.na(k[[2]]))
    if(length(pos_none)>0){
      k[[2]][pos_none]<-"Unknown"
    }
  }
  if(!(names(k)[1]=="Seq")){
    names(k)[1] <- "Seq"
  }

  return(k)
}

.COUNT_ELEMENTS<-function(k,opt,distance_mat,tie_flag,weight,random.seed,alg_opt){
  if(opt=="isotype"){
    pre_data<-.PREPROCESSING(k$Seq,k$isotype,distance_mat,NULL)
  }else{
    pre_data<-.PREPROCESSING(k$Seq,k$cluster,distance_mat,NULL)
  }
  if(alg_opt=="two-step"){
    edgeList<-.FIND_EDGES_v1(pre_data$dist_mat,pre_data$countData)
    edgeList<-.DELETE_MULTIPLE_PARENTS(edgeList$node_list,edgeList$edges_final,pre_data$dist_mat,pre_data$countData,tie_flag,weight,random.seed)
  }else{
    germline.index<-which(apply(pre_data$countData, 1, function(r) any(grepl("Germline",r))))
    adj_mat<-.ADJ_MAT(pre_data$dist_mat)
    edgeList<-.EDGES(pre_data$dist_mat,adj_mat,germline.index,weight,random.seed)
  }
  columns <- grep('Ratio', colnames(pre_data$countData), value=TRUE)
  values <- lapply(edgeList$node_list, function(x) unlist(pre_data$countData[x,columns]))
  df_values <- do.call("rbind", values)
  df_values<-as.data.frame(df_values)
  return(nrow(df_values))
}

.SPLITCOMBINATION<-function(k,opt,distance_mat,tie_flag,weight,random.seed,alg_opt){
  list_IGHG<-list()
  list_IGHA<-list()
  list_IGHM<-list()
  list_IGAG<-list()
  list_other<-list()
  list<-list()

  if(opt=="isotype"){
    pre_data<-.PREPROCESSING(k$Seq,k$isotype,distance_mat,NULL)
  }else{
    pre_data<-.PREPROCESSING(k$Seq,k$cluster,distance_mat,NULL)
  }
  if(alg_opt=="two-step"){
    edgeList<-.FIND_EDGES_v1(pre_data$dist_mat,pre_data$countData)
    edgeList<-.DELETE_MULTIPLE_PARENTS(edgeList$node_list,edgeList$edges_final,pre_data$dist_mat,pre_data$countData,tie_flag,weight,random.seed)
  }else{
    germline.index<-which(apply(pre_data$countData, 1, function(r) any(grepl("Germline",r))))
    adj_mat<-.ADJ_MAT(pre_data$dist_mat)
    edgeList<-.EDGES(pre_data$dist_mat,adj_mat,germline.index,weight,random.seed)
  }
  columns <- grep('Ratio', colnames(pre_data$countData), value=TRUE)
  values <- lapply(edgeList$node_list, function(x) unlist(pre_data$countData[x,columns]))
  df_values <- do.call("rbind", values)
  df_values<-as.data.frame(df_values)
  count<-unlist(apply(df_values,1,function(x) which(x>0)))
  df_counts<-table(count)
  names(df_counts)<-c(columns)
  pos_name<-unname(which(sapply(names(df_counts), function(z) any(grepl("germline.*", z)))))#names(df_counts)[which(grepl("germline",names(df_counts)))]
  df_counts<-df_counts[ -pos_name]
  if(any(grepl("[I][G][H][G].*", names(df_counts)))){
    pos_g<-unname(which(sapply(names(df_counts), function(z) any(grepl("[I][G][H][G].*", z)))))
    df_counts <-c(df_counts, sum(df_counts[pos_g]))
    names(df_counts)[length(df_counts)]<-'IGHG'
    df_counts<-df_counts[-pos_g]
  }else{
    pos_name<-unname(which(sapply(columns, function(z) any(grepl("germline.*", z)))))#names(df_counts)[which(grepl("germline",names(df_counts)))]
    columns<-columns[ -pos_name]
    names(df_counts)<-columns
  }
  max_pos<-which(df_counts == max(df_counts[df_counts>0]))
  if(any(grepl("IGHG", names(max_pos))) & any(grepl("IGHM", names(max_pos)))){
    list_IGHG<-append(list_IGHG,k)
  }else if(any(grepl("IGHM", names(max_pos))) & any(grepl("IGHA", names(max_pos)))){
    list_IGHA<-append(list_IGHA,k)
  }else if(any(grepl("IGHA", names(max_pos))) & any(grepl("IGHG", names(max_pos)))){
    list_IGAG<-append(list_IGAG,k)
  }else if (any(grepl("IGHA", names(max_pos)))){
    list_IGHA<-append(list_IGHA,k)
  }else if(any(grepl("IGHM", names(max_pos)))){
    list_IGHM<-append(list_IGHM,k)
  }else if(any(grepl("IGHG", names(max_pos)))){
    list_IGHG<-append(list_IGHG,k)
  }else{
    list_other<-append(list_other,k)
  }

  list <- list(list_IGHG=list_IGHG,list_IGHA=list_IGHA,list_IGHM=list_IGHM,list_IGAG=list_IGAG,list_other=list_other)

  return(list)
}

.SPLITLINEAGES<-function(k){
  list_IGHG<-list()
  list_IGHA<-list()
  list_IGHM<-list()
  list_IGAG<-list()
  list_other<-list()
  list<-list()

  res_iso<- unname(which(sapply(k, function(x) any(grepl(paste(.REFERENCE_INFO()$type_iso,collapse="|"), x)))))
  count_igg<-sum(stringr::str_count(unlist(k[res_iso]), "IGHG"))
  count_iga<-sum(stringr::str_count(unlist(k[res_iso]), "IGHA"))
  count_igm<-sum(stringr::str_count(unlist(k[res_iso]), "IGHM"))
  if((count_igg>count_iga) &  (count_igm>count_iga) & (count_igg==count_igm)){
    list_IGHG<-append(list_IGHG,k)
  }else if((count_iga>count_igg) &  (count_igm>count_igg) & (count_igm==count_iga)){
    list_IGHA<-append(list_IGHA,k)
  }else if((count_iga>count_igm) &  (count_igg>count_igm) & (count_iga==count_igg)){
    list_IGAG<-append(list_IGAG,k)
  }else if ((count_iga>count_igg) & (count_iga>count_igm)){
    list_IGHA<-append(list_IGHA,k)
  }else if ((count_igm>count_igg) & (count_igm>count_iga)){
    list_IGHM<-append(list_IGHM,k)
  }else if ((count_igg>count_iga) & (count_igg>count_igm)){
    list_IGHG<-append(list_IGHG,k)
  }else{
    list_other<-append(list_other,k)
  }

  list <- list(list_IGHG=list_IGHG,list_IGHA=list_IGHA,list_IGHM=list_IGHM,list_IGAG=list_IGAG,list_other=list_other)

  return(list)

}


.HAS_EMPTY_LIST <- function(x) {
  if(is.list(x)) {
    if (length(x)==0) {
      return(TRUE)
    }else {
      return(any(vapply(x, .HAS_EMPTY_LIST, logical(1))))
    }
  } else {
    return(FALSE)
  }
}
