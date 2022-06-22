#' Function to get shared/public elements across multiple repertoires
#'
#'@description Function to get shared elements across multiple repertoires, specified by the feature.columns parameter (a column of the VDJ matrix). If two columns are specified in feature.columns, the resulting shared features will combine the values from each column (at a per-cell level).
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param feature.columns Character or character vector columns of features to be assayed
#' @param grouping.column string - the repertoire/group-defining column (default to 'sample_id').
#' @param specific.groups vector of strings or NULL - if only the shared elements from specific repertoires should be taken into account. If NULL, will output the shared/public elements across all repertoires.
#' @param find.public.all boolean - if T, will look for the public elements across all repertoires
#' @param find.public.percentage list - the first element denotes the percentage of repertoires to get shared elements for, the second element is the maximum number of repertoire combinations to consider (can be NULL to consider all).
#' @param treat.combined.features string - 'exclude' will exclude combined features with one element missing, 'include' will include and considers them as a new feature value.
#' @param output.format string - 'df' to get a shared element dataframe (with columns = Repertoire and Public), 'list' for a list of shared elements.
#' @return Either a dataframe of public elements across multiple repertoires or a list.
#' @export


#' @examples
#' VDJ_get_public(VDJ = small_vgm[[1]],
#' feature.columns='VDJ_cdr3s_aa', find.public.all=TRUE,
#' output.format='df')
#'



#ADD VENN DIAGRAMS TO PUBLICS


VDJ_public <- function(VDJ,
                       feature.columns,
                       grouping.column,
                       specific.groups,
                       find.public.all,
                       find.public.percentage,
                       treat.combined.features,
                       output.format){

  if(missing(VDJ)) stop('Please input your data as a VDJ matrix')
  if(missing(feature.columns)) feature.columns <- 'CDR3aa'
  if(missing(grouping.column)) grouping.column <- 'sample_id'
  if(missing(specific.groups)) specific.groups <- NULL
  if(missing(find.public.all)) find.public.all <- T
  if(missing(find.public.percentage)) find.public.percentage <- list(0.6, NULL)
  if(missing(treat.combined.features) & length(feature.columns==2)) treat.combined.features <- 'exclude'
  if(missing(output.format)) output.format <- 'plot'

  get_feature_combinations <- function(x, y, split.x, split.y, split.by=';', collapse.by=';', combine.sequences=F){
   if(split.x==T) x <- stringr::str_split(x, split.by ,simplify=T)[1,]
   if(split.y==T) y <- stringr::str_split(y, split.by ,simplify=T)[1,]

   ccombs <- expand.grid(x,y)
   if(!combine.sequences){
     ccombs<-paste0(ccombs[,1], ' ', ccombs[,2])
   }else{
     ccombs<-paste0(ccombs[,1], ccombs[,2])

   }
   ccombs <- paste0(ccombs, collapse=collapse.by)

   return(ccombs)
  }

  if(('CDR3aa' %in% feature.columns) & !('CDR3aa' %in% colnames(VDJ))){
    VDJ$CDR3aa <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T, combine.sequences=T)} else '', VDJ$VDJ_cdr3s_aa, VDJ$VJ_cdr3s_aa)
  }

  for(i in 1:length(feature.columns)){
   if(!feature.columns[i] %in% names(VDJ)){
     stop("Please provide valid feature column name(s) contained within VDJ")
   }
  }

  if(length(feature.columns)==2){
   if(treat.combined.features=='exclude'){
     combined_features <- mapply(function(x,y) if(!is.null(x) & !is.null(y) & !is.na(x) & !is.na(y) & x!='' & y!='') {get_feature_combinations(x,y,split.x=T,split.y=T)} else '', VDJ[,feature.columns[[1]]], VDJ[,feature.columns[[2]]])
   }else{
     combined_features <- mapply(function(x,y) get_feature_combinations(x,y,split.x=T,split.y=T), VDJ[,feature.columns[[1]]], VDJ[,feature.columns[[2]]])
     VDJ$new_feature <- combined_features
     new_feature <- paste0(feature.columns[[1]], '/', feature.columns[[2]])
     names(sample_dfs[[i]])[names(VDJ)=='new_feature'] <- new_feature
   }
   feature.columns <- paste0(feature.columns[[1]], '/', feature.columns[[2]])
  }

  if(is.null(specific.groups)){
    repertoire_numbers <- unique(VDJ[,grouping.column])
  }else{
    repertoire_numbers <- unique(specific.groups)
  }


  if(find.public.all){

    unique_per_rep <- list()
    for(i in 1:length(repertoire_numbers)) {
      VDJ_subset <- VDJ[which(VDJ[grouping.column]==repertoire_numbers[i]),]
      rep_elements <- VDJ_subset[,feature.columns]
      rep_elements <- lapply(rep_elements, function(x) stringr::str_split(x, ';'))
      rep_elements <- unlist(unique(rep_elements))
      rep_elements <- rep_elements[rep_elements!='']
      unique_per_rep[[i]] <- rep_elements
    }
    public_elements <- Reduce(intersect, unique_per_rep)

    if(output.format=='df'){
      reps <- rep(list(repertoire_numbers), length(public_elements))
      public <- unlist(public_elements)
      return(data.frame(Repertoires=matrix(reps), Public=public))
    }else if(output.format == 'list'){
      return(list(public_elements, repertoire_numbers))
    }else{

    }


  }else if(!is.null(find.public.percentage[[1]])){
    max_rep_number <- as.integer(find.public.percentage[[1]] * length(repertoire_numbers))
    combs <- utils::combn(repertoire_numbers, max_rep_number)

    if(!is.null(find.public.percentage[[2]])){
      if(find.public.percentage[[2]] < ncol(combs)) { combs <- combs[,1:find.public.percentage[[2]]] }
    }

    unique_per_comb <- list()
    reps_per_comb <- list()

    for(i in 1:ncol(combs)){
      reps <- combs[,i]

      unique_per_rep <- list()

      for(j in 1:length(reps)){
        VDJ_subset <- VDJ[which(VDJ[grouping.column]==repertoire_numbers[j]),]
        rep_elements <- VDJ_subset[,feature.columns]
        rep_elements <- lapply(rep_elements, function(x) stringr::str_split(x, ';'))
        rep_elements <- unlist(unique(rep_elements))
        rep_elements <- rep_elements[rep_elements!='']
        unique_per_rep[[j]] <- rep_elements
      }

      unique_per_comb[[i]] <- Reduce(intersect, unique_per_rep)
      reps_per_comb[[i]] <- reps
    }

    if(output.format=='df'){
    #  df_list <- list()
    #  for(i in 1:length(unique_per_comb)){
    #  df_list[[i]] <- data.frame(Repertoires=rep(reps_per_comb[[i]], length(unique_per_comb[[i]])), Public=unique_per_comb[[i]])
    #  }
    #  df_final <- do.call('rbind', df_list)
    #  unique_public <- unique(df_final$Public)
    #  unique_per_repertoire <- lapply(unique(df_final$Public), function(x) unique(df_final$Repertoires[which(df_final$Public==x)]))
    #  unique_per_repertoire <- lapply(unique_per_repertoire, function(x) x[order(nchar(x), x)])
    #  return(data.frame(Repertoires=matrix(unique_per_repertoire), Public=unique_public))
      return(data.frame(Repertoire=matrix(reps_per_comb), Public=matrix(unique_per_comb)))

    }else if(output.format == 'list'){
      return(list(matrix(unique_per_comb), reps_per_comb))
    }else{

    }
  }else{
    stop('Method not implemented yet!')
  }
}





mtcars <- datasets::mtcars


set1 <- rownames(subset(mtcars, mpg>18))
set2 <- rownames(subset(mtcars, qsec>18))
set3 <- rownames(subset(mtcars, cyl<5))
set4 <- rownames(subset(mtcars, wt>2.5))
set5 <- rownames(subset(mtcars, gear==4))
carset <- list(mpg=set1, qsec=set2, cyl=set3, wt=set4, gear=set5)

ggplot() +
  geom_venn(carset) +
  theme_void()

oo <- VDJ_public(vdj)
oo

vdj <- read.csv('./data/OVA_vgm_all_annotations.csv')
