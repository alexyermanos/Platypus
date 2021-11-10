#' Plot igraph and ggplot objects

#' @description PlotGraphs takes as input the output of AntibodyForest or ForestMetrics functions and plots all corresponding networks in the single cell immune repertoire or the corresponding ggplot object, the user specifies, from all clone lineages.
#' The function gives the option in the user to store each tree or ggplot object within the repertoire in pdf format.
#' @param graphs a list of networks (Output of AntibodyForest function) or metrics (Output of ForestMetrics function).
#' @param no_arg  element of list the user desires to plot : integer value,if the user desires to plot a metric and NULL, if the user desires to plot the networks.
#' @param topdf logical value, TRUE if user wants to store plots in pdf format (the no_arg element of each list is saved in a separate page of the pdf).
#' @param filename name of saved pdf file based on the user's preferences.
#' @return Empty, output plots are written to file as specified by the user with the parameter filename
#' @export
#' @seealso AntibodyForest, ForestMetrics
#' @examples
#' \dontrun{
#' PlotGraphs(graphs,no_arg=NULL,topdf=TRUE,filename)
#' PlotGraphs(graphs,no_arg5,topdf=TRUE,filename)
#'}

AbForests_PlotGraphs<-function(graphs,no_arg,topdf,filename){

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
      .WRITE_TO_PDF(filename)
      lapply(temp, function(p){
        if(is.null(p$network)){
          p<-p[sapply(p,function(x) "network" %in% names(x))][[1]]
        }
        plot(p$network)
        if (!is.null(p$legend)){
          graphics::legend(p$legend[[1]],bty = p$legend[[2]],legend=p$legend[[3]],fill=p$legend[[4]],border=p$legend[[5]],inset=p$legend[[6]],xpd=p$legend[[7]],title=p$legend$title,cex=p$legend[[10]])

          }
      })

      grDevices::dev.off()
    }else{
      lapply(temp, function(p){
        if(is.null(p$network)){
          p<-p[sapply(p,function(x) "network" %in% names(x))][[1]]
        }
        plot(p$network)
        if (!is.null(p$legend)){
          graphics::legend(p$legend[[1]],bty = p$legend[[2]],legend=p$legend[[3]],fill=p$legend[[4]],border=p$legend[[5]],inset=p$legend[[6]],xpd=p$legend[[7]],title=p$legend$title,cex=p$legend[[10]])
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
          if(all(match( gg_class, unlist(lapply(print_res, class))))){
            .WRITE_TO_PDF(filename)
            invisible(lapply(print_res, print))
            grDevices::dev.off()
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
