GEX_phenotype <- function(seurat.object, cell.state.names, cell.state.markers, default){
  require(do)
  require(useful)
  Cap<-function(x){
    temp<-c()
    for (i in 1:length(x)){
      s <- strsplit(x, ";")[[i]]
      temp[i]<-paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)), sep="", collapse=";")
    } 
    return(temp)
  }
  if(missing(cell.state.markers)&default==T){
    cell.state.markers<-c("CD4+;CD44-",
                          "CD4+;IL7R+;CD44+",
                          "CD4+;CD44+;IL7R-;IFNG+",
                          "CD8A+;TCF7+;CD44-",
                          "CD8A+;CX3CR1+;IL7R-",
                          "CD8A+;IL7R+;CD44+",
                          "PDCD1+;CD8A+",
                          "CD19+;CD27-;CD38-",
                          "FAS+;CD19+",
                          "SDC1+",
                          "CD38+;FAS-")
  }
  if(missing(cell.state.names)&default==T){
    cell.state.names<-c("NaiveCd4",
                        "MemoryCd4",
                        "ActivatedCd4",
                        "NaiveCd8",
                        "EffectorCd8",
                        "MemoryCd8",
                        "ExhaustedCd8",
                        "NaiveBcell",
                        "GerminalcenterBcell",
                        "Plasmacell",
                        "MemoryBcell")
  }
  
  is.hum<-all(useful::find.case(rownames(seurat.object),case="upper"))
  if(is.hum==F){
    cell.state.markers<-Cap(cell.state.markers)
  }
  if(is.hum==T){
    cell.state.markers <- toupper(cell.state.markers) 
  }
  #parse cell state markers
  cell.state.markers<-Replace(cell.state.markers,from=";", to="&")
  cell.state.markers<-Replace(cell.state.markers,from="\\+", to=">0")
  cell.state.markers<-Replace(cell.state.markers,from="-", to="==0")
  #execute cmd
  seurat.object[["previous.ident"]] <- Idents(object = seurat.object)#(clusters ID)
  Idents(seurat.object)<-"Unclassified"
  cmd<-c()
  for(i in 1:length(cell.state.names)){
    
    cmd[i]<-paste0(cell.state.names[i],"<-WhichCells(seurat.object, slot = 'counts', expression =", cell.state.markers[i],")")
    is.exist<-tryCatch(expr=length(eval(parse(text=cmd[i]))), error = function(x){
      x<-F
      return(x)})
    if(is.exist!=F){
      Idents(object = seurat.object, cells = eval(parse(text=cell.state.names[i])) ) <- cell.state.names[i] 
    }
    
  }
  seurat.object[["cell.state"]] <- Idents(object = seurat.object)
  Idents(object = seurat.object)<-seurat.object[["previous.ident"]]
  
  return(seurat.object)
}