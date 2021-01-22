#' Runs a GO term analysis on a submitted list of genes. Works with the output of GEX_topN_DE_genes_per_cluster or a custom list of genes to obtain GOterms.
#' @param GEX.cluster.genes.output Either output of Platypus::GEX_cluster_genes or custom character vector containing gene symbols. A custom gene list will not be further filtered or ordered.
#' @param topNgenes How many of the most significant up or down regulated genes should be considered for GO term analysis. All genes will be used if left empty.
#' @param species The species the genes belong to. Default: "Mm".
#' @param ontology Ontology used for the GO terms. "MF", "BP" or "CC" possible. Default: "BP"
#' @param up.or.down Whether up or downregulated genes should be used for GO term analysis if GEX_cluster_genes output is used. Default: "up"
#' @param MT.Rb.filter logical, if mitochondrial and ribosomal genes should be filtered out.
#' @param kegg logical, if KEGG pathway analysis should be conducted. Requires internet connection. Default: False.
#' @return Returns a list of data frames containing the TopGO and the TopKEGG output containing the significant GO terms.
#' @export
#' @examples
#' \dontrun{
#' test <- GEX_GOterm(DE_genes_cluster,MT.Rb.filter = T, species= "Mm", ontology = "MF")
#' test2 <- GEX_GOterm(rownames(DE_genes_cluster[[1]]),MT.Rb.filter = T, species= "Mm", ontology = "BP")
#'}



GEX_GOterm <- function(GEX.cluster.genes.output, topNgenes, ontology, species,  up.or.down , MT.Rb.filter, kegg){

  require(org.Mm.eg.db)
  require(edgeR)

  if (missing(GEX.cluster.genes.output)) {print("Pleas submit either GEX_cluster_genes output or list of gene Symbols")}
  if (missing(topNgenes)) {}
  if (missing(ontology)) {ontology <- "BP"}
  if (missing(species)) {species <- "Mm"}
  if (missing(up.or.down)) {up.or.down <- "up"}
  if (missing(MT.Rb.filter)) {MT.Rb.filter <- T}
  if (missing(kegg)){kegg <- F}

  gene.list <- list()
  list_topGO <- list()
  list_topKEGG <- list()
  list <- list()

  #Check whether GEX_cluster_genes output was submitted or list as character vector
  if (class(GEX.cluster.genes.output[[1]])=="data.frame"){

    if(MT.Rb.filter==T){
      for (i in 1:length(GEX.cluster.genes.output)){
        #Filter Genes
        GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output[[i]]
        if (nrow(GEX.cluster.genes.output[[i]][grep('^RPS', rownames(GEX.cluster.genes.output[[i]])),])) #Only filter if there is at least one occurence, otherwise there is error
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output[[i]][-grep('^RPS', rownames(GEX.cluster.genes.output[[i]])),] #remove all genen startign with RPS
        }
        if (nrow(GEX.cluster.genes.output[[i]][grep('^RPL', rownames(GEX.cluster.genes.output[[i]])),]))
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output_filtered[-grep('^RPL', rownames(GEX.cluster.genes.output_filtered)),] #remove all genen startign with RPL
        }
        if (nrow(GEX.cluster.genes.output[[i]][grep('^MT', rownames(GEX.cluster.genes.output[[i]])),]))
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output_filtered[-grep('^MT', rownames(GEX.cluster.genes.output_filtered)),] #remove all genen startign with MT
        }
        GEX.cluster.genes.output[[i]] <- GEX.cluster.genes.output_filtered
      }
    }

      for (i in 1:length(GEX.cluster.genes.output)){
        #Add entrez-ID
        rownames(GEX.cluster.genes.output[[i]])<- str_to_title(rownames( GEX.cluster.genes.output[[i]])) #convert Symbol from all uppercase to mixed (eg. GZMK to Gzmk)
        GEX.cluster.genes.output[[i]]$symbol <- select(org.Mm.eg.db,keys=rownames( GEX.cluster.genes.output[[i]]), keytype="SYMBOL" ,columns=c("ENTREZID", "GENENAME"))  #rownames(DE_gene_list_filtered)

        #filter for positive logfoldchanges and arrange for increasing logfoldchanges
        if(up.or.down=="down"){
          GEX.cluster.genes.output[[i]] %>% filter(., GEX.cluster.genes.output[[i]]$avg_logFC < 0) %>% arrange(p_val_adj)-> GEX.cluster.genes.output[[i]]
        }else{
          GEX.cluster.genes.output[[i]] %>% filter(., GEX.cluster.genes.output[[i]]$avg_logFC > 0) %>% arrange(p_val_adj)-> GEX.cluster.genes.output[[i]]
        }

        if(missing(topNgenes)){
          print("missing topNgenes argument: all genes will be used")
          gene.list[[i]] <- GEX.cluster.genes.output[[i]]$symbol$ENTREZID
        }else{
          gene.list[[i]] <- head(GEX.cluster.genes.output[[i]]$symbol$ENTREZID, topNgenes)
        }
      }
  }else{

    GEX.cluster.genes.output <- data.frame(symbol=GEX.cluster.genes.output, dummy=NA) #add dummy, otherwise filter does not recognize df
    rownames(GEX.cluster.genes.output) <- str_to_upper((GEX.cluster.genes.output$symbol))

    if(MT.Rb.filter==T){
        #Filter Genes
        GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output
        if (nrow(GEX.cluster.genes.output[grep('^RPS', rownames(GEX.cluster.genes.output)),])) #Only filter if there is at least one occurence, otherwise there is error
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output[-grep('^RPS', rownames(GEX.cluster.genes.output)),] #remove all genen startign with RPS
        }
        if (nrow(GEX.cluster.genes.output[grep('^RPL', rownames(GEX.cluster.genes.output)),]))
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output_filtered[-grep('^RPL', rownames(GEX.cluster.genes.output_filtered)),] #remove all genen startign with RPL
        }
        if (nrow(GEX.cluster.genes.output[grep('^MT', rownames(GEX.cluster.genes.output)),]))
        {
          GEX.cluster.genes.output_filtered <- GEX.cluster.genes.output_filtered[-grep('^MT', rownames(GEX.cluster.genes.output_filtered)),] #remove all genen startign with MT
        }
        GEX.cluster.genes.output <- GEX.cluster.genes.output_filtered
    }

    rownames(GEX.cluster.genes.output)<- str_to_title((GEX.cluster.genes.output$symbol)) #convert Symbol from all uppercase to mixed (eg. GZMK to Gzmk)
    GEX.cluster.genes.output$symbol <- select(org.Mm.eg.db,keys=rownames( GEX.cluster.genes.output), keytype="SYMBOL" ,columns=c("ENTREZID", "GENENAME"))  #rownames(DE_gene_list_filtered)
    gene.list[[1]] <- GEX.cluster.genes.output$symbol$ENTREZID
  }

  for (i in 1:length(gene.list)){

    #Do GO-Term analysis
    go <- goana(gene.list[[i]], species = species)
    list_topGO[[i]] <- topGO(go, ontology = ontology, number = Inf)

    if (kegg == T){
      #if requested to allso KEGG-pathway analysis
      keg <- kegga(gene.list[[i]], species = species)
      list_topKEGG[[i]] <- topKEGG(keg)
    }
  }
  list[[1]]<-list_topGO

  if (kegg==T){
  list[[2]]<-list_topKEGG
  }

  return(list)

}


