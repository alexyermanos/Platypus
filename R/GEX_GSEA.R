#' Conducts a Gene Set Enrichment Analysis (GSEA) on a set of genes submitted in a data frame with a metric each. Works with the output of GEX_genes_cluster or a custom data frame contianing the gene symbols either in a column "symbols" or as rownames and a metric for each gene. The name of the column containing the metric has to be declared via the input metric_colname.
#' @param GEX.cluster.genes.output Data frame containing the list of gene symbols and a metric. Function works directly with GEX_cluster_genes output.
#' @param MT.Rb.filter Logical, should Mitotic and Ribosomal genes be filtered out of the geneset. True by default.
#' @param filter Character vector containing the identifying symbol sequence for the genes which should be filtered out, if MT.Rb.filter == T. By default set to c("MT-", "RPL", "RPS").
#' @param path_to_pathways Path to gmt file containing the gene sets. Can be downloaded from MSigDB.
#' @param metric_colname Name of column which contains the metric used for the ranking of the submitted genelist. "avg_logFC" is used by default.
#' @param pval_adj_cutoff Only genes with a more significant adjusted pvalue are considered. Default: 0.001
#' @param Enrichment.Plots List of Gene-set names which should be plotted as Enrichment plots in addition to the top 10 Up and Downregulated Genesets.
#' @return Returns a list containing a tibble with the gene sets and their enrichment scores and Enrichment plots. List element [[1]]: Dataframe with Genesets and statistics. [[2]]: Enrichment plots of top10 Up regulated genesets. [[3]]: Enrichment plots of top10 Down regulated genesets. [[4]]: Enrichment plots of submited gene-sets in parameter Enrichment.Plot.
#' @export
#' @examples
#' \dontrun{
#' df <- GEX_cluster_genes(gex_combined[[1]])
#' output <- GEX_GSEA(GEX.cluster.genes.output =  df[[1]], MT.Rb.filter = T, path_to_pathways = "./c5.go.bp.v7.2.symbols.gmt")
#' cowplot::plot_grid(plotlist=output[[2]], ncol=2)
#' View(gex_gsea[[1]])
#'}

GEX_GSEA <- function(GEX.cluster.genes.output, MT.Rb.filter, filter, path_to_pathways, metric_colname, pval_adj_cutoff, Enrichment.Plots){
  if (missing(filter)) {filter <- c("MT-", "RPL", "RPS")}
  if (missing(metric_colname)) {metric_colname <- "avg_logFC"}
  if (missing(pval_adj_cutoff)) {pval_adj_cutoff <- 0.001}
  
  require(dplyr)
  require(fgsea)
  require(tibble)
  require(stringr)
  require(stats)
  require(ggplot2)

  # change metric colname to 'stats' for further downstream analysis
  # print(head(GEX.cluster.genes.output))
  colnames <- names(GEX.cluster.genes.output)[names(GEX.cluster.genes.output) == metric_colname] <- 'stats'
  # if symbols are given as rownames, put them in own column
  if(!is.character(GEX.cluster.genes.output$symbol)){
    GEX.cluster.genes.output$symbol <- rownames(GEX.cluster.genes.output)
  }
  # print(head(GEX.cluster.genes.output))
    # Filter out RBS, RBL and MT genes if desired
    if (MT.Rb.filter==T){
      exclude <- c()
      for (j in filter) {
        exclude <- c(exclude, str_which(GEX.cluster.genes.output$symbol, j))
      }
      if(length(exclude) > 0){
        print(paste0("Filtering ", length(exclude), "genes"))
        df <- GEX.cluster.genes.output[-exclude,]
      }else{
        print("No genes matching filter found")
        df <- GEX.cluster.genes.output
      }
    }else{
      df <- GEX.cluster.genes.output
    }
  # fgsea_res_Tidy %>%
    #   dplyr::select(-leadingEdge, -ES) %>%
    #   arrange(padj) %>%
    #   DT::datatable()-> data.table

    # create ranked list
    df %>% dplyr::filter(., p_val_adj<pval_adj_cutoff)%>% dplyr::select("symbol","stats")%>% na.omit()%>%dplyr::arrange(-stats)%>% distinct(symbol, .keep_all = TRUE)-> df_ranked
    df_ranked <- deframe(df_ranked)
    pathway_MSig <- gmtPathways(path_to_pathways)
    #Run GSEA %>% safe as df
    print("pre-gsea")
    fgsea_res <- fgseaMultilevel(pathways=pathway_MSig, stats=df_ranked, minSize=2, maxSize=500)
    print(fgsea_res)
    print("post-gsea")
    fgsea_res_Tidy <- fgsea_res %>%
      as_tibble() %>%
      arrange(desc(NES))
    print(fgsea_res_Tidy)
    topPathwaysUp <- fgsea_res[ES > 0][head(order(pval), n=10), pathway]
    topPathwaysDown <- fgsea_res[ES < 0][head(order(pval), n=10), pathway]
    topPathways<- c(topPathwaysUp, rev(topPathwaysDown))

    plotsUp <- list()
    if(length(topPathwaysUp)>0){
      for (k in 1:length(topPathwaysUp)){
        plotsUp[[k]]<- plotEnrichment(pathway_MSig[[topPathwaysUp[[k]]]],
                       df_ranked) + labs(title=topPathwaysUp[[k]])
      }
    }
    plotsDown <- list()
    if(length(topPathwaysDown)>0){
      for (k in 1:length(topPathwaysDown)){
        plotsDown[[k]]<- plotEnrichment(pathway_MSig[[topPathwaysDown[[k]]]],
                                      df_ranked) + labs(title=topPathwaysDown[[k]])
      }
    }

    plotsCustom <- list()
    if (!missing(Enrichment.Plots)) {
    for (k in 1:length(Enrichment.Plots)){
      plotsCustom[[k]]<- plotEnrichment(pathway_MSig[[Enrichment.Plots[[k]]]],
                                      df_ranked) + labs(title=Enrichment.Plots[[k]])
    }
    }
    output <- list()
    output[[1]] <- fgsea_res_Tidy
    output[[2]] <- plotsUp
    output[[3]] <- plotsDown
    output[[4]] <- plotsCustom
    return(output)



}
