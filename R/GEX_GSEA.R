#' Conducts a Gene Set Enrichment Analysis (GSEA) on a set of genes submitted in a data frame with a metric each.
#' Works with the output of GEX_genes_cluster or a custom data frame containing the gene symbols either in a column "symbols" or as rownames and a metric for each gene.
#' The name of the column containing the metric has to be declared via the input metric.colname.
#' @param GEX.cluster.genes.output Data frame containing the list of gene symbols and a metric. Function works directly with GEX_cluster_genes output.
#' @param MT.Rb.filter Logical, should Mitotic and Ribosomal genes be filtered out of the geneset. True by default.
#' @param filter Character vector containing the identifying symbol sequence for the genes which should be filtered out, if MT.Rb.filter == T. By default set to c("MT-", "RPL", "RPS").
#' @param path.to.pathways Either a path to gmt file containing the gene sets (can be downloaded from MSigDB) or vector where first element specifies species and second element specifies the MSigDB collection abbreviation. E.g.: c("Homo sapiens", "H"). Mouse C7 (immunologic signature) gene set will be used by default.
#' @param metric.colname Name of column which contains the metric used for the ranking of the submitted genelist. "avg_logFC" is used by default.
#' @param pval.adj.cutoff Only genes with a more significant adjusted pvalue are considered. Default: 0.001
#' @param Enrichment.Plots List of Gene-set names which should be plotted as Enrichment plots in addition to the top 10 Up and Downregulated Genesets.
#' @param my.own.geneset A list, where each element contains a gene list and is named with the corresponding pathway name. Default is set to FALSE, so that gene sets from MSigDB are used. Should not contain ".gmt" in name.
#' @param eps Numeric, specifying boundary for calculating the p value in the GSEA.
#' @param platypus.version Function works with V2 and V3, no need to set this parameter.
#' @param verbose Print run parameters and status to console
#' @return Returns a list containing a tibble with the gene sets and their enrichment scores and Enrichment plots. List element [[1]]: Dataframe with Genesets and statistics. [[2]]: Enrichment plots of top10 Up regulated genesets. [[3]]: Enrichment plots of top10 Down regulated genesets. [[4]]: Enrichment plots of submited gene-sets in parameter Enrichment.Plot.
#' @export
#' @examples
#' \dontrun{
#' df <- GEX_cluster_genes(gex_combined[[1]])
#'
#' #Using gmt file to perform gsea
#' output <- GEX_GSEA(GEX.cluster.genes.output =  df[[1]], MT.Rb.filter = TRUE
#' , path.to.pathways = "./c5.go.bp.v7.2.symbols.gmt")
#' cowplot::plot_grid(plotlist=output[[2]], ncol=2)
#' View(gex_gsea[[1]])
#'
#' #Directly downloading gene set collection from MSigDB to perform gsea
#' output <- GEX_GSEA(GEX.cluster.genes.output =  df[[1]], MT.Rb.filter = TRUE
#' , path.to.pathways = c("Mus musculus", "C7"))
#'
#' #Using your own gene list to perform gsea
#' output <- GEX_GSEA(GEX.cluster.genes.output =  df[[1]], MT.Rb.filter = TRUE
#' , my.own.geneset = my_geneset)
#'}

GEX_GSEA <- function(GEX.cluster.genes.output, MT.Rb.filter, filter, path.to.pathways, metric.colname, pval.adj.cutoff, Enrichment.Plots, my.own.geneset, eps, platypus.version, verbose){

  p_val_adj <- NULL
  stats <- NULL
  symbol <- NULL
  NES <- NULL
  ES <- NULL
  pval <- NULL
  pathway <- NULL

  ###
  if(missing(verbose)) verbose <- F
  if (missing(filter)) {
    filter <- c("MT-", "RPL", "RPS")
    if(verbose) message("filter parameter set to c('MT-', 'RPL', 'RPS'). Mitochondrial and ribosomal genes will be filtered out")}
  if (missing(metric.colname)) {
    metric.colname <- "avg_logFC"
    if(verbose) message("metric.colname parameter set to avg_logFC")}
  if (missing(pval.adj.cutoff)) {
    pval.adj.cutoff <- 0.001
    if(verbose) message("pval.adj.cutoff parameter set to 0.001")}
  if (missing(my.own.geneset)) {my.own.geneset <- F}
  if (missing(path.to.pathways)) {
    path.to.pathways <- c("Mus musculus", "C7")}
  if (class(my.own.geneset) == "logical") {if(verbose) message(paste0("MSigDB collection: ", path.to.pathways))}
  if (class(my.own.geneset) == "list") {if(verbose) message("Own gene set is being used")}
  if (missing(eps)) {
    eps <- 1e-10
    if(verbose) message("eps parameter set to 1e-10")}

  platypus.version <- "Does not matter"

  # change metric colname to 'stats' for further downstream analysis
  colnames <- names(GEX.cluster.genes.output)[names(GEX.cluster.genes.output) == metric.colname] <- 'stats'
  # if symbols are given as rownames, put them in own column
  if(!is.character(GEX.cluster.genes.output$symbol)){
    GEX.cluster.genes.output$symbol <- rownames(GEX.cluster.genes.output)
  }
    # Filter out RBS, RBL and MT genes if desired
    if (MT.Rb.filter==T){
      exclude <- c()
      for (j in filter) {
        exclude <- c(exclude, stringr::str_which(GEX.cluster.genes.output$symbol, j))
      }
      if(length(exclude) > 0){
        print(paste0("Filtering ", length(exclude), "genes"))
        df <- GEX.cluster.genes.output[-exclude,]
      }else{
        warning("No genes matching filter found")
        df <- GEX.cluster.genes.output
      }
    }else{
      df <- GEX.cluster.genes.output
    }

    # create ranked list
    df %>% dplyr::filter(p_val_adj<pval.adj.cutoff)%>% dplyr::select("symbol","stats")%>% stats::na.omit()%>%dplyr::arrange(-stats)%>% dplyr::distinct(symbol, .keep_all = TRUE)-> df_ranked
    df_ranked <- tibble::deframe(df_ranked)

    if (class(my.own.geneset) == "logical"){
      if(sum(grepl(".gmt$", path.to.pathways)) == length(path.to.pathways)){
        pathway_MSig <- fgsea::gmtPathways(path.to.pathways)
      } else{
        #PACKAGE is not loaded! Please revise!
        pathway_MSig <- msigdbr::msigdbr(species = path.to.pathways[[1]], category=path.to.pathways[[2]])
        pathway_MSig <-  split(x = toupper(pathway_MSig$gene_symbol), f = pathway_MSig$gs_name)
      }
    } else{
      pathway_MSig <- my.own.geneset
    }

    #Run GSEA %>% safe as df
    fgsea_res <- fgsea::fgseaMultilevel(pathways=pathway_MSig, stats=df_ranked, minSize=2, maxSize=500, eps = eps)
    print(fgsea_res)
    fgsea_res_Tidy <- fgsea_res %>%
      tidyr::as_tibble() %>%
      dplyr::arrange(IRanges::desc(NES))
    print(fgsea_res_Tidy)
    topPathwaysUp <- fgsea_res[ES > 0][utils::head(order(pval), n=10), pathway]
    topPathwaysDown <- fgsea_res[ES < 0][utils::head(order(pval), n=10), pathway]
    topPathways<- c(topPathwaysUp, rev(topPathwaysDown))

    plotsUp <- list()
    if(length(topPathwaysUp)>0){
      for (k in 1:length(topPathwaysUp)){
        plotsUp[[k]]<- fgsea::plotEnrichment(pathway_MSig[[topPathwaysUp[[k]]]],
                       df_ranked) + ggplot2::labs(title=topPathwaysUp[[k]])
      }
    }
    plotsDown <- list()
    if(length(topPathwaysDown)>0){
      for (k in 1:length(topPathwaysDown)){
        plotsDown[[k]]<- fgsea::plotEnrichment(pathway_MSig[[topPathwaysDown[[k]]]],
                                      df_ranked) + ggplot2::labs(title=topPathwaysDown[[k]])
      }
    }

    plotsCustom <- list()
    if (!missing(Enrichment.Plots)) {
    for (k in 1:length(Enrichment.Plots)){
      plotsCustom[[k]]<- fgsea::plotEnrichment(pathway_MSig[[Enrichment.Plots[[k]]]],
                                      df_ranked) + ggplot2::labs(title=Enrichment.Plots[[k]])
    }
    }
    output <- list()
    output[[1]] <- fgsea_res_Tidy
    output[[2]] <- plotsUp
    output[[3]] <- plotsDown
    output[[4]] <- plotsCustom
    return(output)
}
