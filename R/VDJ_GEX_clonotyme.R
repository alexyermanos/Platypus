#' Pseudotime analysis for scRNA and repertoire sequencing datasets
#'
#' @param method Pseudotime analylysis method to be used. Possible parameters are monocle3 or velocyto. monocle3 is being used as a default method. For velocyto analysis please run on Cluster and it is only available for UNIX based systems.
#' @param version Platypus version to use "v2" or "v3". version 2 used by default.
#' @param top.N.clonotypes How many clonotypes to show per sample in the Ridgeline plots and on the Velocyto UMAP.
#' @param vdj.gex.matrix.output If Platypus v3 is used, the input to this function has to be the output of the VDJ_GEX_matrix function.
#' @param vdj.analyze.output If Platypus v2 is used, the VDJ_analyze output has to be supplied.
#' @param gex.automate.output If Platypus v2 is used, the GEX_automate output has to be supplied here.
#' @param exclude.clusters Please enter a cluster number if you'd like to exclude a certain cluster from analysis. Cluster will be assigned to different partition in Monocle3 analysis and therefore pseudotime distance will be set to infinity. Cells from this cluster will be deleted from the dataset in the Velocyto analysis.
#' @param colors Vector containing custom colors to be used for highlighting the clonotypes. If left empty, default colors will be assigned.
#' @param show.cells Logical, should cells be shown in the Ridgeline plots. True by default.
#' @param highlight.genes Vector containing gene names. The expressionlevels of these genes along pseudotime will be plotted.
#' @param dropest.output.list List containing the cell.counts.matrices.rds from the Dropest alignment for Velocyto analysis.
#' @param velocyto.gex.merged Logical whether samples should be shown in combined UMAP or sepeartely.
#' @param velocyto.file.name String used as file name when saving the output pdf
#' @param velocyto.out.dir Directory to save the output files. By default the current working directory.
#' @param velocyto.save.rds If RDS objects should be saved as well. Default = F.
#' @param velocyto.norm.scale.factor Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.n.variable.features Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.neighbor.dim Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.cluster.resolution Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.mds.dim Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.nCount_spliced Cutoff threshold. cells with less spliced gene counts will be omitted. Filtering of bad quality cells.
#' @param velocyto.percent.mt Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.normalisation.method Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.selection.method Parameter for GEX analysis of the cell.count.matrices.
#' @param velocyto.deltaT Parameter for Velocyto analysis
#' @param velocyto.kCells Parameter for Velocyto analysis
#' @param velocyto.fit.quantile Parameter for Velocyto analysis
#' @param velocyto.kGenes Parameter for Velocyto analysis
#' @return If method=monocle3, the function returns a list element: [[1]] UMAP colored by Pseudotime, [[2]] Ridgeline plots showing the density of each of the top.N.clonotypes per cluster along pseudotime., [[3]] Gene expression plots highlighting the gene expression across pseudotime colored by transcriptional cluster, [[4]] Gene expression plots highlighting the gene expression across pseudotime colored by colotype. If method=velocyto, plots and RDS will be saved to velocyto.out.dir.
#' @export
#' @examples
#' \dontrun{
#' #----Method=monocle3----
#'
#' # Version 2
#'    vdj_repertoire_tcells <- VDJ_analyze(VDJ.out.directory =VDJ.out.directory.list,
#'     filter.1HC.1LC = T)
#'    gex_acute <- Platypus::GEX_automate(GEX.outs.directory.list = dir_to_gex[1:1],
#'     integration.method = "scale.data", mito.filter = 20, cluster.resolution = 0.5,
#'    VDJ.gene.filter = T)
#'
#' clonotyme_output <- VDJ_GEX_clonotyme(vdj.analyze.output = vdj_repertoire_tcells,
#'  gex.analyze.output = gex_acute, version="v2", exclude.clusters=7, highlight.genes="sell",
#'   colors = c("blue", "red", "black", "orange")
#' clonotyme_output[[4]]
#'
#' # Version 3
#'     VGM <-
#'     readRDS("C:/Users/rapha/Downloads/TEMPLATE_VDJ_GEX_mat_Bcells_r2_150521.rds")
#'
#' clonotyme_output <- VDJ_GEX_clonotyme(vdj.gex.matrix.output = VGM, version="v3",
#'  highlight.genes="sell", top.N.clonotypes = 1)
#'
#'
#' #---Method=velocyto----
#'
#' #Dropest Alignment: Run on EULER CLUSTER
#' #env2lmod
#' #module load gcc/4.8.5 python/3.7.4
#' #module load gcc/4.8.5 dropest/0.8.6
#' #module load gcc/4.8.5 r/4.0.2
#' #module load gcc/4.8.5 hdf5/1.10.1
#' #module load gcc/4.8.5 openmpi/4.0.2
#' #module load gcc/4.8.5 r/4.0.2
#' #bsub -W 2880 -R 'rusage[mem=20000]'
#' /cluster/home/rakuhn/dropEst/dropest -V -C 6000 -f -g
#' /cluster/scratch/rakuhn/mm10-2020-A/refdata-gex-mm10-2020-A/genes/genes.gtf
#' -c /cluster/home/rakuhn/dropEst/configs/10x.xml
#' /cluster/scratch/rakuhn/cellranger_v5/g1/outs/possorted_genome_bam.bam
#'
#' #Load required VDJ.analyze.output on EULER CLUSTER
#'
#' vdj_repertoire_tcells
#' <- readRDS("/cluster/home/rakuhn/RPII/vdj_repertoire_tcells.rds")
#' vdj_repertoire_tcells
#' <- head(vdj_repertoire_tcells,2)
#' #Only select the first two repertoires since we only want to analyze these two.
#'
#' # Load the two corresponding Dropest cell.count.matrices.rds
#'
#' dropest.output.list <- list()
#' dropest.output.list[[1]]
#' <- readRDS("/cluster/home/rakuhn/RPII/old_bam/gex1/cell.counts.matrices.rds")
#' dropest.output.list[[2]]
#' <- readRDS("/cluster/home/rakuhn/RPII/old_bam/gex2/cell.counts.matrices.rds")
#'
#' # Run Velocyto using Clonotyme
#'
#' VDJ_GEX_clonotyme(method = "velocyto", version = "v2",
#'  vdj.analyze.output = vdj_repertoire_tcells,
#'  dropest.output.list = dropest.output.list,
#'  top.N.clonotypes = 3, exclude.clusters = 8, highlight.genes = "sell",
#'  velocyto.gex.merged = T, velocyto.out.dir = ".", velocyto.save.rds = F)
#'
#'}


VDJ_GEX_clonotyme <- function(method,
                              version,
                              top.N.clonotypes,
                              vdj.gex.matrix.output,
                              vdj.analyze.output,
                              gex.automate.output,
                              exclude.clusters,
                              colors,
                              show.cells,
                              highlight.genes,
                              dropest.output.list,
                              velocyto.gex.merged,
                              velocyto.file.name,
                              velocyto.out.dir,
                              velocyto.save.rds,
                              velocyto.norm.scale.factor,
                              velocyto.n.variable.features,
                              velocyto.neighbor.dim,
                              velocyto.cluster.resolution,
                              velocyto.mds.dim,
                              velocyto.nCount_spliced,
                              velocyto.percent.mt,
                              velocyto.normalisation.method,
                              velocyto.selection.method,
                              velocyto.deltaT,
                              velocyto.kCells,
                              velocyto.fit.quantile,
                              velocyto.kGenes){ #

  if(missing(method)){method <- "monocle3"}
  if(missing(version)){version <- "v2"}
  if(missing(top.N.clonotypes)){top.N.clonotypes <- 2}

  #adding for variable CRAN consistency
  pseudotime <- NULL
  clonotype <- NULL
  p.size <- NULL
  p.shape <- NULL
  gene_short_name_toupper <- NULL
  expectation <- NULL
  percent.mt <- NULL



  #####################
  # Monocle3 Pipeline #
  #####################

  if(method == "monocle3"){

    print("Running Monocle3 pipeline")
    #Commented out for Platypus v3 namespace / are noted in DESCRIPTION Suggests
    #require(monocle3)
    #require(ggplot2)
    #require(SeuratWrappers)
    #require(ggridges)
    #require(dplyr)

    if(missing(show.cells)){show.cells <- T}

    if(version=="v3"){
      print("Platypus v3 is beeing used.")
      if(missing(vdj.gex.matrix.output)){print("Please provide Output of VDJ_GEX_matrix as an input through the vdj.gex.matrix.output parameter.")}
      if(!missing(vdj.analyze.output)){print("Please use vdj.gex.matrix output as input parameter when using version 3.")}
      if(!missing(gex.automate.output)){print("Please use vdj.gex.matrix output as input parameter when using version 3.")}

      number.of.samples <- length(stats::na.omit(unique(vdj.gex.matrix.output[[2]]$sample_id)))


      print("Reading Input...")
      Cell.Data.Set <- suppressWarnings(SeuratWrappers::as.cell_data_set(vdj.gex.matrix.output[[2]]))
      print(suppressMessages(monocle3::plot_cells(Cell.Data.Set)))
      SummarizedExperiment::colData(Cell.Data.Set)$barcodes <- rownames(SummarizedExperiment::colData(Cell.Data.Set)) #adding barcodes column to Cell.Data.Set

      print("Setting UMAP Partitions...")
      Cell.Data.Set@clusters[["UMAP"]]$partitions <- factor(Cell.Data.Set@clusters[["UMAP"]]$partitions, levels = c(levels(Cell.Data.Set@clusters[["UMAP"]]$partitions), "2"))
      if(!missing(exclude.clusters)){
        Cell.Data.Set@clusters[["UMAP"]]$partitions[which(Cell.Data.Set@clusters[["UMAP"]]$clusters==as.character(exclude.clusters))] <- "2"
      }
      SummarizedExperiment::colData(Cell.Data.Set)$partitions <- Cell.Data.Set@clusters$UMAP$partitions
      names(Cell.Data.Set@clusters[["UMAP"]]$partitions) <- names(Cell.Data.Set@clusters[["UMAP"]]$clusters)
      partition.plot <- monocle3::plot_cells(Cell.Data.Set, color_cells_by = "partitions")
      print(partition.plot)

      #create clonotype column
      print("Integrating Clonotypes...")
      SummarizedExperiment::colData(Cell.Data.Set)$clonotype <- "other"

      level.list <- c()
      for (k in 1:number.of.samples){
        for (i in 1:top.N.clonotypes){

          temp_barcodes <- vdj.gex.matrix.output[[1]][which(vdj.gex.matrix.output[[1]]$sample_id == paste0("s", k)),][which(vdj.gex.matrix.output[[1]]$clonotype_id_10x == paste0("clonotype", k)),]$barcode
          #temp_barcodes <- gsub(temp_barcodes, pattern = "-1",replacement = "")
          temp_barcodes <- temp_barcodes[which(temp_barcodes!="NA")]
          #print(paste0("Sample ", k))
          #print(paste0("Clonotype ", i))
          #print(temp_barcodes)
          index <- list()
          for(j in 1:length(temp_barcodes)){
            if(length(which(SummarizedExperiment::colData(Cell.Data.Set)$barcodes==temp_barcodes[j]))){
              index[j]<- which(SummarizedExperiment::colData(Cell.Data.Set)$barcodes==temp_barcodes[j])
            }
          }

          index <- unlist(index)
          #print(index)

          SummarizedExperiment::colData(Cell.Data.Set)$clonotype[index] <- paste0("Sample ",k," Clonotype ", i)
          level.list <- c(level.list, paste0("Sample ",k," Clonotype ", i))
        }
      }

    } else if(version=="v2"){
      print("Platypus v2 is beeing used.")
      if(missing(vdj.analyze.output)){print("Please provide Output of VDJ_analyze() as an input through the vdj.analyze.output parameter.")}
      if(missing(gex.automate.output)){print("Please provide Output of GEX_analyze() as an input through the gex.automate.output parameter.")}
      if(!missing(vdj.gex.matrix.output)){print("Please use vdj.analyze.output or gex.automate.output as Input parameters")}

      # Make new barcode column
      gex.automate.output[[1]]$barcodes <- colnames(gex.automate.output[[1]])
      print("filter out -1_2")
      gex.automate.output[[1]]$barcodes <- gsub(gex.automate.output[[1]]$barcodes,pattern = "-1_2",replacement = "")
      print("filter out -1_1")
      gex.automate.output[[1]]$barcodes <- gsub(gex.automate.output[[1]]$barcodes,pattern = "-1_1",replacement = "")
      print("filter out -1")
      gex.automate.output[[1]]$barcodes <- gsub(gex.automate.output[[1]]$barcodes,pattern = "-1",replacement = "")
      gex.automate.output[[1]]$barcode <- gex.automate.output[[1]]$barcodes

      number.of.samples <- length(unique(gex.automate.output[[1]]$sample))

      print("Reading Input...")
      Cell.Data.Set <- suppressWarnings(SeuratWrappers::as.cell_data_set(gex.automate.output[[1]]))
      print(suppressMessages(monocle3::plot_cells(Cell.Data.Set)))

      print("Setting UMAP Partitions...")
      Cell.Data.Set@clusters[["UMAP"]]$partitions <- factor(Cell.Data.Set@clusters[["UMAP"]]$partitions, levels = c(levels(Cell.Data.Set@clusters[["UMAP"]]$partitions), "2"))
      if(!missing(exclude.clusters)){
        Cell.Data.Set@clusters[["UMAP"]]$partitions[which(Cell.Data.Set@clusters[["UMAP"]]$clusters==as.character(exclude.clusters))] <- "2"
      }
      SummarizedExperiment::colData(Cell.Data.Set)$partitions <- Cell.Data.Set@clusters$UMAP$partitions
      names(Cell.Data.Set@clusters[["UMAP"]]$partitions) <- names(Cell.Data.Set@clusters[["UMAP"]]$clusters)
      partition.plot <- monocle3::plot_cells(Cell.Data.Set, color_cells_by = "partitions")
      print(partition.plot)

      #create clonotype column
      print("Integrating Clonotypes...")
      SummarizedExperiment::colData(Cell.Data.Set)$clonotype <- "other"
      level.list <- c()
      for (k in 1:number.of.samples){
        for (i in 1:top.N.clonotypes){

          temp_barcodes <- stringr::str_split(vdj.analyze.output[[k]]$barcodes[i],pattern=";")[[1]]
          temp_barcodes <- gsub(temp_barcodes, pattern = "-1",replacement = "")

          index <- list()
          for(j in 1:length(temp_barcodes)){
            if(length(which(SummarizedExperiment::colData(Cell.Data.Set)$barcodes==temp_barcodes[j]))){
              index[j]<- which(SummarizedExperiment::colData(Cell.Data.Set)$barcodes==temp_barcodes[j])
            }
          }
          index <- unlist(index)
          SummarizedExperiment::colData(Cell.Data.Set)$clonotype[index] <- paste0("Sample ",k," Clonotype ", i)
          level.list <- c(level.list, paste0("Sample ",k," Clonotype ", i))
        }
      }
    }else{
      print("Please specify whether you want to use platypus version v2 or v3")
    }

    print("Learning Graph...")
    Cell.Data.Set <- monocle3::learn_graph(Cell.Data.Set)

    print("Setting Trajectory Root...")
    Cell.Data.Set <- monocle3::order_cells(Cell.Data.Set)
    print("Calculating Pseudotime...")
    pseudotime.plot <- monocle3::plot_cells(Cell.Data.Set,
                                            color_cells_by = "pseudotime",
                                            label_cell_groups=FALSE,
                                            label_leaves=FALSE,
                                            label_branch_points=FALSE,
                                            graph_label_size=1.5)

    SummarizedExperiment::colData(Cell.Data.Set)$pseudotime <- pseudotime(Cell.Data.Set, reduction_method = "UMAP")

    print("Preparing Ridgeline Plots...")
    df <- as.data.frame(SummarizedExperiment::colData(Cell.Data.Set))
    df <- df[which(df$pseudotime!=Inf),]
    df2 <- df
    df2$clonotype <- "all"
    df <- df[which(df$clonotype!="other"),]
    df$clonotype <- as.factor(df$clonotype)
    df$clonotype <- ordered(df$clonotype, levels = c(unique(level.list), "all"))
    df <- rbind(df, df2)
    df$p.size <- 0.2
    df$p.size[which(df$clonotype=="all")] <- 0
    df$p.shape <- "circle"
    df$p.shape[which(df$clonotype=="all")] <- "square"

    if(missing(colors)){
      colors <- grDevices::rainbow(top.N.clonotypes*number.of.samples)
      cols <- c(colors, "grey")
    }else{
      if(length(colors)!=top.N.clonotypes*number.of.samples){print("Please provide enough colors for each clonotype in each sample you want to show (top.N.clonotypes*number.of.samples)")}
      cols <- c(colors, "grey")
    }

    print("Plotting Ridgeline Plots...")
    if(show.cells){
      ridgeline.plot <- ggplot2::ggplot(df, ggplot2::aes(x = pseudotime, y =clonotype,fill=clonotype, alpha=0.5))+ ggridges::geom_density_ridges(ggplot2::aes(point_size=p.size, point_shape=p.shape),scale=2, jittered_points = TRUE, point_alpha=1)+ggplot2::theme_classic()+ggplot2::scale_fill_manual(values=cols)+ggridges::scale_point_size_continuous(range = c(0, 2), guide = "none")+ggplot2::scale_discrete_manual(aesthetics = "point_shape", values = c(21, NA))
    }else{
      ridgeline.plot <- ggplot2::ggplot(df, ggplot2::aes(x = pseudotime, y =clonotype,fill=clonotype))+ ggridges::geom_density_ridges(scale=2, jittered_points = F, point_alpha=1)+ggplot2::theme_classic()+ggplot2::scale_fill_manual(values=cols)
    }
    print(ridgeline.plot)
    # #nolegend
    # ggplot(df, aes(x = pseudotime, y =clonotype,fill=clonotype, alpha=0.5))+ geom_density_ridges(aes(point_size=p.size, point_shape=p.shape),scale=2, jittered_points = TRUE, point_alpha=1)+theme_classic()+scale_fill_manual(values=cols)+scale_point_size_continuous(range = c(0, 2), guide = "none")+scale_discrete_manual(aesthetics = "point_shape", values = c(21, NA))+theme(legend.position = "none")

    print("---DONE---")

    # Do Gene-plots
    print("Preparing Gene Plots...")
    SummarizedExperiment::rowData(Cell.Data.Set)$gene_short_name <- rownames(SummarizedExperiment::rowData(Cell.Data.Set))
    Cell.Data.Set <- monocle3::estimate_size_factors(Cell.Data.Set)

    print(paste0("Selecting genes: ", unlist(highlight.genes)))
    fda.usc::fdata(Cell.Data.Set)$gene_short_name_toupper <- toupper(fda.usc::fdata(Cell.Data.Set)$gene_short_name)
    genes <- row.names(subset(fda.usc::fdata(Cell.Data.Set), gene_short_name_toupper %in% toupper(highlight.genes)))
    Cell.Data.Set_subset <- Cell.Data.Set[genes,]

    print("Plotting Genes in Pseudotime...")
    gene.plot.cluster <- monocle3::plot_genes_in_pseudotime(Cell.Data.Set_subset, trend_formula = "~ splines::ns(pseudotime, df=3)", min_expr = NULL, vertical_jitter = 0.2, color_cells_by = "seurat_clusters")
    print(gene.plot.cluster)

    gene.plot.cluster$data <- gene.plot.cluster$data[order(gene.plot.cluster$data$clonotype, decreasing=F),]

    gene.plot.clonotype <- ggplot2::ggplot(ggplot2::aes(pseudotime, expression), data = gene.plot.cluster$data) + ggplot2::geom_point(ggplot2::aes(color = clonotype, size=clonotype, alpha=clonotype))#, position = position_jitter(horizontal_jitter, vertical_jitter)
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::geom_line(ggplot2::aes(x = pseudotime, y = expectation),
                                                           data = gene.plot.cluster$data)
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::scale_y_log10() + ggplot2::facet_wrap(~feature_label, nrow = length(highlight.genes),
                                                                              ncol = 1, scales = "free_y")
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::expand_limits(y = c(0, 1))
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::ylab("Expression")
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::xlab("Pseudotime")
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::theme_minimal()
    gene.plot.clonotype <- gene.plot.clonotype + ggplot2::scale_color_manual(breaks = c(level.list, "other"), values=cols) + ggplot2::scale_size_manual(breaks = c(level.list, "other"), values=c(2,2,2,2,1.5)) + ggplot2::scale_alpha_manual(breaks = c(level.list, "other"), values=c(0.7,0.7,0.7,0.7,0.4))#+theme(legend.position = "bottom", legend.box="vertical", legend.margin=margin())
    print(gene.plot.clonotype)

    print("---DONE---")
    # gene.plot.clonotype + theme(legend.position = "none",strip.text.x = element_text(size = 14))


    output <- list()
    output[[1]] <- print(pseudotime.plot)
    output[[2]] <- print(ridgeline.plot)
    output[[3]] <- print(gene.plot.cluster)
    output[[4]] <- print(gene.plot.clonotype)
    return(output)
  }


  #####################
  # Velocyto Pipeline #
  # Run on Cluster    #
  #####################

  if(method=="velocyto"){
    print("Velocyto.R chosen as Pseudotime method. Please use DropEst aligned Intron and Exton reads as Input. Provide Cell.counts.matrices.rds as Input through dropest.output parameter.")

    if(missing(dropest.output.list)){print("Please provide Exon and Intron alignments by DropEst: cell.counts.matrices.rds required as Input.")}
    if(missing(velocyto.out.dir)){velocyto.out.dir <- "."}
    #require(Seurat)
    #require(SeuratWrappers)
    #require(ggplot2)
    #require(stringr)
    #require(scales)

    if(missing(velocyto.norm.scale.factor)){velocyto.norm.scale.factor <- 10000}
    if(missing(velocyto.n.variable.features)){velocyto.n.variable.features <- 2000}
    if(missing(velocyto.neighbor.dim)){velocyto.neighbor.dim <- 1:10}
    if(missing(velocyto.cluster.resolution)){velocyto.cluster.resolution <- 0.5}
    if(missing(velocyto.mds.dim)){velocyto.mds.dim <- 1:10}
    if(missing(velocyto.nCount_spliced)){velocyto.nCount_spliced <- 1000}
    if(missing(velocyto.percent.mt)){velocyto.percent.mt <- 2}
    if(missing(velocyto.normalisation.method)){velocyto.normalisation.method <- "LogNormalize"}
    if(missing(velocyto.selection.method)){velocyto.selection.method <- "vst"}
    if(missing(velocyto.deltaT)){velocyto.deltaT <- 1}
    if(missing(velocyto.kCells)){velocyto.kCells <- 25}
    if(missing(velocyto.fit.quantile)){velocyto.fit.quantile <- 0.02}
    if(missing(velocyto.kGenes)){velocyto.kGenes <- 1}

    if(missing(velocyto.file.name)){velocyto.file.name <- "Velocyto"}
    file.name <- velocyto.file.name


    switch(Sys.info()[['sysname']],
           Windows= {print("Windows system detected")
             operating.system <- "Windows"},
           Linux  = {print("Linux system detected")
             operating.system <- "Linux"},
           Darwin = {print("MAC system detected")
             operating.system <- "Darwin"})

    if(operating.system=="Windows"){
      print("Package Velocyto.R is now available for Windows. Please switch to UNIX based System.")
    }else{
      print("Loading Velocyto.R...")
      #require(velocyto.R) #Not possible, as we are not requiring this package in the DESCRIPTION file !

      print("Loading Functions...")
      # Defining Functions

      make_Seurat <- function(cell.count.mat, sample){
        gex <- SeuratObject::as.Seurat(x=cell.count.mat)
        gex[["spliced"]]<- SeuratObject::CreateAssayObject( counts = cell.count.mat$exon)
        gex[["unspliced"]] <- SeuratObject::CreateAssayObject( counts = cell.count.mat$intron)
        gex[["spanning"]] <- SeuratObject::CreateAssayObject( counts = cell.count.mat$spanning)
        gex$sample <- sample
        return(gex)
      }

      #####

      VDJ.filter <- function(S.object){
        holding_upper_gene_names <- toupper(rownames(S.object))
        antibody_gene_indices <- which(grepl((holding_upper_gene_names),pattern = "^IGHA")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHG")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHM")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHE")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGK")==F &
                                         grepl((holding_upper_gene_names),pattern = "^IGHV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^JCHAIN")==F&
                                         grepl((holding_upper_gene_names),pattern = "^IGL")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRAV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRAC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDC")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRGJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDD")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRDJ")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRBV")==F &
                                         grepl((holding_upper_gene_names),pattern = "^TRAJ")==F)
        S.object <- S.object[antibody_gene_indices,]
        return(S.object)
      }

      #####

      Velocyto_Gex_anaylze_prep <- function(obj, velocyto.norm.scale.factor, velocyto.n.variable.features, velocyto.neighbor.dim,velocyto.cluster.resolution, velocyto.mds.dim, velocyto.nCount_spliced, velocyto.percent.mt, velocyto.normalisation.method, velocyto.selection.method){

        nCount_spliced <- NULL

        obj <- VDJ.filter(obj)
        obj[["barcode"]]<- gsub(colnames(obj),pattern = "-1_1",replacement = "")
        obj[["barcode"]]<- gsub(obj[["barcode"]],pattern = "-1_1",replacement = "")
        obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-") + Seurat::PercentageFeatureSet(obj, pattern = "^mt-")

        obj <- subset(obj, subset=nCount_spliced>velocyto.nCount_spliced) # filter out bad quality cells according to: https://github.com/velocyto-team/velocyto.R/issues/131
        obj <- subset(obj, subset=percent.mt<velocyto.percent.mt)

        obj <- Seurat::NormalizeData(obj, normalization.method = velocyto.normalisation.method, scale.factor = velocyto.norm.scale.factor)
        obj <- Seurat::FindVariableFeatures(obj, selection.method = velocyto.selection.method, nfeatures = velocyto.n.variable.features)
        all.genes <- rownames(obj)
        obj <- Seurat::ScaleData(obj, features = all.genes)
        obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
        obj <- Seurat::FindNeighbors(obj, dims = velocyto.neighbor.dim)
        obj <- Seurat::FindClusters(obj, resolution = velocyto.cluster.resolution)
        obj <- Seurat::RunUMAP(obj, dims = velocyto.mds.dim)
        return(obj)
      }

      #####

      Velocyto_sample_run <- function(sample, sample.name, file.name, velocyto.out.dir, velocyto.save.rds, velocyto.deltaT, velocyto.kCells, velocyto.fit.quantile){
        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sample)))
        names(x = ident.colors) <- levels(x = sample)
        cell.colors <- ident.colors[Seurat::Idents(object = sample)]
        names(x = cell.colors) <- colnames(x = sample)
        sample <- SeuratWrappers::RunVelocity(object = sample, deltaT = velocyto.deltaT, kCells = velocyto.kCells, fit.quantile = velocyto.fit.quantile)

        if(velocyto.save.rds){saveRDS(gex, paste0(velocyto.out.dir,"/",file.name,"_",sample.name,"_post_velocyto.rds"))}
        return(sample)
      }

      #####

      Velocyto_sample_plots <- function(sample, sample.name, file.name, velocyto.out.dir){
        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sample)))
        names(x = ident.colors) <- levels(x = sample)
        cell.colors <- ident.colors[Seurat::Idents(object = sample)]
        names(x = cell.colors) <- colnames(x = sample)

        grDevices::pdf(paste0(velocyto.out.dir,"/",file.name,'.plot.pdf'), height=5, width=6)
        velocyto::show.velocity.on.embedding.cor(emb = SeuratObject::Embeddings(object = sample, reduction = "umap"), vel = SeuratObject::Tool(object = sample, slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = velocyto::ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
        grDevices::dev.off()

        Seurat::FeaturePlot(sample, highlight.genes, order=T)
        ggplot2::ggsave(paste0(velocyto.out.dir,"/",file.name,"_",sample.name,'.feature.plot.pdf'), height=5, width=6)

        Seurat::DimPlot(sample, reduction="umap")
        ggplot2::ggsave(paste0(velocyto.out.dir,"/",file.name,"_",sample.name,'.umap.pdf'), height=5, width=6)

        return(sample)
      }

      #####

      Velocyto_pseudotime_plots <- function(sample, file.name, sample.name, gene, velocyto.out.dir, velocyto.deltaT, velocyto.kCells, velocyto.fit.quantile, velocyto.kGenes){
        #https://nbisweden.github.io/single-cell_sib_scilifelab/session-trajectories/4_velocity.html


        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sample)))
        names(x = ident.colors) <- levels(x = sample)
        cell.colors <- ident.colors[Seurat::Idents(object = sample)]

        grDevices::pdf(paste0(velocyto.out.dir,"/",file.name,"_",sample.name,'_', gene ,'.gene.pdf'), height =7, width=25)
        fit.quantile<-0.02
        velocyto::gene.relative.velocity.estimates(
          SeuratObject::GetAssayData(sample, slot = "data", assay = "spliced"),
          SeuratObject::GetAssayData(sample, slot = "data", assay = "unspliced"),
          cell.emb = SeuratObject::Embeddings(sample, "umap"),
          deltaT=velocyto.deltaT,kCells = velocyto.kCells,kGenes=velocyto.kGenes,
          show.gene = gene,
          old.fit = SeuratObject::Tool(sample, slot = "RunVelocity"),
          do.par=T,
          #cell.colors = cell.colors
        )
        grDevices::dev.off()
        return(sample)
      }

      Velocyto_highlight_clonotypes <- function(sample, barcodes_vdj, clonotype.number, file.name, sample.name, colors){

        #Show both clonotypes on same UMAP

        # Make new barcode column
        sample$barcodes <- colnames(sample)
        print("filter out -1_2")
        sample$barcodes <- gsub(sample$barcodes,pattern = "-1_2",replacement = "")
        print("filter out -1_1")
        sample$barcodes <- gsub(sample$barcodes,pattern = "-1_1",replacement = "")
        print("filter out -1")
        sample$barcodes <- gsub(sample$barcodes,pattern = "-1",replacement = "")

        # Name all "Other" clonotype
        sample$clonotype <- rep("Other",length(colnames(sample)))


        print(length(barcodes_vdj))
        print(barcodes_vdj)
        for(k in 1:length(barcodes_vdj)){
          for(m in 1:length(barcodes_vdj[[k]])){
            barcodes[[k]][[m]] <- stringr::str_split(barcodes_vdj[[k]][m], ";")
            barcodes[[k]][[m]] <- unlist(barcodes[[k]][[m]])
            barcodes[[k]][[m]] <- gsub(barcodes[[k]][[m]],pattern = "-1",replacement = "")
          }
        }

        for(k in 1:length(barcodes)){
          for(m in 1:length(barcodes[[k]])){
            for(i in 1:length(barcodes[[k]][[m]])){
              sample$clonotype[which(sample$barcodes==barcodes[[k]][[m]][i])] <- paste0("Clonotype_",m ," Mouse_",k)
            }
          }
        }

        Seurat::Idents(sample) <- sample$clonotype

        ident.colors <- (scales::hue_pal())(n = length(x = levels(x = sample)))
        names(x = ident.colors) <- levels(x = sample)
        ident.colors[which(names(ident.colors)=="Other")]<- "#C0C0C0" #grey

        for(k in 1:length(barcodes)){
          for(m in 1:length(barcodes[[k]])){
            ident.colors[which(names(ident.colors)==paste0("Clonotype_", m," Mouse_", k))]<- colors[((k-1)*length(barcodes[[k]])+m)]
          }
        }

        t1 <- utils::head(ident.colors,1)
        t1 <- velocyto::ac(t1, alpha=0.1)
        t2 <- utils::tail(ident.colors,(length(ident.colors)-1))
        t2 <- velocyto::ac(t2, alpha=1)
        ident.colors <- append(t1,t2)
        cell.colors <- ident.colors[Seurat::Idents(object = sample)]

        names(x = cell.colors) <- colnames(x = sample)

        grDevices::pdf(paste0(velocyto.out.dir,"/",file.name,"_",sample.name,'.clonotypes.on.velocyto.plot.pdf'), height=5, width=6)
        velocyto::show.velocity.on.embedding.cor(emb = SeuratObject::Embeddings(object = sample, reduction = "umap"), vel = SeuratObject::Tool(object = sample, slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = cell.colors, cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0)
        grDevices::dev.off()

      }

      ####

      Velocyto_Wrapper <- function(gex, colors,  sample.name, clonotype.number, highlight.genes, file.name, velocyto.save.rds, exclude.clusters, velocyto.out.dir, vdj_repertoire, velocyto.norm.scale.factor, velocyto.n.variable.features, velocyto.neighbor.dim,velocyto.cluster.resolution, velocyto.mds.dim, velocyto.nCount_spliced, velocyto.percent.mt, velocyto.normalisation.method, velocyto.selection.method, velocyto.deltaT, velocyto.kCells, velocyto.fit.quantile, velocyto.kGenes){


        seurat_clusters <- NULL

        gex <- VDJ.filter(S.object = gex)
        gex <- Velocyto_Gex_anaylze_prep(obj=gex, velocyto.norm.scale.factor=velocyto.norm.scale.factor, velocyto.n.variable.features=velocyto.n.variable.features, velocyto.neighbor.dim=velocyto.neighbor.dim,velocyto.cluster.resolution=velocyto.cluster.resolution, velocyto.mds.dim=velocyto.mds.dim, velocyto.nCount_spliced=velocyto.nCount_spliced, velocyto.percent.mt=velocyto.percent.mt, velocyto.normalisation.method=velocyto.normalisation.method, velocyto.selection.method=velocyto.selection.method)
        if(velocyto.save.rds){saveRDS(gex, paste0(velocyto.out.dir,"/",file.name,"_",sample.name,"_pre_velocyto.rds"))}

        print("Remove Bcells")
        for(i in exclude.clusters){
          gex <- subset(gex, seurat_clusters!=i)
        }

        if(velocyto.save.rds){saveRDS(gex, paste0(velocyto.out.dir,"/",file.name,"_",sample.name,"_pre_velocyto_removedclusters.rds"))}

        gex <- Velocyto_sample_run(sample = gex, sample.name=sample.name, file.name=file.name, velocyto.out.dir = velocyto.out.dir, velocyto.save.rds = velocyto.save.rds, velocyto.deltaT = velocyto.deltaT, velocyto.kCells = velocyto.kCells, velocyto.fit.quantile=velocyto.fit.quantile)
        gex <- Velocyto_sample_plots(sample = gex, sample.name=sample.name, file.name=file.name, velocyto.out.dir = velocyto.out.dir)
        gex <- Velocyto_pseudotime_plots(sample = gex, file.name=file.name, sample.name=sample.name, gene = highlight.genes, velocyto.out.dir = velocyto.out.dir, velocyto.deltaT = velocyto.deltaT, velocyto.kCells = velocyto.kCells, velocyto.fit.quantile=velocyto.fit.quantile, velocyto.kGenes = velocyto.kGenes)

        gex <- Velocyto_highlight_clonotypes(colors= colors, sample = gex, barcodes_vdj = vdj_repertoire, clonotype.number = clonotype.number, file.name = file.name, sample.name = sample.name)

        return(gex)
      }

      #####
      #
      #_____________________START VELOCYTO__________________________
      #
      # Read VDJ-repertoire and isolate barcodes to highlight
      barcodes <- list()

      if(version == "v2"){
        print("Reading in V(D)J repertoire from vdj.analyze.output...")
        vdj_repertoire_tcells <- vdj.analyze.output

        #get Barcodes to highlight
        number.of.samples <- length(vdj_repertoire_tcells)
        for(i in 1:length(vdj_repertoire_tcells)){
          barcodes[[i]] <- utils::head(vdj_repertoire_tcells[[i]]$barcodes, top.N.clonotypes)
        }



      }else if(version == "v3"){
        print("Reading in V(D)J repertoire from vdj.gex.matrix.output...")
        vdj_repertoire_tcells <- vdj.gex.matrix.output[[1]]

        #get Barcodes to highlight
        number.of.samples <-length(dropest.output.list)
        for(i in 1:length(dropest.output.list)){
          barcodes[[i]] <- utils::head(dplyr::filter(top.N.clonotypes, vdj_repertoire_tcells$sample_id == paste0('s', i))$barcodes)
          #Does not work yet...
        }
      }else{
        print("Please specify Platypus version you want to use as either v2 or v3.")
      }


      #define colors
      if(missing(colors)){
        colors <- grDevices::rainbow(top.N.clonotypes*number.of.samples)
        cols <- c(colors, "grey")
      }else{
        if(length(colors)!=top.N.clonotypes*number.of.samples){print("Please provide enough colors for each clonotype in each sample you want to show (top.N.clonotypes*number.of.samples)")}
        cols <- c(colors, "grey")
      }


      print("Reading in DropEst cell.count.matrices.rds...")

      sample_num <- length(dropest.output.list)

      gex <- list()
      for(i in 1:sample_num){
        gex[[i]] <- make_Seurat(cell.count.mat = dropest.output.list[[i]], sample = paste0("s", i))
      }

      if(velocyto.gex.merged){
        merged_gex <- gex[[1]]
        for(i in 2:length(gex)){
          merged_gex <- merge(merged_gex, gex[[i]]) #Merge GEX together to plot Clonotypes on same UMAP.

          Velocyto_Wrapper(colors=cols, gex = merged_gex, clonotype.number = top.N.clonotypes, highlight.genes = highlight.genes, file.name = file.name, velocyto.save.rds = velocyto.save.rds, exclude.clusters = exclude.clusters, velocyto.out.dir = velocyto.out.dir, vdj_repertoire = barcodes, sample.name="", velocyto.norm.scale.factor=velocyto.norm.scale.factor, velocyto.n.variable.features=velocyto.n.variable.features, velocyto.neighbor.dim=velocyto.neighbor.dim,velocyto.cluster.resolution=velocyto.cluster.resolution, velocyto.mds.dim=velocyto.mds.dim, velocyto.nCount_spliced=velocyto.nCount_spliced, velocyto.percent.mt=velocyto.percent.mt, velocyto.normalisation.method=velocyto.normalisation.method, velocyto.selection.method=velocyto.selection.method, velocyto.deltaT=velocyto.deltaT, velocyto.kCells=velocyto.kCells, velocyto.fit.quantile = velocyto.fit.quantile, velocyto.kGenes= velocyto.kGenes)
        }

      }else{
        for(i in length(gex)){
          Velocyto_Wrapper(colors = cols, gex = gex[[i]], clonotype.number = top.N.clonotypes, highlight.genes = highlight.genes, file.name = file.name, velocyto.save.rds = velocyto.save.rds, exclude.clusters = exclude.clusters, velocyto.out.dir = velocyto.out.dir, vdj_repertoire = barcodes[[i]], sample.name= paste0("S", i), velocyto.norm.scale.factor=velocyto.norm.scale.factor, velocyto.n.variable.features=velocyto.n.variable.features, velocyto.neighbor.dim=velocyto.neighbor.dim,velocyto.cluster.resolution=velocyto.cluster.resolution, velocyto.mds.dim=velocyto.mds.dim, velocyto.nCount_spliced=velocyto.nCount_spliced, velocyto.percent.mt=velocyto.percent.mt, velocyto.normalisation.method=velocyto.normalisation.method, velocyto.selection.method=velocyto.selection.method, velocyto.deltaT=velocyto.deltaT, velocyto.kCells=velocyto.kCells, velocyto.fit.quantile = velocyto.fit.quantile, velocyto.kGenes = velocyto.kGenes)
        }
      }

    }
  }

}

