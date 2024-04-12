#'Makes a Circos plot from the VDJ_analyze output. Connects the V gene with the corresponding J gene for each clonotype.
#' @param VGM The output of the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]]) has to be supplied. For Platypus v2: The output of the VDJ_GEX_integrate function (Platypus platypus.version v2). A list of data frames for each sample containing the clonotype information and cluster membership information.
#' @param VDJ.or.VJ Determines whether to plot the V J gene pairing of the alpha or beta chain. "VDJ", "VJ" or "both" as possible inputs. Default: "both".
#' @param label.threshold Genes are only labeled if the count is larger then the label.threshold. By default all label.threshold = 0 (all genes are labeled).
#' @param c.count.label Boolean, lets the user decide if the gene and count labels should be plotted or not. Default = T.
#' @param c.count.label.size Determines the font size of the gene labels. By default the font size for count labels is 0.6.
#' @param gene.label Boolean, lets the user decide if the gene labels should be plotted or not.
#' @param gene.label.size Determines the font size of the gene labels. By default the labelsize is automatically adjusted to 0.7 for labels with two or less digits, 0.6 for labels between 2 and 6 digits, and 0.4 for all longer labels. A manually defined font size will be the same for all labels!
#' @param platy.theme Allows plotting in the new "pretty" theme or the older "spiky" theme without group labels and radial arrangement of gene.labels. Default = "pretty".
#' @param arr.col Data.frame with three columns where the first two indicate the names of genes, clonotypes or clusters to be connected, and the third corresponds to the color of the arrow. Default set to data.frame(c("dummy.clonotype"), c("dummy.cluster"), c("dummy.color")), so no arrow is drawn.
#' @param arr.direction Either 1 or -1 and determines the direction of the arrow. Default=1.
#' @param c.threshold Only clonotypes are considered with a frequency higher then c.threshold. Allows to filter for only highly expanded clonotypes.
#' @param topX Filters for the top X clonotypes and only plots the respective gene combinations or cluster memberships.
#' @param clonotype.per.gene.threshold How many clonotypes are required to plot a sector for a gene. Filters the rows and colums of the final adjacency matrix.
#' @param filter1H1L Whether to filter the input VGM in "v3" to only include cells with 1 VDJ and 1 VJ chain. Defaults to TRUE
#' @param cell.level Logical, defines whether weight of connection should be based on number of clonotypes of number of cells. Default: number of clonotypes.
#' @param clonotype.column Which column in VGM contains the clonotyping information? Default="clonotype_id_10X".
#' @param platypus.version Which platypus.version of platypus is being used. Default = v3. Set to v3 if VDJ_GEX_matrix.output[[1]] is used
#' @return Returns a circos plot and a list object with the following elememts for N samples: [[1 to N]] The first N listelements corresponds to the recorded circos plots for N beeing the number or samples in the VGM. Since Circlize uses the R base plotting funciton, this is not a ggplot object but can still be replotted by calling the first list element. [[N+1]] Adjacency matrix forwarded to VDJ_circos(). This Matrix contains the counts and can be used for manual replotting using VDJ_circos directly. [[N+2]] Contains a named list with colors for each connection drawn and can be used for manual replotting using VDJ_circos directly. [[N+3]] Contains a named list with grouping information and can be used for manual replotting using VDJ_circos directly.
#' @export
#' @examples
#' \donttest{
#'  usage_circos_VDJVJ <- VDJ_VJ_usage_circos(Platypus::small_vgm[[1]])
#'  usage_circos_VDJVJ[[1]]
#'}
#'


VDJ_VJ_usage_circos <- function(VGM,
                                VDJ.or.VJ,
                                label.threshold,
                                cell.level,
                                c.threshold,
                                clonotype.per.gene.threshold,
                                c.count.label,
                                c.count.label.size,
                                platypus.version,
                                filter1H1L,
                                gene.label,
                                gene.label.size,
                                arr.col,
                                arr.direction,
                                topX,
                                platy.theme,
                                clonotype.column){

  if(missing(label.threshold)){label.threshold <- 0}
  if(missing(VDJ.or.VJ)){VDJ.or.VJ <- "both"}
  if(missing(cell.level)){cell.level <- FALSE}
  if(missing(c.threshold)){c.threshold <- 0}
  if(missing(clonotype.per.gene.threshold)){clonotype.per.gene.threshold <- 0}
  if(missing(c.count.label)){c.count.label <- TRUE}
  if(missing(c.count.label.size)){c.count.label.size <- 0.6}
  if(missing(platypus.version)){platypus.version <- "v3"}
  if(missing(filter1H1L)){filter1H1L <- TRUE}
  if(missing(gene.label)){gene.label <- TRUE}
  if(missing(gene.label.size)){gene.label.size <- "undef"}
  if(missing(arr.col)){arr.col <- data.frame(c("dummy1"), c("dummy2"), c(""))}
  if(missing(arr.direction)){arr.direction <- 1}
  if(missing(topX)){topX <- "all"}
  if(missing(platy.theme)){platy.theme <- "pretty"}
  if(missing(clonotype.column)){clonotype.column <- "clonotype_id_10x"}
  clonotype.frequency <- paste0("clonotype_frequency_", stringr::str_split(clonotype.column, pattern="_")[[1]][3])


  clonotypes_topX <- NULL
  VDJ.GEX.matrix <- NULL
  VDJ.GEX_list <- NULL
  clonotype <- NULL
  TRAV <- NULL
  TRBV <- NULL
  TRAJ <- NULL
  TRBJ <- NULL
  TRV <- NULL
  TRJ <- NULL
  c <- NULL
  d <- NULL
  e <- NULL
  f <- NULL
  Vgene_usage_matrix <- NULL
  dummy_alpha_df <- NULL
  dummy_beta_df <- NULL
  dummy <- NULL
  nm <- NULL
  group <- NULL
  grid.col <- NULL
  plot <- NULL
  circos.recorded <- NULL
  clonotype.frequency <- NULL


  #naming compatibility
  VDJ.GEX.matrix <- list()
  VDJ.GEX.matrix[[1]] <- VGM
  VGM <- NULL

  # If new version with VDJ_GEX_matric output should be used
  if(platypus.version=="v3"){

      if(cell.level == FALSE){
        message("WARNING: If clonotype strategy is not based on unique V or J genes per clonotype, this setting [cell.level=F] might be questionable. One clonotype might then be represented in several Circos connections between V or J genes. The names of genes of simulatneously used chains will be pasted together.")
        message("WARNING: If Circos plotting error occurs: Maybe your `gap.degree` is too large so that there is no space to allocate sectors -> You might want to increase clonotype.per.gene.threshold to reduce number of sectors in your Circos plots")

      }

      message(paste("Chosen clonotype column: ", clonotype.column))
      clonotype.frequency <- paste0("clonotype_frequency_", stringr::str_split(clonotype.column, pattern="_")[[1]][3])
      message(paste("Chosen clonotype.frequency column: ", clonotype.frequency))

      #filter for 1H1L
      if(filter1H1L==TRUE){
        VDJ.GEX.matrix[[1]]<-VDJ.GEX.matrix[[1]][which((VDJ.GEX.matrix[[1]]$VDJ_chain_count==1)&(VDJ.GEX.matrix[[1]]$VJ_chain_count==1)),]
      }

      #filter out clonotypes with less then c.threshold cells
      if(clonotype.frequency %in% colnames(VDJ.GEX.matrix[[1]])){ #check if 10x frequency already exists... will only be created after reclonotyping.
        VDJ.GEX.matrix[[1]] <-VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]][[clonotype.frequency]] >= c.threshold),]
        message(paste("Chosen clonotype.frequency column: ", clonotype.frequency))
      }else{
        VDJ.GEX.matrix[[1]] <-VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]][["clonotype_frequency"]] >= c.threshold),]
        message(paste("Chosen clonotype.frequency column: ", "clonotype_frequency"))
      }

      #replace empty fields "" with "None"
        VDJ.GEX.matrix[[1]]$VJ_vgene[which(VDJ.GEX.matrix[[1]]$VJ_vgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VJ_jgene[which(VDJ.GEX.matrix[[1]]$VJ_jgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VDJ_vgene[which(VDJ.GEX.matrix[[1]]$VDJ_vgene == "")] <- "None"
        VDJ.GEX.matrix[[1]]$VDJ_jgene[which(VDJ.GEX.matrix[[1]]$VDJ_jgene == "")] <- "None"

      TRBV <- c()
      TRAV <- c()
      TRBJ <- c()
      TRAJ <- c()

      #split up into samples to plot individually
      VDJ.GEX_list <- list()
      for (i in 1:length(unique(VDJ.GEX.matrix[[1]]$sample_id))){
        VDJ.GEX_list[[i]] <- VDJ.GEX.matrix[[1]][which(VDJ.GEX.matrix[[1]]$sample_id== unique(VDJ.GEX.matrix[[1]]$sample_id)[i]),]
      }

      # filter topX clonotypes

      if(topX != "all"){
        for(k in 1:length(VDJ.GEX_list)){
          clonotypes_topX <- names(utils::head(sort(table(VDJ.GEX_list[[k]][[clonotype.column]]),decreasing = TRUE),topX))

          #filter and keep only cells of topX clonotypes
          VDJ.GEX_list[[k]] <- VDJ.GEX_list[[k]][which(VDJ.GEX_list[[k]][[clonotype.column]] %in% clonotypes_topX),]
        }
      }


      for (i in 1:length(VDJ.GEX_list)){

        VDJ.GEX_list[[i]]$alpha_VJ_gene <- paste(VDJ.GEX_list[[i]]$VJ_vgene, VDJ.GEX_list[[i]]$VJ_jgene, sep = "_")
        VDJ.GEX_list[[i]]$beta_VJ_gene <- paste(VDJ.GEX_list[[i]]$VDJ_vgene, VDJ.GEX_list[[i]]$VDJ_jgene, sep = "_")

        #get all V and J genes across all samples
        c <- unique(VDJ.GEX_list[[i]]$VJ_vgene)
        TRBV <- append(TRBV,c)
        d <- unique(VDJ.GEX_list[[i]]$VDJ_vgene)
        TRAV <- append(TRAV,d)
        e <- unique(VDJ.GEX_list[[i]]$VJ_jgene)
        TRBJ<- append(TRBJ,e)
        f <- unique(VDJ.GEX_list[[i]]$VDJ_jgene)
        TRAJ <- append(TRAJ,f)
      }
      #make two branches for matrix V vs J genes
      if(VDJ.or.VJ=="both"){
        TRV <- append(TRAV, TRBV)
        TRJ <- append(TRAJ, TRBJ)
      }else if(VDJ.or.VJ=="VDJ"){
        TRV <- TRAV
        TRJ <- TRAJ
      }else if(VDJ.or.VJ=="VJ"){
        TRV <- TRBV
        TRJ <- TRBJ
      }else{
        warning("Please specify VDJ.or.VJ or leave empty to plot both.")
      }

      # create matrix Vgenes vs Jgenes
      Vgene_usage_matrix <- list()
      dummy_alpha_df <- list()
      dummy_beta_df <- list()

      for (k in 1:length(VDJ.GEX_list)){
        #Create a matrix with rows being heavy chain and columns being light chain v genes
        Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRV)), ncol = length(unique(TRJ)))

        #give the row names and col names
        rownames(Vgene_usage_matrix[[k]]) <- unique(TRV)
        colnames(Vgene_usage_matrix[[k]]) <- unique(TRJ)
      }

      for (k in 1:length(VDJ.GEX_list)){
        message(paste0("Processing sample ", k))
        #create dummy df which will contain the counts for each combination
        if(cell.level==TRUE){
          dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$alpha_VJ_gene))
          colnames(dummy_alpha_df[[k]]) <- c("vjgene", "count")
          dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX_list[[k]]$beta_VJ_gene))
          colnames(dummy_beta_df[[k]]) <- c("vjgene", "count")

        }else{

          dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype.column]],VDJ.GEX_list[[k]]$alpha_VJ_gene, sep="/and/")))
          colnames(dummy) <- c("pasted")
          dummy$clonotype <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)
          dummy$gene <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)[,2]
          dummy_alpha_df[[k]] <- as.data.frame(table(dummy$gene))
          colnames(dummy_alpha_df[[k]]) <- c("vjgene", "count")

          dummy <- as.data.frame(unique(paste(VDJ.GEX_list[[k]][[clonotype.column]],VDJ.GEX_list[[k]]$beta_VJ_gene, sep="/and/")))
          colnames(dummy) <- c("pasted")
          dummy$clonotype <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)
          dummy$gene <- stringr::str_split_fixed(dummy$pasted, "/and/", 2)[,2]
          dummy_beta_df[[k]] <- as.data.frame(table(dummy$gene))
          colnames(dummy_beta_df[[k]]) <- c("vjgene", "count")

        }



        #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$vgenes column
        for (i in 1:nrow(Vgene_usage_matrix[[k]])){
          for (j in 1:ncol(Vgene_usage_matrix[[k]])){
            if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_alpha_df[[k]]$vjgene){
              Vgene_usage_matrix[[k]][i,j] <- dummy_alpha_df[[k]][which(dummy_alpha_df[[k]]$vjgene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
            } else {
              Vgene_usage_matrix[[k]][i,j] <- 0
            }

          }
        }
        #go elementwise in the matrix and count the occurancies of each combination in the VDJ.clonotype.output$jgenes column
        for (i in 1:nrow(Vgene_usage_matrix[[k]])){
          for (j in 1:ncol(Vgene_usage_matrix[[k]])){
            if (paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j]) %in% dummy_beta_df[[k]]$vjgene){
              Vgene_usage_matrix[[k]][i,j] <- dummy_beta_df[[k]][which(dummy_beta_df[[k]]$vjgene == paste0(rownames(Vgene_usage_matrix[[k]])[i], "_",colnames(Vgene_usage_matrix[[k]])[j])),"count"]
            }
            #no else clause. would overwrite alpha chain genes in matrix    }
          }
        }

      }



      # group by TRAV;TRBV; TRAJ and TRBJ
      nm = unique(unlist(dimnames(Vgene_usage_matrix[[1]])))
      group = structure(stringr::str_sub(nm, 1,4), names = nm)
      group = factor(group)
      #set colors: same color across all samples
      grid.col <- stats::setNames(grDevices::rainbow(length(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]])))),sample(union(rownames(Vgene_usage_matrix[[1]]), colnames(Vgene_usage_matrix[[1]]))))

      # filter out genes which have fewer clonotypes then clonotype.per.gene.threshold
      for (k in 1:length(Vgene_usage_matrix)){
        Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][which(rowSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold),]
        Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][,which(colSums(Vgene_usage_matrix[[k]])>=clonotype.per.gene.threshold)]

      }

  }else if (platypus.version=="v2"){

  # In Case old platypus.version of Platypus is being used: VDJ_analyze output is the input

    for (i in 1:length(VDJ.GEX.matrix)) {
      VDJ.GEX.matrix[[i]] <- VDJ.GEX.matrix[[i]][which(VDJ.GEX.matrix[[i]]$frequency >=
                                                                 c.threshold), ]
    }
    TRBV <- c()
    TRAV <- c()
    TRBJ <- c()
    TRAJ <- c()
    for (i in 1:length(VDJ.GEX.matrix)) {
      VDJ.GEX.matrix[[i]]$alpha_VJ_gene <- paste(VDJ.GEX.matrix[[i]]$HC_vgene,
                                                     VDJ.GEX.matrix[[i]]$HC_jgene, sep = "_")
      VDJ.GEX.matrix[[i]]$beta_VJ_gene <- paste(VDJ.GEX.matrix[[i]]$LC_vgene,
                                                    VDJ.GEX.matrix[[i]]$LC_jgene, sep = "_")
      c <- unique(VDJ.GEX.matrix[[i]]$HC_vgene)
      TRBV <- append(TRBV, c)
      d <- unique(VDJ.GEX.matrix[[i]]$LC_vgene)
      TRAV <- append(TRAV, d)
      e <- unique(VDJ.GEX.matrix[[i]]$HC_jgene)
      TRBJ <- append(TRBJ, e)
      f <- unique(VDJ.GEX.matrix[[i]]$LC_jgene)
      TRAJ <- append(TRAJ, f)
    }
    if (VDJ.or.VJ == "both") {
      TRV <- append(TRAV, TRBV)
      TRJ <- append(TRAJ, TRBJ)
    }
    else if (VDJ.or.VJ == "VDJ") {
      TRV <- TRAV
      TRJ <- TRAJ
    }
    else if (VDJ.or.VJ == "VJ") {
      TRV <- TRBV
      TRJ <- TRBJ
    }
    else {
      warning("Please specify VDJ.or.VJ or leave empty to plot both.")
    }
    Vgene_usage_matrix <- list()
    dummy_alpha_df <- list()
    dummy_beta_df <- list()
    for (k in 1:length(VDJ.GEX.matrix)) {
      Vgene_usage_matrix[[k]] <- matrix(nrow = length(unique(TRV)),
                                        ncol = length(unique(TRJ)))
      rownames(Vgene_usage_matrix[[k]]) <- unique(TRV)
      colnames(Vgene_usage_matrix[[k]]) <- unique(TRJ)
    }
    for (k in 1:length(VDJ.GEX.matrix)) {
      if (cell.level == FALSE) {
        dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_VJ_gene))
        colnames(dummy_alpha_df[[k]]) <- c("vjgene",
                                           "count")
        dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$beta_VJ_gene))
        colnames(dummy_beta_df[[k]]) <- c("vjgene",
                                          "count")
      }
      else {
        dummy_alpha_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$alpha_VJ_gene))
        colnames(dummy_alpha_df[[k]]) <- c("vjgene",
                                           "clonotype_counts")
        for (i in 1:nrow(dummy_alpha_df[[k]])) {
          dummy_alpha_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$alpha_VJ_gene ==
                                                                              dummy_alpha_df[[k]]$vjgene[i]), "frequency"])
        }
        dummy_beta_df[[k]] <- as.data.frame(table(VDJ.GEX.matrix[[k]]$beta_VJ_gene))
        colnames(dummy_beta_df[[k]]) <- c("vjgene",
                                          "clonotype_counts")
        for (i in 1:nrow(dummy_beta_df[[k]])) {
          dummy_beta_df[[k]]$count[i] <- sum(VDJ.GEX.matrix[[k]][which(VDJ.GEX.matrix[[k]]$beta_VJ_gene ==
                                                                             dummy_beta_df[[k]]$vjgene[i]), "frequency"])
        }
      }
      for (i in 1:nrow(Vgene_usage_matrix[[k]])) {
        for (j in 1:ncol(Vgene_usage_matrix[[k]])) {
          if (paste0(rownames(Vgene_usage_matrix[[k]])[i],
                     "_", colnames(Vgene_usage_matrix[[k]])[j]) %in%
              dummy_alpha_df[[k]]$vjgene) {
            Vgene_usage_matrix[[k]][i, j] <- dummy_alpha_df[[k]][which(dummy_alpha_df[[k]]$vjgene ==
                                                                         paste0(rownames(Vgene_usage_matrix[[k]])[i],
                                                                                "_", colnames(Vgene_usage_matrix[[k]])[j])),
                                                                 "count"]
          }
          else {
            Vgene_usage_matrix[[k]][i, j] <- 0
          }
        }
      }
      for (i in 1:nrow(Vgene_usage_matrix[[k]])) {
        for (j in 1:ncol(Vgene_usage_matrix[[k]])) {
          if (paste0(rownames(Vgene_usage_matrix[[k]])[i],
                     "_", colnames(Vgene_usage_matrix[[k]])[j]) %in%
              dummy_beta_df[[k]]$vjgene) {
            Vgene_usage_matrix[[k]][i, j] <- dummy_beta_df[[k]][which(dummy_beta_df[[k]]$vjgene ==
                                                                        paste0(rownames(Vgene_usage_matrix[[k]])[i],
                                                                               "_", colnames(Vgene_usage_matrix[[k]])[j])),
                                                                "count"]
          }
        }
      }
    }
    nm = unique(unlist(dimnames(Vgene_usage_matrix[[1]])))
    group = structure(stringr::str_sub(nm, 1, 4), names = nm)
    grid.col <- stats::setNames(grDevices::rainbow(length(union(rownames(Vgene_usage_matrix[[1]]),
                                              colnames(Vgene_usage_matrix[[1]])))), sample(union(rownames(Vgene_usage_matrix[[1]]),
                                                                                                 colnames(Vgene_usage_matrix[[1]]))))
    for (k in 1:length(Vgene_usage_matrix)) {
      Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][which(rowSums(Vgene_usage_matrix[[k]]) >=
                                                                 clonotype.per.gene.threshold), ]
      Vgene_usage_matrix[[k]] <- Vgene_usage_matrix[[k]][,
                                                         which(colSums(Vgene_usage_matrix[[k]]) >= clonotype.per.gene.threshold)]
    }


    #--- Forward Vgene_usage_matrix to Plotting function
  }else{
    warning("Please specify platypus platypus.version as either v2 or v3.")
  }
  plot <- list()
  for (i in 1:length(Vgene_usage_matrix)){
    VDJ_circos(Vgene_usage_matrix[[i]], group = group, grid.col=grid.col, label.threshold = label.threshold, c.count.label=c.count.label, gene.label = gene.label, gene.label.size = gene.label.size, c.count.label.size=c.count.label.size, arr.col = arr.col, arr.direction = arr.direction, platy.theme = platy.theme)
    circos.recorded <- grDevices::recordPlot()
    plot[[i]] <- circos.recorded
  }

  plot[[i+1]] <- Vgene_usage_matrix
  plot[[i+2]] <- grid.col
  plot[[i+3]] <- group
  return(plot)
}
