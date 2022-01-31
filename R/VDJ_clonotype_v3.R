#' Updated clonotyping function based on implications for cells with different chain numbers than 1 VDJ 1 VJ chains.
#'@description This function offers two types of hierarchical clonotyping. The hierarchical option "single.chains" only merges cell with a single chain into clonotypes composed of cells with 1 VDJ 1 VJ chain. This is based on the assumption, that during mRNA capture and RT-PCR in GEMs, not all transcripts are captured and therefore cells may result missing a VDJ or VJ chain.
#' The hierarchical option "double.and.single.chains" is based on the assumption, that cells with 1 VDJ and 2 VJ chains exist. For a review of the work concerning such cells as well as 2 VDJ 1 VJ cells please consult: https://doi.org/10.4049/jimmunol.1800904. The user may set a threshold of occurrence number above which cells with 1 VDJ 2 VJ chains are considered to be true and other cells with 1 VDJ 1 VJ, 1 VDJ 0 VJ and 0 VDJ 1 VDJ may be merged into the same clonotype by the strategy provided by the user. Cells with 2 VDJ chains are currently not considered in this process, as these are reported to be much rarer and, if appearing in the dataset are more likely to be doublets.
#' We advice the user to carefully examine the output after hierarchical clonotyping before proceeding with further analysis.
#' @param VDJ For platypus v2 output from VDJ_analyze function. This should be a list of clonotype dataframes, with each list element corresponding to a single VDJ repertoire. For platypus v3 VDJ output from the VDJ_GEX_matrix function (VDJ_GEX_matrix.output[[1]])
#' @param clone.strategy (Updated keywords, previous format is also functional) String describing the clonotyping strategy. Possible options are 10x.default, cdr3.nt, cdr3.aa, VDJJ.VJJ, VDJJ.VJJ.cdr3length, VDJJ.VJJ.cdr3length.cdr3.homology, VDJJ.VJJ.cdr3length.VDJcdr3.homology, cdr3.homology, VDJcdr3.homology. cdr3.aa will convert the default cell ranger clonotyping to amino acid based. 'VDJJ.VJJ' groups B cells with identical germline genes (V and J segments for both heavy chain and light chain. Those arguments including 'cdr3length' will group all sequences with identical VDJ and VJ CDR3 sequence lengths. Those arguments including 'cdr3.homology' will additionally impose a homology requirement for CDRH3 and CDRL3 sequences.'CDR3.homology',or 'CDRH3.homology' will group sequences based on homology only (either of the whole CDR3 sequence or of the VDJ CDR3 sequence respectively).
#' All homology calculations are performed on the amino acid level.
#' @param homology.threshold Numeric value between 0 and 1 corresponding to the homology threshold forn the clone.strategy arguments that require a homology threshold. Default value is set to 70 percent sequence homology. For 70 percent homology, 0.3 should be supplied as input.
#' @param hierarchical Character. Defaults to "none". This is an extention specifically for cells with aberrant numbers of chains (i.e. 0VDJ 1VJ, 1VDJ 0VJ, 0VDJ 2VJ, 2VDJ 0VJ). Cells with 2VDJ 2VJ are filtered out as these are most likely doublets.
#' If set to "none" aberrant cells are assigned to their own clonotypes.
#' If set to "single.chains" the function will proceed in two steps: 0. Prefiltering: cells with 2 VDJ 2 VJ chains as well as cells with 2 VDJ and any number of VJ chains are filtered out. 1. define clonotypes classically with all cells containing exactly 1VDJ 1VJ chains. 2. For cells with only a single chain (either VDJ or VJ), check if any clone exists, which matches the clonotyping criteria for this chain. If true, add this cell to that clone. If false, create a new clone containing that cell. In case that more than 1 existing clone matches the aberrant cell, the cell is assigned to the most frequent existing clone. Two reasons are behind this decision: 2.1. The aberrant cells is numerically more likely to be a part of the more frequent existing clone. 2.2 In case of a wrong assignment, the effect of the error is lower, if an already expanded clone is increase by one count, rather than a existing non-expanded clone being assigned a second entry and thereby resulting as expanded. Cells with three chains are assigned to their own clonotypes
#' If set to "double.and.single.chains" the function will proceed as if set to "single.chains" but include two more steps
#' 3. Check the frequency of each cell 1 VDJ 2 VJ chain exact clone (by exact nucleotide CDR3 matching). Only if this count exceeds the triple.chain.count.threshold, the clone is used as a "hub clone". This protects from merging clonotypes on the basis of rare doublets.
#' 4. Merge existing clonotypes into the 1 VDJ 2 VJ clonotypes as they match with the assumption that e.g. a cell with 1 VDJ 1 VJ is part of that same clonotype, but missing a VJ chain due to stochastical sampling
#' @param triple.chain.count.threshold Minimal occurrance frequency for any cell with more than 2 of either VDJ or VJ chain (e.g. 2 VDJ 1 VJ) for it to be considered as a trustworthy clone for hierarchical clonotyping ONLY when hierarchical is set to "double.and.single.chains". Defaults to 3, meaning that, an exact combination of three chains needs to appear in the dataset at least 3 times for it to be considered as a clone, into which other cells are merged. (For the counting of exact combination of chains CDR3 nucleotide string matching is used, even if clonotyping by homology)
#' @param global.clonotype Logical specifying whether clonotyping should occur across samples or only within a single sample (grouping via sample_id column).
#' @param VDJ.VJ.1chain Logical specifying whether cells other than once with 1 VDJ and 1 VJ chains should be considered.
#' @param output.format Parameter output.format is deprecated. If non VGM-style output is required please refer to the function VDJ_clonotype. Output is VGM style VDJ by cell dataframe
#' @param platypus.version Only "v3" available
#' @return Returns a VGM[[1]]-type dataframe. All columns of the input VDJ dataframe are maintained and extra columns are added with clonotyping information.New columns are named by clonotyping strategy so to allow for multiple clonotyping identifiers to be present in the same VDJ dataframe and make comparisons between these straighforward.
#' @export
#' @examples
#' reclonotyped_vgm <- VDJ_clonotype_v3(VDJ=Platypus::small_vgm[[1]],
#' clone.strategy="cdr3.nt",
#' hierarchical = "none", global.clonotype = TRUE)
#'
#' reclonotyped_vgm <- VDJ_clonotype_v3(VDJ=Platypus::small_vgm[[1]],
#' clone.strategy="cdr3.homology", homology.threshold = 0.5,
#' hierarchical = "single.chains", global.clonotype = TRUE)
#'
VDJ_clonotype_v3 <- function(VDJ,
                             clone.strategy,
                             homology.threshold,
                             hierarchical,
                             triple.chain.count.threshold,
                             global.clonotype,
                             VDJ.VJ.1chain,
                             output.format,
                             platypus.version){

  #for CRAN checks
  Nr_of_VDJ_chains <- NULL
  Nr_of_VJ_chains <- NULL
  ccombs1 <- NULL
  unique_id_internal <- NULL
  sample_id <- NULL
  new_clonal_feature <- NULL


  #### Functions definition ####
  match_ab_single <- function(sample_dfs, sample_aberrant){ #for matching cells with only a single chain in total to cells with anything more than that. Needs a new_clonal_feature column to compare
    if(nrow(sample_aberrant) > 0){
      curr_ab <- sample_aberrant

      message(paste0("Attempting to merge in ", length(unique(sample_aberrant$unique_id_internal)), " aberrant cells"))

      for(cel in 1:nrow(curr_ab)){
        clone_matches <- which(stringr::str_detect(sample_dfs$new_clonal_feature, curr_ab$new_clonal_feature[cel]))
        if(sample_aberrant$Nr_of_VDJ_chains[cel] == 0){
        }
        if(length(unique(sample_dfs$new_clonal_feature[clone_matches])) > 1){ #This returns TRUE if multiple defined 1VDJ 1VJ clones match the pattern of the aberrant query clone

          #update: prioritisation of merging cells into clonotypes of the same sample
          same_sample <- subset(sample_dfs[clone_matches,], sample_id == curr_ab$sample_id[cel])
          #are there any existing clones in the same sample as the current aberrant cell which match this cells pattern?
          if(nrow(same_sample) > 0){ #-> yes
            #Assigning the aberrant query clone to the most frequent matching clone of the same sample
            curr_ab$new_clonal_feature[cel] <- names(which.max(table(same_sample$new_clonal_feature)))
          } else { #-> no. Proceed as before update
            #Assigning the aberrant query clone to the most frequent matching clone
            curr_ab$new_clonal_feature[cel] <- names(which.max(table(sample_dfs$new_clonal_feature[clone_matches])))
          }
        } else if(length(unique(sample_dfs$new_clonal_feature[clone_matches])) == 1){#This returns TRUE if exactly one predefined clone matches the pattern of the abberant query clone
          curr_ab$new_clonal_feature[cel] <- unique(sample_dfs$new_clonal_feature[clone_matches]) #Assigning the aberrant query clone to the only matching clone
        }  #ELSE: no clone found with the light chain of this cell => open a new clone by leaving the clonal feature as is
      }
      sample_aberrant <- curr_ab #re transpose into list
    }
    return(sample_aberrant)
  }

  match_ab_single_homology <- function(sample_dfs, sample_aberrant, clone.strategy, homology.threshold){ #for matching cells with only a single chain in total to cells with anything more than that. Needs a new_clonal_feature column to compare
    if(nrow(sample_aberrant) > 0){
      curr_ab <- sample_aberrant

      message(paste0("Attempting to merge in ", length(unique(sample_aberrant$unique_id_internal)), " aberrant cells"))

      for(cel in 1:nrow(curr_ab)){

        #check if the chain matches any already existing clone
        clone_matches <- which(stringr::str_detect(sample_dfs$new_clonal_feature, curr_ab$new_clonal_feature[cel]))
        #now check for homology
        if(length(clone_matches) > 0){
          #Get normalized string distances between the current aberrant cell and the clones with matching v gene usage and cdr3 length
          if(nchar(curr_ab$VDJ_cdr3s_aa[cel]) != 0){
            dists_VDJ <- stringdist::stringdist(sample_dfs$VDJ_cdr3s_aa[clone_matches], curr_ab$VDJ_cdr3s_aa[cel]) / nchar(curr_ab$VDJ_cdr3s_nt[cel])
          } else if(nchar(curr_ab$VDJ_cdr3s_aa[cel]) == 0){
            if(clone.strategy %in% c("Hvj.Lvj.CDR3length.CDR3.homology", "CDR3.homology")){
              dists_VDJ <- rep(0, length(clone_matches))
            } else if(clone.strategy %in% c("Hvj.Lvj.CDR3length.CDRH3.homology", "CDRH3.homology")){ #problem case: looking for VDJ homology only, but cell has no VDJ chain. Will be its own clonotype
              dists_VDJ <- rep(100, length(clone_matches))
            }
          }
          if(nchar(curr_ab$VJ_cdr3s_nt[cel]) != 0 & clone.strategy %in% c("Hvj.Lvj.CDR3length.CDR3.homology", "CDR3.homology")){
            dists_VJ <- stringdist::stringdist(sample_dfs$VJ_cdr3s_aa[clone_matches], curr_ab$VJ_cdr3s_aa[cel]) / nchar(curr_ab$VJ_cdr3s_aa[cel])
          } else {dists_VJ <- rep(0, length(clone_matches))}
          dists <- dists_VDJ + dists_VJ #getting combined distances and testing vs. threshold

          #Directo  homology criterion
          if(any(dists <= homology.threshold)){
            curr_ab$new_clonal_feature[cel] <- sample_dfs$new_clonal_feature[clone_matches][which.min(dists)] #Assigning the aberrant query clone to the only matching clone
            #ELSE: no clone found with the light chain of this cell => open a new clone
          } else { curr_ab$new_clonal_feature[cel] <- paste0("hclust_nomatch_", curr_ab$unique_id_internal)}
        }
      }
      sample_aberrant <- curr_ab #re transpose into list
    }
    return(sample_aberrant)
  }
  #### Parameter setup ####
  platypus.version <- "v3"
  if(missing(global.clonotype)) global.clonotype <- F
  if(missing(clone.strategy)) clone.strategy <- "cdr3.aa"
  if(missing(VDJ.VJ.1chain)) VDJ.VJ.1chain <- T
  if(missing(VDJ)) stop("Please provide input data as VDJ")
  if(missing(hierarchical)) hierarchical <- "none"
  if(missing(homology.threshold)) homology.threshold <- 0.3
  if(missing(triple.chain.count.threshold)) triple.chain.count.threshold <- 5

  if(!missing(output.format)){
    message("Parameter output.format is deprecated. Standard output is VGM style VDJ by cell dataframe")
  }
  if(hierarchical == F){
    hierarchical <- "none"
    message("After function updates, please set hierachical to either 'none', 'single.chains', 'double.and.single.chains'. hierarchical was set to 'none' based on your current input for backwards compatibility")
  }
  if(hierarchical == T){
    hierarchical <- "single.chains"
    message("After function updates, please set hierachical to either 'none', 'single.chains', 'double.and.single.chains'. hierarchical was set to 'single.chains' based on your current input for backwards compatibility")
  }
  if(!hierarchical %in% c("none", "single.chains", "double.and.single.chains")){
    stop("Please set hierachical to either 'none', 'single.chains', 'double.and.single.chains'.")
  }

  #Making cloning stategy fitting with VDJ / VJ naming scheme
  #This way the old keyworks will still work and this update should not break any old code
  #remember for renaming later
  clone.strategy.as.input <- clone.strategy
  switch(clone.strategy,
         VDJJ.VJJ.cdr3length.cdr3.homology = {clone.strategy <- 'Hvj.Lvj.CDR3length.CDR3.homology'},
         VDJJ.VJJ.cdr3length.VDJcdr3.homology  = {clone.strategy <- 'Hvj.Lvj.CDR3length.CDRH3.homology'},
         cdr3.homology = {clone.strategy <- 'CDR3.homology'},
         VDJcdr3.homology = {clone.strategy <- 'CDRH3.homology'},
         VDJJ.VJJ = {clone.strategy <- "hvj.lvj"},
         VDJJ.VJJ.cdr3length = {clone.strategy <- "hvj.lvj.cdr3lengths"})

  if(!clone.strategy %in% c("10x.default", "cdr3.nt", "cdr3.aa", "hvj.lvj", "hvj.lvj.cdr3lengths", "Hvj.Lvj.CDR3length.CDR3.homology", "Hvj.Lvj.CDR3length.CDRH3.homology", "CDR3.homology", "CDRH3.homology")){
    stop("Please enter one of the following clonotyping stategies: \n  10x.default, cdr3.nt, cdr3.aa, VDJJ.VJJ, VDJJ.VJJ.cdr3length, VDJJ.VJJ.cdr3length.cdr3.homology, VDJJ.VJJ.cdr3length.VDJcdr3.homology, cdr3.homology, VDJcdr3.homology")
  } else {
    message(paste0("Clonotyping strategy: ", clone.strategy.as.input))
  }
  #### prep for hierarchical ####
  #adding a unique identifier column for later
  VDJ$unique_id_internal <- c(paste0("id", 1:nrow(VDJ)))
  VDJ$new_clonal_feature <- NA
  #adding a column for nchars and setting to empty strings if length 0 (makes pasting and comparing strings later easier)
  VDJ$nchar_VDJ_cdr3s_aa <- nchar(VDJ$VDJ_cdr3s_aa)
  VDJ$nchar_VDJ_cdr3s_aa[VDJ$nchar_VDJ_cdr3s_aa == 0] <- ""
  VDJ$nchar_VJ_cdr3s_aa <- nchar(VDJ$VJ_cdr3s_aa)
  VDJ$nchar_VJ_cdr3s_aa[VDJ$nchar_VJ_cdr3s_aa == 0] <- ""

  #run through of the columns so to make them all characters and empty strings if neccessary
  for(i in 1:ncol(VDJ)){
    if(names(VDJ)[i] %in% c("VDJ_cdr3s_aa", "VJ_cdr3s_aa", "VDJ_cdr3s_nt", "VJ_cdr3s_nt", "VDJ_vgene", "VDJ_jgene", "VJ_vgene", "VJ_jgene")){
      VDJ[,i] <- as.character(VDJ[,i])
      VDJ[is.na(VDJ[,i]),i] <- ""
      VDJ[VDJ[,i] == "NA",i] <- ""
    }
  }

  if(!"sample_id" %in% names(VDJ) & global.clonotype == F){
    stop("sample_id column needed for clonotyping by sample was not found in the input dataframe")
  }
  if(!"sample_id" %in% names(VDJ) & hierarchical != "none" ){
    warning("sample_id column not found in input dataframe. This is used to merge single chain cells into existing clonotypes of the same sample with a higher priority than other samples. This feature is not available without a valid sample_id column")
    VDJ$sample_id <- "no_info_provided"
  }

  #redefining count so to make sure
  VDJ$Nr_of_VDJ_chains <- 0
  VDJ$Nr_of_VDJ_chains[nchar(VDJ$VDJ_cdr3s_aa) > 0] <- 1
  VDJ$Nr_of_VDJ_chains[stringr::str_detect(VDJ$VDJ_cdr3s_aa, ";")] <- 2
  VDJ$Nr_of_VJ_chains <- 0
  VDJ$Nr_of_VJ_chains[nchar(VDJ$VJ_cdr3s_aa) > 0] <- 1
  VDJ$Nr_of_VJ_chains[stringr::str_detect(VDJ$VJ_cdr3s_aa, ";")] <- 2

  if(hierarchical == "none"){ #NOT hierarchical
    #only include clones with one heavy and one light chain
    if(VDJ.VJ.1chain){
      VDJ<- VDJ[which(VDJ$Nr_of_VDJ_chains==1 & VDJ$Nr_of_VJ_chains==1),]
    }####STOP strict
  } else { # hierarchical (either all or only single chain)
    #only include clones with one heavy and one light chain
    if(VDJ.VJ.1chain){
      message("Hierarchical clonotyping is specifically designed to incorporate cells with abberand numbers of chains. Filtering for 1VDJ 1VJ chain thereby defeats its purpose. Function will continue without filtering.")
    }
    prior_filtering <- nrow(VDJ)
    #filtering cells with 2 VDJ and 2 VJ chains as but not cells with e.g. 2 VJ and 1 VDJ chain. Cells with only 1 chain in total are also kept
    VDJ <- subset(VDJ,  (VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains == 0) | (VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains == 1) | (VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains == 2) | (VDJ$Nr_of_VDJ_chains == 0 & VDJ$Nr_of_VJ_chains == 2) | (VDJ$Nr_of_VDJ_chains == 0 & VDJ$Nr_of_VJ_chains == 1))
    if(nrow(VDJ) > 0){
      message(paste0("Filtered out ", prior_filtering - nrow(VDJ), " cells containing more than one VDJ AND VJ chain or two VDJ chains, as these likely correspond to doublets"))
    } else {stop("After filtering for doublets, no cells remained")}

    #here we pick only single chain cells as aberrant cells. All others are treated as if hierarchical was set to none
    aberrant <-  subset(VDJ, (VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains == 0) | (VDJ$Nr_of_VDJ_chains == 0 & VDJ$Nr_of_VJ_chains == 1))

    if(nrow(aberrant) > 0){ #check if there are aberrant cells

      #removing these cells by the added unique id from the main table
      VDJ <- VDJ[!VDJ$unique_id_internal %in% aberrant$unique_id_internal, ]

      if(hierarchical == "double.and.single.chains"){
        #here we pick out two sets: single chain cells and cells with double chains of either VDJ or VJ chains (including those with 3 chains in total)
        double_aberrant <- subset(VDJ, (VDJ$Nr_of_VDJ_chains > 1 & VDJ$Nr_of_VJ_chains == 0) | (VDJ$Nr_of_VDJ_chains > 1 & VDJ$Nr_of_VJ_chains == 1) | (VDJ$Nr_of_VDJ_chains == 1 & VDJ$Nr_of_VJ_chains > 1) | (VDJ$Nr_of_VDJ_chains == 0 & VDJ$Nr_of_VJ_chains > 1))

        VDJ <- VDJ[!VDJ$unique_id_internal %in% double_aberrant$unique_id_internal, ]

        if(nrow(double_aberrant) > 0){ #checking that we have cells

          #to be able to remap original cells later we keep a copy before making pseudocells
          backup_double_aberrant <- double_aberrant

          #### Checking minimal count to use triple aberrant cells as hub cells for clonotypes ####
          triple_aberrant <- subset(double_aberrant, (double_aberrant$Nr_of_VDJ_chains > 1 & double_aberrant$Nr_of_VJ_chains == 1) | (double_aberrant$Nr_of_VDJ_chains == 1 & double_aberrant$Nr_of_VJ_chains > 1))

          #to be able to remap original cells later we keep a copy before making pseudocells
          backup_triple_aberrant <- triple_aberrant

          #order CDR3s alphabetically
          for(z in 1:nrow(triple_aberrant)){
            triple_aberrant$VDJ_cdr3s_nt_check[z] <- paste0(sort(stringr::str_split(triple_aberrant$VDJ_cdr3s_nt[z], ";", simplify = T)[1,]), collapse = "")
            triple_aberrant$VJ_cdr3s_nt_check[z] <- paste0(sort(stringr::str_split(triple_aberrant$VJ_cdr3s_nt[z], ";", simplify = T)[1,]), collapse = "")
          }
          #add new checking feature
          triple_aberrant$thresh_check_feature <- paste0(triple_aberrant$VDJ_cdr3s_nt_check, triple_aberrant$VJ_cdr3s_nt_check)
          #getting those cells which are higher in frequency then the set threshold
          high_conf_unique_features <- names(table(triple_aberrant$thresh_check_feature)[table(triple_aberrant$thresh_check_feature) > triple.chain.count.threshold])


          message(paste0("Found ", length(high_conf_unique_features), " exact matching clones with 3 chains and a frequency of at least ", triple.chain.count.threshold,". These will be used as high confidence clonotypes."))

          #keep track of the ID of these cells
          high_conf_triple_chains <- triple_aberrant$unique_id_internal[triple_aberrant$thresh_check_feature %in% high_conf_unique_features]

          #not needed anymore
          triple_aberrant <- NULL

          #removing all double and triple aberant cells by the added unique id from the main table
          VDJ <- VDJ[!VDJ$unique_id_internal %in% double_aberrant$unique_id_internal, ]

          #now the key part +> use the unique identifer to split up double chain cells in to 2 normal cells each and pass them trough as normal
          for(i in 1:ncol(double_aberrant)){
            if(names(double_aberrant)[i] %in% c("VDJ_cdr3s_aa", "VJ_cdr3s_aa", "VDJ_cdr3s_nt", "VJ_cdr3s_nt", "VDJ_vgene", "VDJ_jgene", "VJ_vgene", "VJ_jgene")){
              double_aberrant[!stringr::str_detect(double_aberrant[,i], ";"),i] <- paste0(double_aberrant[!stringr::str_detect(double_aberrant[,i], ";"),i], ";")
            }
          }
          #making pseudocells
          pseudo_cells <- list()
          for(i in 1:nrow(double_aberrant)){
            curr <- rbind(double_aberrant[i,], double_aberrant[i,]) #make dataframe with two identical rows
            for(j in 1:ncol(curr)){
              if(names(curr)[j] %in% c("VDJ_cdr3s_aa", "VJ_cdr3s_aa", "VDJ_cdr3s_nt", "VJ_cdr3s_nt", "VDJ_vgene", "VDJ_jgene", "VJ_vgene", "VJ_jgene")){
                #add a ; to chains that do not have one to to make splitting easier
                #remove ;from the end
                curr[,j] <- gsub(";$", "", curr[,j])
                curr[!stringr::str_detect(curr[,j], ";"),j] <- paste0(curr[!stringr::str_detect(curr[,j], ";"),j], ";",curr[!stringr::str_detect(curr[,j], ";"),j])
                #For cells in row 1 = make those equal to the first "chain" of the aberrant cell
                curr[1,j] <- stringr::str_split(curr[1,j], ";", simplify = T)[1,1]
                #For cells in row 2 = make those equal to the second "chain" of the aberrant cell
                curr[2,j] <- stringr::str_split(curr[2,j], ";", simplify = T)[1,2]
              }
            }
            #updating nchar values
            curr$nchar_VDJ_cdr3s_aa <- nchar(curr$VDJ_cdr3s_aa)
            curr$nchar_VDJ_cdr3s_aa[curr$nchar_VDJ_cdr3s_aa == 0] <- ""
            curr$nchar_VJ_cdr3s_aa <- nchar(curr$VJ_cdr3s_aa)
            curr$nchar_VJ_cdr3s_aa[curr$nchar_VJ_cdr3s_aa == 0] <- ""
            #adding to list
            pseudo_cells[[i]] <- curr
          }
          #coerce into one dataframe
          pseudo_cells <- dplyr::bind_rows(pseudo_cells)
          #now to distribute these cells into either the normal VDJ or the aberrant cells

          normal_pseudo <- subset(pseudo_cells, pseudo_cells$Nr_of_VDJ_chains != 0 & pseudo_cells$Nr_of_VJ_chains != 0)
          aberrant_pseudo <- subset(pseudo_cells, pseudo_cells$Nr_of_VDJ_chains == 0 | pseudo_cells$Nr_of_VJ_chains == 0)

          #bind with corresponding dataframes
          VDJ <- rbind(VDJ, normal_pseudo)
          aberrant <- subset(aberrant, !unique_id_internal %in% aberrant_pseudo$unique_id_internal)
          aberrant <- rbind(aberrant, aberrant_pseudo)

        } else {#end checking if we have double aberrant cells
          clone.strategy <- "single.chains" #if there are not double aberrant chains, we switch back to normal hierarchical clonotyping
        }

      }
    } else {#end checking if we have  aberrant cells
      clone.strategy <- "none"
    }
  }
  #### prep for global clonotying ####
  if(global.clonotype==F){ # loop through each repertoire individually
    repertoire.number <- unique(as.character(VDJ$sample_id))
    sample_dfs <- list()
    sample_aberrant <- list()
    for(j in 1:length(repertoire.number)){ #Add samples to list
      sample_dfs[[j]] <- subset(VDJ, sample_id==repertoire.number[j])
      if(hierarchical != "none"){
        #subsetting the aberrant cells two. The chance is quite high that in small repertoires some of these resulting dataframes will contain 0 rows.
        sample_aberrant[[j]] <- subset(aberrant, sample_id==repertoire.number[j])
      }
    }
  } else { #global clonotyping. Convert to list so to be able to run it through the same loop as for non global clonotyping
    sample_dfs <- list(VDJ)
    if(hierarchical != "none"){
      sample_aberrant <- list(aberrant)
    }
  }
  #### Start sample loop ####
  for(i in 1:length(sample_dfs)){
    ####Clonotyping strategies

    message(paste0("Processing sample ", i, " of ", length(sample_dfs)))

    if(clone.strategy=="10x.default"){ ####START 10x default
      sample_dfs[[i]]$new_clonal_feature <- sample_dfs[[i]]$clonotype_id_10x
      if(hierarchical != "none"){ ####START Hierarchical 10x default
        sample_aberrant[[i]] <- lapply(sample_aberrant[[i]], function(x){
          x$new_clonal_feature <- x$clonotype_id_10x
          return(x)})
      } ####STOP Hierarchical 10x default
    } #### STOP default 10x
    else if(clone.strategy=="cdr3.nt"){ ####START cdr3.nt

      sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_cdr3s_nt,
                                                  sample_dfs[[i]]$VJ_cdr3s_nt, sep="")

      if(hierarchical != "none"){ ####START Hierarchical single cdr3.nt

        #define clonal features
        sample_aberrant[[i]]$new_clonal_feature <- paste(sample_aberrant[[i]]$VDJ_cdr3s_nt,
                                                         sample_aberrant[[i]]$VJ_cdr3s_nt, sep="")

        #match to main dataframe and update clonal features accordingly
        sample_aberrant[[i]] <- match_ab_single(sample_dfs[[i]], sample_aberrant[[i]])
      } ####STOP Hierarchical cdr3.nt
    } ####STOP cdr3.nt
    else if(clone.strategy=="cdr3.aa"){ ####START cdr3.aa

      sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_cdr3s_aa,
                                                  sample_dfs[[i]]$VJ_cdr3s_aa, sep="")

      if(hierarchical != "none"){ ####START Hierarchical single cdr3.aa

        #define clonal features
        sample_aberrant[[i]]$new_clonal_feature <- paste(sample_aberrant[[i]]$VDJ_cdr3s_aa,
                                                         sample_aberrant[[i]]$VJ_cdr3s_aa, sep="")

        #match to main dataframe and update clonal features accordingly
        sample_aberrant[[i]] <- match_ab_single(sample_dfs[[i]], sample_aberrant[[i]])
      } ####STOP Hierarchical cdr3.aa
    } ####STOP cdr3.aa
    else if(clone.strategy=="hvj.lvj"){ ####START hvj.lvj
      sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                  sample_dfs[[i]]$VDJ_jgene,
                                                  sample_dfs[[i]]$VJ_vgene,
                                                  sample_dfs[[i]]$VJ_jgene,sep="")

      if(hierarchical != "none"){ ####START Hierarchical single cdr3.nt

        sample_aberrant[[i]]$new_clonal_feature <- paste(sample_aberrant[[i]]$VDJ_vgene,
                                                         sample_aberrant[[i]]$VDJ_jgene,
                                                         sample_aberrant[[i]]$VJ_vgene,
                                                         sample_aberrant[[i]]$VJ_jgene,sep="")

        #match to main dataframe and update clonal features accordingly
        sample_aberrant[[i]] <- match_ab_single(sample_dfs[[i]], sample_aberrant[[i]])
      } ####STOP Hierarchical hvj.lvj
    } ####STOP hvj.lvj
    else if(clone.strategy=="hvj.lvj.cdr3lengths"){ ####START hvj.lvj.cdr3lengths
      sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                  sample_dfs[[i]]$VDJ_jgene,
                                                  sample_dfs[[i]]$nchar_VDJ_cdr3s_aa,
                                                  sample_dfs[[i]]$VJ_vgene,
                                                  sample_dfs[[i]]$VJ_jgene,
                                                  sample_dfs[[i]]$nchar_VJ_cdr3s_aa,sep="")

      if(hierarchical != "none"){ ####START Hierarchical single hvj.lvj.cdr3lengths

        sample_aberrant[[i]]$new_clonal_feature <- paste(sample_aberrant[[i]]$VDJ_vgene,
                                                         sample_aberrant[[i]]$VDJ_jgene,
                                                         sample_aberrant[[i]]$nchar_VDJ_cdr3s_aa,
                                                         sample_aberrant[[i]]$VJ_jgene,
                                                         sample_aberrant[[i]]$VJ_jgene,
                                                         sample_aberrant[[i]]$nchar_VJ_cdr3s_aa,sep="")

        sample_aberrant[[i]] <- match_ab_single(sample_dfs[[i]], sample_aberrant[[i]]) #match to main dataframe and update clonal features accordingly
      } #### Stop hierarchical single hvj.lvj.cdr3lengths
    } ####STOP hvj.lvj.cdr3lengths
    else if(clone.strategy=="Hvj.Lvj.CDR3length.CDR3.homology" | clone.strategy=="Hvj.Lvj.CDR3length.CDRH3.homology"){   ####START hvj.lvj.cdr3lengths homology

      sample_dfs[[i]]$new_clonal_feature <- paste(sample_dfs[[i]]$VDJ_vgene,
                                                  sample_dfs[[i]]$VDJ_jgene,
                                                  sample_dfs[[i]]$nchar_VDJ_cdr3s_aa,
                                                  sample_dfs[[i]]$VJ_vgene,
                                                  sample_dfs[[i]]$VJ_jgene,
                                                  sample_dfs[[i]]$nchar_VJ_cdr3s_aa,sep="")
      #loop over unique new clonal features
      curr_cl <- list()
      for(k in 1:length(unique(sample_dfs[[i]]$new_clonal_feature))){
        curr_cl[[k]] <- subset(sample_dfs[[i]], new_clonal_feature == unique(sample_dfs[[i]]$new_clonal_feature)[k])
        if(nrow(curr_cl[[k]]) > 1){ #only passing ones to hclust which have at least two entries
          nchars_VDJ <- nchar(curr_cl[[k]]$VDJ_cdr3s_aa)
          nchars_VDJ[nchars_VDJ == 0] <- 1 #account for missing chains
          #open stringdist matrix
          VDJ_distance <- stringdist::stringdistmatrix(curr_cl[[k]]$VDJ_cdr3s_aa, curr_cl[[k]]$VDJ_cdr3s_aa,method = "lv")/nchars_VDJ
          if (clone.strategy=="Hvj.Lvj.CDR3length.CDR3.homology"){
            nchars_VJ <- nchar(curr_cl[[k]]$VJ_cdr3s_aa)
            nchars_VJ[nchars_VJ == 0] <- 1 #account for missing chains
            #open stringdist matrix
            VJ_distance <- stringdist::stringdistmatrix(curr_cl[[k]]$VJ_cdr3s_aa,curr_cl[[k]]$VJ_cdr3s_aa,method = "lv")/nchars_VJ
          } else {VJ_distance <- 0} #in case of Hvj.Lvj.CDR3length.CDRH3.homology clonotyping
          #draw clusters
          combined_distance <- VDJ_distance + VJ_distance
          diag(combined_distance) <- NA
          hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
          hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
          #add to clonal features
          curr_cl[[k]]$new_clonal_feature <- paste(curr_cl[[k]]$new_clonal_feature,"_",k, "_",hclust_combined_cut)
        }
      }
      sample_dfs[[i]] <- dplyr::bind_rows(curr_cl) #reassemble the dataframe

      if(hierarchical != "none"){ ####START Hierarchical single hvj.lvj.cdr3lengths

        sample_aberrant[[i]]$new_clonal_feature <- paste(sample_aberrant[[i]]$VDJ_vgene,
                                                         sample_aberrant[[i]]$VDJ_jgene,
                                                         sample_aberrant[[i]]$nchar_VDJ_cdr3s_nt,
                                                         sample_aberrant[[i]]$VJ_vgene,
                                                         sample_aberrant[[i]]$VJ_jgene,
                                                         sample_aberrant[[i]]$nchar_VDJ_cdr3s_nt,sep="")

        sample_aberrant[[i]] <- match_ab_single_homology(sample_dfs[[i]], sample_aberrant[[i]], clone.strategy, homology.threshold) #match to main dataframe and update clonal features accordingly
      } #### STOP hierarchical single hvj.lvj.cdr3lengths
    } ####STOP hvj.lvj.cdr3lengths homology
    else if (clone.strategy=="CDR3.homology" | clone.strategy=="CDRH3.homology"){####START homology

      if(nrow(sample_dfs[[i]]) > 1){ #only passing ones to hclust which have at least two entries
        nchars_VDJ <- nchar(sample_dfs[[i]]$VDJ_cdr3s_aa)
        nchars_VDJ[nchars_VDJ == 0| is.na(nchars_VDJ)] <- 1 #account for missing chains
        cdr3_VDJ <- sample_dfs[[i]]$VDJ_cdr3s_aa

        #open stringdist matrix
        VDJ_distance <- stringdist::stringdistmatrix(cdr3_VDJ, cdr3_VDJ,method = "lv")/nchars_VDJ
        if (clone.strategy=="CDR3.homology"){
          nchars_VJ <- nchar(sample_dfs[[i]]$VJ_cdr3s_aa)
          nchars_VJ[nchars_VJ == 0 | is.na(nchars_VJ)] <- 1 #account for missing chains
          cdr3_VJ <- sample_dfs[[i]]$VJ_cdr3s_aa

          #open stringdist matrix
          VJ_distance <- stringdist::stringdistmatrix(cdr3_VJ,cdr3_VJ,method = "lv")/nchars_VJ
        } else {VJ_distance <- 0} #in case of Hvj.Lvj.CDR3length.CDRH3.homology clonotyping
        #draw clusters
        combined_distance <- VDJ_distance + VJ_distance
        diag(combined_distance) <- NA

        hclust_combined <- stats::hclust(stats::as.dist(combined_distance)) #convert combined_distance to a distance object
        hclust_combined_cut <- stats::cutree(hclust_combined, h = homology.threshold)
        #add to clonal features

        sample_dfs[[i]]$new_clonal_feature <- paste0("hclust", hclust_combined_cut)
      }
      if(hierarchical != "none"){ ####START Hierarchical single homology

        sample_aberrant[[i]]$new_clonal_feature <- "hclust" #The match_ab_single_homology function runs a str_detect to find clonal_features in the main dataframe which contain the clonal feature...
        #of the aberrant cell in question. Here we set all clonal features of the aberrant cells to hclust. This string is also contained in all clonal features of the main table...
        #(see above). Therefore, for each aberrant cell all cells from the main dataframe are taken into consideration and cells are merged in only based on homology
        sample_aberrant[[i]] <- match_ab_single_homology(sample_dfs[[i]], sample_aberrant[[i]], clone.strategy, homology.threshold) #match to main dataframe and update clonal features accordingly
      }####STOP Hierarchical single homology
    }####STOP homology

    #### regrouping ####
    #grouping and merging again to get accurate frequencies
    if(hierarchical != "none"){
      aggr_aberrant <- sample_aberrant[[i]]
      if(any(names(aggr_aberrant) != names(sample_dfs[[i]]))){
        stop("Not able to join main and merging dataframes due to column mismatch")
      } else {
        sample_dfs[[i]] <- rbind(sample_dfs[[i]], aggr_aberrant)
      }
      #recondensing pseudocells into original ones and matching clonotypes accordingly

      if(hierarchical == "double.and.single.chains"){

        #input are the dataframe of normal cells and the dataframe of aberrant cells. Both contain pseudocells
        #we first have to pick out these pseudocells by getting only the rows which have duplicated unique internal ids
        pseudo_cells <- sample_dfs[[i]][sample_dfs[[i]]$unique_id_internal %in% names(table(sample_dfs[[i]]$unique_id_internal)[table(sample_dfs[[i]]$unique_id_internal) > 1]), ]

        #now remove these from the sample_dfs[[i]] dataframe
        sample_dfs[[i]] <- subset(sample_dfs[[i]], !sample_dfs[[i]]$unique_id_internal %in% pseudo_cells$unique_id_internal)

        if(length(unique(pseudo_cells$unique_id_internal)) > 1){
          original_aberrant <- list()
          unique_ids_covered <- c()
          for(j in 1:length(unique(pseudo_cells$unique_id_internal))){
            if(!unique(pseudo_cells$unique_id_internal)[j] %in% unique_ids_covered){ #check if this id has already been processed as part of a clonotype
              curr_pseudo_cell <- subset(pseudo_cells, unique_id_internal == unique(pseudo_cells$unique_id_internal)[j])
              #getting the associated clonotypes for each cell
              curr_clonotypes <- curr_pseudo_cell$new_clonal_feature #this vector will always be length 2
              if(unique(pseudo_cells$unique_id_internal)[j] %in% high_conf_triple_chains){ #checking if this was passed the triple.chain.count.threshold
                #create an identifer to raise awarness that this is a joined clonotype
                update_clonotype_feature <- paste0("joint_", curr_clonotypes[1], "_", curr_clonotypes[2])

                #adding the new clonotype feature to ALL cells that matched at least one of the two pseudocells of that clonotype
                matches <- which(sample_dfs[[i]]$new_clonal_feature %in% curr_pseudo_cell$new_clonal_feature)

                if(length(matches) > 0){
                  sample_dfs[[i]]$new_clonal_feature[matches] <- update_clonotype_feature
                }

                #add a new clonal feature to the original cell. This will make it distinct from any other cells
                #now find the original cells that made up the pseudocells with this clonal feature and add the updated clonal feature
                all_matching_orig_ids <- pseudo_cells$unique_id_internal[pseudo_cells$new_clonal_feature %in% curr_clonotypes]
                #additional filtering step required: we only keep those original cells of which we have two rows aka BOTH pseudocells in that table. This way we only group exact matches but not extra pseudocells which may match this clonotype, but of which the original cell contained a non matching chain
                full_matching_orig_ids <- all_matching_orig_ids[duplicated(all_matching_orig_ids) == T]

                original_aberrant[[j]] <- subset(backup_double_aberrant, unique_id_internal %in% full_matching_orig_ids)
                original_aberrant[[j]]$new_clonal_feature <- update_clonotype_feature
                #keep track so we dont have to repeat the loop
                unique_ids_covered <- c( unique_ids_covered,original_aberrant[[j]]$unique_id_internal)

              } else { # did not pass the triple.chain.count.threshold

                #create an identifer to raise awarness that this is a low conf clonotype of only aberrant cells
                update_clonotype_feature <- paste0("double_cell_", "_matching_", curr_clonotypes[1], "_", curr_clonotypes[2])

                #add a new clonal feature to the original cell. This will make it distinct from any other cells
                #now find the original cells that made up the pseudocells with this clonal feature and add the updated clonal feature
                all_matching_orig_ids <- pseudo_cells$unique_id_internal[pseudo_cells$new_clonal_feature %in% curr_clonotypes]
                #additional filtering step required: we only keep those original cells of which we have two rows aka BOTH pseudocells in that table. This way we only group exact matches but not extra pseudocells which may match this clonotype, but of which the original cell contained a non matching chain
                full_matching_orig_ids <- all_matching_orig_ids[duplicated(all_matching_orig_ids) == T]

                original_aberrant[[j]] <- subset(backup_double_aberrant, unique_id_internal %in% full_matching_orig_ids)
                original_aberrant[[j]]$new_clonal_feature <- update_clonotype_feature

                #remove this from the sample_df dataframe
                sample_dfs[[i]]  <- subset(sample_dfs[[i]], !unique_id_internal %in% unique(pseudo_cells$unique_id_internal[pseudo_cells$new_clonal_feature %in% curr_clonotypes]))

                #keep track so we dont have to repeat the loop
                unique_ids_covered <- c(unique_ids_covered, original_aberrant[[j]]$unique_id_internal)
              }
            }
          } #end loop over pseudo cells clonal features
          #finally bind the sample_dfs to the aberrant cells
          sample_dfs[[i]] <- rbind(sample_dfs[[i]], dplyr::bind_rows(original_aberrant))
        }
      }
    }

    clones_g <- as.data.frame(sample_dfs[[i]] %>% dplyr::group_by(new_clonal_feature) %>% dplyr::summarise(new_clonal_frequency = dplyr::n()))
    clones_g <- clones_g[order(clones_g$new_clonal_frequency, decreasing = T),]
    clones_g$new_clonotype_id <- paste0("clonotype", 1:nrow(clones_g))
    sample_dfs[[i]] <- merge(sample_dfs[[i]], clones_g, by = "new_clonal_feature")
    #Order by frequency
    sample_dfs[[i]] <-sample_dfs[[i]][with(sample_dfs[[i]], order(-new_clonal_frequency)), ]
  }####STOP sample loop

  #### Output ####
  VDJ.GEX.matrix <- dplyr::bind_rows(sample_dfs)

  #cleanup
  VDJ.GEX.matrix <- VDJ.GEX.matrix[,-c(which(names(VDJ.GEX.matrix) %in% c("nchar_VDJ_cdr3s_aa", "nchar_VJ_cdr3s_aa","unique_id_internal")))]
  VDJ.GEX.matrix <- VDJ.GEX.matrix[,c(2:ncol(VDJ.GEX.matrix), 1)]

  if(!global.clonotype){
    if(any(stringr::str_detect(names(VDJ.GEX.matrix), paste0("clonotype_id_",clone.strategy.as.input)))){
      message("Removing columns from VDJ.GEX.matrix stemming from a previous clonotyping run with same setting")
      VDJ.GEX.matrix <- VDJ.GEX.matrix[,!names(VDJ.GEX.matrix) %in% c(paste0("clonotype_id_",clone.strategy.as.input), paste0("clonal_feature_",clone.strategy.as.input), paste0("clonotype_frequency_",clone.strategy.as.input))]
    }
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonotype_id")] <- paste0("clonotype_id_",clone.strategy.as.input)
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_feature")] <- paste0("clonal_feature_",clone.strategy.as.input)
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_frequency")] <- paste0("clonotype_frequency_",clone.strategy.as.input)
  } else {
    if(any(stringr::str_detect(names(VDJ.GEX.matrix), paste0("global_clonotype_id_",clone.strategy.as.input)))){
      VDJ.GEX.matrix <- VDJ.GEX.matrix[,!names(VDJ.GEX.matrix) %in% c(paste0("global_clonotype_id_",clone.strategy.as.input), paste0("global_clonal_feature_",clone.strategy.as.input), paste0("global_clonotype_frequency_",clone.strategy.as.input))]
      message("Removing columns from VDJ.GEX.matrix stemming from a previous clonotyping run with same setting")
    }
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonotype_id")] <- paste0("global_clonotype_id_",clone.strategy.as.input)
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_feature")] <- paste0("global_clonal_feature_",clone.strategy.as.input)
    names(VDJ.GEX.matrix)[which(names(VDJ.GEX.matrix) == "new_clonal_frequency")] <- paste0("global_clonotype_frequency_",clone.strategy.as.input)
  }
  return(VDJ.GEX.matrix)
}
