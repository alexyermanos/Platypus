#find_public_cells function:
#'Function searches for public cells and clones in vgm outputs created by the VDJ_GEX_matrix() function from platypus. The definition of a public cell can be provided in the input: define.public
#'
#'In this function cdr3 sequences from a primary dataset (input: vgm) are compared to cdr3 sequences in a second (public) dataset (input: vdj_public).
#'
#'For example if define.public = "VDJ_cdr3s_aa" is used (default), any cell from the primary dataset with a cdr3 heavy chain amino acid sequence matching any public cdrh3 aa sequence is marked as a public cell.
#'Also, all other cells which belong to the same clone (as the public cell) are marked as public clones.
#'
#'This function returns a list containing in it's first element the new vgm with the added public status information and the second element containing a data frame of all clonotypes and their public status.
#'
#'
#'This function can be used to simply find public cells and clones within a dataset on the vdj (vmg[[1]]) level, or it can also integrate the "public" information into the gex (vgm[[2]]).
#'Integration into vgm[[2]] should only be done if gene expression of the public cells or public clones is analyzed afterwards, due to much longer computational time.
#'
#'
#'Definitions for define.public = "VDJ_cdr3s_aa":
#'If a cell in the first dataset and the public data (2nd dataset) have the exact same cdr3 heavy chain aa sequence, this cell is defined as a "public cell".
#'If any cell within a clonotype (defined by 10x default) is a public cell, this clonotype is defined as a "public clone". This information is stored additionally in the second list element of the output.
#'
#'
#'The public status information is stored in 6-7 newly created columns:
#'1. "overlap_type" - Character, indicating whether the cell has any CDR3 amino acid overlaps with the public data. Only created if check.overlap.type = TRUE
#'2. "public_definition" - Character, stores the information which was used to define a cell as public (e.g. the CDRH3 sequence)
#'3. "public" - Logical, if TRUE the corresponding cell in the df is defined as a public cell
#'4. "counts" - Integer, counts the number of cells in the "public" data
#'5. "nr_unique_sample_ids" - Integer, indicates in how many different samples from the public data cdrh3 matches are occurring. Indicates how "public" this cdr3h sequence really is. (In how many different mice or organs of the public data is a cell with the exact same cdrh3 aa seq. present)
#'6. "matches_barcode" - List, contains  barcodes of all cells within the public data which have the same cdrh3 sequence
#'7."public_clonotype" - Logical, if TRUE this cell belongs to a clonotype which has at least one public cell
#'
#'
#'@param vgm output of VDJ_GEX_matrix function in platypus, first dataset (public cell/clone information will only be stored in vgm_1)
#'@param vdj_public vgm[[1]] of all public data used for comparison joined together (multiple datasets can be joined, no information is stored here)
#'@param define.public Character, defines which column(s) in the vgm[[1]] are used to define a public cell. Default is set to "VDJ_cdr3s_aa". Other options are: "VJ_cdr3s_aa", "VJ_cdr3s_aa.VDJ_cdr3s_aa", "VDJ_cdr3s_nt"
#'@param check.overlap.type Logical, defines if overlap type information should be stored in the new vgm. Default is set to FALSE, setting this to TRUE will increase computational time. (Information on cd3 overlaps on the aa level will be stored)
#'@param integrate.into.gex Logical, defines if the public cell/clone information should be integrated into the vgm[[[2]]. Default is to FALSE, if set to TRUE computational time increases significantly.
#'@param vdj.only Logical, if set to TRUE a vdj dataframe can be used as input (vgm = vdj_dataframe). A vdj dataframe with integrated public information is returned instead of a new vgm. Default is FALSE
#'
#'Note: Aberrant numbers of V(D)J chains are filtered in the vgm



find_public_cells <- function(vgm,
                              vdj_public,
                              define.public,
                              check.overlap.type,
                              integrate.into.gex,
                              vdj.only){


  if(missing(define.public)) define.public <- VDJ_cdr3s_aa
  if(missing(integrate.into.gex)) integrate.into.gex <- FALSE
  if(missing(check.overlap.type)) check.overlap.type <- FALSE
  if(missing(vdj.only)) vdj.only <- FALSE

  #Extract and filter our vdj information

  if(vdj.only == FALSE){
    vdj_our_data <- vgm[[1]]
    #vdj_our_data <- subset(vdj_our_data, vdj_our_data$Nr_of_VDJ_chains == 1 & vdj_our_data$Nr_of_VJ_chains == 1)
  }

  if(vdj.only == TRUE){
    vdj_our_data <- vgm
  }


  if(check.overlap.type == TRUE){

    # 1. An overlap_type column is created which indicates if a cell has any cdr3 overlap with any of the cells in the public data (2nd dataset)
    vdj_our_data$overlap_type <- NA

    # 2. Different types of cdr3 aa overlap information is stored in a vdj copy, only integrate overlap_type information into vdj_our_data
    tmp <- vdj_our_data
    tmp$CDRH3_aa <- NA
    tmp$CDRL3_aa <- NA
    tmp$CDRL3_CDRH3_aa <- FALSE

    #look at CDRH3_aa overlaps in public data
    for (k in 1:nrow(tmp)){
      if(tmp$VDJ_cdr3s_aa[k] != ""){ #check that there's actually a sequence present
        if(tmp$VDJ_cdr3s_aa[k] %in% vdj_public$VDJ_cdr3s_aa){
          tmp$CDRH3_aa[k] <- TRUE #there's a match
        }
        else{
          tmp$CDRH3_aa[k] <- FALSE
        }
      }else{
        tmp$CDRH3_aa[k] <- FALSE
      }
    }


    #look at CDRL3_aa overlaps
    for (k in 1:nrow(tmp)){
      if(tmp$VJ_cdr3s_aa[k] != ""){
        if(tmp$VJ_cdr3s_aa[k] %in% vdj_public$VJ_cdr3s_aa){
          tmp$CDRL3_aa[k] <- TRUE
        }
        else{
          tmp$CDRL3_aa[k] <- FALSE
        }
      }else{
        tmp$CDRL3_aa[k] <- FALSE
      }
    }


    #look if any clone has both
    for (k in 1:nrow(tmp)){
      if(tmp$VJ_cdr3s_aa[k] != "" & tmp$VDJ_cdr3s_aa[k] != ""){ #check for both
        if(tmp$VDJ_cdr3s_aa[k] %in% vdj_public$VDJ_cdr3s_aa){ #first check if clone even has a heavy chain overlap
          tmp_subset <- subset(vdj_public, vdj_public$VJ_cdr3s_aa == tmp$VJ_cdr3s_aa[k] &  vdj_public$VDJ_cdr3s_aa == tmp$VDJ_cdr3s_aa[k] )
          if(nrow(tmp_subset)>0){
            tmp$CDRL3_CDRH3_aa[k] <- TRUE
          }
          else{
            tmp$CDRL3_CDRH3_aa[k] <- FALSE
          }
        }else{
          tmp$CDRL3_CDRH3_aa[k] <- FALSE
        }
      }
    }


    #Store only overlap_type information in vdj_our_data

    for(i in 1:nrow(vdj_our_data)){
      if(tmp$CDRL3_aa[i] == TRUE){ #all with overlaps are set to true --> all other will still be FALSE
        vdj_our_data$overlap_type[i] <- "CDRL3_aa"
      }
    }

    for(i in 1:nrow(vdj_our_data)){
      if(tmp$CDRH3_aa[i] == TRUE){
        vdj_our_data$overlap_type[i] <- "CDRH3_aa"
      }
    }

    for(i in 1:nrow(vdj_our_data)){
      if(tmp$CDRL3_CDRH3_aa[i] == TRUE){
        vdj_our_data$overlap_type[i] <- "CDRL3_CDRH3_aa"
      }
    }

    vdj_our_data$overlap_type[is.na(vdj_our_data$overlap_type)] <- "No overlap"

  }




  #Which definition of a public cell is used?
  public_definition <- NULL


  if(define.public == "VDJ_cdr3s_aa"){ #cdr3 aa heavy chain match
    #Copy the column chosen for the public clone definition (cdr3 heavy chain aa sequence: VDJ_cdr3s_aa) in our and in public vdj data
    vdj_our_data$public_definition <- vdj_our_data$VDJ_cdr3s_aa
    vdj_public$public_definition <- vdj_public$VDJ_cdr3s_aa
  }

  if(define.public == "VDJ_cdr3s_nt"){ #cdr3 nt heavy chain match
    vdj_our_data$public_definition <- vdj_our_data$VDJ_cdr3s_nt
    vdj_public$public_definition <- vdj_public$VDJ_cdr3s_nt
  }

  if(define.public == "VJ_cdr3s_aa"){ #cdr3 aa light chain match
    vdj_our_data$public_definition <- vdj_our_data$VJ_cdr3s_aa
    vdj_public$public_definition <- vdj_public$VJ_cdr3s_aa
  }

  if(define.public == "VJ_cdr3s_aa.VDJ_cdr3s_aa"){ #cdr3 aa light+heavy chain match
    vdj_our_data$public_definition <- paste(vdj_our_data$VJ_cdr3s_aa, vdj_our_data$VDJ_cdr3s_aa, sep = ".")
    vdj_public$public_definition <- paste(vdj_public$VJ_cdr3s_aa, vdj_public$VDJ_cdr3s_aa, sep = ".")
  }



  #Create new columns. Descriptions below are for the case that cdrh3 is chosen to define a public cell.
  vdj_our_data$public <- NA #Column indicating whether this cells is public or not (TRUE/FALSE)
  vdj_our_data$counts <- NA #Column indicating how many cells from the public data have the same cdrh3 sequence
  vdj_our_data$nr_unique_sample_ids <- NA #Column indicating in how many different samples from the public data cdrh3 matches are occurring. Indicates how "public" this cdr3h sequence really is.
  vdj_our_data$matches_barcode <- NA #Saves barcodes of all cells within the public data which have the same cdrh3 sequence

  vdj_our_data$public_clonotype <- NA #Is any cell within this clonotype defined as public? --> TRUE


  #look at overlaps in public data
  for (k in 1:nrow(vdj_our_data)){
    if(vdj_our_data$public_definition[k] != ""){ #check if there's actually a sequence present (should actually not be necessary if cells with aberrant nr of V(D)J chains have been removed )
      if(vdj_our_data$public_definition[k] %in% vdj_public$public_definition){

        vdj_our_data$public[k] <- TRUE #there's a match

        tmp <- subset(vdj_public, vdj_public$public_definition == vdj_our_data$public_definition[k]) #create a df of all cells in the public data with a cdrh3 match to the current cell
        vdj_our_data$counts[k] <- nrow(tmp) #Number of rows in this df indicate how many cdr3h matches are found in the public dataset
        vdj_our_data$nr_unique_sample_ids[k] <- length(unique(tmp$sample_id)) #In how many different samples of the public data are the matching cells found?

        barcodes <- tmp$barcode
        vdj_our_data$matches_barcode[k] <- list(barcodes)
      }
      else{
        vdj_our_data$public[k] <- FALSE
      }
    }else{
      vdj_our_data$public[k] <- FALSE
    }
  }


  if(vdj.only == TRUE){
    return(vdj_our_data)
  }


  #Add "public_clone" information to every cell. If one cell within a clone is defined as public, the "public_clone" status of every cell within this clone is set to TRUE.

  #Additionally: A new df containing information on clonotype public status is created. Contains information on sample_id, clonotype_id_10x and the public status of this clonotype.

  public_clonotypes <- NULL

  #To do this: Go through every clonotype within your data. (This needs to be done for every sample individually since clonotype names overlap between samples)

  sample_ids <- unique(vdj_our_data$sample_id)


  for (i in 1:length(sample_ids)){ #Go through every sample separately
    sample <- subset(vdj_our_data, vdj_our_data$sample_id == sample_ids[i])
    clonotypes <- unique(sample$clonotype_id_10x)

    public_clonotypes_sample <- NULL
    public_clonotypes_sample <- data.frame(matrix(ncol = 3, nrow = length(clonotypes)))
    x <- c("sample_id", "clonotype_id_10x", "public")
    names(public_clonotypes_sample) <- x

    public_clonotypes_sample$sample_id <- sample_ids[i]

    for(j in 1:length(clonotypes)){ #Look at each clonotype within the current sample. Are there any public cells within this clonotype?
      curr_clone <- subset(sample, sample$clonotype_id_10x == clonotypes[j])

      public_clonotypes_sample$clonotype_id_10x[j] <- clonotypes[j]

      if(any(curr_clone$public == TRUE)){
        vdj_our_data$public_clonotype[vdj_our_data$sample_id == sample_ids[i] & vdj_our_data$clonotype_id_10x == clonotypes[j]] <- TRUE
        public_clonotypes_sample$public[j] <- TRUE
      }else{
        vdj_our_data$public_clonotype[vdj_our_data$sample_id == sample_ids[i] & vdj_our_data$clonotype_id_10x == clonotypes[j]] <- FALSE
        public_clonotypes_sample$public[j] <- FALSE
      }
    }

    public_clonotypes <- rbind(public_clonotypes, public_clonotypes_sample)
  }



  ##Integration of public clone information into the gene expression data (vgm[[2]])

  #reintegration into vgm[[1]]
  vgm[[1]] <- vdj_our_data

  if(integrate.into.gex == TRUE){

    vgm[[2]]$public <- rep("",ncol(vgm[[2]]))
    vgm[[2]]$nr_unique_sample_ids <- rep("",ncol(vgm[[2]])) #Also add information on how "public" this public cell is. In how many different samples (mice/organs) does this cells cdr3 sequence occur (e.g. cdrh3 aa sequence, depending on user input for public status definition)
    vgm[[2]]$public_clonotype <- rep("",ncol(vgm[[2]]))

    bc_match_index <- NULL


    for (i in 1:ncol(vgm[[2]])){
      #if you have a match vgm[[1]]$orig_barcode %in% vgm[[2]]$orig_barcode[i] will return one TRUE and the rest FALSE --> so length(unique()) will be 2
      if (length(unique(vgm[[1]]$orig_barcode %in% vgm[[2]]$orig_barcode[i])) == 1) { #if no match is found only FALSE is returned for each barcode --> then there's only false in the unique() list
        bc_match_index[i] <- NA
      }else{
        bc_match_index[i] <- which(vgm[[1]]$orig_barcode == vgm[[2]]$orig_barcode[i]) #which() returns index of the vgm[[1]] original b
        vgm[[2]]$public[i] <- vgm[[1]]$public[bc_match_index[i]]
        vgm[[2]]$public_clonotype[i] <- vgm[[1]]$public_clonotype[bc_match_index[i]]
        vgm[[2]]$nr_unique_sample_ids[i] <- vgm[[1]]$nr_unique_sample_ids[bc_match_index[i]]
      }
    }
  }

  tmp <- list(vgm, public_clonotypes)

  return(tmp)

}
