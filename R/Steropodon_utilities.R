#' Annotates a bio3d pdb structure's atoms or residues

#' @description Annotates a Steropodon structure's atoms or residues
#' @param pdb bio3d pdb object to be annotated
#' @param res.idx list or vector of integers - residue IDs to be annotated
#' @param atom.idx list or vector of integers - atom IDs to be annotated
#' @param feature list or vector of strings - features/labels to be annotated
#' @return an annotated bio3d object
#' @export
#' @examples
#' \dontrun{
#' annotate_structure(pdb = pdb, res.idx = res.idx, feature = feature)
#'}
annotate_structure <- function(pdb, res.idx, atom.idx, feature){
  if(missing(pdb)) stop('Please input your PDB structure')
  if(missing(res.idx)) res.idx <- NULL
  if(missing(atom.idx)) atom.idx <- NULL
  if(missing(feature)) feature <- 'new_feature'

  if(!is.null(res.idx)){
    pdb$atom$resno <- as.character(pdb$atom$resno)
    pdb$atom[[feature]] <- res.idx[pdb$atom$resno]
    pdb$atom$resno <- as.integer(pdb$atom$resno)

  }else if(!is.null(atom.idx)){
    pdb$atom$eleno <- as.character(pdb$atom$eleno)
    pdb$atom[[feature]] <- atom.idx[pdb$atom$eleno]
    pdb$atom$eleno <- as.integer(pdb$atom$eleno)
  }

  return(pdb)
}


#' Trims a bio3d pdb object given a grouping and a specific value

#' @description Splits a bio3d pdb object given a grouping and a specific value
#' @param pdb bio3d pdb object to split
#' @param grouping string or vector of strings - grouping factor for the split (e.g., 'chain' will only split chains such as VDJ and VJ/ heavy and light).
#' @param specific.values string or vector of strings - the specific regions to kept in the structure, depending on the 'grouping' parameter (e.g., grouping = 'chain' and specific.values = 'VDJ' will remove all keep chains).
#' @param combine.groupings bool - if TRUE, the regions inserted in specific.values will all be combined in the same final structure.
#' @param combine.values bool - if TRUE, will combine the 'grouping' features for a more specific trimming (e.g., 'chain' and 'region' to select exactly which hypervariable regions from any chain should be kept - 'VDJ_CDR1', 'VJ_CDR1' etc.).
#' @return a trimmed bio3d object, with regions kept as dictated by the specific.values parameter.
#' @export
#' @examples
#' \dontrun{
#' only_heavy <- split_structure(pdb = pdb,
#'   grouping = 'chain',
#'   specific.values = 'VDJ')
#'}
split_structure <- function(pdb, grouping, specific.values, combine.groupings, combine.values){
  if(missing(pdb)) stop('Please input your PDB structure')
  if(missing(grouping)) grouping <- 'chain'
  if(missing(specific.values)) specific.values <- NULL
  if(missing(combine.groupings)) combine.groupings <- F
  if(missing(combine.values)) combine.values <- F

  combined <- NULL

  is_null_pdb <- function(pdb){
    out <- F
    if(nrow(pdb$atom) == 0){
      out <- T
    }
    out
  }

  if(combine.groupings & (length(grouping) > 1)){
    pdb$atom$combined <- rep('', nrow(pdb$atom))
    for(group in grouping){
      pdb$atom$combined <- paste0(pdb$atom$combined, '_', pdb$atom[[group]])
    }
    pdb$atom <- pdb$atom %>%
                dplyr::mutate(combined = substring(combined,2))

    grouping <- 'combined'
  }

  if(is.null(specific.values)){
    specific.values <- unique(pdb$atom[[grouping]])
  }

  if(length(specific.values) > 1 & !combine.values){
    out_pdbs <- vector(mode = 'list', length = length(specific.values))
    for(i in 1:length(specific.values)){
      res_idx <- pdb$atom$resno[pdb$atom[[grouping]] %in% specific.values[i]]
      idx <- bio3d::atom.select(pdb, resno = res_idx)
      trimmed_struct <- bio3d::trim(pdb, inds = idx)
      trimmed_struct$atom$combined <- NULL
      out_pdbs[[i]] <- trimmed_struct
    }
    names(out_pdbs) <- specific.values
    out_pdbs <- out_pdbs[!sapply(out_pdbs, is_null_pdb)]

  }else{
    res_idx <- pdb$atom$resno[pdb$atom[[grouping]] %in% specific.values]
    idx <- bio3d::atom.select(pdb, resno = res_idx)
    trimmed_struct <- bio3d::trim(pdb, inds = idx)
    trimmed_struct$atom$combined <- NULL
    out_pdbs <- trimmed_struct

    if(is_null_pdb(out_pdbs)){
      return(NULL)
    }
  }

  return(out_pdbs)
}


#' Turns a Steropodon nested list into a single list

#' @description Turns a Steropodon nested list into a single list
#' @param steropodon.object  a nested list of predicted structure objects (per sample, per clonotype)
#' @return a single list of Steropodon objects
#' @export
#' @examples
#' \dontrun{
#' steropodon_list <- unnest_steropodon(steropodon_object)
#'}
unnest_steropodon <- function(steropodon.object){
  custom_flatten <- function (x, use.names = TRUE, classes = "ANY"){
    len <- sum(rapply(x, function(x) 1L, classes = classes))
    y <- vector("list", len)
    i <- 0L
    items <- rapply(x, function(x) {
      i <<- i + 1L
      y[[i]] <<- x
      TRUE
    }, classes = classes)
    if (use.names && !is.null(nm <- names(items)))
      names(y) <- nm
    y
  }

  if(inherits(steropodon.object[[1]], 'list')){
    steropodon.object <- custom_flatten(steropodon.object)
  }

  return(steropodon.object)
}


#' Converts a single list of Steropodon objects into a nested list (per sample, per clonotype)

#' @description Converts a single list of Steropodon objects into a nested list (per sample, per clonotype)
#' @param steropodon.list named list - Steropodon objects unnested using unnest_steropodon. The list should be named, including the sample and clonptype IDs for easier nesting.

#' @return a nested list of Steropodon objects (per sample, per clonotype).
#' @export
#' @examples
#' \dontrun{
#' steropodon_object <- nest_steropodon(steropodon_list)
#'}
nest_steropodon <- function(steropodon.list){
  get_barcode_dataframe <- function(steropodon.list){
    unique_names <- names(steropodon.list)
    unique_names <- lapply(unique_names, function(x) stringr::str_split(x, '\\.')[[1]])
    df <- data.frame(do.call('rbind', unique_names))
    colnames(df) <- c('sample', 'clonotype', 'id')

    return(df)
  }

  df <- get_barcode_dataframe(steropodon.list)

  samples <- unique(df$sample)
  steropodon_nested <- vector(mode = 'list', length = length(samples))

  for(i in 1:length(samples)){
    df_subset <- df[df$sample == samples[i],]
    clonotypes <- unique(df_subset$clonotype)

    steropodon_nested[[i]] <- vector(mode = 'list', length = length(clonotypes))

    for(j in 1:length(clonotypes)){
      df_subset2 <- df_subset[df_subset$clonotype == clonotypes[j],]
      ids <- unique(df_subset2$id)

      steropodon_nested[[i]][[j]] <- vector(mode = 'list', length = length(ids))

      for(k in 1:length(ids)){
        steropodon_nested[[i]][[j]][[k]] <- steropodon.list[[paste0(c(samples[i], clonotypes[j], ids[k]), collapse = '.')]]
      }
      names(steropodon_nested[[i]][[j]]) <- ids
    }
    names(steropodon_nested[[i]]) <- clonotypes
  }
  names(steropodon_nested) <- samples


  return(steropodon_nested)
}


#' Selects a single bio3d pdb object from a given Steropodon object slot/attribute
#' @description Selects a single bio3d pdb object from a given Steropodon object slot/attribute
#' @param steropodon.object a single Steropdon object, as obtained from Steropodon_model.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @return a bio3d pdb object
#' @export
#' @examples
#' \dontrun{
#' pdb <- select_structure(steropodon_igfold[[1]][[1]], structure = 'CDRH3')
#'}
select_structure <- function(steropodon.object, structure){

    if(structure == 'structure'){
      pdb <- steropodon.object@structure
    }else if(structure == 'H'){
      pdb <- steropodon.object@H
    }else if(structure == 'L'){
      pdb <- steropodon.object@L
    }else if(structure == 'CDRH3'){
      pdb <- steropodon.object@CDRH3
    }else if(structure == 'CDRL3'){
      pdb <- steropodon.object@CDRL3
    }else if(structure == 'complex'){
      pdb <- steropodon.object@complex
    }else if(structure == 'paratope'){
      pdb <- steropodon.object@paratope
    }else if(structure == 'epitope'){
      pdb <- steropodon.object@epitope
    }else if(structure == 'core'){
      pdb <- steropodon.object@core
    }else if (structure == 'pdbs'){
      pdb <- steropodon.object@pdbs
    }else{
      stop(paste0('Could not find structure ', structure, ' in your Steropodon object!'))
    }

    #if(is.null(structure)){
    #  stop(paste0('Could not find structure ', structure, ' in your Steropodon object!'))
    #}

    return(pdb)
}



#' Modifies an attribute/slot of a single Steropodon object

#' @description Modifies an attribute/slot of a single Steropodon object
#' @param steropodon.object a single Steropdon object, as obtained from Steropodon_model.
#' @param pdb bio3d pdb object - will replace the object present in the attribute determined by the 'structure' parameter.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @return a Steropodon object with the attribute/slot named in 'structure' replaced with the bio3d pdb structure from 'pdb'.
#' @export
#' @examples
#' \dontrun{
#' modified_steropodon <- modify_structure(steropodon_object,
#'                                         pdb = pdb,
#'                                         structure = 'CDRH3')
#'}
modify_structure <- function(steropodon.object, pdb, structure){

    if(structure == 'structure'){
      steropodon.object@structure <- pdb
    }else if(structure == 'H'){
      steropodon.object@H <- pdb
    }else if(structure == 'L'){
      steropodon.object@L <- pdb
    }else if(structure == 'CDRH3'){
      steropodon.object@CDRH3 <- pdb
    }else if(structure == 'CDRL3'){
      steropodon.object@CDRL3 <- pdb
    }else if(structure == 'complex'){
      steropodon.object@complex <- pdb
    }else if(structure == 'paratope'){
      steropodon.object@paratope <- pdb
    }else if(structure == 'epitope'){
      steropodon.object@epitope <- pdb
    }else if(structure == 'core'){
      steropodon.object@core <- pdb
    }else if(structure == 'pdbs'){
      steropodon.object@pdbs <- pdb
    }else{
      stop(paste0('Could not find structure ', structure, ' in your Steropodon object!'))
    }

    #if(is.null(structure)){
    #  stop(paste0('Could not find structure ', structure, ' in your Steropodon object!'))
    #}
    return(steropodon.object)
}


#' Performs sequence and iterative structural alignment

#' @description Sequence and structural alignment of bio3d xyz coordinates, similar to PymMOL's 'align' command.
#' Will first perform sequence alignment, then structural superposition on the aligned coordinates for n cycles ('max.cycles' parameter), removing outliers based on a cutoff ('cutoff' parameter)
#' @param mobile bio3d xyz coordinates - coordinates to be aligned
#' @param fixed bio3d xyz coordinates - template coordinates
#' @param mobile.inds vector of integers - C-alpha atom indices for the mobile coordinates.
#' @param fixed.inds vector of integers - C-alpha atom indices for the fixed coordinates.
#' @param aln.method string - sequence alignment algorithm. Currently, only MAFFT is implemented (aln.method = 'mafft').
#' @param max.cycles integer - the maximum number of iterations (superposition followed by outlier rejection) to be done in the sequence alignment and iterative structural superposition algorithm.
#' @param cutoff float - the distance cutoff at which outliers will be rejected in the sequence alignment and iterative structural superposition algorithm.
#' @param return.mobile boolean - if TRUE, will only return the mobile coordinates. Else, will return both fixed and nobile as bio3d pdb objects.
#' @return a nested list of Steropodon structures with all regions not included in 'specific.values' removed. For more information, see the Steropodon vignette (https://alexyermanos.github.io/Platypus/articles/Steropodon.html).
#' @export
#' @examples
#' \dontrun{
#' aligned <- sequence_structure_superpose(mobile_pdb,
#' fixed_pdb,
#' max.cycles = 10,
#' cutoff = 0.5)
#'}
sequence_structure_superpose <- function(mobile,
                                         fixed,
                                         mobile.inds  = NULL,
                                         fixed.inds = NULL,
                                         aln.method = 'mafft',
                                         max.cycles = 10,
                                         cutoff = 0.5,
                                         return.mobile = F
                                        ){

    #Utilities
    xyz.dist <- function(v){
      a <- v[1:3]
      b <- v[4:6]

      return(sqrt(sum((a-b)**2)))
    }

    resi.dev <- function(xyz.a, xyz.b, cycle=1, cutoff = 0.5){
      k <- matrix(xyz.a, ncol = 3, byrow = T)
      l <- matrix(xyz.b, ncol = 3, byrow = T)

      devs <- apply(cbind(k,l), 1, "xyz.dist")
      m <- stats::median(devs)
      std <- stats::sd(devs)

      cut <- m + (2*std)
      inds <- which(devs > cut)

      if((std < cutoff) | (length(inds)==0)){
        return(NULL)
      }else{
        cat(" Cycle ", i, ": ", length(inds), " atoms rejected", "\n", sep="")
        cat(" Mean: ", round(m, 1),
            " Std: ", round(std, 1),
            " Cut: ", round(cut, 1), "\n", sep="" )
        return(inds)
      }
    }

    remap.inds <- function(pdb.init, inds.init, inds.trunc.atom){
      ## Map back to indices for the entire PDB given
      inds.full <- NULL
      inds.full$atom <- inds.init$atom[inds.trunc.atom]
      inds.full$xyz <- bio3d::atom2xyz(inds.full$atom)
      inds.full$logical <- bio3d::atom2xyz(1:nrow(pdb.init$atom)) %in% inds.full$xyz
      return(inds.full)
    }

    parse.pdb <- function(pdb, gaps, s, i){
      pdbseq <- bio3d::aa321(pdb$atom[pdb$calpha, "resid"])
      aliseq <- toupper(s$ali[i, ])
      tomatch <- gsub("X", "[A-Z]", aliseq[!bio3d::is.gap(aliseq)])
      start.num <- regexpr(pattern = paste(c(stats::na.omit(tomatch[1:15])), collapse = ""), text = paste(pdbseq, collapse = ""))[1]

      nseq <- rep(NA, length(aliseq))
      ali.res.ind <- which(!bio3d::is.gap(aliseq))

      ali.res.ind <- ali.res.ind[1:length(pdbseq)]
      nseq[ali.res.ind] = start.num:((start.num - 1) + length(tomatch))

      pdb$atom <- cbind(pdb$atom, index=seq(1, nrow(pdb$atom)))
      ca.ali <- pdb$atom[pdb$calpha, ][nseq, ]
      at.inds <- ca.ali[, "index"]
      return(at.inds)
    }


    #Main function
    if(!is.null(fixed.inds)){
      if(length(fixed.inds$atom) < 2){
        stop("align: insufficent atom indices for fitting")
      }

      a <- bio3d::trim.pdb(fixed, fixed.inds)

    }else{
      a <- fixed
      fixed.inds <- bio3d::atom.select(fixed, 'all', verbose=FALSE)
    }

    if(!is.null(mobile.inds)){
      if(length(mobile.inds$atom) < 2){
        stop("align: insufficent atom indices for fitting")
      }

      b <- bio3d::trim.pdb(mobile, mobile.inds)

    }else{
      b <- mobile
      mobile.inds <- bio3d::atom.select(mobile, 'all', verbose=FALSE)
    }

    ## PDB list for sequence alignment
    pdb.list <- NULL
    pdb.list[[1]] <- a
    pdb.list[[2]] <- b

    ## Sequence alignment - MODIFIED AS MUSCLE COULD NOT BE FOUND/ ERRORS W SEQALN - ALSO MAFFT IS NOTICEABLY FASTER
    ## WILL ADD MORE ALIGNMENT METHODS IN THE FUTURE
    s <- lapply(pdb.list, bio3d::pdbseq)
    temp_file <- tempfile()
    fasta_name <- paste0(temp_file, '.fasta')
    aligned_name <- paste0(temp_file, '_aligned.fasta')
    seqinr::write.fasta(as.list(s), names=c('fixed', 'mobile'), file.out = fasta_name)

    if(aln.method == 'mafft'){
      system(paste0('mafft ', fasta_name, ' > ', aligned_name))
    }else{
      stop('Alignment method not implemented yet!')
    }

    s <- bio3d::read.fasta(aligned_name)

    unlink(fasta_name)
    unlink(aligned_name)

    gaps <- bio3d::gap.inspect(s$ali)

    ## Parse truncated PDBs
    at.inds.a <- parse.pdb(a, gaps, s, 1)
    at.inds.b <- parse.pdb(b, gaps, s, 2)

    ## Fetch indices for fitting (truncated pdb)
    at.a <- as.numeric(at.inds.a[gaps$f.inds])
    at.b <- as.numeric(at.inds.b[gaps$f.inds])

    ## Indices for full pdb - done with the truncated ones
    a.inds.full <- remap.inds(fixed, fixed.inds, at.a)
    b.inds.full <- remap.inds(mobile, mobile.inds, at.b)

    ## Perform the initial fitting
    fit <- bio3d::rot.lsq(mobile$xyz, fixed$xyz, xfit = b.inds.full$logical, yfit = a.inds.full$logical)

    rmsd.init <- bio3d::rmsd(as.vector(fixed$xyz), fit, a.inds = a.inds.full$xyz, b.inds = b.inds.full$xyz)
    cat("\n")
    cat(" Initial RMSD (", length(gaps$f.inds), " atoms): ", rmsd.init, "\n", sep="")


    ## Refinement process
    rmsd.all <- c(rmsd.init)
    for(i in 1:max.cycles){
      if(i > max.cycles){
        break
      }
      ## Find residues with largest structural deviation
      exc <- resi.dev(fixed$xyz[a.inds.full$xyz], fit[b.inds.full$xyz], cycle = i, cutoff = cutoff)

      if(is.null(exc)){
        break
      }else{
        ## Remove atoms for new round of fitting
        exc <- bio3d::atom2xyz(exc)

        tmp <- 1:length(a.inds.full$logical)
        exc.a <- tmp[which(a.inds.full$logical)][exc]
        a.inds.full$logical[exc.a] <- FALSE

        tmp <- 1:length(b.inds.full$logical)
        exc.b <- tmp[which(b.inds.full$logical)][exc]
        b.inds.full$logical[exc.b] <- FALSE

        ## Build new xyz and atom indices
        a.inds.full$xyz <- which(a.inds.full$logical)
        b.inds.full$xyz <- which(b.inds.full$logical)

        a.inds.full$atom <- bio3d::xyz2atom(a.inds.full$xyz)
        b.inds.full$atom <- bio3d::xyz2atom(b.inds.full$xyz)

        ## Fit based on new indices
        fit <- bio3d::rot.lsq(mobile$xyz, fixed$xyz, xfit = b.inds.full$logical, yfit = a.inds.full$logical)

        ## Calculate RMSD
        tmp.rmsd <- bio3d::rmsd(as.vector(fixed$xyz), fit, a.inds = a.inds.full$xyz, b.inds = b.inds.full$xyz)
        rmsd.all <- c(rmsd.all, tmp.rmsd)
        num.resi <- length(which(a.inds.full$logical))/3

        cat("  RMSD (", num.resi, " of ", length(gaps$f.inds), " atoms): ", tmp.rmsd, "\n", sep="")
      }
    }

    a.inds.full$logical <- NULL
    b.inds.full$logical <- NULL

    mobile$xyz <- bio3d::as.xyz(fit)
    mobile$atom[c('x', 'y', 'z')] <- matrix(mobile$xyz, ncol = 3, byrow = T)
    mobile$superposition.indices <- b.inds.full

    fixed$superposition.indices <- a.inds.full

    if(return.mobile){
      return(mobile)
    }else{
      return(list(mobile = mobile, fixed = fixed))
    }
}


#' Performs a single iteration of structural alignment (Kabsch algorithm)

#' @description Performs a single iteration of structural alignment (Kabsch algorithm)
#' @param fixed bio3d xyz coordinates - template coordinates
#' @param mobile bio3d xyz coordinates - coordinates to be aligned
#' @param mobile.inds vector of integers - C-alpha atom indices for the mobile coordinates.
#' @param fixed.inds vector of integers - C-alpha atom indices for the fixed coordinates.
#' @return mobile coordinates aligned to the fixed ones, as a bio3d xyz object.
#' @export
#' @examples
#' \dontrun{
#' aligned <- sequence_structure_superpose(mobile_pdb,
#' fixed_pdb)
#'}
structure_superpose <- function(fixed,
                                mobile,
                                mobile.inds = NULL,
                                fixed.inds = NULL){

  xyz <- suppressWarnings(bio3d::fit.xyz(fixed = fixed$xyz,
                                         mobile = mobile$xyz,
                                         mobile.inds = mobile.inds,
                                         fixed.inds = fixed.inds))

  mobile$xyz <- xyz
  mobile$atom[c('x', 'y', 'z')] <- matrix(mobile$xyz, ncol = 3, byrow = T)
  return(mobile)
}


#' Extract framework/hypervariable regions annotated with VDJ_call_MIXCR

#' @description Extract framework/hypervariable regions annotated with VDJ_call_MIXCR
#' @param VDJ.matrix VDJ matrix, VGM[[1]].
#' @param chain.to.extract string - 'VDJ' for heavy, 'VJ' for light, 'VDJ.VJ' for both.
#' @param as.nucleotide bool - TRUE for nucleotides, FALSE for amino acids.
#' @param regions.to.extract vector of strings - framework/hypervariable region to extract and paste together ('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4').

#' @return VDJ dataframe with a new column - 'mixcr_assembled' with the extracted MIXCR regions.
#' @export
#' @examples
#' \dontrun{
#' sequences <- extract_MIXCR(VDJ, chain.to.extract = 'VDJ')
#'}
extract_MIXCR <- function(VDJ.matrix, chain.to.extract, as.nucleotide, regions.to.extract){
      if(missing(VDJ.matrix)) stop('Input the VDJ dataframe obtained after calling VDJ_call_MIXCR')
      if(missing(chain.to.extract)) chain.to.extract <- 'VDJ.VJ'
      if(missing(as.nucleotide)) as.nucleotide <- F
      if(missing(regions.to.extract)) regions.to.extract <- list('FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4')


      VDJ.matrix$mixcr_assembled <- rep('', nrow(VDJ.matrix))

      if(chain.to.extract!='VDJ.VJ'){
        if(as.nucleotide){
          col_name <- paste0(chain.to.extract, '_nSeq')
        }else{
          col_name <- paste0(chain.to.extract, '_aaSeq')
        }

        for(region in regions.to.extract){
          VDJ.matrix$mixcr_assembled <- paste0(VDJ.matrix$mixcr_assembled, gsub('_', '', VDJ.matrix[, paste0(col_name, region)]))
        }
      }else if(chain.to.extract=='VDJ.VJ'){
        if(as.nucleotide==T){
          col_name <- '_nSeq'
        }else{
          col_name <- '_aaSeq'
        }
        extracted_VDJ <- rep('', nrow(VDJ.matrix))
        extracted_VJ <- rep('', nrow(VDJ.matrix))

        for(region in regions.to.extract){
          extracted_VDJ <- paste0(extracted_VDJ, gsub('_', '', VDJ.matrix[,paste0('VDJ',col_name, region)]))
          extracted_VJ <- paste0(extracted_VJ, gsub('_', '', VDJ.matrix[,paste0('VJ',col_name, region)]))
        }

        VDJ.matrix$mixcr_assembled <- paste0(extracted_VDJ, ':', extracted_VJ)

      }
      return(VDJ.matrix)
}


#' Sequence aligment using MAFFT

#' @description Sequence aligment using MAFFT
#' @param structure.list list - list of Steropodon structures as bio3d pdb objects
#' @param alignment.method string - alignment method. Currently, only MAFFT is supported (alignment.method = 'mafft').
#' @return a bio3d pdbs object of aligned sequences and their corresponding coordinates
#' @export
#' @examples
#' \dontrun{
#' aligned <- sequence_alignment(structures, alignment.method = 'mafft')
#'}
sequence_alignment <- function(structure.list,
                               alignment.method
                               ){

  sequences <- lapply(structure.list, function(x) paste0(bio3d::pdbseq(x), collapse = ''))
  seqinr::write.fasta(as.list(sequences), names = names(sequences), file.out='temp.fasta')

  if(alignment.method == 'mafft'){
    system(paste0('mafft ', 'temp.fasta', ' > ', 'aligned.fasta'))
  }else{
    stop('Method not implemented yet!')
  }

  aln <- bio3d::read.fasta('aligned.fasta')
  unlink('temp.fasta')
  unlink('aligned.fasta')

  pdbs <- bio3d::read.fasta.pdb(aln = aln, pdblist = structure.list)
  pdbs$xyz <- pdbs$xyz[]

  return(pdbs)
}


#' Saves a Steropodon object's structure as a PDB file
#' @description Saves a Steropodon object's structure as a PDB file

#' @param steropodon.object a single Steropodon object or a nested list (per sample, per clonotype).
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param dir string - path to the directory for saving the PDB files
#' @return no returns. Will save the PDBs in the specified directory.
#' @export
#' @examples
#' \dontrun{
#' write_steropodon_pdbs(steropodon_igfold,
#' structure = 'CDRH3',
#' dir = './pdbs')
#'}
write_steropodon_pdbs <- function(steropodon.object, structure, dir){

  if(!is.null(dir)){
    if(!dir.exists(dir)) dir.create(dir)
  }

  if(!inherits(steropodon.object, 'Steropodon')){
    steropodon_list <- unnest_steropodon(steropodon.object)
    structure_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
    file_names <- unlist(lapply(names(steropodon_list), function(x) stringr::str_replace_all(x, '\\.', '_')))
    if(!is.null(dir)){
      file_list <- paste0(dir, '/', file_names, '.pdb')
    }else{
      file_list <- paste0(file_names, '.pdb')
    }
    chains <- unique(structure_list[[1]]$atom$chain)
    chain_dict <- toupper(letters)[1:length(chains)]
    names(chain_dict) <- chains
    mapply(function(x,y) bio3d::write.pdb(xyz = x$xyz,
                                          type = x$atom$type,
                                          resno = x$atom$resno,
                                          resid = x$atom$resid,
                                          eleno = x$atom$eleno,
                                          elety = x$atom$elety,
                                          chain = chain_dict[x$atom$chain],
                                          o = x$atom$o,
                                          b = x$atom$b,
                                          segid = x$atom$segid,
                                          elesy = x$atom$elesy,
                                          charge = x$atom$charge,
                                          file = y),
                                            structure_list, file_list)

  }else{
    pdb <- select_structure(steropodon.object, structure = structure)
    file_names <- paste0(steropodon.object@structure_id, collapse = '_')

    if(!is.null(dir)){
      file_list <- paste0(dir, '/', file_names, '.pdb')
    }else{
      file_list <- paste0(file_names, '.pdb')
    }

    chains <- unique(pdb$atom$chain)
    chain_dict <- toupper(letters)[1:length(chains)]
    names(chain_dict) <- chains
    bio3d::write.pdb(xyz = pdb$xyz,
                        type = pdb$atom$type,
                        resno = pdb$atom$resno,
                        resid = pdb$atom$resid,
                        eleno = pdb$atom$eleno,
                        elety = pdb$atom$elety,
                        chain = chain_dict[pdb$atom$chain],
                        o = pdb$atom$o,
                        b = pdb$atom$b,
                        segid = pdb$atom$segid,
                        elesy = pdb$atom$elesy,
                        charge = pdb$atom$charge,
                        file = file_list)
  }

  return(list(chain_dict = chain_dict, file_list = file_list))
}


#' Removes NA values from a Steropodon structure
#' @description Removes NA values from a Steropodon structure
#' @param steropodon.object a single Steropodon object or a nested list (per sample, per clonotype).
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @return the Steropodon object with the structure in the slot specified by 'structure' without any NA values.
#' @export
#' @examples
#' \dontrun{
#' no_na <- remove_na_steropodon(steropodon.object,
#' structure = 'CDRH3')
#'}
remove_na_steropodon <- function(steropodon.object, structure){
  remove_na_pdb <- function(pdb){
    na_residues <- unique(pdb$atom$resno[is.na(pdb$atom$chain)])
    if(length(na_residues) > 0){
      pdb <- bio3d::trim(pdb, inds = bio3d::atom.select(pdb, resno = na_residues))
    }

    return(pdb)
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  pdbs <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
  pdbs <- lapply(pdbs, function(x) remove_na_pdb(x))
  steropodon_list <- mapply(function(x,y) modify_structure(x, structure = structure, pdb = y), steropodon_list, pdbs)

  steropodon_object <- nest_steropodon(steropodon_list)
}
