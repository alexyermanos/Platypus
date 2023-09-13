#' Calculates RMSD/LDDT/TM-score/custom scores between all structures in a Steropodon nested list


#' @description Computes a distance matrix using the structural distance scores defined in the distance.metric parameter.
#' Options include: 'rmsd' for calculating the root mean square deviation between all coordinates, 'lddt' for the local distance difference test, 'tmscore' for the template modelling score, 'custom' for a weighted RMSD score (weighted by chains or framework/variable regions).

#' @param steropodon.object a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param distance.metric string - the distance metric to be used. Options include: 'rmsd' for calculating the root mean square deviation between all coordinates, 'lddt' for the local distance difference test, 'tmscore' for the template modelling score, 'custom' for a weighted RMSD score (weighted by chains or framework/variable regions).
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param seq.struct.superpose bool - if TRUE, will perform a sequence alignment followed by an iterative structural superposition (removing outlier atoms in the fit). This is similar to the 'align' command in PyMOL.
#' @param struct.superpose bool - if TRUE, will perform a single structure superposition/ Kabsch algorithm iteration.
#' @param compare.cores bool - if TRUE, will compare the core/ shared structural regions across all structures as obtained from the Steropoodon_find_cores function.
#' @param cutoff flaot - the distance cutoff at which outliers will be rejected in the sequence alignment and iterative structural superposition algorithm (seq.struct.superpose = TRUE).
#' @param max.cycles integer - the maximum number of iterations (superposition followed by outlier rejection) to be done in the sequence alignment and iterative structural superposition algorithm (seq.struct.superpose = TRUE).
#' @param alignment.method string - sequence alignment method to be used when seq.struct.superpose = TRUE. Currently only MAFFT is implemented (alignment.method = 'mafft').
#' @param custom.specific.values vector of strings - VDJ or VJ regions to be considered in the custom distance algorithm (distance.metric = 'custom'). If only specific regions for a given chain are to be picked, the 'features' parameter should include both 'chain' and 'region'.
#' @param custom.features vector of strings or string - the sequence features to be considered for the distance calculation. 'region' only for FR1/CDR1/etc regions (across all chains), 'chain' for either 'VDJ' or 'VJ' chains (inserted in specific.values).
#' @param custom.weights vector of integer - weight to be assigned in the custom score for each region/chain included in the custom.specific.values parameter.
#' @param custom.normalize.sequence.length bool - if TRUE, will normalize the score per region by that region's number of amino acids.
#' @param custom.commutative bool - if TRUE, will normalize the distance by the length of both sequence lengths of the target and template structures (pdb1 and pdb2) used in the distance calculation.
#' @param plot.results bool - if TRUE, will output a heatmap of the distance matrix.


#' @return the inter-structure distance matrix or a heatmap of the distance values.
#' @export
#' @examples
#' \dontrun{
#' superposed_seq_struct <- Steropodon_superposition(steropodon_igfold,
#' sequence.structure.superpose = T,
#' structure.superpose = F,
#' max.cycles = 10, cutoff = 0.5)
#'
#' distance_matrix <- superposed_seq_struct %>%
#' Steropodon_distances(distance.metric = 'rmsd', plot.results = F)
#'}


Steropodon_distances <- function(steropodon.object,
                                 distance.metric,
                                 structure,
                                 seq.struct.superpose,
                                 struct.superpose,
                                 compare.cores,
                                 cutoff,
                                 max.cycles,
                                 alignment.method,
                                 custom.specific.values,
                                 custom.features,
                                 custom.weights,
                                 custom.normalize.sequence.length,
                                 custom.commutative,
                                 plot.results
                                ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object!')
  if(missing(distance.metric)) distance.metric <- 'custom'
  if(missing(structure)) structure <- 'structure'
  if(missing(seq.struct.superpose)) seq.struct.superpose <- T
  if(missing(struct.superpose)) struct.superpose <- F
  if(missing(compare.cores)) compare.cores <- F
  if(missing(cutoff)) cutoff <- 0.5
  if(missing(max.cycles)) max.cycles <- 10
  if(missing(alignment.method)) alignment.method <- 'mafft'
  if(missing(custom.specific.values)) custom.specific.values <- c('VDJ_CDR1', 'VDJ_CDR2','VDJ_CDR3','VJ_CDR1','VJ_CDR2','VJ_CDR3')
  if(missing(custom.features)) custom.features <- c('chain', 'region')
  if(missing(custom.weights)) custom.weights <- rep(1, length(custom.specific.values))
  if(missing(custom.normalize.sequence.length)) custom.normalize.sequence.length <- T
  if(missing(custom.commutative)) custom.commutative <- T
  if(missing(plot.results)) plot.results <- T

  #Calls seq/struct alignment and the distance metric

  compute_parallel_distance <- function(steropodon.pdb.list,
                                        distance.metric,
                                        seq.struct.superpose,
                                        struct.superpose,
                                        compare.cores,
                                        cutoff,
                                        max.cycles,
                                        alignment.method,
                                        custom.specific.values,
                                        custom.features,
                                        custom.weights,
                                        custom.normalize.sequence.length,
                                        custom.commutative
                                        ){


    #Custom weighted RMSD score, normalized by sequence length, for various regions of a structure.
    #Modified from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009675
    compute_custom_score <- function(pdb1, pdb2,
                                     specific.values,
                                     features,
                                     weights,
                                     normalize.sequence.length,
                                     commutative
                                   ){

      if(is.null(pdb1$indices)){
        stop('Could not find the superposed c-alpha indices for custom calculation. Ensure seq.struct.superpose is set to T when distance.metric = "custom" or that Steropodon_core was run before if compare.cores is set to T')
      }
      pdb1$atom$pasted <- apply( pdb1$atom[ , features ] , 1 , paste , collapse = "_" )
      pdb2$atom$pasted <- apply( pdb2$atom[ , features ] , 1 , paste , collapse = "_" )

      lengths1 <- unlist(lapply(specific.values, function(x) length(unique(pdb1$atom$resno[pdb1$atom$pasted == x]))))
      lengths2 <- unlist(lapply(specific.values, function(x) length(unique(pdb2$atom$resno[pdb2$atom$pasted == x]))))

      specific.values <- specific.values[lengths1!=0]
      weights <- weights[lengths1!=0]
      lengths1 <- lengths1[lengths1!=0]
      lengths2 <- lengths2[lengths2!=0]

      #weights <- weights / sum(weights)

      indices1 <- lapply(specific.values, function(x) bio3d::atom2xyz(intersect(pdb1$atom$eleno[pdb1$atom$pasted == x], pdb1$indices$atom)))
      indices2 <- lapply(specific.values, function(x) bio3d::atom2xyz(intersect(pdb2$atom$eleno[pdb2$atom$pasted == x], pdb2$indices$atom)))

      dist <- mapply(function(x,y) bio3d::rmsd(pdb1$xyz[x], pdb2$xyz[y]), indices1, indices2) * weights
      dist <- (dist^2) * weights

      if(normalize.sequence.length){
        if(commutative){
          dist <- dist * lengths1 * lengths2
          dist <- sqrt(sum(dist) / sum(weights * lengths1 * lengths2 ))
        }else{
          dist <- dist * lengths1
          dist <- sqrt(sum(dist) / sum(weights * lengths1 ))
        }

      }else{
        dist <- sqrt(sum(dist) / sum(weights))
      }

      return(dist)
    }

    #Computes AlphaFold's modified version of the lDDT score on superposed c-alpha coordinates
    #pdb1 = true, pdb2 = predicted
    #Large cutoff so all residue pairs are considered and score is commutative
    compute_af_lddt <- function(pdb1, pdb2, cutoff = 100, mask = rep(TRUE, length(pdb1$indices$atom))){
      if(is.null(pdb1$indices)){
        stop('Could not find the superposed c-alpha indices for lDDT calculation. Ensure seq.struct.superpose is set to T when distance.metric = "lddt"')
      }

      coords1 <- matrix(pdb1$xyz, ncol = 3, byrow = T)[pdb1$indices$atom,]
      coords2 <- matrix(pdb2$xyz, ncol = 3, byrow = T)[pdb2$indices$atom,]

      dmat1 <- as.matrix(dist(coords1, diag = T, upper = T)) + 1e-10
      dmat2 <- as.matrix(dist(coords2, diag = T, upper = T)) + 1e-10

      mask = as.matrix(mask) * 1
      mask = mask %*% t(mask)
      dists_to_score <- ((dmat1 < cutoff) * 1) * mask * (1 - diag(nrow(dmat1)))

      dist_l1 = abs(dmat1 - dmat2)
      score <- 0.25 * ((dist_l1 < 0.5) * 1 + (dist_l1 < 1) * 1 + (dist_l1 < 2) * 1 + (dist_l1 < 4) * 1)
      norm <- 1 / (1e-10 + sum(dists_to_score))
      score <- norm * (1e-10 + sum(dists_to_score * score))

      return(score)
    }

    #Computes the Template Modeling Score on superposed c-alpha coordinates
    compute_tmscore <- function(pdb1, pdb2){
      if(is.null(pdb1$indices)){
        stop('Could not find the superposed c-alpha indices for Template Modeling Score calculation. Ensure seq.struct.superpose is set to T when distance.metric = "tmscore" or that Steropodon_core was run before if compare.cores is set to T')
      }

      L_target <- length(unique(pdb2$atom$resno))
      L_common <- length(pdb2$indices$atom)

      coords1 <- matrix(pdb1$xyz, ncol = 3, byrow = T)[pdb1$indices$atom,]
      coords2 <- matrix(pdb2$xyz, ncol = 3, byrow = T)[pdb2$indices$atom,]
      d_0 <- 1.24 * ((L_target - 15)^(1/3)) - 1.8
      d_i <- ( sqrt(rowSums((coords1 - coords2) ^ 2)) / d_0 ) ^ 2 + 1
      tmscore <- sum( 1 / d_i) / L_target

      return(tmscore)
    }


    struct_seq_align <- function(mobile,
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
          m <- median(devs)
          std <- sd(devs)

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
          start.num <- regexpr(pattern = paste(c(na.omit(tomatch[1:15])), collapse = ""), text = paste(pdbseq, collapse = ""))[1]

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

    struct_align <- function(pdb1, pdb2, mobile.inds = NULL, fixed.inds = NULL){
      xyz <- suppressWarnings(bio3d::fit.xyz(pdb1$xyz, pdb2$xyz, mobile.inds = mobile.inds, fixed.inds = fixed.inds))
      pdb2$xyz <- xyz
      pdb2$atom[c('x', 'y', 'z')] <- matrix(pdb2$xyz, ncol = 3, byrow = T)
      return(pdb2)
    }

    compute_distance <- function(pdb1, pdb2,
                                 distance.metric,
                                 seq.struct.superpose,
                                 struct.superpose,
                                 compare.cores,
                                 cutoff,
                                 max.cycles,
                                 alignment.method){

      if(compare.cores){
        if(is.null(pdb1$core.resno)){
          stop('Could not find the core indices. Please call Steropodon_find_core before this function!')
        }
        pdb1 <- bio3d::trim(pdb1, inds = bio3d::atom.select(pdb1, resno = pdb1$core.resno))
        pdb2 <- bio3d::trim(pdb2, inds = bio3d::atom.select(pdb2, resno = pdb2$core.resno))
        pdb1 <- bio3d::trim(pdb1, inds = bio3d::atom.select(pdb1, elety = 'CA'))
        pdb2 <- bio3d::trim(pdb2, inds = bio3d::atom.select(pdb2, elety = 'CA'))
      }

      if(seq.struct.superpose){
        aligned <- struct_seq_align(fixed = pdb1,
                                    mobile = pdb2,
                                    cutoff = cutoff,
                                    max.cycles = max.cycles,
                                    aln.method = alignment.method)
        pdb1 <- aligned$fixed
        pdb2 <- aligned$mobile

        pdb1$indices <- pdb1$superposition.indices
        pdb2$indices <- pdb2$superposition.indices
      }

      if(struct.superpose){
        pdb2 <- struct_align(pdb1, pdb2)
      }

      return(distance.metric(pdb1, pdb2))
    }

    if(!is.function(distance.metric)){
      if(distance.metric == 'rmsd'){
        distance.metric <- bio3d::rmsd

      }else if(distance.metric == 'tmscore'){
        distance.metric <- compute_tmscore

      }else if(distance.metric == 'lddt'){
        distance.metric <- compute_af_lddt

      }else if(distance.metric == 'custom'){
        distance.metric <- function(x,y) compute_custom_score(pdb1 = x, pdb2 = y,
                                                              specific.values = custom.specific.values,
                                                              features = custom.features,
                                                              weights = custom.weights,
                                                              normalize.sequence.length = custom.normalize.sequence.length,
                                                              commutative = custom.commutative)

      }else{
        stop('Distance metric not implemented yet!')
      }
    }

    compute_distances <- function(x,y) compute_distance(x, y,
                                                       distance.metric = distance.metric,
                                                       seq.struct.superpose = seq.struct.superpose,
                                                       struct.superpose = struct.superpose,
                                                       compare.cores = compare.cores,
                                                       cutoff = cutoff,
                                                       max.cycles = max.cycles,
                                                       alignment.method = alignment.method)

    #Will modify for non-commutative metrics and for diag != NA
    m <- bigstatsr::FBM(length(steropodon.pdb.list), length(steropodon.pdb.list))
    combs <- combn(1:length(steropodon.pdb.list), 2)
    ident <- rbind(1:length(steropodon.pdb.list), 1:length(steropodon.pdb.list))
    combs <- cbind(combs, ident)
    cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    out <- foreach::foreach(d = 1:ncol(combs), .combine = 'c') %dopar% {

      dist <- compute_distances(steropodon.pdb.list[[combs[1,d]]], steropodon.pdb.list[[combs[2,d]]])

      m[combs[1,d], combs[2,d]] <- dist
      m[combs[2,d], combs[1,d]] <- dist
      NULL
    }
    parallel::stopCluster(cl)

    m <- m[]
    return(m)
  }

  steropodon_list <- unnest_steropodon(steropodon.object)
  seq_ids <- names(steropodon_list)
  pdb_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
  out_matrix <- compute_parallel_distance(pdb_list,
                                           distance.metric = distance.metric,
                                           seq.struct.superpose = seq.struct.superpose,
                                           struct.superpose = struct.superpose,
                                           compare.cores = compare.cores,
                                           cutoff = cutoff,
                                           max.cycles = max.cycles,
                                           alignment.method = alignment.method,
                                           custom.specific.values = custom.specific.values,
                                           custom.features = custom.features,
                                           custom.weights = custom.weights,
                                           custom.normalize.sequence.length = custom.normalize.sequence.length,
                                           custom.commutative = custom.commutative
                                          )

  colnames(out_matrix) <- seq_ids
  rownames(out_matrix) <- seq_ids


  if(plot.results){
    out_matrix <- pheatmap::pheatmap(out_matrix, color = viridis::viridis(10), cluster_rows=F, cluster_cols=F, main=paste0(distance.metric, ' structure distances'))
  }
  return(out_matrix)
}
