#' Computes a weighted RMSD score for two (superposed) structures


#'@description Custom weighted RMSD score, normalized by sequence length, for various regions of a structure.
#' The 'specific.values' parameter defines a vector of the VDJ/VJ regions to be considered (will only get the RMSD across these regions, rest will be masked).
#' The 'weights' parameter defines a vector of weights for each region from the 'features' parameter.
#' Modified from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009675

#' @param pdb1 bio3d::pdb object of the first structure.
#' @param pdb2 bio3d::pdb object of the second structure.
#' @param specific.values vector of strings - VDJ or VJ regions to be considered. If only specific regions for a given chain are to be picked, the 'features' parameter should include both 'chain' and 'region'.
#' @param features vector of strings or string - the sequence features to be considered for the distance calculation. 'region' only for FR1/CDR1/etc regions (across all chains), 'chain' for either 'VDJ' or 'VJ' chains (inserted in specific.values).
#' @param weights vector of integer - weight to be assigned in the custom score for each region/chain included in the specific.values parameter.
#' @param normalize.sequence.length bool - if TRUE, will normalize the score per region by that region's number of amino acids.
#' @param commutative bool - if TRUE, will normalize the distance by the length of both sequence lengths of the target and template structures (pdb1 and pdb2) used in the distance calculation.

#' @return the score value (float).
#' @export
#' @examples
#' \dontrun{
#' compute_custom_score(pdb1, pdb2,
#' specific.values = c('VDJ_CDR1', 'VDJ_CDR2', 'VDJ_CDR3'),
#' features = c('region', 'chain'),
#' weights = c(1,1,3),
#' normalize.sequence_length = TRUE)
#'}

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

#' Computes the LDDT score for two (superposed) structures

#'@description Computes the LDDT score for two (superposed) structures. Atoms can be masked from the score calculation using the 'mask' parameter.
#' Modified from https://github.com/google-deepmind/alphafold/blob/7c9114c8423ac9db981d8365168464bab09b3e54/alphafold/model/lddt.py.
#' Coordinates should be superposed beforehand.


#' @param pdb1 bio3d::pdb object of the first structure.
#' @param pdb2 bio3d::pdb object of the second structure.
#' @param cutoff integer - maximum distance for a pair of points to be included.
#' @param mask vector of boolean - which atoms should be masked/not considered in the score calculation (if FALSE).


#' @return the score value (float).
#' @export
#' @examples
#' \dontrun{
#' compute_af_lddt(pdb1, pdb2)
#'}
#'
compute_af_lddt <- function(pdb1, pdb2, cutoff = 15, mask = rep(TRUE, length(pdb1$indices$atom))){
  if(is.null(pdb1$indices)){
    stop('Could not find the superposed c-alpha indices for lDDT calculation. Ensure seq.struct.superpose is set to T when distance.metric = "lddt"')
  }

  coords1 <- matrix(pdb1$xyz, ncol = 3, byrow = T)[pdb1$indices$atom,]
  coords2 <- matrix(pdb2$xyz, ncol = 3, byrow = T)[pdb2$indices$atom,]

  dmat1 <- as.matrix(stats::dist(coords1, diag = T, upper = T)) + 1e-10
  dmat2 <- as.matrix(stats::dist(coords2, diag = T, upper = T)) + 1e-10

  mask = as.matrix(mask) * 1
  mask = mask %*% t(mask)
  dists_to_score <- ((dmat1 < cutoff) * 1) * mask * (1 - diag(nrow(dmat1)))

  dist_l1 = abs(dmat1 - dmat2)
  score <- 0.25 * ((dist_l1 < 0.5) * 1 + (dist_l1 < 1) * 1 + (dist_l1 < 2) * 1 + (dist_l1 < 4) * 1)
  norm <- 1 / (1e-10 + sum(dists_to_score))
  score <- norm * (1e-10 + sum(dists_to_score * score))

  return(score)
}

#' Computes the LDDT score for two (superposed) structures

#' @description Computes the TM-score between two structures.
#' Modified from: https://www.blopig.com/blog/2017/01/tm-score/

#' @param pdb1 bio3d::pdb object of the first structure.
#' @param pdb2 bio3d::pdb object of the second structure.

#' @return the score value (float).
#' @export
#' @examples
#' \dontrun{
#' compute_af_lddt(pdb1, pdb2)
#'}
#'
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
