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
