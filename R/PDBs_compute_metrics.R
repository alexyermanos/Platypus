#' Compute prediction metrics for antibody-antigen complexes

#' @description Computes the prediction metrics for antibody-antigen complexes obtained from ESMFold or AlphaFold, by using Platypus functions such as AlphaFold_prediction().
#' Metrics are computed at the binding interface level (e.g., interface pLDDT).
#' @param model string - either 'alphafold' or 'esmfold', depending on the protein folding model used.
#' @param output.folder string - path to the output folder with modelled complexes (antibody - antigen)
#' @param cutoff float - minimnum distance (Angstroms) for protein-protein interfaces. Defaults to 5.

#' @return a dataframe with metrics per modelled antibody-antigen complex.
#' @export
#' @examples
#' \dontrun{
#' metrics <- PDBs_compute_metrics(model = 'esmfold',
#' output.folder = './predicted_structures',
#' cutoff = 5)
#'}


PDBs_compute_metrics <- function(model, output.folder, cutoff=5) {

  #Creating the output object based on the model
  if(model == "alphafold") {
    output.dirs <- list.dirs(output.folder, full.names=FALSE, recursive=FALSE)
    output <- list()
    sample <- c()
    for(i in 1:length(output.dirs)) {
      output[[length(output)+1]] <- bio3d::read.pdb(paste0(output.folder,"/", output.dirs[i],
                                                    "/", "ranked_0.pdb"))
      sample <- c(sample,paste0(readr::parse_number(output.dirs[i])))
    }
    names(output) <- sample
    output <- output[order(as.integer(names(output)))]
  }

  if(model == "esmfold") {
    files <- list.files(path = output.folder)
    sample <- c()
    output <- list()
    for(i in 1:length(files)) {
      output[[length(output)+1]] <- bio3d::read.pdb(paste0(output.folder, "/", files[i]))
      sample <- c(sample,paste0(readr::parse_number(files[i])))
    }
    names(output) <- sample
    output <- output[order(as.integer(names(output)))]
  }

  #Sub-functions to calculate metrics
  PDB_compute_mmd <- function(pdb, cutoff) {
    try({
      # Initializing interface column to keep track of interacting atoms
      pdb$atom$interface <- rep(0,nrow(pdb$atom))

      # Antigen - C
      antigen <- bio3d::atom.select(pdb, chain="C")
      mab <- bio3d::atom.select(pdb, chain=c("A", "B"))
      bs.C <- bio3d::binding.site(pdb, a.inds=antigen, b.inds=mab, cutoff=cutoff)

      pdb$atom$interface[bs.C$inds$atom] <- 1
      pdb$atom$interface[-bs.C$inds$atom] <- 0

      res.C <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      # mab VDJ - A
      mab.vdj <- bio3d::atom.select(pdb, chain="A")
      antigen <- bio3d::atom.select(pdb, chain="C")
      bs.A <- bio3d::binding.site(pdb, a.inds=mab.vdj, b.inds=antigen, cutoff=cutoff)

      pdb$atom$interface[bs.A$inds$atom] <- 1
      pdb$atom$interface[-bs.A$inds$atom] <- 0

      res.A <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      # mab VJ - B
      mab.vj <- bio3d::atom.select(pdb, chain="B")
      antigen <- bio3d::atom.select(pdb, chain="C")
      bs.B <- bio3d::binding.site(pdb, a.inds=mab.vj, b.inds=antigen, cutoff=cutoff)

      pdb$atom$interface[bs.B$inds$atom] <- 1
      pdb$atom$interface[-bs.B$inds$atom] <- 0

      res.B <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      ##Distance Measure
      b.antibody.A <- pdb$atom  %>% dplyr::filter(chain == "A") %>%
        dplyr::filter(resno %in% res.A) %>% dplyr::filter(elety == "CA")

      b.antibody.B <- pdb$atom  %>% dplyr::filter(chain == "B") %>%
        dplyr::filter(resno %in% res.B) %>% dplyr::filter(elety == "CA")

      b.antibody <- rbind(b.antibody.A,b.antibody.B)

      antibody.xyz <- c()
      for(ii in 1:nrow(b.antibody)){
        antibody.xyz <- c(antibody.xyz,
                          b.antibody$x[[ii]], b.antibody$y[[ii]], b.antibody$z[[ii]])}

      b.antigen <- pdb$atom %>% dplyr::filter(chain == "C") %>%
        dplyr::filter(resno %in% res.C) %>% dplyr::filter(elety == "CA")

      antigen.xyz <- c()
      for(ii in 1:nrow(b.antigen)){
        antigen.xyz <- c(antigen.xyz,
                         b.antigen$x[[ii]],b.antigen$y[[ii]],b.antigen$z[[ii]])}

      dist.mat.list <- bio3d::dist.xyz(antibody.xyz,antigen.xyz)
      m.dist <- apply(dist.mat.list,1, min) %>% mean() %>% round(2)
    })

    if (exists("m.dist")) {return(m.dist)} else {return(NA)}
  }

  PDB_compute_mmd_all <- function(pdb) {
    try({
      # Initializing interface column to keep track of interacting atoms
      pdb$atom$interface <- rep(0,nrow(pdb$atom))

      ##Distance Measure
      b.antibody.A <- pdb$atom  %>% dplyr::filter(chain == "A") %>%
        dplyr::filter(elety == "CA")

      b.antibody.B <- pdb$atom  %>% dplyr::filter(chain == "B") %>%
        dplyr::filter(elety == "CA")

      b.antibody <- rbind(b.antibody.A,b.antibody.B)

      antibody.xyz <- c()
      for(ii in 1:nrow(b.antibody)){
        antibody.xyz <- c(antibody.xyz,
                          b.antibody$x[[ii]], b.antibody$y[[ii]], b.antibody$z[[ii]])}

      b.antigen <- pdb$atom %>% dplyr::filter(chain == "C") %>%
        dplyr::filter(elety == "CA")

      antigen.xyz <- c()
      for(ii in 1:nrow(b.antigen)){
        antigen.xyz <- c(antigen.xyz,
                         b.antigen$x[[ii]],b.antigen$y[[ii]],b.antigen$z[[ii]])}

      dist.mat.list <- bio3d::dist.xyz(antibody.xyz,antigen.xyz)
      m.dist <- apply(dist.mat.list,1, min) %>% mean() %>% round(2)
    })

    if (exists("m.dist")) {return(m.dist)} else {return(NA)}
  }

  PDB_compute_plddt <- function(pdb, cutoff) {
    try({
      # Initializing interface column to keep track of interacting atoms
      pdb$atom$interface <- rep(0,nrow(pdb$atom))

      # Antigen - C
      antigen <- bio3d::atom.select(pdb, chain="C")
      mab <- bio3d::atom.select(pdb, chain=c("A", "B"))
      bs.C <- bio3d::binding.site(pdb, a.inds=antigen, b.inds=mab, cutoff=cutoff)

      pdb$atom$interface[bs.C$inds$atom] <- 1
      pdb$atom$interface[-bs.C$inds$atom] <- 0

      res.C <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      # mab VDJ - A
      mab.vdj <- bio3d::atom.select(pdb, chain="A")
      antigen <- bio3d::atom.select(pdb, chain="C")
      bs.A <- bio3d::binding.site(pdb, a.inds=mab.vdj, b.inds=antigen, cutoff=cutoff)

      pdb$atom$interface[bs.A$inds$atom] <- 1
      pdb$atom$interface[-bs.A$inds$atom] <- 0

      res.A <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      # mab VJ - B
      mab.vj <- bio3d::atom.select(pdb, chain="B")
      antigen <- bio3d::atom.select(pdb, chain="C")
      bs.B <- bio3d::binding.site(pdb, a.inds=mab.vj, b.inds=antigen, cutoff=cutoff)

      pdb$atom$interface[bs.B$inds$atom] <- 1
      pdb$atom$interface[-bs.B$inds$atom] <- 0

      res.B <- pdb$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()

      ##pLDDT
      m.A <- pdb$atom %>% dplyr::filter(chain == "A") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.A) %>% .$b %>% mean()

      m.B <- pdb$atom %>% dplyr::filter(chain == "B") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.B) %>% .$b %>% mean()

      m.C <- pdb$atom %>% dplyr::filter(chain == "C") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.C) %>% .$b %>% mean()

      m.bind <- mean(m.A,m.B,m.C) %>% round(.,2)
    })

    if (exists("m.bind")) {return(m.bind)} else {return(NA)}
  }

  PDB_compute_plddt_all <- function(pdb) {
    try({
      # Initializing interface column to keep track of interacting atoms
      pdb$atom$interface <- rep(0,nrow(pdb$atom))

      # Residues
      res.C <- pdb$atom %>% .$resno %>% unique()
      res.A <- pdb$atom %>% .$resno %>% unique()
      res.B <- pdb$atom %>% .$resno %>% unique()

      ##pLDDT
      m.A <- pdb$atom %>% dplyr::filter(chain == "A") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.A) %>% .$b %>% mean()

      m.B <- pdb$atom %>% dplyr::filter(chain == "B") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.B) %>% .$b %>% mean()

      m.C <- pdb$atom %>% dplyr::filter(chain == "C") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res.C) %>% .$b %>% mean()

      m.bind <- mean(m.A,m.B,m.C) %>% round(.,2)
    })

    if (exists("m.bind")) {return(m.bind)} else {return(NA)}
  }

  #Calculating mmd and plddt, both entire structures and metrics based on cutoff
  output.mmd.cutoff <- c()
  output.mmd.all <- c()
  output.plddt.cutoff <- c()
  output.plddt.all <- c()

  for (pdb in output) {
    mmd.cutoff <- PDB_compute_mmd(pdb, cutoff)
    output.mmd.cutoff <- c(output.mmd.cutoff, mmd.cutoff)

    mmd.all <- PDB_compute_mmd_all(pdb)
    output.mmd.all <- c(output.mmd.all, mmd.all)

    plddt.cutoff <- PDB_compute_plddt(pdb, cutoff)
    output.plddt.cutoff <- c(output.plddt.cutoff, plddt.cutoff)

    plddt.all <- PDB_compute_plddt_all(pdb)
    output.plddt.all <- c(output.plddt.all, plddt.all)
  }

  #Storing metrics on dataframe
  pdb.id <- names(output)
  df <- data.frame(pdb.id)
  df$mmd <- output.mmd.all
  df$mmd.cutoff <- output.mmd.cutoff
  df$plddt <- output.plddt.all
  df$plddt.cutoff <- output.plddt.cutoff

  #Returning dataframe
  return(df)
}
