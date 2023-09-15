#' Splits a bio3d pdb object given a grouping and a specific value

#' @description Computes the root mean square deviation (RMSD) between antibody (-antigen) structures predicted with AlphaFold and ESMFold. Outputs a dataframe of RMSD per predicted structure pairs.
#' @param alphafold.output.folder string - path to the directory containing structures predicted with AlphaFold.
#' @param esmfold.output.folder string - path to the directory containing structures predicted with ESMFold

#' @return a dataframe of RMSD per pair of predicted structures (using AlphaFold and ESMFold).
#' @export
#' @examples
#' \dontrun{
#' rmsds <- AF_ESM_compute_RMSD(alphafold.output.folder = './af_structures',
#' esmfold.output.folder = './esm_structures')
#'}


AF_ESM_compute_RMSD <- function(alphafold.output.folder, esmfold.output.folder) {
  #First creating the alphafold_output object
  output.dirs <- list.dirs(alphafold.output.folder, full.names = FALSE, recursive = FALSE)
  alphafold.output <- list()
  sample <- c()
  for(i in 1:length(output.dirs)) {
    alphafold.output[[length(alphafold.output)+1]] <- biop3d::read.pdb(paste0(
      alphafold.output.folder,"/", output.dirs[i],"/", "ranked_0.pdb"))
    sample <- c(sample,paste0(parse_number(output.dirs[i])))
  }
  names(alphafold.output) <- sample
  alphafold.output <- alphafold.output[order(as.integer(names(alphafold.output)))]

  #Then creating the esmfold.output object
  files <- list.files(path = esmfold.output.folder)
  sample <- c()
  esmfold.output <- list()
  for(i in 1:length(files)) {
    esmfold.output[[length(esmfold.output)+1]] <- bio3d::read.pdb(paste0(esmfold.output.folder,
                                                                  "/", files[i]))
    sample <- c(sample,paste0(parse_number(files[i])))
  }
  names(esmfold.output) <- sample
  esmfold.output <- esmfold.output[order(as.integer(names(esmfold.output)))]

  #Calculating pairwise RMSD
  rmsd <- c()
  for (i in seq(length(alphafold.output))) {
    alphafold <- alphafold.output[[paste(i)]]
    esmfold <- esmfold.output[[paste(i)]]
    align <- struct.aln(alphafold, esmfold, exefile='msa')

    rmsd <- c(rmsd, bio3d::rmsd(alphafold, esmfold, fit=TRUE,
                         a.inds=align$a.inds$xyz, b.inds=align$b.inds$xyz))
  }

  #Storing RMSD on dataframe
  pdb.id <- names(alphafold.output)
  df <- data.frame(pdb.id)
  df$rmsd <- rmsd

  #Returning dataframe
  return(df)
}
