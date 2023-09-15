#' Predicts binding interface between the epitope bins and antibodies

#' @description Predicts binding interface between the epitope bins and antibodies. Requires function modelled using ESMFold or AlphaFold, via the AlphaFold_prediction function.
#' @param model string - either 'alphafold' or 'esmfold', depending on the protein folding model used.
#' @param output.folder string - path to the output folder with modelled complexes (antibody - antigen)
#' @param bins vector - residue ids of the epitope bins for the antigen modelled in the antibody-antigen complex.
#' @param cutoff float - minimnum distance (Angstroms) for protein-protein interfaces. Defaults to 5.

#' @return an dataframe of epitope bins and inferred binding interface residues.
#' @export
#' @examples
#' \dontrun{
#' df <- PDBs_predict_bin()
#'}

PDBs_predict_bin <- function(model, output.folder, bins, cutoff=5) {

  #Creating the output object based on the model
  if(model == "alphafold") {
    output.dirs <- base::list.dirs(output.folder, full.names = FALSE, recursive = FALSE)
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

  # Calculating intersections
  preds <- c()
  for(id in names(output)) {
    struc <- output[[paste(id)]]
    binder <- bio3d::atom.select(struc, chain="C")
    mab <- bio3d::atom.select(struc, chain=c("A", "B"))
    bs <- bio3d::binding.site(struc, a.inds=binder, b.inds=mab, cutoff=cutoff)

    intersects <- purrr::map(bins, ~intersect(.x, bs$resno))

    ifelse(lapply(intersects, length)==0, predict <- NA,
           predict <- names(intersects[which.max(lengths(intersects))]))
    preds <- c(preds, predict)
  }

  #Storing predictions on dataframe
  pdb.id <- names(output)
  df <- data.frame(pdb.id)
  df$predicted.bin <- preds

  #Returning dataframe
  return(df)
}
