#'Analysis of antibody structures
#'
#'@description VDJ_structure_analysis is a Platypus function for the analysis of 3d structures of antibodies.
#'The function is designed in a way to be a follow up of the alphafold_prediction function of the Platypus package.
#'The input of this function is a list that has the VDJ MIXCR output in its first element and a list of structures in its second element.
#'This is the output format of the alphafold_prediction function. The function can also be used to visualize structures from PDB files directly.
#'
#'The function can visualize the structures of proteins and is especially designed for visualization of antibodies, where the
#'framework and CDR regions are automatically annotated. If the antibody was predicted together with a antigen, the function can visualize
#'binding interaction and detect the binding site residues. Furthermore, it can determine binding site metrics as average distance and model accuracy.
#'The function has a variety of arguments to create the desired visualization.
#'
#'@param VDJ.structure The VDJ.structure is a list object with the VDJ MIXCR out data frame as its first element containing the germline reference sequences. The second list element contains a list for every predicted structures with the top ranked predictions in pdb format read in by the bio3d:read_pdb() function.
#'@param cells.to.vis  In the cells.to.vis argument, the barcodes of the cells that should be analyzed can be specified. It can be a single barcode or a list of barcodes. If all the elements of the VDJ.structure object shall be included, the cells.to.vis argument can be set to "ALL" (cells.to.vis = "ALL).
#'@param rank AlphaFold predicts multiple models for a given structure which are then ranked according to model's confidence with ranked_0.pdb as the best fit. By default the function uses the ranked_0 model but if a different rank is desired this can be specified in teh rank argument (rank = "ranked_5.pdb").
#'@param rank.list The function can also be used to analyze multiple ranked models per structure. This can be done by specifying the rank.list argument. rank.list = list("ranked_0.pdb", "ranked_1.pdb", "ranked_2.pdb") would use the top three most confident models per structure that are specified in the cells.to.vis argument.
#'@param overlay Multiple antibody structures are by default visualized separately. In case an overlay is desired this can be done by setting overlay to TRUE (overlay = T).
#'@param PDB.file Instead of visualizing the structure form the VDJ.structure object, the function can accept a path to a pdb.file as well.
#'@param spin.speed Protein animations can spin around an internal axis. By default the spin speed is set to 0. To start rotation the spin.speed argument can be set to the desired rotation speed.
#'@param VDJ.anno When using a VDJ structure object as an input, the regions of the antibody are by default annotated based on the MIXCR columns. If annotation is not desired the VDJ.anno argument can be set to FALSE.
#'@param color.cdr3 The color of the cdr3 can be changed by specifying the color.cdr3 argument with a the desired color in HEX format (color.cdr3 = "#eb4034")
#'@param color.cdr2 The color of the cdr2 can be changed by specifying the color.cdr2 argument with a the desired color in HEX format (color.cdr2 = "#eb4034")
#'@param color.cdr1 The color of the cdr1 can be changed by specifying the color.cdr1 argument with a the desired color in HEX format (color.cdr1 = "#eb4034")
#'@param color.frameworks The color of the frameworks can be changed by specifying the color.frameworks argument with a the desired color in HEX format (color.frameworks = "#eb4034")
#'@param label By default the annotated regions are labeled. The label can be disabled by setting the label argument to FASLE.
#'@param label.size The label size can be adjusted by specifying the label.size argument. By default the label size is set to 12.
#'@param bk.opac The label background opacity can be defined with the bk.opac argument. The default opacity is 0.8.
#'@param font.opac The opacity of the label's font can be set by specifying the font.opac argument. The default opacity is 1.
#'@param font.col The color of the font can be set by the font.col argument. It has to be in a HEX format.
#'@param anno.seq By the anno.seq argument any residues of the structure can be annotated. Every domain of the structure is handled as a separate chain, named alphabetically from A-Z, according to the order in the FASTA file.  So for an anybody the HC is Chain A, the LC Chain B and the Antigen Chain C. The format of the anno.seq argument is a vector with the index of the starting residue as its first element, the index of the end residue as second element, the chain in the third element and the color in HEX format as the fourth element. Optionally a label can be added by specifying the text in the fifth element anno.seq = c(4,12,"C","#9900ff","label"). For annotating multiple sequences at once, a list of vectors can be used anno.seq = list(c(4,12,"C","#9900ff","label1"),c(4,12,"B","#9350ff","label2"),...)
#'@param color.molecule The color of the non annotated residues can be set by the color.molecule argument in HEX format.
#'@param color.sheets The color of beta sheets can be set by specifying the color.sheets argument with the desired color in HEX format.
#'@param color.helix The color of alpha helices can be set by specifying the color.helix argument with the desired color in HEX format.
#'@param angle.x The molecule can be rotated around the x-axis by setting the angle.x argument in degree.
#'@param angle.y The molecule can be rotated around the y-axis by setting the angle.y argument in degree.
#'@param angle.z The molecule can be rotated around the z-axis by setting the angle.z argument in degree.
#'@param r3dmol.code Additional r3dmol code to modify appearance of the final structure. Defaults to "". Character input is run via eval(parse(x))
#'@param plddt.plot AlphaFold has a pLDDT confidence score for every residue of the structure. An additional plot of the structure colored according to the pLDDT score will be returned by setting the plddt.plot to TRUE.
#'@param structure.plot By default the structure is visualized. However, if the function is only used to get binding site metrics, the structure.plot argument can be set to FALSE. Like this nos structures are plotted.
#'@param antigen.interaction If the antibody is predicted together with the antigen the antigen.interaction argument can be set to TRUE in order to get some binding site metrics. The function will determine the binding site residues based on the bio3d::binding.site() function as well as the average minimal distance of every binding site residue from the antibody to the antigen. Furthermore, the mean confidence (pLDDT) of the binding site residues is calculated. All the results are summarized in a data frame.
#'@param binding.site.cutoff Cutoff for the bio3d::binding.site() function. Default is 5.
#'@param dist.mat If set to TRUE, the binding residues distance matrix for every structure will be returned as a list. Default = FALSE, so only the minimal distances in the summary data frame is returned.
#'@param SASA if set to TRUE the Solvent Accessible Surface Area will be calculated for every Structure
#'@param hydrophobicity If set to TRUE, the per residue hydrophobicity will be calculated.
#'@param charge If set to TRUE, the per residue charge will be calculated.
#'@param metrics.plot If set to TRUE, a structure plot colored according to the metrics will be shown.
#'@param BindingResidues.plot If the antigen.interaction is enabled, the binding site residues can be visualized on the structure by setting the BindingResidues.plot argument to TRUE.
#'@param binding.residue.barplot The binding.residue.bar plot can be set to TRUE to get a bar plot visualizing to which regions of the antibody the binding site residues belong too. A bar plot is produced for every structure separately as well as one summary bar plot over all analyzed structures.
#'@param binding.residue.barplot.style There are two styles available for the binding site residue bar plot. By default all regions are listed separately. By setting the binding.residue.barplot.style to "FR_CDR", only framework and CDR is distinguished.
#' @return ADD DESCRIPTION OF RETURN VALUE HERE
#' @export
#' @examples
#' \dontrun{
#'
#' ADD EXAMPLE CODE HERE
#'
#' }



VDJ_structure_analysis <- function(VDJ.structure,
                               cells.to.vis,
                               rank,
                               rank.list,
                               overlay,
                               PDB.file,
                               spin.speed,
                               color.frameworks,
                               color.cdr3,
                               color.cdr2,
                               color.cdr1,
                               VDJ.anno,
                               label,
                               label.size,
                               bk.opac,
                               font.opac,
                               font.col,
                               anno.seq,
                               color.molecule,
                               color.sheets,
                               color.helix,
                               angle.x,
                               angle.y,
                               angle.z,
                               r3dmol.code,
                               plddt.plot,
                               structure.plot,
                               antigen.interaction,
                               binding.site.cutoff,
                               dist.mat,
                               SASA,
                               hydrophobicity,
                               charge,
                               metrics.plot,
                               BindingResidues.plot,
                               binding.residue.barplot,
                               binding.residue.barplot.style
                               ) {


  if(missing(VDJ.structure)) {VDJ.structure.input <- F}
  else{VDJ.structure.input <- T}

  if(missing(rank) & missing(rank.list)) {rank <- "ranked_0.pdb"}
  if(missing(rank.list)) {rank.list.input <- F
  rank.list <- F}
  else{rank.list.input <- T}
  if(missing(overlay)) {overlay <- F}
  if(missing(PDB.file)) {PDB.file.input <- F}
  else{PDB.file.input <- T}

  if(missing(spin.speed)) {spin.speed <- 0}
  if(missing(label)) {label <- T}
  if(missing(VDJ.anno)) {VDJ.anno <- T}
  if(missing(label.size)) {label.size <- 12}
  if(missing(bk.opac)) {bk.opac <- 0.8}
  if(missing(font.opac)) {font.opac <- 1}
  if(missing(font.col)) {font.col <- "balck"}
  if(missing(anno.seq)) {anno.seq.input <- F}
  else(anno.seq.input <- T)

  if(missing(color.cdr3)) {color.cdr3 <- "#ff00ff"}
  if(missing(color.cdr2)) {color.cdr2 <- "#00cc00"}
  if(missing(color.cdr1)) {color.cdr1 <- "#66ccff"}
  if(missing(color.frameworks)) {color.frameworks <- "#ff9900"}
  if(missing(color.molecule)) {color.molecule <- "#00cc96"}
  if(missing(color.helix)) {color.helix <- "#ff7f0e"}
  if(missing(color.sheets)) {color.sheets <- "#636efa"}
  if(missing(angle.x)) {angle.x <- 0}
  if(missing(angle.y)) {angle.y <- 90}
  if(missing(angle.z)) {angle.z <- 0}
  if(missing(r3dmol.code)) {r3dmol.code <- ''}
  if(missing(plddt.plot)) {plddt.plot <- F}
  if(missing(structure.plot)) {structure.plot <- T}
  if(missing(antigen.interaction)) {antigen.interaction <- F}
  if(missing(binding.site.cutoff)) {binding.site.cutoff <- 5}
  if(missing(dist.mat)) {dist.mat <- F}
  if(missing(SASA)) {SASA <- F}
  if(missing(hydrophobicity)) {hydrophobicity <- F}
  if(missing(charge)) {charge <- F}
  if(missing(metrics.plot)) {metrics.plot <- F}
  if(missing(BindingResidues.plot)) {BindingResidues.plot <- F}
  if(missing(binding.residue.barplot)) {binding.residue.barplot <- F}
  if(missing(binding.residue.barplot.style)) {binding.residue.barplot.style <- "all_regions"}


  if(!label) {
    bk.opac <- 0
    font.opac <- 0
  }

  vis.structure <- NULL
  models <- NULL
  anno <- NULL
  out.plot <- NULL
  . <- NULL
  barcode <- NULL
  interface <- NULL
  chain <- NULL
  resno <- NULL
  elety <- NULL
  region <- NULL
  Nr_bind_res <- NULL
  Sum_nr_bind_res_pct <- NULL

  ##Read in the structure(s) to be analyzed / visualized

  if(VDJ.structure.input){
    #Check if cells.to.vis is specified
    if(missing(cells.to.vis)) stop("Please specify the entries of the VDJ.structure object to be visualiezed in the 'cells.to.vis' argument")

    cells.to.vis <- as.list(cells.to.vis)
    if(cells.to.vis[[1]] == "ALL"){
      cells.to.vis <- VDJ.structure[[1]]$barcode %>% as.list()

      ## Check if there are structures for all entries and otherwise just use the ones with structures
      empty_entries <- 1
      for(i in VDJ.structure[[2]]){empty_entries <- empty_entries * !purrr::is_empty(i)}

      if(!empty_entries){
        message("Not for every entry a structure is defined! All the defined structures are analysed.")
        idx <- c()
        ii <- 1
        for(i in VDJ.structure[[2]]){
          if(!purrr::is_empty(i)){idx <- c(idx,ii)}
          ii <- ii + 1
        }
        cells.to.vis <- names(VDJ.structure[[2]])[idx] %>% as.list()
      }
    }



    if(rank.list.input){

      vis.structure <- list()
      for(i in 1:length(cells.to.vis)){
        for(n in 1:length(rank.list)){
          vis.structure[[length(vis.structure)+1]] <- VDJ.structure[[2]][[cells.to.vis[[i]]]][[rank.list[[n]]]]
        }
      }

      rank <- rep(as.vector(rank.list),length(cells.to.vis))

    }
    else{

      #Set rank
      rank <- as.list(rank)
      if(length(rank) == 1 & length(cells.to.vis) > 1){
        rank <- rep(rank[[1]],length(cells.to.vis))
      }


      #Define the structure that should be visualized
      vis.structure <- list()
      for(i in 1:length(cells.to.vis)){
        vis.structure[[i]] <- VDJ.structure[[2]][[cells.to.vis[[i]]]][[rank[[i]]]]
      }

    }
  }

  if(PDB.file.input){
    #Check if input is a list of PDB files or just one PDB file
    if(length(names(PDB.file) == 4)){
      PDB.file <- list(PDB.file)
    }

    #Define the structure that should be visualized
    vis.structure <- list()
    names.model <- list()
    for(i in 1:length(PDB.file)){
      vis.structure[[i]] <- PDB.file[[i]]

      names.model[[i]] <- PDB.file[[i]]$call %>% toString()
    }
  }


  #Define a list of return elements
  return.list <- list()
  BindingResidues.list <- list()
  dist.mat.list <- list()





  if(antigen.interaction){

    metrics.list <- list()


    if(binding.residue.barplot){

      bin.res.list <- list()
      #Define the location of the given elements
      VDJ_FR1 <- list()
      VDJ_CDR1 <- list()
      VDJ_FR2 <- list()
      VDJ_CDR2 <- list()
      VDJ_FR3 <- list()
      VDJ_CDR3 <- list()
      VDJ_FR4 <- list()

      VJ_FR1 <- list()
      VJ_CDR1 <- list()
      VJ_FR2 <- list()
      VJ_CDR2 <- list()
      VJ_FR3 <- list()
      VJ_CDR3 <- list()
      VJ_FR4 <- list()


      for(i in 1:length(cells.to.vis)){
        for(n in 1:(length(vis.structure)/length(cells.to.vis))){
          sequnce.to.vis <- VDJ.structure[[1]] %>% dplyr::filter(.,barcode == cells.to.vis[[i]])
          VDJ_aa <- sequnce.to.vis$VDJ_aa_mixcr

          VDJ_FR1[[length(VDJ_FR1)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR1)
          VDJ_CDR1[[length(VDJ_CDR1)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR1)
          VDJ_FR2[[length(VDJ_FR2)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR2)
          VDJ_CDR2[[length(VDJ_CDR2)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR2)
          VDJ_FR3[[length(VDJ_FR3)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR3)
          VDJ_CDR3[[length(VDJ_CDR3)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR3)
          VDJ_FR4[[length(VDJ_FR4)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR4)

          VJ_aa <- sequnce.to.vis$VJ_aa_mixcr

          VJ_FR1[[length(VJ_FR1)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR1)
          VJ_CDR1[[length(VJ_CDR1)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR1)
          VJ_FR2[[length(VJ_FR2)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR2)
          VJ_CDR2[[length(VJ_CDR2)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR2)
          VJ_FR3[[length(VJ_FR3)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR3)
          VJ_CDR3[[length(VJ_CDR3)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR3)
          VJ_FR4[[length(VJ_FR4)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR4)

        }
      }
     }





    out.list <- list()
    for(i in 1:length(vis.structure)){
      ## Calculate binding site metrics
      ############################################################################
      struc <- vis.structure[[i]]
      struc$atom$interface <- rep(0,nrow(struc$atom))

      # atom selection
      a.inds <- bio3d::atom.select(struc, chain="C")
      b.inds <- bio3d::atom.select(struc, chain=c("A", "B"))

      # identify interface residues
      bs <- bio3d::binding.site(struc, a.inds=a.inds, b.inds=b.inds, cutoff = binding.site.cutoff)

      # use interface column to store interface in PDB file
      struc$atom$interface[ bs$inds$atom ] <- 1
      struc$atom$interface[ -bs$inds$atom ] <- 0

      res_C <- struc$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()
      ###################################################################

      struc <- vis.structure[[i]]
      struc$atom$interface <- rep(0,nrow(struc$atom))

      # atom selection
      a.inds <- bio3d::atom.select(struc, chain="A")
      b.inds <- bio3d::atom.select(struc, chain="C")

      # identify interface residues
      bs <- bio3d::binding.site(struc, a.inds=a.inds, b.inds=b.inds, cutoff = binding.site.cutoff)

      # use interface column to store interface in PDB file
      struc$atom$interface[ bs$inds$atom ] <- 1
      struc$atom$interface[ -bs$inds$atom ] <- 0

      res_A <- struc$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()
      ###################################################################

      struc <- vis.structure[[i]]
      struc$atom$interface <- rep(0,nrow(struc$atom))

      # atom selection
      a.inds <- bio3d::atom.select(struc, chain="B")
      b.inds <- bio3d::atom.select(struc, chain="C")

      # identify interface residues
      bs <- bio3d::binding.site(struc, a.inds=a.inds, b.inds=b.inds, cutoff = binding.site.cutoff)

      # use interface column to store interface in PDB file
      struc$atom$interface[ bs$inds$atom ] <- 1
      struc$atom$interface[ -bs$inds$atom ] <- 0

      res_B <- struc$atom %>% dplyr::filter(.,interface == 1) %>% .$resno %>% unique()


      m_A <- struc$atom %>% dplyr::filter(chain == "A") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res_A) %>% .$b %>% mean() %>% round(.,2)

      m_B <- struc$atom %>% dplyr::filter(chain == "B") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res_B) %>% .$b %>% mean() %>% round(.,2)

      m_C <- struc$atom %>% dplyr::filter(chain == "C") %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(resno %in% res_C) %>% .$b %>% mean() %>% round(.,2)

      m_rest <- struc$atom %>% dplyr::distinct(resno, .keep_all = T) %>% dplyr::filter(!resno %in% c(res_A,res_B,res_C)) %>% .$b %>% mean() %>% round(.,2)



      ##Distance Measure
      b_antibody_A <- struc$atom  %>% dplyr::filter(chain == "A") %>% dplyr::filter(resno %in% res_A) %>% dplyr::filter(elety == "CA")

      b_antibody_B <- struc$atom  %>% dplyr::filter(chain == "B") %>% dplyr::filter(resno %in% res_B) %>% dplyr::filter(elety == "CA")

      b_antibody <- rbind(b_antibody_A,b_antibody_B)

      antibody_xyz <- c()
      for(ii in 1:nrow(b_antibody)){
        antibody_xyz <- c(antibody_xyz,b_antibody$x[[ii]],b_antibody$y[[ii]],b_antibody$z[[ii]])
      }

      b_antigen <- struc$atom %>% dplyr::filter(chain == "C") %>% dplyr::filter(resno %in% res_C) %>% dplyr::filter(elety == "CA")
      antigen_xyz <- c()
      for(ii in 1:nrow(b_antigen)){
        antigen_xyz <- c(antigen_xyz,b_antigen$x[[ii]],b_antigen$y[[ii]],b_antigen$z[[ii]])
      }

      dist.mat.list[[i]] <- bio3d::dist.xyz(antibody_xyz,antigen_xyz)
      m_dist <- apply(dist.mat.list[[i]],1, min) %>% mean() %>% round(2)




      ##Metrics
      metrics_df <- data.frame(struc$atom)
      metrics_df$charge <- NULL


      ##Solvent Accessible Surface Area of Binding site residues
      if(SASA){

        SASA_df <- vanddraabe::FreeSASA.diff(struc$atom)
        SASA_exp_df <- SASA_df$uniq.atom.ids %>% stringr::str_split("_") %>% do.call(rbind, .) %>% data.frame()
        names(SASA_exp_df) <- c("resid", "resno", "chain", "elety","eleno")
        SASA_df <- cbind(SASA_df,SASA_exp_df)
        SASA_df$resno <- as.integer(SASA_df$resno)
        SASA_df$eleno <- as.integer(SASA_df$eleno)
        SASA_df$uniq.atom.ids <- NULL

        metrics_df <- dplyr::left_join(metrics_df,SASA_df,by = c("resno","chain","elety","eleno","resid"))



        m_SASA_A <- SASA_df %>% dplyr::filter(chain == "A" & resno %in% res_A) %>% .$SASA.prot %>% mean()
        m_SASA_B <- SASA_df %>% dplyr::filter(chain == "B" & resno %in% res_B) %>% .$SASA.prot %>% mean()
        m_SASA_C <- SASA_df %>% dplyr::filter(chain == "C" & resno %in% res_C) %>% .$SASA.prot %>% mean()
        m_SASA <- SASA_df$SASA.prot %>% mean()

        SASA_out_list <- list(
          "mean_SASA_HC" = m_SASA_A,
          "mean_SASA_LC" = m_SASA_B,
          "mean_SASA_antigen" = m_SASA_C,
          "mean_SASA_overall" = m_SASA
        )

      }




      #Hydrophobicity
      if(hydrophobicity){
        df_hydph <- data.frame()

        for(CHAIN in unique(struc$atom$chain)){
          df_chain <- struc$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
          df_chain$hydrophobicity <- struc$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::hydrophobicity()

          df_hydph <- rbind(df_hydph,df_chain)
        }

        metrics_df <- dplyr::left_join(metrics_df,df_hydph,by = c("resno","chain"))


        m_Hydph_A <- df_hydph %>% dplyr::filter(chain == "A" & resno %in% res_A) %>% .$hydrophobicity %>% mean()
        m_Hydph_B <- df_hydph %>% dplyr::filter(chain == "B" & resno %in% res_B) %>% .$hydrophobicity %>% mean()
        m_Hydph_C <- df_hydph %>% dplyr::filter(chain == "C" & resno %in% res_C) %>% .$hydrophobicity %>% mean()
        m_Hydph <- df_hydph$hydrophobicity %>% mean()

        Hydph_out_list <- list(
          "mean_Hydrophobicity_HC" = m_Hydph_A,
          "mean_Hydrophobicity_LC" = m_Hydph_B,
          "mean_Hydrophobicity_antigen" = m_Hydph_C,
          "mean_Hydrophobicity_overall" = m_Hydph
        )
      }


      #Charge
      if(charge){
        df_charge <- dplyr::data_frame()

        for(CHAIN in unique(struc$atom$chain)){
          df_chain <- struc$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
          df_chain$charge <- struc$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::charge()

          df_charge <- rbind(df_charge,df_chain)
        }

        metrics_df <- dplyr::left_join(metrics_df,df_charge,by = c("resno","chain"))

        m_charge_A <- df_charge %>% dplyr::filter(chain == "A" & resno %in% res_A) %>% .$charge %>% mean()
        m_charge_B <- df_charge %>% dplyr::filter(chain == "B" & resno %in% res_B) %>% .$charge %>% mean()
        m_charge_C <- df_charge %>% dplyr::filter(chain == "C" & resno %in% res_C) %>% .$charge %>% mean()
        m_charge <- df_charge$charge %>% mean()

        charge_out_list <- list(
          "mean_charge_HC" = m_charge_A,
          "mean_charge_LC" = m_charge_B,
          "mean_charge_antigen" = m_charge_C,
          "mean_charge_overall" = m_charge
        )

      }

      #Add the metrics data frame to the list
      metrics.list[[length(metrics.list)+1]] <- metrics_df











      ## Write output summary
      if(VDJ.structure.input){
        out.list[[i]] <- list(
             "barcode"=cells.to.vis[[ceiling(i/length(rank.list))]],
             "rank" = rank[[i]] ,
             "Mean_plddt_non_bind_resi:" = m_rest,
             "Mean_plddt_bind_resi_HC:" = m_A,
             "Mean_plddt_bind_resi_LC:" = m_B,
             "Mean_plddt_bind_resi_antigen:" = m_C,
             "Mean_dist_bind_resi" = m_dist,
             "Bind_res_HC" = paste0(res_A,collapse = ","),
             "Bind_res_LC" = paste0(res_B,collapse = ","),
             "Bind_res_antigen" = paste0(res_C,collapse = ",")
        )

      }

      if(PDB.file.input){
        out.list[[i]] <- list(
          "barcode"=names.model[[i]],
          "Mean_plddt_non_bind_resi" = m_rest,
          "Mean_plddt_bind_resi_chain_A" = m_A,
          "Mean_plddt_bind_resi_chain_B" = m_B,
          "Mean_plddt_bind_resi_chain_C" = m_C,
          "Mean_dist_bind_resi" = m_dist,
          "Bind_res_A" = paste0(res_A,collapse = ","),
          "Bind_res_B" = paste0(res_B,collapse = ","),
          "Bind_res_C" = paste0(res_C,collapse = ",")
        )
      }

      if(SASA){out.list[[i]] <- append(out.list[[i]],SASA_out_list)}
      if(hydrophobicity){out.list[[i]] <- append(out.list[[i]],Hydph_out_list)}
      if(charge){out.list[[i]] <- append(out.list[[i]],charge_out_list)}



      if(BindingResidues.plot){

        # Set up the initial viewer
        BindingResidues.list[[i]] <- r3dmol::r3dmol(
          viewer_spec = r3dmol::m_viewer_spec(
            cartoonQuality = 10,
            lowerZoomLimit = 50,
            upperZoomLimit = 350
          )

        ) %>%
          # Add model to scene
          r3dmol::m_add_model(data =  r3dmol::m_bio3d(struc)) %>%
          # Zoom to encompass the whole scene
          r3dmol::m_zoom_to() %>%
          # Set style of structures
          r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "#b3b3b3")) %>%

          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = res_A, chain = "A"), # Style alpha helix
            style = r3dmol::m_style_cartoon(color = "#ff00ff")
          ) %>%

          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = res_B, chain = "B"), # Style alpha helix
            style = r3dmol::m_style_cartoon(color = "#ff00ff")
          ) %>%

          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = res_C, chain = "C"), # Style alpha helix
            style = r3dmol::m_style_cartoon(color = "#0000ff")
          ) %>%



          # Rotate the scene by given angle on given axis
          r3dmol::m_rotate(angle = 90, axis = "y") %>%

        # Animate the scene by spinning it
         r3dmol::m_spin(speed = spin.speed)

      }


      if(binding.residue.barplot){

        VDJ_FR1_c <- 0
        VDJ_CDR1_c <- 0
        VDJ_FR2_c <- 0
        VDJ_CDR2_c <- 0
        VDJ_FR3_c <- 0
        VDJ_CDR3_c <- 0
        VDJ_FR4_c <- 0
        for(n in 1:length(res_A)){
          if(res_A[[n]] %in% VDJ_FR1[[i]][[1]]:VDJ_FR1[[i]][[2]]){VDJ_FR1_c <- VDJ_FR1_c + 1}
          if(res_A[[n]] %in% VDJ_CDR1[[i]][[1]]:VDJ_CDR1[[i]][[2]]){VDJ_CDR1_c <- VDJ_CDR1_c + 1}
          if(res_A[[n]] %in% VDJ_FR2[[i]][[1]]:VDJ_FR2[[i]][[2]]){VDJ_FR2_c <- VDJ_FR2_c + 1}
          if(res_A[[n]] %in% VDJ_CDR2[[i]][[1]]:VDJ_CDR2[[i]][[2]]){VDJ_CDR2_c <- VDJ_CDR2_c + 1}
          if(res_A[[n]] %in% VDJ_FR3[[i]][[1]]:VDJ_FR3[[i]][[2]]){VDJ_FR3_c <- VDJ_FR3_c + 1}
          if(res_A[[n]] %in% VDJ_CDR3[[i]][[1]]:VDJ_CDR3[[i]][[2]]){VDJ_CDR3_c <- VDJ_CDR3_c + 1}
          if(res_A[[n]] %in% VDJ_FR4[[i]][[1]]:VDJ_FR4[[i]][[2]]){VDJ_FR4_c <- VDJ_FR4_c + 1}
        }


        VJ_FR1_c <- 0
        VJ_CDR1_c <- 0
        VJ_FR2_c <- 0
        VJ_CDR2_c <- 0
        VJ_FR3_c <- 0
        VJ_CDR3_c <- 0
        VJ_FR4_c <- 0
        for(n in 1:length(res_B)){
          if(res_B[[n]] %in% VJ_FR1[[i]][[1]]:VJ_FR1[[i]][[2]]){VJ_FR1_c <- VJ_FR1_c + 1}
          if(res_B[[n]] %in% VJ_CDR1[[i]][[1]]:VJ_CDR1[[i]][[2]]){VJ_CDR1_c <- VJ_CDR1_c + 1}
          if(res_B[[n]] %in% VJ_FR2[[i]][[1]]:VJ_FR2[[i]][[2]]){VJ_FR2_c <- VJ_FR2_c + 1}
          if(res_B[[n]] %in% VJ_CDR2[[i]][[1]]:VJ_CDR2[[i]][[2]]){VJ_CDR2_c <- VJ_CDR2_c + 1}
          if(res_B[[n]] %in% VJ_FR3[[i]][[1]]:VJ_FR3[[i]][[2]]){VJ_FR3_c <- VJ_FR3_c + 1}
          if(res_B[[n]] %in% VJ_CDR3[[i]][[1]]:VJ_CDR3[[i]][[2]]){VJ_CDR3_c <- VJ_CDR3_c + 1}
          if(res_B[[n]] %in% VJ_FR4[[i]][[1]]:VJ_FR4[[i]][[2]]){VJ_FR4_c <- VJ_FR4_c + 1}
        }

        bin.res.list[[i]] <- list("VDJ_FR1" = VDJ_FR1_c,
          "VDJ_CDR1" = VDJ_CDR1_c,
          "VDJ_FR2" = VDJ_FR2_c,
          "VDJ_CDR2" =VDJ_CDR2_c,
          "VDJ_FR3" = VDJ_FR3_c,
          "VDJ_CDR3" = VDJ_CDR3_c,
          "VDJ_FR4" = VDJ_FR4_c,
          "VJ_FR1" = VJ_FR1_c,
          "VJ_CDR1" = VJ_CDR1_c,
          "VJ_FR2" = VJ_FR2_c,
          "VJ_CDR2" = VJ_CDR2_c,
          "VJ_FR3" = VJ_FR3_c,
          "VJ_CDR3" = VJ_CDR3_c,
          "VJ_FR4" = VJ_FR4_c)

      } # End if binding residue barplot




    } #End of for loop over vis.structure

    #Creating a data from the out list
    df_out <- do.call(rbind.data.frame,out.list)




    #Adding to the return list
    return.list[[1]] <- df_out
    return.list[[length(return.list) + 1]] <- BindingResidues.list
    if(dist.mat) {return.list[[length(return.list) + 1]] <- dist.mat.list}





    if(binding.residue.barplot){
      bin.res.df <- do.call(rbind.data.frame,bin.res.list)
      bin.res.df$FrameWork <- bin.res.df %>% dplyr::select(dplyr::matches("*FR*")) %>% rowSums()
      bin.res.df$CDR <- bin.res.df %>% dplyr::select(dplyr::matches("*CDR*"))  %>% rowSums()


      if(binding.residue.barplot.style == "all_regions" | binding.residue.barplot.style == "both"){
        bar_plot <- list()
        for(pp in 1:nrow(bin.res.df)){

        bin.res.pp <- bin.res.df %>% dplyr::select(c("VDJ_CDR1","VDJ_FR2","VDJ_CDR2","VDJ_FR3","VDJ_CDR3","VDJ_FR4","VJ_FR1","VJ_FR2","VJ_CDR2","VJ_FR3","VJ_CDR3","VJ_FR4"))
        bin.res.pp <- bin.res.pp[pp,] %>% t() %>% data.frame()
        colnames(bin.res.pp) <- c("Nr_bind_res")
        bin.res.pp$region <- bin.res.pp %>% row.names()
        bin.res.pp$region <- factor(bin.res.pp$region, levels =  bin.res.pp$region)

        # Basic barplot
        bar_plot[[pp]] <-ggplot2::ggplot(data=bin.res.pp, ggplot2::aes(x=region, y=Nr_bind_res, fill=region)) +
          ggplot2::geom_bar(stat="identity") +
          ggplot2::coord_flip() +
          ggplot2::ylab("Nr. of residues in binding site") +
          ggplot2::xlab("Region") +
          ggplot2::ggtitle(cells.to.vis[[ceiling(pp/length(rank.list))]]) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none") +
          ggplot2::scale_fill_manual(values = c("VDJ_FR1" = color.frameworks,
                                       "VDJ_CDR1" = color.cdr1,
                                       "VDJ_FR2" = color.frameworks,
                                       "VDJ_CDR2" = color.cdr2,
                                       "VDJ_FR3" = color.frameworks,
                                       "VDJ_CDR3" = color.cdr3,
                                       "VDJ_FR4" = color.frameworks,
                                       "VJ_FR1" = color.frameworks,
                                       "VJ_CDR1" = color.cdr1,
                                       "VJ_FR2" = color.frameworks,
                                       "VJ_CDR2" = color.cdr2,
                                       "VJ_FR3" = color.frameworks,
                                       "VJ_CDR3" = color.cdr3,
                                       "VJ_FR4" = color.frameworks))



        }
        return.list[[length(return.list)+1]] <- bar_plot


        ### Plot overall usage
        overall_df <- bin.res.df %>% dplyr::select(c("VDJ_CDR1","VDJ_FR2","VDJ_CDR2","VDJ_FR3","VDJ_CDR3","VDJ_FR4","VJ_FR1","VJ_FR2","VJ_CDR2","VJ_FR3","VJ_CDR3","VJ_FR4"))
        overall_df <- mapply(sum,overall_df) %>% data.frame()
        names(overall_df) <- c("Sum_nr_bind_res")
        overall_df$region <- overall_df %>% row.names()
        overall_df$region <- factor(overall_df$region, levels =  overall_df$region)
        overall_df$Sum_nr_bind_res_pct <- overall_df$Sum_nr_bind_res / sum(overall_df$Sum_nr_bind_res)

        # Basic barplot
        bar_plot <-ggplot2::ggplot(data=overall_df, ggplot2::aes(x=region, y=Sum_nr_bind_res_pct, fill=region)) +
          ggplot2::geom_bar(stat="identity") +
          ggplot2::scale_y_continuous(labels=scales::percent) +
          ggplot2::coord_flip() +
          ggplot2::ylab("Percentage of binding site residues") +
          ggplot2::xlab("Region") +
          ggplot2::ggtitle("Position of binding site residues") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none") +
          ggplot2::scale_fill_manual(values = c("VDJ_FR1" = color.frameworks,
                                                "VDJ_CDR1" = color.cdr1,
                                                "VDJ_FR2" = color.frameworks,
                                                "VDJ_CDR2" = color.cdr2,
                                                "VDJ_FR3" = color.frameworks,
                                                "VDJ_CDR3" = color.cdr3,
                                                "VDJ_FR4" = color.frameworks,
                                                "VJ_FR1" = color.frameworks,
                                                "VJ_CDR1" = color.cdr1,
                                                "VJ_FR2" = color.frameworks,
                                                "VJ_CDR2" = color.cdr2,
                                                "VJ_FR3" = color.frameworks,
                                                "VJ_CDR3" = color.cdr3,
                                                "VJ_FR4" = color.frameworks))



        return.list[[length(return.list)+1]] <- bar_plot
        return.list[[length(return.list)+1]] <- overall_df
      } # if style == "all_regions"





      if(binding.residue.barplot.style == "FR_CDR" | binding.residue.barplot.style == "both"){
        bar_plot <- list()
        for(pp in 1:nrow(bin.res.df)){

          bin.res.pp <- bin.res.df %>% dplyr::select(c("FrameWork", "CDR"))
          bin.res.pp <- bin.res.pp[pp,] %>% t() %>% data.frame()
          colnames(bin.res.pp) <- c("Nr_bind_res")
          bin.res.pp$region <- bin.res.pp %>% row.names()
          bin.res.pp$region <- factor(bin.res.pp$region, levels =  bin.res.pp$region)

          # Basic barplot
          bar_plot[[pp]] <-ggplot2::ggplot(data=bin.res.pp, ggplot2::aes(x=region, y=Nr_bind_res, fill=region)) +
            ggplot2::geom_bar(stat="identity", width = 0.5) +
            ggplot2::ylab("Nr. of residues in binding site") +
            ggplot2::xlab("Region") +
            ggplot2::ggtitle(cells.to.vis[[ceiling(pp/length(rank.list))]]) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none") +
            ggplot2::scale_fill_manual(values = c("FrameWork" = color.frameworks,
                                                  "CDR" = color.cdr3))



        }
        return.list[[length(return.list)+1]] <- bar_plot

        ### Plot overall usage
        overall_df <- mapply(sum,bin.res.df %>% dplyr::select(c("FrameWork", "CDR"))) %>% data.frame()
        names(overall_df) <- c("Sum_nr_bind_res")
        overall_df$region <- overall_df %>% row.names()
        overall_df$region <- factor(overall_df$region, levels =  overall_df$region)
        overall_df$Sum_nr_bind_res_pct <- overall_df$Sum_nr_bind_res / sum(overall_df$Sum_nr_bind_res)

        # Basic barplot
        bar_plot <-ggplot2::ggplot(data=overall_df, ggplot2::aes(x=region, y=Sum_nr_bind_res_pct, fill=region)) +
          ggplot2::geom_bar(stat="identity",width = 0.5) +
          ggplot2::scale_y_continuous(labels=scales::percent) +
          ggplot2::ylab("Percentage of binding site residues") +
          ggplot2::xlab("Region") +
          ggplot2::ggtitle("Position of binding site residues") +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "none") +
          ggplot2::scale_fill_manual(values = c("FrameWork" = color.frameworks,
                                                "CDR" = color.cdr3))



        return.list[[length(return.list)+1]] <- bar_plot
        return.list[[length(return.list)+1]] <- overall_df



      }#End if style equals FR_CDR




    } # End of if binding.residue.plot

  } # End of if antigen.interaction



  #If antigen interaction is disabled
  if(antigen.interaction == F){

    metrics.list <- list()

    #Loop over the vis.structure
    for(STRUC in vis.structure){

      metrics_df <- data.frame(STRUC$atom)
      metrics_df$charge <- NULL

      #SASA
      if(SASA & !antigen.interaction){

        SASA_df <- vanddraabe::FreeSASA.diff(STRUC$atom)
        SASA_exp_df <- SASA_df$uniq.atom.ids %>% stringr::str_split("_") %>% do.call(rbind, .) %>% data.frame()
        names(SASA_exp_df) <- c("resid", "resno", "chain", "elety","eleno")
        SASA_df <- cbind(SASA_df,SASA_exp_df)
        SASA_df$resno <- as.integer(SASA_df$resno)
        SASA_df$eleno <- as.integer(SASA_df$eleno)
        SASA_df$uniq.atom.ids <- NULL

        metrics_df <- dplyr::left_join(metrics_df,SASA_df,by = c("resno","chain","elety","eleno","resid"))
      }

        #Hydrophobicity
        if(hydrophobicity){
          df_hydph <- data.frame()

          for(CHAIN in unique(STRUC$atom$chain)){
            df_chain <- STRUC$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
            df_chain$hydrophobicity <- STRUC$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::hydrophobicity()

            df_hydph <- rbind(df_hydph,df_chain)
          }

          metrics_df <- dplyr::left_join(metrics_df,df_hydph,by = c("resno","chain"))

        }

      #Charge
      if(charge){
        df_charge <- dplyr::data_frame()

        for(CHAIN in unique(STRUC$atom$chain)){
          df_chain <- STRUC$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% dplyr::select(c("resno","chain"))
          df_chain$charge <- STRUC$atom %>% dplyr::filter(chain == CHAIN) %>% dplyr::distinct(.,resno, .keep_all = T) %>% .$resid %>% stringr::str_to_title() %>% seqinr::a() %>% Peptides::charge()

          df_charge <- rbind(df_charge,df_chain)
        }

        metrics_df <- dplyr::left_join(metrics_df,df_charge,by = c("resno","chain"))


      }

      ## Add metric dataframe to list
      metrics.list[[length(metrics.list)+1]] <- metrics_df


    }#End of Vis.structure loop

  }#End of if not antigen.interaction




  ##Add to return list
  if(SASA | charge | hydrophobicity) {return.list[[length(return.list)+1]] <- metrics.list}

  if(metrics.plot){

    metrics.plot <- list()

    for(i in 1:length(vis.structure)){
      Plot.Struc <- vis.structure[[i]]

      Plot.Struc$atom <- metrics.list[[i]]

      if(SASA){

        #For loop to get rid of NA values
        out.vec <- c()
        CurVar <- 0
        for(i in metrics.list[[i]]$SASA.prot){
          if(!is.na(i)){CurVar <- i}
          out.vec <- c(out.vec,CurVar)
        }

        Plot.Struc$atom$b <- out.vec

        cartoon_styles <- r3dmol::m_style_cartoon()
        cartoon_styles$cartoon$colorscheme <- list(prop = "b", gradient = "roygb", min = 0, max = 100)

        SASA_plot <- r3dmol::r3dmol() %>%
          r3dmol::m_add_model(r3dmol::m_bio3d(Plot.Struc)) %>%
          r3dmol::m_set_style(style = cartoon_styles) %>%
          r3dmol::m_zoom_to() %>%
          r3dmol::m_spin(speed = spin.speed)

        metrics.plot[[length(metrics.plot)+1]] <- SASA_plot
      }

      if(hydrophobicity){

        Plot.Struc$atom$b <- Plot.Struc$atom$hydrophobicity

        cartoon_styles <- r3dmol::m_style_cartoon()
        cartoon_styles$cartoon$colorscheme <- list(prop = "b", gradient = "roygb", min = -4.5, max = 4.5)

        Hydph_plot <- r3dmol::r3dmol() %>%
          r3dmol::m_add_model(r3dmol::m_bio3d(Plot.Struc)) %>%
          r3dmol::m_set_style(style = cartoon_styles) %>%
          r3dmol::m_zoom_to() %>%
          r3dmol::m_spin(speed = spin.speed)

        metrics.plot[[length(metrics.plot)+1]] <- Hydph_plot

      }


      if(charge){

        Plot.Struc$atom$b <- Plot.Struc$atom$charge

        cartoon_styles <- r3dmol::m_style_cartoon()
        cartoon_styles$cartoon$colorscheme <- list(prop = "b", gradient = "roygb", min = -1, max = 1)

        Charge_plot <- r3dmol::r3dmol() %>%
          r3dmol::m_add_model(r3dmol::m_bio3d(Plot.Struc)) %>%
          r3dmol::m_set_style(style = cartoon_styles) %>%
          r3dmol::m_zoom_to() %>%
          r3dmol::m_spin(speed = spin.speed)

        metrics.plot[[length(metrics.plot)+1]] <- Charge_plot

      }





      return.list[[length(return.list)+1]] <- metrics.plot
    }


  }











  ###Plots
  if(structure.plot){

    if(VDJ.structure.input){

      #If VDJ annotation shall be visualized (requires mixcr output columns )
      if(VDJ.anno){

        if(anno.seq.input != F){

          SeqAnno <- list()
          for(i in 1:length(anno.seq)){

            if(length(anno.seq[[i]]) > 4){

              anno.label <- paste0(
                'r3dmol::m_add_label(
            text = "',anno.seq[[i]][[5]],'",
            sel = r3dmol::m_sel(resi = ceiling(',((as.integer(anno.seq[[i]][[1]]) + as.integer(anno.seq[[i]][[2]]))/2),'), chain = "',anno.seq[[i]][[3]],'"),
            style = r3dmol::m_style_label(
            backgroundColor = "',anno.seq[[i]][[4]],'",
            backgroundOpacity = 0.9,
            fontSize = label.size,
            fontOpacity = 1,
            fontColor = font.col,
            alignment = "center")
            ) %>%' )

            } else{anno.label <- ""}

            SeqAnno[[i]] <- paste0('r3dmol::m_set_style(
          sel = r3dmol::m_sel(resi = ',anno.seq[[i]][[1]],':',anno.seq[[i]][[2]],',chain = "',anno.seq[[i]][[3]],'"),
          style = r3dmol::m_style_cartoon(color = "',anno.seq[[i]][[4]],'", arrows = TRUE)
          ) %>%',
                                   anno.label)
          }
        }
        else{
          SeqAnno <- ""
        }





        #Define the location of the given elements
        VDJ_FR1 <- list()
        VDJ_CDR1 <- list()
        VDJ_FR2 <- list()
        VDJ_CDR2 <- list()
        VDJ_FR3 <- list()
        VDJ_CDR3 <- list()
        VDJ_FR4 <- list()

        VJ_FR1 <- list()
        VJ_CDR1 <- list()
        VJ_FR2 <- list()
        VJ_CDR2 <- list()
        VJ_FR3 <- list()
        VJ_CDR3 <- list()
        VJ_FR4 <- list()


        for(i in 1:length(cells.to.vis)){
          for(n in 1:(length(vis.structure)/length(cells.to.vis))){
            sequnce.to.vis <- VDJ.structure[[1]] %>% dplyr::filter(.,barcode == cells.to.vis[[i]])
            VDJ_aa <- sequnce.to.vis$VDJ_aa_mixcr

            VDJ_FR1[[length(VDJ_FR1)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR1)
            VDJ_CDR1[[length(VDJ_CDR1)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR1)
            VDJ_FR2[[length(VDJ_FR2)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR2)
            VDJ_CDR2[[length(VDJ_CDR2)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR2)
            VDJ_FR3[[length(VDJ_FR3)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR3)
            VDJ_CDR3[[length(VDJ_CDR3)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqCDR3)
            VDJ_FR4[[length(VDJ_FR4)+1]] <- stringr::str_locate(VDJ_aa,sequnce.to.vis$VDJ_aaSeqFR4)

            VJ_aa <- sequnce.to.vis$VJ_aa_mixcr

            VJ_FR1[[length(VJ_FR1)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR1)
            VJ_CDR1[[length(VJ_CDR1)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR1)
            VJ_FR2[[length(VJ_FR2)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR2)
            VJ_CDR2[[length(VJ_CDR2)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR2)
            VJ_FR3[[length(VJ_FR3)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR3)
            VJ_CDR3[[length(VJ_CDR3)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqCDR3)
            VJ_FR4[[length(VJ_FR4)+1]] <- stringr::str_locate(VJ_aa,sequnce.to.vis$VJ_aaSeqFR4)

          }
        }

        models <- list()
        for(i in 1:length(vis.structure)){
          models[[i]] <- paste0('# Add model to scene
          r3dmol::m_add_model(data = r3dmol::m_bio3d(vis.structure[[',i,']])) %>% ')
        }


        anno <- list()
        for(i in 1:length(vis.structure)){

          anno[[i]] <- paste0(
            '
          ## Annotate VDJ

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_FR1[[',i,']][[1]]:VDJ_FR1[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%


          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_CDR1[[',i,']][[1]]:VDJ_CDR1[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.cdr1, arrows = TRUE)
          ) %>%

          # Add a label
          r3dmol::m_add_label(
            text = "VDJ_CDR1",
            sel = r3dmol::m_sel(resi = ceiling((VDJ_CDR1[[',i,']][[1]] + VDJ_CDR1[[',i,']][[2]])/2), chain = "A"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr1,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")) %>%


          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_FR2[[',i,']][[1]]:VDJ_FR2[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%



          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_CDR2[[',i,']][[1]]:VDJ_CDR2[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.cdr2, arrows = TRUE)
          ) %>%

          # Add a label
          r3dmol::m_add_label(
            text = "VDJ_CDR2",
            sel = r3dmol::m_sel(resi = ceiling((VDJ_CDR2[[',i,']][[1]] + VDJ_CDR2[[',i,']][[2]])/2), chain = "A"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr2,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")
          ) %>%

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_FR3[[',i,']][[1]]:VDJ_FR3[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%


          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_CDR3[[',i,']][[1]]:VDJ_CDR3[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.cdr3, arrows = TRUE)
          ) %>%

          # Add a label
          r3dmol::m_add_label(
            text = "VDJ_CDR3",
            sel = r3dmol::m_sel(resi = ceiling((VDJ_CDR3[[',i,']][[1]] + VDJ_CDR3[[',i,']][[2]])/2), chain = "A"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr3,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")
          ) %>%

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VDJ_FR4[[',i,']][[1]]:VDJ_FR4[[',i,']][[2]], chain = "A"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%





          ## Annotate VJ

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_FR1[[',i,']][[1]]:VJ_FR1[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%



          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_CDR1[[',i,']][[1]]:VJ_CDR1[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.cdr1, arrows = TRUE)
          ) %>%

          r3dmol::m_add_label(
            text = "VJ_CDR1",
            sel = r3dmol::m_sel(resi = ceiling((VJ_CDR1[[',i,']][[1]] + VJ_CDR1[[',i,']][[2]])/2), chain = "B"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr1,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")
          ) %>%

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_FR2[[',i,']][[1]]:VJ_FR2[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%



          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_CDR2[[',i,']][[1]]:VJ_CDR2[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.cdr2, arrows = TRUE)
          ) %>%

          r3dmol::m_add_label(
            text = "VJ_CDR2",
            sel = r3dmol::m_sel(resi = ceiling((VJ_CDR2[[',i,']][[1]] + VJ_CDR2[[',i,']][[2]])/2), chain = "B"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr2,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")
          ) %>%

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_FR3[[',i,']][[1]]:VJ_FR3[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%



          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_CDR3[[',i,']][[1]]:VJ_CDR3[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.cdr3, arrows = TRUE)
          ) %>%

          r3dmol::m_add_label(
            text = "VJ_CDR3",
            sel = r3dmol::m_sel(resi = ceiling((VJ_CDR3[[',i,']][[1]] + VJ_CDR3[[',i,']][[2]])/2), chain = "B"),
            style = r3dmol::m_style_label(
              backgroundColor = color.cdr3,
              backgroundOpacity = bk.opac,
              fontSize = label.size,
              fontOpacity = font.opac,
              fontColor = font.col,
              alignment = "center")
          ) %>%

          # Annotate
          r3dmol::m_set_style(
            sel = r3dmol::m_sel(resi = VJ_FR4[[',i,']][[1]]:VJ_FR4[[',i,']][[2]], chain = "B"),
            style = r3dmol::m_style_cartoon(color = color.frameworks, arrows = TRUE)
          ) %>%

        ')
        }




      if(overlay){
        out.plot <- eval(parse(text = paste0(
            '# Set up the initial viewer
          r3dmol::r3dmol(
            viewer_spec = r3dmol::m_viewer_spec(
              cartoonQuality = 10,
              lowerZoomLimit = 50,
              upperZoomLimit = 350
            )

          ) %>%

            ',paste0(models,collapse = ""),'

            # Zoom to encompass the whole scene
            r3dmol::m_zoom_to() %>%
            # Set style of structures
            r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
            # Set style of specific selection (selecting by secondary)
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "s"),
              style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
            ) %>%
            # Style the alpha helix
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
              style = r3dmol::m_style_cartoon(color = "',color.helix,'")
            ) %>%

            ',paste0(anno, collapse = ""),'

            ',paste0(SeqAnno,collapse = ""),'


            # Rotate the scene by given angle on given axis
            r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
            r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
            r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
            # Animate the scene by spinning it
            r3dmol::m_spin(speed = spin.speed)',r3dmol.code)
          ))

        return.list[[length(return.list)+1]] <- out.plot
      }




      if(!overlay){
        for(ii in 1:length(vis.structure)){
          out.plot <- eval(parse(text = paste0(
            '# Set up the initial viewer
            r3dmol::r3dmol(
              viewer_spec = r3dmol::m_viewer_spec(
                cartoonQuality = 10,
                lowerZoomLimit = 50,
                upperZoomLimit = 350
              )

            ) %>%

              ',models[[ii]],'

              # Zoom to encompass the whole scene
              r3dmol::m_zoom_to() %>%
              # Set style of structures
              r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
              # Set style of specific selection (selecting by secondary)
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "s"),
                style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
              ) %>%
              # Style the alpha helix
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
                style = r3dmol::m_style_cartoon(color = "',color.helix,'")
              ) %>%

              ',anno[[ii]],'

              ',paste0(SeqAnno,collapse = ""),'


              # Rotate the scene by given angle on given axis
              r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
              r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
              r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
              # Animate the scene by spinning it
              r3dmol::m_spin(speed = spin.speed)',r3dmol.code)
          ))

          return.list[[length(return.list)+1]] <- out.plot
        }
      }















      }#End if VDj_anno


      if(VDJ.anno == F){

        if(anno.seq.input != F){

          SeqAnno <- list()
          for(i in 1:length(anno.seq)){

            if(length(anno.seq[[i]]) > 4){

              anno.label <- paste0(
                'r3dmol::m_add_label(
            text = "',anno.seq[[i]][[5]],'",
            sel = r3dmol::m_sel(resi = ceiling(',((as.integer(anno.seq[[i]][[1]]) + as.integer(anno.seq[[i]][[2]]))/2),'), chain = "',anno.seq[[i]][[3]],'"),
            style = r3dmol::m_style_label(
            backgroundColor = "',anno.seq[[i]][[4]],'",
            backgroundOpacity = 0.9,
            fontSize = label.size,
            fontOpacity = 1,
            fontColor = font.col,
            alignment = "center")
            ) %>%' )

            } else{anno.label <- ""}

            SeqAnno[[i]] <- paste0('r3dmol::m_set_style(
          sel = r3dmol::m_sel(resi = ',anno.seq[[i]][[1]],':',anno.seq[[i]][[2]],',chain = "',anno.seq[[i]][[3]],'"),
          style = r3dmol::m_style_cartoon(color = "',anno.seq[[i]][[4]],'", arrows = TRUE)
          ) %>%',
                                   anno.label)
          }
        }
        else{
          SeqAnno <- ""
        }





        models <- list()
        for(i in 1:length(cells.to.vis)){
          models[[i]] <- paste0('# Add model to scene
          r3dmol::m_add_model(data = r3dmol::m_bio3d(vis.structure[[',i,']])) %>% ')
        }




        if(overlay){
          out.plot <- eval(parse(text = paste0(
            '# Set up the initial viewer
          r3dmol::r3dmol(
            viewer_spec = r3dmol::m_viewer_spec(
              cartoonQuality = 10,
              lowerZoomLimit = 50,
              upperZoomLimit = 350
            )

          ) %>%

            ',paste0(models,collapse = ""),'

            # Zoom to encompass the whole scene
            r3dmol::m_zoom_to() %>%
            # Set style of structures
            r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
            # Set style of specific selection (selecting by secondary)
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "s"),
              style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
            ) %>%
            # Style the alpha helix
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
              style = r3dmol::m_style_cartoon(color = "',color.helix,'")
            ) %>%

            ',paste0(SeqAnno,collapse = ""),'

            # Rotate the scene by given angle on given axis
            r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
            r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
            r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
            # Animate the scene by spinning it
            r3dmol::m_spin(speed = spin.speed)')
          ))

          return.list[[length(return.list)+1]] <- out.plot
        }


        if(!overlay){
          for(ii in 1:length(cells.to.vis)){
            out.plot <- eval(parse(text = paste0(
              '# Set up the initial viewer
            r3dmol::r3dmol(
              viewer_spec = r3dmol::m_viewer_spec(
                cartoonQuality = 10,
                lowerZoomLimit = 50,
                upperZoomLimit = 350
              )

            ) %>%

              ',models[[ii]],'

              # Zoom to encompass the whole scene
              r3dmol::m_zoom_to() %>%
              # Set style of structures
              r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
              # Set style of specific selection (selecting by secondary)
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "s"),
                style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
              ) %>%
              # Style the alpha helix
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
                style = r3dmol::m_style_cartoon(color = "',color.helix,'")
              ) %>%

              ',paste0(SeqAnno,collapse = ""),'

              # Rotate the scene by given angle on given axis
              r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
              r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
              r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
              # Animate the scene by spinning it
              r3dmol::m_spin(speed = spin.speed)',r3dmol.code)
            ))

            return.list[[length(return.list)+1]] <- out.plot
          }
        }












      } #End of if VDJ_anno == F
    } #End of VDJ.structure != F

    if(PDB.file.input){


      if(anno.seq.input != F){

        SeqAnno <- list()
        for(i in 1:length(anno.seq)){

          if(length(anno.seq[[i]]) > 4){

            anno.label <- paste0(
            'r3dmol::m_add_label(
            text = "',anno.seq[[i]][[5]],'",
            sel = r3dmol::m_sel(resi = ceiling(',((as.integer(anno.seq[[i]][[1]]) + as.integer(anno.seq[[i]][[2]]))/2),'), chain = "',anno.seq[[i]][[3]],'"),
            style = r3dmol::m_style_label(
            backgroundColor = "',anno.seq[[i]][[4]],'",
            backgroundOpacity = 0.9,
            fontSize = label.size,
            fontOpacity = 1,
            fontColor = font.col,
            alignment = "center")
            ) %>%' )

          } else{anno.label <- ""}

          SeqAnno[[i]] <- paste0('r3dmol::m_set_style(
          sel = r3dmol::m_sel(resi = ',anno.seq[[i]][[1]],':',anno.seq[[i]][[2]],',chain = "',anno.seq[[i]][[3]],'"),
          style = r3dmol::m_style_cartoon(color = "',anno.seq[[i]][[4]],'", arrows = TRUE)
          ) %>%',
          anno.label)
        }
      }
      else{
        SeqAnno <- ""
      }





      models <- list()
      for(i in 1:length(PDB.file)){
        models[[i]] <- paste0('# Add model to scene
          r3dmol::m_add_model(data = r3dmol::m_bio3d(vis.structure[[',i,']])) %>% ')
      }






      if(overlay){
        out.plot <- eval(parse(text = paste0(
          '# Set up the initial viewer
          r3dmol::r3dmol(
            viewer_spec = r3dmol::m_viewer_spec(
              cartoonQuality = 10,
              lowerZoomLimit = 50,
              upperZoomLimit = 350
            )

          ) %>%

            ',paste0(models,collapse = ""),'

            # Zoom to encompass the whole scene
            r3dmol::m_zoom_to() %>%
            # Set style of structures
            r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
            # Set style of specific selection (selecting by secondary)
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "s"),
              style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
            ) %>%
            # Style the alpha helix
            r3dmol::m_set_style(
              sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
              style = r3dmol::m_style_cartoon(color = "',color.helix,'")
            ) %>%

            ',paste0(SeqAnno,collapse = ""),'

            # Rotate the scene by given angle on given axis
            r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
            r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
            r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
            # Animate the scene by spinning it
            r3dmol::m_spin(speed = spin.speed)',r3dmol.code)
        ))

        return.list[[length(return.list)+1]] <- out.plot
      }


      if(!overlay){
        for(ii in 1:length(PDB.file)){
          out.plot <- eval(parse(text = paste0(
            '# Set up the initial viewer
            r3dmol::r3dmol(
              viewer_spec = r3dmol::m_viewer_spec(
                cartoonQuality = 10,
                lowerZoomLimit = 50,
                upperZoomLimit = 350
              )

            ) %>%

              ',models[ii],'

              # Zoom to encompass the whole scene
              r3dmol::m_zoom_to() %>%
              # Set style of structures
              r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "',color.molecule,'")) %>%
              # Set style of specific selection (selecting by secondary)
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "s"),
                style = r3dmol::m_style_cartoon(color = "',color.sheets,'", arrows = TRUE)
              ) %>%
              # Style the alpha helix
              r3dmol::m_set_style(
                sel = r3dmol::m_sel(ss = "h"), # Style alpha helix
                style = r3dmol::m_style_cartoon(color = "',color.helix,'")
              ) %>%

              ',paste0(SeqAnno,collapse = ""),'

              # Rotate the scene by given angle on given axis
              r3dmol::m_rotate(angle = angle.x, axis = "x") %>%
              r3dmol::m_rotate(angle = angle.y, axis = "y") %>%
              r3dmol::m_rotate(angle = angle.z, axis = "z") %>%
              # Animate the scene by spinning it
              r3dmol::m_spin(speed = spin.speed)',r3dmol.code)
          ))

          return.list[[length(return.list)+1]] <- out.plot
        }
      }
    } # End of PDB_file
  }



  if(plddt.plot){

      for(i in 1:length(vis.structure)){
        cartoon_styles <- r3dmol::m_style_cartoon()
        cartoon_styles$cartoon$colorscheme <- list(prop = "b", gradient = "roygb", min = 50, max = 100)

        plddt_plot <- r3dmol::r3dmol() %>%
          r3dmol::m_add_model(r3dmol::m_bio3d(vis.structure[[i]])) %>%
          r3dmol::m_set_style(style = cartoon_styles) %>%
          r3dmol::m_zoom_to() %>%
          r3dmol::m_spin(speed = spin.speed)


        return.list[[length(return.list)+1]] <- plddt_plot
      }


  }#End if plddt plot




  return(return.list)


} # End of function











