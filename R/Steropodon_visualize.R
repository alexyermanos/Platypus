#' Plots a Steropodon modelled structure using r3dmol

#' @description Plots a structure modelled via Steropodon_model using r3dmol.
#' Structures can be colored by physicochemical properties, framework/hypervariable regions, chains, epitopes/paratopes, and invariant core regions, as determined in the color.by parameter.

#' @param steropodon.object  a nested list of predicted structure objects (per sample, per clonotype) or a single Steropodon object.
#' @param structure string - the structure saved inside the Steropodon object to be chosen: 'structure' for the whole receptor structure (VDJ and VJ chains),'H' for the heavy chain, 'L' for the light chain,
#' 'CDRH3' for the CDR3 region of the heavy chain, 'CDRL3' for the CDR3 region in the light chain, 'paratope' for the paratope structure (after using Steropodon_dock), 'epitope' for the epitope structure (after using Steropodon_dock),
#' 'core' for the core/structurally non-variable region across all structures in the Steropodon nested list (after using the Steropodon_find_core function), 'complex' for the modelled antibody-antigen complex (after using Steropodon_dock).
#' @param color.by string - feature in the Steropodon object to color the structure by.
#' @param single.plot boolean - if TRUE, will pool all structures into a single plot.
#' @param specific.color.values named list - list of feature names and colors (list('#4287f5'='VDJ_CDR1') to color only VDJ_CDR1 regions ).
#' @param additional.r3dmol.code named list - additional r3dmol code for the final visualizations/plots.
#' @param show.label boolean - if TRUE, will also show a label with the feature name as specified by color.by.
#' @param animate boolean - if TRUE, will create an animation of all structures in a given (nested) list of Steropodon objects.

#' @return No returns. Will output an r3dmol visualization of the structure(s) in the plots tab on RStudio.
#' @export
#' @examples
#' \dontrun{
#'steropodon_igfold$s1$clonotype1$`1` %>%
#'  Steropodon_trim(structure = 'structure',
#'                  grouping = c('chain', 'region'),
#'                  specific.values = c('VDJ_CDR1',
#'                  'VDJ_CDR2',
#'                  'VDJ_CDR3',
#'                  'VJ_CDR1',
#'                  'VJ_CDR2',
#'                  'VJ_CDR3'),
#'                  combine.values = T,
#'                  combine.groupings = T) %>%
#'  Steropodon_visualize(structure = 'structure',
#'                       color.by = c('chain', 'region'))
#'}


Steropodon_visualize <- function(steropodon.object,
                                 structure,
                                 color.by,
                                 single.plot,
                                 specific.color.values,
                                 additional.r3dmol.code,
                                 show.label,
                                 animate
                               ){

  if(missing(steropodon.object)) stop('Please input your Steropodon object/ a (nested) list of objects!')
  if(missing(structure)) structure <- 'structure'
  if(missing(color.by)) color.by <- c('chain','region')
  if(missing(single.plot)) single.plot <- 'overlap'
  if(missing(specific.color.values)) specific.color.values <- c(
                                                                '#4287f5'='VDJ_CDR1','#4287f5'='VJ_CDR1',
                                                                '#923fd1'='VDJ_CDR2','#923fd1'='VJ_CDR2',
                                                                '#e31948'='VDJ_CDR3','#e31948'='VJ_CDR3'
                                                                )

  if(missing(additional.r3dmol.code)) additional.r3dmol.code <- NULL
  if(missing(show.label)) show.label <- T
  if(missing(animate)) animate <- F

  if(animate){
    single.plot <- 'overlap'
  }

  combined <- NULL

  preprocess_pdbs <- function(pdb, color.by, numeric.feature){
    chains <- unique(pdb$atom$chain)
    chain_dict <- toupper(letters)[1:length(chains)]
    names(chain_dict) <- chains
    pdb$atom$chain <- chain_dict[pdb$atom$chain]

    if(numeric.feature){
      pdb$atom$b <- pdb$atom[[color.by]]
    }

    return(pdb)
  }

  get_features_position_list <- function(pdb, color.by, specific.color.values, numeric.feature){

    if(!numeric.feature){
      if(length(color.by) > 1){
        new_feature <- paste0(color.by, collapse = '_')
        pdb$atom$combined <- rep('', nrow(pdb$atom))

        for(col in color.by){
          pdb$atom$combined <- paste0(pdb$atom$combined, '_', pdb$atom[[col]])
        }

        pdb$atom <- pdb$atom %>%
          dplyr::mutate(combined = substring(combined,2))

        colnames(pdb$atom)[colnames(pdb$atom) == 'combined'] <- new_feature
        color.by <- new_feature
      }

      if(is.null(specific.color.values)){
        specific.color.values <- unique(pdb$atom[[color.by]])
        names(specific.color.values) <- grDevices::rainbow(length(specific.color.values))
      }

      features_positions_list <- list()
      features_positions_list$resid <- list()
      features_positions_list$color <- list()
      features_positions_list$feature <- list()
      features_positions_list$feature_name <- color.by

      i <- 1
      while(length(specific.color.values) > 0){
        if(length(which(pdb$atom[[color.by]] == specific.color.values[1])) == 0){
          specific.color.values <- specific.color.values[specific.color.values!=specific.color.values[1]]
          next
        }

        features_positions_list$resid[[i]] <- unique(pdb$atom$resno[pdb$atom[[color.by]] == specific.color.values[1]])
        features_positions_list$color[[i]] <- names(specific.color.values)[1]
        features_positions_list$feature[[i]] <- specific.color.values[1]
        specific.color.values <- specific.color.values[specific.color.values!=specific.color.values[1]]
        i <- i + 1
      }
    }else{
      if(color.by == 'hydrophobicity'){
        features_positions_list <- list(prop = "b", gradient = "roygb", min = -4.5, max = 4.5)

      }else if(color.by == 'charge'){
        features_positions_list <- list(prop = "b", gradient = "roygb", min = -1, max = 1)

      }else if(color.by == 'SASA'){
        features_positions_list <- list(prop = "b", gradient = "roygb", min = 0, max = 100)

      }else if(color.by == 'b'){
        features_positions_list <- list(prop = "b", gradient = "roygb", min = 50, max = 100)

      }else{
        features_positions_list <- list(prop = "b", gradient = "roygb", min = min(pdb$atom[[color.by]]), max = max(pdb$atom[[color.by]]))
      }
    }

    return(features_positions_list)
  }

  create_plot <- function(pdb.list,
                          position.list,
                          numeric.feature,
                          additional.r3dmol.code,
                          single.plot,
                          show.label,
                          animate
                          ){

    #TO DO: color each structure w a different shade (to better notice structure differences)
    #colfunc_heavy <- grDevices::colorRampPalette(c("#fabea2", "#f56f31"))
    #colfunc_light <- grDevices::colorRampPalette(c("#f5e598", "#fcda35"))
    #col_list_heavy <- colfunc_heavy(length(pdb.list))
    #col_list_light <- colfunc_light(length(pdb.list))

    if(numeric.feature){
      cartoon_styles <- r3dmol::m_style_cartoon()
      cartoon_styles$cartoon$colorscheme <- position.list[[1]]
    }

    plot_string <- list()
    plot_grid <- list()
    for(i in 1:length(pdb.list)){
      plot_string[[i]] <- paste0('r3dmol::m_add_model(data = r3dmol::m_bio3d(pdb.list[[',i,']])) %>% r3dmol::m_zoom_to() %>% r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = "#a09fa1"))')

      if(numeric.feature){
        feat_string <- 'r3dmol::m_set_style(style = cartoon_styles)'
        plot_string[[i]] <- paste0(plot_string[[i]], ' %>% ', feat_string)

      }else{
        for(j in 1:length(position.list[[i]]$feature)){
          resid <- paste0(position.list[[i]]$resid[[j]], collapse = ',')
          resid <- paste0('c(', resid, ')')
          color <- position.list[[i]]$color[[j]]
          label <- position.list[[i]]$feature[[j]]

          if(show.label){
            feat_string <- paste0('r3dmol::m_set_style(sel = r3dmol::m_sel(resi = ', resid, ', model = -1), style = r3dmol::m_style_cartoon(color = "', color, '")) %>%
                                 r3dmol::m_add_label(text = "', label, '", sel = r3dmol::m_sel(resi = ', resid, ', model = -1), style = r3dmol::m_style_label(backgroundColor = "', color,'", backgroundOpacity = 0.8, fontSize = 12, fontOpacity = 1, fontColor = "black", alignment = "center"))')
          }else{
            feat_string <- paste0('r3dmol::m_set_style(sel = r3dmol::m_sel(resi = ', resid, ', model = -1), style = r3dmol::m_style_cartoon(color = "', color, '"))')
          }


          plot_string[[i]] <- paste0(plot_string[[i]], ' %>% ', feat_string)
        }
      }

      if(!is.null(additional.r3dmol.code)){
        plot_string[[i]] <- paste0(plot_string[[i]], ' %>% ', additional.r3dmol.code)
      }

      if(single.plot == 'grid'){
        plot_grid[[i]] <- paste0('r3dmol::r3dmol() %>% ', plot_string[[i]])
        plot_grid[[i]] <- parse(text = plot_grid[[i]])
      }
    }

    if(single.plot == 'grid'){
      out <- r3dmol::m_grid(
        viewer = plot_grid,
        control_all = TRUE,
        viewer_config = r3dmol::m_viewer_spec(
          backgroundColor = "black"
        )
      )

    }else if(single.plot == 'overlap'){
      plot_string <- paste0(plot_string, collapse = ' %>% ')
      plot_string <- paste0('r3dmol::r3dmol() %>% ', plot_string)
      if(animate){
        plot_string <- paste0(plot_string, 'r3dmol::m_animate(list(loop = "forward")', sep = ' %>% ')
      }
      eval(parse(text = plot_string))

    }else{
      stop('Unrecognized method for combining structure plots!')
    }
  }

  if(inherits(steropodon.object, 'Steropodon')){
    pdb_list <- list()
    pdb_list[[1]] <- select_structure(steropodon.object, structure = structure)
    names(pdb_list) <- paste0(steropodon.object@structure_id, collapse = '.')

  }else if(inherits(steropodon.object, 'list')){
    if(inherits(steropodon.object[[1]], 'Steropodon')){
      pdb_list <- lapply(steropodon.object, function(x) select_structure(x, structure = structure))
      pdb_names <- lapply(steropodon.object, function(x) paste0(x@structure_id, collapse = '.'))
      names(pdb_list) <- pdb_names

    }else{
      steropodon_list <- unnest_steropodon(steropodon.object)
      pdb_list <- lapply(steropodon_list, function(x) select_structure(x, structure = structure))
      names(pdb_list) <- names(steropodon_list)
    }

  }else{
    stop('Unrecognized input: please input a single Steropodon object, a list of objects, or a nested list of Steropodon objects (per sample, per clonotype)')
  }

  numeric.feature <- F
  if(length(color.by) == 1){
    if(is.numeric(pdb_list[[1]]$atom[[color.by]])){
      numeric.feature <- T
    }
  }

  feature_positions_list <- lapply(pdb_list, function(x) get_features_position_list(x, color.by = color.by, specific.color.values = specific.color.values, numeric.feature = numeric.feature))
  pdb_list <- lapply(pdb_list, function(x) preprocess_pdbs(x, color.by = color.by, numeric.feature = numeric.feature))

  create_plot(pdb.list = pdb_list,
              position.list = feature_positions_list,
              numeric.feature = numeric.feature,
              additional.r3dmol.code = additional.r3dmol.code,
              single.plot = single.plot,
              show.label = show.label,
              animate = animate)
}

