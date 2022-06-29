#' Wrapper function of VDJ_antigen_integrate function
#'
#'@description Wraps the VDJ_antigen_integrate function and uses it to annotate a VDJ dataframe with antigen information. Needs to VDJ_db_load to be executed first, with preprocess=T and vgm.names=T to obtain the same column names as in the VDJ (to allow for sequence matching).
#' @param VDJ VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.
#' @param db.list list of database dataframes or csv file paths, obtained from VDJ_db_load with .
#' @param database.features list of features/column names to be integrated from the databases.
#' @param match string - sequences by which to match and integrate the antigen information. Currently, only 'cdr3.aa' and 'cdrh3.aa' are supported, as all databases have these two sequence types ('VJ_cdr3s_aa','VDJ_cdr3s_aa').
#' @param homology string - 'exact' for exact sequence matchings, 'homology' for homology matching.
#' @param lv.distance integer - maximum Levehnstein distance threshold for the homology matchings.

#' @return VDJ with new columns - antigen information integrated from the antigen databases.
#' @export
#' @examples
#' \dontrun{
#' VDJ_db_annotate(VDJ=VDJ,db.list=db.list,database.features='Epitope',match='cdr3.aa',homology=FALSE)
#'}

VDJ_db_annotate <- function(VDJ, db.list, database.features, match, homology, lv.distance){
  if(missing(VDJ)) stop('Please input your data as a VDJ')
  if(missing(db.list)) stop('Please input a named list of the databases you want annotated')
  if(missing(database.features)) stop('Please input the feature columns you want your VDJ to be annotated with. Make sure the names are present in all databases - use vgm.names=T or keep.only.common=T in VDJ_db_load')
  if(missing(match)) match <- 'cdr3.aa'
  if(missing(homology)) homology <- T
  if(missing(lv.distance) & homology==T) lv.distance <- 16

  if(homology){
    matching_type <- 'homology'
  }else{
    matching_type <- 'exact'
  }

  output_VDJ <- VDJ_antigen_integrate(VDJ,
                                       db.list,
                                       antigen.features=database.features,
                                       match.by=match,
                                       matching.type=matching_type,
                                       distance.threshold=lv.distance,
                                       sample.id=F,
                                       VDJ.VJ.1chain=T,
                                       aberrant.chosen.sequences=F,
                                       output.format='vgm')

  return(output_VDJ)
}
