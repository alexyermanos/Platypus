#' Infer germline from the desired software/caller

#'@description Function to infer the germline from the tree

#' @param VDJ VDJ dataframe obtained after calling VDJ_call_MIXCR or any other germline you want to use
#' @param germlines.from MIXCR or any other germline caller - default: MIXCR
#' @param VDJ.only boolean - if T, only Heavy Chain (VDJ) germline will be inferred

#' @return VDJ with the updated germline
#' @export
#' @examples
#' \dontrun{
#' VDJ_germline(VDJ, germlines.from='MIXCR',
#' VDJ.only=T)
#'}



VDJ_germline <- function(VDJ, germlines.from, VDJ.only){

  if(missing(VDJ))
    stop('Input the VDJ dataframe obtained after calling VDJ_call_MIXCR or any')
  if(missing(germlines.from))
    germlines.from <- "MIXRC"
  if(missing(VDJ.only))
    VDJ.only <- F

  sample_id <- NULL
  clonotype_id <- NULL
  VDJ_nSeqVRegion <- NULL
  VDJ_nSeqDRegion <- NULL
  VDJ_nSeqJRegion <- NULL
  VDJ_nSeqFR4 <- NULL
  VJ_nSeqVRegion <- NULL
  VJ_nSeqJRegion <- NULL
  VJ_nSeqFR4 <- NULL


  vgm_sc[[1]][["VDJ_germline"]] <- ""
  vgm_sc[[1]][["VJ_germline"]] <- ""

  if(germlines.from == "MIXCR" | germlines.from == "mixcr") {

    for (si in unique(VDJ$sample_id)){
      for (ci in unique(VDJ[VDJ$sample_id==si,]$clonotype_id)) {

        this_clone = subset(VDJ, sample_id==si & clonotype_id==ci)

        # cat(si, " ", ci, "\n")
        # print(table(this_clone$VDJ_nSeqVRegion))

        ##### VDJ - HEAVY CHAIN #####

        if (!(all(this_clone$VDJ_nSeqVRegion ==""))) {# | all(this_clone$VDJ_nSeqDRegion=="")

          VDJ[VDJ$sample_id==si & VDJ$clonotype_id==ci,]$VDJ_germline <-
            paste0(
              ### V GENE GERMLINE
              names(which.max(table(subset(this_clone, VDJ_nSeqVRegion!="")$VDJ_nSeqVRegion))),
              ### D GENE GERMLINE
              ifelse(!(all(this_clone$VDJ_nSeqDRegion =="")), names(which.max(table(subset(this_clone, VDJ_nSeqDRegion!="")$VDJ_nSeqDRegion))),""), # TODO might want to use CDR3 from MIXCR instead of empty when D gene is missing (but part of CDR3 is in J), discuss in meeting
              ### J GENE GERMLINE
              ifelse(!(all(this_clone$VDJ_nSeqJRegion =="")), names(which.max(table(subset(this_clone, VDJ_nSeqJRegion!="")$VDJ_nSeqJRegion))),
                     ifelse(!(all(this_clone$VDJ_nSeqFR4 =="")), names(which.max(table(subset(this_clone, VDJ_nSeqFR4!="")$VDJ_nSeqFR4))),"")))
        }
        else { # if we are missing V we leave the germline empty
          #TODO try using which.max(table(VDJ_nt_MIXCR)) and see if it works better; or just another sequence taken from VDJ_nt_MIXCR)
          VDJ[VDJ$sample_id==si & VDJ$clonotype_id==ci,]$VDJ_germline <- ""
          # <- VDJ[["VDJ_nt_mixcr"]][VDJ$sample_id==si & VDJ$clonotype_id==ci,]
        }
        if(!VDJ.only) {
          ##### VJ - LIGHT CHAIN #####

          if (!(all(this_clone$VJ_nSeqVRegion ==""))) {

            VDJ[VDJ$sample_id==si & VDJ$clonotype_id==ci,]$VJ_germline <-
              paste0(
                ### V GENE GERMLINE
                names(which.max(table(subset(this_clone, VJ_nSeqVRegion!="")$VJ_nSeqVRegion))),
                ### J GENE GERMLINE
                ifelse(!(all(this_clone$VJ_nSeqJRegion =="")), names(which.max(table(subset(this_clone, VJ_nSeqJRegion!="")$VJ_nSeqJRegion))),
                       ifelse(!(all(this_clone$VJ_nSeqFR4 =="")), names(which.max(table(subset(this_clone, VJ_nSeqFR4!="")$VJ_nSeqFR4))),"")))
            }
            else {
              VDJ[VDJ$sample_id==si & VDJ$clonotype_id==ci,]$VJ_germline <- ""
              }
          }
        }

      }

    VDJ$VDJ_trimmed_ref <- VDJ[["VDJ_germline"]]
    if(!VDJ.only)
      VDJ$VJ_trimmed_ref <- VDJ[["VJ_germline"]]
    }

  return(VDJ)
}
