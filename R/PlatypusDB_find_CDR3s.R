#' Loads and saves RData objects from the PlatypusDB
#' @param VDJ.cdr3s.aa Character A VDJ CDR3s amino acid sequence to search for
#' @param VJ.cdr3s.aa Character A VJ CDR3s amino acid sequence to search for
#' @param projects.to.search Optional character vector. Defaults to "ALL". Names of projects to search within.
#' @return A list of subsets of VDJ matrices from projects containing the query VDJ CDR3 (out[[1]]), the VJ CDR3 (out[[2]]) and cells containing both the query VDJ and VJ CDR3s (out[[3]])
#' @export
#' @examples
#' \dontrun{
#' public_clones <- PlatypusDB_find_CDR3s <- function(VDJ.cdr3s.aa = "CMRYGNYWYFDVW"
#' , VJ.cdr3s.aa = "CLQHGESPFTF", projects.to.search = "ALL")
#' }
#'
PlatypusDB_find_CDR3s <- function(VDJ.cdr3s.aa,
                                  VJ.cdr3s.aa,
                                  projects.to.search){
  platypus_url_lookup <- NULL


  if(missing(VDJ.cdr3s.aa)) VDJ.cdr3s.aa <- c("none")
  if(missing(VJ.cdr3s.aa)) VJ.cdr3s.aa <- c("none")
  if(missing(projects.to.search)) projects.to.search <- "ALL"

  print("Getting lookup table...")
  #Get the lookup table
  tryCatch({
  load(url("https://storage.googleapis.com/platypusdb_lookup/platypus_url_lookup.RData"))
    #FOR DEV
    #platypusdb_lookup <- new_lookup

    }, error=function(e){
    print(e)
    print(paste0("Failed to load lookup table. Please verify internet connection"))})

  platypusdb_lookup <- platypus_url_lookup

  print("Got lookup table")

  if(projects.to.search == "ALL"){
    to_download <- paste0(platypusdb_lookup$project_id, "__CDR3s")
  } else {
    to_download <- paste0(platypusdb_lookup$project_id[which(platypusdb_lookup$project_id %in% projects.to.search)], "__CDR3s")
  }

  out.list <- list()
  for(i in 1:length(to_download)){
    url_ <- paste0("https://storage.googleapis.com/platypusdb_cdr3s/", to_download[i], ".RData")
    load(url(url_), envir = .GlobalEnv)
    out.list[[i]] <- get(to_download[i])
    out.list[[i]]$project_id <- stringr::str_split(to_download[i], "__", simplify = T)[,1]
    rm(list = ls(pattern = to_download[i]), envir = .GlobalEnv)
    closeAllConnections()
  }

  out.list <- do.call(rbind, out.list)

  res <- list()

  out.list$VDJ_cdr3s_aa[out.list$VDJ_cdr3s_aa == ""] <- "NOCDR3"
  out.list$VJ_cdr3s_aa[out.list$VJ_cdr3s_aa == ""] <- "NOCDR3"

  for(i in 1:length(VDJ.cdr3s.aa)){
  res[[1]] <- subset(out.list, stringr::str_detect(VDJ.cdr3s.aa[i], out.list$VDJ_cdr3s_aa))
  print(paste0("Query VDJ CDR3 ", i  ," was found a total of ", nrow(res[[1]]), " times in the selected projects"))
  }

  for(i in 1:length(VJ.cdr3s.aa)){
  res[[2]] <- subset(out.list, stringr::str_detect(VJ.cdr3s.aa[i], out.list$VJ_cdr3s_aa))
  print(paste0("Query VJ CDR3 ", i ," was found a total of ", nrow(res[[2]]), " times in the selected projects"))
  }

  if(length(VDJ.cdr3s.aa) == length(VJ.cdr3s.aa)){
    res.clones <- list()
    for(i in 1:length(VDJ.cdr3s.aa)){
      res.clones[[i]] <- subset(out.list, stringr::str_detect(VJ.cdr3s.aa[i], out.list$VJ_cdr3s_aa) & stringr::str_detect(VDJ.cdr3s.aa[i], out.list$VDJ_cdr3s_aa))
      res.clones$query_clone <- i
      print(paste0("Query clone ", paste0(VJ.cdr3s.aa[i],"/",VDJ.cdr3s.aa[i]), "was found a total of ", nrow(res.clones[[i]]), " times in the selected projects"))
    }
    res[[3]] <- do.call(rbind, res.clones)
    names(res) <- c("VDJ CDR3 occurences", "VJ CDR3 occurences", "Clone occurences")
  } else {
    print("Not searching for clones (VDJ and VJ CDR3s), because input vectors are of different lengths")
    names(res) <- c("VDJ CDR3 occurences")
  }
  print("Done")
  return(res)
}
