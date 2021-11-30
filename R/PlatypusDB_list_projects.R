#' Loads and saves RData objects from the PlatypusDB
#' @param keyword Character. Keyword by which to search project ids (First Author, Year) in the database. Defaults to an empty string ("") which will list all projects currently available
#' @return A list of metadata tables by project. List element names correspond to project ids to use in the PlatypusDB_fetch function
#' @export
#' @examples
#' \dontrun{
#'
#' #Get list of all available projects and metadata.
#' PlatypusDB_projects <- PlatypusDB_list_projects()
#'
#' #Names of list are project ids to use in PlatypusDB_fetch function
#' names(PlatpyusDB_projects)
#' #Common format: first author, date, letter a-z (all lowercase)
#'
#' View metadata of a specific project
#' View(PlatypusDB_projects[["Yermanos2021a"]])
#'
#' }
#'
PlatypusDB_list_projects <- function(keyword){

  platypus_url_lookup <- NULL

  if(missing(keyword)) keyword <- ""

  if(keyword == ""){
    print("Fetching list of all projects...")
  } else {
    print(paste0("Searching for keyword ", keyword, " and returning list of hits"))
  }

  #Get the lookup table
  tryCatch({
    load(url("https://storage.googleapis.com/platypusdb_lookup/platypus_url_lookup.RData"))
    #FOR DEV
    #platypusdb_lookup <- new_lookup

    }, error=function(e){
    print(e)
    print(paste0("Failed to load lookup table. Please verify internet connection"))})

  platypusdb_lookup <- platypus_url_lookup #Reassignment

  platypusdb_meta <- subset(platypusdb_lookup, stringr::str_detect(platypusdb_lookup$filetype, "metadata"))

  if(keyword != ""){
    platypusdb_meta <- subset(platypusdb_meta, stringr::str_detect(platypusdb_meta$name, keyword))
    print(paste0("Found ", nrow(platypusdb_meta), " project ids containing keyword"))
    if(nrow(platypusdb_meta) == 0){
      stop("Please retry with a different keyword or provide an empty string to the keyword argument to list all projects")
    }
  }

  print("Downloading metadata list...")
  out.list <- list()
  for(i in 1:nrow(platypusdb_meta)){
    tryCatch({

      curr_download_name <- gsub("\\.RData","",platypusdb_meta$name[i]) #have a name ready to use for objects in the r enviroment without the .RData extension

      load(url(platypusdb_meta$url[i]), envir = .GlobalEnv) #download and load to global enviroment

        out.list[[i]] <- get(curr_download_name) #add the just loaded object into a list
        names(out.list)[i] <- platypusdb_meta$project_id[i]
        rm(list = ls(pattern = curr_download_name, envir = .GlobalEnv), envir = .GlobalEnv)

    }, error=function(e){
      print(e)
      print(paste0("Failed to load",  platypusdb_meta$url[i]))})
  }

    print("Done")
    return(out.list) #return loaded files
  }

