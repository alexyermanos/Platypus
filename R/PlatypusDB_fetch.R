#' Loads and saves RData objects from the PlatypusDB
#' @param PlatypusDB.links Character vector. One or more links to files in the PlatypusDB. Links are constructed as follows: "\%Project id\%/\%sample_id\%/\%filetype\%". Any of the three can be "ALL", to download all files fitting the other link elements. If \%filetype\% is GEXmatrix or VDJmatrix or metadata, \%sample_id\% needs to be "", as these are elements which are not divided by sample. See examples for clarification. See last example on how to download AIRR compliant data. Feature Barcode (FB) data will be downloaded both for GEX and VDJ if present and does not need to be specified in the path. For sample_id entries the the metadata table for a given project via the function PlatypusDB_list_projects()
#'@param save.to.disk Boolean. Defaults to FALSE. Whether to save downloaded files individually to the directory specified in path.to.save
#'@param load.to.enviroment Boolean. Defaults to TRUE. Whether to load objects directly into the current .GlobalEnv. An array of the names of the loaded objects will be returned. !Be aware of RAM limitations of your machine when downloading multiple large files.
#'@param load.to.list Boolean. Defaults to FALSE. Whether to return loaded objects as a list. !Be aware of RAM limitations of your machine when downloading multiple large files.
#'@param path.to.save System path to save files to.
#'@param combine.objects Boolean. Defaults to TRUE. Whether to combine objects if appropriate. e.g. VDJ and GEX RData objects for a sample are saved as two independent objects and downloaded as such, to allow for flexibility. If combine.objects is set to TRUE, the function will coerce RData objects of each loaded sample or of each loaded VDJ_GEX_matrix appropriately. Combined input of VDJ and GEX Rdata objects can be directly supplied to the VDJ_GEX_matrix function.
#' @return A list of loaded project files as R objects if load.to.list = T or a name of these object loaded to the enviroment if load.to.enviroment = T.
#' @export
#' @examples
#' \dontrun{
#'
#'#Get a list of available projects by name
#'names(PlatypusDB_list_projects())
#'
#' #Load the VDJ_GEX_matrix of a project as an object and
#' #also save it to disk for later.
#' #This will download the VDJ and GEX part of the VDJ_GEX_matrix and combine
#' PlatypusDB_fetch(PlatypusDB.links = c("Kuhn2021a//ALL")
#' ,save.to.disk = FALSE,load.to.enviroment = TRUE, load.to.list = FALSE
#' , combine.object = TRUE,path.to.save = "/Downloads")
#'
#' #Load VDJ dataframe of the VDJ GEX matrix for all samples of one project
#' loaded_list <- PlatypusDB_fetch(PlatypusDB.links = c("Kuhn2021a//VDJmatrix")
#' ,save.to.disk = FALSE,load.to.enviroment = FALSE, load.to.list = TRUE)
#'
#' #Load the VDJ and GEX RData of 2 samples from
#' #2 different projects which can be directly passed
#' #on to the VDJ_GEX_matrix function to integrate
#'#downloaded_objects <- PlatypusDB_fetch(
#'#PlatypusDB.links = c("Project1/s1/ALL", "Project1/s2/ALL")
#'#,save.to.disk = FALSE,load.to.enviroment = FALSE, load.to.list = TRUE
#'#, combine.objects = TRUE)
#'
#' #integrated_samples <- VDJ_GEX_matrix_DB(data.in = downloaded_objects)
#'
#' #Download metadata objects for projects
#' list_of_metadata_tables <- PlatypusDB_fetch(
#' PlatypusDB.links = c("Kuhn2021a//metadata")
#' ,save.to.disk = FALSE,load.to.enviroment = FALSE, load.to.list = TRUE)
#'
#' #Dowload of airr_rearrangement.tsv
#' #Load VDJ.RData into a list
#'#downloaded_objects <- PlatypusDB_fetch(
#'#PlatypusDB.links = c("Project1/ALL/VDJ.RData"),save.to.disk = FALSE
#'#,load.to.enviroment = FALSE, load.to.list = TRUE)
#'
#' #Extract airr_rearrangement table for sample 1
#' #airr_rearrangement <- downloaded_objects[[1]][[1]][[6]]
#' #Index hierarchy: Sample, VDJ or GEX, VDJ element
#'
#' #Save for import to AIRR compatible pipeline
#' #write.table(airr_rearrangement, file = "airr_rearrangement_s1.tsv", sep='\t',
#' #row.names = FALSE, quote=FALSE)
#'
#' }
#'
PlatypusDB_fetch <- function(PlatypusDB.links,
                             save.to.disk,
                             load.to.enviroment,
                             load.to.list,
                             path.to.save,
                             combine.objects){

  project_id <- NULL
  sample_id <- NULL
  filetype <- NULL
  group <- NULL
  platypus_url_lookup <- NULL


  if(missing(save.to.disk)) save.to.disk <- F
  if(missing(load.to.enviroment)) load.to.enviroment <- T
  if(missing(load.to.list)) load.to.list <- F
  if(missing(combine.objects)) combine.objects <- T

  #if not loaded into R the files needs to go somewhere
  if(load.to.enviroment == F & load.to.list == F) save.to.disk <- T

  if(load.to.enviroment == T & load.to.list == T){
    warning("Downloaded objects will be loaded to .GlobalEnv and additionally returned as list. This doubles the storage usage. Having both load.to.enviroment and load.to.list set to TRUE is not recommended")}

  if(save.to.disk == T & missing(path.to.save)){
    message("No path.to.save specified. Saving files to current working directory")
    path.to.save <- getwd()
  }
  if(save.to.disk == F & missing(path.to.save)){
    path.to.save <- "none"}
  #making sure we have a / at the end of the path
  if(substr(path.to.save, nchar(path.to.save), nchar(path.to.save)) != "/"){
    path.to.save <- paste0(path.to.save, "/")
  }

  if(save.to.disk == T){
    combine.objects <- F
    load.to.list <- F
    load.to.enviroment <- F
  }

  for(i in 1:length(PlatypusDB.links)){

    if(stringr::str_detect(PlatypusDB.links[i], "\\.RData") == F){
      PlatypusDB.links[i] <- paste0(PlatypusDB.links[i], ".RData")
    }

    if(stringr::str_detect(PlatypusDB.links[i], "\\.zip")){stop(".zip files are not downloadable via R. Please refer to the PlatypusDB website")}

    if(stringr::str_detect(PlatypusDB.links[i], ".*/.*/.*\\.RData") == F & stringr::str_detect(PlatypusDB.links[i], ".*//.*\\.RData") == F){
      stop(paste0("Input link ", i, " does not fit the neccessary format. Please input a link with the format '[Project id]/[Sample id]/[filetype].RData' or '[Project id]//[filetype].RData' to download VDG_GEX_matrix or project metadata files. To download .zip files please refer to the PlatypusDB website. Any link element can be replaced with ALL or ALL.RData, for the final element, to download all files of that catagory"))
    }
  }

  #Get the lookup table
  tryCatch({
    load(url("https://storage.googleapis.com/platypusdb_lookup/platypus_url_lookup.RData"))
    #FOR DEV
    #platypusdb_lookup <- new_lookup

    }, error=function(e){
    message(paste0("Failed to load lookup table. Please verify internet connection \n", e))})

  platypusdb_lookup <- platypus_url_lookup #Reassignment

  #FORDEV
  #lookup <- data.frame("project_id" = c("A","A","B","B","B","C","C"), "sample_id" = c("s1","s1","","","","s1","s1"), "filetype" = c("VDJ.RData", "GEX.RData", "gexVGM.RData", "vdjVGM.RData", "metadata.RData", "VDJ.RData","gex.RData"), "url" = c("AA","AB","BA","BB","BC"))

  #add group column for combining to lookup
  platypusdb_lookup$group <- 1:nrow(platypusdb_lookup)
  gcount <- nrow(platypusdb_lookup)+1

  #make sense of the present link
  link_parts <- as.data.frame(stringr::str_split(PlatypusDB.links, "/", simplify = T))

  link_parts$full <- PlatypusDB.links #This is link_parts[,4]

  #check for patterns in more than one link that could be combined as files
  if(nrow(link_parts) > 1){
    if(any(duplicated(paste0(link_parts[,1], link_parts[,2])))){ #input like c("agedcnsbcells/S1/VDJ.RData", "agedcnsbcells/S1/GEX.RData")
      stop("Detected two links pointing to the same sample. To download VDJ and GEX Rdata of one samples please use path 'project_id/[sample_id or ALL]/ALL'. To download a processed VDJ GEX matrix please use path 'project_id//ALL'")
    }
  }


  if(any(!link_parts[,1] %in% unique(platypusdb_lookup$project_id))){ #Stop if project ids not found in lookup
    stop("One or more selected project ids were not present on the PlatypusDB. Please revise project id inputs")
  }

  #iterate over links
  to_download_list <- list()
  for(i in 1:nrow(link_parts)){
    #look for patterns were objects could be combined
    if(stringr::str_detect(link_parts[i,4], "/ALL.RData") & link_parts[i,2] != "" & link_parts[i,2] != "ALL"){
      to_download <- subset(platypusdb_lookup, project_id == link_parts[i,1] & sample_id == link_parts[i,2] & stringr::str_detect(name, ".zip") == F)
      if(nrow(to_download) == 2){
      to_download$group <- gcount
      gcount <- gcount + 1
      to_download_list[[i]] <- to_download
      } else{ #do not group if unexpectately there are more of different entries present
      to_download_list[[i]] <- to_download
      }

    } else if(stringr::str_detect(link_parts[i,4], "ALL/ALL.RData")){
      to_download <- subset(platypusdb_lookup, project_id == link_parts[i,1] & sample_id != "" & stringr::str_detect(name, ".zip") == F)

      for(j in 1:length(unique(to_download$sample_id))){
        to_download_s <- subset(to_download, sample_id == unique(to_download$sample_id)[j])

        if(all(c("VDJ.RData", "GEX.RData") %in% to_download_s$filetype) & nrow(to_download_s == 2)){
        to_download_s$group <- gcount
        gcount <- gcount + 1
        to_download_list[[length(to_download_list)+1]] <- to_download_s

        } else{
        to_download_list[[length(to_download_list)+1]] <- to_download_s
        }
      }

    } else if(stringr::str_detect(link_parts[i,4], "//ALL.RData")){
      to_download <- subset(platypusdb_lookup, project_id == link_parts[i,1] & filetype %in% c("VDJmatrix.RData", "GEXmatrix.RData"))
      if(nrow(to_download) == 2){

        to_download$group <- gcount
        gcount <- gcount + 1

        to_download_list[[i]] <- to_download
      } else {
        stop("More than 2 matrix files detected for this project #174")
      }
    } else {

    #subset hierarchically
    if(link_parts[i,1] == "ALL"){
      to_download <- platypusdb_lookup
    }else{
      to_download <- subset(platypusdb_lookup, project_id == link_parts[i,1])
    }
    if(link_parts[i,2] == "ALL"){
    } else { to_download <- subset(to_download, sample_id == link_parts[i,2])
    }
    if(link_parts[i,3] == "ALL"){
    }else{
      to_download <- subset(to_download, filetype == link_parts[i,3])
      to_download$group <- gcount
      gcount <- gcount + 1
    }
      to_download_list[[i]] <- to_download
    }

    to_download$size <- as.numeric(stringr::str_extract(to_download$size, "\\d+"))

  }


  to_download <- do.call(rbind, to_download_list) #combine into dataframe

  to_download <- to_download[duplicated(to_download$url) == F,] #make unique in case one file is entered more than once
  to_download$group <- as.character(to_download$group) #make character so that the order of rows is not disturbed by the unique() in the following

  if(nrow(to_download) == 0) stop("No files found matching PlatypusDB.links. Please revise input")

  out.list <- list()
  if(combine.objects == F){ #IF NO OBJECTS ARE TO BE COMBINED. group column will not be used
  for(i in 1:nrow(to_download)){
    tryCatch({

      curr_download_name <- gsub("\\.RData","",to_download$name[i]) #have a name ready to use for objects in the r enviroment without the .RData extension

      message(paste0(Sys.time(), ": Starting download of ", to_download$name[i],"..."))

      if(save.to.disk == T){ #if objects are to be saved to disk, this happens here.
        utils::download.file(to_download$url[i], destfile = paste0(path.to.save,curr_download_name, ".RData")) #Saving directly to disk to avoid RAM usage

      } else { #Save to disk == F

      load(url(to_download$url[i]), envir = .GlobalEnv) #download and load to global enviroment

      if(load.to.list == T){
        out.list[[i]] <- get(curr_download_name) #if a list of objects is to be returned, add the just loaded object into a list
        names(out.list)[i] <- curr_download_name
      } else {
        out.list[[i]] <- curr_download_name #if not list of objects is to be returned, only the name of the loaded object is appended and will be returned for reference. (see also last few lines of the function)
      }

      if(load.to.enviroment == F){ #if objects should not present in the global enviroment at the end of the function run, we delete them from the .GlobalEnv again.
        rm(list = ls(pattern = curr_download_name, envir = .GlobalEnv), envir = .GlobalEnv)
      }
      }

    }, error=function(e){
      message(paste0("Failed to load",  to_download$url[i], "\n", e))})
  }


  } else if(combine.objects == T ){ #IF COMBINING: iterating over download groups

    for(j in 1:length(unique(to_download$group))){
      curr_to_download <- subset(to_download, group == unique(to_download$group)[j])
      #get the order right => VDJ first, GEX second
      if(nrow(curr_to_download) == 2 & (curr_to_download$filetype[1] %in% c("VDJ.RData", "GEX.RData") | curr_to_download$filetype[1] %in% c("VDJmatrix.RData", "GEXmatrix.RData"))){
       curr_to_download <- curr_to_download[c(2,1),] #swap rows
      }

      tryCatch({

      if(nrow(curr_to_download) == 2 & (curr_to_download$filetype[1] %in% c("VDJ.RData", "GEX.RData")| curr_to_download$filetype[1] %in% c("VDJmatrix.RData", "GEXmatrix.RData"))){
        cur.out.list <- list()

        curr_download_names <- gsub("\\.RData","",curr_to_download$name) #getting the names of the objects. [1] is VDJ [2] is GEX

        for(i in 1:nrow(curr_to_download)){ #load both objects!
        message(paste0(Sys.time(), ": Starting download of ", curr_to_download$name[i],"..."))
        load(url(curr_to_download$url[i]), envir = .GlobalEnv) #download
        } #End of loop. Now both files which are to be combined (are part of one predefined group) should be present in the R enviroment

        #Combine the objects by parsing. This saves ram because we do not have to copy the GEX object into a new list
        cmd <- paste0(curr_download_names[2],"[[1]] <- ", curr_download_names[1], "[[1]]")
        try(eval(parse(text=cmd)))

        #remove the VDJ object
        rm(list = ls(pattern = curr_download_names[1], envir = .GlobalEnv), envir = .GlobalEnv)

        #rename the R object
        if(stringr::str_detect(curr_to_download$filetype[1], "matrix")){ #rename if a full VGM is being downloaded
          comb_download_name <- paste0(curr_to_download$project_id[1], "__VDJGEXmatrix")
          eval(parse(text = "assign(comb_download_name, get(curr_download_names[2]), envir = .GlobalEnv)"))
          rm(list = ls(pattern = curr_download_names[2], envir = .GlobalEnv), envir = .GlobalEnv)
        } else { #rename if RData is being downloaded
          comb_download_name <- paste0(curr_to_download$project_id[1], "_", curr_to_download$sample_id[1] , "_VDJGEXdata")
          eval(parse(text = "assign(comb_download_name, get(curr_download_names[2]), envir = .GlobalEnv)"))
          rm(list = ls(pattern = curr_download_names[2], envir = .GlobalEnv), envir = .GlobalEnv)
        }

        #We now have the combined download objects named comb_download_name in the global enviroment. The other objects have been deleted.
        #Therefore we can proceed with the returning / saving modes as above.

        if(load.to.list == T){
          out.list[[j]] <- get(comb_download_name)
          names(out.list)[j] <- comb_download_name
        } else {
          out.list[[j]] <- comb_download_name
        }

        if(save.to.disk == T){
          utils::download.file(to_download$url[i], destfile = paste0(path.to.save,curr_download_name, ".RData")) #Saving directly to disk to avoid RAM usage
        }

        if(load.to.enviroment == F){
          rm(list = ls(pattern = comb_download_name, envir = .GlobalEnv), envir = .GlobalEnv)
        }

        #This else if is there to catch an instance were the input PlatypusDB paths asked for both two files which should be combined, as well as extra files which are not to be combined (e.g. metadata). Here we can proceed as if combine.objects == F
      } else if(nrow(curr_to_download) == 1){ #if the group number of that entry was unique

        curr_download_name <- gsub("\\.RData","",curr_to_download$name[1])

        message(paste0(Sys.time(), ": Starting download of ", curr_to_download$name[1],"..."))
        load(url(curr_to_download$url[1]), envir = .GlobalEnv) #download and load to global enviroment

        if(load.to.list == T){
          out.list[[j]] <- get(curr_download_name) #if a list of objects is to be returned, add the just loaded object into a list
          names(out.list)[j] <- curr_download_name
        } else {
          out.list[[j]] <- curr_download_name
        }

        if(save.to.disk == T){
          utils::download.file(to_download$url[i], destfile = paste0(path.to.save,curr_download_name, ".RData")) #Saving directly to disk to avoid RAM usage
        }

        if(load.to.enviroment == F){
          rm(list = ls(pattern = curr_download_name, envir = .GlobalEnv), envir = .GlobalEnv)
        }
        #end of else if(nrow(curr_to_download) == 1)
      } else if(nrow(curr_to_download) > 2){ #looks like wrong grouping or no grouping
        message("Downloads of different samples are not combined. Please set load.to.list = T to return a single list containing info of all downloaded samples which can be used as input to the VDJ_GEX_matrix function")

        for(k in 1:nrow(to_download)){
          tryCatch({

            curr_download_name <- gsub("\\.RData","",to_download$name[k]) #have a name ready to use for objects in the r enviroment without the .RData extension

            message(paste0(Sys.time(), ": Starting download of ", to_download$name[k],"..."))

            if(save.to.disk == T){ #if objects are to be saved to disk, this happens here.
              utils::download.file(to_download$url[k], destfile = paste0(path.to.save,curr_download_name, ".RData")) #Saving directly to disk to avoid RAM usage

            } else { #Save to disk == F

              load(url(to_download$url[k]), envir = .GlobalEnv) #download and load to global enviroment

              if(load.to.list == T){
                out.list[[k]] <- get(curr_download_name) #if a list of objects is to be returned, add the just loaded object into a list
                names(out.list)[k] <- curr_download_name
              } else {
                out.list[[k]] <- curr_download_name #if not list of objects is to be returned, only the name of the loaded object is appended and will be returned for reference. (see also last few lines of the function)
              }

              if(load.to.enviroment == F){ #if objects should not present in the global enviroment at the end of the function run, we delete them from the .GlobalEnv again.
                rm(list = ls(pattern = curr_download_name, envir = .GlobalEnv), envir = .GlobalEnv)
              }
            }

          }, error=function(e){
            message(paste0("Failed to load",  to_download$url[i], "\n", e))})


        }
      }
      }, error=function(e){
        message(e)
        message(paste0("Failed to load",  names(out.list)[j]))})
    } #end of loop over files
  } #end of if(combine.object == T)

  if(load.to.list == F){
    return(unlist(out.list)) #unlist to return a simple array of names of loaded files
  } else{
    return(out.list) #return loaded files
  }
}
