#' Load and preprocess a list of antigen-specific databases

#'@description Preprocessing function for several antigen databases for both TCRs (VDJdb, McPAS-TCR, TBAdb) and BCRs (TBAdb), saving them either at a specified path, or loading them as a database list for downstream integration/analyses.
#' @param databases list of databases to be processed and saved. Currently supported ones include: VDJdb(='vdjdb'), McPAS-TCR(='mcpas'), TBAdb(='tbdadb_tcr' or 'tbadb_bcr').
#' @param file.paths list of file paths for the specified databases (in the database parameter). If NULL, will try to locally download the databases from the archived download links.
#' @param preprocess boolean - if T, will preprocess each database individually.
#' @param species string - either 'Human' or 'Mouse', the species for the processed database. Needs preprocess=T.
#' @param filter.sequences string - 'VDJ' to remove rows with NA VDJ sequences, 'VJ' to remove rows with NA VJ sequences, 'VDJ.VJ' to remove rows with both VDJ and VJ sequences missing. Needs preprocess=T.
#' @param remove.na string or NULL - 'all' will remove all rows with missing values from the database, 'common' will remove only rows with missing values for the shared columns among all databases ('VJ_cdr3s_aa','VDJ_cdr3s_aa','Species','Epitope','Antigen species'), 'vgm' will remove missing values for columns shared with the VDJ object (specific to each database). Needs preprocess=T.
#' @param vgm.names boolean - if T, will change all column names of the shared columns (with VDJ) to match those from VDJ. Use this to integrate the antigen data into VDJ using VDJ_antigen_integrate or VDJ_db_annotate. Needs preprocess=T.
#' @param keep.only.common boolean - if T, will only keep the columns shared between all databases ('VJ_cdr3s_aa','VDJ_cdr3s_aa','Species','Epitope','Antigen species') for each processed database. Needs preprocess=T.
#' @param output.format string - 'df.list' to save all databases as a list, 'save' to save them as csv files.
#' @param saving.path string - directory where the processed databases should be locally saved if output.format='save'.

#' @return Processed antigen-specific databases for both TCRs and BCRs.
#' @export

#' @examples
#' \dontrun{
#' VDJ_db_load(databases=list('vdjdb'),file.paths=NULL,
#' preprocess=TRUE,species='Mouse',filter.sequences='VDJ.VJ',
#' remove.na='vgm', vgm.names=TRUE, keep.only.common=TRUE,
#' output.format='df.list')
#'}

VDJ_db_load <- function(databases,
                        file.paths,
                        preprocess,
                        species,
                        filter.sequences,
                        remove.na,
                        vgm.names,
                        keep.only.common,
                        output.format,
                        saving.path
                        ){


  #### Internal function definition ####

  .single_db_preprocess <- function(database.df,
                                    database,
                                    filter.species,
                                    filter.na.sequences,#all, only VDJ, only VJ
                                    use.vgm.names,
                                    filter.na.features,
                                    keep.only.common.cols
  ){

    if(missing(database.df)) stop('Input your loaded database')
    if(missing(database)) database <- 'vdjdb'
    if(missing(filter.species)) filter.species <- NULL
    if(missing(filter.na.sequences)) filter.na.sequences <- 'VDJ.VJ'
    if(missing(use.vgm.names)) use.vgm.names <- T
    if(missing(filter.na.features)) filter.na.features <- NULL
    if(missing(keep.only.common.cols)) keep.only.common.cols <- T

    output_db <- database.df


    vdjdb_common_cols <- list('cdr3.alpha', 'cdr3.beta', 'v.alpha', 'j.alpha', 'v.beta', 'd.beta', 'j.beta', 'species', 'antigen.epitope', 'antigen.gene', 'antigen.species', 'reference.id')
    names(vdjdb_common_cols) <- list('VJ_cdr3s_aa', 'VDJ_cdr3s_aa', 'VJ_vgne', 'VJ_jgene', 'VDJ_vgene', 'VDJ_dgene', 'VDJ_jgene', 'Species', 'Epitope', 'Antigen gene', 'Antigen species', 'Antigen ID')

    mcpas_common_cols <- list('CDR3.alpha.aa', 'CDR3.beta.aa', 'Species', 'Pathology', 'Epitope.peptide', 'Protein.ID')
    names(mcpas_common_cols) <- list('VJ_cdr3s_aa', 'VDJ_cdr3s_aa', 'Species', 'Antigen species', 'Epitope', 'Antigen ID')

    tbadb_bcr_common_cols <- list('CDR3.light.aa', 'CDR3.heavy.aa', 'Vlight', 'Jlight', 'Vheavy', 'Dheavy', 'Jheavy', 'Species', 'Antigen', 'Antigen.sequence')
    names(tbadb_bcr_common_cols) <- list('VJ_cdr3s_aa', 'VDJ_cdr3s_aa', 'VJ_vgene', 'VJ_jgene', 'VDJ_vgene', 'VDJ_dgene', 'VDJ_jgene', 'Species', 'Antigen species', 'Epitope')

    tbadb_tcr_common_cols <- list('CDR3.alpha.aa', 'CDR3.beta.aa', 'Valpha', 'Jalpha', 'Vbeta', 'Dbeta', 'Jbeta', 'Species', 'Antigen', 'Antigen.sequence')
    names(tbadb_tcr_common_cols) <- list('VJ_cdr3s_aa', 'VDJ_cdr3s_aa', 'VJ_vgene', 'VJ_jgene', 'VDJ_vgene', 'VDJ_dgene', 'VDJ_jgene', 'Species', 'Antigen species', 'Epitope')

    names <-  list(names(vdjdb_common_cols), names(mcpas_common_cols), names(tbadb_bcr_common_cols), names(tbadb_tcr_common_cols))
    common_names <- Reduce(intersect, names)
    if(database=='vdjdb'){
      db_names <- vdjdb_common_cols
    }else if(database=='mcpas'){
      db_names <- mcpas_common_cols
    }else if(database=='tbadb_bcr'){
      db_names <- tbadb_bcr_common_cols
    }else if(database=='tbdadb_tcr'){
      db_names <- tbadb_tcr_common_cols
    }else{
      stop('Database not implemented yet!')
    }

    if(('homosapiens' %in% tolower(output_db[,db_names$Species])) | ('homo sapiens' %in% tolower(output_db[,db_names$Species]))){
      output_db[,db_names$Species][which(tolower(output_db[,db_names$Species])== 'homosapiens')] <- 'Human'
      output_db[,db_names$Species][which(tolower(output_db[,db_names$Species])== 'homo sapiens')] <- 'Human'
      output_db[,db_names$Species][which(tolower(output_db[,db_names$Species])== 'musmusculus')] <- 'Mouse'
      output_db[,db_names$Species][which(tolower(output_db[,db_names$Species])== 'mus musculus')] <- 'Mouse'
    }

    if(!is.null(filter.species)){
      output_db <- output_db[which(output_db[db_names$Species]==filter.species),]
    }

    if(!is.null(filter.na.sequences)){
      if(filter.na.sequences=='VDJ'){
        output_db <- output_db[which(!is.na(output_db[,db_names$VDJ_cdr3s_aa])),]
      }else if(filter.na.sequences=='VJ'){
        output_db <- output_db[which(!is.na(output_db[,db_names$VJ_cdr3s_aa])),]
      }else if(filter.na.sequences=='VDJ.VJ'){
        output_db <- output_db[which(!is.na(output_db[,db_names$VJ_cdr3s_aa]) & !is.na(output_db[,db_names$VDJ_cdr3s_aa])),]
      }
    }

    if(!is.null(filter.na.features)){

      if(filter.na.features=='all'){
        output_db <- stats::na.omit(output_db)

      }else if(filter.na.features=='vgm'){
        db_names_temp <- db_names[names(db_names)!='VDJ_cdr3s_aa' & names(db_names)!='VJ_cdr3s_aa']
        for(feature in db_names_temp){
          output_db <- output_db[which(!is.na(output_db[,feature])),]
        }
      }else if(filter.na.features=='common'){
        db_names_temp <- common_names[names(common_names)!='VDJ_cdr3s_aa' & names(common_names)!='VJ_cdr3s_aa']
        for(feature in db_names_temp){
          output_db <- output_db[which(!is.na(output_db[,feature])),]
        }
      }else{
        for(feature in filter.na.features){
          output_db <- output_db[which(!is.na(output_db[,feature])),]
        }
      }
    }

    if(keep.only.common.cols){
      output_db <- dplyr::select(output_db, tidyselect::all_of(unlist(unname(db_names[common_names]))))
    }

    if(use.vgm.names){
      colnames(output_db)[colnames(output_db) %in% db_names] <- names(db_names)[db_names %in% colnames(output_db)]
    }

    return(output_db)
  }

  #### Main function body ####


  if(missing(databases)) databases <- list('vdjdb')
  if(missing(file.paths)) file.paths <- NULL
  if(missing(preprocess)) preprocess <- T
  if(missing(species)) species <- 'Mouse'
  if(missing(filter.sequences)) filter.sequences <- 'VDJ.VJ'
  if(missing(remove.na)) remove.na <- 'common'
  if(missing(vgm.names)) vgm.names <- T
  if(missing(keep.only.common)) keep.only.common <- T
  if(missing(output.format)) output.format <- 'df.list'
  if(missing(saving.path)) saving.path <- './data'


  db_list <- list('vdjdb', 'mcpas', 'tbadb_tcr', 'tbadb_bcr')

  vdjdb_master_url <- 'https://raw.githubusercontent.com/antigenomics/vdjdb-db/master/latest-version.txt'
  mcpas_url <- 'http://friedmanlab.weizmann.ac.il/McPAS-TCR/session/7f617e999fa2e89c384a41f5c68cbf5e/download/downloadDB?w='
  tbadb_url <- 'https://db.cngb.org/dc_assets/data/dc_pird/data_template/TBAdb.xlsx'

  db_output <- list()
  for(i in 1:length(databases)){
    if(databases[i]=='vdjdb'){
      if(!is.null(file.paths)){
        if(file.exists(file.paths[i])){
          db_output[[i]] <- utils::read.csv(file.paths[i], sep='\t', header=T)
        }

      }else{
        latest_vdjdb_url <- utils::read.table(vdjdb_master_url)[1,]
        temp <- tempfile()
        downloaded_file <- utils::download.file(latest_vdjdb_url, temp)
        db_output[[i]] <- utils::read.csv(unz(temp, 'vdjdb_full.txt'), sep='\t', header=T)
        unlink(temp)
      }

      db_output[[i]][db_output[[i]]==''] <- NA

      if(preprocess){
        db_output[[i]] <- .single_db_preprocess(db_output[[i]], database='vdjdb', filter.species=species, filter.na.sequences=filter.sequences, use.vgm.names=vgm.names, filter.na.features=remove.na, keep.only.common.cols=keep.only.common)
      }


    }else if(databases[i]=='mcpas'){
      if(!is.null(file.paths)){
        if(file.exists(file.paths[i])){
          db_output[[i]] <- utils::read.csv(file.paths[i], sep=',')
        }

      }else{
        db_output[[i]]<-utils::read.table(mcpas_url, sep=',')
        colnames(db_output[[i]]) <- db_output[[i]][1,]
        db_output[[i]] <- db_output[[i]][-1,]
      }


      if(preprocess){
        db_output[[i]] <- .single_db_preprocess(db_output[[i]], database='mcpas', filter.species=species, filter.na.sequences=filter.sequences, use.vgm.names=vgm.names, filter.na.features=remove.na, keep.only.common.cols=keep.only.common)
      }


    }else if(databases[i]=='tbadb_tcr'){

      if(!is.null(file.paths)){
        if(file.exists(file.paths[i])){
          db_output[[i]] <- readxl::read_excel(file.paths[i],sheet='TCR-AB')
        }

      }else{
        temp <- tempfile()
        downloaded_file <- utils::download.file(tbadb_url, temp)
        db_output[[i]] <- readxl::read_excel(temp, sheet='TCR-AB')
      }
      db_output[[i]][db_output[[i]]=='-'] <- NA

      if(preprocess){
        db_output[[i]] <- .single_db_preprocess(db_output[[i]], database='tbadb_tcr', filter.species=species, filter.na.sequences=filter.sequences, use.vgm.names=vgm.names, filter.na.features=remove.na, keep.only.common.cols=keep.only.common)
      }


    }else if(databases[i]=='tbadb_bcr'){

      if(!is.null(file.paths)){
        if(file.exists(file.paths[i])){
          db_output[[i]] <- readxl::read_excel(file.paths[i],sheet='BCR')
        }

      }else{
        temp <- tempfile()
        downloaded_file <- utils::download.file(tbadb_url, temp)
        db_output[[i]] <- readxl::read_excel(temp, sheet='BCR')
      }

      if(preprocess){
        db_output[[i]] <- .single_db_preprocess(db_output[[i]], database='tbadb_bcr', filter.species=species, filter.na.sequences=filter.sequences, use.vgm.names=vgm.names, filter.na.features=remove.na, keep.only.common.cols=keep.only.common)
      }
      db_output[[i]][db_output[[i]]=='-'] <- NA

    }else if(!(databases[i] %in% db_list)){
      warning(paste0('Database ', databases[i], ' not yet supported!'))
      db_output[[i]] <- 'empty' #Because you can't name NULLs........................
    }
  }

  names(db_output) <- databases
  db_output <- db_output[db_output!='empty']


  if(output.format=='df.list'){
    return(db_output)
  }else{
    names <- names(db_output)
    for(i in 1:length(db_output)){
      save = paste0(saving.path, '/', names[i], '.csv')

      utils::write.csv(db_output[[i]], file=save)
    }
  }
}

