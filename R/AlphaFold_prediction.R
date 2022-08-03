#'Structure prediction of Mixcr wrapper output with Alpha Fold
#'
#'
#'@description This function takes the output from the VDJ_call_MIXCR function as input in the VDJ.mixcr.out argument and predicts the structure with Alpha Fold.
#'From the VDJ.mixcr.out object the full length VDJ & VJ sequence containing all the frameworks and CDR's is used to predict the structure of
#'the variable part with Alpha Fold multi. If the user has no access to the Euler function, the fucntion just returns a fasta file with the VDJ and VJ sequence,
#'that can be used for running Alpha Fold on a Cluster. 
#'For users that have a login to the Euler cluster, this function will automatically connect to Euler and start Alpha Fold for all the indicated sequences.
#'After the prediction is finished the same function can be used to import the predicted structures as a pdb file and add it to the input as a list object
#'
#'@note For running Alpha Fold on Euler, the user needs to have access to GPU usage. This is automatically activated if one is part of the Reddy Euler Group.
#'
#'@param VDJ.mixcr.out Contains the output from the VDJ_call_MIXCR function with VJ_aa_mixcr and VDJ_aa_mixcr columns containing the full length amino acid sequence of Framework 1 - 4.
#'@param cells.to.predict Here you can specify 10x barcodes for the cells of the VDJ.mixcr.out that should be used for structure prediction. It can be set to "ALL" if the antibody structure of all cells shall be predicted. 
#'@param max.template.date This is a parameter for running Alpha Fold and a date can be specified in the following format: "yyyy-mm-dd" This tells Aplhpa Fold which state of the databases it shall use.
#'@param fasta.storage.path Here you can specify where the function saves the fasta files needed as an Alpha Fold input. By default files an 'AlphaFold_Fasta' directory with all the fasta files is created in the same directory as the R script runs.
#'@param dir.name By default the function creates a directory named AlphaFold_Fasta, where the FASTA files created for prediction. The name of this directory can be changed by specifying the dir.name argument.
#'@param euler.user.name If running Alpha Fold on Euler is requested, the user name needs to be specified in this parameter. Make sure that you have access to GPU usage on the Cluster. You will be prompted to enter you password by the "ssh" package which handles your credentials in a safe manner.
#'@note If running prediction on Euler, the function will create a "AlphaFold_Fasta" directory in sour scratch on the cluster where all the fasta files are uploaded. The output files will be saved as well in this directory.
#'@param import This argument is for telling the function to import predicted structures. It is by default set to FALSE, which will initiate prediction not import. There are two options for importing predicted structures: Import = "euler" will start a connection to Euler and import the pdb files from the "AlphaFold_Fasta/output" directory. Import = "local" will import the pdb files form a local directory.
#'@param euler.dirname If import = "euler" is used the name of the directory containing the Alpha Fold output directory can be specified in euler.dirname. It is set to "AlphaFold_Fasta" by default and is expected to be on your scratch.
#'@param euler.path If import = "euler" is used the path to the directory containing the Alpha Fold output folder can be specified in euler.path. By default the function expects the output in the AlphaFold_Fasta directory on your scratch. In case you wanna import the data from a different location you can specify the path here. The function expects a sub directory named output which contains sub directories named after the specific barcodes. (../scratch/AlphaFold_Fasta/output/s4_AGCCTAATCCCTTGCA-1/) 
#'@param import.local.dirnames If import = "local" is used the function expects a directory named 'Output_AlphaFold' in the same directory as the script runs. In case you do not wanna import all the pdb files off all samples in the 'Output_AlphaFold directory you can specify a sub directories in the import.local.dirnames parameter. (import.local.dirnames = c(s4_AGCCTAATCCCTTGCA-1_ranked,s4_CCCATACCACGTTGGC-1_ranked,...))
#'@param import.local.path If import = "local" is used you can specify the path to the AlphaFold_Output directory here. By default it is expected in the same directory as the r script runs.
#'@param n.ranked Alpha Fold returns 21 predictions for each sequence which are ranked for 0:20. The ranked_0 is the most accurate according to the model. Here you can specify how many of the top ranked structures are added to the output object. By default only the most accurate structure 'ranked_0.pdb' is integrated.
#'@param platypus.version This function is not directly depended on other Platypus functions but was developed to be compatible with v3.
#'@param rm.local.fasta Here you can specify if the local AlphaFold_Fasta directory shall be deleted from your local computer after uploaded to the scratch on Euler. By default it is set to TRUE, to keep your environment clean. If the function is not used in the Euler modus it is set to FALSE, so you will have the fasta files as an output.
#'@param rm.euler.files Here you can specify if the files on Euler shall be deleted after importing them. It is set to FASLE by default to reduce the risk of unintentionally deleting the predictions. However, make sure to keep you scratch environment clean.
#'@param rm.local.output Here you an specify if the downloaded output folder from Euler shall be deleted after the import. It is set to true by default to keep you environment clean.
#'@param output.path If the data is downloaded from the cluster it is by default stored in a sub folder in the current directory. If the data should be downloaded at a different location this can be specified in the output.path. 
#'@param antigen.fasta.path It can be of interest to predict the antibody structure together with the antigen to see interaction. For this purpose a path to a FASTA file containing the amino acid sequence of the antigen can be specified in the antigen.fasta.path argument. This will add the antigen sequence to every antibody prediction.
#'@param fasta.directory.path The prediction function can also be used to predict structure directly from amino acid FASTA files without specifying a the VDJ.mixcr.out argument. For this the path to a directory, congaing all the FASTA files of interest can be specified in the fasta.directory.path argument. The files just need to have the .fasta extension. If multiple FASTA files are in the directory, the function will predict all of them separately. 
#'
#'
#'
#'@return This function returns a list with the VDJ.mixcr.out in the first element and a list of pdb files as a second element
#'

AlphaFold_prediction <- function(VDJ.mixcr.out,
                                 cells.to.predict,
                                 max.template.date,
                                 dir.name,
                                 fasta.storage.path,
                                 euler.user.name,
                                 rm.local.fasta,
                                 import,
                                 import.local.path,
                                 import.local.dirnames,
                                 euler.dirname,
                                 euler.dirpath,
                                 n.ranked,
                                 rm.euler.files,
                                 rm.local.output,
                                 output.path,
                                 antigen.fasta.path,
                                 fasta.directory.path,
                                 platypus.version
                                 
    ){
  
  
  
  
  if(missing(fasta.storage.path)) {CurDir <- getwd()}
  else {CurDir <- fasta.storage.path}
  
  if(missing(euler.user.name)) {euler.user.name <- FALSE}
  if(missing(import)) {import <- FALSE}
  if(missing(platypus.version)) {platypus.version <- "v3"}
  else if(platypus.version != "v3") {warning("This function was developed for Platypus version 3 only.")}
  
  if(missing(rm.local.fasta)) {rm.local.fasta <- TRUE }
  if(missing(rm.euler.files)) {rm.euler.files <- FALSE }
  if(missing(rm.local.output)) {rm.local.output <- TRUE }
  if(missing(antigen.fasta.path)) {antigen.fasta.path <- FALSE}
  if(missing(fasta.directory.path)) {fasta.directory.path <- FALSE}
  
  
  
  antigen.fasta <- NULL
  
  
  
  if(import == FALSE){
    if(fasta.directory.path == FALSE){
  
      if(missing(dir.name)){dir.name <- "AlphaFold_Fasta"}
      if(missing(max.template.date)) {max.template.date <- NULL}
      if(missing(cells.to.predict)) {stop("Missing argument for cells.to.predict : Provide a list with barcodes or set it to 'ALL' to predict the structure of all cells")}
      if(cells.to.predict[1] != "ALL") {VDJ.mixcr.out <- dplyr::filter(VDJ.mixcr.out, barcode %in% cells.to.predict)}
      
      
      #Filter for cells that have exact 1 VDJ and 1 or 2 VJ sequences
      VDJ_predict <- VDJ.mixcr.out %>% dplyr::filter(Nr_of_VDJ_chains == 1 & between(Nr_of_VJ_chains, 1, 2))
      
      
      #For Cells that have two VJ sequences only select the first one
      VDJ_predict <- VDJ_predict %>% dplyr::mutate(VJ_aa_mixcr = ifelse(grepl(";",VJ_aa_mixcr),
                                                                 paste0(
                                                                   strsplit(VDJ_predict$VJ_aaSeqFR1[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqCDR1[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqFR2[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqCDR2[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqFR3[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqCDR3[1],";")[[1]][1],
                                                                   strsplit(VDJ_predict$VJ_aaSeqFR4[1],";")[[1]][1]
                                                                 ),
                                                                 VJ_aa_mixcr))
      
      #Remove the "_" from the end of the aa sequence
      VDJ_predict <- tidyr::separate(VDJ_predict, VDJ_aa_mixcr,c("VDJ_aa_mixcr", NA) ,sep = "_$",)
      VDJ_predict <- tidyr::separate(VDJ_predict, VJ_aa_mixcr,c("VJ_aa_mixcr", NA) ,sep = "_$",)
      
      ## Filter all amino acid sequences that contain a gap "*" or "_"
      VDJ_predict <- dplyr::filter(VDJ_predict, str_detect(VDJ_aa_mixcr, "[*_]",negate = TRUE))
      
      #Reduce the dataframe to the important columns for structure prediction
      VDJ_predict <- dplyr::select(VDJ_predict,c("barcode","VDJ_aa_mixcr","VJ_aa_mixcr"))
      
      
      
      
      
      
      #Create a directory for the fasta files
      OutDir <- dir.name
      
      #When the directory is already existing a new directory with a number increment at the end will be created
      if (dir.exists(file.path(CurDir, OutDir))){
        
        BaseDir <- OutDir
        
        i = 1
        while (dir.exists(file.path(CurDir, OutDir))){
              OutDir <- paste0(BaseDir,i)
              i <- i + 1
        }
         
        dir.create(file.path(CurDir,OutDir))
        
      }
      else {
        dir.create(file.path(CurDir,OutDir))
      }
      
     
      #Add an antigen structure to the prediction
      if(antigen.fasta.path != FALSE){
        antigen.fasta <- readr::read_file(antigen.fasta.path)
      }
      else{
        antigen.fasta <- ""
      }
      
      #Write the Fasta files for every cell and store it in the directory
      for(i in 1:nrow(VDJ_predict)) {
        cell <- VDJ_predict[i, ]
        
        file_name <- paste0(cell$barcode,".fasta")
        file_name_VDJ <- paste0(cell$barcode,"_","VDJ")
        file_name_VJ <- paste0(cell$barcode,"_","VJ")
        
        
        Fastafile <- file(file.path(CurDir,OutDir,file_name))
        writeLines(c(paste0(">",file_name_VDJ),cell$VDJ_aa_mixcr,'',paste0(">",file_name_VJ),cell$VJ_aa_mixcr,"",antigen.fasta), Fastafile)
        close(Fastafile)
        
        #Create a txt file that contains all the barcodes
        write(cell$barcode,file = file.path(CurDir,OutDir,"barcodes.txt"), append = TRUE)
        
        
      }
      
      
      if(!is.null(max.template.date)) {
        write(max.template.date,file = file.path(CurDir,OutDir,"max.template.date.txt"))
      }
      
      
      
    
    
    
    
      ## Euler
    
      if(euler.user.name != FALSE) {
      
        
        ##Dowload setup and run script to the Fasta directory
        #download.file("https://gitlab.com/lucas.stalder1/ab-structure-prediction/-/raw/main/setup_alphafold_platypus/setup_alphafold_platypus.sh", file.path(CurDir,OutDir,"setup_alphafold_platypus.sh"))
        #download.file("https://gitlab.com/lucas.stalder1/ab-structure-prediction/-/raw/main/run_alphafold/run_alphafold.sh", file.path(CurDir,OutDir,"run_alphafold.sh"))
        utils::download.file("https://storage.googleapis.com/platypus_private/stalder2022a_setup_alphafold_platypus.sh", file.path(CurDir,OutDir,"setup_alphafold_platypus.sh"))
        utils::download.file("https://storage.googleapis.com/platypus_private/stalder2022a_run_alphafold.sh", file.path(CurDir,OutDir,"run_alphafold.sh"))
        
        
        
        ## Establish connection to server; This will ask for a password
        session <- ssh::ssh_connect(paste0(euler.user.name,"@euler.ethz.ch"))
        file_path <- file.path(CurDir, OutDir)
        
        print(session)
        
        # Copy the files to euler on scratch
        ssh::scp_upload(session, file_path, to = paste0("/cluster/scratch/",euler.user.name))
        
        
        #go to scratch on euler
        ssh::ssh_exec_wait(session, command = c(paste0('cd ',"/cluster/scratch/",euler.user.name,"/",OutDir),
                                          "chmod 755 setup_alphafold_platypus.sh",
                                          paste0("bash ","/cluster/scratch/",euler.user.name,"/",OutDir,"/run_alphafold.sh")))
        
      
        #Dissconnect
        ssh::ssh_disconnect(session)
        
        #Remove local files
        if(rm.local.fasta) {
          unlink(file.path(CurDir,OutDir,recursive = TRUE))
        }
      
      
      }
      
    }
  }
  
  
  
  ### Predict the structure of a protein based on a fasta file
  if(fasta.directory.path != FALSE) {
    if(euler.user.name == FALSE) {stop("In order to run alphafold on Euler the 'euler.user.name' argument must be specified")}
    
    OutDir <- basename(fasta.directory.path)
    
    #Write a barcode.txt file
    Fasta.Files <- list.files(fasta.directory.path, pattern = "*.fasta")
    Fasta.Files <- stringr::str_split(Fasta.Files,".fasta")
    Fasta.Files <- sapply(Fasta.Files,"[[",1)
    write(Fasta.Files,file = file.path(fasta.directory.path,"barcodes.txt"), append = TRUE)
    
    ##Dowload setup and run script to the Fasta directory
    utils::download.file("https://gitlab.com/lucas.stalder1/ab-structure-prediction/-/raw/main/setup_alphafold_platypus/setup_alphafold_platypus.sh", file.path(fasta.directory.path,"setup_alphafold_platypus.sh"))
    utils::download.file("https://gitlab.com/lucas.stalder1/ab-structure-prediction/-/raw/main/run_alphafold/run_alphafold.sh", file.path(fasta.directory.path,"run_alphafold.sh"))
    
    ## Establish connection to server; This will ask for a password
    session <- ssh::ssh_connect(paste0(euler.user.name,"@euler.ethz.ch"))
    file_path <- fasta.directory.path
    
    print(session)
    
    # Copy the files to euler on scratch
    ssh::scp_upload(session, file_path, to = paste0("/cluster/scratch/",euler.user.name))
    
    
    #go to scratch on euler
    ssh::ssh_exec_wait(session, command = c(paste0('cd ',"/cluster/scratch/",euler.user.name,"/",OutDir),
                                            "chmod 755 setup_alphafold_platypus.sh",
                                            paste0("bash ","/cluster/scratch/",euler.user.name,"/",OutDir,"/run_alphafold.sh")))
    
    
    #Dissconnect
    ssh::ssh_disconnect(session)
    
  }
  
  
  
  
  
  
  
  
  
  ### IMPORT the pdb files of the highst ranked outputs
  if(import != FALSE) {
    
    if(missing(output.path)) {CurDir <- getwd()}
    else {CurDir <- output.path}
    
    if(missing(n.ranked)) {n.ranked <- 1}
    
    ## if import set to "euler" the function connects to the server imports the pdb files directly from the Alpha Fold output
    if(import == "euler"){
      
      
      if(euler.user.name == FALSE) {stop("The username for euler must be specified to import the structure predictions")}
      
      
      ## Establish connection to server; This will ask for a password
      session <- ssh::ssh_connect(paste0(euler.user.name,"@euler.ethz.ch"))
      
      #Search for the AlphaFold directory
      if(missing(euler.dirname)) {
        euler.dirname <- "AlphaFold_Fasta"
        out <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("cd ","/cluster/scratch/",euler.user.name),'ls | grep "AlphaFold_Fasta*" | wc -l')))
        if(out[1] > 1) {stop("There are multiple AlphaFold_Fasta* files on euler. Please specify euler.dirname")}
      }
      
      ##Check if the path to the directory is present
      if(missing(euler.dirpath)){
        out <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("[ -d ","/cluster/scratch/",euler.user.name,"/",euler.dirname," ] && echo TRUE"))))
        if(out[1] != "TRUE") {stop(paste("The AlphaFold_Fasta* folder could not be found at:",paste0("/cluster/scratch/",euler.user.name,"/",euler.dirname),"Please sepcify 'euler.dirpath'"))}
        euler.dirpath <- paste0("/cluster/scratch/",euler.user.name,"/",euler.dirname)
      }
      else{
        out <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("[ -d ",euler.dirpath," ] && echo TRUE"))))
        if(out[1] != "TRUE") {stop(paste("The AlphaFold_Fasta* folder could not be found at the indicated path:",euler.dirpath))}
      }
      
      
      ##Check if all the structures are in the output directory
      out1 <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("cd ",euler.dirpath),'ls | grep " *.bsub" | wc -l')))
      out2 <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("cd ",euler.dirpath,"/output"),'ls | wc -l')))
      
      if(out1[1] > out2[1]) {
        if(askYesNo("There were more jobs sumitted than present in the output directory. The prediction might not yet be finsihed. Do you wanna proceed") != TRUE) {stop("The import was stopped")}
      }
      
      
      ##Download the output files to Output_Alphafold Directory
      
      #When the directory is already existing a new directory with a number increment at the end will be created
      OutDir <- "Output_AlphaFold"
      if (dir.exists(file.path(CurDir, OutDir))){
        BaseDir <- OutDir
        
        i = 1
        while (dir.exists(file.path(CurDir, OutDir))){
          OutDir <- paste0(BaseDir,i)
          i <- i + 1
        }
        
        dir.create(file.path(CurDir,OutDir))
        
      }
      else {
        dir.create(file.path(CurDir,OutDir))
      }
      
      ##Download all the highest ranked outputs to the Output folder
      pdb_list <- list()
      out <- utils::capture.output(ssh::ssh_exec_wait(session, command = c(paste0("ls ",euler.dirpath,"/output"))))
      for(i in 2:length(out)-1 ){
        
        
        
        ssh::ssh_exec_wait(session, command = c(paste0("cd ",euler.dirpath,"/output/",out[i],"/",out[i],"/"),paste0("mkdir ",out[i],"_ranked"),paste0("mv ranked* ",out[i],"_ranked")))
        ssh::scp_download(session,paste0(euler.dirpath,"/output/",out[i],"/",out[i],"/",out[i],"_ranked/"), to = paste0(CurDir,"/",OutDir))
        
        # Read in the pdb files to a list object with "n.ranked" number of highest ranked structure predictions
        n_list <- list()
        n_name <- c()
        for(n in 1:n.ranked) {
          n_list[[length(n_list)+1]] <- bio3d::read.pdb(paste0(CurDir,"/",OutDir,"/",out[i],"_ranked/","ranked_",n-1,".pdb"))
          n_name <- c(n_name,paste0("ranked_",n-1,".pdb"))
        }
        names(n_list) <- n_name
        pdb_list[[length(pdb_list)+1]] <- n_list
        
      }
      names(pdb_list) <- out[2:length(out)-1]
      
      #Delete the Euler Directory
      if(rm.euler.files) {
        ssh:ssh_exec_wait(session,command = c("rm -r ",euler.dir.path))
      }
      
      #Dissconnect
      ssh::ssh_disconnect(session)
      
      #Delete the local output directory
      if(rm.local.output) {
        unlink(file.path(CurDir,OutDir),recursive = TRUE)
      }
      
    } # End of import Euler
    
    
    
    ## Import the pdb files form an already downloaded euler output with the right structure
    if(import == "local") {
      
      if(missing(import.local.path)) {stop("For loacal import of pdb files the 'import.local.path' must be specified")}
      if(missing(import.local.dirnames)) {
        import.local.dirnames <- list.dirs(import.local.path, full.names = FALSE, recursive = FALSE)
      }
      
      
      pdb_list <- list()
      for(i in 1:length(import.local.dirnames)) {
      
        # Read in the pdb files to a list object with "n.ranked" number of highest ranked structure predictions
        n_list <- list()
        n_name <- c()
        for(n in 1:n.ranked) {
          n_list[[length(n_list)+1]] <- bio3d::read.pdb(paste0(import.local.path,"/",import.local.dirnames[i],"/","ranked_",n-1,".pdb"))
          n_name <- c(n_name,paste0("ranked_",n-1,".pdb"))
        }
        names(n_list) <- n_name
        pdb_list[[length(pdb_list)+1]] <- n_list
        
      }
    
      import.local.barcodes <- c()
      for(DirName in import.local.dirnames){
        import.local.barcodes <- c(import.local.barcodes,strsplit(DirName,"_ranked")[[1]])
      }
      names(pdb_list) <- import.local.barcodes
      
      
    } # End local imput
    

    if(missing(VDJ.mixcr.out)){
      Function_Output <- list("PDB_AF_structure" = pdb_list)
    }
    
    else{
    ### Add Null lines to the list for missing / filtered barcodes
    missing_barcode <- dplyr::filter(VDJ.mixcr.out, !barcode %in% names(pdb_list))$barcode
    empty_list <- vector(mode = "list", length = length(missing_barcode))
    names(empty_list) <- missing_barcode
    pdb_list <- append(pdb_list,empty_list)
    
    ##Create an output list object with the VDJ.mixcr.out in position 1 and the structure list object in position 2
    Function_Output <- list("VDJ.mixcr.out" = VDJ.mixcr.out,
         "PDB_AF_structure" = pdb_list)
    }

    return(Function_Output)

    
  
  } # End of import
} #End of function  





