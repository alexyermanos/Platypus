#' Function to set up the data so that it can be read and processed by CellPhoneDB, which saves the results of the analysis in a directory "out" and adds them in a new vgm (output of VDJ_GEX_matrix function) list item (CellPhoneDB). Needs Python to be installed in the computer. Running time can take some minutes. 
#' @param vgm.input Output of the VDJ_GEX_matrix function 
#' @param column Character vector. Mandatory. Column name of VDJ_GEX_matrix[[2]] where the groups to be tested for interactions are located
#' @param groups Strings vector. Mandatory. Vector with groups of the indicated column to be compared.
#' @param organism Character vector. Defaults to "human". If == "mouse" the function converts the gene's mouse names into the human ones. 
#' @param gene.id Character vector. Defaults to "ensembl". Possible arguments: "ensembl" , "hgnc_symbol", "gene_name". Indicates the gene ID used by CellPhoneDB during the analysis. CellPhoneDB specific method argument.
#' @param  analysis.method Character vector. Defaults to "statistical_analysis" Possible arguments: "statistical_analysis". CellPhoneDB is developing also "degs_analysis" method, which will be included among the possible arguments once released. Indicates the analysis method used by CellPhoneDB. CellPhoneDB specific method argument.
#' @param project.name Character vector. Defaults to NULL. If a name is given by the user, a subfolder with this name is created in the output folder. CellPhoneDB specific method argument.
#' @param iterations Numerical. Defaults to 1000. Number of iterations for the statistical analysis. CellPhoneDB specific method argument.
#' @param threshold Numerical. By defaults not specified. % of cells expressing the specific ligand/receptor. Range: 0<= threshold <=1. 
#' @param result.precision Numerical. Defaults to 3. Number of decimal digits in results. CellPhoneDB specific method argument.
#' @param subsampling Logical. Defaults to FALSE. If set to TRUE it enables subsampling. CellPhoneDB specific method argument.
#' @param subsampling.num.pc Numerical. Defaults to 100, if subsampling == TRUE. Number of PCs to use. CellPhoneDB specific method argument.
#' @param subsampling.num.cells Numerical. Defaults to 1/3 of cells, if subsampling == TRUE. Number of cells to subsample. CellPhoneDB specific method argument.
#' @param subsampling.log Logical. No default, mandatory when subsampling. Enables subsampling log1p for non log-transformed data inputs. CellPhoneDB specific method argument.
#' @param pvalue Numerical. Defaults to 0.05. P-value threshold. CellPhoneDB specific statistical method argument.
#' @param debug.seed Numerical. Deafults to -1. Debug random seed. To disable it use a value >=0. CellPhoneDB specific statistical method argument
#' @param threads Numerical. Defaults to -1. Number of threads to use (needs to be >=1). CellPhoneDB specific statistical method argument. 
#' @param install.cellphonedb Logical. Defaults to TRUE. Installs the CellPhoneDB Python package if set ==TRUE. 
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter
#' @return VDJ_GEX_matrix object with additional list item (VDJ_GEX_matrix[[10]]), containing results and plots of CellPhoneDB analysis. Saves in the directory the input files, results and plots of CellPhoneDB analysis. 
#' @export
#' @examples
#' \dontrun{
#' vgm_cellphonedb<-CellPhoneDB_analyse(vgm.input=vgm_m, organism="mouse", groups=c(3,6,9), gene.id="ensembl", analysis.method= "statistical_analysis", install.cellphonedb = FALSE, subsampling= TRUE, subsampling.num.pc=100, subsampling.num.cells=70, subsampling.log=TRUE, project.name = "test")
#'}


#Depending on the state of the connection to ensembl and whether this is down, the function might not work if it needs to convert the genes identity. In these cases, try at some other moment and the connection should hopefully be back. 
#Needs Python to be installed in the computer. 

CellPhoneDB_analyse<-function(vgm.input, #mandatory  #the function takes vgm as an input
                              column, #mandatory #column name of vgm[[2]] where the groups to be tested for interactions are located
                              groups, #mandatory  # vector with groups of interest in the indicated column
                              organism, #default "human", if == "mouse" the function needs to convert the gene's name 
                              ##the next argumemts are the CellPhoneDB specific method arguments 
                              gene.id, #possible arguments: "ensembl" , "hgnc_symbol", "gene_name". By deafult "ensembl"
                              analysis.method,  #possible arguments: "statistical_analysis", in the future also "degs_analysis". By deafult "statistical_analysis"
                              project.name, #Name of the project. A subfolder with this name is created in the output folder
                              iterations, #Number of iterations for the statistical analysis. By default 1000
                              threshold, #% of cells expressing the specific ligand/receptor. By default not specified
                              result.precision, #Number of decimal digits in results. By default 3
                              
                              subsampling, # If set to TRUE it enables subsampling. By default == FALSE 
                              subsampling.num.pc, #Subsampling NumPC argument (number of PCs to use). By default 100, if subsampling == TRUE
                              subsampling.num.cells, #Number of cells to subsample to. By default 1/3 of cells, if subsampling == TRUE
                              subsampling.log, #Enable subsampling log1p for non log-transformed data inputs. No default. Takes boolean values  
                              
                              ## statistical method arguments 
                              pvalue, #P-value threshold. By default 0.05.
                              debug.seed, #Debug random seed -1. To disable it please use a value >=0. By deafult -1
                              threads, #Number of threads to use. >=1. By default -1
                              install.cellphonedb, #Boolean argument that installs the CellPhoneDB Python package if set ==TRUE. By default TRUE. 
                              platypus.version
){
  
  platypus.version <- "v3"
  
  #set default value for column
  if(missing(column)){
    column="seurat_clusters"
  }
  
  #set default value for organism
  if(missing(organism)){
    organism = "human"
  }
  
  #set default value for gene.id
  if(missing(gene.id)){
    gene.id = "ensembl"
  }
  
  #set default value for analysis.method
  if(missing(analysis.method)){
    analysis.method = "statistical_analysis"
  }
  
  #set default value for project.name
  if(missing(project.name)){
    project.name = NULL
  }
  
  #set default value for iterations
  if(missing(iterations)){
    iterations = 1000
  }
  
  #set default value for threshold
  if(missing(threshold)){
    threshold = NULL
  }
  
  #set default value for result.precision
  if(missing(result.precision)){
    result.precision = 3
  }
  
  #set default value for subsampling
  if(missing(subsampling)){
    subsampling = FALSE
  }
  
  #set default value for subsampling specifications
  if(subsampling==TRUE){
    if(missing(subsampling.num.pc)){
      subsampling.num.pc = 100
    }
    
    if(missing(subsampling.num.cells)){
      subsampling.num.cells = "temp"
    }
    
    if(missing(subsampling.log)){
      print("subsampling.log value needs to be specified")
    }
    
  }else{
    
    if(missing(subsampling.num.pc)){
      subsampling.num.pc = NULL
    }
    
    if(missing(subsampling.num.cells)){
      subsampling.num.cells = NULL
    }
    
    if(missing(subsampling.log)){
      subsampling.log = NULL
    }
    
  }
  
  #set default value for pvalue
  if(missing(pvalue)){
    pvalue = 0.05
  }
  
  #set default value for debug.seed
  if(missing(debug.seed)){
    debug.seed = -1
  }
  
  #set default value for threads
  if(missing(threads)){
    threads = -1
  }
  
  #set default value for install.cellphonedb
  if(missing(install.cellphonedb)){
    install.cellphonedb=TRUE
  }
  
  
  ##Format vgm[[2]] data according to CellPhoneDB input
  
  #take subset of vgm[[2]] according to the input column and input groups given by the User 
  
  where_groups<-c()
  for(i in groups){
    #get the index of the groups 
    where_groups<-append(where_groups, which(vgm.input[[2]][[column]]==as.character(i)))
  }
  input.counts<-as.matrix(vgm.input[[2]]@assays$RNA@data[,sort(where_groups)])
  
  #COUNTS file
  
  if(organism == "mouse"){ 
    
    #need to translate gene names
    
    # Basic function to convert mouse to human gene names
    convertMouseGene <- function(x){
      
      #usual code, currently not working because of new release of Ensembl and some bugs in getLDS function related to that 
      #human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      #mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
      
      #use different host mirrors
      human = biomaRt::useMart(biomart ="ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
      mouse = biomaRt::useMart(biomart ="ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
      
      #convert the gene names to the input gene.id
      if(gene.id == "ensembl"){
        hgenes = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=TRUE)
        return(unique(hgenes))
      }else if(gene.id =="hgnc_symbol"){
        hgenes = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_id"), martL = human, uniqueRows=TRUE)
        return(unique(hgenes))
      }else if(gene.id == "gene_name"){
        hgenes = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("external_gene_name"), martL = human, uniqueRows=TRUE)
        return(unique(hgenes))
      }
    }
    
    Mouse_to_human<-convertMouseGene(dimnames(input.counts)[[1]])
    
    #if there are genes which were not converted, don't include them in input.counts file, otherwise substitute the mouse gene names with the human converted ones  
    if (length(which(!(toupper(Mouse_to_human[,1]) %in% dimnames(input.counts)[[1]] )))>0 ){
      dimnames(input.counts)[[1]][which(!(toupper(Mouse_to_human[,1]) %in% dimnames(input.counts)[[1]]))] <-NA
    }
    
    dimnames(input.counts)[[1]][match(toupper(Mouse_to_human[,1]), dimnames(input.counts)[[1]])]<-Mouse_to_human[,2]
    
  }
  
  #save matrix in directory as a .txt file
  
  #write files locally changing directory in console 
  
  # get the current DIRECTORY
  
  my_directory<-getwd()
  dir_counts <-paste(my_directory, "/input.counts.txt", sep="")
  
  write.table(input.counts, dir_counts, quote=FALSE, sep='\t')
  
  #META file
  
  #take subset of vgm[[2]] according to the input column and input groups given by the User
  input.meta <- cbind( 
    #first column of meta file, cells barcodes
    rownames(vgm.input[[2]]@meta.data[ sort( 
      where_groups),] ), 
    #second column of meta file, cells annotations
    as.character(vgm.input[[2]]@meta.data[[column]][sort( 
      where_groups)]))%>%
    as.matrix()
  
  #rewrite the cell annotation including the interacting group that cell belongs to
  n=1
  for(i in groups){
    input.meta[which(input.meta[,2] == as.character(i)), 2]<-paste(c(paste(c("group", n), sep=""), as.character(i)), collapse ="_")
    n<-n+1
  }
  
  
  #save matrix in directory as a .txt file in the current directory 
  
  dir_meta <-paste(my_directory, "/input.meta.txt", sep="")
  write.table(input.meta, dir_meta, quote=FALSE, sep='\t', row.names=FALSE, col.names = c("cell", "cell type"))
  
  
  #set default value for subsampling specification number of cells (that depends on the total numver of cells)
  if(!(is.null(subsampling.num.cells))){
    if(subsampling.num.cells == "temp"){
      subsampling.num.cells = ceiling((1/3)*length(input.meta[,1]))
    }
  }
  
  ##Terminal commands are called from R in order to run CellPhoneDB analysis and plot dotplot and heatmaps.
  
  #The result files and plots are saved in the save directory of the function.
  
  #Avoid installing CellPhoneDB Python package if install.cellphonedb is set to FALSE
  if(install.cellphonedb == TRUE){
    term1<-rstudioapi::terminalExecute("pip3 install cellphonedb")
    #wait until this command has been executed, then close the terminal 
    while (is.null(rstudioapi::terminalExitCode(term1))) {
      Sys.sleep(0.1)
    }
    rstudioapi::terminalKill(term1)
  }
  
  #commands to make sure everything runs smoothly and without errors 
  term2<-rstudioapi::terminalExecute("export LC_ALL=en_US.UTF-8")
  while (is.null(rstudioapi::terminalExitCode(term2))) {
    Sys.sleep(0.1)
  }
  rstudioapi::terminalKill(term2)
  
  term3<-rstudioapi::terminalExecute("export LANG=en_US.UTF-8")
  while (is.null(rstudioapi::terminalExitCode(term3))) {
    Sys.sleep(0.1)
  }
  rstudioapi::terminalKill(term3)
  
  #pasting togther the user-dependent arguments of cellphonedb 
  cellphoneDB_analysis<-paste("cellphonedb method", analysis.method, "input.meta.txt input.counts.txt --counts-data", gene.id, "--iterations", iterations, "--result-precision", result.precision)
  
  #include only the arguments that were specified by the user 
  
  if(!(is.null(project.name))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, " --project-name", project.name)
  }
  
  if(!(is.null(threshold))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--threshold", threshold )
  }
  
  if(subsampling == TRUE){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--subsampling" )
  }
  
  if(!(is.null(subsampling.num.pc))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--subsampling-num-pc", subsampling.num.pc )
  }
  
  if(!(is.null(subsampling.num.cells))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--subsampling-num-cells", subsampling.num.cells )
  }
  
  if(!(is.null(subsampling.log))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, " --subsampling-log", subsampling.log)
  }
  
  if(!(is.null(pvalue))){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--pvalue", pvalue)
  }
  
  if(debug.seed >=0){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--debug-seed", debug.seed)
  }
  
  if(threads >=1 ){
    cellphoneDB_analysis<-paste(cellphoneDB_analysis, "--threads",threads)
  }
  
  term4<-rstudioapi::terminalExecute(cellphoneDB_analysis)
  
  while (is.null(rstudioapi::terminalExitCode(term4))) {
    Sys.sleep(0.1)
  }
  #rstudioapi::terminalKill(term4)
  
  #this makes sure here are no errors due to markupsafe not being available on the last version of Python
  term5<-rstudioapi::terminalExecute("sudo pip3 install markupsafe==2.0.1")
  while (is.null(rstudioapi::terminalExitCode(term5))) {
    Sys.sleep(0.1)
  }
  rstudioapi::terminalKill(term5)
  
  ##plotting the results 
  
  #Specify the files directories
  #Make sure to indicate the right path where to find the file (according to whether an additional project is added in the out folder)
  
  if(!is.null(project.name)){
    deconvoluted_path<-paste(my_directory, "/out/",project.name,"/deconvoluted.txt", sep ="")
    means_path<-paste(my_directory, "/out/",project.name,"/means.txt", sep ="") 
    significant_means_path<-paste(my_directory, "/out/",project.name,"/significant_means.txt", sep ="")
    pvalues_path<-paste(my_directory, "/out/",project.name,"/pvalues.txt", sep ="")
    interaction_count_path<-paste(my_directory, "/out/",project.name,"/interaction_count.txt", sep ="") 
    count_network_path<-paste(my_directory, "/out/",project.name,"/count_network.txt", sep ="")
    
  } else {
    deconvoluted_path<-paste(my_directory, "/out/deconvoluted.txt", sep ="")
    means_path<-paste(my_directory, "/out/means.txt", sep ="")
    significant_means_path<-paste(my_directory, "/out/significant_means.txt", sep ="")
    pvalues_path<-paste(my_directory, "/out/pvalues.txt", sep ="")
    interaction_count_path<-paste(my_directory, "/out/interaction_count.txt", sep ="")
    count_network_path<-paste(my_directory, "/out/count_network.txt", sep ="")
  }
  
  dot_plot <- "cellphonedb plot dot_plot"
  
  #in case a subfolder for the results is specified, add that path to the arguments of the method
  
  if(!(is.null(project.name))){
    dot_plot<-paste(dot_plot, "--output-path", my_directory)
    dot_plot<-paste(dot_plot, "/out/", project.name, sep ="")
    dot_plot<-paste(dot_plot, "--means-path", means_path)
    dot_plot<-paste(dot_plot, "--pvalues-path", pvalues_path)
  }else{
    dot_plot<-paste(dot_plot, "--output-path", my_directory)
    dot_plot<-paste(dot_plot, "/out", sep ="")
    dot_plot<-paste(dot_plot, "--means-path", means_path)
    dot_plot<-paste(dot_plot, "--pvalues-path", pvalues_path)
  }
  
  term6<-rstudioapi::terminalExecute(dot_plot)
  while (is.null(rstudioapi::terminalExitCode(term6))) {
    Sys.sleep(0.1)
  }
  rstudioapi::terminalKill(term6)
  
  heatmap_plot<- paste("cellphonedb plot heatmap_plot", my_directory)
  heatmap_plot<- paste(heatmap_plot, "/input.meta.txt", sep = "")
  
  if(!(is.null(project.name))){
    heatmap_plot<-paste(heatmap_plot, "--output-path", my_directory)
    heatmap_plot<-paste(heatmap_plot, "/out/", project.name, sep ="")
    heatmap_plot<-paste(heatmap_plot, "--pvalues-path", pvalues_path)
  }else{
    heatmap_plot<-paste(heatmap_plot, "--output-path", my_directory)
    heatmap_plot<-paste(heatmap_plot, "/out", sep ="")
    heatmap_plot<-paste(heatmap_plot, "--pvalues-path", pvalues_path)
  }
  
  term7<-rstudioapi::terminalExecute(heatmap_plot)
  while (is.null(rstudioapi::terminalExitCode(term7))) {
    Sys.sleep(0.1)
  }
  rstudioapi::terminalKill(term7)
  
  ## Read txt data into R as data frames and combine them in a list to be included in the vgm 
  
  #Write file directories of images. Make sure to indicate the right path where to find the file (according to whether an additional project is added in the out folder) 
  
  if(!is.null(project.name)){
    heatmap_log_path<-paste(my_directory, "/out/",project.name,"/heatmap_log_count.pdf", sep ="")
    heatmap_path<-paste(my_directory, "/out/",project.name,"/heatmap_count.pdf", sep ="") 
    plot_path<-paste(my_directory, "/out/",project.name,"/plot.pdf", sep ="")
  } else {
    heatmap_log_path<-paste(my_directory, "/out/heatmap_log_count.pdf", sep ="")
    heatmap_path<-paste(my_directory, "/out/heatmap_count.pdf", sep ="")
    plot_path<-paste(my_directory, "/out/plot.pdf", sep ="")
  }
  
  #Create the list reading the .txt files as data.frames and the .pdf files as images 
  cellphonedb.list<-list(
    deconvoluted = read.table(deconvoluted_path, header = TRUE, sep = "\t", dec = "."),
    means=read.table(means_path, header = TRUE, sep = "\t", dec = "."),
    significant.means=read.table(significant_means_path, header = TRUE, sep = "\t", dec = "."),
    pvalues=read.table(pvalues_path, header = TRUE, sep = "\t", dec = "."),
    interaction.count=read.table(interaction_count_path, header = TRUE, sep = "\t", dec = "."),
    count.network=read.table(count_network_path, header = TRUE, sep = "\t", dec = "."),
    heatmap_log_count=knitr::include_graphics(heatmap_log_path),
    heatmap_count=knitr::include_graphics(heatmap_path),
    dot_plot=knitr::include_graphics(plot_path)
  )
  
  #integrate the list in vgm as the 10th element of vgm, keeping the current list elements
  
  vgm.input[[10]] <- cellphonedb.list
  names(vgm.input)[10]<- "CellPhoneDB"
  return(vgm.input)
  
}


#' Function to cutomise the Dot Plot of CellPhoneDB analysis results. 
#' @param vgm.input Output of the VDJ_GEX_matrix function. Mandatory. Object where to save the dotplot.
#' @param selected.rows Strings vector. Defaults to NULL. Vector of rows to plot (interacting genes pair), one per line.
#' @param selected.columns Strings vector. Defaults to  NULL. Vector of columns (interacting groups) to plot, one per line
#' @param threshold.type Character vector. Defaults to NULL. Possible arguments: "pvalue", "log2means", "pvalue_topn", "log2means_topn". Which thresholding system the user wants to use. 
#' @param threshold.value Numerical. Defaults to NULL. Value below/above (depending on whether it's pvalue or log2means) which genes to plot are selected. 
#' @param project.name Character vector. Defaults to NULL. Subfolder where to find the output of the CellPhoneDB analysis and where to save the dot_plot output plot.
#' @param filename Character vector. Defaults to "selection_plot.pdf". Name of the file where the dot_plot output plot will be saved. 
#' @param width Numerical. Defaults to 8. Width of the plot. 
#' @param height Numerical. Defaults to 10. Height of the plot. 
#' @param text.size Numerical. Defaults to 12. Font size of the plot axes
#' @param return.vector Logical. Defaults to FALSE. If set to TRUE, it includes the vector of genes_pairs present in the dot_plot in the VDJ_GEX_matrix[[10]] list. 
#' @param platypus.version This function works with "v3" only, there is no need to set this parameter
#' @return VDJ_GEX_matrix object with output dot plot added to VDJ_GEX_matrix[[10]] list. 
#' @export
#' @examples
#' \dontrun{
#' vgm_cellphonedb<-dot_plot(vgm.input=vgm_cellphonedb, selected.columns = c("group_1_3|group_2_6", "group_1_3|group_3_9", "group_2_6|group_1_3", "group_2_6|group_3_9", "group_3_9|group_1_3", "group_3_9|group_2_6"), threshold.type="pvalue_topn", threshold.value=50, project.name = "test", height = 12, width=6, text.size=10, return.vector=TRUE)
#'}

dot_plot = function(vgm.input, #vgm where to save the dotplot
                    
                    selected.rows, #Defaults to NULL #vector of rows to plot (interacting genes pair), one per line 
                    selected.columns , #Defaults to  NULL #vector of columns (interacting groups) to plot, one per line
                    threshold.type, #Defaults to NULL. Possible arguments: "pvalue", "log2means", "pvalue_topn", "log2means_topn". Which thresholding system the user wants to use. 
                    threshold.value, #Defaults to NULL #value below/above (depending on whether it's pvalue or log2means) which select the genes to plot
                    project.name , #Defaults to  NULL
                    filename, #Defaults to "selection_plot.pdf"
                    width, #Defaults to 8
                    height, #Defaults to 10
                    text.size, #Defaults to 12
                    return.vector, #Defaults to FALSE. If set to TRUE, it includes the vector of genes_pairs present in the dot_plot in the vgm CellPhoneDB list item. 
                    platypus.version #Defaults to v3"
                    
){
  
  platypus.version <- "v3"
  
  #set default value for selected.rows
  if(missing(selected.rows)){
    selected.rows=NULL
  }
  
  #set default value for selected.columns
  if(missing(selected.columns)){
    selected.columns = NULL
  }
  
  #set default value for pvalue.threshold
  if(missing(threshold.type)){
    threshold.type = NULL
  }
  
  #set default value for pvalue.threshold
  if(missing(threshold.value)){
    threshold.value = NULL
  }
  
  #set default value for filename
  if(missing(filename)){
    filename = "selection_plot.pdf"
  }
  
  #set default value for project.name
  if(missing(project.name)){
    project.name = NULL
  }
  
  #set default value for width
  if(missing(width)){
    width = 8
  }
  
  #set default value for height
  if(missing(height)){
    height = 10
  }
  
  #set default value for text.size
  if(missing(text.size)){
    text.size = 12
  }
  
  #set default value for return.vector
  if(missing(return.vector)){
    return.vector = FALSE
  }
  
  means_separator = '\t'
  pvalues_separator = '\t'
  
  #get files directories
  my_directory<-getwd()
  
  if(!is.null(project.name)){
    means_path<-paste(my_directory, "/out/",project.name,"/means.txt", sep ="") 
    pvalues_path<-paste(my_directory, "/out/",project.name,"/pvalues.txt", sep ="")
    plot_path<-paste(my_directory, "/out/",project.name, sep ="")
  } else {
    means_path<-paste(my_directory, "/out/means.txt", sep ="")
    pvalues_path<-paste(my_directory, "/out/pvalues.txt", sep ="")
    plot_path<-paste(my_directory, "/out", sep ="")
  }
  
  all_pval = read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)
  all_means = read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected.rows)){
    selected.rows = intr_pairs
  }
  
  if(is.null(selected.columns)){
    selected.columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected.rows, intr_pairs), selected.columns]
  sel_means = all_means[match(selected.rows, intr_pairs), selected.columns]
  intr_pairs<-selected.rows
  #select only rows where at least one interacting pair has p-value velow the threshold 
  if(!(is.null(threshold.value))&!(is.null(threshold.type))){
    
    #according to threshold.type values are selected in different ways 
    if(threshold.type == "pvalue"){
      sel_col<-NULL
      sel_pval$interacting_pairs<-selected.rows
      #create a column counting how many times there's a pvalue below the threshold in the row 
      sel_pval$count<-rep(0, length(selected.rows))
      for(g in 1:length(selected.rows)){
        for(c in 1:length(selected.columns)){
          if(!(is.na(sel_pval[g,c])) & sel_pval[g,c]<=threshold.value){
            sel_pval$count[g]<-sel_pval$count[g]+1
          }
        }
      } 
      #only rows with at least one pvalue below the threshold are selected
      sel_pval<-filter(sel_pval, count>0) 
      sel_pval$count<-NULL
      selected.rows<-NULL
      selected.rows<-sel_pval$interacting_pairs
      sel_pval$interacting_pairs<-NULL
      sel_means = sel_means[match(selected.rows, intr_pairs), ]
    }
    
    #take the same approach to select the means 
    if(threshold.type == "log2means"){
      sel_col<-NULL
      sel_means$interacting_pairs<-selected.rows
      sel_means$count<-rep(0, length(selected.rows))
      for(g in 1:length(selected.rows)){
        for(c in 1:length(selected.columns)){
          if(!(is.na(sel_means[g,c])) & sel_means[g,c]>=threshold.value){
            sel_means$count[g]<-sel_means$count[g]+1
          }
        }
      } 
      sel_means<-filter(sel_means, count>0) 
      sel_means$count<-NULL
      selected.rows<-NULL
      selected.rows<-sel_means$interacting_pairs
      sel_means$interacting_pairs<-NULL
      sel_pval = sel_pval[match(selected.rows, intr_pairs), ]
    }
    
    if(threshold.type == "pvalue_topn"){
      sel_col<-NULL
      sel_pval$interacting_pairs<-selected.rows
      
      #pivot the data.frame so that we have all the pvalues in a single column with the respective interacting groups pair
      sel_piv<-pivot_longer(sel_pval, all_of(selected.columns), names_to = "group_pair", values_to = "temp")
      sel_piv<-cbind(pvalue=sel_piv$temp, sel_piv)
      sel_piv$temp<-NULL
      
      #order data.frame accoridng to the pvalue values
      sel_piv<-sel_piv[do.call(order,sel_piv),]
      
      # take only the first n rows
      sel_piv<-sel_piv[1:threshold.value,]
      sel_pval<-sel_pval[match(unique(sel_piv$interacting_pairs), sel_pval$interacting_pairs),] 
      
      selected.rows<-NULL
      selected.rows<-sel_pval$interacting_pairs
      sel_pval$interacting_pairs<-NULL
      sel_means = sel_means[match(selected.rows, intr_pairs), ]
    }
    
    if(threshold.type == "log2means_topn"){
      sel_col<-NULL
      sel_means$interacting_pairs<-selected.rows
      
      #pivot the data.frame so that we have all the pvalues in a single column with the respective interacting groups pair
      sel_piv<-pivot_longer(sel_means, all_of(selected.columns), names_to = "group_pair", values_to = "temp")
      sel_piv<-cbind(mean=sel_piv$temp, sel_piv)
      sel_piv$temp<-NULL
      
      #order data.frame accoridng to the pvalue values
      sel_piv<-sel_piv[do.call("order",c(sel_piv, list(decreasing=TRUE))),]
      
      # take only the first n rows
      sel_piv<-sel_piv[1:threshold.value,]
      sel_means<-sel_means[match(unique(sel_piv$interacting_pairs), sel_means$interacting_pairs),] 
      
      selected.rows<-NULL
      selected.rows<-sel_means$interacting_pairs
      sel_means$interacting_pairs<-NULL
      sel_pval = sel_pval[match(selected.rows, intr_pairs), ]
    }
  }
  
  df_names = expand.grid(selected.rows, selected.columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- grDevices::colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  #make the dot plot, including only the selcted rows and columns 
  selection_dot_plot<-ggplot2::ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=text.size, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size=text.size, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  #save the plot in the output directory
  ggplot2::ggsave(filename, path=plot_path, width = width, height = height, limitsize=FALSE)
  
  #include the plot in the vgm 
  vgm.input[[10]]$selection_dot_plot<-selection_dot_plot
  
  #if argument return.vector==TRUE, include the vector of interacting gene pairs in the vgm 
  if(return.vector==TRUE){
    vgm.input[[10]]$interacting_pairs<-selected.rows
  }
  return(vgm.input)
}

