#' Function to import the annotations and alignments from IgBLAST output into the VDJ dataframe.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description  Imports the IgBLAST annotations and alignments from IgBLAST output files, stored in the output folders of Cell Ranger, into a VDJ dataframe obtained from the minimal_VDJ() function in Platypus.
#' @param VDJ dataframe - VDJ object as obtained from the minimal_VDJ() function in Platypus.
#' @param VDJ.directory string - path to parent directory containing the output folders (one folder for each sample) of Cell Ranger. This pipeline assumes that the sample IDs and contigs IDs have not been modified and that the IgBLAST output file names have not been changed from the default changeo settings. Each sample directory should contain a 'filtered_contig_igblast_db-pass.tsv' file.
#' @param file.path.list list - list containing the paths to the 'filtered_contig_igblast_db-pass.tsv' files, in which the names of each item should refer to an sample ID.
#' @param method string - denotes the way the IgBLAST germline annotations from the 'filtered_contig_igblast_db-pass.tsv' files should be appended to the VDJ dataframe. Options: 'replace' or 'attach'. Defaults to 'append'.
#' 'replace'  : The original annotation columns in the VDJ dataframe are replaced with the IgBLAST annotations. The original columns are kept with the suffix '_10x'.
#' 'append'   : The IgBLAST annotation columns are stored in columns with the suffix '_IgBLAST'.
#' @return The VDJ dataframe with the appended IgBLAST annotations and alignments.
#' @examples
#' \dontrun{
#' VDJ <- import_igblast_annotations <- function(VDJ,
#'                                               VDJ.directory)
#'}


VDJ_import_igblast_annotations <- function(VDJ,
                                           VDJ.directory,
                                           file.path.list,
                                           method){
  
  # If no 'VDJ' dataframe is provided, a message is returned and execution is stopped
  if(missing(VDJ)){stop("ERROR: Please provide a VDJ dataframe as obtained from the 'minimal_VDJ()' function in Platypus.")}
  # If not all default columns are present in the input 'VDJ' dataframe, a message is returned and execution is stopped
  default_columns <- c("sample_id", "barcode", "celltype", "isotype", "VDJ_contig_id", "VJ_contig_id", "VDJ_consensus_id", "VJ_consensus_id", "clonotype_id", "clonotype_frequency", "VDJ_chain", "VDJ_chain_count", "VDJ_umis", "VDJ_vgene", "VDJ_dgene", "VDJ_jgene", "VDJ_cgene", "VDJ_fwr1_nt", "VDJ_fwr1_aa", "VDJ_cdr1_nt", "VDJ_cdr1_aa", "VDJ_fwr2_nt", "VDJ_fwr2_aa", "VDJ_cdr2_nt", "VDJ_cdr2_aa", "VDJ_fwr3_nt", "VDJ_fwr3_aa", "VDJ_cdr3_nt", "VDJ_cdr3_aa", "VDJ_fwr4_nt", "VDJ_fwr4_aa", "VDJ_sequence_nt_raw", "VDJ_sequence_nt_trimmed", "VDJ_sequence_aa_trimmed", "VDJ_consensus_nt_raw", "VDJ_consensus_nt_trimmed", "VDJ_consensus_aa_trimmed", "VDJ_germline_nt_raw", "VDJ_germline_nt_trimmed", "VDJ_germline_aa_trimmed", "VJ_chain", "VJ_chain_count", "VJ_umis", "VJ_vgene", "VJ_jgene", "VJ_cgene", "VJ_fwr1_nt", "VJ_fwr1_aa", "VJ_cdr1_nt", "VJ_cdr1_aa", "VJ_fwr2_nt", "VJ_fwr2_aa", "VJ_cdr2_nt", "VJ_cdr2_aa", "VJ_fwr3_nt", "VJ_fwr3_aa", "VJ_cdr3_nt", "VJ_cdr3_aa", "VJ_fwr4_nt", "VJ_fwr4_aa", "VJ_sequence_nt_raw", "VJ_sequence_nt_trimmed", "VJ_sequence_aa_trimmed", "VJ_consensus_nt_raw", "VJ_consensus_nt_trimmed", "VJ_consensus_aa_trimmed", "VJ_germline_nt_raw", "VJ_germline_nt_trimmed", "VJ_germline_aa_trimmed")
  if(!all(default_columns %in% colnames(VDJ))){stop(paste0("ERROR: The provided VDJ dataframe lacks the following columns: ", paste(default_columns[!default_columns %in% colnames(VDJ)], collapse = " ")))}
  # If both a 'VDJ.directory' and 'file.path.list' are specified, a message is returned and execution is stopped
  if(!missing(VDJ.directory) & !missing(file.path.list)){stop('ERROR: Both a path to the parent directory containing the output folders (one for each sample) from Cell Ranger and a list of file paths are given as input. Please provide one parent direcotry or one list of file paths.')}
  # Create list with paths to TSV files stored in Cell Ranger output directories (one path for each sample) within the parent directory 'VDJ.directory'
  if(missing(file.path.list)){file.path.list <- list.dirs(path=VDJ.directory,full.names=T,recursive=F); file.path.list <- paste0(file.path.list, "/filtered_contig_igblast_db-pass.tsv"); names(file.path.list) <- basename(dirname(file.path.list))}
  # Check whether all sample IDs in the VDJ dataframe are present in the 'file.path.list'
  if(!all(unique(VDJ$sample_id) %in% names(file.path.list))){stop("ERROR: The VDJ dataframe contains sample IDs that are not found in the VDJ directory.")}
  # Check whether all files in the 'file.path.list' exist
  if(!all(file.exists(file.path.list))){stop(paste0("ERROR: The IgBLAST output files could not be found of sample", names(file.path.list)[!file.exists(file.path.list)], "."))}
  # If the 'method' parameter is missing, it is set to 'append'
  if(missing(method)){method <- "append"}
  # If the 'method' parameter is not recognized, a message is returned and execution is stopped
  if(!method %in% c("replace", "append")){stop("ERROR: The specified method is not recognized. Please choose from the following options: 'replace' or 'append'.")}
  
  # Read in all IgBLAST output files and store in 'igblast_annotations' list, while removing the '-1' from the 'sequence_id' column
  igblast_annotations <- lapply(file.path.list, function(x){
    df <- read.table(x, header = TRUE, sep = "\t")
    df$sequence_id <- gsub(pattern = "-1", replacement = "", df$sequence_id)
    return(df)})
  
  
  # Make selection of columns that contain the gene and region annotations
  annotation_columns <- c("VDJ_vgene", "VDJ_vgene_cigar", "VDJ_vgene_start", "VDJ_vgene_end", "VDJ_vgene_germline_start", "VDJ_vgene_germline_end",
                          "VDJ_dgene", "VDJ_dgene_cigar", "VDJ_dgene_start", "VDJ_dgene_end", "VDJ_dgene_germline_start", "VDJ_dgene_germline_end",
                          "VDJ_jgene", "VDJ_jgene_cigar", "VDJ_jgene_start", "VDJ_jgene_end", "VDJ_jgene_germline_start", "VDJ_jgene_germline_end",
                          "VDJ_cgene", 
                          "VDJ_fwr1_nt", "VDJ_cdr1_nt", "VDJ_fwr2_nt", "VDJ_cdr2_nt", "VDJ_fwr3_nt", "VDJ_cdr3_nt", "VDJ_fwr4_nt", "VDJ_junction_nt",
                          "VDJ_fwr1_aa", "VDJ_cdr1_aa", "VDJ_fwr2_aa", "VDJ_cdr2_aa", "VDJ_fwr3_aa", "VDJ_cdr3_aa", "VDJ_fwr4_aa", "VDJ_junction_aa",
                               
                          "VJ_vgene", "VJ_vgene_cigar", "VJ_vgene_start", "VJ_vgene_end", "VJ_vgene_germline_start", "VJ_vgene_germline_end",
                          "VJ_jgene", "VJ_jgene_cigar", "VJ_jgene_start", "VJ_jgene_end", "VJ_jgene_germline_start", "VJ_jgene_germline_end", 
                          "VJ_cgene", 
                          "VJ_fwr1_nt", "VJ_cdr1_nt", "VJ_fwr2_nt", "VJ_cdr2_nt", "VJ_fwr3_nt", "VJ_cdr3_nt", "VJ_fwr4_nt", "VJ_junction_nt",
                          "VJ_fwr1_aa", "VJ_cdr1_aa", "VJ_fwr2_aa", "VJ_cdr2_aa", "VJ_fwr3_aa", "VJ_cdr3_aa", "VJ_fwr4_aa", "VJ_junction_aa")
  
  # Make selection of columns that contain the alignments
  alignment_columns <- c("VDJ_sequence_alignment_nt", "VDJ_sequence_alignment_aa", "VDJ_germline_alignment_nt", "VDJ_germline_alignment_aa",
                         "VJ_sequence_alignment_nt", "VJ_sequence_alignment_aa", "VJ_germline_alignment_nt", "VJ_germline_alignment_aa")
  
  
  translate_IMGT_sequence <- function(seq){
    
    # Converts DNA sequence in IMGT numbering format into protein sequence in IMGT numbering format
    # Arguments:
    # - seq: DNA sequence in IMGT numbering format (in which gaps / unoccupied locations are indicated by a '.')
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # If the input sequence is not a multiple of 3, the sequence is first trimmed before translation
    if((nchar(seq) %% 3) != 0){seq <- substr(seq, 1, nchar(seq) - (nchar(seq)%%3))}
    
    # Define the codon table
    codon_table <- list(
      "AAA" = "K", "AAC" = "N", "AAG" = "K", "AAT" = "N",
      "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACT" = "T",
      "AGA" = "R", "AGC" = "S", "AGG" = "R", "AGT" = "S",
      "ATA" = "I", "ATC" = "I", "ATG" = "M", "ATT" = "I",
      "CAA" = "Q", "CAC" = "H", "CAG" = "Q", "CAT" = "H",
      "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCT" = "P",
      "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGT" = "R",
      "CTA" = "L", "CTC" = "L", "CTG" = "L", "CTT" = "L",
      "GAA" = "E", "GAC" = "D", "GAG" = "E", "GAT" = "D",
      "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCT" = "A",
      "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGT" = "G",
      "GTA" = "V", "GTC" = "V", "GTG" = "V", "GTT" = "V",
      "TAA" = "*", "TAC" = "Y", "TAG" = "*", "TAT" = "Y",
      "TCA" = "S", "TCC" = "S", "TCG" = "S", "TCT" = "S",
      "TGA" = "*", "TGC" = "C", "TGG" = "W", "TGT" = "C",
      "TTA" = "L", "TTC" = "F", "TTG" = "L", "TTT" = "F",
      "..." = "."
    )
    
    # Extract codons from input 'seq'
    codons <- sapply(1:(nchar(seq)/3), function(x) substr(seq, (x*3-2), (x*3)))
    
    # Translate each codon separately into an amino acid using the defined 'codon_table'
    aa_seq <- paste(codon_table[codons], collapse = "")
    
    # Return protein sequence
    return(aa_seq)
  }
  
  
  retrieve_cell_annotations <- function(igblast_annotations,
                                        cell_info,
                                        translate.IMGT.sequences){
    
    # Retrieves the annotations from the IgBLAST output files for one cell using the sample ID and contig IDs.
    # Arguments:
    # - igblast_annotations: list of dataframes containing the IgBLAST output (one dataframe for each sample)
    # - cell_info: list of 'sample_id', 'VDJ_contig_id', and 'VJ_contig_ID'
    # - translate.IMGT.sequences: bool indicating whether the nucleotide sequences (with IMGT numbering) should be translated into protein sequences
    # Authors: Valentijn Tromp, Daphne van Ginneken
    
    # Extract the 'sample_id', the 'VDJ_contig_id', and the 'VJ_contig_id' from the 'cell_info'
    sample_id <- cell_info[["sample_id"]]
    VDJ_contig_id <- cell_info[["VDJ_contig_id"]]
    VJ_contig_id <- cell_info[["VJ_contig_id"]]
    
    # Make selection of columns to be selected from IgBLAST output file
    columns <- c("sequence_id", "sequence_alignment", "germline_alignment",
                 "v_call", "v_cigar", "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
                 "d_call", "d_cigar", "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
                 "j_call", "j_cigar", "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end",
                 "c_call",
                 "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4",
                 "junction", "junction_aa")
    
    # Extract the annotations from the IgBLAST output using the sample ID and contig IDs
    subset <- igblast_annotations[[sample_id]][igblast_annotations[[sample_id]]$sequence_id %in% c(VDJ_contig_id, VJ_contig_id), columns]
    
    # Remove the rownames from the 'subset' dataframe
    rownames(subset) <- NULL
    
    # Rename the columns containing the gene symbols
    colnames(subset)[colnames(subset) %in% c("v_call", "d_call", "j_call", "c_call")] <- c("vgene", "dgene", "jgene", "cgene")
    
    # Rename the columns containing the CIGAR strings
    colnames(subset)[colnames(subset) %in% c("v_cigar", "d_cigar", "j_cigar")] <- c("vgene_cigar", "dgene_cigar", "jgene_cigar")
    
    # Rename the columns containing the start and end position of the genes in the recovered sequence
    colnames(subset) <- gsub(pattern = "_sequence_(start|end)$", replacement = "gene_\\1", colnames(subset))
    
    # Rename the columns containing the start and end position of the genes in the germline sequence
    colnames(subset) <- gsub(pattern = "_germline_(start|end)$", replacement = "gene_germline_\\1", colnames(subset))
    
    # Append '_nt' to the column names of the columns containing DNA sequences
    colnames(subset)[colnames(subset) %in% c("sequence_alignment", "germline_alignment",  "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4", "junction")] <- paste0(c("sequence_alignment", "germline_alignment", "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4", "junction"), "_nt")
    
    # Translate the columns with DNA sequences into protein sequences using the 'translate_IMGT_sequence' function (if 'translate.IMGT.sequences' is set to TRUE)
    subset[, c("sequence_alignment_aa", "germline_alignment_aa", "fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa", "cdr3_aa", "fwr4_aa")] <- sapply(c("sequence_alignment_nt", "germline_alignment_nt", "fwr1_nt", "cdr1_nt", "fwr2_nt", "cdr2_nt", "fwr3_nt", "cdr3_nt", "fwr4_nt"), function(x) sapply(1:nrow(subset), function(y) translate_IMGT_sequence(subset[y, x])))
    
    # Transform the 'subset' dataframe into a 1x46 dataframe and append 'VDJ_' and 'VJ_' as prefixes to the column names (if one of the contig IDs could not be found in the IgBLAST output, the corresponding columns will be set to NA)
    subset_VDJ <- subset[subset$sequence_id == VDJ_contig_id,]; colnames(subset_VDJ) <- paste0("VDJ_", colnames(subset_VDJ)); if(nrow(subset_VDJ) == 0){subset_VDJ[1,] <- NA}
    subset_VJ <- subset[subset$sequence_id == VJ_contig_id,]; colnames(subset_VJ) <- paste0("VJ_", colnames(subset_VJ)); if(nrow(subset_VJ) == 0){subset_VJ[1,] <- NA}
    cell_annotations <- cbind(subset_VDJ, subset_VJ)
    
    # Remove the 'VJ_dgene' column
    cell_annotations[, grep(pattern = "VJ_dgene", colnames(cell_annotations), value = TRUE)] <- NULL
    
    # Return the 'cell_annotations' dataframe
    return(cell_annotations)
  }
  
  
  # Create list of 'sample_id', 'VDJ_contig_id', and 'VJ_contig_id' for each row/cell in the VDJ dataframe
  cell_info_list <- lapply(1:nrow(VDJ), function(x) list(sample_id = VDJ[x, "sample_id"], VDJ_contig_id = VDJ[x, "VDJ_contig_id"], VJ_contig_id = VDJ[x, "VJ_contig_id"]))
  
  # Retrieve all annotations for each cell using the 'cell_info_list' and the 'retrieve_cell_annotations' function, and store in the 'all_IgBLAST_annotations' dataframe
  all_IgBLAST_annotations <- do.call(rbind, lapply(cell_info_list, function(x) retrieve_cell_annotations(igblast_annotations = igblast_annotations, cell_info = x, translate.IMGT.sequences = translate.IMGT.sequences)))
  
  
  # If the 'method' parameter is set to 'replace', the original annotation columns in the VDJ dataframe will get the suffix '_10x'
  if(method == "replace"){
    
    # The original columns in 'annotation_columns' will get the suffix '_10x'
    colnames(VDJ)[colnames(VDJ) %in% annotation_columns] <- paste0(colnames(VDJ)[colnames(VDJ) %in% annotation_columns], "_10x")
    
    # Append the IgBLAST annotation columns to the VDJ dataframe
    VDJ[, annotation_columns] <- all_IgBLAST_annotations[, annotation_columns]
    
    # Append the alignment columns to the VDJ dataframe
    VDJ[, alignment_columns] <- all_IgBLAST_annotations[, alignment_columns]
    
    # Redefine the order of the columns in the VDJ dataframe
    new_order <- c("sample_id", "barcode", "celltype", "isotype",
                   "VDJ_contig_id", "VJ_contig_id",
                   "VDJ_consensus_id", "VJ_consensus_id",
                   "clonotype_id", "clonotype_frequency",
                   
                   "VDJ_chain", "VDJ_chain_count", "VDJ_umis",
                   "VDJ_vgene", "VDJ_vgene_10x", "VDJ_vgene_cigar", "VDJ_vgene_start", "VDJ_vgene_end", "VDJ_vgene_germline_start", "VDJ_vgene_germline_end",
                   "VDJ_dgene", "VDJ_dgene_10x", "VDJ_dgene_cigar", "VDJ_dgene_start", "VDJ_dgene_end", "VDJ_dgene_germline_start", "VDJ_dgene_germline_end",
                   "VDJ_jgene", "VDJ_jgene_10x", "VDJ_jgene_cigar", "VDJ_jgene_start", "VDJ_jgene_end", "VDJ_jgene_germline_start", "VDJ_jgene_germline_end",
                   "VDJ_cgene", "VDJ_cgene_10x",
                   "VDJ_fwr1_nt", "VDJ_fwr1_nt_10x", "VDJ_fwr1_aa", "VDJ_fwr1_aa_10x",
                   "VDJ_cdr1_nt", "VDJ_cdr1_nt_10x", "VDJ_cdr1_aa", "VDJ_cdr1_aa_10x", 
                   "VDJ_fwr2_nt", "VDJ_fwr2_nt_10x", "VDJ_fwr2_aa", "VDJ_fwr2_aa_10x",
                   "VDJ_cdr2_nt", "VDJ_cdr2_nt_10x", "VDJ_cdr2_aa", "VDJ_cdr2_aa_10x",
                   "VDJ_fwr3_nt", "VDJ_fwr3_nt_10x", "VDJ_fwr3_aa", "VDJ_fwr3_aa_10x",
                   "VDJ_cdr3_nt", "VDJ_cdr3_nt_10x", "VDJ_cdr3_aa", "VDJ_cdr3_aa_10x",
                   "VDJ_fwr4_nt", "VDJ_fwr4_nt_10x", "VDJ_fwr4_aa", "VDJ_fwr4_aa_10x",
                   "VDJ_junction_nt", "VDJ_junction_aa",
                   "VDJ_sequence_nt_raw", "VDJ_sequence_nt_trimmed", "VDJ_sequence_aa_trimmed",
                   "VDJ_consensus_nt_raw", "VDJ_consensus_nt_trimmed", "VDJ_consensus_aa_trimmed",
                   "VDJ_germline_nt_raw", "VDJ_germline_nt_trimmed", "VDJ_germline_aa_trimmed", "VDJ_germline_nt_cdr3_filled", "VDJ_germline_aa_cdr3_filled",
                   "VDJ_sequence_alignment_nt", "VDJ_sequence_alignment_aa", "VDJ_germline_alignment_nt", "VDJ_germline_alignment_aa",
                   
                   "VJ_chain", "VJ_chain_count", "VJ_umis",
                   "VJ_vgene", "VJ_vgene_10x", "VJ_vgene_cigar", "VJ_vgene_start", "VJ_vgene_end", "VJ_vgene_germline_start", "VJ_vgene_germline_end",
                   "VJ_jgene", "VJ_jgene_10x", "VJ_jgene_cigar", "VJ_jgene_start", "VJ_jgene_end", "VJ_jgene_germline_start", "VJ_jgene_germline_end", 
                   "VJ_cgene", "VJ_cgene_10x",
                   "VJ_fwr1_nt", "VJ_fwr1_nt_10x", "VJ_fwr1_aa", "VJ_fwr1_aa_10x",
                   "VJ_cdr1_nt", "VJ_cdr1_nt_10x", "VJ_cdr1_aa", "VJ_cdr1_aa_10x", 
                   "VJ_fwr2_nt", "VJ_fwr2_nt_10x", "VJ_fwr2_aa", "VJ_fwr2_aa_10x",
                   "VJ_cdr2_nt", "VJ_cdr2_nt_10x", "VJ_cdr2_aa", "VJ_cdr2_aa_10x",
                   "VJ_fwr3_nt", "VJ_fwr3_nt_10x", "VJ_fwr3_aa", "VJ_fwr3_aa_10x",
                   "VJ_cdr3_nt", "VJ_cdr3_nt_10x", "VJ_cdr3_aa", "VJ_cdr3_aa_10x",
                   "VJ_fwr4_nt", "VJ_fwr4_nt_10x", "VJ_fwr4_aa", "VJ_fwr4_aa_10x",
                   "VJ_junction_nt", "VJ_junction_aa",
                   "VJ_sequence_nt_raw", "VJ_sequence_nt_trimmed", "VJ_sequence_aa_trimmed",
                   "VJ_consensus_nt_raw", "VJ_consensus_nt_trimmed", "VJ_consensus_aa_trimmed",
                   "VJ_germline_nt_raw", "VJ_germline_nt_trimmed", "VJ_germline_aa_trimmed", "VJ_germline_nt_cdr3_filled", "VJ_germline_aa_cdr3_filled",
                   "VJ_sequence_alignment_nt", "VJ_sequence_alignment_aa", "VJ_germline_alignment_nt", "VJ_germline_alignment_aa")
    
    # Retrieve the column names from the 'VDJ' dataframe that are not yet included in the 'new_order' vector and add these column names at the end of the vector
    unrecognized_columns <- colnames(VDJ)[!colnames(VDJ) %in% new_order]
    new_order <- c(new_order, unrecognized_columns)
    
    # Reorder the VDJ dataframe
    VDJ <- VDJ[, new_order]
  }
  
  # If the 'method' parameter is set to 'append', the IgBLAST annotations are stored in columns with the suffix '_IgBLAST'
  if(method == "append"){
    
    # Append the IgBLAST annotations columns with the suffix '_IgBLAST'
    VDJ[, paste0(annotation_columns, "_IgBLAST")] <- all_IgBLAST_annotations[, annotation_columns]
    
    # Append the alignment columns to the VDJ dataframe
    VDJ[, paste0(alignment_columns, "_IgBLAST")] <- all_IgBLAST_annotations[, alignment_columns]
    
    # Redefine the order of the columns in the VDJ dataframe
    new_order <- c("sample_id", "barcode", "celltype", "isotype",
                   "VDJ_contig_id", "VJ_contig_id",
                   "VDJ_consensus_id", "VJ_consensus_id",
                   "clonotype_id", "clonotype_frequency",
                   
                   "VDJ_chain", "VDJ_chain_count", "VDJ_umis",
                   "VDJ_vgene", "VDJ_vgene_IgBLAST", "VDJ_vgene_cigar_IgBLAST", "VDJ_vgene_start_IgBLAST", "VDJ_vgene_end_IgBLAST", "VDJ_vgene_germline_start_IgBLAST", "VDJ_vgene_germline_end_IgBLAST",
                   "VDJ_dgene", "VDJ_dgene_IgBLAST", "VDJ_dgene_cigar_IgBLAST", "VDJ_dgene_start_IgBLAST", "VDJ_dgene_end_IgBLAST", "VDJ_dgene_germline_start_IgBLAST", "VDJ_dgene_germline_end_IgBLAST",
                   "VDJ_jgene", "VDJ_jgene_IgBLAST", "VDJ_jgene_cigar_IgBLAST", "VDJ_jgene_start_IgBLAST", "VDJ_jgene_end_IgBLAST", "VDJ_jgene_germline_start_IgBLAST", "VDJ_jgene_germline_end_IgBLAST",
                   "VDJ_cgene", "VDJ_cgene_IgBLAST",
                   "VDJ_fwr1_nt", "VDJ_fwr1_nt_IgBLAST", "VDJ_fwr1_aa", "VDJ_fwr1_aa_IgBLAST",
                   "VDJ_cdr1_nt", "VDJ_cdr1_nt_IgBLAST", "VDJ_cdr1_aa", "VDJ_cdr1_aa_IgBLAST", 
                   "VDJ_fwr2_nt", "VDJ_fwr2_nt_IgBLAST", "VDJ_fwr2_aa", "VDJ_fwr2_aa_IgBLAST",
                   "VDJ_cdr2_nt", "VDJ_cdr2_nt_IgBLAST", "VDJ_cdr2_aa", "VDJ_cdr2_aa_IgBLAST",
                   "VDJ_fwr3_nt", "VDJ_fwr3_nt_IgBLAST", "VDJ_fwr3_aa", "VDJ_fwr3_aa_IgBLAST",
                   "VDJ_cdr3_nt", "VDJ_cdr3_nt_IgBLAST", "VDJ_cdr3_aa", "VDJ_cdr3_aa_IgBLAST",
                   "VDJ_fwr4_nt", "VDJ_fwr4_nt_IgBLAST", "VDJ_fwr4_aa", "VDJ_fwr4_aa_IgBLAST",
                   "VDJ_junction_nt_IgBLAST", "VDJ_junction_aa_IgBLAST",
                   "VDJ_sequence_nt_raw", "VDJ_sequence_nt_trimmed", "VDJ_sequence_aa_trimmed",
                   "VDJ_consensus_nt_raw", "VDJ_consensus_nt_trimmed", "VDJ_consensus_aa_trimmed",
                   "VDJ_germline_nt_raw", "VDJ_germline_nt_trimmed", "VDJ_germline_aa_trimmed", "VDJ_germline_nt_cdr3_filled", "VDJ_germline_aa_cdr3_filled",
                   "VDJ_sequence_alignment_nt_IgBLAST", "VDJ_sequence_alignment_aa_IgBLAST", "VDJ_germline_alignment_nt_IgBLAST", "VDJ_germline_alignment_aa_IgBLAST",
                   
                   "VJ_chain", "VJ_chain_count", "VJ_umis",
                   "VJ_vgene", "VJ_vgene_IgBLAST", "VJ_vgene_cigar_IgBLAST", "VJ_vgene_start_IgBLAST", "VJ_vgene_end_IgBLAST", "VJ_vgene_germline_start_IgBLAST", "VJ_vgene_germline_end_IgBLAST",
                   "VJ_jgene", "VJ_jgene_IgBLAST", "VJ_jgene_cigar_IgBLAST", "VJ_jgene_start_IgBLAST", "VJ_jgene_end_IgBLAST", "VJ_jgene_germline_start_IgBLAST", "VJ_jgene_germline_end_IgBLAST",
                   "VJ_cgene", "VJ_cgene_IgBLAST",
                   "VJ_fwr1_nt", "VJ_fwr1_nt_IgBLAST", "VJ_fwr1_aa", "VJ_fwr1_aa_IgBLAST",
                   "VJ_cdr1_nt", "VJ_cdr1_nt_IgBLAST", "VJ_cdr1_aa", "VJ_cdr1_aa_IgBLAST", 
                   "VJ_fwr2_nt", "VJ_fwr2_nt_IgBLAST", "VJ_fwr2_aa", "VJ_fwr2_aa_IgBLAST",
                   "VJ_cdr2_nt", "VJ_cdr2_nt_IgBLAST", "VJ_cdr2_aa", "VJ_cdr2_aa_IgBLAST",
                   "VJ_fwr3_nt", "VJ_fwr3_nt_IgBLAST", "VJ_fwr3_aa", "VJ_fwr3_aa_IgBLAST",
                   "VJ_cdr3_nt", "VJ_cdr3_nt_IgBLAST", "VJ_cdr3_aa", "VJ_cdr3_aa_IgBLAST",
                   "VJ_fwr4_nt", "VJ_fwr4_nt_IgBLAST", "VJ_fwr4_aa", "VJ_fwr4_aa_IgBLAST",
                   "VJ_junction_nt_IgBLAST", "VJ_junction_aa_IgBLAST",
                   "VJ_sequence_nt_raw", "VJ_sequence_nt_trimmed", "VJ_sequence_aa_trimmed",
                   "VJ_consensus_nt_raw", "VJ_consensus_nt_trimmed", "VJ_consensus_aa_trimmed",
                   "VJ_germline_nt_raw", "VJ_germline_nt_trimmed", "VJ_germline_aa_trimmed", "VJ_germline_nt_cdr3_filled", "VJ_germline_aa_cdr3_filled",
                   "VJ_sequence_alignment_nt_IgBLAST", "VJ_sequence_alignment_aa_IgBLAST", "VJ_germline_alignment_nt_IgBLAST", "VJ_germline_alignment_aa_IgBLAST")
    
    # Retrieve the column names from the 'VDJ' dataframe that are not yet included in the 'new_order' vector and add these column names at the end of the vector
    unrecognized_columns <- colnames(VDJ)[!colnames(VDJ) %in% new_order]
    print(unrecognized_columns)
    new_order <- c(new_order, unrecognized_columns)
    
    # Reorder the VDJ dataframe 
    VDJ <- VDJ[, new_order]
  }
  
  # Return the VDJ dataframe with the appended IgBLAST annotations and alignments
  return(VDJ)
}