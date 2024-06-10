#' Minimal version of the VDJ building part from VDJ_GEX_matrix() function. Optimized for for Cell Ranger v7 and suitable for older Cell Ranger versions.
#' Authors: Valentijn Tromp, Daphne van Ginneken, Tudor-Stefan Cotet, Victor Kreiner, Aurora Desideri Perea, Evgenios Kladis, Anamay Samant
#' @description  Imports Cell Ranger output into R dataframe for downstream analyses. Minimal version of the VDJ building part from the 'VDJ_GEX_matrix()' function of Platypus package. Adapted for Cell Ranger v7 and older versions as well. Seurat objects can be integrated by matching barcodes from the Seurat object's metadata with the barcodes in the 'barcode' column of the VDJ dataframe.
#' @param VDJ.directory string - path to parent directory containing the output folders (one folder for each sample) of Cell Ranger. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires. ! Neccessary 5 files within this folder: 'filtered_contig_annotations.csv', 'filtered_contig.fasta', 'consensus_annotations.csv', 'consensus.fasta', and 'concat_ref.fasta'.
#' @param VDJ.sample.list list - list of paths to the output folders (one folder for each sample) of Cell Ranger. The names of the items in the list will be used as sample names in the VDJ dataframe (instead of the name of the sample/run directory). This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires. ! Neccessary 5 files within this folder: 'filtered_contig_annotations.csv', 'filtered_contig.fasta', 'consensus_annotations.csv', 'consensus.fasta', and 'concat_ref.fasta'.
#' @param remove.divergent.cells bool - if TRUE, cells with more than one VDJ transcript or more than one VJ transcript will be excluded. This could be due to multiple cells being trapped in one droplet or due to light chain dual expression (concerns ~0.2-0.5% of B cells, see DOI:10.1084/jem.181.3.1245). Defaults to FALSE.
#' @param complete.cells.only bool - if TRUE, only cells with both a VDJ transcripts and a VJ transcript are included in the VDJ dataframe. Keeping only cells with 1 VDJ and 1 VJ transcript could be preferable for downstream analysis. Defaults to FALSE.
#' @param trim.germlines bool - if TRUE, the raw germline sequences of each clone will be trimmed using the the consensus sequences of that clone as reference seqeunces (using Biostrings::pairwiseAlignment with the option "global-local" and a gap opening cost = gap.opening.cost). Defaults to FALSE.
#' @param fill.germline.CDR3 bool - if TRUE, the CDR3 region in the trimmed germline sequences are replaced by the most frequently observed CDR3 sequences within the clone in order to obtain germline sequences that are more likely to encode producible and productive antibodies.
#' @param gap.opening.cost float or Inf - the cost for opening a gap in 'Biostrings::pairwiseAlignment()' when aligning and trimming germline sequences. If set to Inf, a gapless alignment will be performed. Defaults to 10.
#' @param gap.extension.cost float or Inf - the cost for extending a gap in 'Biostrings::pairwiseAlignment()' when aligning and trimming germline sequences. If set to Inf, a gapless alignment will be performed. Defaults to 4.
#' @param parallel bool - if TRUE, the per-sample VDJ building is executed in parallel (parallelized across samples). Defaults to FALSE.
#' @param num.cores integer - number of cores to be used when parallel = TRUE. Defaults to all available cores - 1 or the number of sample folders in 'VDJ.directory' (depending which number is smaller).
#' @return Returns the VDJ dataframe / VGM[[1]] object. Each row in this dataframe represent one cell, or one unique cell barcode.
#' @export
#' @examples
#' \dontrun{
#' VDJ <- VDJ_build(VDJ_directory)
#'}


VDJ_build <- function(VDJ.directory,
                      VDJ.sample.list,
                      remove.divergent.cells,
                      complete.cells.only,
                      trim.germlines,
                      gap.opening.cost,
                      gap.extension.cost,
                      fill.germline.CDR3,
                      parallel,
                      num.cores){


  # If both the 'VDJ.directory' and 'VDJ.sample.list' parameters are not specified, a message is returned and execution is stopped
  if(missing(VDJ.directory) & missing(VDJ.sample.list)){stop('The path to the parent directory containing the output folders (one for each sample) from Cell Ranger, or the list of these sample directories are missing.')}
  # If both a 'VDJ.directory' and 'VDJ.sample.list' are specified, a message is returned and execution is stopped
  if(!missing(VDJ.directory) & !missing(VDJ.sample.list)){stop('Both a path to the parent directory containing the output folders (one for each sample) from Cell Ranger and a list of sample directories are given as input. Please provide one parent direcotry or one list of sample directories.')}
  # Create list with paths to output folders of Cell Ranger (one path for each sample) in parent dirctory 'VDJ.directory'
  if(missing(VDJ.sample.list)){VDJ.sample.list <- list.dirs(path = VDJ.directory, full.names = TRUE, recursive = FALSE); names(VDJ.sample.list) <- basename(VDJ.sample.list)}
  # Add a '/' at the end of each directory (if not already present)
  VDJ.sample.list <- sub("/?$", "/", VDJ.sample.list)
  # If not specified, the 'complete.cells.only' parameter is set to FALSE
  if(missing(complete.cells.only)){complete.cells.only <- FALSE}
  # If not specified, the 'remove.divergent.cells' parameter is set to FALSE
  if(missing(remove.divergent.cells)){remove.divergent.cells <- FALSE}
  # If not specified, the 'trim.germlines' parameter is set to FALSE
  if(missing(trim.germlines)){trim.germlines <- FALSE}
  # If the 'trim.germlines' parameter is set to FALSE, but the 'gap.opening.cost', 'gap.extension.cost'. or 'fill.germline.CDR3' parameters are specified, a message is returned and execution is stopped
  if(!trim.germlines && (!missing(gap.opening.cost) | !missing(gap.extension.cost) | !missing(fill.germline.CDR3))){stop("When 'trim.germlines' is set to FALSE, the following parameters should not be specified: 'gap.opening.cost', 'gap.extension.cost', and 'fill.germline.CDR3'.")}
  # If not specified, the 'fill.germline.CDR3' parameter is set to FALSE
  if(missing(fill.germline.CDR3)){fill.germline.CDR3 <- FALSE}
  # If not specified, the 'gap.opening.cost' parameter is set to 10
  if(trim.germlines && missing(gap.opening.cost)){gap.opening.cost <- 10}
  # If not specified, the 'gap.extension.cost' parameter is set to 4
  if(trim.germlines && missing(gap.extension.cost)){gap.extension.cost <- 4}
  # If not specified, the 'fill.germline.CDR3' parameter is set to FALSE
  if(trim.germlines && missing(fill.germline.CDR3)){fill.germline.CDR3 <- FALSE}
  # If not specified, the 'parallel' parameter is set to FALSE
  if(missing(parallel)) parallel <- FALSE
  # If the number of cores is not specified, while 'parallel' is set to TRUE, the number is set to the number of CPU cores on the current host, or the number of samples in 'VDJ.sample.list', if this number is lower
  if(parallel && missing(num.cores)){num.cores <- min(c((parallel::detectCores() - 1), length(VDJ.sample.list)))}
  # If 'parallel' is set to FALSE, the 'num.cores' parameter is set to FALSE
  if(!parallel && missing(num.cores)){num.cores <- NULL}
  # Return a message informing that performing pairwise alignment for all germline sequences significantly increases computation time (if 'trim.germlines' is set to TRUE)
  if(trim.germlines){message("WARNING: Setting 'trim.germlines' to TRUE will significantly increase computation time, please be patient!\n")}

  # Set all global variables to NULL - CRAN fixes
  sequence_trimmed <- NULL
  chain <- NULL
  barcode <- NULL
  VDJ_clonotype_id <- NULL
  VJ_clonotype_id <- NULL
  VDJ_chain_count <- NULL
  VJ_chain_count <- NULL
  clonotype_id <- NULL
  
  
  get_trim_seqs_from_json <- function(annotation_file){
    
    # Obtains trimmed sequences from raw sequences using the start and end indices of the annotated 'L-REGION+V-REGION' and 'J-REGION' using 'all_contig_annotations.json' or 'consensus_annotations.json' file from Cell Ranger
    # Arguments:
    # - annotations_JSON: path to JSON file containing annotations of contig or consensus sequences ('all_contig_annotations.json' or 'consensus_annotations.json')
    # Authors: Victor Kreiner, Tudor-Stefan Cotet, Valentijn Tromp

    # Read JSON file and store data  object in 'annotations_list'
    # The annotations list, containing a sublist for each contig/consensus sequence, is eventually transformed into a dataframe in which each row represents one trimmed contig/consensus sequence
    # A sequence may contain the following annotation regions: '5-UTR', 'L-REGION+V-REGION', 'D-REGION', 'J-REGION', 'C-REGION' (which are listed in 'annotations')
    annotations_list <- jsonlite::read_json(annotation_file)
    
    # Define function to extract the following info for each contig/consensus sequence:
    # - id: contig or consensus ID depending on the JSON file provided ('all_contig_annotations.json' or 'consensus_annotations.json', respectively)
    # - start: the start index of the 'L-region+V-region' in the raw sequence
    # - stop: the stop index of the 'L-region+V-region' in the raw sequence
    # - sequence trimmed: the trimmed sequence (given the start and stop indices)
    extract_info <- function(single_annotation){
      
      # Check if both 'L-REGION+V-REGION' and 'J-REGION' annotations are present
      if (all(c("L-REGION+V-REGION", "J-REGION") %in% sapply(single_annotation$annotations, function(x) x$feature$region_type))){
        
        # Identify 'L-REGION+V-REGION' annotation (if present)
        lv_region <- single_annotation$annotations[[which(sapply(single_annotation$annotations, function(x) x$feature$region_type == "L-REGION+V-REGION"))]]
        # Identify 'J-REGION' annotation (if present)
        j_region <- single_annotation$annotations[[which(sapply(single_annotation$annotations, function(x) x$feature$region_type == "J-REGION"))]]
        
        # Create dataframe 'trimmed_sequence' that contains the contig/consensus ID and the trimmed sequence
        trimmed_sequence <- data.frame(
          # Retrieve contig/consensus ID
          id = single_annotation$contig_name,
          # Trim raw sequence using 'L-REGION+V-REGION' start and 'J-REGION' end indices
          sequence_trimmed = substr(single_annotation$sequence,
                                    start = as.numeric(lv_region$contig_match_start)+1,
                                    stop = as.numeric(j_region$contig_match_end))
        )
        
        # Return 'trimmed_sequence' dataframe
        return(trimmed_sequence)
      }
    } 
    
    # Use 'lapply' to apply the 'extract_info' function to each element of 'annotations_list' and store in 'trimmed_sequences' list
    trimmed_sequences <- lapply(annotations_list, extract_info)
    # Transform the 'trimmed_sequences' list into a dataframe by 'rbind'
    trimmed_sequences <- do.call(rbind, trimmed_sequences[!sapply(trimmed_sequences, is.null)])
    # Remove "-1" from contig/consensus IDs
    trimmed_sequences$id <- gsub(trimmed_sequences$id, pattern = "-1", replacement = "")
    
    # Return the 'trimmed_sequences' dataframe
    return(trimmed_sequences)
  }
  
  
  translate_DNA <- function(sequence){
    
    # Translates a DNA sequence with gaps ('-') into a protein sequence while maintaining the appropriate reading frame
    # Arguments:
    # - sequence: nucleotide sequence to be translated
    # Authors: Anamay Samant, Daphne van Ginneken, Valentijn Tromp
    
    # If an empty string is provided, return an empty string
    if(sequence == ""){
      return("")
    }
    
    # Import codon table
    genetic_code <- list(
      "AAA"="K", "AAC"="N", "AAG"="K", "AAT"="N",
      "ACA"="T", "ACC"="T", "ACG"="T", "ACT"="T",
      "AGA"="R", "AGC"="S", "AGG"="R", "AGT"="S",
      "ATA"="I", "ATC"="I", "ATG"="M", "ATT"="I",
      "CAA"="Q", "CAC"="H", "CAG"="Q", "CAT"="H",
      "CCA"="P", "CCC"="P", "CCG"="P", "CCT"="P",
      "CGA"="R", "CGC"="R", "CGG"="R", "CGT"="R",
      "CTA"="L", "CTC"="L", "CTG"="L", "CTT"="L",
      "GAA"="E", "GAC"="D", "GAG"="E", "GAT"="D",
      "GCA"="A", "GCC"="A", "GCG"="A", "GCT"="A",
      "GGA"="G", "GGC"="G", "GGG"="G", "GGT"="G",
      "GTA"="V", "GTC"="V", "GTG"="V", "GTT"="V",
      "TAA"="*", "TAC"="Y", "TAG"="*", "TAT"="Y",
      "TCA"="S", "TCC"="S", "TCG"="S", "TCT"="S",
      "TGA"="*", "TGC"="C", "TGG"="W", "TGT"="C",
      "TTA"="L", "TTC"="F", "TTG"="L", "TTT"="F"
    )
    
    # Split the nucleotide sequence into codons
    codons <- sapply(1:(round(nchar(sequence)/3)), function(x) substr(sequence, x*3-2, x*3))
    
    # Create an empty string to store the protein sequence
    aa_seq <- ""
    
    # Iterate through the codons in the 'codons' vector
    for(i in codons){

      # Codons that contain a gap ("-") are translated into a gap ("-") in the protein sequence
      if(grepl(pattern = "-", i, fixed = TRUE)){aa_seq <- paste0(aa_seq, "-")}
      
      # Otherwise, codons are translated using the 'genetic_code' list
      else{aa_seq <- paste0(aa_seq, genetic_code[[i]])}
    }
    
    # Return the protein sequence
    return(aa_seq)
  }


  trim_seq <- function(seq1, 
                       seq2, 
                       gap.opening.cost, 
                       gap.extension.cost,
                       sequence.type){
    
    # Perform global-local alignment with seq1 as pattern and seq2 as subject and thereby trim seq2; this type of alignment aims to find the best alignment over the entire length of seq1 on seq2
    # Arguments:
    # - seq1: sequence used as reference during pairwise alignment
    # - seq2: sequence that will be trimmed during pairwise alignment
    # - gap.opening.cost: cost for opening a gap in 'Biostrings::pairwiseAlignment' 
    # - gap.extension.cost: cost for extending a gap in 'Biostrings::pairwiseAlignment'
    # - sequence.type: sequence type of the sequences that are being aligned and trimmed ('DNA' or 'AA')
    # Authors: Evgenios Kladis, Valentijn Tromp

    # If 'seq1' or 'seq2' is missing, an empty string is returned
    if(is.na(seq1) | is.na(seq2)){
      return ("")
    } else{
      
    # When the sequence.type is set to 'DNA', perform pairwise alignment with two objects of class 'DNAString'
    if(sequence.type == "DNA"){
      alignment <- Biostrings::pairwiseAlignment(pattern = Biostrings::DNAString(seq1),
                                                 subject = Biostrings::DNAString(seq2),
                                                 type = "global-local",
                                                 gapOpening = gap.opening.cost,
                                                 gapExtension = gap.extension.cost)
    }
    
    # When the sequence.type is set to 'DNA', perform pairwise alignment with two objects of class 'AAString'
    if(sequence.type == "AA"){
      alignment <- Biostrings::pairwiseAlignment(pattern = Biostrings::AAString(seq1),
                                                 subject = Biostrings::AAString(seq2),
                                                 type = "global-local",
                                                 gapOpening = gap.opening.cost,
                                                 gapExtension = gap.extension.cost)
    }
    
    # Return the subject of the alignment, which is trimmed seq2
    return(as.character(alignment@subject))
    }
  }


  make_VDJ_sample <- function(sample.ID,
                              VDJ.sample.list,
                              remove.divergent.cells,
                              complete.cells.only,
                              trim.germlines,
                              gap.opening.cost,
                              gap.extension.cost,
                              fill.germline.CDR3){

    # Obtains all necessary data from a VDJ directory to construct the VDJ dataframe for one sample
    # Arguments:
    # - sample.ID: string -
    # - VDJ.sample.list: list - 
    # - remove.divergent.cells: bool - if TRUE, cells with more than one VDJ transcript or more than one VJ transcript will be excluded. This could be due to multiple cells being trapped in one droplet or due to light chain dual expression (concerns ~2-5% of B cells, see DOI:10.1084/jem.181.3.1245). Defaults to FALSE.
    # - complete.cells.only: bool - if TRUE, only cells with both a VDJ transcripts and a VJ transcript are included in the VDJ dataframe. Keeping only cells with 1 VDJ and 1 VJ transcript could be preferable for downstream analysis. Defaults to FALSE.
    # - trim.germlines: bool - if TRUE, the raw germline sequences of each clone will be trimmed using the the consensus sequences of that clone as reference seqeunces (using BIostrings::pairwiseAlignment with the option "global-local" and a gap opening cost = gap.opening.cost). Defaults to FALSE.
    # - gap.opening.cost: float or Inf - the cost for opening a gap in 'Biostrings::pairwiseAlignment' when aligning and trimming germline sequences. 
    # - gap.extension.cost: float or Inf - the cost for extending a gap in 'Biostrings::pairwiseAlignment' when aligning and trimming germline sequences.
    # - fill.germline.CDR3: bool - if TRUE, the CDR3 regions in the trimmed germline sequences are replaced by the most frequently observed CDR3 sequences within the clone
    # Authors: Tudor-Stefan Cotet, Aurora Desideri Perea, Anamay Samant, Valentijn Tromp
    
    # Retrieve the path to the sample/run directory from the 'VDJ.sample.list' using the 'sample.ID'
    VDJ.sample.directory <- VDJ.sample.list[[sample.ID]]
    
    # 1. Read in the following files from Cell Ranger output and store in list 'out':
    # [[1]]: 'filtered_contig_annotations.csv' contains high-level annotations of each high-confidence contig from cell-associated barcodes
    # [[2]]: 'filtered_contig.fasta' contains contig sequences from barcodes that passed the filtering steps and are annotated in 'filtered_contig_annotations.csv'
    # [[3]]: 'consensus_annotations.csv' contains high-level annotations of each clonotype consensus sequence
    # [[4]]: 'consensus.fasta' contains the consensus sequences of each assembled contig, which is identical to the sequence of the most frequent exact subclonotype
    # [[5]]: 'concat_ref.fasta' contains the concatenated V(D)J reference segments for the segments detected on each consensus sequence and serves as an approximate reference for each consensus sequence
    
    # Check if the required files are present in the sample directory, otherwise, a message is returned and execution is stopped
    if(!file.exists(paste0(VDJ.sample.directory,"filtered_contig_annotations.csv"))){stop(paste0("Could not find the 'filtered_contig_annotations.csv' file in ", VDJ.sample.directory))}
    if(!file.exists(paste0(VDJ.sample.directory,"filtered_contig.fasta"))){stop(paste0("Could not find the 'filtered_contig.fasta' file in ", VDJ.sample.directory))}
    if(!file.exists(paste0(VDJ.sample.directory,"consensus_annotations.csv"))){stop(paste0("Could not find the 'consensus_annotations.csv' file in ", VDJ.sample.directory))}
    if(!file.exists(paste0(VDJ.sample.directory,"consensus.fasta"))){stop(paste0("Could not find the 'consensus.fasta' file in ", VDJ.sample.directory))}
    if(!file.exists(paste0(VDJ.sample.directory,"concat_ref.fasta"))){stop(paste0("Could not find the 'concat_ref.fasta' file in ", VDJ.sample.directory))}
    
    # Create empty list to store warnings during data retrieval out of single sample directory
    warnings <- list()
    
    # Create dataframe to store number of cells/barcodes filtered when 'remove.divergent.cells' and/or 'complete.cells.only' is/are set to TRUE (all counts default to 0)
    filtering_counts <- data.frame(divergent.cells = 0, incomplete.cells = 0, row.names = basename(VDJ.sample.directory))
    
    # Create empty list ('out') with length 5
    out <- vector("list", 5)

    # Read 'filtered_contig_annotations.csv' file and store dataframe in out[[1]]
    out[[1]] <- utils::read.csv(paste(VDJ.sample.directory, "filtered_contig_annotations.csv", sep = ""), sep = ",", header = TRUE)
    # Add "_aa" to the names of columns containing protein sequences for clarity
    colnames(out[[1]]) <- gsub(colnames(out[[1]]), pattern = "^(fwr\\d|cdr\\d)$", replacement = "\\1_aa")
    # Remove "raw_" from column names for consistency
    colnames(out[[1]]) <- gsub(colnames(out[[1]]), pattern = "raw_", replacement = "")
    # Remove "-1" from cell barcodes and contig IDs
    out[[1]]$barcode <- gsub(out[[1]]$barcode, pattern = "-1", replacement = "")
    out[[1]]$contig_id <- gsub(out[[1]]$contig_id, pattern = "-1", replacement = "")
    # Remove "_" from 'v_gene', 'd_gene', 'j_gene', and 'c_gene' column
    colnames(out[[1]]) <- gsub(colnames(out[[1]]), pattern = "_gene", replacement = "gene")
    # Concatenate nt and aa sequences of FWR1-4 and CDR1-3 regions and place concatenated sequence in new column in dataframe in out[[1]] (if the FWR1-4 and CDR1-3 regions are annotated in the 'filtered_contig_annotations.csv' file)
    if(all(c("fwr1_nt","cdr1_nt","fwr2_nt","cdr2_nt","fwr3_nt","cdr3_nt","fwr4_nt","fwr1_aa","cdr1_aa","fwr2_aa","cdr2_aa","fwr3_aa","cdr3_aa","fwr4_aa") %in% colnames(out[[1]]))){
      out[[1]]$sequence_nt_trimmed <- paste0(out[[1]][["fwr1_nt"]], out[[1]][["cdr1_nt"]],
                                             out[[1]][["fwr2_nt"]], out[[1]][["cdr2_nt"]],
                                             out[[1]][["fwr3_nt"]], out[[1]][["cdr3_nt"]],
                                             out[[1]][["fwr4_nt"]])
      out[[1]]$sequence_aa_trimmed <- paste0(out[[1]][["fwr1_aa"]], out[[1]][["cdr1_aa"]],
                                             out[[1]][["fwr2_aa"]], out[[1]][["cdr2_aa"]],
                                             out[[1]][["fwr3_aa"]], out[[1]][["cdr3_aa"]],
                                             out[[1]][["fwr4_aa"]])
    # If the FWR1-4 and CDR1-3 regions are not annotated (e.g., in older Cell Ranger versions where only the CDR3 sequence is annotated in 'filtered_contig_annotations.csv'), trimmed sequences are retrieved from the 'all_contig_annotations.json' file (if present)
    }else{
      # Check if the 'all_contig_annotations.json' file is present in the sample directory
      if(file.exists(paste0(VDJ.sample.directory, "all_contig_annotations.json"))){
        # Raw sequences are trimmed using the annotated 'L-REGION+V-REGION' and 'J-REGION' and stored in 'trimmed_contigs' dataframe
        trimmed_contigs <- get_trim_seqs_from_json(paste0(VDJ.sample.directory, "all_contig_annotations.json"))
        # Trimmed sequences are appended to dataframe in out[[1]] by merging by 'contig_id'
        out[[1]] <- out[[1]] |>
          dplyr::left_join(trimmed_contigs, by = c("contig_id" = "id")) |>
          dplyr::rename(sequence_nt_trimmed = sequence_trimmed)
        # Translate trimmed nucleotide sequences into protein sequences using the 'translate_DNA()' function
        out[[1]]$sequence_aa_trimmed <- sapply(out[[1]]$sequence_nt_trimmed, function(x) if(!is.na(x)){return(translate_DNA(x))}else{return(NA)})
      # If the 'all_contig_annotations.json' file is absent in the sample directory, a warning message is returned and the 'sequence_nt_trimmed' and 'sequence_aa_trimmed' columns are set to NA
      } else{
        warnings <- c(warnings, paste0("Could not find the 'all_contig_annotations.json' file in ", VDJ.sample.directory, " so 'sequence_nt_trimmed' and 'sequence_aa_trimmed' columns of sample ", basename(VDJ.sample.directory), "  remain empty."))
        out[[1]]$sequence_nt_trimmed <- NA
        out[[1]]$sequence_aa_trimmed <- NA
      }
    }

    # Read 'filtered.contig.fasta' file and store sequences in dataframe in out[[2]]
    # Since additional filtration steps are performed by enclone during the clonotype grouping step, some contigs present in 'filtered_contig.fasta' may be absent in 'filtered_contig_annotations.csv'
    out[[2]] <- data.frame(Biostrings::readDNAStringSet(paste(VDJ.sample.directory, "filtered_contig.fasta", sep = "")))
    # Rename column with raw nucleotide sequences of dataframe in out[[2]] to 'sequence_nt_raw'
    colnames(out[[2]]) <- 'sequence_nt_raw'
    # Add contig IDs from FASTA file in separate column in dataframe in out[[2]] and remove "-1" from and contig IDs
    out[[2]]$contig_id <- gsub(rownames(out[[2]]), pattern = "-1", replacement = "")

    # Read 'consensus_annotations.csv' file and store in out[[3]]
    out[[3]] <- utils::read.csv(paste0(VDJ.sample.directory, "consensus_annotations.csv"), sep = ",", header = TRUE)
    # Add "_aa" to column names with aa sequences for clarity
    colnames(out[[3]]) <- gsub(colnames(out[[3]]), pattern = "^(fwr\\d|cdr\\d)$", replacement = "\\1_aa")
    # Modify consensus IDs by placing an underscore between 'consensus' and the subsequent number for consistency
    out[[3]]$consensus_id <- unlist(lapply(out[[3]]$consensus_id, function(x) gsub(x, pattern = "(consensus)(.*)", replacement = "\\1_\\2")))
    # Remove double underscore (if present)
    out[[3]]$consensus_id <- unlist(lapply(out[[3]]$consensus_id, function(x) gsub(x, pattern = "__", replacement = "_")))
    # Remove "_" from 'v_gene', 'd_gene', 'j_gene', and 'c_gene' column
    colnames(out[[3]]) <- gsub(colnames(out[[3]]), pattern = "_gene", replacement = "gene")
    # Concatenate nt and aa sequences of the FWR1-4 and CDR1-3 regions and place concatenated sequence in new column in dataframe in out[[3]] (if the FWR1-4 and CDR1-3 regions are annotated in the 'consensus_annotations.csv' file)
    if(all(c("fwr1_nt","cdr1_nt","fwr2_nt","cdr2_nt","fwr3_nt","cdr3_nt","fwr4_nt","fwr1_aa","cdr1_aa","fwr2_aa","cdr2_aa","fwr3_aa","cdr3_aa","fwr4_aa") %in% colnames(out[[3]]))){
      out[[3]]$consensus_nt_trimmed <- paste0(out[[3]][["fwr1_nt"]], out[[3]][["cdr1_nt"]],
                                              out[[3]][["fwr2_nt"]], out[[3]][["cdr2_nt"]],
                                              out[[3]][["fwr3_nt"]], out[[3]][["cdr3_nt"]],
                                              out[[3]][["fwr4_nt"]])
      out[[3]]$consensus_aa_trimmed <- paste0(out[[3]][["fwr1_aa"]], out[[3]][["cdr1_aa"]],
                                              out[[3]][["fwr2_aa"]], out[[3]][["cdr2_aa"]],
                                              out[[3]][["fwr3_aa"]], out[[3]][["cdr3_aa"]],
                                              out[[3]][["fwr4_aa"]])
    # If the FWR1-4 and CDR1-3 regions are not annotated (e.g., in older Cell Ranger versions where only the CDR3 sequence is annotated in 'consensus_annotations.csv'), trimmed sequences are retrieved from the 'consensus_annotations.json' file (if present)
    }else{
      # Check if the 'consensus_annotations.json' file is present in the sample directory
      if(file.exists(paste0(VDJ.sample.directory, "consensus_annotations.json"))){
        # Raw sequences are trimmed using the annotated 'L-REGION+V-REGION' and 'J-REGION' and stored in 'trimmed_consensuses' dataframe
        trimmed_consensuses <- get_trim_seqs_from_json(paste0(VDJ.sample.directory, "consensus_annotations.json"))
        # Trimmed sequences are appended to dataframe in out[[3]] by merging by 'consensus_id'
        out[[3]] <- out[[3]] |>
          dplyr::left_join(trimmed_consensuses, by = c("consensus_id" = "id")) |>
          dplyr::rename(consensus_nt_trimmed = sequence_trimmed)
        # Translate trimmed nucleotide sequences and store protein sequences in the 'consensus_aa_trimmed' column
        out[[3]]$consensus_aa_trimmed <- sapply(out[[3]]$consensus_nt_trimmed, function(x) if(!is.na(x)){return(translate_DNA(x))}else{return(NA)})
        # If the 'consensus_annotations.json' file is absent in the sample directory, a warning message is returned and the 'consensus_nt_trimmed' and 'consensus_aa_trimmed' columns are set to NA
      } else{
        warnings <- c(warnings, paste0("Could not find the 'consensus_annotations.json' file in ", VDJ.sample.directory, " so 'consensus_nt_trimmed' and 'consensus_aa_trimmed' columns of sample ", basename(VDJ.sample.directory), " remain empty."))
        out[[3]]$consensus_nt_trimmed <- NA
        out[[3]]$consensus_aa_trimmed <- NA
      }
    }

    # Read 'consensus.fasta' file and store sequences in dataframe in out[[4]]
    out[[4]] <- data.frame(Biostrings::readDNAStringSet(paste0(VDJ.sample.directory, "consensus.fasta")))
    # Rename column with raw nucleotide consensus sequences of dataframe in out[[4]] to 'consensus_nt_raw'
    colnames(out[[4]]) <- 'consensus_nt_raw'
    # Add consensus IDs from the 'consensus.fasta' file in a separate column in the dataframe in out[[4]]
    out[[4]]$consensus_id <- rownames(out[[4]])
    
    # Read 'concat_ref.fasta' file and store sequences in dataframe in out[[5]]
    out[[5]] <- data.frame(Biostrings::readDNAStringSet(paste0(VDJ.sample.directory, "concat_ref.fasta")))
    # Rename column with raw nucleotide germline sequences in 'out[[5]]' to 'germline_nt_raw
    colnames(out[[5]]) <- 'germline_nt_raw'
    # Add consensus IDs from FASTA file in separate column in dataframe in out[[5]]
    out[[5]]$consensus_id <- rownames(out[[5]])
    # Modify raw consensus IDs by replacing '_concact_ref_' with '_consensus_'
    out[[5]]$consensus_id <- unlist(lapply(out[[5]]$consensus_id, function(x) gsub(x, pattern = "_concat_ref_", replacement = "_consensus_")))

    
    # 2. Merge data frames in list 'out' and rename items as follows:
    # [[1]] and [[2]] will be merged in 'filtered_contigs'
    # [[3]], [[4]], and [[5]] will be merged in 'consensuses'

    # Merge 'filtered_contig_annotations.csv' in out[[1]] with 'filtered_contig.fasta' in out[[2]] by 'contig_id'
    filtered_contigs <- out[[1]] |>
      dplyr::left_join(out[[2]], by = 'contig_id')
    # Merge 'consensus_annotations.csv' in out[[3]] with 'consensus.fasta' in out[[4]] and with 'concat_ref.fasta' in out[[5]] by 'consensus_id'
    consensuses <- out[[3]] |>
            dplyr::left_join(out[[4]], by="consensus_id") |>
            dplyr::left_join(out[[5]], by="consensus_id")
    
    
    # 3. Add 'germline_nt_trimmed' and 'germline_aa_trimmed' columns to 'consensuses' dataframe (if the 'trim.germlines' parameter is set to TRUE)

    # If 'trim.germlines' is set to TRUE, perform pairwise alignment with the raw germline sequences and the consensuses sequences using 'trim_seq()' function in parallel with 'mapply()' to obtain trimmed germline sequences
    if(trim.germlines){
      # Trim raw nucleotide germline sequences by performing an alignment of type 'global-local' with the trimmed nucleotide consensus sequence and store the trimmed sequences in the 'germline_nt_trimmed' column
      consensuses$germline_nt_trimmed <- mapply(trim_seq, seq1 = consensuses$consensus_nt_trimmed, seq2 = consensuses$germline_nt_raw, gap.opening.cost = gap.opening.cost, gap.extension.cost = gap.extension.cost, sequence.type = "DNA")
      # Translate the trimmed nucleotide germline sequences into protein sequences and store the translated germline sequences in the 'germline_aa_trimmed' column
      consensuses$germline_aa_trimmed <- sapply(consensuses$germline_nt_trimmed, function(x) if(!is.na(x)){return(translate_DNA(x))}else{return(NA)})
    }
    # If 'trim.germlines' is set to FALSE (default), the 'germline_nt_trimmed' and 'germline_aa_trimmed' columns are set to NA
    else{
      consensuses$germline_nt_trimmed <- NA
      consensuses$germline_aa_trimmed <- NA
    }
    
    
    # 4. Add 'germline_nt_predicted' and 'germline_aa_predicted' columns to 'consensuses' dataframe (if the 'fill.germline.CDR3' parameter is set to TRUE)
    
    # If 'fill.germline.CDR3' is set to TRUE, perform pairwise alignment with the trimmed germline sequences and the most frequently observed CDR3 sequences using the 'trim_seq()' function in parallel with 'sapply()' to obtain the CDR3 region in the germline 
    if(fill.germline.CDR3){
      consensuses$germline_nt_cdr3_filled <- sapply(consensuses$consensus_id, function(x){
        # Retrieve the most frequently observed sequence for the CDR3 region from the 'filtered_contigs' dataframe
        most_abundant_cdr3_seq <- names(table(filtered_contigs[filtered_contigs$consensus_id == x, "cdr3_nt"])[which.max(table(filtered_contigs[filtered_contigs$consensus_id == x, "cdr3_nt"]))])
        # Retrieve the trimmed germline sequence from the 'consensuses' dataframe using the consensus ID
        germline_nt_trimmed <- consensuses[consensuses$consensus_id == x, "germline_nt_trimmed"]
        # Perform pairwise alignment to select CDR3 region from 'germline_nt_trimmed' sequence using the 'most_abundant_cdr3_seq' as reference
        germline_cdr3 <- trim_seq(seq1 = most_abundant_cdr3_seq, seq2 = germline_nt_trimmed, gap.opening.cost = gap.opening.cost, gap.extension.cost = gap.extension.cost, sequence.type = "DNA")
        # If the selected germline CDR3 sequence ('germline_cdr3') occurs only once in the 'germline_nt_trimmed' sequence, this sequence is replaced by the most frequently observed CDR3 sequence ('most_abundant_cdr3_seq') to get the 'germline_nt_cdr3_filled' sequence
        if(stringr::str_count(germline_nt_trimmed, pattern = germline_cdr3) == 1){
          germline_nt_cdr3_filled <- gsub(pattern = germline_cdr3, replacement = most_abundant_cdr3_seq, germline_nt_trimmed)
        }
        # If the selected germline CDR3 sequence ('germline_cdr3') occurs more than once in the 'germline_nt_trimmed' sequence, the 'germline_nt_cdr3_filled' sequence is set to NA
        if(stringr::str_count(germline_nt_trimmed, pattern = germline_cdr3) != 1){
          germline_nt_cdr3_filled <- NA
        }
        # Return the 'germline_nt_cdr3_filled' sequence
        return(germline_nt_cdr3_filled)
      })
      # Translate the adapted germline sequences into protein sequences with the 'translate_DNA()' function and store the translated sequences in the 'germline_aa_cdr3_filled' column
      consensuses$germline_aa_cdr3_filled <- sapply(consensuses$germline_nt_cdr3_filled, function(x) if(!is.na(x)){translate_DNA(x)}else{return(NA)})
    }
    # If 'fill.germline.CDR3' is set to FALSE (default), the 'germline_nt_cdr3_filled' and 'germline_aa_cdr3_filled' columns are set to NA
    else{
      consensuses$germline_nt_cdr3_filled <- NA
      consensuses$germline_aa_cdr3_filled <- NA
    }
    
    
    # 4. Select the columns out of 'filtered_contigs' and 'consensuses' for each transcript and perform VDJ and VJ pairing
    # Simultaneously perform the filtering of cells/barcodes which contain multiple VDJ or VJ transcripts, if 'remove.divergent.cells' is set to TRUE (default)

    # Make selection of columns to be selected out of 'filtered_contigs' for each VDJ (heavy/beta chain) and VJ (light/alpha chain) transcript
    columns_contigs <- c("barcode",
                         "contig_id",
                         "clonotype_id",
                         "consensus_id",
                         "chain",
                         "umis",
                         "vgene","dgene","jgene","cgene",
                         "fwr1_nt","fwr1_aa",
                         "cdr1_nt","cdr1_aa",
                         "fwr2_nt","fwr2_aa",
                         "cdr2_nt","cdr2_aa",
                         "fwr3_nt","fwr3_aa",
                         "cdr3_nt","cdr3_aa",
                         "fwr4_nt","fwr4_aa",
                         "sequence_nt_raw",
                         "sequence_nt_trimmed",
                         "sequence_aa_trimmed")

    # Make selection of columns to be selected out of 'consensuses' for each VDJ (heavy/beta chain) and VJ (light/alpha chain) transcript
    columns_consensuses <- c("consensus_id",
                             "consensus_nt_raw",
                             "consensus_aa_raw",
                             "consensus_nt_trimmed",
                             "consensus_aa_trimmed",
                             "germline_nt_raw",
                             "germline_nt_trimmed",
                             "germline_aa_trimmed",
                             "germline_nt_cdr3_filled",
                             "germline_aa_cdr3_filled")

    # All column names
    all_columns <- c(columns_consensuses, columns_contigs)

    # If columns are not present in 'filtered_contigs', add them with NA values
    for(col in columns_contigs){
      if(!(col %in% colnames(filtered_contigs))){
        filtered_contigs[[col]] <- NA
      }
    }

    # If columns are not present in 'consensuses', add them with NA values
    for(col in columns_consensuses){
      if(!(col %in% colnames(consensuses))){
        consensuses[[col]] <- NA
      }
    }

    # Pre-merge the filtered contigs and consensuses (before splitting into the VDJ and VJ subsets)
    contigs_merged <- subset(filtered_contigs, select = columns_contigs) |>
      dplyr::left_join(subset(consensuses, select = columns_consensuses), by = "consensus_id")

    # Create a VDJ subset by selecting contigs encoding a TRB, TRG, or IGH (heavy) chain
    vdj_subset <- subset(contigs_merged, chain %in% c("TRB","TRG","IGH"))
    # Add "VDJ_" to column names of 'vdj_subset' (except the 'barcode' column)
    colnames(vdj_subset) <- ifelse(colnames(vdj_subset)  == "barcode", yes = colnames(vdj_subset), no = paste0("VDJ_", colnames(vdj_subset)))

    # Create a VJ subset by selecting contigs encoding a TRA, TRD, IGK, or IGL (light) chain
    vj_subset <- subset(contigs_merged, chain %in% c("TRA","TRD","IGK","IGL")) |>
      # Select all column except 'dgene' (absent in VJ transcripts)
      dplyr::select(all_columns[all_columns != "dgene"])
    # Add "VJ_" to column names of 'vj_subset' (except the 'barcode' column)
    colnames(vj_subset) <- ifelse(colnames(vj_subset) == "barcode", yes = colnames(vj_subset), no = paste0("VJ_", colnames(vj_subset)))

    # If 'remove.divergent.cells' is set to TRUE (default), filter out cells with multiple VDJ transcripts or multiple VJ transcripts by filtering out non-unique barcodes
    if(remove.divergent.cells){
      # Retrieve non-unique barcodes in 'vdj_subset' and 'vd_subset' and store in 'non_singlets_vdj' and 'non_singlets_vj' vector, respectively
      non_singlets_vdj <- unique(vdj_subset$barcode[duplicated(vdj_subset$barcode)])
      non_singlets_vj <- unique(vj_subset$barcode[duplicated(vj_subset$barcode)])
      non_singlets <- unique(c(non_singlets_vdj, non_singlets_vj))
      # Subset by the unique barcodes per VDJ and VJ dataframe
      vdj_subset_processed <- subset(vdj_subset, !(barcode %in% non_singlets))
      vj_subset_processed <- subset(vj_subset, !(barcode %in% non_singlets))
      # Save number of cells excluded due to suspicion of being divergent singlets in 'filtering_counts' dataframe
      filtering_counts[basename(VDJ.sample.directory), "divergent.cells"] <- length(non_singlets)
      }

    # If 'remove.divergent.cells' is set to FALSE, merge rows with identical 'barcode' in 'VDJ_subset' and VJ_subset'
    if(!remove.divergent.cells){
      # Concatenate non-missing values in each column into a comma-separated string for 'vdj_subset' and 'vj_subset'
      vdj_subset_processed <- stats::aggregate(data = vdj_subset, x = . ~ barcode, FUN = \(x)toString(stats::na.omit(x)), na.action = identity)
      vj_subset_processed <- stats::aggregate(data = vj_subset, x = . ~ barcode, FUN = \(x)toString(stats::na.omit(x)), na.action = identity)
    }
    
    # Merge 'vdj_subset' with 'vj_subset' by 'barcode'
    VDJ_df <- dplyr::full_join(vdj_subset_processed, vj_subset_processed, by = "barcode")
    
    # Add number of contigs per barcode for both VDJ and VJ transcript in columns 'VDJ_chain_count' and 'VJ_chain_count'
    VDJ_df$VDJ_chain_count <- sapply(VDJ_df$barcode, function(x) nrow(subset(vdj_subset, barcode == x)))
    VDJ_df$VJ_chain_count <- sapply(VDJ_df$barcode, function(x) nrow(subset(vj_subset, barcode == x)))
    
    # Combine 'VDJ_clonotype_id' and 'VDJ_clonotype_id'
    VDJ_df <- VDJ_df |>
      dplyr::mutate(clonotype_id = ifelse(!is.na(VDJ_clonotype_id) & !is.na(VJ_clonotype_id), 
                                          ifelse(VDJ_clonotype_id == VJ_clonotype_id, VDJ_clonotype_id, paste(VDJ_clonotype_id, VJ_clonotype_id, sep = ", ")),
                                          ifelse(is.na(VDJ_clonotype_id), VJ_clonotype_id, VDJ_clonotype_id))) |>
      dplyr::select(-VDJ_clonotype_id, -VJ_clonotype_id)
    
    # If 'complete.cells.only' is set to TRUE, filter out cells that have a 'VDJ_chain_count' or 'VJ_chain_count' of 0
    if(complete.cells.only){
      number_cells_before <- nrow(VDJ_df)
      VDJ_df <- subset(VDJ_df, VDJ_chain_count > 0 & VJ_chain_count > 0)
      number_cells_after <- nrow(VDJ_df)
      # Save number of cells excluded due to missing VDJ or VJ transcript in 'filtering_counts' dataframe
      filtering_counts[basename(VDJ.sample.directory), "incomplete.cells"] <- number_cells_before - number_cells_after
    }
    
    
    # 5. Finalize VDJ dataframe

    # Add column 'sample_id' to 'VDJ' dataframe and use name of sample directory folder
    VDJ_df$sample_id <- sample.ID
    
    # Remove duplicated values from 'clonotype_id' column
    VDJ_df$clonotype_id <- sapply(VDJ_df$clonotype_id, function(x) unique(strsplit(x, split = ", ")[[1]]))

    # Add column 'celltype' based on chains in 'VDJ_chain' and 'VJ_chain' column ("TRA", "TRB", "TRD", and TRG" --> "T cell"; "IGH", "IGK", and "IGL" --> "B cell")
    VDJ_df$celltype[stringr::str_detect(paste0(VDJ_df$VDJ_chain,VDJ_df$VJ_chain), "TR")] <- "T cell"
    VDJ_df$celltype[stringr::str_detect(paste0(VDJ_df$VDJ_chain,VDJ_df$VJ_chain), "IG")] <- "B cell"
    
    # Add the 'isotype' column based on the 'VDJ_cgene' column 
    VDJ_df <- VDJ_df |>
      dplyr::mutate(
        isotype = ifelse(celltype == "B cell",
                         paste(sub("^IGH([A-Z])([0-9]*)([A-Z]*)$", "Ig\\1\\2", strsplit(VDJ_cgene, split = ", ")[[1]]), collapse = ", "),
                         NA))

    # Add column 'clonotype_frequencies'
    VDJ_df$clonotype_frequency <- unlist(lapply(VDJ_df$clonotype_id, function(x) length(which(VDJ_df$clonotype_id == x))))

    # Specify order of columns in 'VDJ_df'
    new_order <- c("sample_id", "barcode", "celltype", "isotype",
                   "VDJ_contig_id", "VJ_contig_id",
                   "VDJ_consensus_id", "VJ_consensus_id",
                   "clonotype_id", "clonotype_frequency",

                   "VDJ_chain", "VDJ_chain_count", "VDJ_umis",
                   "VDJ_vgene", "VDJ_dgene", "VDJ_jgene", "VDJ_cgene",
                   "VDJ_fwr1_nt", "VDJ_fwr1_aa",
                   "VDJ_cdr1_nt", "VDJ_cdr1_aa",
                   "VDJ_fwr2_nt", "VDJ_fwr2_aa",
                   "VDJ_cdr2_nt", "VDJ_cdr2_aa",
                   "VDJ_fwr3_nt", "VDJ_fwr3_aa",
                   "VDJ_cdr3_nt", "VDJ_cdr3_aa",
                   "VDJ_fwr4_nt", "VDJ_fwr4_aa",
                   "VDJ_sequence_nt_raw", "VDJ_sequence_nt_trimmed", "VDJ_sequence_aa_trimmed",
                   "VDJ_consensus_nt_raw", "VDJ_consensus_nt_trimmed", "VDJ_consensus_aa_trimmed",
                   "VDJ_germline_nt_raw", "VDJ_germline_nt_trimmed", "VDJ_germline_aa_trimmed", "VDJ_germline_nt_cdr3_filled", "VDJ_germline_aa_cdr3_filled",

                   "VJ_chain", "VJ_chain_count", "VJ_umis",
                   "VJ_vgene", "VJ_jgene", "VJ_cgene",
                   "VJ_fwr1_nt", "VJ_fwr1_aa",
                   "VJ_cdr1_nt", "VJ_cdr1_aa",
                   "VJ_fwr2_nt", "VJ_fwr2_aa",
                   "VJ_cdr2_nt", "VJ_cdr2_aa",
                   "VJ_fwr3_nt", "VJ_fwr3_aa",
                   "VJ_cdr3_nt", "VJ_cdr3_aa",
                   "VJ_fwr4_nt", "VJ_fwr4_aa",
                   "VJ_sequence_nt_raw", "VJ_sequence_nt_trimmed", "VJ_sequence_aa_trimmed",
                   "VJ_consensus_nt_raw", "VJ_consensus_nt_trimmed", "VJ_consensus_aa_trimmed",
                   "VJ_germline_nt_raw", "VJ_germline_nt_trimmed", "VJ_germline_aa_trimmed", "VJ_germline_nt_cdr3_filled", "VJ_germline_aa_cdr3_filled")
    VDJ_df <- VDJ_df[, new_order]

    # Reorder rows in 'VDJ_df' dataframe by 'clonotype_id'
    VDJ_df <- VDJ_df |>
      dplyr::arrange(as.integer(stringr::str_extract(clonotype_id, "\\d+")))
    
    # Combine 'VDJ_df' dataframe, 'warnings' list, and 'filtering_counts' dataframe in 'output' list
    output <- list(VDJ_df = VDJ_df, warnings = warnings, filtering = filtering_counts)
    
    # Return 'output' list
    return(output)
  }
  
  
  # Define partial function for each sample (parallelized or non-parallelized):
  partial_function <- function(sample.ID){
    
    # Execute make_VDJ_sample() function
    make_VDJ_sample(sample.ID,
                    VDJ.sample.list = VDJ.sample.list,
                    remove.divergent.cells = remove.divergent.cells,
                    complete.cells.only = complete.cells.only,
                    trim.germlines = trim.germlines,
                    gap.opening.cost = gap.opening.cost,
                    gap.extension.cost = gap.opening.cost,
                    fill.germline.CDR3 = fill.germline.CDR3)
  }

  # If 'parallel' is set to TRUE, the VDJ dataframe will be constructed for each sample in parallel (if the number of samples is not exceeding the number of cores specified or available)
  if(parallel){

    # Retrieve the operating system
    operating_system <- Sys.info()[['sysname']]

    # If the operating system is Linux or Darwin, 'mclapply' is used for parallelization
    if(operating_system %in% c('Linux', 'Darwin')){
      # Construct a VDJ dataframe for each sample output folder in 'VDJ.directory' and store output in the list 'output_list'
      output_list <- parallel::mclapply(mc.cores = num.cores, X = names(VDJ.sample.list), FUN = partial_function)
      # Name items in 'output_list' according to sample names in 'VDJ.sample.list'
      names(output_list) <- basename(VDJ.sample.list)
      }

    # If the operating system is Windows, "parLapply" is used for parallelization
    if(operating_system == "Windows"){
      # Create cluster
      cluster <- parallel::makeCluster(num.cores)
      # Construct a VDJ dataframe for each sample output folder in 'VDJ.directory' and store output in the list 'output_list'
      output_list <- parallel::parLapply(cluster, X = names(VDJ.sample.list), fun = partial_function)
      # Name items in 'output_list' according to sample names in 'VDJ.sample.list'
      names(output_list) <- basename(VDJ.sample.list)
      # Stop cluster
      parallel::stopCluster(cluster)
      }
  }

  # If 'parallel' is set to FALSE, the per-sample VDJ building is not executed in parallel
  if(!parallel){
    # Construct a VDJ dataframe for each sample output folder in 'VDJ.directory' and store output in the list 'output_list'
    output_list <- lapply(names(VDJ.sample.list), partial_function)
    # Name items in 'output_list' according to sample names in 'VDJ.sample.list'
    names(output_list) <- basename(VDJ.sample.list)
    }

  # Select and print warning messages in 'output_list'
  warnings_output <- unlist(sapply(output_list, function(x) x$warnings))
  for(warning in warnings_output){message(paste(c("WARNING: ", i, "\n"), collapse = ""))}
  
  # Select and row bind filtering counts in 'output_list' 
  filtering_output <- do.call(rbind, lapply(output_list, function(x) x$filtering))
  # Print 'filtering_output' to show number of cells/barcodes excluded (if 'remove.divergent.cells' and/or 'complete.cells.only' is/are set to TRUE)
  message("WARNING: Please be aware of the number of cells excluded, when 'remove.divergent.cells' or 'complete.cells.only' is set to TRUE:\n")
  print(knitr::kable(filtering_output))
  
  # Select and row bind VDJ dataframes in 'output_list'
  VDJ_output <- do.call(rbind, lapply(output_list, function(x) x$VDJ_df))
  
  # Remove rownames from 'VDJ_output'
  rownames(VDJ_output) <- NULL
  
  # Ensure every missing value and every empty string are set to NA
  VDJ_output[VDJ_output == 'None'] <- NA
  VDJ_output[VDJ_output == 'NA'] <- NA
  VDJ_output[VDJ_output == 'NULL'] <- NA
  VDJ_output[VDJ_output == ''] <- NA
  VDJ_output[is.null(VDJ_output)] <- NA
  
  # Return 'VDJ_output'
  return(VDJ_output)
}