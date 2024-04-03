#' Function to convert VDJ dataframe into an AIRR-formatted TSV file.
#' Authors: Valentijn Tromp, Daphne van Ginneken
#' @description  Takes a VDJ dataframe, typically obtained from the 'minimal_VDJ()' function in Platypus, along with the imported IgBLAST annotations and alignments, as obtained from the 'import_IgBLAST_annotations()' function in Platypus, and converts it into a tab-separated values (TSV) file formatted according to the AIRR (Adaptive Immune Receptor Repertoire) guidelines.
#' @param VDJ dataframe - VDJ object as obtained from the 'minimal_VDJ()' function in Platypus, together with the imported IgBLAST annotations and alignments, as obtained from the 'import_IgBLAST_annotations' function in Platypus.
#' @param include list - a nested list specifying the samples and their associated clonotypes to include in the output TSV file. Each sublist represents a sample, where the sublist name is the sample name and the elements within the sublist are the clonotypes of that sample. If not provided, all samples and clonotypes are included.
#' @param columns list - a list specifying the columns to include in the output TSV file. At minimum, the following columns must be specified: 'sequence_id', 'clone_id', 'sequence', 'sequence_alignment', 'germline_alignment', 'v_call', 'v_sequence_start', 'v_sequence_end', 'v_germline_start', 'v_germline_end', 'j_call', 'j_sequence_start', 'j_sequence_end', 'j_germline_start', and 'j_germline_end'. The items in this list should correspond to the column names in the VDJ dataframe, while the names of the items in this list should refer to the column names of the output TSV file.
#' @param complete.rows.only bool - if TRUE, only complete rows (without any missing values) are included in the output TSV file. If FALSE, rows with missing values are retained in the output. Defaults to TRUE.
#' @param filter.rows.with.stop.codons bool - if TRUE, rows containing sequences with stop codons (TAA, TAG, TGA) in the 'sequence_alignment' and 'germline_alignment' columns are filtered out from the output TSV file. Defaults to TRUE.
#' @param output.file string - string specifying the path to the output file. If no path is specified, the output is written to 'airr_rearrengement.tsv' in the current working directory.
#' @return None
#' @examples
#' \dontrun{
#' VDJ_to_AIRR(VDJ)
#'}


VDJ_to_AIRR <- function(VDJ,
                        include,
                        columns,
                        complete.rows.only,
                        filter.rows.with.stop.codons,
                        output.file){
  
  # If no 'VDJ' dataframe is provided, a message is returned and execution is stopped
  if(missing(VDJ)){stop("ERROR: Please provide a VDJ dataframe as obtained from the 'minimal_VDJ()' function in Platypus.")}
  # If the 'include' parameter is not specified, all samples and clonotypes are appended to the 'include' list
  if(missing(include)){include <- lapply(unique(VDJ$sample_id), function(x) unique(VDJ[VDJ$sample_id == x, "clonotype_id"])); names(include) <- unique(VDJ$sample_id)}
  # Set default columns for each required column in the AIRR-formatted TSV file
  default_columns <- list(sequence_id = "barcode", 
                          clone_id = "clonotype_id",
                          sequence = "VDJ_sequence_nt_raw", 
                          sequence_alignment = "VDJ_sequence_alignment_nt_IgBLAST", 
                          germline_alignment = "VDJ_germline_alignment_nt_IgBLAST",
                          v_call = "VDJ_vgene_IgBLAST",
                          v_sequence_start = "VDJ_vgene_start_IgBLAST",
                          v_sequence_end = "VDJ_vgene_end_IgBLAST",
                          v_germline_start = "VDJ_vgene_germline_start_IgBLAST",
                          v_germline_end = "VDJ_vgene_germline_end_IgBLAST",
                          j_call = "VDJ_jgene_IgBLAST",
                          j_sequence_start = "VDJ_jgene_start_IgBLAST",
                          j_sequence_end = "VDJ_jgene_end_IgBLAST",
                          j_germline_start = "VDJ_jgene_germline_start_IgBLAST",
                          j_germline_end = "VDJ_jgene_germline_end_IgBLAST")
  # If the 'column' parameter is not specified, the defaults columns will be selected, and if the default columns could not be found, columns without the suffix '_IgBLAST' are selected (if present)
  if(missing(columns)){columns <- default_columns; for(i in names(columns)){if(!columns[[i]] %in% colnames(VDJ)){if(gsub(pattern = "_IgBLAST", replacement = "", columns[[i]]) %in% colnames(VDJ)){columns[[i]] <- gsub(pattern = "_IgBLAST", replacement = "", columns[[i]])}}}}
  # Check whether all necessary columns are specified in the 'columns' list and whether these columns are also present in the VDJ dataframe
  if(!all(c("sequence_id", "clone_id", "sequence", "sequence_alignment", "germline_alignment", "v_call", "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end", "j_call", "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end") %in% names(columns)) | !all(columns %in% colnames(VDJ))){stop("ERROR: Not all necessary columns could be found in the VDJ dataframe.")}
  # If the 'complete.rows.only' parameter is missing, it is set to TRUE
  if(missing(complete.rows.only)){complete.rows.only <- TRUE}
  # If the 'filter.rows.with.stop.codons' parameter is missing, it is set to TRUE
  if(missing(filter.rows.with.stop.codons)){filter.rows.with.stop.codons <- TRUE}
  # If the 'output' parameter is not specified, the TSV is written to the working directory as 'airr_rearrangement.tsv'
  if(missing(output.file)){output.file <- paste0(getwd(), "/airr_rearrangement.tsv")}
  
  # Make a subset of the VDJ dataframe with the specified samples and clonotypes in the 'include' parameter
  VDJ_subset <- do.call(rbind, lapply(names(include), function(x) VDJ[VDJ$sample_id == x & VDJ$clonotype_id %in% include[[x]], ]))
  
  # If there are multiple samples present in the 'VDJ_subset' dataframe, the sample ID is pasted as a prefix to the barcodes and the clonotype IDs
  if(length(unique(VDJ_subset$sample_id)) > 1){VDJ_subset$barcode <- paste(VDJ_subset$sample_id, VDJ_subset$barcode, sep = "_"); VDJ_subset$clonotype_id <- paste(VDJ_subset$sample_id, VDJ_subset$clonotype_id, sep = "_")}
  
  # Create the AIRR-formatted dataframe
  output_tsv <- VDJ_subset[, unlist(columns)]
  colnames(output_tsv) <- names(columns)
  
  # If the 'complete.rows.only' parameter is set to TRUE, remove rows with NA values in the 'output_tsv' dataframe
  if(complete.rows.only){output_tsv <- output_tsv[complete.cases(output_tsv), ]}
  
  # If the 'filter.rows.with.stop.codons' parameter is set to TRUE, replace sequences with stop codons (TAA/TAG/TGA) in the 'sequence_alignment' and the 'germline_alignment' column with 'STOP'
  if(filter.rows.with.stop.codons){
    for(alignment in c("sequence_alignment", "germline_alignment")){
      output_tsv[, alignment] <- unlist(lapply(output_tsv[, alignment], function(seq){
        codons <- sapply(1:(nchar(seq)/3), function(x) substr(seq, (x*3-2), (x*3)))
        if("TAA" %in% codons | "TAG" %in% codons | "TGA" %in% codons){return("STOP")}
        if(!("TAA" %in% codons) && !("TAG" %in% codons) && !("TGA" %in% codons)){return(seq)}
      }))
    }
  }
  
  # Remove rows with 'STOP' values in the 'sequence_alignment' column from the 'output_tsv' dataframe
  output_tsv <- output_tsv[output_tsv$sequence_alignment != "STOP", ]
  
  # Iterate through the clones present in the 'output_tsv' dataframe
  for(clone in unique(output_tsv$clone_id)){
    # Iterate through the columns containing germline annotations
    for(column in c("germline_alignment", "v_germline_start", "v_germline_end", "j_germline_start", "j_germline_end")){
      # If no germline alignment is present anymore due to the filtering of sequences containing stop codons, remove the clone from the 'output_tsv' dataframe
      if(nrow(output_tsv[output_tsv$clone_id == clone & output_tsv$germline_alignment != "STOP", ]) == 0){
        output_tsv <- output_tsv[output_tsv$clone_id != clone,]
      }
      # If different germline annotations are present within the same clonotype, replace them with the most abundant annotation, of the rows that do not contain 'STOP' in the 'germline_alignment'
      if(nrow(output_tsv[output_tsv$clone_id == clone & output_tsv$germline_alignment != "STOP", ]) != 0){
        output_tsv[output_tsv$clone_id == clone, column] <- names(table(output_tsv[output_tsv$clone_id == clone & output_tsv$germline_alignment != "STOP", column]))[which.max(table(output_tsv[output_tsv$clone_id == clone & output_tsv$germline_alignment != "STOP", column]))]
      }
    }
  }
  
  # Write the dataframe to the specified 'output.file'
  write.table(output_tsv, file = output.file, sep = "\t", row.names = FALSE, quote = FALSE)
}