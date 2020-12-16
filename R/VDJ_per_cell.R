#' Analyzes and processes the repertoire sequencing data from cellranger vdj. This provides information on the single-cell level for each clone, as opposed to the output from VDJ_analyze.
#' @title VDJ_per_cell
#' @param VDJ.directory Character vector with each element containing the path to the output of cellranger vdj runs. This corresponds to the same object used for the VDJ_analyze function. Multiple repertoires to be integrated in a single transcriptome should be supplied as multiple elements of the character vector. This can be left blank if supplying the clonotypes and contig files directly as input. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @param clonotype.list Output from either VDJ_analyze or VDJ_clonotype functions. This list should correspond to a single GEX.list object, in which each list element in clonotype.list is found in the GEX.object. Furthermore, the i'th entry in the directory supplied to GEX.list should correspond to the i'th element in the clonotype.list object.
#' @param contig.list List of dataframe based on the all_contigs.csv file from cellranger vdj output. If 10x sequencing was not used then this object should be formatted with the same columns as the 10x object.
#' @param fasta.list Contains the full-length sequence information in the same format as filtered_contig.fasta file from the output of cellranger.
#' @param reference.list Contains the reference sequence information in the same format as concat_ref.fasta file from the output of cellranger.
#' @param filtered.contigs Logical indicating if the filtered contigs file should be used. TRUE will read VDJ information from only the filtered output of cellranger. FALSE will read the all contigs file from cellranger. Default set to TRUE (filtered output)
#' @return Returns a list of dataframes containing
#' @export
#' @examples
#' \dontrun{
#' check_VDJ_per_cell <- VDJ_per_cell(clonotype.list = output.from.VDJ_analyze, VDJ.directory = "path/to/cellranger/outs/")
#' }
VDJ_per_cell <- function(clonotype.list,
                         VDJ.directory,
                         contig.list,
                         fasta.list,
                         reference.list,
                         filtered.contigs,
                         annotations.json) {
  require(seqinr)
  require(jsonlite)
  require(dplyr)
  require(msa)
  require(stringr)
  require(Biostrings)
  library(parallel)
  if(.Platform$OS.type == "unix") options(mc.cores=detectCores() -1 , print("Using parallel package")) else options(mc.cores=1)
  if(missing(VDJ.directory)) print("No output directory supplied. Assuming clonotype and contig are provided as list objects")
  if(missing(contig.list)) print("No contig.list supplied. Assuming contigs should be extracted from working directory")
  if(missing(filtered.contigs)) filtered.contigs <- TRUE

  # Extracting mouse_number from dir name
  if (all(str_detect(VDJ.directory, "mouse"))) {mouse_number <- substr(sub('.*mouse', '', VDJ.directory), 1, 1)}

  print("Reading in input files")

  # Reading files -----------------------------------------------------------

  if(missing(VDJ.directory)==F) {
    if (str_sub(VDJ.directory, -1) == "/") {slash <- ""} else {slash = "/"}
    if(filtered.contigs==T) VDJ.directory_contigs <- paste0(VDJ.directory, slash, "/filtered_contig_annotations.csv")
    if(filtered.contigs==F) VDJ.directory_contigs <- paste0(VDJ.directory, slash, "/all_contig_annotations.csv")
    VDJ.directory_fasta <- paste0(VDJ.directory, slash, "/filtered_contig.fasta")
    VDJ.directory_reference <- paste0(VDJ.directory, slash, "/concat_ref.fasta")
    VDJ.directory_annotations <- paste0(VDJ.directory, slash, "/all_contig_annotations.json")

    contig.list <- lapply(VDJ.directory_contigs, function(x) utils::read.csv(x, stringsAsFactors = FALSE,sep=",",header=T))
    fasta.list <- lapply(VDJ.directory_fasta, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    reference.list <- lapply(VDJ.directory_reference, function(x) seqinr::read.fasta(x, as.string = T,seqonly = F,forceDNAtolower = F))
    annotations.list <- lapply(VDJ.directory_annotations, read_json)
  }

  VDJ.per.cell <- list()

  # each element in VDJ.per.cell is a dataframe containing per cell  --------

  for(i in 1:length(clonotype.list)) {

  # Reference and barcodes --------------------------------------------------

    contig.list[[i]]$chain[which(contig.list[[i]]$chain == "TRB")] <- "IGH"
    contig.list[[i]]$chain[which(contig.list[[i]]$chain == "TRA")] <- "IGK"

    # references
    references_HC <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {
      temp_index <- gsub("_consensus_", "_concat_ref_",
                         contig.list[[i]]$raw_consensus_id[contig.list[[i]]$raw_clonotype_id == x
                       & contig.list[[i]]$chain == "IGH" & contig.list[[i]]$is_cell == "True"])

      references_HC[[x]] <- as.character(reference.list[[i]][unique(temp_index[temp_index != "None"])])
      if (length(references_HC[[x]]) == 0) {references_HC[[x]] <- ""}
    }

    references_LC <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {
      temp_index <- gsub("_consensus_", "_concat_ref_",
                         contig.list[[i]]$raw_consensus_id[contig.list[[i]]$raw_clonotype_id == x
                                                           & contig.list[[i]]$chain == "IGK" & contig.list[[i]]$is_cell == "True"])
      references_LC[[x]] <- as.character(reference.list[[i]][unique(temp_index[temp_index != "None"])])
      if (length(references_LC[[x]]) == 0) {references_LC[[x]] <- ""}
    }

    # Barcodes
    barcode <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {
      barcode[[x]] <- contig.list[[i]]$barcode[contig.list[[i]]$raw_clonotype_id==x]
    }


  # Trimmed sequences -------------------------------------------------------

    # Selecting only sequences with start and end annotations
    filtered <- mclapply(annotations.list[[i]], function(x)
                any(sapply(x$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION"))
              & any(sapply(x$annotations, function(x) x$feature$region_type=="J-REGION")))
    annotations.list[[i]] <- annotations.list[[i]][unlist(filtered)]

    # Returns list of annotations for each cell in each clonotype
    seq_per_clonotype <- function(annotations, barcodes) {
      annotations_barcodes <- unlist(sapply(annotations, function(x) x$barcode))
      mcmapply(function(x) annotations[annotations_barcodes %in% x], barcodes)
    }

    # Extracts trimmed sequence for each cell in the clonotype
    extract_trimmed <- function(seq_per_clono) {
      mclapply(seq_per_clono,
               function(clono) {
                 temp_start <- sapply(clono, function(x) x$annotations[sapply(x$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION")][[1]]$contig_match_start)
                 temp_end <- sapply(clono, function(x) x$annotations[sapply(x$annotations, function(x) x$feature$region_type=="J-REGION")][[1]]$contig_match_end)
                 sequence <- sapply(clono, function(x) x$sequence)
                 contig_id <-  sapply(clono, function(x) x$contig_name)
                 chain <- sapply(clono, function(x) x$annotations[[1]]$feature$chain)
                 data.frame(contig_id = contig_id, chain = chain, sequence = substr(sequence, temp_start + 1, temp_end))
               }
      )
    }

    seq_trimmed <- extract_trimmed(seq_per_clonotype(annotations.list[[i]], barcode))


  # Aligns trimmed sequences against reference to get trimmed reference --------

    get_aligned_ref <- function(trimmed, reference, chain) {
      if (length(reference) > 1) {reference = reference[which.min(str_length(reference))]}
      if (reference == "") {return("")}
      else {
        alignments <- pairwiseAlignment(trimmed$sequence[trimmed$chain == chain], reference, type = "local")
        as.character(subject(alignments[which.max(score(alignments))]))
      }
    }

    print("Alignment of trimmed sequences to get trimmed germline. This will take ~2 minutes")
    ref_trimmed_HC <- mcmapply(get_aligned_ref, seq_trimmed, references_HC, MoreArgs = list("IGH"))
    ref_trimmed_LC <- mcmapply(get_aligned_ref, seq_trimmed, references_LC, MoreArgs = list("IGK"))


  # Pasted cdr3 -------------------------------------------------------------

    pasted_cdr3 <- mclapply(barcode, function(x) {

      cdr3_HC <- filter(contig.list[[i]][contig.list[[i]]$barcode %in% x,],
                        is_cell=="True", productive=="True", chain=="IGH") %>% select (cdr3, barcode, chain)
      cdr3_LC <- filter(contig.list[[i]][contig.list[[i]]$barcode %in% x,],
                        is_cell=="True", productive=="True", chain %in% c("IGK", "IGL")) %>% select (cdr3, barcode, chain)

      cdr3 <- inner_join(cdr3_HC, cdr3_LC, by="barcode", suffix = c("_HC", "_LC"))

       data.frame(pasted_cdr3 = paste0(cdr3$chain_HC, ":", cdr3$cdr3_HC, ";", cdr3$chain_LC, ":",cdr3$cdr3_LC),
                        CDR3H = cdr3$cdr3_HC, CDR3L = cdr3$cdr3_LC, barcode = cdr3$barcode)
    })

  # Extracting gene usage features from contig.list -------------------------

    extract_feature <- function(barcodes, trimmed_sequence, trimmed_reference, feature, chain) {
      matching <- contig.list[[i]]$barcode%in%barcodes & contig.list[[i]]$chain == chain & contig.list[[i]]$is_cell == "True"
      if (trimmed_reference == "") {trimmed_reference <- rep("", length(contig.list[[i]]$barcode[matching]))}

      data <- data.frame(barcode = contig.list[[i]]$barcode[matching], contig.list[[i]][matching, feature],
                         full_seq = as.character(fasta.list[[i]][contig.list[[i]]$contig_id[matching]]),
                         trimmed_ref = trimmed_reference)
      data <- inner_join(data, as.data.frame(trimmed_sequence, col.names = c("contig_id", "chain", "sequence")), by = "contig_id")
    }

    data_HC <- mcMap(extract_feature, barcode, seq_trimmed, ref_trimmed_HC,
                     MoreArgs = list(c("umis", "contig_id", "c_gene", "v_gene", "j_gene", "raw_clonotype_id"), "IGH"))

    data_LC <- mcMap(extract_feature, barcode, seq_trimmed, ref_trimmed_LC,
                     MoreArgs = list(c("umis", "contig_id", "c_gene", "v_gene", "j_gene"), "IGK"))

    VDJ.per.cell[[i]] <- mcMap(inner_join, data_HC, data_LC, by="barcode", suffix=list(c("_HC", "_LC")))
    VDJ.per.cell[[i]] <- mcMap(left_join, VDJ.per.cell[[i]], pasted_cdr3, by="barcode", suffix=list(c("", "")))

  }
  return(VDJ.per.cell)
}
