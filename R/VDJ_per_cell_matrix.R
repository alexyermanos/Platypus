#' Processes and organizes the repertoire sequening data from cellranger vdj and returns a single matrix where the rows are individual cells and the columns are repertoire features.
#' @param VDJ.out.directory Character vector with each element containing the path to the output of cellranger vdj runs. Multiple repertoires to be integrated in a single transcriptome should be supplied as multiple elements of the character vector. This can be left blank if supplying the clonotypes and contig files directly as input. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder. This is compatible with B and T cell repertoires (both separately and simultaneously).
#' @return returns a single matrix where the rows are individual cells and the columns are repertoire features.
#' @export
#' @examples
#' \dontrun{
#' example.vdj.analyze <- VDJ_per_cell_matrix(VDJ.out.directory = "~/path/to/cellranger/vdj/outs/")
#' }
#'
VDJ_per_cell_matrix <- function(VDJ.out.directory) {
  require(seqinr)
  require(jsonlite)
  require(tidyverse)
  require(msa)
  require(stringr)
  require(Biostrings)
  require(parallel)
  if(.Platform$OS.type == "unix") options(mc.cores=detectCores()-1, print("Using parallel package"))
  else options(mc.cores=1)

  # Accounting for directory name with or without "/"
  if (substr(VDJ.out.directory, nchar(VDJ.out.directory), nchar(VDJ.out.directory)) == "/") {slash <- ""} else {slash <- "/"}
  vdj <- VDJ_analyze(VDJ.out.directory)[[1]]

  print("Reading in input files")
  VDJ.out.directory_reference <- paste0(VDJ.out.directory, slash, "concat_ref.fasta")
  VDJ.out.directory_annotations <- paste0(VDJ.out.directory, slash, "all_contig_annotations.json")
  VDJ.out.directory_contigs <- paste0(VDJ.out.directory, slash, "filtered_contig_annotations.csv")

  contig.list <- read.csv(VDJ.out.directory_contigs, stringsAsFactors = FALSE,sep=",",header=T)
  reference.list <- read.fasta(VDJ.out.directory_reference, as.string = T,seqonly = F,forceDNAtolower = F)
  annotations.list <- read_json(VDJ.out.directory_annotations)

  # references
  references_HC <- list()
  for (x in vdj$clonotype_id) {
    temp_index <- gsub("_consensus_", "_concat_ref_",
                       contig.list$raw_consensus_id[contig.list$raw_clonotype_id == x & contig.list$chain == "IGH" & contig.list$is_cell == "True"])
    references_HC[[x]] <- as.character(reference.list[unique(temp_index[temp_index != "None"])])
    if (length(references_HC[[x]]) == 0) {references_HC[[x]] <- ""}
  }

  references_LC <- list()
  for (x in vdj$clonotype_id) {
    temp_index <- gsub("_consensus_", "_concat_ref_",
                       contig.list$raw_consensus_id[contig.list$raw_clonotype_id == x & contig.list$chain == "IGK" & contig.list$is_cell == "True"])
    references_LC[[x]] <- as.character(reference.list[unique(temp_index[temp_index != "None"])])
    if (length(references_LC[[x]]) == 0) {references_LC[[x]] <- ""}
  }

  # Barcodes
  barcode <- list()
  for (x in vdj$clonotype_id) {
    barcode[[x]] <- contig.list$barcode[contig.list$raw_clonotype_id==x]
  }
  # Selecting only sequences with start and end annotations
  filtered <- mclapply(annotations.list, function(x)
    any(sapply(x$annotations, function(x) x$feature$region_type=="L-REGION+V-REGION"))
    & any(sapply(x$annotations, function(x) x$feature$region_type=="J-REGION")))
  annotations.list <- annotations.list[unlist(filtered)]

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
               data.frame(contig_id = contig_id, chain = chain, sequence = substr(sequence, temp_start+1, temp_end-1))
             }
    )
  }

  seq_trimmed <- extract_trimmed(seq_per_clonotype(annotations.list, barcode))

  # Aligns trimmed sequences against reference to get trimmed reference sequence
  get_alinged_ref <- function(trimmed, reference, chain) {
    if (reference == "") {return(data.frame(""))} else {
      alignments <- pairwiseAlignment(trimmed$sequence[trimmed$chain == chain], reference, type = "local")
      data.frame(trimmed, ref_trimmed = as.character(subject(alignments[which.max(score(alignments))])))
    }
  }
  print("Alignment of trimmed sequences to get trimmed germline. This will take ~2 minutes")
  ref_trimmed_HC <- mcmapply(get_alinged_ref, seq_trimmed, references_HC, MoreArgs = list("IGH"), SIMPLIFY = F)
  ref_trimmed_LC <- mcmapply(get_alinged_ref, seq_trimmed, references_LC, MoreArgs = list("IGK"), SIMPLIFY = F)

  temp_contig <- contig.list %>%
    filter(is_cell == "True", productive == "True", high_confidence == "True") %>%
    select(!c(is_cell, productive, high_confidence, full_length)) %>%
    rename(., "cdr3"="cdr3_aa","raw_clonotype_id" = "clonotype_id")

  temp_contig <- temp_contig %>% inner_join(., select(bind_rows(seq_trimmed), contig_id, sequence), by="contig_id", suffix=c("", "_trimmed")) %>%
    inner_join(., select(bind_rows(ref_trimmed_HC), contig_id, ref_trimmed), by="contig_id", suffix=c("","")) %>%
    inner_join(., select(bind_rows(ref_trimmed_LC), contig_id, ref_trimmed), by="contig_id", suffix=c("",""))

  HC <- filter(temp_contig, chain == "IGH") %>% select(!c(chain))
  LC <- filter(temp_contig, chain %in% c("IGK", "IGL")) %>% select(!c(c_gene, chain))

  VDJ.per.cell <- inner_join(HC, LC, by ="barcode", suffix = c("_HC", "_LC")) %>%

    relocate(., clonotype_id_HC, cdr3_nt_HC, cdr3_nt_LC, cdr3_aa_HC, cdr3_aa_LC, c_gene, barcode,
             sequence_HC, sequence_LC, ref_trimmed_HC, ref_trimmed_LC) %>%
    select(!c(contig_id_LC, contig_id_HC, clonotype_id_LC)) %>%
    rename("clonotype_id_HC" = "clonotype_id") %>%
    mutate(sequence_HC_aa = as.character(Biostrings::translate(DNAStringSet(sequence_HC)))) %>%
    mutate(sequence_LC_aa = as.character(Biostrings::translate(DNAStringSet(sequence_LC)))) %>%
    select(-matches("reads|length|raw"))

  return(VDJ.per.cell)
}
