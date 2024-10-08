% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VDJ_call_enclone.R
\name{VDJ_call_enclone}
\alias{VDJ_call_enclone}
\title{(Re)clonotype a VDJ object using cellranger's enclone tool}
\usage{
VDJ_call_enclone(
  VDJ,
  VDJ.directory,
  global.clonotype,
  samples.to.clonotype,
  samples.to.combine,
  same.origin,
  output.format,
  operating.system,
  parallel
)
}
\arguments{
\item{VDJ}{VDJ or VDJ.GEX.matrix[[1]] object, as obtained from the VDJ_GEX_matrix function in Platypus.}

\item{VDJ.directory}{string - directory for the VDJ data, should be the main folder which includes the individual sample folders (each with the all_contig_annotations.json file that is used by enclone)}

\item{global.clonotype}{bool - if T, will use clonotype definitions irrespective of samples. Must also be T is you wish to merge clonotypes from two specific (which should be specified in the samples.to.combine parameter)}

\item{samples.to.clonotype}{- vector - lists the samples names which should be clonotyped. The unspecified samples will keep their old clonotype defintions.}

\item{samples.to.combine}{- vector or list of vectors - lists the samples which you wish to have their clonotypes merged (e.g., c('s1','s2') to only merge the first 2 samples, or list(c('s1','s3'), c('s2', 's4')) to merge the first and third, second and fourth, respectively). global.clonotype must be set to T!}

\item{same.origin}{bool - if the merged samples come from the same donor, with the same or with different origins. If two datasets come from the same origin, enclone will filter to remove certain artifacts.}

\item{output.format}{string - 'vgm' to output a VGM-specific VDJ dataframe (all samples in the same dataframe).}

\item{operating.system}{string - operating system on which enclone will be run. 'Windows' for Windows, 'Linux' for Linux, 'Darwin' for MacOS.}

\item{parallel}{bool - if T, the program will be executed in parallel, on no. cores = max. available cores - 1.}
}
\value{
Reclonotyped VDJ object using the enclone software and 10x-specific clonotype deifinition.
}
\description{
Calls recon to clonotype a VDJ object given a VDJ.directory (with sample folders which should include the all_contig_annotations.json file) - outputs a new VDJ with updated clonotype_id, clonotype_id_10x, and clonotype_frequency columns
}
\examples{
\dontrun{
VDJ_call_enclone(vdj, VDJ.directory, samples.to.combine = c('s1', 's2', 's3'), global.clonotype = T)
}
}
