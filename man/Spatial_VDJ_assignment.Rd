% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spatial_VDJ_assignment.R
\name{Spatial_VDJ_assignment}
\alias{Spatial_VDJ_assignment}
\title{Assign simulated immune repertoire sequences (BCR or TCR) simulated by Echidna to transcriptome and location in a spatial image in function of cell type.}
\usage{
Spatial_VDJ_assignment(
  GEX_matrix,
  vgm,
  vgm_VDJ,
  celltype,
  simulated_VDJ,
  method = c("random", "density", "germline")
)
}
\arguments{
\item{GEX_matrix}{Dataframe containing barcode, imagecol and imagerow from bcs_merge.}

\item{vgm}{Output of VDJ_GEX_matrix function with already the simulated VDJ data.}

\item{vgm_VDJ}{Dataframe from VDJ_GEX_matrix output (vgm[[1]]).}

\item{celltype}{Character designating the cell type that we want to study either "B" or "T".}

\item{simulated_VDJ}{Large list, output of Echidna simulate_repertoire function. Only needed if we use simulated data.}

\item{method}{Character to chose the assignment method of BCR or TCR to transcriptomic information, it can be "random", "density" or "germline".}
}
\value{
A dataframe corresponding to the VDJ (VGM[[1]]) with GEX_barcode and x, y coordinates column (allowing to localise each BCR or TCR on the spatial image).
}
\description{
Assign simulated immune repertoire sequences (BCR or TCR) simulated by Echidna to transcriptome and location in a spatial image in function of cell type.
}
\examples{
\dontrun{
#1)Assignment random to GEX
random_BCR_assignment <- Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,
vgm = vgm_with_simulated_VDJ,
vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B",
simulated_VDJ = simulated_B_cells_VDJ, method = "random")

#2)Assignment density-based
density_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,
vgm = vgm_with_simulated_VDJ,
vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B",
simulated_VDJ = simulated_B_cells_VDJ, method = "density")

#3)Assignment germline-based
germline_BCR_assignment<-Spatial_VDJ_assignment(GEX_matrix = GEX_matrix,
vgm = vgm_with_simulated_VDJ,
vgm_VDJ = vgm_with_simulated_VDJ$VDJ, celltype = "B",
simulated_VDJ = simulated_B_cells_VDJ, method = "germline")
}
}
