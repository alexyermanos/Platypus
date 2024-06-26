% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VDJ_GEX_clonotyme.R
\name{VDJ_GEX_clonotyme}
\alias{VDJ_GEX_clonotyme}
\title{Pseudotime analysis for scRNA and repertoire sequencing datasets}
\usage{
VDJ_GEX_clonotyme(
  method,
  version,
  top.N.clonotypes,
  vdj.gex.matrix.output,
  vdj.analyze.output,
  gex.automate.output,
  exclude.clusters,
  colors,
  show.cells,
  highlight.genes,
  dropest.output.list,
  velocyto.gex.merged,
  velocyto.file.name,
  velocyto.out.dir,
  velocyto.save.rds,
  velocyto.norm.scale.factor,
  velocyto.n.variable.features,
  velocyto.neighbor.dim,
  velocyto.cluster.resolution,
  velocyto.mds.dim,
  velocyto.nCount_spliced,
  velocyto.percent.mt,
  velocyto.normalisation.method,
  velocyto.selection.method,
  velocyto.deltaT,
  velocyto.kCells,
  velocyto.fit.quantile,
  velocyto.kGenes,
  root.selection,
  root.marker,
  ridgeline.separator,
  genes.for.module.score,
  root.nodes,
  color.cells
)
}
\arguments{
\item{method}{Pseudotime analylysis method to be used. Possible parameters are monocle3 or velocyto. monocle3 is being used as a default method. For velocyto analysis please run on Cluster and it is only available for UNIX based systems.}

\item{version}{Platypus version to use "v2" or "v3". version 2 used by default.}

\item{top.N.clonotypes}{How many clonotypes to show per sample in the Ridgeline plots and on the Velocyto UMAP.}

\item{vdj.gex.matrix.output}{If Platypus v3 is used, the input to this function has to be the output of the VDJ_GEX_matrix function.}

\item{vdj.analyze.output}{If Platypus v2 is used, the VDJ_analyze output has to be supplied.}

\item{gex.automate.output}{If Platypus v2 is used, the GEX_automate output has to be supplied here.}

\item{exclude.clusters}{Please enter a cluster number if you'd like to exclude a certain cluster from analysis. Cluster will be assigned to different partition in Monocle3 analysis and therefore pseudotime distance will be set to infinity. Cells from this cluster will be deleted from the dataset in the Velocyto analysis.}

\item{colors}{Vector containing custom colors to be used for highlighting the clonotypes. If left empty, default colors will be assigned.}

\item{show.cells}{Logical, should cells be shown in the Ridgeline plots. True by default.}

\item{highlight.genes}{Vector containing gene names. The expressionlevels of these genes along pseudotime will be plotted.}

\item{dropest.output.list}{List containing the cell.counts.matrices.rds from the Dropest alignment for Velocyto analysis.}

\item{velocyto.gex.merged}{Logical whether samples should be shown in combined UMAP or sepeartely.}

\item{velocyto.file.name}{String used as file name when saving the output pdf}

\item{velocyto.out.dir}{Directory to save the output files. By default the current working directory.}

\item{velocyto.save.rds}{If RDS objects should be saved as well. Default = F.}

\item{velocyto.norm.scale.factor}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.n.variable.features}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.neighbor.dim}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.cluster.resolution}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.mds.dim}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.nCount_spliced}{Cutoff threshold. cells with less spliced gene counts will be omitted. Filtering of bad quality cells.}

\item{velocyto.percent.mt}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.normalisation.method}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.selection.method}{Parameter for GEX analysis of the cell.count.matrices.}

\item{velocyto.deltaT}{Parameter for Velocyto analysis}

\item{velocyto.kCells}{Parameter for Velocyto analysis}

\item{velocyto.fit.quantile}{Parameter for Velocyto analysis}

\item{velocyto.kGenes}{Parameter for Velocyto analysis}

\item{root.selection}{Character. Method for root selection. Defaults to "manual"}

\item{root.marker}{Character. Marker to use as Root. Defaults to "SELL"}

\item{ridgeline.separator}{Character. Variable to group ridgeline plots by. Defaults to "clonotype"}

\item{genes.for.module.score}{List of vectors of genes. With module scores inferred via Seurat::AddModuleScore(). Default is set to NULL.}

\item{root.nodes}{Labeled node from trajectory plot to specify root nodes root nodes for pseudotime trajectory.}

\item{color.cells}{For the module score plot decide how the cells should be coloured (based on e.g. group_id, seurat_clusters etc.). Default = 'seurat_clusters'}
}
\value{
If method=monocle3, the function returns a list element: [[1]] UMAP colored by Pseudotime, [[2]] Ridgeline plots showing the density of each of the top.N.clonotypes per cluster along pseudotime., [[3]] Gene expression plots highlighting the gene expression across pseudotime colored by transcriptional cluster, [[4]] Gene expression plots highlighting the gene expression across pseudotime colored by colotype. If method=velocyto, plots and RDS will be saved to velocyto.out.dir.
}
\description{
Pseudotime analysis for scRNA and repertoire sequencing datasets
}
\examples{
\dontrun{
#----Method=monocle3----

# Version 2
   vdj_repertoire_tcells <- VDJ_analyze(VDJ.out.directory =VDJ.out.directory.list,
    filter.1HC.1LC = T)
   gex_acute <- Platypus::GEX_automate(GEX.outs.directory.list = dir_to_gex[1:1],
    integration.method = "scale.data", mito.filter = 20, cluster.resolution = 0.5,
   VDJ.gene.filter = T)

clonotyme_output <- VDJ_GEX_clonotyme(vdj.analyze.output = vdj_repertoire_tcells,
 gex.analyze.output = gex_acute, version="v2", exclude.clusters=7, highlight.genes="sell",
  colors = c("blue", "red", "black", "orange")
clonotyme_output[[4]]

# Version 3
    VGM <-
    readRDS("C:/Users/rapha/Downloads/TEMPLATE_VDJ_GEX_mat_Bcells_r2_150521.rds")

clonotyme_output <- VDJ_GEX_clonotyme(vdj.gex.matrix.output = VGM, version="v3",
 highlight.genes="sell", top.N.clonotypes = 1)


#---Method=velocyto----

#Dropest Alignment: Run on EULER CLUSTER
#env2lmod
#module load gcc/4.8.5 python/3.7.4
#module load gcc/4.8.5 dropest/0.8.6
#module load gcc/4.8.5 r/4.0.2
#module load gcc/4.8.5 hdf5/1.10.1
#module load gcc/4.8.5 openmpi/4.0.2
#module load gcc/4.8.5 r/4.0.2
#bsub -W 2880 -R 'rusage[mem=20000]'
/cluster/home/rakuhn/dropEst/dropest -V -C 6000 -f -g
/cluster/scratch/rakuhn/mm10-2020-A/refdata-gex-mm10-2020-A/genes/genes.gtf
-c /cluster/home/rakuhn/dropEst/configs/10x.xml
/cluster/scratch/rakuhn/cellranger_v5/g1/outs/possorted_genome_bam.bam

#Load required VDJ.analyze.output on EULER CLUSTER

vdj_repertoire_tcells
<- readRDS("/cluster/home/rakuhn/RPII/vdj_repertoire_tcells.rds")
vdj_repertoire_tcells
<- head(vdj_repertoire_tcells,2)
#Only select the first two repertoires since we only want to analyze these two.

# Load the two corresponding Dropest cell.count.matrices.rds

dropest.output.list <- list()
dropest.output.list[[1]]
<- readRDS("/cluster/home/rakuhn/RPII/old_bam/gex1/cell.counts.matrices.rds")
dropest.output.list[[2]]
<- readRDS("/cluster/home/rakuhn/RPII/old_bam/gex2/cell.counts.matrices.rds")

# Run Velocyto using Clonotyme

VDJ_GEX_clonotyme(method = "velocyto", version = "v2",
 vdj.analyze.output = vdj_repertoire_tcells,
 dropest.output.list = dropest.output.list,
 top.N.clonotypes = 3, exclude.clusters = 8, highlight.genes = "sell",
 velocyto.gex.merged = T, velocyto.out.dir = ".", velocyto.save.rds = F)

}
}
