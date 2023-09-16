#' S4 class for the Steropodon object.
#'
#' @slot structure bio3d pdb object of main modelled structure.
#' @slot sequence string - structure's corresponding amino acid sequence.
#' @slot barcodes vector of strings - cell barcodes corresponding to a given sequence/structure.
#' @slot pdbs bio3d pdbs object of aligned sequences and structures.
#' @slot H bio3d pdb object - heavy chain structure.
#' @slot L bio3d pdb object - light chain structure.
#' @slot CDRH3 bio3d pdb object - CDRH3.
#' @slot CDRL3 bio3d pdb object - CDRL3
#' @slot antigen bio3d pdb object - antigen
#' @slot complex bio3d pdb object - antibody-antigen complex from Steropodon_dock
#' @slot properties dataframe - structure physicochemical properties.
#' @slot paratope bio3d pdb object - paratope of complex from Steropodon_interface
#' @slot epitope bio3d pdb object - epitope of complex from Steropodon_interface
#' @slot core bio3d pdb object - invariant core from Steropodon_find_core
#' @slot structure_id string - unique structure id.
methods::setClass('Steropodon',
    slots = c(
      structure = 'ANY',
      sequence = 'ANY',
      barcodes = 'ANY',
      pdbs = 'ANY',
      H = 'ANY',
      L = 'ANY',
      CDRH3 = 'ANY',
      CDRL3 = 'ANY',
      antigen = 'ANY',
      complex = 'ANY',
      properties = 'ANY',
      paratope = 'ANY',
      epitope = 'ANY',
      core = 'ANY',
      structure_id = 'ANY'
    )
)
