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



methods::setMethod(f='show', signature='Steropodon',
 definition=function(object) {
   cat('Steropodon object', '\n')

   cat('Sample: ', object@structure_id$sample, '  ', 'Clonotype: ', object@structure_id$clonotype, '  ', 'Rank: ', object@structure_id$rank)
   cat('\n')

   cat('Number of cells: ', length(object@barcodes))
   cat('\n')

   structures <- c()

   if(!is.null(object@structure)){
     structures <- c(structures, 'main')
   }

   if(!is.null(object@complex)){
     structures <- c(structures, 'complex')
   }

   if(!is.null(object@antigen)){
     structures <- c(structures, 'antigen')
   }

   if(!is.null(object@paratope)){
     structures <- c(structures, 'paratope')
   }

   if(!is.null(object@epitope)){
     structures <- c(structures, 'epitope')
   }

   if(!is.null(object@pdbs)){
     structures <- c(structures, 'PDBs')
   }

   if(!is.null(object@core)){
     structures <- c(structures, 'core')
   }

   if(!is.null(object@properties)){
     structures <- c(structures, 'properties')
   }

   if(!is.null(object@H)){
     structures <- c(structures, 'H')
   }

   if(!is.null(object@L)){
     structures <- c(structures, 'L')
   }

   if(!is.null(object@CDRH3)){
     structures <- c(structures, 'CDRH3')
   }

   if(!is.null(object@CDRL3)){
     structures <- c(structures, 'CDRL3')
   }

   cat('Structures available: ', paste0(structures, collapse = ', '))

   cat('\n')

   cat('Sequence: ', object@sequence)

   cat('\n')

   cat('\n')

   cat('\n')

 }
)
